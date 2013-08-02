#!/usr/bin/env python

from optparse import OptionParser
import simplejson as json
import os
import sys
import logging
import traceback
from collections import namedtuple
import itertools
import pickle
import math
import random
from rtree import index
import multiprocessing
from multiprocessing import Pool, Manager, Lock, cpu_count
from functools import partial

logging.basicConfig(format="[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s")
logging.getLogger().setLevel(logging.INFO)

Location = namedtuple('Location', ['pos', 'locationid', 'pop'])

def make_spatial_index(rtree_file, data):

    logging.info('making spatial index')
    rtree = None

    if not os.path.exists(rtree_file + ".dat"):
        logging.info("making rtree {}".format(rtree_file))
        rtree = index.Rtree(rtree_file)
        i = 0
        for feat in data:
            # coords = feat.pos.coords
            coords = feat.pos
            rtree.insert(i, (coords[0], coords[1], coords[0], coords[1]), feat)
            i += 1
        logging.info('finished making rtree')
    else:
        logging.info('path exists %s' % (rtree_file))
        rtree = index.Rtree(rtree_file)

    return rtree

def make_spatial_index_memory(data):

    logging.info('making spatial index')

    rtree = index.Index()
    for i, feat in enumerate(data):
        coords = feat.pos
        rtree.insert(i, (coords[0], coords[1], coords[0], coords[1]), feat)

    logging.info('finished making rtree')
    return rtree


"""
The input is a file with 

    location_id, x, y, location_str, population
    
the id is [0,N-1] and the location_str can be anything.

The x,y are expected to be lat,lon with srid=4326, and are transformed to srid=900913 before being
added to the spatial index.
"""
def read_locations(fname):

    pkl = fname + ".pkl"
    if os.path.exists(pkl):
        data = pickle.load(open(pkl, 'r'))
    else:
        import django.contrib.gis.geos as geos

        data = []
        with open(fname, 'r') as f:
            for l in f.readlines():
                location, x, y, locationid, pop = l.strip().split(',')
                # we create a geos geometry object and transform it to 900913 - this allows us to do simple distance calculations in metres
                pos = geos.Point([float(y), float(x)])
                pos.srid = 4326
                pos.transform(900913)
                # data.append(Location(locationid=locationid, pos=pos, pop=float(pop)))
                data.append(Location(locationid=locationid, pos=pos.coords, pop=float(pop)))
        pickle.dump(data, open(pkl, 'w'), -1)
    return data


"""
starting a node: i, find the radiation probability t_ij to each other node j.
Computing t_ij for a pair i_j involves computing the sum s_ij of populations
in a circle centered on i. This is the part of the computation which takes all
the time, and I use a spatial index to speed it up. Of course, if you had enough
memory you could compute a distance matrix, but you have the cost of making the matrix
and its not clear it would be any faster in computing s_ij.

Rtree provides rtree.nearest(rect) which give a list of objects in the index nearest to the
current point (rect) and sorted in order of distance. I build this list just once, then loop
through it for each j, and stop when I encounter the object j, which is the boundary of circle
of radius d_ij centred on i. There are no distance computations involved then, and the 'nearest'
list is computed just once, so I think this is probably as fast as the method can be.

s_ij is not symmetric- each node i has a different environment, so the calculation must be done
for each ij pair in the network.

Most of the t_ij's will be extrememly small, so impose a cutoff- I don't yet have a good estimate
for what this cutoff should be.

I/O with multiprocessing is a littler trickier than I had thought, so each task creates its own
output file. This is actually better, since it avoids locking issues- at the end of the program
just cat the files together:

  cat t_ij_*.dat | sort -n -k 1 -k 2 > all_t_ij.dat
  
  
Any ideas welcome.
"""
def radiation_task(i, data, norm, min_tij, lock, counter):
    cp = multiprocessing.current_process()
    data_i = data[i]
    pos_i, pop_i = data_i.pos, data_i.pop
    logging.info('working on {} : {}% @{}'.format(i, 1.0 * counter.value / len(data), cp.name))
    
    out_file = open('t_ij_{}.dat'.format(i), 'w')
    # rtree doesn't seem to be threadsafe, so lock the access
    # with lock:
    
    # this returns a list of objects in order of distance to pos_i
    hits_i = [h for h in idx.nearest((pos_i[0], pos_i[1], pos_i[0], pos_i[1]), len(data), objects=False)]

    # get t_ij for each other j
    for j, data_j in enumerate(data):
        if i == j: continue            
        pos_j, pop_j = data_j.pos, data_j.pop
        # distance ordered iterator to objects nearest to location 'i'- only need to consider objects closer than d_ij
        s_ij = 0
        for hit in hits_i:
            if hit == j:
                break
            else:
                s_ij += data[hit].pop

        # the spatial index returns the locations i as the first element so exclude it from count
        # s_ij = s_ij - (pop_i + pop_j)
        s_ij -= pop_i
        mob_i_j = (1.0 * pop_i * norm)
            
        # logging.info('s_ij: {} {} -> {}'.format(pop_i, pop_j, s_ij))
        t_ij = (1.0 * mob_i_j * pop_i * pop_j) / ((pop_i + s_ij) * (pop_i + pop_j + s_ij))            
        if t_ij > min_tij:
            # with lock:
            out_file.write("%d\t%d\t%d\t%d\t%f\t%f\n" % (i, j, pop_i, pop_j, s_ij, t_ij))
                # logging.info("i->j: {} -> {} t_ij({})".format(i, j, t_ij))

    counter.value += 1
    return None

# have to make this a global to prevent it getting garbage collected- could be a bug with partial stuff
# TODO:: also figure out how to do I/O with multiprocessing
idx = None

def radiation_mp(opts):
    global idx
    # global mp_out_file
      
    data = read_locations(opts.location_file)
    idx = make_spatial_index_memory(data)
    data_len = len(data)
    
    norm = opts.norm  # some normalisation depending on the total flow (or something)
    num_cpus = cpu_count() if opts.nprocs == 0 else opts.nprocs
    
    logging.info('num_cpus: {}'.format(num_cpus))
    logging.info('min_tij: {}'.format(opts.min_tij))
    
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    
    # create a counter to monitor progress of all tasks
    counter = manager.Value('counter', 0)
    
    partial_task = partial(radiation_task, data=data, norm=opts.norm, min_tij=opts.min_tij, lock=lock, counter=counter)
    pool = Pool(processes=num_cpus)
    
    # shuffling the order so that the progress is more or less fair
    r = range(len(data))    
    random.shuffle(r)
    pool.map(partial_task, r)
    
    # finish the pool
    pool.close()
    pool.join()
    
    logging.info('done')
    return
    
def radiation(opts):
    data = read_locations(opts.location_file)

    idx = make_spatial_index(opts.rtree_file, data)
    data_len = len(data)
    
    norm = opts.norm  # some normalisation depending on the total flow (or something)
    
    out_file = open('tij_ED.csv', 'w')
    
    for i, data_i in enumerate(data):
        pos_i, pop_i = data_i.pos, data_i.pop
        logging.info('working on {} / {}'.format(i, data_len))
        hits_i = [h for h in idx.nearest((pos_i[0], pos_i[1], pos_i[0], pos_i[1]), data_len, objects=False)]

        for j, data_j in enumerate(data):
            if i == j: continue            
            pos_j, pop_j = data_j.pos, data_j.pop
            # d_ij = euclidean_distance(pos_i, pos_j)

            # logging.info('working on {} {} / {}'.format(i, j, data_len))
            # logging.info('i,j {} -> {}'.format(pos_i, pos_j))

            # distance ordered iterator to objects nearest to location 'i'- only need to consider objects closer than d_ij
            s_ij = 0
            for hit in hits_i:
                if hit == j:
                    break
                else:
                    s_ij += data[hit].pop

            # the spatial index returns the locations i so exclude it from count
            # s_ij = s_ij - (pop_i + pop_j)
            s_ij = s_ij - (pop_i)
            mob_i_j = (1.0 * pop_i * norm)
            
            # logging.info('s_ij: {} {} -> {}'.format(pop_i, pop_j, s_ij))
            t_ij = (1.0 * mob_i_j * pop_i * pop_j) / ((pop_i + s_ij) * (pop_i + pop_j + s_ij))            
            if t_ij > 0.5:
              out_file.write("%d\t%d\t%d\t%d\t%f\t%f\n" % (i, j, pop_i, pop_j, s_ij, t_ij))
              # logging.info("i->j: %d -> %d : %f" % (i, j, t_ij))
        
        # just for testing
        #if i >= 10: break
           
    out_file.close()    
    
    return

def main():
    parser = OptionParser()
    parser.add_option("--rtree_file", dest="rtree_file", help="rtree file", type='str', default="location_rtree")
    parser.add_option("--location_file", dest="location_file", help="location file", type='str', default="radiation_ED.csv")
    parser.add_option("--outfile", dest="outfile", help="output file", default="tij_ED.csv", type='str')
    parser.add_option("--norm", dest="norm", help="normalisation factor", type='float', default=1.0)    
    parser.add_option("--mp", dest="mp", help="use multiprocessing", type='int', default=1)    
    parser.add_option("--nprocs", dest="nprocs", help="number of processes", type='int', default=0)    
    parser.add_option("--min_tij", dest="min_tij", help="min value for t_ij", type='float', default=0.000002)
    (opts, args) = parser.parse_args()

    if opts.mp == 1:
        radiation_mp(opts)
    else:
        radiation(opts)
    return

if __name__ == '__main__':
    main()
