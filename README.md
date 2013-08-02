spatialradiation
================

The radiation model of spatial attraction was proposed by Simini et al [1]. 

The basic result is the average flux $$$T\_{ij}$$$ from location $$$i$$$ to location $$$j$$$:

$$ \left\< T\_{ij} \right\> = T_i \frac{m_i n_j}{(m_i + s\_{ij}) (m_i + n_j + s\_{ij})}$$

where $$$T_i = \sum\_{j \neq i} T\_{ij}$$$ is the total number of commuters that start their journey at location $$$i$$$, and $$$m_i$$$, $$$n_j$$$ are the populations at locations $$$i$$$ and $$$j$$$.

The term $$$s\_{ij}$$$ is the sum of populations in a circle centered on location $$$i$$$, and with radius $$$r\_{ij}$$$ the distance between $$$i$$$ and $$$j$$$.

Given the locations and populations, the flux is fairly straightforward to compute. Here, we use an rtree index to improve the performance of the computation of $$$s\_{ij}$$$. The rtree index can efficiently return the locations within a fixed distance without computing any distances.

Since the flux for each pair of locations can be independently computed, we use multiprocessing tasks to parallelise the computation.

The input should be a list of comma separated values

    location_id, x, y, location_str, population
    0,-7.045541,54.379834,34066,219.000000
    1,-7.144696,54.305961,34067,648.000000

Each multiprocessing job $$$k$$$ writes to its own file "t_ij_k.dat", which can be easily joined together at the end:

    cat t_ij_*.dat | sort -n -k 1 -k 2 > all_t_ij.dat



Improvements
======
There are many possible improvements, and ideas are welcome. For the moment the use of the rtree means we are restricted to a 'flat earth' approximation of distance.





[1]: http://www.nature.com/nature/journal/v484/n7392/fig_tab/nature10856_F1.html "A universal model for mobility and migration patterns"