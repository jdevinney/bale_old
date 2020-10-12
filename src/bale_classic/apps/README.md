# apps

The applications in bale are meant to showcase:

- the challenges in writing interesting, distributed, parallel applications that are efficient and high performing at scale
- the challenges of getting these codes to use aggregated communications. Including:
  - getting the code right the first time
  - rapidly experimenting at scale during algorithm development
  - reading and understanding the code
  - implementing algorithms with different tolerances for latency
- our quest for the From the Book (FTB) implementation of each app

### Implementations

Each of the applications in bale are implemented in multiple ways to showcase the pros and cons of each. In general the models used are: AGP : standard PGAS model that uses Atomics, Gets, and Puts (AGP), [exstack](../exstack/README.md), [exstack2](../exstack/README.md), [convey](../convey/README.md). In some applications we have included other variants in the "alternatives" directory.

### List of apps

- [histogram](histo_src/README.md)
- [indexgather](ig_src/README.md)
- [transpose_matrix](transpose_matrix_src/README.md) 
- [randperm](randperm_src/README.md)
- [permute_matrix](permute_matrix_src/README.md)
- [triangle counting](triangle_src/README.md)
- [toposort](topo_src/README.md) 
- [sparse_matrix_io](sparse_matrix_io_src/README.md)
- [Single Source Shortest Path (sssp)](sssp_src/README.md)

