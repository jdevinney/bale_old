# apps

Each of the applications in bale are implemented in multiple ways to showcase the pros and cons of each. In general the models used are:

1. agi      : standard PGAS model that uses puts, gets, and atomics
2. exstack  : a bulk synchronous buffering model (see [exstack](../exstack/README.md) library)
3. exstack2 : a asynchronous variant of exstack (see [exstack](../exstack/README.md) library)
4. conveyors: a more mature and sophisticated asynchronous model that is independent of exstack and exstack2.  (see [convey](../convey/README.md) library)
5. In some applications we have included other variants in the "alternatives" directory.   

The apps are:

- [histogram](histo_src/README.md)
- [indexgather](ig_src/README.md)
- [toposort](topo_src/README.md) 
- [transpose_matrix](transpose_matrix_src/README.md) 
- [triangle counting](triangle_src/README.md)
- [randperm](randperm_src/README.md)
- [permute_matrix](permute_matrix_src/README.md)
- [write_sparse_matrix](write_sparse_matrix_src/README.md)

