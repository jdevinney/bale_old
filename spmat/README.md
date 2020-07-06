# spmat (Sparse Matrix library)

The main data structure in this library is a distributed sparse matrix (implemented as Compressed Sparse Row (CSR)). The rows are distributed to PEs in a round-robin fashion and all nonzeros of any given row have affinity to a single PE. Besides their uses in physical sciences, sparse matrices are also useful in representing the adjacency matrix of a graph. This library is mostly a collection of functions that act on sparse matrices or graphs. There are also a few functions that act on or create permutations of the numbers {0,..., n-1}. 


Several functions in this library are instructive enough that they are implemented in a variety of ways (AGI, exstack, exstack2, and conveyors). Those functions are transpose_matrix, permute_matrix, and randpermp (create a random permutation of {0,...,n-1} in parallel).

## Matrix/Graph generation
* Flat uniform (uses Erdos-Renyi model)
* Geometric Random Graph (see [wikipedia/RandomGeometricGraph](https://en.wikipedia.org/wiki/Random_geometric_graph))
* Kronecker Product Graphs
* I/O using Matrix Market format