# spmat (Sparse Matrix library)

The main data structure in this library is a distributed sparse matrix
(implemented as Compressed Sparse Row (CSR)). The rows are distributed
to PEs in a round-robin fashion and all nonzeros of any given row have
affinity to a single PE. Besides their uses in physical sciences,
sparse matrices are also useful in representing the adjacency matrix
of a graph. This library is mostly a collection of functions that act
on sparse matrices or graphs. There are also a few functions that act
on or create permutations of the numbers {0,..., n-1}.

Several functions in this library are instructive enough that they are
implemented in a variety of ways (AGI, exstack, exstack2, and
conveyors). Those functions are transpose_matrix, permute_matrix, and
randpermp (create a random permutation of {0,...,n-1} in parallel).

## Matrix/Graph generation
* Flat uniform (uses Erdos-Renyi model)
* Geometric Random Graph (see [wikipedia/RandomGeometricGraph](https://en.wikipedia.org/wiki/Random_geometric_graph))
* Kronecker Product Graphs
* I/O using Matrix Market format

### Erdos-Renyi
The Erdos-Renyi random graph model with parameter p for
n vertices flips a weighted coin (heads with probability p) for every
potential edge in the graph. The edge is inserted into the graph if
the coin flip results in a heads, and left out otherwise. In
matrix-speak, we are generating the lower-half of the adjacency matrix
of a graph and inserting a 1 for each heads and a 0 otherwise. The
pseudo-code :

    for i = 0...n
      for j = 0...i
        if(random() < p)
	   A[i][j] = 1
	else
	   A[i][j] = 0

Generating a large graph according to this psuedo-code is inefficient
as it requires O(n^2) random numbers to be generated. A better way is
found in the paper, "Efficient Generation of Large Random Networks" by
Bategeli and Brandes. This uses the fact that the geometric
distribution models the number of tails between two heads. We have
implemented both the naive and the more efficient algorithms in
bale. Also note that one can create a directed graph with this model
but flipping a coin for each potential directed edge between two
vertices.

### Geometric Random Graphs

This random graph model is rather simple. Similar to the Erdos Renyi
model, it has a single parameter r (between 0 and 1). For each vertex,
a point is randomly placed in the unit square. An edge is placed
between vertices i and j if the points corresponding to these vertices
are within distance r of each other. To generate a graph under this
model we break the unit square into square sectors (usually of length
and width r). This reduces the number of interpoint distances we need
to calculate since edges can only exist between points in the same
sector or in neighboring sectors (including diagonally neighboring).

These graphs present an interesting alternate to the Erdos-Renyi
random graph. Their generation in parallel is also interesting in its
own right. The assignment of points to PEs during the edge generation
to reduce the amount of communication is an interesting discussion, as
is the assignment of rows of the adjacency matrix once the edges are
determined. These need not be the same and in fact it would be
difficult to make them so given the sparse matrix data structure's
requirment that the rows of a matrix must be evenly distributed to the
PEs (with any remainder rows going to the lowest indexed PEs). Our
current implementation assigns points to sectors starting with the top
left sector and working left-to-right and wrapping to the next row
down at the end of a row of sectors. Point numbers increase with
sector number. We assign sectors in round-robin order to PEs for edge
generation. Then, to create the matrix, points are distributed in
block order (the first k points go to PE 0, the next k go to PE 1,
etc). This allows us to evenly distribute the points to PEs while
maintaining some of the locality inherent in these graphs.

While we know there exist communication-free algorithms for generating
geometric random graphs, we decided to implement an algorithm that
requires communication. Why? Bale is primarily about improving the
lives of distributed parallel programmers. The communication free
algorithm is not representative of an algorithm a researcher would
think of on his/her first attempt. We want to see how different
parallel programming models behave under the kinds of algorithms that
people write as they evolve their algorithms.

### Kronecker Product Graphs

We chose to implement Kronecker product graphs in bale to test out our Triangle counting implementations.
For more details see. The parallel generation of these graphs is not particularly challenging or interesting. See
"Design, Generation, and Validation of Extreme Scale Power-Law Graphs"
by Kepner et. al. for more details.

### Matrix Market I/O

The spmat library has the ability to read
matrices in Matrix Market format. This function reads the matrices in serial
using one PE and then distributes the resulting matrix to all PEs. For
this reason it is not meant to scale to large matrices, but it is
useful for debugging or sanity checking on known examples.