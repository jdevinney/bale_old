## triangle
### Definition
Find the number of triangles in a given simple unweighted graph. 

A triangle is a set of three vertices {u,w,v} where edges {u,w}, {w,v} and {u,v} are in the graph.

### Algorithm
This uses matrix algebra approach to counting triangles in a graph.
The simple graph is presented as a lower triangular {0,1}-matrix.

We preform the matrix computation that counts the number of nonzeros in the matrix (LL^t .& L).
### Discussion
This is here to shadow the algorithms in bale_classic.

### References
See the book, "Graph Algorithms in the Language of Linear Algebra",
edited by Gilbert, and Kepner for more details on our approach to this problem.

