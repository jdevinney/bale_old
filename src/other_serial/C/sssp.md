## sssp Single Source Shortest Path
### Definition
We are given the adjacency matrix for a graph with non-negative edge weights *c(v,w)*
and an given source vertex, *v_0*.
We wish to find the lightest weighted path from *v_0* to all other vertices,
where the weight of a path is the sum of the weights of the edges in the path.
Note, if the graph is undirected we are given the full (symmetric) adjacency matrix.

### Algorithms
We consider three algorithms: Dijsktra's, Delta-Stepping, and Bellman-Ford.
Dijsktra's algorithm is not in bale_classic because it is a serial algorithm.
Delta-Stepping and Bellman-Ford are here as shadows of the parallel versions.

//TODO I don't want to repeat all the stuff in bale_classic README
### Discussion

// TODO move to README

### References
"Delta-stepping: a parallelizable shortest path algorithm" by U. Meyer and P. Sanders.
