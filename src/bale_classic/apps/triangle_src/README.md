# Triangle counting

This uses matrix algebra approach to counting triangles in a graph.
See the book, "Graph Algorithms in the Language of Linear Algebra",
edited by Gilbert, and Kepner for more details on our approach to this problem.

We have implemented two algorithms to count triangles in a graph. The
first computes (L & L * U) and the second computes (L & U * L), where
'&' means element-wise AND and '*' is ordinary matrix
multiplication. This application is interesting because it is quite
amenable to aggregation and there are ways to reduce communication if
one wanted to speed up the calculation.
