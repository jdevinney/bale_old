# Triangle counting

This uses matrix algebra approach to counting triangles in a graph.

We have implemented two algorithms to count triangles in a graph. The
first computes (L & L * U) and the second computes (L & U * L), where
'&' means element-wise AND and '*' is ordinary matrix multiplication. This application is intereting because it is quite amenable to aggregation and there are ways to reduce communication if one wanted to speed up the calculation.

We have several ways to create input for the triangle counting application. The first is our ordinary random graphs from the spmat library. The second is more tailored to this application: Kronecker product graphs. These graphs can be created with a known number of triangles which helps in the validation of the code. See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al. for more information. Finally, we can read a matrix in matrix-market format.