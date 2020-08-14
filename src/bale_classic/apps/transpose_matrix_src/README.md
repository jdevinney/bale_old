# transpose_matrix

In this application we transpose a distributed sparse matrix. 
See the [spmat](../../spmat/README.md) library for this data structure. 
The algorithm we use to transpose a matrix has two parts. In the first
part, PEs use a histogram pattern to calculate the number of entries
in each column of the matrix. This allows us to allocate the exact
space needed to store the transpose matrix. In the second phase, we
again use the histogram pattern, but this time we send the nonzeros of
the transpose matrix to the correct PE. 
For example, if *A[i,j]* lives on PE k in the original matrix, 
 PE *k* sends a *(i,j)* to PE *m* to create *AT[j,i]* (where *AT* is the transpose of *A*). 
Both phases of the transpose function are well-suited for aggregation.
