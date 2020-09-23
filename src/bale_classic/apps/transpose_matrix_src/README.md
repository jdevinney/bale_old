# transpose_matrix

## Description

In this application we transpose a distributed sparse matrix. See the [spmat](../../spmat/README.md) library for this data structure. 

## Discussion

The algorithm we use to transpose a matrix has two parts. In the first part, PEs use a histogram pattern to calculate the number of entries in each column of the matrix. This allows us to allocate the exact space needed to store the transpose matrix. In the second phase, we
again use the histogram pattern, but this time we send the nonzeros of the transpose matrix to the correct PE. For example, if *A[i,j]* lives on PE k in the original matrix, PE *k* sends a *(i,j)* to PE *m* to create *AT[j,i]* (where *AT* is the transpose of *A*). Both phases of the transpose function are well-suited for aggregation.

An interesting difference between the AGP and aggregated versions comes in the second phase. In the AGP version, PEs are placing nonzeros in AT directly into the distributed sparse matrix data structure via remote writes. To do that in parallel, PEs must be able to atomically reserve a spot in the "nonzero" array for their writes. Since we don't know in what order these nonzeros will arrive. the PEs are contending for these writes. In the aggregated versions, PEs are being sent nonzeros (via an aggregation library) and process them in serial. So there is no need for atomic operations. This phenomenon occurs in other bale apps (histogram for example) and is fairly common when going from AGP style paralell code to aggregated code.