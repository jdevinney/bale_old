## permute_matrix

### Definition
We apply given row and column permutations to a sparse matrix.

### Algorithm
We produce a new sparse matrix data structure by copying the nonzeros and value entries
to their new positions based on the row offsets of the permuted matrix.  As we are copying
the entries we replace the nonzeros (the column indices) with their new column indices 
given by the column permutation.

### Discussion
This app is really just a wrapper that calls the routine in the sparse matrix library. 

It is a app in bale_classic because the communication pattern is interesting.

### References

