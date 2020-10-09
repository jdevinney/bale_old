## transpose_matrix

### Definition
Compute the transpose of a given sparse matrix.

### Algorithm
We start by computing columns counts. These become row counts in the transpose. 
With these we can allocate the memory and set the row offsets for the transpose.
Then we go through the `nonzero[ ]` and `value[ ]` arrays one row at a time.
We write the given row and value to the location in the transpose matrix
given by the of the nonzero (column number) of the original matrix

### Discussion
This apps is a timer wrapper for the routine in the sparse matrix library.

This is C_bale to shadow the app in bale_classic.
It is an app in bale_classic because the communication pattern is interesting.

### References
