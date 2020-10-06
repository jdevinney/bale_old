## transpose_matrix
### Definition
Given a sparse matrix, produce the sparsemat_t data structure for the transpose of the matrix.
This apps is a timer wrapper for the routine in the sparse matrix library.

### Algorithm
We start by computing columns counts. These become row counts in the transpose. With these
we can allocate the memory and set the row offsets in the transpose.
Then we go through the `nonzero[ ]` and `value[ ]` arrays one row at a time.
We write the given row and value to the appropriate location 
of the `nonzero[ ]` and `value[ ]` arrays of the transpose sparsemat_t.

### Discussion

### References
