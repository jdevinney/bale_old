## toposort
### Definition
We are given a matrix that is a random row and column permutation 
of an upper triangular matrix (with ones on the diagonal).
Such a matrix has been called a morally upper triangular matrix.
This algorithm finds a row and column permutation that, when applied,
returns it to an upper triangular form.

### Algorithms
First we generate the problem by reading in or generating 
an upper triangular matrix with ones on the diagonal 
and then applying random row and column permutations.
This is the input to the toposort algorithm.
The output of toposort is a row and a column permutation that, 
if applied, would result in an upper triangular matrix.

We find these permutations by finding pivot and removing 
elements from the matrix.  A pivot element is the (row,col) pair 
for the nonzero in a row that has only one nonzero.
We put the row and col indices in the next available entries
in the row and column permutations and "delete" the row and column
from the matrix. We then repeat the process.
 
#### enqueuing pivots
In the process of "deleting" a row and column we will make rows that
only have a single non-zero.

   We set the row and column permutations,  rperm and cperm, one pivot at a time.

   N = number of rows
   for( pos=N-1; pos > 0; pos-- ) {
     pick a row, r, with a single nonzero, c.
     say (r,c) is the pivot and set rperm[pos] = r and cprem[pos] = c
     Note: a pivot always exists because the matrix is morally upper tri.

     cross out that row r and col c 
   }

   Meaning of cross out:
   Rather than changing the matrix by deleting rows and column and then searching the 
   new matrix for the next pivot.  We do the obvious thing of keeping row counts, where
   rowcnt[i] is the number of non-zeros in row i and we use a really cool trick 
   of keeping the sum of the live column indices for the non-zeros in each row.
   That is, rowsum[i] is the sum of the column indices, not the sum of the non-zero elements,
   for the non-zeros in row i.  To "delete a column" one decrements the rowcnt by one and 
   the rowsum by the corrsponding column index. 
   The cool trick is that, when the rowcnt gets to one, the rowsum is the column that is left.
#### loop to find pivots


### Discussion

### References
