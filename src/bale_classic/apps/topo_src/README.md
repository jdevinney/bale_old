# toposort

## Description

The input to toposort is a sparse matrix A, where A is "morally
unit-upper-triangular". That is, we start we an upper-triangular
matrix *T* with ones on the diagonal (that's the "unit" part of the
name). We don't care about the values of the non-zeros in the matrix,
only their position. Next, we randomly (and separately) permute the
rows and columns of *T* to get a matrix *M*.  The matrix *M* has been
called a morally (unit) upper-triangular matrix.

Given a morally unit-upper-triangular matrix, *M*, the goal of
toposort is to create row and column permutations such that when these
permutations are applied to *M*, the result is an unit-upper
triangular matrix. The answer need not be unique but since *M* is a
row and column permutation of *T*, we know T is a solution.

If you observe the set of the column labels of nonzeros
in a particular row; a column permutation might change the labels of
the elements in the set, but doesn't change the cardinality of the
set. Likewise for columns. Hence, there must be a row in *M* with a
single non-zero.  If we remove that row and the column it intersects,
we are left with smaller, morally unit upper triangular matrix. This is the
motivation behind a simple algorithm (and the outline of an induction
proof of its correctness).

The outline of the algorithm is as follows:

* For all rows r with a single non-zero in a column c, put the pair (r,c) onto a queue.
* While the queue is not empty: 
  * pop an (r,c) pair from the queue
  * claim a new position as the last row and column of the permutations being created 
  * remove all the non-zeros in column c
  * if any row now has a single non-zero, enqueue that (r,c) pair

Rather than changing the matrix by deleting rows and column and then searching the 
new matrix for the next row.  We do the obvious thing of keeping an array of row counts.
*rowcnt[i]* is the number of undeleted non-zeros in *row i* and
we use a cool trick to find the last surviving column of a row.
We initialize an array, *rowsum[i]*, to be the sum of the column indices in *row i*.
When we "delete" a column we decrement *rowcnt[i]* by 1 and *rowsum[i]* by that column index.
Hence, when the *rowcnt[i]* = 1, *rowsum[i]* contains the column that is left.

In parallel there are three race conditions or synchronization issues to address.

The first is reading and writing the queue of rows to be processed.
One way to handle it is to introduce the notion of a levels.
Within a level all threads process the all the rows on their queues 
and by doing so create new degree one rows. These rows are placed on the 
appropriate queues for the next level. There is a barrier between levels.

Threads race to pick their position in *rperm* and *cperm*. 
One could handle this race for the pivots with a fetch_and_add,
instead we use parallel prefix to claim enough room for the pivots 
in the current level on each thread then assign them in order per thread.

Threads race to update the *rowcnt* and *rowsum* arrays. 
We handle this with levels and atomic memory operations.

## Discussion

The toposort application is a step up in complexity from [histogram](../histo_src/README.md) and
[indexgather](../ig_src/README.md).  Depending on its input, it typically enjoys a
significant amount of parallelism, but it is not completely order and
latency tolerant. One of the most intereting things about toposort is
how well it works with aggregation techniques (like exstack or
conveyors). It is not obvious from the description of the algorithm
that it should be so well suited for aggregation. Over several rounds
of brainstorming, we have found new ways to think about implementing
toposort that make it even more amenable to aggregation.  In this way,
toposort is a poster-child for aggregation. 


