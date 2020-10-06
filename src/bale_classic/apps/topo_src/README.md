# toposort

#### Glossary

*unit-upper-triangular matrix*: A matrix that is upper-triangular and where every diagonal entry is equal to 1.

*morally unit-upper-triangular*: A matrix M is morally unit-upper-triagular if there are square permutation matrices P and Q and a unit-upper-triangular matrix T, such that P M Q = T. 

## Definition

The input to toposort is a morally unit-upper-triangular sparse matrix M.  To obtain M we create a unit-upper-triangular matrix and apply random permutations to the rows and columns of that matrix. Note that we don't care about the values of the non-zeros in any matrix in toposort, only their position.  The goal of
toposort is to create row and column permutations such that when these permutations are applied to *M*, the result is an unit-upper triangular matrix. The answer need not be unique. A solution always exists since *M* is a row and column permutation of *T*.

#### Base Algorithm

If you observe the set of the column labels of nonzeros
in a particular row; a column permutation might change the labels of
the elements in the set, but doesn't change the cardinality of the
set. Likewise for columns. Hence, there must be a row in *M* with a
single non-zero.  If we remove that row and the column it intersects,
we are left with smaller, morally unit upper triangular matrix. This is the
motivation behind a simple algorithm (and the outline of an induction
proof of its correctness).

The outline of an algorithm is as follows:

Let M be an N by N morally unit-upper-triangular matrix. 

* For all rows r with a single non-zero in a column c, put the pair (r,c) onto a queue.
* pos = N-1
* While the queue is not empty: 
  * pop an (r,c) pair from the queue
  * rperm[r] = pos; cperm[c] = pos; pos--;
  * remove all the non-zeros in column c
  * if any row now has a single non-zero, enqueue that (r,c) pair

Rather than changing the matrix by deleting rows and column and then searching the 
new matrix for the next row.  We do the obvious thing of keeping an array of row counts.
*rowcnt[i]* is the number of undeleted non-zeros in *row i*. We use a cool trick to find the last surviving column of a row. We initialize an array, *rowsum[i]*, to be the sum of the column indices in *row i*.
When we "delete" a column we decrement *rowcnt[i]* by 1 and *rowsum[i]* by that column index.
Hence, when the *rowcnt[i]* = 1, *rowsum[i]* contains the column that is left. 

## Discussion

#### Parallel Considerations

In parallel there are three race conditions or synchronization issues to address.

1. The first is reading and writing the queue of rows to be processed.
   One way to handle it is to introduce the notion of a levels.
   Within a level all threads process the all the rows on their queues 
   and by doing so create new degree one rows. These rows are placed on the 
   appropriate queues for the next level. There is a barrier between levels.

2. Threads race to pick their position in *rperm* and *cperm*. 
   One could handle this race for the pivots with a fetch_and_inc for each new pivot. An improvement on this idea is to use one fetch_and_add to claim enough room for all the local pivots 
   in the current level on each thread then assign them in order per thread.

3. Threads race to update the *rowcnt* and *rowsum* arrays. 

#### Why is it in bale?	

The toposort application is a major step up in complexity from [histogram](../histo_src/README.md) and
[indexgather](../ig_src/README.md). Successfully implementing toposort (including all of the pre-computation of its input) requires 3 other bale apps (transpose_matrix, randperm, permute_matrix). In this way, toposort represents a crucible for any new parallel programming model.

Depending on its input, toposort typically enjoys a significant amount of parallelism, but it is not completely order and latency tolerant. One of the most interesting things about toposort is how well it works with aggregation techniques (like exstack or conveyors). However, it is not obvious from the description of the algorithm that it should be so well suited for aggregation. In this way, toposort is a poster-child for the wide applicability of aggregation. 

Over several rounds of brainstorming, we have found new ways to think about implementing toposort that make it even more amenable to aggregation. See the toposort_cooler.upc code in the alternates directory for instance. 

#### From the Book?

We are not really satisfied with any of the versions of toposort in bale_classic. A more modern language (like Chapel or Rust) helps make toposort more readable. But we still haven't seen the one FTB.



The toposort algorithm is more complicated than histogram and indexgather.
It typically enjoys a significant amount of parallelism, but it 
is not completely order and latency tolerant.

To prepare the input to the toposort algorithm, 
we start we an upper-triangular matrix <b>T</b> with no zeros on the diagonal. 
We don't care about the values of the non-zeros in the matrix, only their position.
Next we randomly permute the rows and columns of <b>T</b> to get a matrix <b>M</b>. 
The matrix <b>M</b> has been called a morally triangular matrix.

Given a morally triangular matrix, <b>M</b>, the goal of toposort is 
to create row and column permutations such that when these 
permutations are applied to <b>M</b>,
the result is an upper triangular matrix with no zeros on the diagonal.
Note, the answer need not be unique and since <b>M</b> is a row and column permutation
of <b>T</b>, there must be a solution.

We use a breadth first search algorithm based on the following observations:
    - The rows (and columns) of a sparse matrix partition the set of non-zeros in the matrix.
    - Row and column permutations preserve both partitions.

For example, if you create a set of the nonzeros in a particular row;
a column permutation might change the labels of the elements in the set,
but doesn't change the cardinality of the set. Likewise for columns.
Hence, there must be a row in <b>M</b> with a single non-zero.
If we remove that row and column, 
we are left with smaller, morally upper triangular matrix. This is the motivation behind a simple algorithm.

The outline of the algorithm is as follows:
\verbatim
For all rows with a single non-zero, put its non-zero onto a queue.
While the queue is not empty: 
   pop a non-zero from the queue (this represents row r and column c)
   claim its new position as the last row and column of the permutations being created 
   remove all the non-zeros in column c
   if any row now has a single non-zero, enqueue that non-zero
\endverbatim

Rather than changing the matrix by deleting rows and column and then searching the 
new matrix for the next row.  We do the obvious thing of keeping and array of row counts,
<b>rowcnt[i]</b> is the number of non-zeros in <b>row i</b> and
we use a cool trick to find the column of a row with <b>rowcnt[i]</b> equal 1.
We initialize an array, <b>rowsum[i]</b>, to be the sum of the column indices in <b> row i</b>.
When we "delete" a column we decrement <b>rowcnt[i]</b> and <b>rowsum[i]</b> by that column index.
Hence, when the <b>rowcnt[i]</b> gets down to one, the <b>rowsum[i]</b> is the column that is left.

In parallel there are three race conditions or synchronization issues to address..

The first is reading and writing the queue of rows to be processed.
One way to handle it is to introduce the notion of a levels.
Within a level all threads process the all the rows on their queues 
and by doing so create new degree one rows. These rows are placed on the 
appropriate queues for the next level. There is a barrier between levels.

Threads race to pick their position in <b>rperm</b> and <b>cperm</b>. 
One could handle this race for the pivots with a fetch_and_add,
instead we use parallel prefix to claim enough room for the pivots 
in the current level on each thread then assign them in order per thread.
   
Threads race to update the <b>rowcnt</b> and <b>rowsum</b> arrays. 
We handle this with levels and atomic memory operations.
