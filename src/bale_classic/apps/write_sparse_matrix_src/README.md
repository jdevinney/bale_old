# Write Sparse Matrix

## Description

The goal of this app is pretty self-evident. Write a distributed sparse matrix to disk. 

## Discussion

It sounds like a very simple task. However, we throw in a twist. bale sparse matrices are stored so that the rows are distributed to PEs in CYCLIC order (that means rows i, i + THREADS, i + 2*THREADS... are local to PE i; this is a contrast to BLOCK layout where PEs get continguous sets of rows). There are pros and cons to both CYCLIC and BLOCK layouts, One of the nice features of CYCLIC layout is that nonzeros tend to be better load balanced to PEs. 

We would like to be able to write a sparse matrix in a way such that we can read it back in on any number of PEs easily. This is one of the pros to BLOCK layout. It is very easy to write a matrix in BLOCK layout just by having each PE dump their data to disk. This data is also very easy to read back into BLOCK layout, using any number of PEs.

We would also like the data to be laid out in row-major order on disk. That would make it more convenient to explore the data on disk if necessary. Both of these goals are automatically and easily achieved using BLOCK layout. However, as we mentioned, CYCLIC layout has its perks too. We would like to be able to read and write matrices from/to CYCLIC layout with any numbers of PEs.

### Writing

One way to do this task would be to shuffle rows from CYCLIC layout to BLOCK layout and then just dump the data. That motivates how the bale implementations of write sparse matrix work. As described this sounds just like the [permute_matrix](../permute_matrix_src/README.md) app. However, we what if we do not want to allocate the space to store an entire second copy of the matrix just to write the matrix to disk? We need to send the rows in batches to the PE that would own them if the matrix were in BLOCK layout and write small buffers of data. This is where the challenge of this application comes in. We need to send just the right data so that each PE is able to fill a write buffer with the data it is supposed to write. This makes this code especially challenging for aggregation, especially asynchronous aggregation. In fact, we only have an AGP and exstack version of write_sparse_matrix currently. We don't think an exstack2 version would be reasonable. We hope to soon write an efficient conveyor version.

### Reading

Though we do not have a read_sparse_matrix in bale yet, such an app would be very similar to the write_sparse_matrix app. PEs would read some amount of the on disk data and farm it out to the PEs who own that data in CYCLIC layout. 

### Other thoughts

Another idea for the write is to always just dump the data in CYCLIC layout. Then when reading back in, the PEs need to figure out which rows they read and where those rows go. This can be tricky since the rows are laid out in a way that depends on the number of writing PEs and the communication pattern during the read depends on the number of PEs used to write and how many are reading. 




Demo program that runs the variants of write_sparse_matrix kernel. It first generates 
a random matrix in FLAT mode and then it writes this matrix to disk
in a directory called 'write_sparse_test'.

We define a sparse matrix dataset to be the following:
 - It lives in a directory of its own
 - It contains one ASCII file called 'metadata' which contains the number of rows, columns and nonzeros in the matrix.
 - It contains N binary 'rowcnt' files. These files contain the number of nonzeros in each row for rows 0..A->numrows
 - It contains N binary 'nonzero' files. These files contain the nonzeros in each row and are ordered by row.

This application is interesting because, as implemented, it requires
us to shuffle the rows of the matrix from cyclic to block. That is,
the nonzeros for row i is stored on PE (i % THREADS) in the
sparsemat_t data structure.  However, we wish PE 0 to write out the
first block of (approx) A->numrows/THREADS rows. One way to do this
would be to call permute_sparse_matrix to get a copy of the matrix
whose rows are distributed in the block layout. However, that would
require 2x the storage space of the matrix. We don't want the
write_sparse_matrix routine to have this requirement. That is where
things get interesting. We have a fixed buffer for writing on each PE
and each PE collects or is sent nonzero data to write in their current
buffer. In AGP (where the PEs just get the data) or synchronous
exstack, this is easy. We don't have a nice way of doing this with
exstack2 or asynchronous conveyors yet. The reason this is a challenge
for asynchronous methods is that PEs can get into a deadlock waiting
for the records they need to complete a write buffer.


