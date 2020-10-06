# Sparse Matrix IO

## Description

The goal of this app to write to disk and then read from disk a distributed sparse matrix dataset.

## Discussion

It sounds like a very simple task. However, we throw in a twist. bale sparse matrices are stored so that the rows are distributed to PEs in CYCLIC order (that means rows i, i + THREADS, i + 2*THREADS... are local to PE i; this is a contrast to BLOCK layout where PEs get continguous sets of rows). There are pros and cons to both CYCLIC and BLOCK layouts, One of the nice features of CYCLIC layout is that nonzeros tend to be better load balanced to PEs. 

We would like to be able to write a sparse matrix in a way such that we can read it back in on any number of PEs easily. This is one of the pros to BLOCK layout. It is very easy to write a matrix in BLOCK layout just by having each PE dump their data to disk. 

We would also like the data to be laid out in row-major order on disk. That would make it more convenient to explore the data on disk if necessary. Both of these goals are automatically and easily achieved using BLOCK layout. However, as we mentioned, CYCLIC layout has its perks too. We would like to be able to read and write matrices from/to CYCLIC layout with any numbers of PEs.

#### Sparse Matrix Dataset

### Writing

One way to do this task would be to shuffle rows from CYCLIC layout to BLOCK layout and then just dump the data. That motivates how the bale implementations of write sparse matrix work. As described this sounds just like the [permute_matrix](../permute_matrix_src/README.md) app. However, what if we do not want to allocate the space to store an entire second copy of the matrix just to write the matrix to disk? We need to send the rows in batches to the PE that would own them if the matrix were in BLOCK layout and write small buffers of data. This is where an additional challenge of this application comes in. We need to send just the right data so that each PE is able to fill a write buffer with the data it is supposed to write. This makes this code especially challenging for aggregation, especially asynchronous aggregation. In fact, we only have an AGP and exstack version of write_sparse_matrix currently. We don't think an exstack2 version would be reasonable. We hope to soon write an efficient conveyor version.

### Reading

Our read_sparse_matrix routines (again, we currently only have AGP and exstack versions) have an extra parameter: the number of readers. We have seen evidence that as the number of cores on multi-core processors grows, I/O could suffer when the number of cores concurrently reading or writing gets too large. For this reason, it seems logical that we should limit the number of PEs taking part in the reading or writing. We have elected to implement this twist in the read_sparse_matrix apps. 

### Why is this in bale?

Another idea for the write is to always just dump the data in CYCLIC layout. Then when reading back in, the PEs need to figure out which rows they read and where those rows go. This can be tricky since the rows are laid out in a way that depends on the number of writing PEs and the communication pattern during the read depends on the number of PEs used to write and how many are reading. 