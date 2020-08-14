# libgetput (A basic parallel library)

This library is meant to act as a wrapper for some of the
functionality that is common between UPC and SHMEM. It allows most of
bale to be easily compiled against UPC or SHMEM. libgetput, as its
name implies, provides basic remote get and put functions:
lgp_get_int64, lgp_put_int64 and lgp_getmem and lgp_putmem (which are
similar to SHMEM functions for single word gets and puts or more
general gets and puts of memory). One key distinction is that the
indexing in libgetput is UPC style indexing. That is we consider the
distributed array to be indexed in round-robin fashion: the first
element has affinity to PE 0, the next element has affinity to PE 1
and so on.

libgetput also supplies a variety of atomic functions, both fetching and non-fetching. Finally, libgetput provides some fundamental collectives: value-based reductions and barriers.