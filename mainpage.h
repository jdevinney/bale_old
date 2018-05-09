/*!
 * \defgroup libgetputgrp Libgetput library (libgetput)
 * This is a simple parallel macro library that exposes very few
 * functions to the users. The most important are the 'get' and 'put' functions.
 *
 * \defgroup exstackgrp exstack library (exstack)
 * Our first attempts at buffered communications librarys. Contains both synchronous
 * (exstack) and asynchronous (exstack2) modes.
 *  
 * \defgroup spmatgrp Sparse Matrix library (spmat)
 * This is a simple parallel sparse matrix library built on top of libgetput. It
 * implements a distributed compressed sparse row data structure.

 */

/*! \mainpage 
 *
 ******************************************************************************************
 
 The bale archive is a summary of work that has been done on programming models
  that help the programmers trade latency for bandwidth.
  
  The bale archive contains some library packages and an applications package.

  Libraries
  - \ref libgetputgrp simple parallel programming language that can be built on UPC or SHMEM. Supports
  basic PGAS functions (like get and put) as well as some atomics.
  - \ref exstackgrp Our first attempt at a buffered communications library was called exstack. exstack required bulk-synchronous programming. Our second attempt was exstack2, which allowed for more asynchronous programming. Both models are included in this library.
  - \b convey: The successor to the exstack libraries. Conveyors should offer better performance and more flexibility.
  - \ref spmatgrp A parallel sparse matrix library

  Apps
  - \subpage histogram_page
  - \subpage indexgather_page
  - \subpage toposort_page 
  
  ******************************************************************************************
  \page histogram_page Histogram
  
  The histogram application is a simple parallel histogram. Each
  processor has a list of updates which are actually just indicies 0 -
  M-1. The job of histogram is to count the number of occurances of
  each index across all PEs. We can accomplish this task by creating a
  distributed array T with length M and initializing it to zero. Then
  each PE can run through its list and increment the appropriate count
  in the array T. As described, this technique could lead to a race
  condition where two or more processors simultaneaously update the
  same entry in T. This can be by forcing all updates to be done with
  atomic increments.

  As described, histogram is extremely latency tolerant. That means we
  could buffer up our remote updates and send whole buffers instead of
  one update at a time. This is exactly what exstacks and conveyors
  help us do (so we don't have to go to the trouble of creating
  buffers and figuring out the right addresses every time we want to
  do something like this).

  ******************************************************************************************
  \page indexgather_page Indexgather
  
  If histogram is like a "put", then indexgather is like a "get". Each
  processor has a list of read indices into a distributed table T. The
  would like to read those indices from T and write the results
  somewhere locally. We can do this quite simply with libgetput's \ref
  lgp_get_int64() function. As described, the indexgather problem is
  also very latency tolerant and a fine candidate for buffered
  communications.

  ******************************************************************************************
  \page toposort_page Toposort
  
  The toposort algorithm is more complicated than histogram and
  indexgather. There will be a paper the more fully describes toposort. For now...
  
  Consider an upper-triangular matrix T with nonzero diagonal. Now
  randomly permute the rows and columns of T to get a matrix T'. The
  goal of toposort is to create row and column permutations (not
  necessarily the original permutations) such that when applied to T',
  forms an upper triangular matrix with nonzero diagonal.

 */

