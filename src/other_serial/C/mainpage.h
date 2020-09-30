/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*!
 * \defgroup Bale_C
 * This is a version of the bale apps written in simple serial C code.
 * There is no need for a infrastucture code in the parallel version 
 * of Bale. 
 * We wrote this version in the hope that it might start discussions
 * among those who are unfamiliar with UPC or SHMEM.
 *
 * \defgroup spmatgrp Sparse Matrix library (spmat)
 * We also have a parts of the sparse matrix library, again,
 * the serial version is much less complicated than the parallel version.
 */

/*! \mainpage
 *

 \section Bale apps in C.
   Why bother with is?

 \section Whats here:
  This contains a serial version of each of the apps in bale.
  The parallel programming library <tt>libgetput</tt> is not needed.
  The parts of the <tt> spmat </tt> library that are needed just 
  live a single file <tt> spmat_utils.c</tt>.

  The apps:
  - \subpage histo_page histogram a large number of items into a large table
  - \subpage ig_page simple random loads from memory
  - \subpage randperm_page generate a random sparse matrix
  - \subpage permute_matrix_page permute a sparse matrix
  - \subpage transpose_matrix_page transpose a sparse matrix
  - \subpage triangle_page count triangles in a graph (given by a lower triangular matrix)
  - \subpage toposort_page a bfs algorithm to re-order the rows and columns of  "morally" upper triangular matrix
  - \subpage unionfind_page application of the union-find data structure
  - \subpage spmat_utils_page stuff

  Also:
  - <tt>runall.sh</tt> a bash script to run all the apps with trivial parameters.

*/

