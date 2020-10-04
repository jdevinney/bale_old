/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/


/*! 
\mainpage Serial C version of Bale
\section intro Introduction 

Why bother with is?

\section contents Contents 
This directory contains a serial version of the bale apps written in C.
There is a serial version of the ``spmat`` library and a few of the utility
functions that would be in the ``libgetput`` library.
These are combined in the files ``spmat_utils.h`` and ``spmat_utils.c``.
We also use the same argument parsing routines as we do in bale_classic
(except we don't use the options needed to handle parallel threads).

The apps Serial C Bale are:
  - **histo**  tests random writes to memory: ``histo.c`` 
  - <b>ig</b>  tests random loads from memory: ``ig.c``
  - <b>randperm</b>  generate a random permutation: ``randperm.c``
  - <b>permute_matrix</b> permute rows and columns of a sparse matrix: ``permute_matrix.c``
  - <b>transpose_matrix</b> transpose a sparse matrix: ``transpose_matrix.c``
  - <b>triangle</b> count triangles in a simple graph: ``triangle.c``
  - <b>toposort</b> a bfs algorithm verify that a matrix is "morally" upper triangular: ``toposort.c``
  - <b>sssp solve</b> Single Source Shortest Path problem on a weighted graph: ``sssp.c``
  - <b>unionfind</b> application of the disjoint union data structure: ``unionfind.c``

Other:
  - **spmat_utils** sparse matrix library and utility functions: ``spmat_utils.c`` and ``spmat_utils.h``
  - the use of argp to parse command line ``std_options.h`` and ``std_options.c``
  - the use of pytest to run some unit tests [pytest](pytest.md)
*/

/*
\defgroup bale_c C_Bale
 * This is a version of the bale apps has been written as simple serial C code.
 * There is no need for a much of the infrastucture code that is in the 
 * "bale classic" parallel version.
 * We wrote this version in the hope that it might start discussions
 * among those who are unfamiliar with UPC or SHMEM.
 */

/* \page pytest
Description of how to run pytest

#include "pytest.md"

*/
