/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
*****************************************************************/ 
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

