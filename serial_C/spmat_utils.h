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

/*! \file spmat_utils.h
 * \brief The header file for spmat library.
 */ 
#ifndef spmat_utils_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <fcntl.h>
#include <assert.h>
#include <stddef.h>
#include <getopt.h>


/*! \struct sparsemat_t 
 * \brief A structure to hold a sparse matrix 
 * \ingroup spmatgrp
 We essentially use the standard Compressed Row Format
 Since we only care about the access patterns, we don't 
 to keep track of the values themselves for the nonzeros
 in the matrix only the position of the nonzeros.
 This cuts the volume to traffic (not the pattern) by some amount.
 This also saves a lot local sorting and combining that would
 be required to actually handle values (other than one).
 
 We store the nonzeros with affinity by rows.
 Ie, if row i has affinity to thread t, then all the nonzeros
 in row i have affinity to thread t.    
 */
typedef struct sparsemat_t{
  int64_t numrows;       //!< the total number of rows in the matrix
  int64_t numcols;       //!< the nonzeros have values between 0 and numcols
  int64_t nnz;           //!< total number of nonzeros in the matrix
  int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * nonzero;     //!< the global array of column indices for nonzeros
  double * value;       //!< the global array of nonzero values
}sparsemat_t;

/*! \brief struct to hold the state used to iterate across the row of a sparse matrix.
 * \ingroup spmatgrp
 */
typedef struct next_nz {
  sparsemat_t * mat;      //!<  the matrix
  int64_t start;
  int64_t idx;
  int64_t stop;
  int64_t col;
} next_nz_t;

next_nz_t * new_nxt_nz( sparsemat_t *mat );
void init_nxt_l_nz(next_nz_t * nxtnz, int64_t row);
int has_nxt_l_nz(next_nz_t * nxtnz);
void incr_nxt_l_nz(next_nz_t * nxtnz);

int64_t * rand_perm(int64_t N, int64_t seed);

sparsemat_t * permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv);
sparsemat_t * transpose_matrix(sparsemat_t *A);

/* misc utility functions */
double wall_seconds();
int64_t dump_array(int64_t *A, int64_t len, int64_t maxdisp, char * name);
int64_t dump_matrix(sparsemat_t * A, int64_t maxrows, char * name);
int64_t write_matrix_mm(sparsemat_t * A, char * name);
sparsemat_t * read_matrix_mm(char * name);

void spmat_stats(sparsemat_t *mat);

int64_t is_upper_triangular(sparsemat_t *A);
int64_t is_lower_triangular(sparsemat_t *A);
int64_t is_perm(int64_t * perm, int64_t N);

int nz_comp(const void *a, const void *b);
int64_t sort_nonzeros( sparsemat_t *mat);
int64_t compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);

sparsemat_t * copy_matrix(sparsemat_t *srcmat);

sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz);
enum ER_TRIANGLE {ER_TRI_L = 0, ER_TRI_U = 1, ER_TRI_LWD = 2, ER_TRI_UWD = 3};

sparsemat_t * random_graph(int64_t n, graph_gen mode, edge_type type, self_loops loops,
			   double edge_density, int64_t seed);

sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type type, self_loops loops, int64_t seed);
sparsemat_t * erdos_renyi_random_graph(int n, double p, edge_type type, self_loops loops, int64_t seed);
sparsemat_t * naive_erdos_random_graph(int n, double p, edge_type type, self_loops loops, int64_t seed);


void clear_matrix(sparsemat_t * mat);

#define spmat_INCLUDED
#endif

