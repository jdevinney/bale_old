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
#define spmat_utils_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <fcntl.h>
#include <getopt.h>
#include <assert.h>


/*! \struct sparsemat_t 
 * \brief A structure to hold a sparse matrix 
 * \ingroup spmatgrp
 * We use the standard Compressed Sparse Row format.
 * The offset array is an array with length numrows + 1. offset[i]
 * is the index where the ith row's data starts in the nonzero and value arrays. 
 * The nonzero array holds the column index for each nonzero in row-major order.
 * The value array holds the nonzero value for each nonzero in row-major order 
 * (it is aligned with the nonzero array.)
 *
*/
typedef struct sparsemat_t{
  int64_t numrows;       //!< the total number of rows in the matrix
  int64_t numcols;       //!< the nonzeros have values between 0 and numcols
  int64_t nnz;           //!< total number of nonzeros in the matrix
  int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * nonzero;     //!< the global array of column indices for nonzeros
  double * value;       //!< the global array of nonzero values (optional)
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

/*! \struct element_t 
 * \brief structure used while reading and writing the MatrixMarket format.
 * We only handle {0,1} matrices, so we don't need triples.
 */

typedef struct edge_t{
  int64_t row;
  int64_t col;
}edge_t;

typedef struct edge_list_t{
  edge_t * edges;
  int64_t nalloc;
  int64_t num;
}edge_list_t;

// struct to sort rows in a matrix with values
typedef struct col_val_t{
   int64_t col;
   double value;
 }col_val_t;

typedef struct triple_t{
  int64_t row;
  int64_t col;
  double val;
}triple_t;

typedef struct kron_args_t{
  char str[256];
  int64_t mode;
  int64_t num_stars;      // 
  int64_t star_size[64];  // can't be too many stars, else the graph would be huge
  int64_t numrows;
} kron_args_t;

typedef enum graph_model {FLAT, GEOMETRIC, KRONECKER} graph_model;
typedef enum edge_type {DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED} edge_type;
typedef enum self_loops {LOOPS, NOLOOPS} self_loops;

next_nz_t * new_nxt_nz( sparsemat_t *mat );
void init_nxt_l_nz(next_nz_t * nxtnz, int64_t row);
int has_nxt_l_nz(next_nz_t * nxtnz);
void incr_nxt_l_nz(next_nz_t * nxtnz);



void             clear_matrix(sparsemat_t * mat);
void             clear_kron_args(kron_args_t * kron_args);
int64_t          compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);
sparsemat_t *    copy_matrix(sparsemat_t *srcmat);

int64_t          dump_array(int64_t *A, int64_t len, int64_t maxdisp, char * name);
int64_t          dump_matrix(sparsemat_t * A, int64_t maxrows, char * name);

sparsemat_t *    erdos_renyi_random_graph(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed);
sparsemat_t *    erdos_renyi_random_graph_naive(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed);

sparsemat_t *    geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed);
kron_args_t *    kron_args_init(char *str);
sparsemat_t *    kronecker_product_graph(kron_args_t * K);
int64_t          tri_count_kron_graph(kron_args_t *K);

sparsemat_t *    init_matrix(int64_t numrows, int64_t numcols, int64_t nnz, int values);

int64_t          is_upper_triangular(sparsemat_t *A);
int64_t          is_lower_triangular(sparsemat_t *A);
int64_t          is_perm(int64_t * perm, int64_t N);

int              nz_comp(const void *a, const void *b);
sparsemat_t *    permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv);

int64_t *        rand_perm(int64_t N, int64_t seed);

sparsemat_t *    random_graph(int64_t n, graph_model model, edge_type edge_type, self_loops loops, double edge_density, int64_t seed);
sparsemat_t *    random_sparse_matrix(int64_t nrows, int64_t ncols, double density, int values, int64_t seed);

sparsemat_t *    read_matrix_mm(char * name);
void             resolve_edge_prob_and_nz_per_row(double * edge_prob, double * nz_per_row,
                                                  int64_t numrows, self_loops loops);
int64_t          sort_nonzeros( sparsemat_t *mat);
void             spmat_stats(sparsemat_t *mat);

double           sssp_dijsktra_linear(sparsemat_t * mat, double *dist, int64_t v0);
double           sssp_dijsktra_heap(sparsemat_t * mat, double *dist, int64_t r0);
double           sssp_delta_stepping(sparsemat_t * mat, double *dist, int64_t r0);

sparsemat_t *    transpose_matrix(sparsemat_t *A);
sparsemat_t *    make_symmetric_from_lower(sparsemat_t * L);
int64_t          write_matrix_mm(sparsemat_t * A, char * name);



double wall_seconds();
#define DEBUG 0
#define Dprintf if(DEBUG) printf
#endif

