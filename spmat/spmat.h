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
/*! \file spmat.h
 * \brief The header file for spmat library.
 * \ingroup spmatgrp
 */ 
#ifndef spmat_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <convey.h>
#include <libgetput.h>


/*! \struct sparsemat_t spmat.h
 * \brief A structure to hold a sparse matrix.
 *
 * We use a distributed version of the standard Compressed Sparse Row
 * (CSR) format.  Since the applications in bale (so far) only need to
 * know whether a matrix entry is zero or nonzero, we don't keep track
 * of the values themselves for the nonzeros in the matrix only the
 * position of the nonzeros.  This reduces the amount of memory needed
 * to store the matrices.  This also saves a lot local sorting and
 * combining that would be required to actually handle values (other
 * than one).
 
 * We store the nonzeros with affinity by rows.  That is, if row i has
 * affinity to PE t, then all the nonzeros in row i have affinity to
 * PE t. The affinity rule for a row is as follows. Row i is
 * assigned to PE (i % NPES).
 *
 * We call a matrix "tidy" if the column indices of nonzeros in every
 * row are sorted in ascending order.
 *
 *  \ingroup spmatgrp
 */
typedef struct sparsemat_t {
  int64_t local;                //!< 0/1 flag specifies whether this is a local or distributed matrix
  int64_t numrows;              //!< the total number of rows in the matrix
  int64_t lnumrows;             //!< the number of rows on this PE 
                                // note lnumrows = (numrows / NPES) + {0 or 1} 
                                //    depending on numrows % NPES
  int64_t numcols;              //!< the nonzeros have values between 0 and numcols
  int64_t nnz;                  //!< total number of nonzeros in the matrix
  int64_t lnnz;                 //!< the number of nonzeros on this PE
  SHARED int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * loffset;            //!< the row offsets for the row with affinity to this PE
  SHARED int64_t * nonzero;     //!< the global array of nonzeros
  int64_t * lnonzero;           //!< the nonzeros with affinity to this PE
}sparsemat_t;


/*! \struct nxnz_t spmat.h
 * \brief A structure to experiment with an iterator that walks across a row of a sparsemat.
 * \ingroup spmatgrp
 */
typedef struct nxnz_t {  //!< next nonzero struct
  sparsemat_t * mat;     //!< the matrix in question
  int64_t row;           //!< the row of the matrix, used to init and to check for consistent usage.
  int64_t first;         //!< the index of the first nonzero in the row
  int64_t idx;           //!< the current idx
  int64_t stop;          //!< the index of the first nonzero in the next row
  int64_t col;           //!< a place to put the nonzero (ie the column of the nonzero)
  //int64_t val;         //   a place to put the value of the nonzero (if we used them).
}nxnz_t;

nxnz_t * init_nxnz(sparsemat_t * mat);
void first_l_nxnz( nxnz_t *nxz, int64_t l_row );
bool has_l_nxnz( nxnz_t *nxz, int64_t l_row );
void incr_l_nxnz( nxnz_t *nxz, int64_t l_row );

void first_S_nxnz( nxnz_t *nxz, int64_t S_row );
bool has_S_nxnz( nxnz_t *nxz, int64_t S_row );
void incr_S_nxnz( nxnz_t *nxz, int64_t S_row );

int64_t rowcount_l( sparsemat_t *mat, int64_t l_row );
int64_t rowcount_S( sparsemat_t *mat, int64_t S_row );

SHARED int64_t * rand_permp_conveyor(int64_t N, int seed);
SHARED int64_t * rand_permp_exstack2(int64_t N, int seed, int64_t buf_cnt);
SHARED int64_t * rand_permp_exstack(int64_t N, int seed, int64_t buf_cnt);
SHARED int64_t * rand_permp_agi(int64_t N, int seed);

sparsemat_t * permute_matrix_conveyor(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);
sparsemat_t * permute_matrix_exstack2(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt);
sparsemat_t * permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt);
sparsemat_t * permute_matrix_agi(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);

sparsemat_t * transpose_matrix_conveyor(sparsemat_t * A);
sparsemat_t * transpose_matrix_exstack2(sparsemat_t * A, int64_t buf_cnt);
sparsemat_t * transpose_matrix_exstack(sparsemat_t * A, int64_t buf_cnt);
sparsemat_t * transpose_matrix_agi(sparsemat_t * A);

int64_t write_sparse_matrix_agi( char * datadir, sparsemat_t * mat);
int64_t write_sparse_matrix_exstack( char * datadir, sparsemat_t * mat, int64_t buf_cnt);

/* wrapper functions */
SHARED int64_t * rand_permp(int64_t N, int seed);
sparsemat_t * permute_matrix(sparsemat_t *omat, SHARED int64_t *rperminv, SHARED int64_t *cperminv);
sparsemat_t * transpose_matrix(sparsemat_t *omat);


/* misc utility functions */
//int write_matrix(sparsemat_t * A, int maxrows, char * name);
int write_matrix_mm(sparsemat_t * A, char * name);
sparsemat_t * read_matrix_mm_to_dist(char * name);
int64_t write_sparse_matrix_metadata(char * dirname, sparsemat_t * A);
int64_t read_sparse_matrix_metadata(char * dirname, int64_t * nr, int64_t * nc, int64_t * nnz, int64_t *nwriters);
sparsemat_t * read_sparse_matrix_agi(char * datadir);


int64_t tril(sparsemat_t * A, int64_t k);
int64_t tril(sparsemat_t * A, int64_t k);
int is_upper_triangular(sparsemat_t *A, int64_t unit_diagonal);
int is_lower_triangular(sparsemat_t *A, int64_t unit_diagonal);
int is_perm(SHARED int64_t * perm, int64_t N);

int nz_comp(const void *a, const void *b);
sparsemat_t * gen_erdos_renyi_graph_dist_naive(int n, double p, int64_t unit_diag, int64_t mode, int64_t seed);
sparsemat_t * gen_erdos_renyi_graph_dist(int n, double p, int64_t unit_diag, int64_t mode, int64_t seed);
sparsemat_t * gen_erdos_renyi_graph_triangle_dist(int n, double p, int64_t unit_diag, int64_t lower, int64_t seed);

sparsemat_t * kron_prod_dist(sparsemat_t * B, sparsemat_t * C, int64_t lower);
sparsemat_t * kron_prod(sparsemat_t * B, sparsemat_t * C);
sparsemat_t * gen_star(int64_t m, int mode);
sparsemat_t * gen_local_mat_from_stars(int64_t M, int64_t * m, int mode);

int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);
sparsemat_t * copy_matrix(sparsemat_t *srcmat);
int sort_nonzeros( sparsemat_t *mat);

sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread);
sparsemat_t * init_local_matrix(int64_t numrows, int64_t numcols, int64_t nnz);

void clear_matrix(sparsemat_t * mat);

#define spmat_INCLUDED
#endif

