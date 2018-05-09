/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
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
  int64_t numrows;              //!< the total number of rows in the matrix
  int64_t lnumrows;             //!< the number of rows on this thread 
                                // note lnumrows = (numrows / THREADS) + {0 or 1} 
                                //    depending on numrows % THREADS
  int64_t numcols;              //!< the nonzeros have values between 0 and numcols
  int64_t nnz;                  //!< total number of nonzeros in the matrix
  int64_t lnnz;                 //!< the number of nonzeros on this thread
  SHARED int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * loffset;            //!< the row offsets for the row with affinity to this thread
  SHARED int64_t * nonzero;     //!< the global array of nonzeros
  int64_t * lnonzero;           //!< the nonzeros with affinity to this thread

  int64_t l_enum_row;           // state to play with row iterator stuff
  int64_t l_enum_idx;
  int64_t l_enum_nstop;

  int64_t S_enum_row;
  int64_t S_enum_idx;
  int64_t S_enum_nstop;

}sparsemat_t;

SHARED int64_t * rand_permp_conveyor(int64_t N, int seed);
SHARED int64_t * rand_permp_exstack2(int64_t N, int seed);
SHARED int64_t * rand_permp_exstack(int64_t N, int seed);
SHARED int64_t * rand_permp_atomic(int64_t N, int seed);

sparsemat_t * permute_matrix_conveyor(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);
sparsemat_t * permute_matrix_exstack2(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);
sparsemat_t * permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);
sparsemat_t * permute_matrix_atomic(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);

sparsemat_t * transpose_matrix_conveyor(sparsemat_t * A);
sparsemat_t * transpose_matrix_exstack2(sparsemat_t * A);
sparsemat_t * transpose_matrix_exstack(sparsemat_t * A);
sparsemat_t * transpose_matrix_atomic(sparsemat_t * A);


/* wrapper functions */
SHARED int64_t * rand_permp(int64_t N, int seed, int64_t model);
sparsemat_t * permute_matrix(sparsemat_t *omat, SHARED int64_t *rperminv, SHARED int64_t *cperminv, int64_t model);
sparsemat_t * transpose_matrix(sparsemat_t *omat, int64_t model );


/* misc utility functions */
int write_matrix(sparsemat_t * A, int maxrows, char * name);

int is_upper_triangular(sparsemat_t *A);
int is_lower_triangular(sparsemat_t *A);
int is_perm(SHARED int64_t * perm, int64_t N);

int nz_comp(const void *a, const void *b);
sparsemat_t * gen_uniform_sparse(int64_t rowsperthread, int64_t max_row_cnt);

int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);
sparsemat_t * copy_matrix(sparsemat_t *srcmat);
int sort_nonzeros( sparsemat_t *mat);

sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread);
void clear_matrix(sparsemat_t * mat);

#define spmat_INCLUDED
#endif

