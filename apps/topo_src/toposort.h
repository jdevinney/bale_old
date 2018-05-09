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

/*! \file toposort.h
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <spmat.h>
#include <locale.h>

double toposort_matrix_exstack(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_exstack_b(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_exstack2(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_upc(SHARED int64_t *tri_rperm, SHARED int64_t *tri_cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_convey(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);

