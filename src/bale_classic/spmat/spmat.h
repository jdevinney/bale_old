/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//
//  All rights reserved.
//  
//  This file is a part of Bale.  For license information see the
//  LICENSE file in the top level directory of the distribution.
//  
 *****************************************************************/ 
/*! \file spmat.h
 * \brief The header file for spmat library.
 * \ingroup spmatgrp
 */ 
#ifndef spmat_INCLUDED
#define spmat_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <convey.h>
#include <libgetput.h>
#include "spmat_enums.h"


/*! \struct sparsemat_t spmat.h
 * \brief A structure to hold a sparse matrix.
 *
 * We use a distributed version of the standard Compressed Sparse Row (CSR) format.  
 * In most cases we don't have values in the sparse matrix.
 * For example, in toposort we only track the positions of the nonzeros,
 * in permute matrix or transpose matrix having values would only double the size
 * of messages. In some apps we store simple graphs in the matrix.
 * In these cases, the value pointer will be NULL.
 * For sssp, we store weighted directed graph in the matrix, so we do need values.
 * We only support values that are doubles.
 *
 * We store the nonzeros with affinity by rows.  That is, if row i has
 * affinity to PE t, then all the nonzeros in row i also have affinity 
 * to PE t. 
 * We used a striped layout for affinity.  
 * That is, row i has affinity to PE (i % NPES).
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
  SHARED int64_t * nonzero;     //!< the global array of column indices of nonzeros
  int64_t * lnonzero;           //!< local array of column indices (for rows on this PE).
  SHARED double * value;        //!< the global array of values of nonzeros. Optional.
  double * lvalue;              //!< local array of values (values for rows on this PE)
}sparsemat_t;


typedef struct d_array_t {
  int64_t num;                  //!< the total number of entries in the array
  int64_t lnum;                 //!< the number of entries on this PE. lnum = (num / NPES) + {0 or 1}
  SHARED double * entry;        //!< the shared array, striped across PE's
  double * lentry;              //!< the localized part of the shared array
} d_array_t;

typedef struct triples_t {
  int64_t * row;
  int64_t * col;
  double * val;
  int64_t numrows;
  int64_t lnumrows;
  int64_t numcols;
  int64_t lnnz;
  int64_t nalloc;
} triples_t;


typedef struct w_edge_t{
  int64_t row;
  int64_t col;
  double val;
}w_edge_t;


typedef struct edge_t{
  int64_t row;
  int64_t col;
}edge_t;

typedef struct edge_list_t{
  edge_t * edges;
  w_edge_t * wedges;
  int64_t nalloc;
  int64_t num;
}edge_list_t;

// struct to sort rows in a matrix with values
typedef struct col_val_t{
   int64_t col;
   double value;
}col_val_t;

// struct to represent a point on the plane. (for geometric graphs)
typedef struct point_t{
  double x;
  double y;
}point_t;



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

int64_t             append_edge(edge_list_t * el, int64_t row, int64_t col);
int64_t             append_weighted_edge(edge_list_t * el, int64_t row, int64_t col, double val);
int64_t             append_triple(triples_t * T, int64_t row, int64_t col, double val);

int64_t             calculate_num_triangles(int kron_mode, int * kron_spec, int kron_num);
void                clear_edge_list(edge_list_t * el);
void                clear_matrix(sparsemat_t * mat);
void                clear_triples(triples_t * T);

int                 compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);
sparsemat_t *       copy_matrix(sparsemat_t *srcmat);

sparsemat_t *       direct_undirected_graph(sparsemat_t * L);

sparsemat_t *       erdos_renyi_random_graph(int64_t n, double p, edge_type edge_type, self_loops loops, uint64_t seed);
sparsemat_t *       gen_star_graph(int64_t m, int mode);
sparsemat_t *       generate_kronecker_graph_from_spec(int mode, int * spec, int num);
sparsemat_t *       geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed, SHARED point_t ** out_points);


edge_list_t *       init_edge_list(int64_t nalloc, int weighted);
sparsemat_t *       init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread, int weighted);
triples_t *         init_triples(int64_t numrows, int64_t numcols, int64_t lnnz, int weighted);
sparsemat_t *       init_local_matrix(int64_t numrows, int64_t numcols, int64_t nnz);

int                 is_upper_triangular(sparsemat_t *A, int64_t unit_diagonal);
int                 is_lower_triangular(sparsemat_t *A, int64_t unit_diagonal);
int                 is_perm(SHARED int64_t * perm, int64_t N);

sparsemat_t *       kronecker_product_of_stars(int64_t M, int64_t * m, int mode);
sparsemat_t *       kronecker_product_graph_local(sparsemat_t * B, sparsemat_t * C);
sparsemat_t *       kronecker_product_graph_dist(sparsemat_t * B, sparsemat_t * C);

sparsemat_t *       permute_matrix(sparsemat_t * A, SHARED int64_t *rperminv, SHARED int64_t *cperminv);
sparsemat_t *       permute_matrix_conveyor(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);
sparsemat_t *       permute_matrix_exstack2(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt);
sparsemat_t *       permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt);
sparsemat_t *       permute_matrix_agp(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv);

void                print_matrix(sparsemat_t * A);

SHARED int64_t *    rand_permp(int64_t N, int seed);
SHARED int64_t *    rand_permp_conveyor(int64_t N, int seed);
SHARED int64_t *    rand_permp_exstack2(int64_t N, int seed, int64_t buf_cnt);
SHARED int64_t *    rand_permp_exstack(int64_t N, int seed, int64_t buf_cnt);
SHARED int64_t *    rand_permp_agp(int64_t N, int seed);
sparsemat_t *       random_graph(int64_t n, graph_model model, edge_type edge_type, self_loops loops,
                                 double edge_density, int64_t seed);
void                resolve_edge_prob_and_nz_per_row(double * edge_prob, double * nz_per_row, int64_t numrows, edge_type edge_type, self_loops loops);

int                 spmat_compare_doubles(double a, double b);
  
sparsemat_t *       transpose_matrix(sparsemat_t * A);
sparsemat_t *       transpose_matrix_conveyor(sparsemat_t * A);
sparsemat_t *       transpose_matrix_exstack2(sparsemat_t * A, int64_t buf_cnt);
sparsemat_t *       transpose_matrix_exstack(sparsemat_t * A, int64_t buf_cnt);
sparsemat_t *       transpose_matrix_agp(sparsemat_t * A);
sparsemat_t *       triples_to_sparsemat(triples_t * T);

int64_t             write_sparse_matrix_agp( char * datadir, sparsemat_t * mat);
int64_t             write_sparse_matrix_exstack( char * datadir, sparsemat_t * mat, int64_t buf_cnt);

/* misc utility functions */
//int write_matrix(sparsemat_t * A, int maxrows, char * name);
int                 write_matrix_mm(sparsemat_t * A, char * name);
sparsemat_t *       read_matrix_mm_to_dist(char * name);
int64_t             write_sparse_matrix_metadata(char * dirname, sparsemat_t * A);
int64_t             read_sparse_matrix_metadata(char * dirname, int64_t * nr, int64_t * nc, int64_t * nnz, int64_t *nwriters);
sparsemat_t *       read_sparse_matrix_agp(char * datadir);


int64_t tril(sparsemat_t * A, int64_t k);
int64_t triu(sparsemat_t * A, int64_t k);


sparsemat_t * gen_erdos_renyi_graph_dist_naive(int n, double p, int64_t unit_diag, int64_t mode, uint64_t seed);
sparsemat_t * gen_erdos_renyi_graph_dist(int n, double p, int64_t unit_diag, int64_t mode, uint64_t seed);
sparsemat_t * gen_erdos_renyi_graph_triangle_dist(int n, double p, int64_t unit_diag, int64_t lower, uint64_t seed);


int sort_nonzeros( sparsemat_t *mat);
int nz_comp(const void *a, const void *b);
int point_comp(const void *a, const void *b);
int col_val_comp(const void *a, const void *b);
//int dbl_comp(const void *a, const void *b);
int edge_comp(const void *a, const void *b);
int w_edge_comp(const void *a, const void *b);

d_array_t * init_d_array(int64_t num); 
d_array_t * copy_d_array(d_array_t *S);
void        set_d_array(d_array_t * A, double v);
int64_t     replace_d_array(d_array_t * D, d_array_t *S);
void        clear_d_array(d_array_t *A);

#endif

