/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

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
#include <time.h>
#include <stdarg.h>
#include <fcntl.h>
#include <getopt.h>
#include <assert.h>

/*! \struct sparsemat_t 
 * \brief The data structure to hold a sparse matrix 
 * \ingroup spmatgrp
 * We use one of the standard Compressed Sparse Row formats.
 * The offset array is an array with length numrows + 1. The offset[] array
 * holds the indices into the arrays nonzero[] and value[].
 * The index offset[i] is where the ith row's data starts in the nonzero[] and value[] arrays. 
 * The nonzero[] array holds the columns number for a particular non-zero in the matrix.
 * The value[] array holds the non-zero value of the corresponding non-zero in the matrix.
 * In case the matrix is a {0,1}-matrix, we don't need values and we set *value=NULL.
*/
typedef struct sparsemat_t{
  int64_t numrows;       //!< the total number of rows in the matrix
  int64_t numcols;       //!< the nonzeros have values between 0 and numcols
  int64_t nnz;           //!< total number of nonzeros in the matrix
  int64_t * offset;      //!< the row offsets into the array of nonzeros
  int64_t * nonzero;     //!< the global array of column indices for nonzeros
  double * value;       //!< the global array of nonzero values (optional)
}sparsemat_t;

/*! \struct w_edge_t  
 * \brief structure to work with weighted edges as a triple.
 */
typedef struct w_edge_t{
  int64_t row;          //!< row
  int64_t col;          //!< col
  double val;           //!< value of M[row,col]
}w_edge_t;

/*! \struct edge_t  
 * \brief structure to work with unweighted edges as a tuple.
 * We take the value of M[row,col] to be 1.
 */
typedef struct edge_t{
  int64_t row;         //!< row
  int64_t col;         //!< col
}edge_t;

/*! \struct edge_list_t  
 * \brief structure to work the non-zeros as a list of edges.
 * This is a convenient format for reading and writing to files
 */
typedef struct edge_list_t{
  edge_t * edges;       //! pointer to a array of edges
  w_edge_t *wedges;     //! pointer to a array of weighted edges
  int64_t nalloc;       //! number of elements allocated for the array
  int64_t num;          //! number of elements in the array
}edge_list_t;

/*! \struct col_val_t  
 * \brief struct use while sorting a row in a matrix with values
 */
typedef struct col_val_t{
   int64_t col;         //!< col
   double value;        //!< val
}col_val_t;

/*! \struct point_  
 * \brief struct to represent a point on the plane. (for geometric graphs)
 */
typedef struct point_t{
  double x;         //!< x
  double y;         //!< y
}point_t;

#if 0  // TODO delete
typedef struct triple_t{
  int64_t row;
  int64_t col;
  double val;
}triple_t;

typedef struct kron_args_t{
  char str[256];
  int mode;
  int num_stars;      // 
  int star_size[64];  // can't be too many stars, else the graph would be huge
  int64_t numrows;
} kron_args_t;
void             clear_kron_args(kron_args_t * kron_args);
kron_args_t *    kron_args_init(char *str);
#endif

typedef enum graph_model {FLAT, GEOMETRIC, KRONECKER} graph_model;
typedef enum edge_type {DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED} edge_type;
typedef enum self_loops {LOOPS, NOLOOPS} self_loops;
//char *graphmodelstr[] = {"FLAT","GEOMETRIC","KRONECKER"};
//char *edgetypestr[] = {"DIRECTED", "UNDIRECTED", "DIRECTED_WEIGHTED", "UNDIRECTED_WEIGHTED"};
//char *loopstr[] = {"LOOPS", "NOLOOPS"};


void             clear_matrix(sparsemat_t * mat);
int64_t          compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat);
sparsemat_t *    copy_matrix(sparsemat_t *srcmat);

int64_t          dump_array(int64_t *A, int64_t len, int64_t maxdisp, char * name);
int64_t          dump_matrix(sparsemat_t * A, int64_t maxrows, char * name);

sparsemat_t *    erdos_renyi_random_graph(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed);
sparsemat_t *    erdos_renyi_random_graph_naive(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed);

sparsemat_t *    geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed);

sparsemat_t *    generate_kronecker_graph_from_spec(int mode, int * spec, int num, int weighted);
//sparsemat_t *    kronecker_product_graph(kron_args_t * K);
sparsemat_t *    kronecker_product_of_stars(int64_t M, int64_t * m, int mode);
sparsemat_t *    kronecker_product_graph_local(sparsemat_t * B, sparsemat_t * C);
sparsemat_t *    kronecker_product_graph_dist(sparsemat_t * B, sparsemat_t * C);
int64_t          tri_count_kron_graph(int kron_mode, int * kron_spec, int kron_num);

//sparsemat_t *    get_input_graph(std_args_t * sargs, std_graph_args_t * gargs);

sparsemat_t *    init_matrix(int64_t numrows, int64_t numcols, int64_t nnz, int values);

int64_t          tril(sparsemat_t * A, int64_t k);
int64_t          triu(sparsemat_t * A, int64_t k);

int64_t          is_upper_triangular(sparsemat_t *A, int64_t ondiag);
int64_t          is_lower_triangular(sparsemat_t *A, int64_t ondiag);
int64_t          is_perm(int64_t * perm, int64_t N);

int              nz_comp(const void *a, const void *b);
sparsemat_t *    permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv);

int64_t *        rand_perm(int64_t N, int64_t seed);

sparsemat_t *    random_graph(int64_t n, graph_model model, edge_type edge_type, self_loops loops, double edge_density, int64_t seed);
//sparsemat_t *    random_sparse_matrix(int64_t nrows, int64_t ncols, double density, int values, int64_t seed);  // TODO delete

sparsemat_t *    read_matrix_mm(char * name);
void             resolve_edge_prob_and_nz_per_row(double * edge_prob, double * nz_per_row,
                                                  int64_t numrows, edge_type edge_type, self_loops loops);
int64_t          sort_nonzeros( sparsemat_t *mat);
void             spmat_stats(sparsemat_t *mat);


sparsemat_t *    transpose_matrix(sparsemat_t *A);
sparsemat_t *    make_symmetric_from_lower(sparsemat_t * L);
int64_t          write_matrix_mm(sparsemat_t * A, char * name);


/*! \struct d_array_t  
 * \brief struct to hold the lenght and the elements in an array of doubles
 * More of a convenience in the serial case, but it is more useful in the parallel case
 */
typedef struct d_array_t {
  int64_t num;                 //!< the total number of entries in the array
  double * entry;              //!< the array of doubles
} d_array_t;

d_array_t * init_d_array(int64_t num); 
d_array_t * read_d_array(char *name);
int64_t     write_d_array(d_array_t *A, char * name);
void        set_d_array(d_array_t * A, double v);
void        clear_d_array(d_array_t *A);
d_array_t  *copy_d_array(d_array_t *src);
int64_t     replace_d_array(d_array_t *dest, d_array_t *src);


double wall_seconds();

#define USE_KNUTH
#ifdef USE_KNUTH
#define CBALE_RAND_MAX 2251799813685248
#include "knuth_rng_double_2019.h"
#else
#define CBALE_RAND_MAX 281474976710656
#endif
double  rand_double();
int64_t rand_int64(int64_t N);
void    rand_seed(int64_t seed);


#define DEBUG 0
#define Dprintf if(DEBUG) printf
#endif

