/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
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

/*! \file toposort.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"
#include <std_options.h>

/*!
\page toposort_page Toposort

The toposort algorithm is more complicated than histogram and indexgather.
It typically enjoys a significant amount of parallelism, but it 
is not completely order and latency tolerant.

To prepare the input to the toposort algorithm, 
we start we an upper-triangular matrix <b>T</b> with no zeros on the diagonal. 
We don't care about the values of the non-zeros in the matrix, only their position.
Next we randomly permute the rows and columns of <b>T</b> to get a matrix <b>M</b>. 
The matrix <b>M</b> has been called a morally triangular matrix.

Given a morally triangular matrix, <b>M</b>, the goal of toposort is 
to create row and column permutations such that when these 
permutations are applied to <b>M</b>,
the result is an upper triangular matrix with no zeros on the diagonal.
Note, the answer need not be unique and since <b>M</b> is a row and column permutation
of <b>T</b>, there must be a solution.

We use a breadth first search algorithm based on the following observations:
    - The rows (and columns) of a sparse matrix partition the set of non-zeros in the matrix.
    - Row and column permutations preserve both partitions.

For example, if you create a set of the nonzeros in a particular row;
a column permutation might change the labels of the elements in the set,
but doesn't change the cardinality of the set. Likewise for columns.
Hence, there must be a row in <b>M</b> with a single non-zero.
If we remove that row and column, 
we are left with smaller, morally upper triangular matrix. This is the motivation behind a simple algorithm.

The outline of the algorithm is as follows:
\verbatim
For all rows with a single non-zero, put its non-zero onto a queue.
While the queue is not empty: 
   pop a non-zero from the queue (this represents row r and column c)
   claim its new position as the last row and column of the permutations being created 
   remove all the non-zeros in column c
   if any row now has a single non-zero, enqueue that non-zero
\endverbatim

Rather than changing the matrix by deleting rows and column and then searching the 
new matrix for the next row.  We do the obvious thing of keeping and array of row counts,
<b>rowcnt[i]</b> is the number of non-zeros in <b>row i</b> and
we use a cool trick to find the column of a row with <b>rowcnt[i]</b> equal 1.
We initialize an array, <b>rowsum[i]</b>, to be the sum of the column indices in <b> row i</b>.
When we "delete" a column we decrement <b>rowcnt[i]</b> and <b>rowsum[i]</b> by that column index.
Hence, when the <b>rowcnt[i]</b> gets down to one, the <b>rowsum[i]</b> is the column that is left.

In parallel there are three race conditions or synchronization issues to address..

The first is reading and writing the queue of rows to be processed.
One way to handle it is to introduce the notion of a levels.
Within a level all threads process the all the rows on their queues 
and by doing so create new degree one rows. These rows are placed on the 
appropriate queues for the next level. There is a barrier between levels.

Threads race to pick their position in <b>rperm</b> and <b>cperm</b>. 
One could handle this race for the pivots with a fetch_and_add,
instead we use parallel prefix to claim enough room for the pivots 
in the current level on each thread then assign them in order per thread.
   
Threads race to update the <b>rowcnt</b> and <b>rowsum</b> arrays. 
We handle this with levels and atomic memory operations.

Run with the --help, -?, or --usage flags for run details.
*/


/*! \brief check the result toposort 
 *
 * check that the permutations are in fact permutations and the check that applying
 * them to the original matrix yields an upper triangular matrix
 * \param mat the original matrix
 * \param rperm the row permutation
 * \param cperm the column permutation
 * \return 0 on success, 1 otherwise
 */
int check_is_triangle(sparsemat_t * mat, SHARED int64_t * rperm, SHARED int64_t * cperm) {
  int ret = 0;

  int rf = is_perm(rperm, mat->numrows);
  int cf = is_perm(cperm, mat->numrows);
  if(!rf || !cf){
    T0_fprintf(stderr,"ERROR: check_is_triangle is_perm(rperm) = %d is_perm(cperm) = %d\n",rf,cf);
    return(1);
  }
  sparsemat_t * mat2 = permute_matrix(mat, rperm, cperm);
  if(!mat2){
    T0_fprintf(stderr,"ERROR: check_is_triangle mat2 is NULL\n");
    return(1);
  }
  if(!is_upper_triangular(mat2, 1)) {
    T0_fprintf(stderr,"ERROR: check_is_triangle fails\n");
    ret = 1;
  }
  clear_matrix(mat2);
  free(mat2);
  return(ret);
}

/*! \brief Generates an input matrix for the toposort algorithm from a unit-triangulare matrix.
 * \param tri_mat An upper (or lower) triangular matrix with unit diagonal.
 * \param rand_seed the seed for random number generator.
 * \return a permuted upper triangular matrix
 */
sparsemat_t * generate_toposort_input(sparsemat_t * tri_mat, uint64_t rand_seed) {
  sparsemat_t * mat = NULL;
  double t;
  //write_matrix_mm(tri_mat, "tri_mat");
  if(!is_upper_triangular(tri_mat, 1)){
    if(is_lower_triangular(tri_mat, 1)){
      mat = transpose_matrix(tri_mat);
    }else{
      T0_fprintf(stderr, "ERROR: toposort: input matrix is not triangular with unit diagonal!\n");
      return(NULL);
    }
  }else{
    mat = tri_mat;
  }

  if(!mat){T0_printf("ERROR: mat is NULL!\n"); return(NULL);}

  // get row and column permutations
  t = wall_seconds();
  SHARED int64_t * rperm = rand_permp(mat->numrows, rand_seed);
  SHARED int64_t * cperm = rand_permp(mat->numrows, rand_seed + 12345);
  //T0_printf("generate perms time %lf\n", wall_seconds() - t);
  lgp_barrier();

  
  if(!rperm || !cperm){
    T0_printf("ERROR: topo_rand_permp returns NULL!\n");fflush(0);
    return(NULL);
  }
  
  lgp_barrier();
  t = wall_seconds();
  sparsemat_t * pmat = permute_matrix(mat, rperm, cperm);
  if(!pmat) {
    T0_printf("ERROR: permute_matrix returned NULL");fflush(0);
    return(NULL);
  }
  //T0_printf("permute matrix time %lf\n", wall_seconds() - t);

  lgp_barrier();
  if(mat != tri_mat){
    clear_matrix( mat );free(mat);
  }
  
  lgp_all_free(rperm);
  lgp_all_free(cperm);
  
  return( pmat );
}

typedef struct args_t{
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };


int main(int argc, char * argv[]) {

  /* process command line */
  args_t args = {0}; // initialize args struct to all zero
  struct argp argp = {NULL, parse_opt, 0,
                      "Parallel topological sort.", children_parsers};

  args.gstd.loops = 1; // force loops into graph
  args.gstd.l_numrows = 500000;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  /* force input graph to be undirected and to have loops, no matter what the options */
  if(args.gstd.directed == 1){
    T0_fprintf(stderr, "toposort needs undirected graph input, overriding -d flag.\n");
    args.gstd.directed = 0;
  }
    
  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }
  
  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_printf("Error! toposort: inmat is NULL");lgp_global_exit(-1);}
  
  // permate the rows and columns of the matrix randomly
  sparsemat_t * mat = generate_toposort_input(inmat, args.std.seed + 23456);
  if(!mat){T0_printf("Error! toposort: mat is NULL");lgp_global_exit(-1);}

  // get the transpose of mat (needed for toposort implemmentations)
  sparsemat_t * tmat = transpose_matrix(mat);
  if(!tmat){T0_printf("Error! toposort: tmat is NULL");lgp_global_exit(-1);}


  if(args.std.dump_files){
    write_matrix_mm(inmat, "topo_inmat");
    write_matrix_mm(mat, "topo_permuted_mat");
  }

  lgp_barrier();
  clear_matrix(inmat); free(inmat);


  // arrays to hold the row and col permutations
  SHARED int64_t *rperm2 = lgp_all_alloc(mat->numrows, sizeof(int64_t));
  SHARED int64_t *cperm2 = lgp_all_alloc(mat->numrows, sizeof(int64_t));

  int64_t use_model;
  double laptime = 0.0;
  char model_str[32];
  for( use_model=1L; use_model < 32; use_model *=2 ) {

    switch( use_model & args.std.models_mask ) {
    case AGP_Model:
      sprintf(model_str, "AGP");
      laptime = toposort_matrix_agp(rperm2, cperm2, mat, tmat);
      break;

    case EXSTACK_Model:
      sprintf(model_str, "Exstack");
      laptime = toposort_matrix_exstack(rperm2, cperm2, mat, tmat, args.std.buffer_size);
      break;

    case EXSTACK2_Model:
      sprintf(model_str, "Exstack2");
      laptime = toposort_matrix_exstack2(rperm2, cperm2, mat, tmat, args.std.buffer_size);
      break;

    case CONVEYOR_Model:
      sprintf(model_str, "Conveyor");
      laptime = toposort_matrix_convey(rperm2, cperm2, mat, tmat);
      break;

    case ALTERNATE_Model:
      //T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      sprintf(model_str, "Alternate");
      laptime = toposort_matrix_cooler(rperm2, cperm2, mat, tmat);
      break;
    
    default:
      continue;
      
    }
    lgp_barrier();

    bale_app_write_time(&args.std, model_str, laptime);

    if( check_is_triangle(mat, rperm2, cperm2) ) {
      T0_fprintf(stderr,"\nERROR: After toposort_matrix_upc: mat2 is not upper-triangular!\n");
    }
  }

  lgp_barrier();

  bale_app_finish(&args.std);
  
  return(0);  
}

