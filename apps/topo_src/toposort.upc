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

/*! \file toposort.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"

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

Usage:
topo [-h][-b count][-M mask][-n num][-f filename][-Z num][-e prob][-D]
- -h prints this help message
- -b count is the number of packages in an exstack(2) buffer
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate
- -n num is the number of rows per thread
- -f filename read the input matrix from filename (in Matrix Market format)
- -Z num use an Erdos Renyi matrix with num being the expected number of nonzeros in a row 
- -e prob use an Erdos Renyi matrix where prob is the probability of an entry in matrix being non-zero 
- -D debugging flag that dumps out input and output files.
*/

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
topo [-h][-b count][-M mask][-n num][-f filename][-Z num][-e prob][-D]\n\
 -h prints this help message\n\
 -b count is the number of packages in an exstack(2) buffer\n\
 -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate\n\
 -n num is the number of rows per thread\n\
 -f filename read the input matrix from filename (in Matrix Market format)\n\
 -F Set the graph model to FLAT\n\
 -G Set the graph model to GEOMETRIC\n\
 -Z num use an Erdos Renyi matrix with num being the expected number of nonzeros in a row \n\
 -e prob use an Erdos Renyi matrix where prob is the probability of an entry in matrix being non-zero \n\
 -D debugging flag that dumps out input and output files.\n\
\n");
  lgp_finalize();
  lgp_global_exit(0);
}


/*! \brief check the result toposort 
 *
 * check that the permutations are in fact permutations and the check that applying
 * them to the original matrix yields an upper triangular matrix
 * \param mat the original matrix
 * \param rperminv the row permutation
 * \param cperminv the column permutation
 * \param dump_files debugging flag
 * \return 0 on success, 1 otherwise
 */
int check_is_triangle(sparsemat_t * mat, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t dump_files) {
  sparsemat_t * mat2;
  int ret = 0;

  int rf = is_perm(rperminv, mat->numrows);
  int cf = is_perm(cperminv, mat->numcols);
  if(!rf || !cf){
    T0_fprintf(stderr,"ERROR: check_is_triangle is_perm(rperminv2) = %d is_perm(cperminv2) = %d\n",rf,cf);
    return(1);
  }
  mat2 = permute_matrix(mat, rperminv, cperminv);
  if(!is_upper_triangular(mat2, 1)) {
    T0_fprintf(stderr,"ERROR: check_is_triangle fails\n");
    ret = 1;
  }
  clear_matrix(mat2);
  free(mat2);
  return(ret);
}

/*! \brief Generates an input matrix for the toposort algorithm
 * \param numrows the number of rows (and columns) in the produced matrix
 * \param prob the probability that there is an edge between two given vertices, i.e. the probability
 *   that a given entry in the matrix is non-zero.
 * \param rand_seed the seed for random number generator that determines the original matrix and the permutations
 * \return the permuted upper triangular matrix
 */
sparsemat_t * generate_toposort_input(int64_t numrows, graph_model model, double edge_prob, uint64_t rand_seed) {
  sparsemat_t * omat;
  int64_t numcols = numrows;
  
  T0_fprintf(stderr,"Creating input matrix for toposort\n");fflush(stderr);
  double t = wall_seconds();
  omat = random_graph(numrows, model, UNDIRECTED, LOOPS, edge_prob, rand_seed);
  T0_printf("generate graph time %lf\n", wall_seconds() - t);
  if(!omat) {printf("Error! omat is NULL");return(NULL);}
  if(!is_lower_triangular(omat, 1)){printf("ERROR: Not LT\n");return(NULL);}

  sparsemat_t * U = transpose_matrix(omat);
  if(!U){printf("Error! U is NULL");return(NULL);}
  if(!is_upper_triangular(U, 1)){printf("ERROR: Not UT\n");return(NULL);}
  
  // get row and column permutations
  t = wall_seconds();
  SHARED int64_t * rperminv = rand_permp(numrows, 1230+MYTHREAD);
  SHARED int64_t * cperminv = rand_permp(numcols, 45+MYTHREAD);
  T0_printf("generate perms time %lf\n", wall_seconds() - t);
  lgp_barrier();

  if(!rperminv || !cperminv){
    T0_printf("ERROR: topo_rand_permp returns NULL!\n");fflush(0);
    return(NULL);
  }
  
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * lcperminv = lgp_local_part(int64_t, cperminv);

  lgp_barrier();
  t = wall_seconds();
  sparsemat_t * mat = permute_matrix(U, rperminv, cperminv);
  T0_printf("permute matrix time %lf\n", wall_seconds() - t);
  
  if(!mat) {
    T0_printf("ERROR: permute_matrix returned NULL");fflush(0);
    return(NULL);
  }

  lgp_barrier();

  clear_matrix( omat );free(omat);
  clear_matrix(U); free(U);
  lgp_all_free(rperminv);
  lgp_all_free(cperminv);

  return( mat );
}

int main(int argc, char * argv[]) {
  lgp_init(argc, argv);

  int64_t i, j, fromth, lnnz, start, end;  
  int64_t pe, row, col, idx;
  double t1;

  int64_t l_numrows = 100000;
  double  nz_per_row = 35;
  int64_t buf_cnt = 1024;
  int64_t rand_seed =  MYTHREAD*MYTHREAD*10000 + 5;
  int64_t numrows, numcols;
  int64_t pos = 0;

  double erdos_renyi_prob = 0.0;
  int64_t models_mask = ALL_Models;
  int64_t printhelp = 0;
  int64_t read_graph = 0;
  char filename[64];
  int64_t dump_files = 0;
  int64_t cores_per_node = 1;
  graph_model model = FLAT;
  
  int opt; 
  while( (opt = getopt(argc, argv, "hb:c:M:n:Ff:GZ:p:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%"PRId64"" , &buf_cnt);  break;
    case 'c': sscanf(optarg,"%"PRId64"" ,&cores_per_node); break;
    case 'M': sscanf(optarg,"%"PRId64"" , &models_mask);  break;
    case 'F': model = FLAT; break;
    case 'G': model = GEOMETRIC; break;
    case 'n': sscanf(optarg,"%"PRId64"" , &l_numrows);  break;
    case 'f': read_graph = 1; sscanf(optarg,"%s", filename); break;
    case 'Z': sscanf(optarg,"%lf" , &nz_per_row);  break;
    case 'e': sscanf(optarg,"%lf" , &erdos_renyi_prob);  break;
    case 'D': dump_files = 1; break;
    default:  break;
    }
  }
  if(printhelp) usage();

  numrows = l_numrows * THREADS;
  numcols = numrows;
  if(erdos_renyi_prob == 0.0){ // use nz_per_row to get erdos_renyi_prob
    erdos_renyi_prob = (2.0*nz_per_row)/(numrows - 1);
    if(erdos_renyi_prob > 1.0)
      erdos_renyi_prob = 1.0;
  } else {                     // use erdos_renyi_prob to get nz_per_row
    nz_per_row = erdos_renyi_prob * numrows;
  }

  T0_fprintf(stderr,"Running toposort on %d threads\n", THREADS);
  T0_fprintf(stderr,"buf_cnt (stack size)           (-b)   %"PRId64"\n", buf_cnt);
  T0_fprintf(stderr,"Number of rows per thread      (-n)   %"PRId64"\n", l_numrows);
  T0_fprintf(stderr,"Avg # of nonzeros per row      (-Z)   %2.2lf\n", nz_per_row);
  T0_fprintf(stderr,"Erdos-Renyi edge probability   (-e)   %lf\n", erdos_renyi_prob);
  T0_fprintf(stderr,"Graph Model              (-F or -G)   %s\n", (model == FLAT ? "FLAT" : "GEOMETRIC"));
  T0_fprintf(stderr,"task mask (M) = %"PRId64" (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);
  
  sparsemat_t * mat = generate_toposort_input(numrows, model, erdos_renyi_prob, rand_seed);
  if(!mat){T0_printf("ERROR: mat is NULL!\n"); lgp_global_exit(1);}

  T0_printf("Input matrix has %"PRId64" rows and %"PRId64" nonzeros\n", mat->numrows, mat->nnz);
  
  sparsemat_t * tmat = transpose_matrix(mat);
  if(!tmat){T0_printf("ERROR: tmat is NULL!\n"); lgp_global_exit(1);}

  lgp_barrier();

  T0_fprintf(stderr,"Run toposort on mat (and tmat) ...\n");
  // arrays to hold the row and col permutations
  SHARED int64_t *rperminv2 = lgp_all_alloc(numrows, sizeof(int64_t));
  SHARED int64_t *cperminv2 = lgp_all_alloc(numcols, sizeof(int64_t));
  double gb_th  = (mat->numrows + mat->numcols*2 + mat->nnz*2)*8;

  int64_t use_model;
  double laptime = 0.0;

  for( use_model=1L; use_model < 32; use_model *=2 ) {

    switch( use_model & models_mask ) {
    case AGI_Model:
      T0_fprintf(stderr,"      AGI: ");
      laptime = toposort_matrix_agi(rperminv2, cperminv2, mat, tmat);
      break;

    case EXSTACK_Model:
      T0_fprintf(stderr,"  Exstack: ");
      laptime = toposort_matrix_exstack(rperminv2, cperminv2, mat, tmat);
      break;

    case EXSTACK2_Model:
      T0_fprintf(stderr," Exstack2: ");
      laptime = toposort_matrix_exstack2(rperminv2, cperminv2, mat, tmat);
      break;

    case CONVEYOR_Model:
      T0_fprintf(stderr," Conveyor: ");
      laptime = toposort_matrix_convey(rperminv2, cperminv2, mat, tmat);
      break;

    case ALTERNATE_Model:
      //T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      T0_fprintf(stderr," ALTERNATE: ");
      laptime = toposort_matrix_cooler(rperminv2, cperminv2, mat, tmat);
      break;
    
    default:
      continue;
      
    }
    lgp_barrier();
    T0_fprintf(stderr,"  %8.3lf seconds\n", laptime);

    if( check_is_triangle(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_upc: mat2 is not upper-triangular!\n");
    }
  }

  lgp_barrier();
  lgp_finalize();
  return(0);  
}

