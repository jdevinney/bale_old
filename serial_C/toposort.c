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

/*! \file toposort.c
 * \brief Demo application that finds an upper triangular form for a matrix.  
 * That is, we are given a matrix that is a random row and column permutation 
 * of a an upper triangular matrix (with ones on the diagonal).
 * This algorithm finds a row and column permutation that would return it
 * to an upper triangular form.
 */

#include "spmat_utils.h"

/*! \page toposort_page Topologically sort a morally upper triangular matrix. 
 
   First we generate the problem by generating an upper triangular matrix
   and applying row and column permutations.
   
   The output of toposort is a row and a column permutation that, if applied,
   would result in an upper triangular matrix.

   We set the row and column permutations,  rperm and cperm, one pivot at a time.

   N = number of rows
   for( pos=N-1; pos > 0; pos-- ) {
     pick a row, r, with a single nonzero, c.
     say (r,c) is the pivot and set rperm[pos] = r and cprem[pos] = c
     Note: a pivot always exists because the matrix is morally upper tri.

     cross out that row r and col c 
   }

   Meaning of cross out:
   Rather than changing the matrix by deleting rows and column and then searching the 
   new matrix for the next pivot.  We do the obvious thing of keeping row counts, where
   rowcnt[i] is the number of non-zeros in row i and we use a really cool trick 
   of keeping the sum of the live column indices for the non-zeros in each row.
   That is, rowsum[i] is the sum of the column indices, not the sum of the non-zero elements,
   for the non-zeros in row i.  To "delete a column" one decrements the rowcnt by one and 
   the rowsum by the corrsponding column index. 
   The cool trick is that, when the rowcnt gets to one, the rowsum is the column that is left.
*/



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
int check_result(sparsemat_t * mat, int64_t * rperminv, int64_t * cperminv, int64_t dump_files) 
{
  sparsemat_t * mat2;
  int ret = 0;
  
  int64_t rf = is_perm(rperminv, mat->numrows);
  int64_t cf = is_perm(cperminv, mat->numcols);
  if(!rf || !cf){
    fprintf(stderr,"ERROR: check_result is_perm(rperminv2) = %ld is_perm(cperminv2) = %ld\n",rf,cf);
    return(1);
  }
  mat2 = permute_matrix(mat, rperminv, cperminv);
  if(is_upper_triangular(mat2) != 0)
    ret = 1;
  if(dump_files) 
    dump_matrix(mat2, 20, "mat2.out");
  clear_matrix(mat2);
  free(mat2);
  return(ret);
}

/*! \brief Generate a matrix that is the a random permutation of a sparse uppper triangular matrix.
 * \param numrows the number of rows (and columns) in the produced matrix
 * \param erdos_renyi_prob the probability that an entry in the matrix is non-zero
 * \param seed the seed for random number generator that determines the original matrix and the permutations
 * \param dump_files is a debugging flag
 * \return the permuted upper triangular matrix
 * 
 * Make the upper triangular matrix is essentially the upper half of an Erdös–Renyi random graph
 * We force the diagonal entry and pick the other edges with probability erdos_renyi_prob
 * then randomly permute the rows and the columns.  
 * The toposort algorithm takes this matrix and finds one of the possibly many row and column permutations 
 *  that would bring the matrix back to an upper triangular form.
 */
sparsemat_t * generate_toposort_input(int64_t numrows, double edge_prob,  uint32_t seed, int64_t dump_files)
{
  int64_t numcols = numrows;
  
  sparsemat_t * tmat = erdos_renyi_tri(numrows, edge_prob, ER_TRI_UWD, seed);
  if(!tmat){
    printf("ERROR: generate_toposort_input: init_matrix failed!\n");
    return(NULL);
  }
  
  if(dump_files) dump_matrix(tmat, 20, "orig_tri.out");
  if(is_upper_triangular(tmat) != 0 ){
    fprintf(stderr,"ERROR: generate_toposort did not start with an upper triangular\n");
    exit(1);
  }
  
  // get random row and column permutations
  int64_t * rperminv = rand_perm(numrows, 1234);
  int64_t * cperminv = rand_perm(numcols, 5678);
  //int64_t * rperminv = rand_perm(numrows, 0);
  //int64_t * cperminv = rand_perm(numcols, 0);
  if(!rperminv || !cperminv){
    printf("ERROR: generate_toposort_input: rand_perm returned NULL!\n");
    exit(1);
  }
  if(dump_files){
    dump_array(rperminv, numrows, 20, "rperm.out");
    dump_array(cperminv, numcols, 20, "cperm.out");
  }
  
  sparsemat_t * mat = permute_matrix(tmat, rperminv, cperminv);
  if(!mat) {
    printf("ERROR: generate_toposort_input: permute_matrix returned NULL");
    exit(1);
  }
  if(dump_files) dump_matrix(mat,20, "perm.out");
  
  clear_matrix( tmat );
  free(tmat);
  free(rperminv);
  free(cperminv);
  
  return(mat);
}


/*!
 * \brief This routine implements the agi variant of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_queue(int64_t *rperm, int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) 
{
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  
  int64_t * queue  = calloc(nr, sizeof(int64_t));
  int64_t * rowtrck = calloc(nr, sizeof(int64_t));
  
  int64_t start, end;
  
  int64_t i, j, row, col, t_row;
   
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  start = end = 0;
  for(i = 0; i < mat->numrows; i++){
    rowtrck[i] = 0L;
    for(j = mat->offset[i]; j < mat->offset[i+1]; j++)
      rowtrck[i] += (1L<<32) + mat->nonzero[j];
    if((rowtrck[i] >> 32) ==  1)
      queue[end++] = i;    
  }
  
  // we a pick a row with a single nonzero = col.
  // setting rperm[pos] = row and cprem[pos] = col
  // in a sense 'moves' the row and col to the bottom
  // right corner of matrix.
  // Next, we cross out that row and col by decrementing 
  //  the rowcnt for any row that contains that col
  // repeat
  //
  double t1 = wall_seconds();
  
  int64_t n_pivots = 0;
  while(start < end){      
    row = queue[start++];
    col = rowtrck[row] & 0xFFFF;  // see cool trick
    
    rperm[row] = nr - 1 - n_pivots;
    cperm[col] = nc - 1 - n_pivots;
    n_pivots++;
  
    // look at this column (tmat's row) to find all the rows that hit it
    for(j=tmat->offset[col]; j < tmat->offset[col+1]; j++) {
      t_row = tmat->nonzero[j];
      assert((t_row) < mat->numrows);
      rowtrck[t_row] -= (1L<<32) + col;
      if( (rowtrck[t_row] >> 32) == 1L ) {
        queue[end++] = t_row;
      }
    }
  }
  
  t1 = wall_seconds() - t1;
  
  if(n_pivots != nr){
    printf("ERROR! toposort_matrix_queue: found %ld pivots but expected %ld!\n", n_pivots, nr);
    exit(1);
  }
  free(queue);
  free(rowtrck);
  return(t1);
}

/*!
 * \brief This routine implements the agi variant of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_loop(int64_t *rperm, int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) 
{
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  
  int64_t * rowtrck = calloc(nr, sizeof(int64_t));
  
  int64_t i, j, col, t_row;
   
  /* initialize rowtrck */
  for(i = 0; i < nr; i++){
    rowtrck[i] = 0L;
    for(j = mat->offset[i]; j < mat->offset[i+1]; j++)
      rowtrck[i] += (1L<<32) + mat->nonzero[j];
  }
  
  // we a pick a row with a single nonzero = col.
  // setting rperm[pos] = row and cprem[pos] = col
  // in a sense 'moves' the row and col to the bottom
  // right corner of matrix.
  // Next, we cross out that row and col by decrementing 
  //  the rowcnt for any row that contains that col
  // repeat
  //
  double t1 = wall_seconds();
  
  int64_t n_pivots = 0;
  while(n_pivots < nr){      
    for(i = 0; i < nr; i++){
			if( (rowtrck[i] >> 32) == 1 ){
			  col = rowtrck[i] & 0xFFFF;  // see cool trick
		   	rperm[i] = nr - 1 - n_pivots;
			  cperm[col] = nc - 1 - n_pivots;
        n_pivots++;
				rowtrck[i] = 0L;
		
			  // look at this column (tmat's row) to find all the rows that hit it
        for(j=tmat->offset[col]; j < tmat->offset[col+1]; j++) {
           t_row = tmat->nonzero[j];
           assert((t_row) < mat->numrows);
           rowtrck[t_row] -= (1L<<32) + col;
				}
			}
		}
	}
  
  t1 = wall_seconds() - t1;
  
  if(n_pivots != nr){
    printf("ERROR! toposort_matrix_queue: found %ld pivots but expected %ld!\n", n_pivots, nr);
    exit(1);
  }
  free(rowtrck);
  return(t1);
}

int main(int argc, char * argv[]) 
{
  #define NUMROWS 20000 
  int64_t numrows=NUMROWS, numcols;
  double erdos_renyi_prob = 0.1;
  uint32_t seed =  123456789;
  double laptime = 0.0;

  enum MODEL {GENERIC_Model=1, LOOP_Model=2, ALL_Models=4};
  uint32_t use_model;
  uint32_t models_mask = ALL_Models - 1;
  uint32_t printhelp = 0;
  uint32_t quiet = 0;
  uint32_t dump_files = 0;

  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:DM:q")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows );  break;
    case 's': sscanf(optarg,"%d" ,&seed );  break;
    case 'e': sscanf(optarg,"%lg", &erdos_renyi_prob); break;
    case 'D': dump_files = 1; break;
    case 'M': sscanf(optarg,"%d", &models_mask); break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }

  numcols = numrows;
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of toposort\n");
    fprintf(stderr,"help                 (-h)\n");
    fprintf(stderr,"number of rows       (-n)= %ld\n", numrows);
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"erdos_renyi_prob     (-e)= %g\n", erdos_renyi_prob);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files           (-D)\n");
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
    if(printhelp)
      return(0);
  }
  

  if(!quiet) printf("Creating input matrix for toposort\n");
  sparsemat_t * mat = generate_toposort_input (numrows, erdos_renyi_prob, seed, dump_files);
  if(!mat){printf("ERROR: topo: generate_toposort_input failed\n"); exit(1);}
  
  if(!quiet){
    printf("Input matrix stats:\n");
    spmat_stats(mat);
  }
  if(dump_files) dump_matrix(mat,20, "mat.out");

  write_matrix_mm(mat, "topo_mat.mm");

  sparsemat_t * tmat = transpose_matrix(mat);
  if(!tmat){printf("ERROR: topo: transpose_matrix failed\n"); exit(1);}

  if(dump_files) dump_matrix(tmat,20, "trans.out");
  write_matrix_mm(tmat, "topo_tmat.mm");

  if(!quiet) printf("Running toposort on mat (and tmat) ...\n");
  // arrays to hold the row and col permutations
  int64_t *rperminv2 = calloc(numrows, sizeof(int64_t));
  int64_t *cperminv2 = calloc(numcols, sizeof(int64_t));
  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC_Model:
      if(!quiet) printf("   using generic toposort: ");
      laptime = toposort_matrix_queue(rperminv2, cperminv2, mat, tmat);
      break;
    case LOOP_Model:
      if(!quiet) printf("   using loop    toposort: ");
      laptime = toposort_matrix_loop(rperminv2, cperminv2, mat, tmat);
      break;

    default:
       continue;
    }
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      fprintf(stderr,"nERROR: After toposort_matrix_queue: mat2 is not upper-triangular!\n");
      exit(1);
    }
    if(!quiet) printf("  %8.3lf seconds \n", laptime);
  }
  return(0);
}

