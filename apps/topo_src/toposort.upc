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
/*! \file toposort.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"

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
int check_result(sparsemat_t * mat, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t dump_files) {
  sparsemat_t * mat2;
  int ret = 0;

  int rf = is_perm(rperminv, mat->numrows);
  int cf = is_perm(cperminv, mat->numcols);
  if(!rf || !cf){
    T0_fprintf(stderr,"ERROR: check_result is_perm(rperminv2) = %d is_perm(cperminv2) = %d\n",rf,cf);
    return(1);
  }
  mat2 = permute_matrix_exstack(mat, rperminv, cperminv);
  if(!is_upper_triangular(mat2))
    ret = 1;
  if(dump_files) write_matrix(mat2, 0, "mat2.out");
  clear_matrix(mat2);
  free(mat2);
  return(ret);
}

/*! \brief Generate a matrix that is the a random permutation of a sparse uppper triangular matrix.
 * \param numrows the number of rows (and columns) in the produced matrix
 * \param nz_per_row the number of non-zeros in the matrix we generate before zero-ing out the lower triangle, 
 *        nz_per_row / numrows is the probability of non-diagonal entry being non-zero
 * \param rand_seed the seed for random number generator that determines the original matrix and the permutations
 * \param model a flag see spmat.h for which communication model to use
 * \param dump_files is a debugging flag
 * \return the permuted upper triangular matrix
 */
sparsemat_t * generate_toposort_input( int64_t numrows, int64_t nz_per_row,  int64_t rand_seed, int64_t model, int64_t dump_files)
{
  minavgmaxD_t mmtmp[1];
  double t1;

  T0_fprintf(stderr,"Creating input matrix for toposort\n");fflush(stderr);
  srand48( rand_seed );
  int64_t numcols = numrows;

  // get a square sparse matrix
  t1 = wall_seconds();
  sparsemat_t * omat = gen_uniform_sparse(numrows, nz_per_row);
  
  // remove entries below the diagonal and add 1's on the diagonal as necessary
  int64_t i, j, col, pivot, pos = 0, start = 0;
  for(i = 0; i < omat->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    if(start == omat->loffset[i+1]){
      printf("ERROR! toposort: empty row! %ld\n", i);
    }
    //printf("TH %d row %ld: ", MYTHREAD, global_row);fflush(0);
    for(j = start; j < omat->loffset[i+1]; j++){
      col = omat->lnonzero[j];
      //printf("row %ld col %ld\n", global_row, col);
      if(col >= global_row){
        if(!pivot){
          if(col != global_row){
            /* if the first entry is not on the diagonal, move it to the diagonal */
            col = global_row;
          }
          pivot = 1;
        }
        omat->lnonzero[pos++] = col;
      }
    }
    if(!pivot){
      /* if we are here, it is because all entries were below diagonal originally */
      /* so we should just shove the pivot in this row */
      //printf("adding a pivot in row %ld\n", global_row);
      omat->lnonzero[pos++] = global_row;
    }
    start = omat->loffset[i+1];
    omat->loffset[i+1] = pos;
  }

  omat->nnz = lgp_reduce_add_l(pos);
  omat->lnnz = pos;

  if(dump_files) write_matrix(omat, 0, "orig.out");
  if(!is_upper_triangular(omat)){
    exit(1);
  }
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( mmtmp, t1, THREADS );  
  T0_fprintf(stderr,"    generate upper triangular %8.3lf (with model %ld)\n", mmtmp->avg, model);

  // get row and column permutations

  t1 = wall_seconds();
  SHARED int64_t * rperminv = rand_permp(numrows, 1230+MYTHREAD, model);
  SHARED int64_t * cperminv = rand_permp(numcols, 45+MYTHREAD, model);
  lgp_barrier();
  if(!rperminv || !cperminv){
    T0_printf("ERROR: topo_rand_permp returns NULL!\n");fflush(0);
    return(NULL);
  }
  if(!is_perm(rperminv, numrows)){
    T0_printf("ERROR: topo_rand_permp rperminv is not a perm\n");fflush(0);
    return(NULL);
  }
  if(!is_perm(cperminv, numcols)){
    T0_printf("ERROR: topo_rand_permp cperminv is not a perm\n");fflush(0);
    return(NULL);
  }
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * lcperminv = lgp_local_part(int64_t, cperminv);
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( mmtmp, t1, THREADS );  
  T0_fprintf(stderr,"    compute row, column perms %8.3lf (with model %ld)\n", mmtmp->avg, model);


  if( dump_files && MYTHREAD == 0){
    FILE * fp = fopen("rperm.out", "w");
    for(i = 0; i < numrows; i++){
      fprintf(fp,"%ld\n", lgp_get_int64(rperminv, i));
    }
    fclose(fp);
    fp = fopen("cperm.out", "w");
    for(i = 0; i < numcols; i++){
      fprintf(fp, "%ld\n", lgp_get_int64(cperminv, i));
    }
    fclose(fp);
  }
  lgp_barrier();

  t1 = wall_seconds();
  sparsemat_t * mat = permute_matrix(omat, rperminv, cperminv, model);
  if(!mat) {
    T0_printf("ERROR: permute_matrix returned NULL");fflush(0);
    return(NULL);
  }
  lgp_barrier();
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( mmtmp, t1, THREADS );  
  T0_fprintf(stderr,"    apply the permutations    %8.3lf (with model %ld)\n", mmtmp->avg, model);

  if(dump_files) write_matrix(mat,0, "perm.out");

  clear_matrix( omat );
  free(omat);
  lgp_all_free(rperminv);
  lgp_all_free(cperminv);

  return( mat );
}

 /*! 
   To make things easy,  we start with an upper triangular matrix and randomly
   permute the rows and the columns.  The toposort algorithm takes this matrix
   and finds one of the possibly many row and column permutations 
   that would bring the matrix back to an upper triangular form.

   We set the row and column permutations,  rperm and cperm, one pivot at a time.

   for( pos=0; pos < number of rows; pos++ ) { 
     { 
       pick a row, r, with a single nonzero, c.
       then (r,c) is the pivot and set rperm[pos] = r and cprem[pos] = c.
       Note: from the above, we know such a row exists.
     } 
     {
       cross out that row r and col c
     }
   }

   Rather than changing the matrix by deleting rows and column and then searching the 
   new matrix for the next row.  We do the obvious thing of keeping row counts,
   rowcnt[i] is the number of non-zeros in row i and
   a really cool trick of keeping the sum of the column indices for the non-zeros
   in the row.
   Note: rowsum[i] is the sum of the column indices, not the non-zero elements,
   for the non-zeros in row i.
   To "delete a column" one decrements the rowcnt by one and the rowsum 
   by the corrsponding column index. 
   Hence, when the rowcnt gets to one, the rowsum is the column that is left.

   In parallel we have two race conditions to handle.

   Threads race to pick their pos in rperm and cperm. 
   We could race for the pivots with a fetch_and_add, 
   instead we use parallel prefix to claim enough room 
   for the pivots in the current level
   on each thread then assign them in order per thread
   
   Threads race to update the rowcnt and rowsum arrays. 
   We handle this with levels and atomic memory operations.

   The notion of a level is that we make a queue of all degree one rows,
   then we process all the rows on this queue and by doing so create
   new degree one rows. These rows are placed on the queue for the next
   level.  This solves the race between one thread processing rows on 
   its queue while other threads are adding rows to that queue.
*/


int main(int argc, char * argv[])
{
//X  typedef struct pkg_rowcnt_t{
//X    int64_t row;
//X    int64_t cnt;
//X  }pkg_rowcnt_t;

  int ret;
  int64_t i, j, fromth, lnnz, start, end;  
  int64_t pe, row, col, idx;
  double t1;
  
  lgp_init();
  /**************************************************************/
  /*                TOPOLOGICAL SORT EXAMPLE                    */
  /**************************************************************/
  /* We feed the toposort algorithm a random permutation of a   */
  /* random upper triangluar sparse matrix.                     */
  /**************************************************************/
//X  pkg_rowcnt_t pkg4;
  int64_t rowsperthread = 1000;
  int64_t nz_per_row = 35;
  int64_t stream_len = 1024;
  int64_t rand_seed =  MYTHREAD*MYTHREAD*MYTHREAD*10000 + 5;
  int64_t numrows, numcols;
  int64_t pos = 0;

  int64_t gen_model = EXSTACK_Model;
  int64_t solve_model = AGI_Model | EXSTACK_Model | EXSTACK2_Model | CONVEYOR_Model;
  int64_t dump_files = 0;
  minavgmaxD_t topo[1];
  setlocale(LC_NUMERIC,"");

  int opt; 
  while( (opt = getopt(argc, argv, "DN:Z:S:M:G:")) != -1 ) {
    switch(opt) {
    case 'N': sscanf(optarg,"%ld" ,&rowsperthread );  break;
    case 'Z': sscanf(optarg,"%ld" ,&nz_per_row    );  break;
    case 'S': sscanf(optarg,"%ld" ,&stream_len    );  break;
    case 'G': sscanf(optarg,"%ld" ,&gen_model);  break;
    case 'M': sscanf(optarg,"%ld" , &solve_model);  break;
    case 'D': dump_files = 1; break;
    default:  break;
    }
  }

  numrows = rowsperthread * THREADS;
  numcols = numrows;
  
  T0_fprintf(stderr,"Running toposort on %ld threads\n", THREADS);
  T0_fprintf(stderr,"Number of rows per thread   (-N)   %ld\n", rowsperthread);
  T0_fprintf(stderr,"Number of nonzeros per row  (-Z)   %ld\n", nz_per_row);
  T0_fprintf(stderr,"stream_len (stack size)     (-S)   %ld\n", stream_len);
  T0_fprintf(stderr,"task mask (M) = %ld (should be 1,2,4,8 for agi, exstack, exstack2, conveyors\n", solve_model);

  if( gen_model != EXSTACK_Model )
  {
     T0_fprintf(stderr,"generate input matrix with model (-G)   %ld\n", gen_model);
  }
  sparsemat_t * mat = generate_toposort_input(numrows, nz_per_row, rand_seed, gen_model, dump_files);
  if(!mat && !MYTHREAD){printf("ERROR: mat is NULL!\n"); exit(1);}

  if(dump_files) write_matrix(mat,0, "mat.out");

  sparsemat_t * tmat = transpose_matrix(mat, gen_model);
  if(!tmat && !MYTHREAD){printf("ERROR: tmat is NULL!\n"); exit(1);}

  if(dump_files) write_matrix(tmat,0, "trans.out");

  lgp_barrier();

  T0_fprintf(stderr,"Run toposort on mat (and tmat) ...\n");
  // arrays to hold the row and col permutations
  SHARED int64_t *rperminv2 = lgp_all_alloc(numrows, sizeof(int64_t));
  SHARED int64_t *cperminv2 = lgp_all_alloc(numcols, sizeof(int64_t));
  double t_atomic, t_exstack, t_exstack2, t_conveyor;
  double gb_th  = (mat->numrows + mat->numcols*2 + mat->nnz*2)*8;

  if( solve_model &  AGI_Model ) {                                         // Do the toposort with UPC and atomics
    t_atomic = toposort_matrix_upc(rperminv2, cperminv2, mat, tmat);
    //T0_fprintf(stderr,"      UPC: toposort matrix:     %8.3lf  %5.3lf GB/th/s\n", t_atomic, gb_th/(t_atomic*(1L<<30)));
    T0_fprintf(stderr,"      UPC: toposort matrix:     %8.3lf seconds\n", t_atomic);
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_upc: mat2 is not upper-triangular!\n");
    }
  }
  lgp_barrier();

  /*
  if( solve_model &  EXSTACK_Model ) {                                     // Do the toposort with exstack
    t_exstack = toposort_matrix_exstack(rperminv2, cperminv2, mat, tmat);
    T0_fprintf(stderr,"  Exstack: toposort matrix:     %8.3lf  %5.3lf GB/th/s\n", t_exstack, gb_th/(t_exstack*(1L<<30)));
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_exstack: mat2 is not upper-triangular!\n");
    }
  }
  lgp_barrier();
  */
  if( solve_model &  EXSTACK_Model ) {                                     // Do the toposort with exstack (no atomics)
    t_exstack = toposort_matrix_exstack_b(rperminv2, cperminv2, mat, tmat);
    //T0_fprintf(stderr,"  Exstack: toposort matrix:     %8.3lf  %5.3lf GB/th/s\n", t_exstack, gb_th/(t_exstack*(1L<<30)));
    T0_fprintf(stderr,"  Exstack: toposort matrix:     %8.3lf seconds\n", t_exstack);
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_exstack: mat2 is not upper-triangular!\n");
    }
  }
  lgp_barrier();

  if( solve_model &  EXSTACK2_Model ) {                                     // Do the toposort with exstack2
    t_exstack2 = toposort_matrix_exstack2(rperminv2, cperminv2, mat, tmat);
    //T0_fprintf(stderr,"  Exstack2 toposort matrix:     %8.3lf  %5.3lf GB/th/s\n", t_exstack2,  gb_th/(t_exstack2*(1L<<30)));
    T0_fprintf(stderr,"  Exstack2 toposort matrix:     %8.3lf seconds\n", t_exstack2);
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_exstack2: mat2 is not upper-triangular!\n");
    }
  }
  lgp_barrier();

  if( solve_model &  CONVEYOR_Model ) {                                     // Do the toposort with conveyors
    t_conveyor = toposort_matrix_convey(rperminv2, cperminv2, mat, tmat);
    //T0_fprintf(stderr,"  Conveyor toposort matrix:     %8.3lf  %5.3lf GB/th/s\n", t_conveyor,  gb_th/(t_conveyor*(1L<<30)));
    T0_fprintf(stderr,"  Conveyor toposort matrix:     %8.3lf seconds\n", t_conveyor);
    if( check_result(mat, rperminv2, cperminv2, dump_files) ) {
      printf("\nERROR: After toposort_matrix_conveyor: mat2 is not upper-triangular!\n");
    }
  }
  lgp_barrier();

  return(0);
  
}

