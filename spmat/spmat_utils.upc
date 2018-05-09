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
/*! \file spmat_utils.upc
 * \brief Utilities to support spmat and wrappers to make it easier to switch 
 * between routines written in the various models.
 */
#include <spmat.h>
#include <assert.h>

/*! \brief Produce a global array the holds a uniform random permutation.
 * \param N the global length of the permutaion
 * \param seed the seed for the random number generator
 * \param mode a flag that sets which model (atomic, exstack, exstack2, conveyor) is used
 * \param mode which buffering model to use
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 *
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp(int64_t N, int seed, int64_t mode)
{
  if( mode == AGI_Model ) 
    return( rand_permp_atomic(N, seed) );
  else if( mode == EXSTACK_Model )
    return( rand_permp_exstack(N, seed) );
  else if( mode == EXSTACK2_Model )
    return( rand_permp_exstack2(N, seed) );
  else if( mode == CONVEYOR_Model )
    return( rand_permp_exstack2(N, seed) );
  return(NULL);
}


/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param omat pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \param mode which buffering model to use 
 * \return a pointer to the matrix that has be computed or NULL on failure
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix(sparsemat_t *omat, SHARED int64_t *rperminv, SHARED int64_t *cperminv, int64_t mode)
{
  if(mode == AGI_Model) 
    return( permute_matrix_atomic(omat, rperminv, cperminv) );
  else if(mode == EXSTACK_Model)
    return( permute_matrix_exstack(omat, rperminv, cperminv) );
  else if(mode == EXSTACK2_Model)
    return( permute_matrix_exstack2(omat, rperminv, cperminv) );
  else if(mode == CONVEYOR_Model)
    return( permute_matrix_exstack2(omat, rperminv, cperminv) );
  return(NULL);
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param omat  pointer to the original matrix
 * \param mode which buffering model to use
 * \return a pointer to the matrix that has be computed or NULL on failure
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix(sparsemat_t *omat, int64_t mode )
{
  if( mode == AGI_Model )
    return( transpose_matrix_atomic(omat) );
  else if( mode == EXSTACK_Model )
    return( transpose_matrix_exstack(omat) );
  else if( mode == EXSTACK2_Model )
    return( transpose_matrix_exstack2(omat) );
  else if( mode == CONVEYOR_Model )
    return( transpose_matrix_conveyor(omat) );
  return(NULL);
}



/*! \brief writes a sparse matrix to a file in a ASCII format
 * \param A pointer to the sparse matrix
 * \param maxrows the number of rows that are written, 0 means everything, 
         otherwise write the first and last maxrows/2 rows
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 * \ingroup spmatgrp
 */
int write_matrix(sparsemat_t *A, int maxrows, char * name)
{
  int64_t i,j, off, nxtoff;
  int64_t stoprow, startrow;
  stoprow = 0;
  startrow = 0;
  if( maxrows > 0 ){
     stoprow = maxrows / 2;
     startrow = A->numrows - maxrows/2;
  }

  lgp_barrier();
  
  if(MYTHREAD == 0){
    FILE * fp = fopen(name, "w");
    fprintf(fp,"\n----- offsets:\n");
    for(i=0; i < stoprow; i++) {
      if( i%THREADS == 0 ) 
        fprintf(fp, "\n");
      fprintf(fp, "%3ld ", lgp_get_int64(A->offset,i));
    }
    if( maxrows > 0 ) fprintf(fp,"\n...\n");
    for(i=startrow; i < A->numrows + 1; i++) {
      if( i%THREADS == 0 ) 
        fprintf(fp, "\n");
      fprintf(fp, "%3ld ", lgp_get_int64(A->offset,i));
    }
    fprintf(fp,"\n--------- nonzeros:\n");
    for(i=0; i < stoprow; i++){
      off    = lgp_get_int64(A->offset, i);
      nxtoff = lgp_get_int64(A->offset, i + THREADS);
      for(j=off; j < nxtoff;j++){
        fprintf(fp, "%9ld %9ld\n",i, lgp_get_int64(A->nonzero, j*THREADS + i%THREADS) );
      }
    }
    if( maxrows > 0 ) fprintf(fp,"\n...\n");
    for(i=startrow; i < A->numrows; i++){
      off    = lgp_get_int64(A->offset, i);
      nxtoff = lgp_get_int64(A->offset, i + THREADS);
      for(j=off; j < nxtoff;j++){
        fprintf(fp, "%9ld %9ld\n",i, lgp_get_int64(A->nonzero, j*THREADS + i%THREADS) );
      }
    }
    fclose(fp);
  }
  lgp_barrier();
  return(0);
}

/*! \brief checks that the sparse matrix is lower triangluar
 * \param A pointer to the sparse matrix
 * \return 0 on success, non-0 on error.
 * kind of a speciality routine to check that toposort might of worked
 * \ingroup spmatgrp
 */
int is_lower_triangular(sparsemat_t *A)
{
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;

  lgp_barrier();

  /* we are hoping for an lower triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j=A->loffset[i]; j < A->loffset[i+1];j++){
      if( A->lnonzero[j] > global_row ) {
        err++;
      }else if( A->lnonzero[j] == global_row){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = lgp_reduce_add_l(err);
  err2 = lgp_reduce_add_l(err2);
  if( err || err2 ){
    if(!MYTHREAD)printf("\nERROR!!! There are %ld nz above diag. and %ld missing pivots in lower.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief checks that the sparse matrix is upper triangluar
 * \param A pointer to the sparse matrix
 * \return 0 on success, non-0 on error.
 * kind of a speciality routine to check that toposort might of worked
 * \ingroup spmatgrp
 */
int is_upper_triangular(sparsemat_t *A)
{
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;
  lgp_barrier();

  /* we are hoping for an upper triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j=A->loffset[i]; j < A->loffset[i+1];j++){
      if( A->lnonzero[j] < global_row ) {
        err++;
      }else if( A->lnonzero[j] == global_row){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = lgp_reduce_add_l(err);
  err2 = lgp_reduce_add_l(err2);
  if( err || err2 ){
    if(!MYTHREAD)printf("\nERROR!!! There are %ld nz below diag. and %ld missing pivots in upper.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief checks that a global array is in fact a permutation
 * \param perm SHARED pointer to the global array
 * \param N the length of the global array
 * \return 1 if it is permutation
 * \ingroup spmatgrp
 */
int is_perm(SHARED int64_t * perm, int64_t N)
{
  int64_t i;
  int64_t * lperm = lgp_local_part(int64_t, perm);
  SHARED int64_t * flag = lgp_all_alloc(N, sizeof(int64_t));
  if( flag == NULL ) return(0);
  int64_t * lflag = lgp_local_part(int64_t, flag);
  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;

  for(i = 0; i < l_N; i++) 
    lflag[i] = 0;
  lgp_barrier();
  for(i = 0; i < l_N; i++) 
   lgp_put_int64(flag, lperm[i], 1);
  lgp_barrier();
  int64_t err = 0L;
  for(i = 0; i < l_N; i++) 
    if(lflag[i] == 0) 
      err++;
  lgp_all_free(flag);
  err = lgp_reduce_add_l(err);
  return(err == 0L);
}

/*! \brief comparison function to support qsort 
 */
int nz_comp(const void *a, const void *b)
{
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief create a random square sparse matrix with N rows and columns 
 * \param N total number of rows and columns
 *  note the total number of rows is a multiple of THREADS
 * \param max_row_cnt max number of columns in a row
 * \return 0 on success, non-0 on error.
 * The definition of uniform here is:
 *  the row counts are uniformly distributed between 1 and max_row_cnt+1
 *  each of the columns are chosen uniformly, with replacement, from 0 to numcols-1
 * Rows can have repeated columns, columns appear in non-decreasing order
 * \ingroup spmatgrp
 */
sparsemat_t * gen_uniform_sparse(int64_t N, int64_t max_row_cnt)
{
  int i,j;
  int64_t lnnz, numrows, numcols, lnumrows, lnumcols;

  numrows = N;
  numcols = numrows;
  lnumrows = (numrows + THREADS - MYTHREAD - 1)/THREADS;
  lnumcols = (numcols + THREADS - MYTHREAD - 1)/THREADS;
  int64_t * loffset = calloc(lnumrows + 1, sizeof(int64_t));

  int64_t rrc;
  lnnz = loffset[0] = 0;
  for(i = 0; i < lnumrows; i++){
     rrc = lrand48() % max_row_cnt;
     rrc += 1;
     loffset[i+1] = loffset[i] + rrc;
     lnnz += rrc;
  }
  
  sparsemat_t * mat = init_matrix(numrows, numcols, lnnz);
  if(!mat){
    printf("ERROR: gen_uniform_sparse: init_matrix failed!\n");
    return(NULL);
  }

  int64_t pos = 0;
  int64_t lastval, start, col;
  mat->loffset[0] = 0;
  for(i = 0; i < lnumrows; i++){

    int64_t cnt = loffset[i+1] - loffset[i];
    
    for(j=0; j < cnt; j++) {
      col = lrand48() % numcols;
      mat->lnonzero[pos + j] = col;
    }
    
    qsort( &(mat->lnonzero[pos]), cnt, sizeof(int64_t), nz_comp );
    
    start = pos;
    lastval = mat->lnonzero[pos];
    mat->lnonzero[pos++] = lastval;
    for(j = 1; j < cnt; j++) {
      if(mat->lnonzero[start + j] != lastval){
        lastval = mat->lnonzero[start + j];
        mat->lnonzero[pos++] = mat->lnonzero[start + j];        
      }
    }
    mat->loffset[i+1] = pos;
  }

  mat->lnnz = pos;
  mat->nnz = lgp_reduce_add_l(pos);

  lgp_barrier();

  free(loffset);
  
  return(mat);
}

/*! \brief sort the non-zeros in each row of a sparse matrix
 * \param mat pointer to the sparse matrix
 */
int sort_nonzeros( sparsemat_t *mat)
{
  int i;
  for(i = 0; i < mat->lnumrows; i++){
     qsort( &(mat->lnonzero[mat->loffset[i]]), mat->loffset[i+1] - mat->loffset[i], sizeof(int64_t), nz_comp );
  }
  lgp_barrier();
  return(0);
}

/*! \brief compare the structs that hold two sparse matrices
 * \param lmat pointer to the left sparse matrix
 * \param rmat pointer to the right sparse matrix
 * \return 0 on success
 * \ingroup spmatgrp
 */
int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat)
{
  int i,j;

  if( lmat->numrows != rmat->numrows ){
    if(!MYTHREAD)printf("(lmat->numrows = %ld)  != (rmat->numrows = %ld)", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->lnumrows != rmat->lnumrows ){
    if(!MYTHREAD)printf("THREAD %03d: (lmat->lnumrows = %ld)  != (rmat->lnumrows = %ld)", 
              MYTHREAD, lmat->lnumrows, rmat->lnumrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    if(!MYTHREAD)printf("(lmat->numcols = %ld)  != (rmat->numcols = %ld)", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    if(!MYTHREAD)printf("(lmat->nnz = %ld)  != (rmat->nnz = %ld)", lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->lnnz != rmat->lnnz ){
    if(!MYTHREAD)printf("THREAD %03d: (lmat->lnnz = %ld)  != (rmat->lnnz = %ld)", 
       MYTHREAD, lmat->lnnz, rmat->lnnz );
    return(1);
  }

  if( lmat->loffset[0] != 0 || rmat->loffset[0] != 0 
    || (lmat->loffset[0] != rmat->loffset[0] ) ){
    if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[0] = %ld)  != (rmat->loffset[0] = %ld)", 
       MYTHREAD, lmat->loffset[0], rmat->loffset[0] );
    return(1);
  }

  
  for(i = 0; i < lmat->lnumrows; i++){
    if( lmat->loffset[i+1] != rmat->loffset[i+1] ){
       if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[%d] = %ld)  != (rmat->loffset[%d] = %ld)", 
          MYTHREAD, i+1, lmat->loffset[i+1], i+1, rmat->loffset[i+1] );
       return(1);
    }
  }
  
  for(j=0; j< lmat->lnnz; j++) {
    if( lmat->lnonzero[j] != rmat->lnonzero[j] ){
      if(!MYTHREAD)printf("THREAD %03d: (lmat->lnonzero[%d] = %ld)  != (rmat->lnonzero[%d] = %ld)", 
                MYTHREAD, j, lmat->lnonzero[j], j, rmat->lnonzero[j] );
      return(1);
    }
  }

  return(0);
}

/*! \brief makes an exact copy of a given sparse matrices
 * \param srcmat pointer to the original sparse matrix
 * \return A pointer to the cloned sparse matrix
 */
sparsemat_t * copy_matrix(sparsemat_t *srcmat)
{
  int i,j;
  int64_t numrows, numcols, lnumrows, lnumcols;
  
  sparsemat_t * destmat = init_matrix(srcmat->numrows, srcmat->numcols, srcmat->lnnz);
  if(!destmat) return(NULL);

  for(i = 0; i < (srcmat->lnumrows)+1; i++){
     destmat->loffset[i] = srcmat->loffset[i];
  }
  for(j=0; j < srcmat->lnnz; j++) {
     destmat->lnonzero[j] = srcmat->lnonzero[j];
  }

  lgp_barrier();
  return(destmat);
}

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total number of rows
 * \param numcols total number of columns
 * \param nnz_this_thread number of nonzero on this thread
 * \return An initialized sparsemat_t or NULL on error.
 * \ingroup spmatgrp
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread){
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->numrows  = numrows;
  mat->lnumrows = (numrows + THREADS - MYTHREAD - 1)/THREADS;
  mat->numcols  = numcols;  
  mat->offset   = lgp_all_alloc(mat->numrows + THREADS, sizeof(int64_t));
  if(mat->offset == NULL) 
    return(NULL);
  mat->loffset  =  lgp_local_part(int64_t, mat->offset);
  int64_t max = lgp_reduce_max_l(nnz_this_thread);
  int64_t total = lgp_reduce_add_l(nnz_this_thread);
  mat->nonzero = lgp_all_alloc(max*THREADS, sizeof(int64_t));
  if(mat->nonzero == NULL) 
    return(NULL);
  mat->lnonzero = lgp_local_part(int64_t, mat->nonzero);
  mat->nnz = total;
  mat->lnnz = nnz_this_thread;

  mat->l_enum_row = -1;  // initialize the row iterator stuff
  mat->l_enum_idx = -1;
  mat->l_enum_nstop = -1;

  mat->S_enum_row = -1;
  mat->S_enum_idx = -1;
  mat->S_enum_nstop = -1;

  return(mat);
}


/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
void clear_matrix(sparsemat_t * mat){
  lgp_all_free(mat->nonzero);
  lgp_all_free(mat->offset);
}
