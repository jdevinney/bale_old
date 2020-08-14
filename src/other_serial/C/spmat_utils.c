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
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FTNESS
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

/*! \file spmat_utils.c
 * \brief Routines to implement basic spmat functions and some
 * utilities for debugging and timing.
 */

#include "spmat_utils.h"

/*! \page spmat_utils_page Sparse matrix routines and few utils */

/*! \brief  A routine to give access to the wall clock timer on most UNIX-like systems.
 * Uses gettimeofday.
 * \return the number of seconds since the beginning of time on UNIX machines.
 */

double wall_seconds() 
{
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1){ perror("gettimeofday"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );

}


/*! \brief create an int64_t array which holds a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation (return identity perm if seed is 0)
 * 
 * This implements the standard serial algorithm, known at least as
 * Fisher-Yates or Knuth shuffle, to generate the uniform permutation.
 * Start with an array holding the identity permutation,
 *   swap he last entry with a random entry in the array,
 *   this determines the last entry,
 *   shorten the array to be all but the last entry 
 *   and repeat
 */
int64_t * rand_perm(int64_t N, int64_t seed) 
{
  int64_t r, i, L, s;

  if( seed != 0 ) srand48( seed );

  int64_t * perm = malloc(N * sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  for(i=0; i<N; i++)
    perm[i] = i;
  if( seed == 0 )
    return(perm);
  for(L=N-1; L>0; L--){
    r = lrand48() % (L+1);  // must be allowed to swap the Lth slot with itself
    s = perm[L];
    perm[L] = perm[r];
    perm[r] = s;
  }
  return(perm);
}


/*! \brief checks that an array is in fact a permutation
 * \param perm pointer to the array
 * \param N the length of the array
 * \return 1 if it is permutation
 *
 * Every element in the flag array will equal 1 iff perm is a permutation.
 */
int64_t is_perm(int64_t * perm, int64_t N) 
{
  int64_t i, err=0L;
  int64_t * flag = calloc(N, sizeof(int64_t));
  if( flag == NULL ) return(0);

  for(i = 0; i < N; i++){
    if( 0 > perm[i] || perm[i] >= N )
      err++;
    else 
      flag[perm[i]]++;
  }
  for(i = 0; i < N; i++) 
    if(flag[i] != 1L) 
      err++;
  free(flag);
  return(err == 0L);
}

/*! \brief writes the first and last part of an array of int64_t's to the specified file
 * \param A pointer to the array
 * \param maxdisp the number of entries that are written, 0 means everything, 
     otherwise write the first and last maxdisplay/2 entries
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 */
int64_t dump_array(int64_t *A, int64_t len, int64_t maxdisp, char * name) 
{
  int64_t i, stoprow, startrow;
  stoprow = 0;
  startrow = 0;
  if((maxdisp > 0) && (maxdisp < len)){
    stoprow = maxdisp / 2;
    startrow = len - maxdisp/2;
  }

  FILE * fp = fopen(name, "w");
  for(i=0; i<len; i = (i<stoprow || i>=startrow)? i+1: startrow){
    fprintf(fp, "%3ld\n", A[i]);
  }
  fclose(fp);
  return(0);
}


/*! \brief dumps a sparse matrix to a file in a ASCII format
 * \param A pointer to the sparse matrix
 * \param maxrows the number of rows that are written, 0 means everything, 
     otherwise write the first and last maxrows/2 rows
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 */
int64_t dump_matrix(sparsemat_t *A, int64_t maxrows, char * name) 
{
  int64_t i,j, off, nxtoff;
  int64_t stoprow, startrow;
  stoprow = A->numrows;
  startrow = A->numrows;
  if( (maxrows > 0) && (maxrows < A->numrows) ){
    stoprow = maxrows / 2;
    startrow = A->numrows - maxrows/2;
  }

  FILE * fp = fopen(name, "w");
  fprintf(fp,"\n----- offsets:\n");
  for(i=0; i<stoprow; i++) 
    fprintf(fp, "%3ld ", A->offset[i]);
  if( stoprow < startrow)
    fprintf(fp, " ... ");
  for(i=startrow; i<A->numrows; i++) 
    fprintf(fp, "%3ld ", A->offset[i]);
  fprintf(fp, "%3ld ", A->offset[A->numrows]);

  fprintf(fp,"\n--------- nonzeros:\n");
  for(i=0; i<stoprow; i++){
    off    = A->offset[i];
    nxtoff = A->offset[i+1];
    for(j=off; j < nxtoff;j++){
      fprintf(fp, "%9ld %9ld",i, A->nonzero[j] );
      if(A->value)
        fprintf(fp," %9.5f\n",A->value[j]);
      else
        fprintf(fp,"\n");
    }
  }
  if( stoprow < startrow)
    fprintf(fp, ".\n.\n.\n");
  for(i=startrow; i<A->numrows; i++){
    off    = A->offset[i];
    nxtoff = A->offset[i+1];
    for(j=off; j < nxtoff;j++){
      fprintf(fp, "%9ld %9ld",i, A->nonzero[j] );
      if(A->value)
        fprintf(fp," %9.5f\n",A->value[j]);
      else
        fprintf(fp,"\n");
    }
  }
  fclose(fp);
  return(0);
}

/*! \brief writes a sparse matrix to a file in a MatrixMarket ASCII format
 * \param A pointer to the sparse matrix
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 */
int64_t write_matrix_mm(sparsemat_t *A, char * name){
  int64_t i,j;

  FILE * fp = fopen(name, "w");       // FIXME error
  if(A->value)
    fprintf(fp,"%%%%MasterMarket matrix coordinate real\n");
  else
    fprintf(fp,"%%%%MasterMarket matrix coordinate pattern\n");
  fprintf(fp, "%ld %ld %ld\n", A->numrows, A->numcols, A->nnz);

  for(i=0; i<A->numrows; i++){
    for(j=A->offset[i]; j < A->offset[i+1];j++)
      if(A->value)
        fprintf(fp, "%ld %ld %lf\n",i+1, A->nonzero[j]+1, A->value[j]);
      else
        fprintf(fp, "%ld %ld\n",i+1, A->nonzero[j]+1);
  }
  fclose(fp);
  return(0);
}


/*! \brief the compare function for qsort called while reading 
 * a MatrixMarket format.
 * NB. We sort on the rows so that we can fill the offset array
 * sequentially in one pass. We sort on the columns so that
 * the matrix will be "tidy"
 */
int elt_comp(const void *a, const void *b) 
{
  element_t * eltA = (element_t *)a;
  element_t * eltB = (element_t *)b;
  //return( eltA->row - eltB->row );
  if( (eltA->row - eltB->row) == 0 )
    return( eltA->col - eltB->col );
  return( eltA->row - eltB->row );
}

/*! \brief the compare function for qsort called while reading 
 * a MatrixMarket format (for matrix with values).
 * NB. We sort on the rows so that we can fill the offset array
 * sequentially in one pass. We sort on the columns so that
 * the matrix will be "tidy"
 */
int triple_comp(const void *a, const void *b) 
{
  triple_t * A = (triple_t *)a;
  triple_t * B = (triple_t *)b;
  if( (A->row - B->row) == 0 )
    return( A->col - B->col );
  return( A->row - B->row );
}

/*! \brief read a sparse matrix from a file in a MasterMarket ASCII format
 * \param name the filename to be read
 * \return a pointer to the sparse matrix or NULL on failure
 * TODO: get this to work for matrices with real values.
 */
sparsemat_t  *read_matrix_mm(char * name) 
{
  int64_t i;
  int64_t nr, nc, nnz;
  char object[24], format[24], field[24];
  int fscanfret;

  // Read the header line of the MatrixMarket format 
  FILE * fp = fopen(name, "r");
  if( fp == NULL ){
    fprintf(stderr,"read_matrix_mm: can't open file %s \n", name);
    exit(1);
  }
    
  fscanfret = fscanf(fp,"%%%%MatrixMarket %s %s %s\n", object, format, field);
  if( (fscanfret != 3 ) || strncmp(object,"matrix",24) || strncmp(format,"coordinate",24) ){
    fprintf(stderr,"read_matrix_mm: Incompatible matrix market format.\n");
    fprintf(stderr,"                First line should be either:\n");
    fprintf(stderr,"                matrix coordinate pattern\n");
    fprintf(stderr,"                OR\n");
    fprintf(stderr,"                matrix coordinate real\n");
    fprintf(stderr,"                OR\n");
    fprintf(stderr,"                matrix coordinate integer\n");
    exit(1);
  }
  
  if(strncmp(field,"pattern",24) && strncmp(field,"real",24) && strncmp(field,"integer",24) ){
    fprintf(stderr,"read_matrix_mm: Incompatible matrix market field.\n");
    fprintf(stderr,"                Last entry on first line should be pattern, real, or integer\n");
    exit(1);
  }
  int value;

  if(strncmp(field,"pattern",24) == 0){
    value = 0; // no values
  }else if(strncmp(field,"real",24) == 0){
    value = 1; // real values
  }else{
    //value = 2; // integer values
    printf("Error: read_matrix_mm: don't yet support integer weights (only doubles)\n");
    return(NULL);
  }

  fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64" %"SCNd64"\n", &nr, &nc, &nnz);
  if( (fscanfret != 3 ) || (nr<=0) || (nc<=0) || (nnz<=0) ){
    fprintf(stderr,"read_matrix_mm: reading nr, nc, nnz\n");
    exit(1);
  }
  
  int64_t pos = 0;
  int64_t row = 0;
  sparsemat_t * ret;
  if(!value){ // no values
    
    // read all the nonzeros into the elts array of (row,col)
    // and sort them before building the sparsemat format
    element_t * elts = calloc(nnz, sizeof(element_t));
    if( elts == NULL ){
      fprintf(stderr,"read_matrix_mm: elts calloc failed\n");
      exit(1);
    }

    for(i=0; i<nnz; i++){
      fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64"\n", &(elts[i].row), &(elts[i].col));
      assert (fscanfret == 2);
      //fprintf(stderr,"--- %ld %ld\n",  elts[i].row, elts[i].col);
      assert ( 0<elts[i].row && elts[i].row <=nr);
      assert ( 0<elts[i].col && elts[i].col <=nc);
      elts[i].row -=1;    // MasterMarket format is 1-up, not 0-up
      elts[i].col -=1;
    }
    
    qsort( elts, nnz, sizeof(element_t), elt_comp);

    ret = init_matrix( nr, nc, nnz, (value > 0));
    if( ret == NULL ){
      fprintf(stderr,"read_matrix_mm: sparsemat calloc failed\n");
      exit(1);
    }

    ret->offset[row] = 0;
    while( pos<nnz ){
      if( elts[pos].row == row ){
    ret->nonzero[pos] = elts[pos].col;
    pos++;
    continue;
      }
      while( row < elts[pos].row ){
    row++;
    ret->offset[row] = pos;
      }
    }
    ret->offset[row+1] = pos;

    free(elts);
    
  }else{ // real 

    triple_t * elts = calloc(nnz, sizeof(triple_t));
    if( elts == NULL ){
      fprintf(stderr,"read_matrix_mm: elts calloc failed\n");
      exit(1);
    }

    for(i=0; i<nnz; i++){
      if(value == 1)
        fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64" %lf\n", &(elts[i].row), &(elts[i].col), &(elts[i].val));
      assert (fscanfret == 3);
      assert ( 0<elts[i].row && elts[i].row <=nr);
      assert ( 0<elts[i].col && elts[i].col <=nc);
      elts[i].row -= 1;    // MasterMarket format is 1-up, not 0-up
      elts[i].col -= 1;
    }

    qsort( elts, nnz, sizeof(triple_t), triple_comp);

    ret = init_matrix( nr, nc, nnz, (value > 0));
    if( ret == NULL ){
      fprintf(stderr,"read_matrix_mm: sparsemat calloc failed\n");
      exit(1);
    }

    ret->offset[row] = 0;
    while( pos<nnz ){
      if( elts[pos].row == row ){
        ret->nonzero[pos] = elts[pos].col;
        ret->value[pos] = elts[pos].val;
        pos++;
        continue;
      }
      while( row < elts[pos].row ){
        row++;
        ret->offset[row] = pos;
      }
    }
    ret->offset[row+1] = pos;
    free(elts);
  }
  fclose(fp);
  return(ret);
}



/*! \brief apply row and column permutations to a sparse matrix 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced 
 */
sparsemat_t * permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv) 
{
  int64_t i, j, row, pos;

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, A->nnz, (A->value!=NULL));
  if(!Ap){printf("ERROR: permute_matrix: init_matrix failed!\n");return(NULL);}

  int64_t * rperm = calloc(A->numrows, sizeof(int64_t));
  if( rperm == NULL ) return(NULL);
  //compute rperm from rperminv 
  for(i=0; i < A->numrows; i++)
    rperm[rperminv[i]] = i;
  
  // fill in permuted rows with permuted columns
  Ap->offset[0] = pos = 0;
  for(i = 0; i < Ap->numrows; i++){
    row = rperm[i];
    for(j = A->offset[row]; j < A->offset[row+1]; j++){
      Ap->nonzero[pos] = cperminv[A->nonzero[j]];
      if(A->value) Ap->value[pos] = A->value[j];
      pos++;      
    }
    Ap->offset[i+1] = pos;
  }
  
  free(rperm);
  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix
 * Since this is serial, we don't have to worry about sorting the nonzero column
 * indicies (the get populated in the correct order).
 * 
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if error.
 */
sparsemat_t * transpose_matrix(sparsemat_t *A) 
{
  int64_t i, j, row, col;
  sparsemat_t * At;
  
  // get column counts
  int64_t * colcnt = calloc( A->numcols, sizeof(int64_t));
  if( colcnt == NULL ) return(NULL);
  int64_t * tmpoffset = colcnt;
  
  // histogram the column counts of A into colcnt
  for( i=0; i< A->nnz; i++){  
    assert( A->nonzero[i] < A->numcols );
    assert( A->nonzero[i] >= 0 ); 
    colcnt[A->nonzero[i]]++;
  }

  At = init_matrix(A->numcols, A->numrows, A->nnz, (A->value!=NULL));
  if(!At){printf("ERROR: transpose_matrix: init_matrix failed!\n");return(NULL);}

  //use the colcnt array to build the offset array and
  //reuse it as we fill in the nonzeros 
  At->offset[0] = 0;
  for(i = 0; i < At->numrows; i++){
    At->offset[i+1] = At->offset[i] + colcnt[i];
    tmpoffset[i] = At->offset[i];
  }

  //redistribute the nonzeros 
  for(row=0; row<A->numrows; row++){
    for(j=A->offset[row]; j<A->offset[row+1]; j++){
      col = A->nonzero[j];
      At->nonzero[ tmpoffset[col] ] = row;
      if(A->value)  // expect (or hope) that the compiler would pull this out of the loop
				At->value[ tmpoffset[col] ] = A->value[j];
      tmpoffset[col]++;
    }
  }

  free(colcnt);
  return(At);
}


/*! \brief checks that the sparse matrix is (strictly, i.e. zero on diagonal) lower triangluar
 * \param A pointer to the sparse matrix
 * \return 0 on success, non-0 on error.
 */
int64_t is_lower_triangular(sparsemat_t *A) 
{
  int64_t j, row;
  int64_t err = 0;

  for(row=0; row < A->numrows; row++){
    for(j=A->offset[row]; j < A->offset[row+1];j++){
      if( A->nonzero[j] > row ) 
        err++;
    }
  }
  return(err);
}

/*! \brief checks that the sparse matrix is (strictly) upper triangluar
 * \param A pointer to the sparse matrix
 * \return 0 on success, non-0 on error.
 */
int64_t is_upper_triangular(sparsemat_t *A) 
{
  int64_t j, row;
  int64_t err = 0;

  for(row=0; row < A->numrows; row++){
    for(j=A->offset[row]; j < A->offset[row+1];j++){
      if( A->nonzero[j] < row ) 
        err++;
    }
  }
  return(err);
}

/*! \brief comparison function for qsort in sort_nonzeros
 */
int nz_comp(const void *a, const void *b) 
{
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief comparison function for qsort in sort_nonzeros in a matrix with values
 */
int cv_comp(const void *a, const void *b)
 {
   return(((int64_t)((col_val_t*)a)->col) - ((int64_t)((col_val_t*)b)->col));
 }
 
/*! \brief sort the  non-zeros in each row of a sparse matrix so that their column indices
 * are in ascending order.
 * \param mat pointer to the sparse matrix
 */
int64_t sort_nonzeros( sparsemat_t *mat)
{
  int64_t i, j;
  if(mat->value){
    // we have to sort the column indicies, but we also have to permute the value array accordingly
    // this is annoying in C
    // we have to create an array of stucts that holds col,val pairs for a row
    // sort that array according to the col keys
    // and then overwrite the row data
    int64_t max = 0;
    for(i = 0; i < mat->numrows; i++)
      if(mat->offset[i+1] - mat->offset[i] > max) max = mat->offset[i+1] - mat->offset[i];

    // allocate a temporary array to hold a row's worth of col, value pairs
    col_val_t * tmparr = calloc(max, sizeof(col_val_t));

    for(i = 0; i < mat->numrows; i++){
      int64_t pos = 0;
      for(j = mat->offset[i]; j < mat->offset[i+1]; j++){
        tmparr[pos].col = mat->nonzero[j];
        tmparr[pos++].value = mat->value[j];
      }
      qsort(tmparr, mat->offset[i+1] - mat->offset[i], sizeof(col_val_t), cv_comp );
      pos = 0;
      for(j = mat->offset[i]; j < mat->offset[i+1]; j++){
        mat->nonzero[j] = tmparr[pos].col;
        mat->value[j] = tmparr[pos++].value;
      }
    }
    free(tmparr);
  }else{
    // no values, just sort the column indicies
    for(i = 0; i < mat->numrows; i++){
      qsort( &(mat->nonzero[mat->offset[i]]), mat->offset[i+1] - mat->offset[i], sizeof(int64_t), nz_comp );
    }
  }
  return(0);
}

/*! \brief compare the structs that hold two sparse matrices
 * \param lmat pointer to the left sparse matrix
 * \param rmat pointer to the right sparse matrix
 * \return 0 on success
 */
int64_t compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat) 
{
  int64_t i,j;
  int values;

  if( lmat->value == NULL && rmat->value == NULL){
    values = 0;
  }else if((lmat->value && !rmat->value) || (rmat->value && !lmat->value)){
    printf("Only one matrix has values!\n");
    return(1);    
  }else{
    values = 1;
  }
  
  if( lmat->numrows != rmat->numrows ){
    printf("(lmat->numrows = %ld)  != (rmat->numrows = %ld)", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->numrows != rmat->numrows ){
    printf("(lmat->numrows = %ld)  != (rmat->numrows = %ld)", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    printf("(lmat->numcols = %ld)  != (rmat->numcols = %ld)", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    printf("(lmat->nnz = %ld)  != (rmat->nnz = %ld)", lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    printf("(lmat->lnnz = %ld) != (rmat->lnnz = %ld)", lmat->nnz, rmat->nnz );
    return(1);
  }

  if( lmat->offset[0] != 0 || rmat->offset[0] != 0 || (lmat->offset[0] != rmat->offset[0] ) ){
    printf("(lmat->offset[0] = %ld)  != (rmat->offset[0] = %ld)", lmat->offset[0], rmat->offset[0] );
    return(1);
  }
  
  for(i = 0; i < lmat->numrows; i++){
    if( lmat->offset[i+1] != rmat->offset[i+1] ){
      printf("(lmat->offset[%ld] = %ld)  != (rmat->offset[%ld] = %ld)", i+1, lmat->offset[i+1], i+1, rmat->offset[i+1] );
      return(1);
    }
  }
  
  for(j=0; j< lmat->nnz; j++){
    if( lmat->nonzero[j] != rmat->nonzero[j] ){
      printf("(lmat->nonzero[%ld] = %ld)  != (rmat->nonzero[%ld] = %ld)", j,
             lmat->nonzero[j], j, rmat->nonzero[j] );      
      return(1);
    }
    if(values){
      if( lmat->value[j] != rmat->value[j] ){
        printf("(lmat->value[%ld] = %lf)  != (rmat->value[%ld] = %lf)", j,
               lmat->value[j], j, rmat->value[j] );
        return(1);
      }
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
  int64_t i,j;
  
  sparsemat_t * destmat = init_matrix(srcmat->numrows, srcmat->numcols, srcmat->nnz, (srcmat->value != NULL));
  if(!destmat) return(NULL);

  for(i = 0; i < (srcmat->numrows)+1; i++){
    destmat->offset[i] = srcmat->offset[i];
  }
  for(j=0; j < srcmat->nnz; j++){
    destmat->nonzero[j] = srcmat->nonzero[j];
    if(srcmat->value) destmat->value[j] = srcmat->value[j];
  }
  return(destmat);
}

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total number of rows
 * \param numcols total number of columns
 * \param nnz number of nonzero
 * \param values flag: 0 = no values, 1 = doubles
 * \return An initialized sparsemat_t (numrows, numcols, nnz are set) or NULL on error.
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz, int values) 
{
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->numrows  = numrows;
  mat->numcols  = numcols;
  mat->nnz  = nnz;  
  mat->offset   = calloc(mat->numrows+1, sizeof(int64_t));
  if(mat->offset == NULL) 
    return(NULL);
  mat->nonzero = calloc(mat->nnz, sizeof(int64_t));
  if(mat->nonzero == NULL) 
    return(NULL);
  if(values){
    mat->value = calloc(mat->nnz, sizeof(double));
    if(mat->value == NULL)
      return(NULL);
  }else
    mat->value = NULL;
  return(mat);
}

 /*************************************************************************************/
 /*                               RANDOM MATRICES                                     */
 /*************************************************************************************/

 /* 
  * What kinds of methods should we have to generate "random" sparse matrices?
  * 1) erdos-renyi - square upper or lower-triangular matrices 
  *    can be used to create symmetric or nonsymmetric square matrices
  * 2) uniform sparse: use ER technology
  * 3) kron_prod_graph (square symmetric or lower-triangular)
  * 4) random geometric? random planar graphs? Social network (preferential attachment, etc)
  * geometric graph - k-nearest neighbor graph, or epsilon ball graph
  *   - has different properties than ER
  *
  */
 // apps that need input matrices:
 // topo: square, upper-triangular, unit-diagonal, random or read in
 // permute_matrix and transpose: any random matrix, or read in matrix
 // triangle: kron_graph (special lower triangular), or any random lower triangular, or read in
 // SSSP: random non-symmetric square with values

 /*! \brief A routine to generate the adjacency matrix of a random graph.
  * If the graph is undirected, this routine only returns a lower-triangular
  * adjancency matrix (since the adjancency matrix would be symmetric and we don't need
  * the redundant entries).
  * 
  * \param n The number of vertices in the graph.
  * \param model FLAT: Erdos-Renyi random, GEOMETRIC: geometric random graph
  * \param edge_type See edge_type enum. Directed, or not, Weighted or not.
  * \param loops see self_loops enum. Does every node have a loop or not.
  * \param edge_density: d in [0, 1), target fraction of edges present.
  * \param seed: RNG seed.
  */
 sparsemat_t * random_graph(int64_t n, graph_model model, edge_type edge_type, self_loops loops,
                            double edge_density, int64_t seed){

   if(model == FLAT){
     
     return(erdos_renyi_random_graph(n, edge_density, edge_type, loops, seed));
     
   }else if(model == GEOMETRIC){
     double r;
     // determine the r that will lead to the desired edge density
     // Expected degree = n*pi*r^2
     // The expected number of edges E = n^2*pi*r^2/2
     // for undirected d = E/(n choose 2)
     // for directed   d = E/(n^2 - n)
     if (edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED){
       r = sqrt(edge_density/M_PI);
     }else{
       printf("ERROR: directed geometric graphs are not supported yet.\n");
       return(NULL);
     }

     return(geometric_random_graph(n, r, edge_type, loops, seed));
     
   }else{
     printf("ERROR: random_graph: Unknown type!\n");
     return(NULL);
   }
   
 }

 /*! \brief Subroutine to create a random sparse matrix.
  * 
  * The matrix is formed using a Erdo-Renyi like method. A[i,j] is nonzero with a fixed probability p.
  * The value of p depends on the density specified.
  *
  * \param nrows The number of rows
  * \param ncols The number of columns
  * \params density d in [0,1) specifying the fraction of entries that are nonzero
  * \params values 0 means all nonzeros are 1, 1 means all nonzeros are [0,1) uniform random.
  * \params seed RNG seed
  */
 
 sparsemat_t * random_sparse_matrix(int64_t nrows, int64_t ncols, double density, int values, int64_t seed){
   int64_t i, r;
   int64_t row, col;
   double lM = log(RAND_MAX);
   double D  = log(1 - density);
   int64_t nnz = 0;
   
   // first loop to count the number of nonzeros
   srand(seed);
   row = 0;
   do { r = rand(); } while(r == RAND_MAX);
   col = 1 + floor((log(RAND_MAX - r) - lM)/D);
   while(row < nrows){
     while(col < ncols){
       nnz++;
       do { r = rand(); } while(r == RAND_MAX);     
       col += 1 + floor((log(RAND_MAX - r) - lM)/D);
     }
     row++;
     col -= ncols;
   }

   sparsemat_t * A = init_matrix(nrows, ncols, nnz, values);

   //second pass: regenerate same random sequence and populate the matrix
   srand(seed);
   row = 0;
   nnz = 0;
   do { r = rand(); } while(r == RAND_MAX);
   col = 1 + floor((log(RAND_MAX - r) - lM)/D);
   while(row < nrows){
     while(col < ncols){
       A->nonzero[nnz++] = col;
       do { r = rand(); } while(r == RAND_MAX);     
       col += 1 + floor((log(RAND_MAX - r) - lM)/D);
     }
     row++;
     A->offset[row] = nnz;
     col -= ncols;
   }

   // fill in the random values
   if(values){
     for(i = 0; i < nnz; i++){
       A->value[i] = (double)rand()/RAND_MAX;
     }
   }
   return(A);
}

//****************** Kronecker Products of Star  Graphs *******************************/

/*! \brief Parses a list describing the input for the Kronecker Product construction.
 * \param list is the ascii string describing which stars are to be used in the construction.
 *   must be close to the form "M: k0 k1 k2 ...k(M-1)"
 *   where M is the mode:
 *    - mode == 0: default, no self loops
 *    - mode == 1: add a self loop to center vertex of star (row zero)
 *    - mode == 2: add a self loop to an outer vertex (the last one)
 *   then a space separated list of sizes.
 * \return a pointer to the argument struct
*/
kron_args_t * kron_args_init(char * list)
{
  int64_t cnt, error=0;
  int i, n;
  int64_t numrows = 1;

  kron_args_t * R = (kron_args_t *) calloc(1,sizeof(kron_args_t));

  strncpy(R->str, list, 256);
  do { 
    cnt =  sscanf(list, "%"PRId64": %n", &(R->mode), &n);
    if(cnt != 1){
      error = 1;     // didn't find the mode or its NULL
      break;
    }
    list += n;

    for (i=0; i<64; i++){
      cnt = sscanf(list, "%"PRId64" %n", &(R->star_size[i]), &n);
      if (cnt != 1)  // not an error, just done
        break;
      list += n;
      numrows *= (R->star_size[i] + 1); // a star is a K{1,m} graph
    }
    if(i < 1 || i > 63){
      error = 2;     // didn't find stars or found too many
      break;
    }
    R->num_stars = i;
    R->numrows = numrows;
  } while(0); // cause goto's are bad

  if( error ){
    printf("ERROR: trouble parsing the kron argument list! (error = %"PRId64")\n", error); 
    return(NULL);
  }
  return(R);
}

/*! \brief Generate the adjacency matrix for the star K_{1,m} (with or without loop edges)
 * \param m  the number of non-hub vertices
 * \param mode
 *  - mode == 0: default, no self loops
 *  - mode == 1: add a self loop to center vertex of star (row zero)
 *  - mode == 2: add a self loop to an outer vertex (the last one)
 * \return the adjacency matrix of the graph
*/
sparsemat_t * gen_star_graph(int64_t m, int mode) 
{
  sparsemat_t * G = init_matrix(m + 1, m + 1, 2*m + (mode > 0 ? 1 : 0), 0);
  
  int64_t i,j;
  int64_t pos = 0;  // counter thru the nonzeros

  G->offset[0] = 0;
  if(mode == 1)     // add a self loop to center vertex 
    G->nonzero[pos++] = 0;

  // Add the top row 
  for(j = 1; j < m+1; j++){
    G->nonzero[pos++] = j;
  }  
  G->offset[1] = pos;

  // rest of the rows 
  for(i = 1; i < m+1; i++){
    G->nonzero[pos++] = 0;
    G->offset[i+1] = pos;
  }
  
  if(mode == 2){  // add a self loop to last vertex and fix the last offset
    G->nonzero[pos++] = m;
    G->offset[m+1] = pos;
  }
  
  return(G);
}

/*! \brief Produce the Kronecker matrix product of two square {0,1}-matrices
 * \param A sparse matrix, A is A_n x A_n
 * \param B sparse matrix, B is B_n x B_n
 * \return $C =  A \bigotimes B$
 * 
 * Key fact is that element wise, using the notation C[row,col], A[r,s] and B[u,v],
 * 
 *  C[row,col] = C[r*B_n + u, s*B_n + v] = A[r,s] * B[u,v]
 *
*/
sparsemat_t * kronecker_mat_product(sparsemat_t * A, sparsemat_t * B) 
{
  int64_t r, s, u, v, j, k;
  int64_t A_n = A->numrows;
  int64_t B_n = B->numrows;
  int64_t C_n = A_n * B_n;
  int64_t C_nnz = A->nnz * B->nnz;

  sparsemat_t * C = init_matrix(C_n, C_n, C_nnz, 0);  // values = 0, for a {0,1}-matrix

  // get the number of nonzeros in each row 
  int64_t * C_rowtmp = calloc(C_n + 1, sizeof(int64_t));
  for(r = 0; r < A_n; r++){
    int64_t da = A->offset[r + 1] - A->offset[r];
    for(s = 0; s < B_n; s++){
      int64_t db = B->offset[s + 1] - B->offset[s]; 
      C_rowtmp[r*B_n + s] = da * db;
    }
  }

  // set the offsets in C 
  // we will reuse C_rowtmp[] to count from C->offset[r] to C->offset[r+1] 
  // as we fill in the nonzeros in row r of C
  C->offset[0] = 0;
  for(r = 0; r < C_n; r++){
    C->offset[r+1] = C->offset[r] + C_rowtmp[r];
    C_rowtmp[r] = C->offset[r];
  }
  

  // fill in the nonzeros of C, we fill them in a "block" order.
  // ie, we write all of B for each nonzero in A 
  for(r = 0; r < A_n; r++){
    for(j = A->offset[r]; j < A->offset[r+1]; j++){
      s = A->nonzero[j];
      for(u = 0; u < B_n; u++){
        for(k = B->offset[u]; k < B->offset[u+1]; k++){
          v = B->nonzero[k];
          C->nonzero[C_rowtmp[r*B_n + u]++] = s * B_n + v;
        }
      }
    }
  }
  
  free(C_rowtmp);

  return(C);
}


/*! \brief Generate the kroncker product of a collection of star graphs.
 * \param K is a struct that holds all the parameters for the construction
 * \param mode 
 *  - mode == 0: default, no self loops
 *  - mode == 1: add a self loop to center vertex of each star 
 *  - mode == 2: add a self loop to an outer vertex (the last vertex) of each star
 * \return the adjacency matrix for graph
*/
sparsemat_t * kronecker_product_graph(kron_args_t * K){
  int64_t i,j;
  int64_t G_n, G_nnz;
  sparsemat_t * star;
  int64_t M = K->num_stars;


  if(M < 2){
    printf("ERROR: kronecker_product_graph requires  number of stars to be >= 2!\n");
    return(NULL);
  }

  sparsemat_t ** mats = calloc(2*M, sizeof(sparsemat_t *));

  mats[0] = gen_star_graph(K->star_size[0], K->mode);

  for(i = 1; i < M; i++){
    star = gen_star_graph(K->star_size[i], K->mode);
    mats[i] = kronecker_mat_product(mats[i-1] , star);
    clear_matrix( star );
  }
  sparsemat_t * R = mats[M-1];  // what we want but is full symmetric and may have loops
  
  // count the number of nonzeros in the lower triangle part of the matrix
  G_n = R->numrows;
  G_nnz = 0;
  for(i=0; i<R->numrows; i++){
    for(j=R->offset[i]; (j<R->offset[i+1]) && (R->nonzero[j] < i ); j++){
      G_nnz++;
    }
  }

  // copy the lower triangle of the matrix (excluding the diagonal) to the returned matrix
  sparsemat_t * G = init_matrix(G_n, G_n, G_nnz, 0);
  G->offset[0] = 0;
  int64_t pos = 0;
  for(i=0; i<R->numrows; i++){
    for(j=R->offset[i]; (j<R->offset[i+1]) && (R->nonzero[j] < i ); j++){
      G->nonzero[pos++] = R->nonzero[j];
    }
    G->offset[i+1] = pos;
  }
  // free mats
  return(G);
}

int64_t tri_count_kron_graph(kron_args_t * K)
{
   double approx, ns;
   int i;

   if( K->mode == 0 ) 
     return 0;
   else if( K->mode == 1 ){
     approx = 1.0;
     for(i = 0; i < K->num_stars; i++)
        approx *= (3*K->star_size[i] + 1);
     return ( round((approx / 6.0) - 0.5 * K->numrows + 1.0/3.0) );
   } else if( K->mode == 2 ){
     ns = (double) K->num_stars;
     return( round( (1.0/6.0)*pow(4,ns) - pow(2.0,ns - 1.0) + 1.0/3.0) );
   } else {
     printf("ERROR: : init_matrix failed!\n");
     return(-1); 
   }
}
 

//****************************** Geometric Graphs *********************************************/
 typedef struct points_t{
   double x;
   double y;
   int64_t index;
 }points_t;
 
 typedef struct sector_t{
   points_t * points;
   int64_t numpoints;
 }sector_t;

// returns the square of the L2 distance
double dist(points_t a, points_t b){
  return((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}


 
 /*! \brief Generates the adjacency matrix for a random geometric graph. 
  * 
  * See https://en.wikipedia.org/wiki/Random_geometric_graph
  * Each vertex corresponds to a point randomly placed in the unit square. Two vertices
  * are adjancent if their corresponding points are within distance r of each other.
  * 
  * \param n The number of vertices
  * \param r The size of the neighborhoods that determine edges.
  * \param type See edge_type. If undirected, this routine returns a lower triangular matrix.
  * \param loops See self_loops. Are there self loops?
  * \param seed A seed for the RNG.
  * \return An adjacency matrix (or lower portion of in the undirected case).
  */
sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed){
  int64_t i, j, k, l;
  double r2 = r*r;
  // We break up the unit square into chunks that are rxr grid.
  // We generate the points uniformly at random over the unit square.
  // We calculate edges by comparing distances between each point and every other point in
  // its own sector and in neighboring sectors.
  
  int64_t nsectors = ceil(1.0/r);
  printf("GEOMETRIC with r = %lf number of sectors = %"PRId64"\n", r, nsectors);
  sector_t ** sectors = calloc(nsectors, sizeof(sector_t*));
  int64_t ** first_index_this_sector = calloc(nsectors, sizeof(int64_t*));
  for(i = 0; i < nsectors; i++){
    sectors[i] = calloc(nsectors, sizeof(sector_t));
    first_index_this_sector[i] = calloc(nsectors, sizeof(int64_t));
    for(j = 0; j < nsectors; j++)
      sectors[i][j].numpoints = 0;
  }


  // First pass, generate the points and count how many will fall in each sector.
  srand(seed);
  for(i = 0; i < n; i++){
    double x = (double)rand()/RAND_MAX;
    double y = (double)rand()/RAND_MAX;
    int64_t row = floor(y/r);
    int64_t col = floor(x/r);
    assert(row < nsectors);
    assert(col < nsectors);
    sectors[row][col].numpoints++;
  }

  // initialize the struct to hold the points
  for(i = 0; i < nsectors; i++){
    for(j = 0; j < nsectors; j++){
      if(j > 0)
        first_index_this_sector[i][j] += first_index_this_sector[i][j-1] + sectors[i][j-1].numpoints;
      else if(i > 0)
        first_index_this_sector[i][j] += first_index_this_sector[i-1][nsectors-1] + sectors[i-1][nsectors-1].numpoints;
      sectors[i][j].points = calloc(sectors[i][j].numpoints, sizeof(points_t));
    }
  }

  // reset numpoints
  for(i = 0; i < nsectors; i++){
    for(j = 0; j < nsectors; j++){
      sectors[i][j].numpoints = 0;
    }
  }
  // Second pass: generate the points and insert them into the struct
  srand(seed);
  for(i = 0; i < n; i++){
    double x = (double)rand()/RAND_MAX;
    double y = (double)rand()/RAND_MAX;
    int64_t row = floor(y/r);
    int64_t col = floor(x/r);
    int64_t li = sectors[row][col].numpoints;
    sectors[row][col].points[li].x = x;
    sectors[row][col].points[li].y = y;
    sectors[row][col].points[li].index = first_index_this_sector[row][col] + li;
    sectors[row][col].numpoints++;
  }
  
  
  // next we will count number of edges
  int64_t node = 0;
  int64_t nedges = 0;
  for(i = 0; i < nsectors; i++){
    for(j = 0; j < nsectors; j++){
      
      sector_t * sec = &sectors[i][j];
      int64_t m = sec->numpoints;
      for(k = 0; k < m; k++){

        // count the edges to lower-indexed nodes within this sector
        for(l = 0; l < k; l++){
          if(dist(sec->points[k], sec->points[l]) < r2)
            nedges++;
        }

        // count the edges to lower-indexed nodes outside the sector
        // to do this, we need to look at sectors to the W, NW, N, and NE.
        // W
        if(j > 0){
          sector_t * sec2 = &sectors[i][j-1];
          for(l = 0; l < sec2->numpoints; l++){
            if(dist(sec->points[k], sec2->points[l]) < r2)
              nedges++;
          }
        } 
        // NW
        if(i > 0 && j > 0){
          sector_t * sec2 = &sectors[i-1][j-1];
          for(l = 0; l < sec2->numpoints; l++){
            if(dist(sec->points[k], sec2->points[l]) < r2)
              nedges++;
          }
        } 
        // N
        if(i > 0){
          sector_t * sec2 = &sectors[i-1][j];
          for(l = 0; l < sec2->numpoints; l++){
            if(dist(sec->points[k], sec2->points[l]) < r2)
              nedges++;
          }
        }
        // NE
        if(i > 0 && j < (nsectors-1)){
          sector_t * sec2 = &sectors[i-1][j+1];
          for(l = 0; l < sec2->numpoints; l++){
            if(dist(sec->points[k], sec2->points[l]) < r2)
              nedges++;
          }
        } 
      }
    }
  }
  
  if(loops == LOOPS)
    nedges += n;
  //printf("nedges = %"PRId64" %lf\n", nedges, nedges/((double)n*(n-1)/2.0));
   
  int weighted = (edge_type == UNDIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, nedges, weighted);

  // go back through the loop and populate the adjacency matrix
  nedges = 0;
  for(i = 0; i < nsectors; i++){
    for(j = 0; j < nsectors; j++){
       
       sector_t * sec = &sectors[i][j];
       int64_t m = sec->numpoints;
       for(k = 0; k < m; k++){
         
         node = sec->points[k].index;

         if(loops == LOOPS){
           A->nonzero[nedges] = node;
           if(weighted) A->value[nedges] = rand()/RAND_MAX;
           nedges++;
         }
         
         // count the edges to lower-indexed nodes within this sector
         for(l = 0; l < k; l++){
           if(dist(sec->points[k], sec->points[l]) < r2){
             A->nonzero[nedges] = sec->points[l].index;
             assert(node >= sec->points[l].index);
             if(weighted) A->value[nedges] = rand()/RAND_MAX;
             nedges++;
           }
         }

         // count the edges to lower-indexed nodes outside the sector
         // to do this, we need to look at sectors to the W, NW, N, and NE.
         // W
         if(j > 0){
           sector_t * sec2 = &sectors[i][j-1];
           for(l = 0; l < sec2->numpoints; l++){
             if(dist(sec->points[k], sec2->points[l]) < r2){
               A->nonzero[nedges] = sec2->points[l].index;
               assert(node >= sec2->points[l].index);
               if(weighted) A->value[nedges] = rand()/RAND_MAX;
               nedges++;
             }
           }
         } 
         // NW
         if(i > 0 && j > 0){
           sector_t * sec2 = &sectors[i-1][j-1];
           for(l = 0; l < sec2->numpoints; l++){
             if(dist(sec->points[k], sec2->points[l]) < r2){
               A->nonzero[nedges] = sec2->points[l].index;
               assert(node >= sec2->points[l].index);
               if(weighted) A->value[nedges] = rand()/RAND_MAX;
               nedges++;
             }
           }
         } 
         // N
         if(i > 0){
           sector_t * sec2 = &sectors[i-1][j];
           for(l = 0; l < sec2->numpoints; l++){
             if(dist(sec->points[k], sec2->points[l]) < r2){
               A->nonzero[nedges] = sec2->points[l].index;
               assert(node >= sec2->points[l].index);
               if(weighted) A->value[nedges] = rand()/RAND_MAX;
               nedges++;
             }
           }
         }
         // NE
         if(i > 0 && j < (nsectors-1)){
           sector_t * sec2 = &sectors[i-1][j+1];
           for(l = 0; l < sec2->numpoints; l++){
             if(dist(sec->points[k], sec2->points[l]) < r2){
               A->nonzero[nedges] = sec2->points[l].index;
               assert(node >= sec2->points[l].index);
               if(weighted) A->value[nedges] = rand()/RAND_MAX;
               nedges++;
             }
           }
         }
         A->offset[node+1] = nedges;
         node++;
       }
              
     }
   }
   sort_nonzeros(A);
   return(A);
 }

 
/*! \brief Generates the lower half of the adjacency matrix for an Erdos-Renyi random graph. 
 * This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
 * this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the 
 * adjancecy matrix using a geometric random variable.
 *
 * \param n The total number of vertices in the graph (rows in the matrix).
 * \param p The probability that each non-loop edge is present.
 * \param edge_type See edge_type. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops.
 * \param seed A random seed.
 * \return A sparsemat_t
 */
sparsemat_t * erdos_renyi_random_graph(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed)
{
  int64_t row, col;
  double lM = log(RAND_MAX);
  double D  = log(1 - p);
  int64_t nnz;
  int64_t r;
  int64_t end = n;
  int64_t ndiag = n;

  assert(n > 0);

  // first loop to count the number of nonzeros
  srand(seed);
  row = 0;
  nnz = 0;
  do { r = rand(); } while(r == RAND_MAX);     
  col = 1 + floor((log(RAND_MAX - r) - lM)/D);
  while(row < n){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) // we hit a diagonal for free!
        ndiag--;
      nnz++;
      do { r = rand(); } while(r == RAND_MAX);     
      col += 1 + floor((log(RAND_MAX - r) - lM)/D);
    }
    row++;
    col -= end;
  }
  if(loops == LOOPS) nnz += ndiag;

  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * mat = init_matrix(n, n, nnz, weighted);

  if(!mat){ printf("ERROR: erdos_renyi_random_graph: init_matrix failed!\n"); return(NULL); }

  // fill in the nonzeros
  srand(seed);
  nnz = 0;
  row = 0;
  mat->offset[0] = 0;
  do { r = rand(); } while(r == RAND_MAX);     
  col = 1 + floor((log(RAND_MAX - r) - lM)/D);
  while(row < n){
    int need_diag = (loops==LOOPS);
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) need_diag = 0;
      mat->nonzero[nnz++] = col;
      do { r = rand(); } while(r == RAND_MAX);     
      col += 1 + floor((log(RAND_MAX - r) - lM)/D);     
    }
    
    if(need_diag){
      mat->nonzero[nnz++] = row;
    }
    row++;
    mat->offset[row] = nnz;
    col -= end;
  }

  // fill in weights
  if(weighted){
    int64_t i;
    for(i = 0; i < nnz; i++){
      mat->value[i] = (double)rand()/RAND_MAX;
    }
  }
  assert(mat->nnz == nnz);
  if(loops == LOOPS && (edge_type == DIRECTED || edge_type == DIRECTED_WEIGHTED))
    sort_nonzeros(mat); // to get the diagonal entry sorted correctly
  return(mat);
}

/*! \brief Generates the upper or lower half of the adjacency matrix for an Erdos-Renyi random graph. 
 * This is here mostly for illustration because it is so slow.
 * It flips a coin for each possible edge.
 * \param n The total number of vertices in the graph (rows in the matrix).
 * \param p The probability that each non-loop edge is present.
 * \param edge_type See edge_type.
 * \param loops See self_loops.
 * \param seed A random seed.
 * \return A sparsemat_t
 */

sparsemat_t * erdos_renyi_random_graph_naive(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed){
  int64_t row, col;
  int64_t P = p*RAND_MAX;
  int64_t pos;
  
  assert(n > 0);
  
  srand(seed);
  int64_t nnz = 0;
  int64_t end = n;
  for(row = 0; row < n; row++){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;    
    for(col = 0; col < end; col++){
      if(col == row) continue;
      if( rand() < P ){
        nnz++;
      }
    }
  }
  if(loops == LOOPS) nnz += n;
  
  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * mat = init_matrix(n, n, nnz, weighted);
  if(!mat){ printf("ERROR: erdos_renyi_random_graph_naive: init_matrix failed!\n"); return(NULL); }

  // fill in the nonzeros
  srand(seed);
  pos = 0;
  mat->offset[0] = 0;
  for(row = 0; row < n; row++){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row + (loops == LOOPS);
    for(col = 0; col < end; col++){
      if(col == row && loops == LOOPS){
        mat->nonzero[pos++] = row;
        continue;
      }
      if( rand() < P ){
        mat->nonzero[pos++] = col;
      }
    }
    mat->offset[row+1] = pos;
  }
  // fill in the weights
  if(weighted){
    int64_t i;
    for(i = 0; i < pos; i++){
      mat->value[i] = (double)rand()/RAND_MAX;
    }
  }
  
  return( mat);
}



/*! \brief prints some stats of a sparse matrix
 * \param mat the sparse matrix
 */
void spmat_stats(sparsemat_t *mat) 
{
  printf("   mat->numrows  = %12"PRId64"\n", mat->numrows);
  printf("   mat->numcols  = %12"PRId64"\n", mat->numcols);
  printf("   mat->nnz      = %12"PRId64"\n", mat->nnz);
  
  int64_t i, d, mindeg, maxdeg, cntdeg;
  double avgdeg;
  mindeg = mat->numcols;
  maxdeg = 0;
  cntdeg = 0;
  for(i=0; i<mat->numrows; i++){
    d = mat->offset[i+1]-mat->offset[i];
    cntdeg += d;
    mindeg = (d < mindeg) ? d : mindeg;
    maxdeg = (d > maxdeg) ? d : maxdeg;
  }
  avgdeg = (double)cntdeg/(double)mat->numrows;
  
  printf("  min, avg, max degree = %"PRId64", %g, %"PRId64"\n", mindeg, avgdeg, maxdeg);
}


/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 */
void clear_matrix(sparsemat_t * mat)
{
  free(mat->nonzero);
  free(mat->offset);
  free(mat->value);
}
