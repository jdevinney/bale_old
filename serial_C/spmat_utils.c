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
  if (retVal == -1) { perror("gettimeofday"); fflush(stderr); }
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

  for(i = 0; i < N; i++) {
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
  for(i=0; i<len; i = (i<stoprow || i>=startrow)? i+1: startrow) {
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
    for(j=off; j < nxtoff;j++)
      fprintf(fp, "%9ld %9ld\n",i, A->nonzero[j] );
  }
  if( stoprow < startrow)
    fprintf(fp, ".\n.\n.\n");
  for(i=startrow; i<A->numrows; i++){
    off    = A->offset[i];
    nxtoff = A->offset[i+1];
    for(j=off; j < nxtoff;j++)
      fprintf(fp, "%9ld %9ld\n",i, A->nonzero[j] );
  }

  fclose(fp);
  return(0);
}

/*! \brief writes a sparse matrix to a file in a MatrixMarket ASCII format
 * \param A pointer to the sparse matrix
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 */
int64_t write_matrix_mm(sparsemat_t *A, char * name) 
{
  int64_t i,j;

  FILE * fp = fopen(name, "w");       // FIXME error
  fprintf(fp,"%%%%MasterMarket matrix coordinate position\n");
  fprintf(fp, "%ld %ld %ld\n", A->numrows, A->numcols, A->nnz);

  for(i=0; i<A->numrows; i++){
    for(j=A->offset[i]; j < A->offset[i+1];j++)
      fprintf(fp, "%ld %ld\n",i+1, A->nonzero[j]+1);
  }
  fclose(fp);
  return(0);
}


/*! \struct element_t 
 * \brief structure used while reading and writing the MatrixMarket format.
 * We only handle {0,1} matrices, so we don't need triples.
 */
typedef struct element_t { 
  int64_t row; //!< row
  int64_t col; //!< col
} element_t;

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

/*! \brief read a sparse matrix from a file in a MasterMarket ASCII format
 * \param name the filename to be read
 * \return a pointer to the sparse matrix or NULL on failure
 */
sparsemat_t  *read_matrix_mm(char * name) 
{
  int64_t i;
  int64_t nr, nc, nnz;
  char object[24], format[24], field[24];
  int fscanfret;

  // Read the header line of the MasterMarket format 
  FILE * fp = fopen(name, "r");
  if( fp == NULL ) {
    fprintf(stderr,"read_matrix_mm: can't open file %s \n", name);
    exit(1);
  }
  //The only format we need to be able to read requires the first line to be:
  //"%%MasterMarket matrix coordinate position"
  fscanfret = fscanf(fp,"%%%%MasterMarket %s %s %s\n", object, format, field);
  if( (fscanfret != 3 ) || strncmp(object,"matrix",24)|| strncmp(format,"coordinate",24) || strncmp(field,"position",24)){
    fprintf(stderr,"read_matrix_mm: can't read format of file %s\n", name);
    exit(1);
  }

  fscanfret = fscanf(fp,"%ld %ld %ld\n", &nr, &nc, &nnz);
  if( (fscanfret != 3 ) || (nr<=0) || (nc<=0) || (nnz<=0) ) {
    fprintf(stderr,"read_matrix_mm: reading nr, nc, nnz\n");
    exit(1);
  }

  // read all the nonzeros into the elts array of (row,col)
  // and sort them before building the sparsemat format
  element_t *elts = calloc(nnz, sizeof(element_t));
  if( elts == NULL ) {
    fprintf(stderr,"read_matrix_mm: elts calloc failed\n");
    exit(1);
  }
  for(i=0; i<nnz; i++) {
    fscanfret = fscanf(fp,"%ld %ld\n", &(elts[i].row), &(elts[i].col));
    assert (fscanfret == 2);
    //fprintf(stderr,"--- %ld %ld\n",  elts[i].row, elts[i].col);
    assert ( 0<elts[i].row && elts[i].row <=nr);
    assert ( 0<elts[i].col && elts[i].col <=nc);
    elts[i].row -=1;    // MasterMarket format is 1-up, not 0-up
    elts[i].col -=1;
  }
  qsort( elts, nnz, sizeof(element_t), elt_comp);

  sparsemat_t *ret = init_matrix( nr, nc, nnz );
  if( ret == NULL ) {
    fprintf(stderr,"read_matrix_mm: sparsemat calloc failed\n");
    exit(1);
  }

  int64_t pos = 0;
  int64_t row = 0;
  ret->offset[row] = 0;
  while( pos<nnz ){
    if( elts[pos].row == row ) {
      ret->nonzero[pos] = elts[pos].col;
      pos++;
      continue;
    }
    while( row < elts[pos].row ) {
      row++;
      ret->offset[row] = pos;
    }
  }
  ret->offset[row+1] = pos;

  free(elts);
  fclose(fp);
  return(ret);
}



/*! \brief apply row and column permutations to a sparse matrix 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv) 
{
  int64_t i, j, row, pos;

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, A->nnz);
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
    for(j = A->offset[row]; j < A->offset[row+1]; j++)
      Ap->nonzero[pos++] = cperminv[A->nonzero[j]];
    Ap->offset[i+1] = pos;
  }
  
  free(rperm);
  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * transpose_matrix(sparsemat_t *A) 
{
  int64_t i, j, row;
  sparsemat_t * At;
  
  // get column counts
  int64_t * tmp = calloc( A->numcols, sizeof(int64_t)); 
  if( tmp == NULL ) return(NULL);
  
  // use tmp to hold the column counts
  // histogram the column counts of A into tmp
  for( i=0; i< A->nnz; i++) {  
    assert( A->nonzero[i] < A->numcols );
    assert( A->nonzero[i] >= 0 ); 
    tmp[A->nonzero[i]]++;
  }

  At = init_matrix(A->numcols, A->numrows, A->nnz);
  if(!At){printf("ERROR: transpose_matrix: init_matrix failed!\n");return(NULL);}

  //use the tmp array to build the offset array and
  //reuse it as we fill in the nonzeros 
  At->offset[0] = 0;
  for(i = 0; i < At->numrows; i++){
    At->offset[i+1] = At->offset[i] + tmp[i];
    tmp[i] = At->offset[i];
  }
  tmp[i] = At->offset[i];

  //redistribute the nonzeros 
  //This is still tidy, right?
  for(row=0; row<A->numrows; row++) {
    for(j=A->offset[row]; j<A->offset[row+1]; j++){
      At->nonzero[ tmp[A->nonzero[j]] ] = row;
      tmp[A->nonzero[j]]++;
    }
  }

  free(tmp);
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

/*! \brief sort the non-zeros in each row of a sparse matrix
 * \param mat pointer to the sparse matrix
 */
int64_t sort_nonzeros( sparsemat_t *mat) 
{
  int64_t i;
  for(i = 0; i < mat->numrows; i++){
    qsort( &(mat->nonzero[mat->offset[i]]), mat->offset[i+1] - mat->offset[i], sizeof(int64_t), nz_comp );
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
  
  for(j=0; j< lmat->nnz; j++) {
    if( lmat->nonzero[j] != rmat->nonzero[j] ){
      printf("(lmat->nonzero[%ld] = %ld)  != (rmat->nonzero[%ld] = %ld)", j, lmat->nonzero[j], j, rmat->nonzero[j] );
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
  int64_t i,j;
  
  sparsemat_t * destmat = init_matrix(srcmat->numrows, srcmat->numcols, srcmat->nnz);
  if(!destmat) return(NULL);

  for(i = 0; i < (srcmat->numrows)+1; i++){
    destmat->offset[i] = srcmat->offset[i];
  }
  for(j=0; j < srcmat->nnz; j++) {
    destmat->nonzero[j] = srcmat->nonzero[j];
  }
  return(destmat);
}

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total number of rows
 * \param numcols total number of columns
 * \param nnz number of nonzero
 * \return An initialized sparsemat_t (numrows, numcols, nnz are set) or NULL on error.
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz) 
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
  return(mat);
}


// TODO: this whole part for generating random matrices needs to be reworked.
// - It really can only make square matrices (cause we only cared about graphs)
// - There is a cute trick to get around the fact that you have
//    to generate the matrix before you can allocate memory, then
//    you have to generate the same matrix to fill it in.
//    (Hence the compile time warning).
//
/*! \brief Generates the upper or lower half of the adjacency matrix for an Erdos-Renyi random graph. 
 * This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
 * this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the 
 * adjancecy matrix using a geometric random variable.
 *
 * \param numrows The total number of vertices in the graph (rows in the matrix).
 * \param p The probability that each non-loop edge is present.
 * \param mode  enum ER_TRIANGLE:
 *              ER_TRI_L, ER_TRI_U, ER_TRI_LWD, ER_TRI_UWD for 
 *              lower triangle, upper triangle, lower with diagonal, upper with diagonal, resp.
 * \param seed A random seed.
 * \return A sparsemat_t
 */
sparsemat_t * erdos_renyi_tri(int numrows, double p, enum ER_TRIANGLE mode, int64_t seed) 
{
  int64_t row, col;
  double lM = log(RAND_MAX);
  double D  = log(1 - p);

  int64_t numcols = numrows;
  int64_t nnz;
  int64_t r;
  int64_t lower, diag; 
  sparsemat_t *mat;

  assert(numrows > 0);
  switch(mode){ 
  case ER_TRI_L   : lower =1; diag = 0; break;
  case ER_TRI_U   : lower =0; diag = 0; break;
  case ER_TRI_LWD : lower =1; diag = 1; break;
  case ER_TRI_UWD : lower =0; diag = 1; break;
  default: 
    fprintf(stderr," ERROR: erdos_renyi_tri unknown mode\n"); return(NULL); 
    break; 
  }
  
  // run this loop twice. the first time to count, then alloc.
  // the second time to store.
  int keep;
  for(keep=0; keep<2; keep++) { 
    srand(seed);
    if(keep) mat->offset[0] = 0;
    nnz = 0;
    row = 0;
    col = (lower) ? 0 : 1;
    if(diag && (lower != 1)) {
      if(keep) mat->nonzero[nnz] = 0;
      nnz++;
    }
    while(row < numrows){
      do { r = rand(); } while(r == RAND_MAX);
  
      col += 1 + floor((log(RAND_MAX - r) - lM)/D);
  
      while((col >= ((lower) ? row : numrows)) && (row < numrows)){
        if(lower){ 
          if(diag) {
            if(keep) mat->nonzero[nnz] = row;
            nnz++;
          }
          col = col - row;
          row ++;
          if(keep) mat->offset[row] = nnz;
        }else{
          row ++;
          col = row + 1 + col - numrows;
          if(keep) mat->offset[row] = nnz;
          if((row < numrows) && diag) {
            if(keep)  mat->nonzero[nnz] = row;
            nnz++;
          }
        }      
      }
      if(row < numrows) {
        if(keep) mat->nonzero[nnz] = col;
        nnz++;
      }
    }
    if(keep) mat->offset[row] = nnz;
    if(keep == 0){
      mat = init_matrix(numrows, numcols, nnz);
      if(!mat){ printf("ERROR: generate_toposort_input: init_matrix failed!\n"); return(NULL); }
    }
  }
  //dump_matrix( mat, 0, "erdos_renyi");
  return( mat);
}

/*! \brief Generates the upper or lower half of the adjacency matrix for an Erdos-Renyi random graph. 
 * This is here mostly for illustration cause it is so slow.
 * It flips a coin for each possible edge
 * \param numrows The total number of vertices in the graph (rows in the matrix).
 * \param p The probability that each non-loop edge is present.
 * \param mode  enum ER_TRIANGLE:
 *              ER_TRI_L, ER_TRI_U, ER_TRI_LWD, ER_TRI_UWD for 
 *              lower triangle, upper triangle, lower with diagonal, upper with diagonal, resp.
 * \param seed A random seed.
 * \return A sparsemat_t
 */

sparsemat_t * naive_erdos_renyi_tri(int numrows, double p, enum ER_TRIANGLE mode, int64_t seed) 
{
  int64_t row, col;
  int64_t P = p*RAND_MAX;
  int64_t numcols = numrows;
  int64_t nnz=0;
  int64_t pos;
  int64_t lower, diag;
  
  assert(numrows > 0);
  switch(mode){ 
  case ER_TRI_L : lower =1; diag = 0; break;
  case ER_TRI_U : lower =0; diag = 0; break;
  case ER_TRI_LWD : lower =1; diag = 1; break;
  case ER_TRI_UWD : lower =0; diag = 1; break;
  default: 
    fprintf(stderr," ERROR: naive_erdos_renyi_tri unknown mode\n"); return(NULL); 
    break; 
  }
  
  srand(seed);
  nnz = 0;
  for(row = 0; row < numrows; row++){
    if( lower ) {
      for(col = 0; col < row; col++){
        if( rand() < P )
          nnz++;
      }
    } else {
      for(col = row+1; col < numcols; col++){
        if( rand() < P )
          nnz++;
      }
    }
    if( diag ) 
      nnz++;
  }
  
  sparsemat_t * mat = init_matrix(numrows, numcols, nnz);
  if(!mat){ printf("ERROR: generate_toposort_input: init_matrix failed!\n"); return(NULL); }
  srand(seed);
  // fill in the nonzeros, including the diagonal entry
  pos = 0;
  mat->offset[0] = 0;
  for(row = 0; row < numrows; row++){
    if( lower ) {
      for(col = 0; col < row; col++){
        if( rand() < P )
          mat->nonzero[pos++] = col;
      }
      if( diag ) 
        mat->nonzero[pos++] = row;
      mat->offset[row+1] = pos;
    } else {
      if( diag ) 
        mat->nonzero[pos++] = row;
      for(col = row+1; col < numcols; col++){
        if( rand() < P )
          mat->nonzero[pos++] = col;
      }
      mat->offset[row+1] = pos;
    }
  }
  return( mat);
}

/*! \brief Generates an Erdos-Renyi random graph
 * Puts together appropriate triangluar matrices created by erdos_renyi_tri
 * \param numrows The total number of vertices in the graph (rows in the matrix).
 * \param p The probability that each non-loop edge is present.
 * \param mode  enum ER_GRAPH:
 *              ER_GRAPH_SIMPLE, ER_GRAPH_DIRECT for a simple (undirected) graph, simple directed graph, resp
 * \param seed A random seed.
 * \return A sparsemat_t
 */

sparsemat_t * erdos_renyi_graph(int64_t numrows, double p, enum ER_GRAPH mode, int64_t seed) 
{
  int64_t row, j;
  int64_t pos;
  sparsemat_t *mat, *U, *L;
  
  assert(numrows > 0);
  switch(mode){ 
  case ER_GRAPH_SIMPLE:
    U = erdos_renyi_tri(numrows, p, ER_TRI_U, seed); 
    L = transpose_matrix(U);
    break;
  case ER_GRAPH_DIRECT: 
    U = erdos_renyi_tri(numrows, p, ER_TRI_U, seed); 
    L = erdos_renyi_tri(numrows, p, ER_TRI_L, seed); 
    break;
  default: 
    fprintf(stderr," ERROR: erdos_renyi_graph unknown mode\n"); return(NULL); 
    break; 
  }
  
  mat = init_matrix(numrows, numrows, U->nnz + L->nnz );
  if(!mat){ printf("ERROR: gen_erdos_renyi_graph: init_matrix failed!\n"); return(NULL); }
  
  pos=0;
  mat->offset[0] = 0;
  for(row=0; row<numrows; row++) {
    for(j=L->offset[row]; j< L->offset[row+1]; j++)
      mat->nonzero[pos++] = L->nonzero[j];     
    for(j=U->offset[row]; j< U->offset[row+1]; j++)
      mat->nonzero[pos++] = U->nonzero[j];     
    mat->offset[row+1] = pos;
  }
  
  return( mat );
}


/*! \brief prints some stats of a sparse matrix
 * \param mat the sparse matrix
 */
void spmat_stats(sparsemat_t *mat) 
{
  printf("   mat->numrows  = %12ld\n", mat->numrows);
  printf("   mat->numcols  = %12ld\n", mat->numcols);
  printf("   mat->nnz      = %12ld\n", mat->nnz);
  
  int64_t i, d, mindeg, maxdeg, cntdeg;
  double avgdeg;
  mindeg = mat->numcols;
  maxdeg = 0;
  cntdeg = 0;
  for(i=0; i<mat->numrows; i++) {
    d = mat->offset[i+1]-mat->offset[i];
    cntdeg += d;
    mindeg = (d < mindeg) ? d : mindeg;
    maxdeg = (d > maxdeg) ? d : maxdeg;
  }
  avgdeg = (double)cntdeg/(double)mat->numrows;
  
  printf("  min, avg, max degree = %ld, %g, %ld\n", mindeg, avgdeg, maxdeg);
}


/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 */
void clear_matrix(sparsemat_t * mat)
{
  free(mat->nonzero);
  free(mat->offset);
}
