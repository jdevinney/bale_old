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
/*! \file spmat_agi.upc
 * \brief Sparse matrix support functions implemented with global addresses and atomics
 */
#include <spmat.h>
#include <assert.h>

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_atomic(int64_t N, int seed)
{  
  int64_t * ltarget, *lperm;
  int64_t r, i, j;
  int64_t pos, numdarts, numtargets, lnumtargets;

  if( seed != 0 ) srand48( seed );

  //T0_printf("Entering rand_permp_atomic...");fflush(0);

  SHARED int64_t * perm = lgp_all_alloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  lperm = lgp_local_part(int64_t, perm);

  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = 2*N;
  int64_t l_M = (M + THREADS - MYTHREAD - 1)/THREADS;

  SHARED int64_t * target = lgp_all_alloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  ltarget = lgp_local_part(int64_t, target);
  
  for(i=0; i<l_M; i++)
    ltarget[i] = -1L;
  lgp_barrier();

  i=0;
  while(i < l_N){                // throw the darts until you get l_N hits
    r = lrand48() % M;
    if( lgp_cmp_and_swap(target, r, -1L, (i*THREADS + MYTHREAD)) == (-1L) ){
      i++;
    }
  }
  lgp_barrier();

  numdarts = 0;
  for(i = 0; i < l_M; i++)    // count how many darts I got
    numdarts += (ltarget[i] != -1L );

  pos = lgp_prior_add_l(numdarts);    // my first index in the perm array is the number 
                                      // of elements produce by the smaller threads
  for(i = 0; i < l_M; i++){
    if(ltarget[i] != -1L ) {
       lgp_put_int64(perm, pos, ltarget[i]);
       pos++;
    }
  }

  lgp_all_free(target);
  lgp_barrier();
  //T0_printf("done!\n");
  return(perm);
}


/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_atomic(sparsemat_t *A, SHARED int64_t *rperminv, SHARED int64_t *cperminv)
{
  //T0_printf("Permuting matrix with single puts\n");
  int64_t i, j, col, row, pos;
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  SHARED int64_t * rperm = lgp_all_alloc(A->numrows, sizeof(int64_t));
  if( rperm == NULL ) return(NULL);
  int64_t *lrperm = lgp_local_part(int64_t, rperm);

  //compute rperm from rperminv 
  for(i=0; i < A->lnumrows; i++){
    lgp_put_int64(rperm, lrperminv[i], i*THREADS + MYTHREAD);
  }

  lgp_barrier();
  
  int64_t cnt = 0, off, nxtoff;
  for(i = 0; i < A->lnumrows; i++){
    row = lrperm[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    cnt += nxtoff - off;
  }
  lgp_barrier();

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, cnt);
  
  // fill in permuted rows
  Ap->loffset[0] = pos = 0;
  for(i = 0; i < Ap->lnumrows; i++){
    row = lrperm[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    for(j = off; j < nxtoff; j++){
      Ap->lnonzero[pos++] = lgp_get_int64(A->nonzero, j*THREADS + row%THREADS);
    }
    Ap->loffset[i+1] = pos;
  }
  
  assert(pos == cnt);
  
  lgp_barrier();
  
  // finally permute column indices
  for(i = 0; i < Ap->lnumrows; i++){
    for(j = Ap->loffset[i]; j < Ap->loffset[i+1]; j++){
      Ap->lnonzero[j] = lgp_get_int64(cperminv, Ap->lnonzero[j]);      
    }
  }
  lgp_barrier();

  lgp_all_free(rperm);

  //T0_printf("done\n");
  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_atomic(sparsemat_t *A)
{
  int64_t counted_nnz_At;
  int64_t lnnz, i, j, col, row, fromth, idx;
  int64_t pos;
  sparsemat_t * At;
  
  //T0_printf("UPC version of matrix transpose...");
  
  // find the number of nnz.s per thread

  SHARED int64_t * shtmp = lgp_all_alloc( A->numcols + THREADS, sizeof(int64_t));
  if( shtmp == NULL ) return(NULL);
  int64_t * l_shtmp = lgp_local_part(int64_t, shtmp);
  int64_t lnc = (A->numcols + THREADS - MYTHREAD - 1)/THREADS;
  for(i=0; i < lnc; i++)
    l_shtmp[i] = 0;
  lgp_barrier();

  for( i=0; i< A->lnnz; i++) {                   // histogram the column counts of A
    assert( A->lnonzero[i] < A->numcols );
    assert( A->lnonzero[i] >= 0 ); 
    pos = lgp_fetch_and_inc(shtmp, A->lnonzero[i]);
  }
  lgp_barrier();


  lnnz = 0;
  for( i = 0; i < lnc; i++) {
    lnnz += l_shtmp[i];
  }
  
  At = init_matrix(A->numcols, A->numrows, lnnz);
  if(!At){printf("ERROR: transpose_matrix_upc: init_matrix failed!\n");return(NULL);}

  int64_t sum = lgp_reduce_add_l(lnnz);      // check the histogram counted everything
  assert( A->nnz == sum ); 

  // compute the local offsets
  At->loffset[0] = 0;
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + l_shtmp[i-1];

  // get the global indices of the start of each row of At
  for(i = 0; i < At->lnumrows; i++)
    l_shtmp[i] = MYTHREAD + THREADS * (At->loffset[i]);
    
  lgp_barrier();

  //redistribute the nonzeros 
  for(row=0; row<A->lnumrows; row++) {
    for(j=A->loffset[row]; j<A->loffset[row+1]; j++){
      pos = lgp_fetch_and_add(shtmp, A->lnonzero[j], (int64_t) THREADS);
      lgp_put_int64(At->nonzero, pos, row*THREADS + MYTHREAD);
    }
  }
  lgp_barrier();


  lgp_barrier();
  //if(!MYTHREAD)printf("done\n");
  lgp_all_free(shtmp);

  return(At);
}


