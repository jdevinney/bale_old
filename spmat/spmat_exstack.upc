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
#include <spmat.h>
#include <exstack.h>
#include <assert.h>

/*! \file spmat_exstack.upc
 * \brief spmat routines written using exstack
 */ 

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * This implements the random dart algorithm to generate the permutation using exstack.
 * Each thread throws its elements of the perm array randomly at large target array by
 * pushing a throw on to the appropriate exstack.
 * The popping pe then claims the requested slot or returns the dart to the sender.
 * 
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 *  
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_exstack(int64_t N, int seed)
{
  int ret;
  int64_t i, j, cnt, pe, pos, fromth, istart, iend;
  int64_t val;
  
  int64_t lN = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = N * 2;
  int64_t lM = (M + THREADS - MYTHREAD - 1)/THREADS;

  typedef struct pkg_t{
    int64_t idx; 
    int64_t val;
  }pkg_t;

  SHARED int64_t * perm = lgp_all_alloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  int64_t * lperm = lgp_local_part(int64_t, perm);

  SHARED int64_t * target = lgp_all_alloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  int64_t * ltarget = lgp_local_part(int64_t, target);
  
  /* initialize perm[i] = i,  the darts*/
  for(i = 0; i < lN; i++)
    lperm[i] = i*THREADS + MYTHREAD;

  /* initialize target[i] = -1 */
  for(i = 0; i < lM; i++)
    ltarget[i] = -1L;

  if( seed != 0 ) srand48( seed );

  lgp_barrier();
  
  double t1 = wall_seconds();
  int64_t rejects = 0;
  pkg_t pkg;
  exstack_t * ex = exstack_init(1024, sizeof(pkg_t));
  if(ex == NULL){return(NULL);}
  
  istart = iend = 0L;
  while(exstack_proceed(ex, (iend == lN))){
    i = iend;
    while(i < lN){
      int64_t r = lrand48() % M;
      pe = r % THREADS;
      pkg.idx = r/THREADS;
      pkg.val = lperm[i];
      if(exstack_push(ex, &pkg, pe) == 0L)
        break;
      i++;
    }
    iend = i;
    
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      if(ltarget[pkg.idx] == -1L){
        ltarget[pkg.idx] = pkg.val;
      }else{
        /* push this pkg back to the sender */
        exstack_push(ex, &pkg, fromth);
        //rejects++;
      }
    }

    lgp_barrier();
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      lperm[--iend] = pkg.val;
    }
    lgp_barrier();
  }
  
  lgp_barrier();
  t1 = wall_seconds() - t1;
  //T0_printf("phase 1 t1 = %8.3lf\n", t1);
  //rejects = lgp_reduce_add_l(rejects);
  //T0_printf("rejects = %ld\n", rejects);

  exstack_clear(ex);
  free(ex);

  /* now locally pack the values you have in target */
  cnt = 0;
  for(i = 0; i < lM; i++){
    if( ltarget[i] != -1L ) {
      ltarget[cnt++] = ltarget[i];
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t total = lgp_reduce_add_l(cnt);
  if(total != N){
    T0_printf("ERROR: rand_permp_exstack: total = %ld should be %ld\n", total, N);
    return(NULL);
  }

  exstack_t * ex1 = exstack_init(1024, sizeof(int64_t));
  if(ex1 == NULL){return(NULL);}

  int64_t offset = lgp_prior_add_l(cnt);
  pe = offset % THREADS;
  i = pos = 0;
  while(exstack_proceed(ex1, (i==cnt))) {
    while(i < cnt){
      if(exstack_push(ex1, &ltarget[i], pe) == 0)
        break;
      i++;
      pe++;
      if(pe == THREADS) pe = 0;
    }
    
    exstack_exchange(ex1);

    while(exstack_pop(ex1, &val, &fromth)){
      lperm[pos++] = val;
    }
  }
  
  pos = lgp_reduce_add_l(pos);
  if(pos != N){
    printf("ERROR! in rand_permp! sum of pos = %ld lN = %ld\n", pos, N);
    return(NULL);
  }
  lgp_barrier();

  exstack_clear(ex1);
  free(ex1);

  return(perm);
}



/*! \brief apply row and column permutations to a sparse matrix using exstack2 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv)
{
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;
  typedef struct pkg_rowcnt_t{
    int64_t row;
    int64_t cnt;
  }pkg_rowcnt_t;
  typedef struct pkg_inonz_t{
    int64_t i;    
    int64_t nonz;
  }pkg_inonz_t;

  sparsemat_t * Ap;
  
  int64_t i, fromth, fromth2, pe, row, lnnz;
  pkg_rowcnt_t pkg_rc;
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * lcperminv = lgp_local_part(int64_t, cperminv);
  
  //if(!MYTHREAD)printf("Permuting matrix with exstack...");
  
  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperminv rowcounts 
  /****************************************************************/
  int64_t * tmprowcnts = calloc(A->lnumrows + 1, sizeof(int64_t));
  
  exstack_t * ex = exstack_init( 1024, sizeof(pkg_rowcnt_t));
  if( ex == NULL ) return(NULL);
  lnnz = row = 0;
  while(exstack_proceed(ex, (row == A->lnumrows))) {
    while(row < A->lnumrows){
      pe = lrperminv[row] % THREADS;
      pkg_rc.row = lrperminv[row] / THREADS;
      pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
      if(exstack_push(ex, &pkg_rc, pe) == 0L)
        break;
      row++;
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg_rc, &fromth)){
      tmprowcnts[pkg_rc.row] = pkg_rc.cnt;
      lnnz += pkg_rc.cnt;
    }
  }
  lgp_barrier();
  exstack_clear(ex);
  free(ex);

  assert(A->nnz = lgp_reduce_add_l(lnnz));

  // allocate pmat to the max of the new number of nonzeros per thread  
  Ap = init_matrix(A->numrows, A->numcols, lnnz);
  if(Ap == NULL) return(NULL);
  lgp_barrier();

  // convert row counts to offsets 
  Ap->loffset[0] = 0;
  for(i = 1; i < Ap->lnumrows+1; i++)
    Ap->loffset[i] = Ap->loffset[i-1] + tmprowcnts[i-1];

  /****************************************************************/
  // re-distribute nonzeros
  // working offset: wrkoff[row] gives the first empty spot on row row  
  /****************************************************************/
  int64_t * wrkoff = calloc(A->lnumrows, sizeof(int64_t)); 
  pkg_rowcol_t pkg_nz;
  
  exstack_t * ex1 = exstack_init( 1024, sizeof(pkg_rowcol_t));
  if( ex1 == NULL )return(NULL);

  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) 
        row++;
      pkg_nz.row = lrperminv[row] / THREADS;
      pkg_nz.col = A->lnonzero[i];
      //printf("th %d: pushing (%ld, %ld) to pe %ld\n", MYTHREAD, lrperminv[row], pkg_nz.col, lrperminv[row] % THREADS);
      if( !exstack_push(ex1, &pkg_nz, lrperminv[row] % THREADS ) )        
        break;
      i++;
    }
    exstack_exchange(ex1);

    while(exstack_pop(ex1, &pkg_nz, &fromth)) {
      //printf("th %d: rcv %ld %ld\n", MYTHREAD, pkg_nz.row, pkg_nz.col);
      Ap->lnonzero[ Ap->loffset[pkg_nz.row] + wrkoff[pkg_nz.row] ] = pkg_nz.col;
      wrkoff[pkg_nz.row]++;
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++){
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%ld] = %ld trc[%ld] = %ld\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  }
  if(error){printf("ERROR! permute_matrix_exstack: error = %ld\n", error);}

  free(wrkoff);
  free(tmprowcnts);
  exstack_clear(ex1);
  free(ex1);


  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg3;
  exstack_t * ex2 = exstack_init( 1024, sizeof(pkg_inonz_t));
  if( ex2 == NULL ) return(NULL);
  i=0;
  while(exstack_proceed(ex2,(i == Ap->lnnz))){
    while(i < Ap->lnnz){     // request the new name for this non-zero
      pkg3.i = i;
      pkg3.nonz = Ap->lnonzero[i] / THREADS;
      pe = Ap->lnonzero[i] % THREADS;
      if( !exstack_push(ex2, &pkg3, pe) )
        break;
      i++;
    }

    exstack_exchange(ex2);

    while(exstack_pop(ex2, &pkg3, &fromth)){ // echo the requests back to requester 
      pkg3.nonz = lcperminv[pkg3.nonz];
      exstack_push(ex2, &pkg3, fromth);
    }

    lgp_barrier();

    exstack_exchange(ex2);

    while(exstack_pop(ex2, &pkg3, &fromth)){ //process the echo to your requests
      Ap->lnonzero[pkg3.i] = pkg3.nonz;
    }
  }

  exstack_clear(ex2);
  free(ex2);

  lgp_barrier();
  
  //if(!MYTHREAD)printf("done\n");

  return(Ap);

}

/*! \brief produce the transpose of a sparse matrix using exstack2
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_exstack(sparsemat_t * A)
{
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, *idxp;

  //if(!MYTHREAD)printf("Exstack version of matrix transpose...");
  
  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();
  
  exstack_t * ex = exstack_init( 1024, sizeof(int64_t));
  if( ex == NULL ) return(NULL);

  lnnz = i = 0;
  while(exstack_proceed(ex, (i == A->lnnz))){
    while(i < A->lnnz){
      col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i] % THREADS;
      if( !exstack_push(ex, &col, pe) )
        break;
      i++;

    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &idx, &fromth)){
      lcounts[idx]++; 
      lnnz++;
    }
  }
  exstack_clear(ex);
  free(ex);

  int64_t sum = lgp_reduce_add_l(lnnz);
  assert( A->nnz == sum ); 
  
  sparsemat_t * At = init_matrix(A->numcols, A->numrows, lnnz);
  if(!At){printf("ERROR: transpose_matrix: init_matrix failed!\n");return(NULL);}

  /* convert colcounts to offsets */
  At->loffset[0] = 0;  
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + lcounts[i-1];
    
  lgp_barrier();
  
  /* redistribute the nonzeros */
  int64_t *wrkoff = calloc(At->lnumrows, sizeof(int64_t));
  if(!wrkoff) {printf("ERROR: transpose_matrix: wrkoff alloc fail!\n"); return(NULL);}

  pkg_rowcol_t pkg_nz;
  
  exstack_t * ex1 = exstack_init( 1024, sizeof(pkg_rowcol_t));
  if( ex1 == NULL ) return(NULL);

  uint64_t numtimespop=0;
  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) 
        row++;
      pkg_nz.row = row * THREADS + MYTHREAD;
      pkg_nz.col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i]%THREADS;
      if( exstack_push(ex1, &pkg_nz, pe) == 0 )
        break;
      i++;
    }

    exstack_exchange(ex1);

    while(exstack_pop(ex1, &pkg_nz, &fromth)){
      numtimespop++;
      At->lnonzero[ At->loffset[pkg_nz.col] + wrkoff[pkg_nz.col] ] = pkg_nz.row;
      wrkoff[pkg_nz.col]++;
    }
  }
  
  lgp_barrier();
  exstack_clear(ex1);
  free(ex1);

  //if(!MYTHREAD)printf("done\n");

  numtimespop = lgp_reduce_add_l(numtimespop);
  if(numtimespop != A->nnz ){
    printf("ERROR: numtimespop %d \n", numtimespop);
    printf("%d wrkoff %d\n", MYTHREAD, wrkoff[0]);
    return(NULL);
  }

  for(i = 0; i < At->lnumrows; i++){
    if(wrkoff[i] != lcounts[i] ) {
      printf("ERROR: %d wrkoff[%d] = %d !=  %d = lcounts[%d]\n", MYTHREAD, i, wrkoff[i],lcounts[i],i);
      return(NULL);
    }
  }
  
  free(wrkoff);
  free(lcounts);
  return(At);
}

