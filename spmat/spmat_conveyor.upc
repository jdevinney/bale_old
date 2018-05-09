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
/*! \file spmat_conveyor.upc
 * \brief spmat routines written with conveyors
 */ 

/*! \brief create a global int64_t array with a uniform random permutation using conveyors
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 *  
 *
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * 
 */
SHARED int64_t * rand_permp_conveyor(int64_t N, int seed)
{
  T0_printf("Don't have an conveyor version for rand_permp, returning rand_permp_atomic\n");
  return( rand_permp_atomic( N, seed) );
}


/*! \brief apply row and column permutations to a sparse matrix using conveyors 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * permute_matrix_conveyor(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv)
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
  int64_t * lrperminv = lgp_local_part(int64_t, rperminv);
  int64_t * lcperminv = lgp_local_part(int64_t, cperminv);

  size_t n_local = 1;
  char *number = getenv("CONVEY_LOCAL_PES");
  if(number)
    n_local = atoi(number);
  //T0_printf("Permuting matrix with conveyors\n");
  
  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperminv rowcounts 
  /****************************************************************/
  int64_t * tmprowcnts = calloc(A->lnumrows + 1, sizeof(int64_t));

  pkg_rowcnt_t pkg_rc;
  pkg_rowcnt_t pkgrc_p;
  convey_t* cnv_rc = convey_new(sizeof(pkg_rowcnt_t), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);

  convey_begin(cnv_rc);
  lnnz = row = 0;
  while(convey_advance(cnv_rc, (row == A->lnumrows))) {
    for( ;row < A->lnumrows; row++){
      pe = lrperminv[row] % THREADS;
      pkg_rc.row = lrperminv[row] / THREADS;
      pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
      if( !convey_push(cnv_rc, &pkg_rc, pe) )
        break;
    }
    while( convey_pull(cnv_rc, &pkgrc_p, NULL) == convey_OK ){
      tmprowcnts[pkgrc_p.row] = pkgrc_p.cnt;
      lnnz += pkgrc_p.cnt;
    }
  }
  lgp_barrier();
  convey_free(cnv_rc);  

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
  pkg_rowcol_t pkg_nz, pkgnz_p;
  
  convey_t* cnv_nz = convey_new(sizeof(pkg_rowcol_t), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  convey_begin(cnv_nz);

  i = row = 0;
  while(convey_advance(cnv_nz, (i == A->lnnz))){
    for( ;i < A->lnnz; i++){
      while( i == A->loffset[row+1] ) // skip empty rows 
        row++; 
      pkg_nz.row = lrperminv[row] / THREADS;
      pkg_nz.col = A->lnonzero[i];
      pe = lrperminv[row] % THREADS;
      if( !convey_push(cnv_nz, &pkg_nz, pe) )
        break;
    }

    while( convey_pull(cnv_nz, &pkgnz_p, NULL) == convey_OK) {
      Ap->lnonzero[ Ap->loffset[pkgnz_p.row] + wrkoff[pkgnz_p.row] ] = pkgnz_p.col;
      wrkoff[pkgnz_p.row]++;
    }
  }
  lgp_barrier();
  convey_free(cnv_nz);

  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++)
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%ld] = %ld trc[%ld] = %ld\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  if(error){printf("ERROR! permute_matrix_conveyor: error = %ld\n", error);}

  free(wrkoff);
  free(tmprowcnts);

  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg_r, pkg_e, pkg_p;
  convey_t* cnv_r = convey_new(sizeof(pkg_inonz_t), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  assert( cnv_r != NULL );
  convey_t* cnv_e = convey_new(sizeof(pkg_inonz_t), SIZE_MAX, n_local, NULL, 0);
  assert( cnv_e != NULL );
  convey_begin(cnv_r);
  convey_begin(cnv_e);
  bool more;
  i=0;
  while( more = convey_advance(cnv_r,(i == Ap->lnnz)), more | convey_advance(cnv_e, !more) ){
    for( ; i < Ap->lnnz; i++){
      pkg_r.i = i;
      pkg_r.nonz = Ap->lnonzero[i] / THREADS;
      pe = Ap->lnonzero[i] % THREADS;
      if( !convey_push(cnv_r, &pkg_r, pe) )
        break;
    }
    while( convey_pull(cnv_r, &pkg_p, &fromth2) == convey_OK ){ 
      pkg_r.i = pkg_p.i;
      pkg_r.nonz = lcperminv[pkg_p.nonz];
      if( !convey_push(cnv_e, &pkg_r, fromth2) ){
        convey_unpull(cnv_r);
        break;
      }
    }
    while( convey_pull(cnv_e, &pkg_p, NULL) == convey_OK ){ 
      Ap->lnonzero[pkg_p.i] = pkg_p.nonz;
    }
  }

  lgp_barrier();
  convey_free(cnv_e);
  convey_free(cnv_r);
  
  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using conveyors
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * transpose_matrix_conveyor(sparsemat_t * A)
{
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  size_t n_local = 1; // set to cores/socket
  char *number = getenv("CONVEY_LOCAL_PES");
  if(number)
    n_local = atoi(number);

  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, idxp;

  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();
  
  convey_t* cnv_cnt = convey_new(sizeof(long), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  convey_begin(cnv_cnt);

  lnnz = i = 0;
  while(convey_advance(cnv_cnt, (i == A->lnnz))){
    for( ;i < A->lnnz; i++){
      col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i] % THREADS;
      if( !convey_push(cnv_cnt, &col, pe) )
        break;
    }
    while( convey_pull(cnv_cnt, &idxp, NULL) == convey_OK ) {
      lcounts[idxp]++; 
      lnnz++;
    }
  }
  convey_free(cnv_cnt);

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

  pkg_rowcol_t pkg_nz, pkg_p;
  
  convey_t* cnv_rd = convey_new(sizeof(pkg_rowcol_t), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  convey_begin(cnv_rd);

  uint64_t numtimespop=0;
  i = row = 0;
  while(convey_advance(cnv_rd, (i == A->lnnz))){
    for( ; i < A->lnnz; i++){
      while( i == A->loffset[row+1] ) 
        row++;
      pkg_nz.row = row * THREADS + MYTHREAD;
      pkg_nz.col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i]%THREADS;
      if( !convey_push(cnv_rd, &pkg_nz, pe) )
        break;
    }
    while(convey_pull(cnv_rd, &pkg_p, NULL ) == convey_OK){
      numtimespop++;
      At->lnonzero[ At->loffset[pkg_p.col] + wrkoff[pkg_p.col] ] = pkg_p.row;
      wrkoff[pkg_p.col]++;
    }
  }
  
  lgp_barrier();
  convey_free(cnv_rd);
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

