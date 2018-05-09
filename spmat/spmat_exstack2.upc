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
/*! \file spmat_exstack2.upc
 * \brief spmat routines written using exstack2
 */ 

/*! \brief create a global int64_t array with a uniform random permutation using a exstack2 implementation
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
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_exstack2(int64_t N, int seed)
{
  T0_printf("Don't have an exstack2 version for rand_permp, returning rand_permp_atomic\n");
  return( rand_permp_atomic( N, seed) );
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
sparsemat_t * permute_matrix_exstack2(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv)
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
  
  //T0_printf("Permuting matrix with exstack2\n");
  
  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperminv rowcounts 
  /****************************************************************/
  int64_t * tmprowcnts = calloc(A->lnumrows + 1, sizeof(int64_t));

  exstack2_t * ex2 = exstack2_init(1024, sizeof(pkg_rowcnt_t));
  if( ex2 == NULL ){return(NULL);}
  lnnz = row = 0;
  while(exstack2_proceed(ex2, (row == A->lnumrows))) {
    while(row < A->lnumrows){
      pe = lrperminv[row] % THREADS;
      pkg_rc.row = lrperminv[row] / THREADS;
      pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
      if( !exstack2_push(ex2, &pkg_rc, pe) )
        break;
      row++;
    }
    while(exstack2_pop(ex2, &pkg_rc, &fromth)){
      tmprowcnts[pkg_rc.row] = pkg_rc.cnt;
      lnnz += pkg_rc.cnt;
    }
  }
  lgp_barrier();
  exstack2_clear(ex2);  
  free(ex2);

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
  
  exstack2_t * exr = exstack2_init(1024, sizeof(pkg_rowcol_t));
  if( exr == NULL )return(NULL);

  i = row = 0;
  while(exstack2_proceed(exr, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) // skip empty rows 
        row++; 
      pkg_nz.row = lrperminv[row] / THREADS;
      pkg_nz.col = A->lnonzero[i];
      //printf("th %d: pushing (%ld, %ld) to pe %ld\n", MYTHREAD, lrperminv[row], pkg_nz.col, lrperminv[row] % THREADS);
      if( !exstack2_push(exr, &pkg_nz, lrperminv[row] % THREADS ) )
        break;
      i++;
    }

    while(exstack2_pop(exr, &pkg_nz, &fromth)) {
      //printf("th %d: rcv %ld %ld\n", MYTHREAD, pkg_nz.row, pkg_nz.col);
      Ap->lnonzero[ Ap->loffset[pkg_nz.row] + wrkoff[pkg_nz.row] ] = pkg_nz.col;
      wrkoff[pkg_nz.row]++;
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++)
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%ld] = %ld trc[%ld] = %ld\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  if(error){printf("ERROR! permute_matrix_exstack: error = %ld\n", error);}

  free(wrkoff);
  free(tmprowcnts);
  exstack2_clear(exr);
  free(exr);

  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg_r, pkg_e;
  exstack2_t *ex2_r = exstack2_init(1024, sizeof(pkg_inonz_t));
  exstack2_t *ex2_e = exstack2_init(1024, sizeof(pkg_inonz_t));
  if( (ex2_r == NULL) || (ex2_e == NULL) ) return(NULL);
  i=0;
  while(exstack2_proceed(ex2_r,(i == Ap->lnnz)) || exstack2_proceed(ex2_e, ex2_r->all_done)){
    while(i < Ap->lnnz){     // request the new name for this non-zero
      pkg_r.i = i;
      pkg_r.nonz = Ap->lnonzero[i] / THREADS;
      pe = Ap->lnonzero[i] % THREADS;
      if(!exstack2_push(ex2_r, &pkg_r, pe))
        break;
      i++;
    }
    while(exstack2_pop(ex2_r, &pkg_e, &fromth2)){ 
      pkg_r.i = pkg_e.i;
      pkg_r.nonz = lcperminv[pkg_e.nonz];
      if( !exstack2_push(ex2_e, &pkg_r, fromth2)){
        exstack2_unpop(ex2_r);
        break;
      }
    }
    while(exstack2_pop(ex2_e, &pkg_r, &fromth)){ 
      Ap->lnonzero[pkg_r.i] = pkg_r.nonz;
    }
  }

  exstack2_clear(ex2_r);
  exstack2_clear(ex2_e);
  free(ex2_r);
  free(ex2_e);

  lgp_barrier();
  
  //if(!MYTHREAD)printf("done\n");

  return(Ap);

}

/*! \brief produce the transpose of a sparse matrix using exstack2
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_exstack2(sparsemat_t * A)
{
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, *idxp;

  //if(!MYTHREAD)printf("Exstack2 version of matrix transpose\n");
  
  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();
  
  exstack2_t * ex2c = exstack2_init(1024, sizeof(int64_t));

  if( ex2c == NULL ) return(NULL);

  lnnz = i = 0;
  while(exstack2_proceed(ex2c, (i == A->lnnz))){
    while(i < A->lnnz){
      col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i] % THREADS;
      if( !exstack2_push(ex2c, &col, pe) )
        break;
      i++;
    }
    while(exstack2_pop(ex2c, &idx, &fromth)){
      lcounts[idx]++; 
      lnnz++;
    }
  }
  exstack2_clear(ex2c);

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
  
  exstack2_t * ex2r = exstack2_init(1024, sizeof(pkg_rowcol_t));
  if( ex2c == NULL ) return(NULL);

  uint64_t numtimespop=0;
  i = row = 0;
  while(exstack2_proceed(ex2r, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) 
        row++;
      pkg_nz.row = row * THREADS + MYTHREAD;
      pkg_nz.col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i]%THREADS;
      if( !exstack2_push(ex2r, &pkg_nz, pe) )
        break;
      i++;
    }
    while(exstack2_pop(ex2r, &pkg_nz, &fromth)){
      numtimespop++;
      At->lnonzero[ At->loffset[pkg_nz.col] + wrkoff[pkg_nz.col] ] = pkg_nz.row;
      wrkoff[pkg_nz.col]++;
    }
  }
  
  lgp_barrier();
  exstack2_clear(ex2r);

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

