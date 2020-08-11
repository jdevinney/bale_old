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

#include <spmat.h>
#include <exstack.h>
/*! \file spmat_exstack2.upc
 * \brief spmat routines written using exstack2
 */ 

/*! \brief create a global int64_t array with a uniform random permutation using a exstack2 implementation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \param buf_cnt the number of packets in an exstack2 buffer
 * \return the permutation
 * 
 * This is a collective call.
 * this implements a variant of the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * 
 *
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * 
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_exstack2(int64_t N, int seed, int64_t buf_cnt) {
  // I needed to use a steady conveyor in the conveyor implementation
  // of the dart throwing algorithm (and exstack2 doesn't have this
  // feature). For now, I just call the exstack version but below you
  // can see my attempt at an exstack2 implementation. It ultimately
  // failed... in the main dart throwing loop, you don't know how many
  // rounds you are going to need, but some of those rounds are quite
  // small (with only a few darts exchanged) and so we need a way to
  // force exstack2 to send even though its buffers are not full.  We
  // worked around this in toposort by having exstack2 resets in every
  // round. I could do that... but I am going to put this off. I am
  // not sure it is worth the effort.
  
  return(rand_permp_exstack(N, seed, buf_cnt));

#if 0
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
  exstack2_t * ex_throw = exstack2_init(buf_cnt, sizeof(pkg_t));
  exstack2_t * ex_reply = exstack2_init(buf_cnt, sizeof(int64_t));
  if(ex_throw == NULL){return(NULL);}
  if(ex_reply == NULL){return(NULL);}

  int64_t more = 1;
  int64_t hits = 0;
  iend = 0;
  while(more || exstack2_proceed(ex_reply, !more)){
    i = iend;
    while(i < lN){
      int64_t r = lrand48() % M;
      pe = r % THREADS;
      pkg.idx = r/THREADS;
      pkg.val = lperm[i];
      if(!exstack2_push(ex_throw, &pkg, pe))
        break;
      i++;
    }
    iend = i;
    
    while(exstack2_pop(ex_throw, &pkg, &fromth)){
      if(ltarget[pkg.idx] == -1L){
        val = pkg.val;
        if(!exstack2_push(ex_reply, &val, fromth)){
          exstack2_unpop(ex_throw);
          break;
        } 
        ltarget[pkg.idx] = pkg.val;
      }else{
        val = -(pkg.val + 1);
        if(!exstack2_push(ex_reply, &val, fromth)){
          exstack2_unpop(ex_throw);
          rejects++;
          break;
        }
      }
    }
    
    while(exstack2_pop(ex_reply, &val, NULL)){
      if(val < 0L)
        lperm[--iend] = -val - 1;
      else
        hits++;
    }
    fprintf(stderr,"MYTHREAD = %"PRId64" hits = %"PRId64"\n", MYTHREAD, hits);
    more = exstack2_proceed(ex_throw, (hits == lN));
  }
  
  lgp_barrier();
  t1 = wall_seconds() - t1;

  exstack2_clear(ex_throw);
  exstack2_reset(ex_reply);
  T0_fprintf(stderr,"Made it!\n");
  
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
    T0_printf("ERROR: rand_permp_exstack2: total = %"PRId64" should be %"PRId64"\n", total, N);
    return(NULL);
  }
  
  int64_t offset = lgp_prior_add_l(cnt);
  pe = offset % THREADS;
  i = pos = 0;
  while(exstack2_proceed(ex_reply, (i==cnt))) {
    while(i < cnt){
      if(!exstack2_push(ex_reply, &ltarget[i], pe))
        break;
      i++;
      pe++;
      if(pe == THREADS) pe = 0;
    }
    
    while(exstack2_pop(ex_reply, &val, &fromth)){
      lperm[pos++] = val;
    }
  }

  pos = lgp_reduce_add_l(pos);
  if(pos != N){
    printf("ERROR! in rand_permp_exstack2! sum of pos = %"PRId64" lN = %"PRId64"\n", pos, N);
    return(NULL);
  }
  lgp_barrier();
  
  exstack2_clear(ex_reply);

  return(perm);
#endif
}



/*! \brief apply row and column permutations to a sparse matrix using exstack2 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \param buf_cnt the number of packets in an exstack2 buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_exstack2(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt) {
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

  exstack2_t * ex2 = exstack2_init(buf_cnt, sizeof(pkg_rowcnt_t));
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
  Ap = init_matrix(A->numrows, A->numcols, lnnz, (A->value != NULL));
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
  
  exstack2_t * exr = exstack2_init(buf_cnt, sizeof(pkg_rowcol_t));
  if( exr == NULL )return(NULL);

  i = row = 0;
  while(exstack2_proceed(exr, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) // skip empty rows 
        row++; 
      pkg_nz.row = lrperminv[row] / THREADS;
      pkg_nz.col = A->lnonzero[i];
      //printf("th %d: pushing (%"PRId64", %"PRId64") to pe %"PRId64"\n", MYTHREAD, lrperminv[row], pkg_nz.col, lrperminv[row] % THREADS);
      if( !exstack2_push(exr, &pkg_nz, lrperminv[row] % THREADS ) )
        break;
      i++;
    }

    while(exstack2_pop(exr, &pkg_nz, &fromth)) {
      //printf("th %d: rcv %"PRId64" %"PRId64"\n", MYTHREAD, pkg_nz.row, pkg_nz.col);
      Ap->lnonzero[ Ap->loffset[pkg_nz.row] + wrkoff[pkg_nz.row] ] = pkg_nz.col;
      wrkoff[pkg_nz.row]++;
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++)
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%"PRId64"] = %"PRId64" trc[%"PRId64"] = %"PRId64"\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  if(error){printf("ERROR! permute_matrix_exstack: error = %"PRId64"\n", error);}

  free(wrkoff);
  free(tmprowcnts);
  exstack2_clear(exr);
  free(exr);

  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg_r, pkg_e;
  exstack2_t *ex2_r = exstack2_init(buf_cnt, sizeof(pkg_inonz_t));
  exstack2_t *ex2_e = exstack2_init(buf_cnt, sizeof(pkg_inonz_t));
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
 * \param buf_cnt the number of packets in an exstack2 buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_exstack2(sparsemat_t * A, int64_t buf_cnt) {
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
  
  exstack2_t * ex2c = exstack2_init(buf_cnt, sizeof(int64_t));

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
  
  sparsemat_t * At = init_matrix(A->numcols, A->numrows, lnnz, (A->value != NULL));
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
  
  exstack2_t * ex2r = exstack2_init(buf_cnt, sizeof(pkg_rowcol_t));
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
    printf("ERROR: numtimespop %"PRId64" \n", numtimespop);
    printf("%d wrkoff %"PRId64"\n", MYTHREAD, wrkoff[0]);
    return(NULL);
  }

  for(i = 0; i < At->lnumrows; i++){
    if(wrkoff[i] != lcounts[i] ) {
      printf("ERROR: %d wrkoff[%"PRId64"] = %"PRId64" !=  %"PRId64" = lcounts[%"PRId64"]\n", MYTHREAD, i, wrkoff[i],lcounts[i],i);
      return(NULL);
    }
  }
  
  
  free(wrkoff);
  free(lcounts);
  return(At);
}

