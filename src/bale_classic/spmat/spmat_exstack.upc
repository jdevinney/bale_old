/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
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
#include <sys/stat.h> // for mkdir()

/*! \file spmat_exstack.upc
 * \brief spmat routines written using exstack
 */ 

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \param buf_cnt number of items in an exstack buffer
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
SHARED int64_t * rand_permp_exstack(int64_t N, int seed, int64_t buf_cnt) {
  int ret;
  int64_t i, j, cnt, pe, pos, fromth, iend;
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
  exstack_t * ex = exstack_init(buf_cnt, sizeof(pkg_t));
  if(ex == NULL){return(NULL);}
  
  iend = 0L;
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
  //T0_printf("rejects = %"PRId64"\n", rejects);

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
    T0_printf("ERROR: rand_permp_exstack: total = %"PRId64" should be %"PRId64"\n", total, N);
    return(NULL);
  }

  exstack_t * ex1 = exstack_init(buf_cnt, sizeof(int64_t));
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
    printf("ERROR! in rand_permp! sum of pos = %"PRId64" lN = %"PRId64"\n", pos, N);
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
 * \param buf_cnt number of items in an exstack buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperminv, SHARED int64_t * cperminv, int64_t buf_cnt) {
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
  
  exstack_t * ex = exstack_init( buf_cnt, sizeof(pkg_rowcnt_t));
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

  assert(A->nnz == lgp_reduce_add_l(lnnz));

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
  
  exstack_t * ex1 = exstack_init( buf_cnt, sizeof(pkg_rowcol_t));
  if( ex1 == NULL )return(NULL);

  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) 
        row++;
      pkg_nz.row = lrperminv[row] / THREADS;
      pkg_nz.col = A->lnonzero[i];
      //printf("th %d: pushing (%"PRId64", %"PRId64") to pe %"PRId64"\n", MYTHREAD, lrperminv[row], pkg_nz.col, lrperminv[row] % THREADS);
      if( !exstack_push(ex1, &pkg_nz, lrperminv[row] % THREADS ) )        
        break;
      i++;
    }
    exstack_exchange(ex1);

    while(exstack_pop(ex1, &pkg_nz, &fromth)) {
      //printf("th %d: rcv %"PRId64" %"PRId64"\n", MYTHREAD, pkg_nz.row, pkg_nz.col);
      Ap->lnonzero[ Ap->loffset[pkg_nz.row] + wrkoff[pkg_nz.row] ] = pkg_nz.col;
      wrkoff[pkg_nz.row]++;
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++){
    if(wrkoff[i] != tmprowcnts[i]){printf("T%d: w[%"PRId64"] = %"PRId64" trc[%"PRId64"] = %"PRId64"\n", MYTHREAD, i, wrkoff[i], i, tmprowcnts[i]);error++;}
  }
  if(error){printf("ERROR! permute_matrix_exstack: error = %"PRId64"\n", error);}

  free(wrkoff);
  free(tmprowcnts);
  exstack_clear(ex1);
  free(ex1);


  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg3;
  exstack_t * ex2 = exstack_init( buf_cnt, sizeof(pkg_inonz_t));
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
 * \param buf_cnt number of items in an exstack buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_exstack(sparsemat_t * A, int64_t buf_cnt) {
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, *idxp;

  //T0_fprintf(stderr,"Exstack version of matrix transpose...");fflush(stderr);
  
  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();
  
  exstack_t * ex = exstack_init( buf_cnt, sizeof(int64_t));
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
  
  exstack_t * ex1 = exstack_init( buf_cnt, sizeof(pkg_rowcol_t));
  if( ex1 == NULL ) return(NULL);

  uint64_t numtimespop=0;
  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ){
        row++;
        assert(row < A->lnumrows);
      }
      pkg_nz.row = row * THREADS + MYTHREAD;
      pkg_nz.col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i]%THREADS;
      if( exstack_push(ex1, &pkg_nz, pe) == 0 )
        break;
      i++;
    }
    //fprintf(stderr,"exstack...i = %ld/%ld\n", i, A->lnnz);fflush(stderr);
    exstack_exchange(ex1);

    while(exstack_pop(ex1, &pkg_nz, &fromth)){
      numtimespop++;
      At->lnonzero[ At->loffset[pkg_nz.col] + wrkoff[pkg_nz.col] ] = pkg_nz.row;
      wrkoff[pkg_nz.col]++;
    }
    //fprintf(stderr,"num_times_pop %ld\n", numtimespop);
  }
  
  lgp_barrier();
  exstack_clear(ex1);
  free(ex1);

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





/****************************************************************************/
/*! \brief Write a sparsemat_t to disk as a sparse matrix dataset.
 *
 * This routine is collective and it's return is single-valued.
 *
 * PE i will collect the ith chunk of rows (and their nonzeros) to write. 
 * 
 * \param dirname The directory where the sparsemat will be written (must be single-valued).
 * \param mat The sparsemat_t to be written.
 * \param buf_cnt number of packets in an exstack buffer
 * \return 0 on SUCCESS nonzero on ERROR.
 */
/****************************************************************************/
int64_t write_sparse_matrix_exstack( char * dirname, sparsemat_t * mat, int64_t buf_cnt){

  int64_t *rowcnt_buf;
  int64_t *nz_buf;
  uint64_t row, col, cnt, max, lrow, pe;
  uint64_t * current_row_to_th, * last_row;

  int64_t i, j, k, rows_per_pass, nr = mat->numrows;
  int64_t pass, num_passes, first_pe;
  int64_t ret, error;
  exstack_t * ex;

  //T0_fprintf(stderr,"\n***** Begining to write sparse matrix %s *****\n", dirname);

  //T0_fprintf(stderr,"Writing out a sparsemat_t with ");  
  //T0_fprintf(stderr,"%lu rows %lu columns and %lu nonzeros\n",
  //         mat->numrows, mat->numcols, mat->nnz);

  lgp_barrier();

  /* get max row density */
  max = 0;
  for(row = 0; row < mat->lnumrows; row++){
    cnt = mat->loffset[row + 1] - mat->loffset[row];
    max = (cnt > max ? cnt : max);
  }

  max = lgp_reduce_max_l(max);
  if(max==0) max = 1;

  rows_per_pass = buf_cnt/(max+1);
  while(rows_per_pass == 0){
    buf_cnt *= 2;
    rows_per_pass = buf_cnt/(max+1);
  }
  //T0_fprintf("rows_per_pass = %"PRId64", max = %lu\n", rows_per_pass, max);
  
  /* allocate space */
  rowcnt_buf         = calloc(THREADS, sizeof(int64_t));  
  current_row_to_th  = calloc(THREADS, sizeof(uint64_t));
  last_row           = calloc(THREADS, sizeof(uint64_t));
  nz_buf             = calloc(max*THREADS, sizeof(int64_t));
  
  ex = exstack_init((max + 1)*rows_per_pass, sizeof(uint64_t));
  if(!ex){
    fprintf(stderr,"write_sparse_matrix_exstack: exstack_init failed\n");lgp_barrier();
    return(-1);
  }
  
  /* first figure out the first row each PE will write. */
  current_row_to_th[0] = 0L;
  for(i = 1; i < THREADS; i++)
    current_row_to_th[i] =  current_row_to_th[i-1] + (nr + THREADS - (i-1) - 1)/THREADS;

  /* now figure out which row on this PE should get sent first to each other PE */
  /* and determine which PE owns the first row you will be writing */
  for(i = 0; i < THREADS; i++){
    if(i == MYTHREAD)
      first_pe = current_row_to_th[i] % THREADS; /* this is the pe who will send your first row to write */
    while((current_row_to_th[i] % THREADS) != MYTHREAD)
      current_row_to_th[i]++;
    if(i > 0)
      last_row[i-1] = current_row_to_th[i];
  }
  last_row[THREADS - 1] = mat->numrows;

  //for(i = 0; i < THREADS; i++)
    //T0_fprintf("%d first_pe = %"PRId64" last_row[%d]= %"PRId64" current_row_to_th[%"PRId64"] = %"PRId64"\n", MYTHREAD, first_pe, i, last_row[i], i, current_row_to_th[i]);
  
  
  /* make the directory */  
  mkdir(dirname, 0770);

  /* write the metadata file */
  write_sparse_matrix_metadata(dirname, mat);
  
  /* open rowcnt file */
  char filename[64];
  sprintf(filename, "%s/rowcnt_%d", dirname, MYTHREAD);
  FILE * rcfp = fopen(filename, "wb");
  
  /* open nonzero file */
  sprintf(filename, "%s/nonzero_%d", dirname, MYTHREAD);
  FILE * nzfp = fopen(filename, "wb");
  
  /*********************************/
  /*        WRITE OUT DATASET      */
  /*********************************/
  //T0_fprintf(stderr, "Writing data...");
  
  int64_t room;
  num_passes = (nr + THREADS - 1)/THREADS;
  num_passes /= (rows_per_pass * THREADS);
  //T0_fprintf(stderr,"num_passes = %"PRId64"\n", num_passes);
  
  error = 0;
  int imdone = 0;
  int64_t recs_written, nnz;
  int64_t rc_pos = 0;
  int64_t nz_pos = 0;
  while(exstack_proceed(ex, imdone)){

    imdone = 1;

    /* add rows_per_pass rows to everyone's buffers */
    for(i = 0; i < THREADS; i++){
      row = current_row_to_th[i];
      lrow = row/THREADS;
      for(k = 0; (row < last_row[i]) && (k < rows_per_pass); row+=THREADS, lrow++,k++){
        cnt = mat->loffset[lrow+1] - mat->loffset[lrow];
        room = exstack_push(ex, &cnt, i);
        if(room < 0L){
          fprintf(stderr, "ERROR: write_sparse_matrix: Trying to push cnt onto full exstack!\n");
          error = 1;
        }
        for(j = mat->loffset[lrow]; j < mat->loffset[lrow+1]; j++){
          col = mat->lnonzero[j];
          room = exstack_push(ex, &col, i);
          if(room < 0L){
            fprintf(stderr,"ERROR: write_sparse_matrix: Trying to push cols onto full exstack!\n");
            error = 1;
          }
        }
      }
      current_row_to_th[i] = row;
      /* you are not done if you did not finish pushing to some PE */
      if(current_row_to_th[i] < last_row[i])
        imdone = 0;
    }
      
    exstack_exchange(ex);

    for(k = 0; k < rows_per_pass; k++){
      pe = first_pe;
      nz_pos = rc_pos = 0L;
      for(i = first_pe; i < first_pe + THREADS; i++, pe++){
        pe = (pe == THREADS ? 0 : pe);
        if(exstack_pop_thread(ex, &cnt, pe)){
          for(j = 0; j < cnt; j++){
            exstack_pop_thread(ex, &col, pe);
            nz_buf[nz_pos++] = col;
          }
          rowcnt_buf[rc_pos++] = cnt;
        }
      }

      /*******************************/
      /*   FIRST WRITE ROW COUNTS    */
      /*******************************/
      /* write row counts */
      recs_written =  fwrite(rowcnt_buf, sizeof(int64_t), rc_pos, rcfp);
      if(recs_written != rc_pos){
        error = 1;
        fprintf(stderr, "write_sparse_matrix_exstack: recs_written != numrows %"PRId64" %"PRId64" on PE %d\n", recs_written, rc_pos, MYTHREAD);
      }
        
      /* count up the number of nonzeros to be written */
      nnz = 0;
      for(i = 0; i < rc_pos; i++)
        nnz += rowcnt_buf[i];
        
      /* write the nonzeros */
      recs_written = fwrite(nz_buf, sizeof(int64_t), nz_pos, nzfp);
      if(recs_written != nz_pos){
        error = 1;
        fprintf(stderr, "write_sparse_matrix_exstack: recs_written != nnz %"PRId64" %"PRId64" on PE %d\n", recs_written, nz_pos, MYTHREAD);
      }
    }
    lgp_barrier();
  }

  
  /* make sure all threads wrote cleanly */
  ret = lgp_reduce_add_l(error);
  if(ret){
    T0_fprintf(stderr,"ERROR: write_sparse_matrix_exstack: error in main loop\n");
    lgp_barrier();
    return(-1);
  }

  /* finished writing */
  lgp_barrier();

  free(rowcnt_buf);
  free(nz_buf);
  free(current_row_to_th);
  free(last_row);
  exstack_clear(ex);

  lgp_barrier();

  //T0_fprintf(stderr,"***** End writing sparse matrix %s *****\n", dirname);
  return(0);
}
