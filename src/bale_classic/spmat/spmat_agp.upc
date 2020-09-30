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
/*! \file spmat_agp.upc
 * \brief Sparse matrix support functions implemented with global addresses and atomics
 */
#include <spmat.h>
#include <sys/stat.h>   // for mkdir()
#include <fcntl.h>

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
SHARED int64_t * rand_permp_agp(int64_t N, int seed) {  
  int64_t * ltarget, *lperm;
  int64_t r, i, j;
  int64_t pos, numdarts, numtargets, lnumtargets;

  lgp_rand_seed(seed);

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
    r = lgp_rand_int64(M);
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
 * \param rperm pointer to the global array holding the row permutation
 * \param cperm pointer to the global array holding the column permutation
 * rperm[i] = j means that row i of A goes to row j in matrix Ap
 * cperm[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_agp(sparsemat_t *A, SHARED int64_t *rperm, SHARED int64_t *cperm) {
  //T0_printf("Permuting matrix with single puts\n");
  int weighted = (A->value != NULL);
  int64_t i, j, col, row, pos;
  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  SHARED int64_t * rperminv = lgp_all_alloc(A->numrows, sizeof(int64_t));
  if( rperminv == NULL ) return(NULL);
  int64_t *lrperminv = lgp_local_part(int64_t, rperminv);

  //compute rperminv from rperm 
  for(i=0; i < A->lnumrows; i++){
    lgp_put_int64(rperminv, lrperm[i], i*THREADS + MYTHREAD);
  }

  lgp_barrier();
  
  int64_t cnt = 0, off, nxtoff;
  for(i = 0; i < A->lnumrows; i++){
    row = lrperminv[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    cnt += nxtoff - off;
  }
  lgp_barrier();

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, cnt, weighted);
  
  // fill in permuted rows
  Ap->loffset[0] = pos = 0;
  for(i = 0; i < Ap->lnumrows; i++){
    row = lrperminv[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    lgp_memget(&Ap->lnonzero[pos], A->nonzero,
               (nxtoff-off)*sizeof(int64_t), off*THREADS + row%THREADS);    
    if(weighted)
      lgp_memget(&Ap->lvalue[pos], A->value,
                 (nxtoff-off)*sizeof(double), off*THREADS + row%THREADS);
    pos += nxtoff - off;
    //for(j = off; j < nxtoff; j++){
    //Ap->lnonzero[pos++] = lgp_get_int64(A->nonzero, j*THREADS + row%THREADS);
    //}
    Ap->loffset[i+1] = pos;
  }
  
  assert(pos == cnt);
  
  lgp_barrier();
  
  // finally permute column indices
  for(i = 0; i < Ap->lnumrows; i++){
    for(j = Ap->loffset[i]; j < Ap->loffset[i+1]; j++){
      Ap->lnonzero[j] = lgp_get_int64(cperm, Ap->lnonzero[j]);      
    }
  }
  lgp_barrier();

  lgp_all_free(rperminv);
  sort_nonzeros(Ap);

  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_agp(sparsemat_t *A) {
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
  int weighted = (A->value != NULL);
  At = init_matrix(A->numcols, A->numrows, lnnz, weighted);
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
      if(weighted) lgp_put_double(At->value, pos, A->lvalue[j]);
    }
  }

  lgp_barrier();
  //if(!MYTHREAD)printf("done\n");
  lgp_all_free(shtmp);

  return(At);
}

/*! \brief Write file called rowcnt_[thread number] into the given directory.
 * \param dirname The directory name
 * \param A  pointer to the matrix
 * \return 0 on success, -1 on failure
 * \ingroup spmatgrp
 */
int64_t write_rowcounts(char * dirname, sparsemat_t * A){
  int64_t i;
  
  SHARED int64_t * rowcnt = lgp_all_alloc(A->numrows, sizeof(int64_t));
  int64_t * lrowcnt = lgp_local_part(int64_t, rowcnt);

  for(i = 0; i < A->lnumrows; i++)
    lrowcnt[i] = A->loffset[i + 1] - A->loffset[i];
  lgp_barrier();

  char fname[64];
  sprintf(fname, "%s/rowcnt_%d", dirname, MYTHREAD);
  FILE * fp = fopen(fname, "wb");

  int64_t global_first_row = 0;
  for(i = 0; i < MYTHREAD; i++)
    global_first_row += (A->numrows + THREADS - i - 1)/THREADS;
  
  for(i = 0; i < A->lnumrows; i++){
    int64_t rc = lgp_get_int64(rowcnt, global_first_row + i);
    if( fwrite(&rc, sizeof(int64_t), 1,  fp) != 1 )
      return(-1);
  }
  fclose(fp);
  lgp_barrier();
  lgp_all_free(rowcnt);
  return(0);
}

/*! \brief Write a sparsemat_t to disk (in the directory dirname).
 * \param dirname The directory to write the sparse matrix data into.
 * \param A  pointer to the matrix
 * \return 0 on success, -1 on error
 * \ingroup spmatgrp
 */
int64_t write_sparse_matrix_agp(char * dirname, sparsemat_t * A){
  int64_t i;
  int values = 0;
  if(A->value)
    values = 1;
  
  /* create the directory  */
  mkdir(dirname, 0770);
  
  /* write metadata ASCII file */
  write_sparse_matrix_metadata(dirname, A);
  
  /* write out the row counts */
  write_rowcounts(dirname, A);
  
  /* figure out max_row_cnt */
  int64_t global_first_row = 0;
  for(i = 0; i < MYTHREAD; i++)
    global_first_row += (A->numrows + THREADS - i - 1)/THREADS;

  int64_t max_row_cnt = 0;
  for(i = 0; i < A->lnumrows; i++)
    max_row_cnt = (max_row_cnt < (A->loffset[i+1] - A->loffset[i]) ? A->loffset[i+1] - A->loffset[i] : max_row_cnt);
  max_row_cnt = lgp_reduce_max_l(max_row_cnt);
  
  /* allocate write buffer */
  int64_t BUFSIZE = 128;
  while(BUFSIZE < max_row_cnt)
    BUFSIZE *= 2;
  int64_t * buf = calloc(BUFSIZE, sizeof(int64_t));
  double * vbuf;
  if(values)
    vbuf = calloc(BUFSIZE, sizeof(double));

  char fname[64];
  sprintf(fname, "%s/nonzero_%d", dirname, MYTHREAD);
  FILE * fp = fopen(fname, "wb");
  FILE * vfp;
  if(values){
    sprintf(fname, "%s/value_%d", dirname, MYTHREAD);
    vfp = fopen(fname, "wb");
  }
  
  /* write out the nonzeros in your block */
  int64_t pos = 0;
  for(i = global_first_row; i < global_first_row + A->lnumrows; i++){
    int64_t row_start = lgp_get_int64(A->offset, i);
    int64_t row_cnt = lgp_get_int64(A->offset, i + THREADS) - row_start;
    if(pos + row_cnt >= BUFSIZE){
      fwrite(buf, sizeof(int64_t), pos, fp); 
      if(values) fwrite(vbuf, sizeof(double), pos, vfp);
      pos = 0;
    }
    int64_t plc = row_start*THREADS + i % THREADS;
    lgp_memget(&buf[pos], A->nonzero, row_cnt*sizeof(int64_t), plc);
    if(values) lgp_memget(&vbuf[pos], A->value, row_cnt*sizeof(int64_t), plc);
    pos += row_cnt;
  }
  fwrite(buf, sizeof(int64_t), pos, fp);
  if(values) fwrite(vbuf, sizeof(int64_t), pos, vfp);

  lgp_barrier();

  fclose(fp);  
  free(buf);

  if(values){
    fclose(vfp);
    free(vbuf);
  }
  
  return(0);  
}

sparsemat_t * read_sparse_matrix_agp(char * dirname){
  int64_t i;
  int64_t nr, nc, nnz, nwriters;
  int64_t values;
  char fname[128];
  
  /* read the metadata file */
  int ret = read_sparse_matrix_metadata(dirname, &nr, &nc, &nnz, &nwriters, &values);
  if(ret){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_agp: read_metadata failed!\n"); return(NULL);
  }

  //T0_fprintf(stderr,"metadata: %ld %ld %ld %ld\n", nc, nc, nnz, nwriters);fflush(stderr);
  
  int64_t rows_to_read = (nr + THREADS - MYTHREAD - 1)/THREADS;
  int64_t first_row_to_read = MYTHREAD * rows_to_read + ((MYTHREAD < nr % THREADS ? 0 : nr % THREADS));
  
  /* read rowcnts in BLOCKED fashion */
  int64_t * read_rc = read_rowcnts(dirname, nr, nwriters);
  if(read_rc == NULL){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_agp: read_rowcnts_failed");
    return(NULL);
  }

  /* distribute row counts in CYCLIC fashion */
  SHARED int64_t * rowcnt = lgp_all_alloc(nr, sizeof(int64_t));
  if(rowcnt == NULL){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_agp: could not allocate rowcnt");
    return(NULL);
  }
  int64_t * lrowcnt = lgp_local_part(int64_t, rowcnt);

  int64_t nnz_to_read = 0;
  for(i = 0; i < rows_to_read; i++){
    lgp_put_int64(rowcnt, first_row_to_read + i, read_rc[i]);
    nnz_to_read += read_rc[i];
  }
  lgp_barrier();
  
  /* calculate lnnz on each PE */
  int64_t lnnz = 0;
  for(i = 0; i < rows_to_read; i++){
    lnnz += lrowcnt[i];
  }
  //fprintf(stderr,"PE %d gets %ld nonzeros\n", MYTHREAD, lnnz);fflush(stderr);
  
  /* initialize sparse matrix */
  sparsemat_t * A = init_matrix(nr, nc, lnnz, values);
  if(!A){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_agp: init_matrix failed!\n");return(NULL);
  }

  /* create A->offset array */
  A->loffset[0] = 0;
  for(i = 0; i < rows_to_read; i++){
    A->loffset[i+1] = A->loffset[i] + lrowcnt[i];
  }

  lgp_barrier();
  lgp_all_free(rowcnt);
  
  /* populate nnz per file */
  SHARED int64_t * nnz_in_file = lgp_all_alloc(nwriters, sizeof(int64_t));
  
  for(i = MYTHREAD; i < nwriters; i+=THREADS){
    sprintf(fname, "%s/nonzero_%d", dirname, i);
    int fd = open(fname, O_RDONLY);
    struct stat buf;
    fstat(fd, &buf);
    lgp_put_int64(nnz_in_file, i, buf.st_size/8);
    //fprintf(stderr,"%ld nnz in file %ld\n", lgp_get_int64(nnz_in_file, i), i);
    close(fd);
  }

  lgp_barrier();
  
  /* find the right spot to seek to in the nonzeros files */
  int64_t skip = lgp_prior_add_l(nnz_to_read);

  int64_t current_file = 0;
  int64_t current_nnz = 0;
  for(current_file = 0; current_file < nwriters; current_file++){
    int64_t nnz_this_file = lgp_get_int64(nnz_in_file, current_file);
    if( current_nnz + nnz_this_file > skip)
      break;
    current_nnz += nnz_this_file;
  }
  
  /* open first file and seek to the right place */
  FILE * fp, *vfp;  
  sprintf(fname, "%s/nonzero_%d", dirname, current_file);
  fp = fopen(fname,"rb");
  if(values){
    sprintf(fname, "%s/value_%d", dirname, current_file);
    vfp = fopen(fname,"rb");
  }

  ret = 0;
  if(current_nnz < skip){
    ret = fseek(fp, skip - current_nnz, SEEK_SET);
    if(ret){
      fprintf(stderr,"ERROR: read_sparse_matrx: seek failed!\n");
    }
    if(values)  fseek(vfp, skip - current_nnz, SEEK_SET);
    current_nnz = skip;
  }
  ret = lgp_reduce_add_l(ret);
  if(ret) return(NULL);
  
  /* read your slice of nonzeros and spray them out to other PEs */
  int64_t current_local_row = 0;
  int new_file = 0;
  int64_t buf_size = 512*512;
  int64_t * buf = calloc(buf_size, sizeof(int64_t));
  double * vbuf;
  if(values) vbuf = calloc(buf_size, sizeof(double));

  while(nnz_to_read){

    if(new_file){
      sprintf(fname, "%s/nonzero_%d", dirname, current_file);
      fp = fopen(fname,"rb");
      if(values){
        sprintf(fname, "%s/value_%d", dirname, current_file);
        vfp = fopen(fname, "rb");
      }
    }
    new_file = 0;
    
    /* read from current file until EOF or you don't need to read any more */
    int64_t num_to_read_now = 0;
    int64_t num_rows_to_read_now = 0;
    while(1){
      num_to_read_now += read_rc[current_local_row + num_rows_to_read_now];
      num_rows_to_read_now++;
      if(num_to_read_now + read_rc[current_local_row + num_rows_to_read_now] >= buf_size)
        break;
      if(current_local_row + num_rows_to_read_now == rows_to_read)
        break;
    }
    
    size_t num_read = fread(buf, sizeof(int64_t), num_to_read_now, fp);
    //fprintf(stderr,"Trying to read %ld from %ld read %ld\n",
    //num_to_read_now, current_file, num_read);
    if(values){
      int64_t num_v_read = fread(vbuf, sizeof(double), num_to_read_now, vfp);
      assert(num_v_read == num_read);
    }
    nnz_to_read -= num_read;
    
    if(num_read != num_to_read_now || nnz_to_read == 0){
      fclose(fp);
      current_file++;
      new_file = 1; 
    }
    
    int64_t pos = 0;
    while(num_read){
      int64_t rc = read_rc[current_local_row];
      int64_t global_row = first_row_to_read + current_local_row;
      int64_t rs = lgp_get_int64(A->offset,global_row)*THREADS + global_row%THREADS;
      lgp_memput(A->nonzero, &buf[pos], rc*sizeof(int64_t), rs);
      if(values)
        lgp_memput(A->value, &vbuf[pos], rc*sizeof(double), rs);
      pos += rc;
      num_read -= rc;
      assert(num_read >= 0);
      current_local_row++;
    }
  }

  lgp_barrier();

  free(buf);
  if(values) free(vbuf);
  free(read_rc);
  lgp_all_free(nnz_in_file);

  //if(!MYTHREAD)
  //print_matrix(A);
  return(A);
}

