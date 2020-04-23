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
/*! \file spmat_utils.upc
 * \brief Utilities to support spmat and wrappers to make it easier to switch 
 * between routines written in the various models.
 */
#include <math.h>
#include <spmat.h>
#include <exstack.h>

/*! \brief Produce a global array the holds a uniform random permutation.
 * \param N the global length of the permutaion
 * \param seed the seed for the random number generator
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework. 
 *
 * This is a collective call.
 *
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp(int64_t N, int seed) {
  SHARED int64_t * p;
  //p = rand_permp_agi(N, seed);
  p = rand_permp_exstack(N, seed, 1024);
  //p = rand_permp_exstack2(N, seed, 1024);
  
  if(!is_perm(p, N)){
    T0_printf("ERROR: rand_permp: not a permutation!\n");fflush(0);
    return(NULL);
  }
  return(p);
}


/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param omat pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \param mode which buffering model to use 
 * \return a pointer to the matrix that has be computed or NULL on failure
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework. 
 *
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix(sparsemat_t *omat, SHARED int64_t *rperminv, SHARED int64_t *cperminv) {
  //return( permute_matrix_agi(omat, rperminv, cperminv) );
    return( permute_matrix_exstack(omat, rperminv, cperminv, 1024) );
  //return( permute_matrix_exstack2(omat, rperminv, cperminv, 1024) );
  //return( permute_matrix_conveyor(omat, rperminv, cperminv) );
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param omat  pointer to the original matrix
 * \return a pointer to the matrix that has be computed or NULL on failure
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it
 * within the apps framework. 
 *
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix(sparsemat_t *omat) {
  sparsemat_t * A;
  //A = transpose_matrix_agi(omat);
  A = transpose_matrix_exstack(omat, 1024);
  //A = transpose_matrix_exstack2(omat, 1024);
  //A = transpose_matrix_conveyor(omat);
  if(!A){return(NULL);}
  
  sort_nonzeros(A);
  return(A);
}


/*! \brief writes a sparse matrix to a file in Matrix Market format
 * \param A pointer to the sparse matrix
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 * \ingroup spmatgrp
 */
int write_matrix_mm(sparsemat_t *A, char * name) {
  if(!MYTHREAD){
    FILE * fp = fopen(name, "w");
    /* write the banner */
    fprintf(fp,"%%%%MatrixMarket matrix coordinate pattern\n");
    fprintf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", A->numrows, A->numcols, A->nnz);
    fclose(fp);
  }
  
  int64_t i, j,k, row;
  for(k = 0; k < THREADS; k++){
    if(k == MYTHREAD){
      FILE * fp = fopen(name, "a");
      for(i = 0; i < A->lnumrows; i++){
        row = i*THREADS + MYTHREAD;
        for(j = A->loffset[i]; j < A->loffset[i+1]; j++){
          fprintf(fp, "%"PRId64" %"PRId64"\n", row + 1, A->lnonzero[j] + 1);
        }
      }
      fclose(fp);
    }
    lgp_barrier();
  }
  
  return(0);
}


//
// this is new code not ready for release in bale2.0
//
int64_t read_sparse_matrix_metadata(char * dirname, int64_t * nr, int64_t * nc, int64_t * nnz, int64_t * nwriters){
  *nr = *nc = *nnz = 0;
  if(MYTHREAD == 0){
    int ret = 0;
    char fname[64];
    sprintf(fname, "%s/metadata", dirname);
    FILE * fp = fopen(fname, "r");
    ret += fscanf(fp, "%"PRId64"\n", nr);
    ret += fscanf(fp, "%"PRId64"\n", nc);
    ret += fscanf(fp, "%"PRId64"\n", nnz);
    ret += fscanf(fp, "%"PRId64"\n", nwriters);
    if(ret != 4){
      fprintf(stderr,"ERROR: read_sparse_matrix_metadata\n");      
    }
    *nr = -1;
    fclose(fp);
  }
  lgp_barrier();
  *nr = lgp_reduce_add_l(*nr);
  *nc = lgp_reduce_add_l(*nc);
  *nnz = lgp_reduce_add_l(*nnz);
  *nwriters = lgp_reduce_add_l(*nwriters);
  if(*nr < 0) return(-1);
  return(0);
}

/*! \brief This function writes out the ASCII metadata file for a sparse matrix dataset.
 * The metadata file is called "metadata" and contains 4 lines
 * - numrows
 * - numcols
 * - number of nonzeros
 * - number of PEs that were writing.
 *
 * \param dirname The name of the directory where the dataset will be written.
 * \param A The sparse matrix that is being written.
 * \ingroup spmatgrp
 */
int64_t write_sparse_matrix_metadata(char * dirname, sparsemat_t * A){
  if(MYTHREAD == 0){
    /* open params file for writing */
    char fname[64];
    sprintf(fname, "%s/metadata", dirname);
    FILE * fp = fopen(fname, "w");
    fprintf(fp, "%"PRId64"\n%"PRId64"\n%"PRId64"\n%d\n", A->numrows, A->numcols, A->nnz,THREADS);
    fclose(fp);
  }
  return(0);
}


/*! \brief Read a sparse matrix in matrix market format on one PE and create a distributed matrix
  from that.
  * Only PE 0 reads the matrix file.
  * 
  * \param name The name of the file.
  * \return The sparsemat_t struct.
  * \ingroup spmatgrp
  */
sparsemat_t * read_matrix_mm_to_dist(char * name) {
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  int64_t nr, nc, nnz = 0, i, pe;
  SHARED int64_t * sh_data;
  sh_data = lgp_all_alloc (THREADS, sizeof(int64_t));

  int64_t * Is;
  int64_t * rowcount;
  int64_t * Js;
  
  while(!MYTHREAD){

    int64_t * nnz_per_th = calloc(THREADS, sizeof(int64_t));

    FILE * fp = fopen(name, "r");
    
    /* read the banner */
    char * object = calloc(64, sizeof(char));
    char * format = calloc(64, sizeof(char));
    char * field = calloc(64, sizeof(char));;
    int ret = fscanf(fp,"%%%%MatrixMarket %s %s %s\n", object, format, field);
    if(ret != 3){
      T0_printf("ERROR: read_matrix_mm: bad read of header!\n");
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
      break;
    }
    if(strcmp(object, "matrix") || strcmp(format,"coordinate")){
      T0_printf("ERROR: read_matrix_mm: unsupported matrix market format! %s %s\n", object, format);
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
      break;
    }
    if(strcmp(field, "integer") && strcmp(field,"pattern")){
      T0_printf("ERROR: read_matrix_mm: unsupported matrix market format %s!\n", field);
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
    }
    
    fscanf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", &nr, &nc, &nnz);
    Is = calloc(nnz, sizeof(int64_t));
    Js = calloc(nnz, sizeof(int64_t));
    rowcount = calloc(nr, sizeof(int64_t));
    if(!rowcount || !Is || !Js){
      T0_printf("ERROR: read_matrix_mm_to_dist: could not allocate arrays\n");
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
      break;
    }
    int64_t row, col, val, pos = 0;
    if(!strcmp(field, "pattern")){
      while(fscanf(fp,"%"PRId64" %"PRId64"\n", &row, &col) != EOF){
        row--;
        col--;
        Is[pos] = row;
        Js[pos++] = col;
        nnz_per_th[row % THREADS]++;
        rowcount[row]++;
      }
    }else{
      while(fscanf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", &row, &col, &val) != EOF){
        row--;
        col--;
        Is[pos] = row;
        Js[pos++] = col;
        nnz_per_th[row % THREADS]++;
        rowcount[row]++;
      }
    }
    
    fclose(fp);
    if(nnz != pos){
      T0_printf("ERROR: read_matrix_mm_to_dist: nnz (%"PRId64") != pos (%"PRId64")\n", nnz, pos);
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
      break;
    }
    for(i = 0; i < THREADS; i++){
      lgp_put_int64(sh_data, i, nnz_per_th[i]);
      lgp_put_int64(sh_data, i+THREADS, nr);
      lgp_put_int64(sh_data, i+2*THREADS, nc);
    }
    free(nnz_per_th);
    break;
  }
  
  lgp_barrier();

  int64_t * lsh_data = lgp_local_part(int64_t, sh_data);
  if(lsh_data[0] == -1)
    return(NULL);
  
  nr = lsh_data[1];
  nc = lsh_data[2];
  sparsemat_t * A = init_matrix(nr, nc, sh_data[MYTHREAD]);
  SHARED int64_t * tmp_offset = lgp_all_alloc(nr + THREADS, sizeof(int64_t));
  
  if(!A || !tmp_offset){
    T0_printf("ERROR: read_matrix_mm_to_dist: failed to init matrix or tmp_offset!\n");
    return(NULL);
  }

  /* set up offset array and tmp_offset */
  lgp_barrier();
  lgp_all_free(sh_data);
  
  if(!MYTHREAD){
    for(i = 0; i < nr; i++)
      lgp_put_int64(tmp_offset, i+THREADS, rowcount[i]);    
    free(rowcount);
  }

  lgp_barrier();

  int64_t * ltmp_offset = lgp_local_part(int64_t, tmp_offset);
  for(i = 1; i <= A->lnumrows; i++){
    A->loffset[i] = ltmp_offset[i] += ltmp_offset[i-1];
  }

  int64_t fromth;
  pkg_rowcol_t pkg;
  exstack_t * ex = exstack_init(256, sizeof(pkg_rowcol_t));
  if( ex == NULL ) return(NULL);

  /* pass around the nonzeros */
  /* this is a strange exstack loop since only PE0 has data to push */
  i = 0;
  while(exstack_proceed(ex, (i == nnz))){
    while(i < nnz){
      pkg.row = Is[i];
      pkg.col = Js[i];
      pe = Is[i] % THREADS;
      if(!exstack_push(ex, &pkg, pe))
        break;
      i++;
    }
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg, &fromth)){
      A->lnonzero[ltmp_offset[pkg.row/THREADS]++] = pkg.col;
    }
  }

  lgp_barrier();
  if(!MYTHREAD){
    free(Is);
    free(Js);
  }
  lgp_all_free(tmp_offset);
  exstack_clear(ex);

  return(A);
}

/*! \brief checks that the sparse matrix is lower triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 * kind of a speciality routine to check that toposort might of worked
 * \ingroup spmatgrp
 */
int is_lower_triangular(sparsemat_t *A, int64_t unit_diagonal) {
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
  err2 = (unit_diagonal ? lgp_reduce_add_l(err2) : 0);
  if( err || err2 ){
    //if(!MYTHREAD)printf("\nThere are %"PRId64" nz above diag. and %"PRId64" missing pivots in lower.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief checks that the sparse matrix is upper triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 *
 * \ingroup spmatgrp
 */
int is_upper_triangular(sparsemat_t *A, int64_t unit_diagonal) {
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
  err2 = (unit_diagonal ? lgp_reduce_add_l(err2) : 0);
  if( err || err2 ){
    //if(!MYTHREAD)printf("\nThere are %"PRId64" nz below diag. and %"PRId64" missing pivots in upper.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief Sets all entries above the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param k Anything above the kth diagonal will be zero (k > 0 is above the main diagonal k < 0 is below)
 * \return 0 for success nonzero if error.
 * \ingroup spmatgrp
 */
int64_t tril(sparsemat_t * A, int64_t k) {
  // remove entries below the diagonal
  int64_t i, j, col, pos = 0, start = 0;
  for(i = 0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j = start; j < A->loffset[i+1]; j++){
      col = A->lnonzero[j];
      if(col - global_row <= k){
        A->lnonzero[pos++] = col;
      }
    }
    start = A->loffset[i+1];
    A->loffset[i+1] = pos;
  }
  A->nnz = lgp_reduce_add_l(pos);
  A->lnnz = pos;
  return(0);
}

/*! \brief Sets all entries below the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param k Anything below the kth diagonal will be zero (k > 0 is above the main diagonal k < 0 is below)
 * \return 0 for success nonzero if error.
 * \ingroup spmatgrp
 */
int64_t triu(sparsemat_t * A, int64_t k) {
  // remove entries below the diagonal
  int64_t i, j, col, pos = 0, start = 0;
  for(i = 0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j = start; j < A->loffset[i+1]; j++){
      col = A->lnonzero[j];
      if(col - global_row >= k){
        A->lnonzero[pos++] = col;
      }
    }
    start = A->loffset[i+1];
    A->loffset[i+1] = pos;
  }
  A->nnz = lgp_reduce_add_l(pos);
  A->lnnz = pos;
  return(0);
}


/*! \brief checks that a global array is in fact a permutation
 * \param perm SHARED pointer to the global array
 * \param N the length of the global array
 * \return 1 if it is permutation
 * \ingroup spmatgrp
 */
int is_perm(SHARED int64_t * perm, int64_t N) {
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
int nz_comp(const void *a, const void *b) {
  return( *(uint64_t *)a - *(uint64_t *)b );
}


/*! \brief Generates the upper or lower half of the adjacency matrix (non-local) for an Erdos-Renyi random
 * graph. This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
 * this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the 
 * adjancecy matrix using a geometric random variable.
 *
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param unit_diag 1 - set the diagonal to all ones, 0's otherwise
 * \param mode 0 : return a symmetric adjacency matrix
               1 : return the lower-triangular portion of the adjacency matrix 
               2 : return the upper-triangular portion of the adjacency matrix
               3 : return an asymmetric random matrix.
 * \param seed A random seed.
 * \return A distributed sparsemat_t 
 */
sparsemat_t * gen_erdos_renyi_graph_dist(int n, double p, int64_t unit_diag, int64_t mode, int64_t seed) {
  //T0_fprintf(stderr,"Entering gen_erdos_renyi_graph_dist...\n");

  sparsemat_t * L, * U;
  switch(mode){
  case 0:
    /* generate the upper triangular portion, then transpose and add */
    U = gen_erdos_renyi_graph_triangle_dist(n, p, unit_diag, 0, seed);
    L = transpose_matrix(U);
    if(!L){T0_printf("ERROR: gen_er_graph_dist: L is NULL!\n"); return(NULL);}
    break;
  case 1:
    return(gen_erdos_renyi_graph_triangle_dist(n, p, unit_diag, 1, seed));
  case 2:
    return(gen_erdos_renyi_graph_triangle_dist(n, p, unit_diag, 0, seed));
  case 3:
    /* generate to separate halves and add together */
    U = gen_erdos_renyi_graph_triangle_dist(n, p, unit_diag, 0, seed);
    L = gen_erdos_renyi_graph_triangle_dist(n, p, unit_diag, 1, rand());    
  }

  int64_t i, j;
  int64_t lnnz = L->lnnz + U->lnnz;
  sparsemat_t * A2 = init_matrix(n, n, lnnz);
  if(!A2){T0_printf("ERROR: gen_er_graph_dist: A2 is NULL!\n"); return(NULL);}
  
  A2->loffset[0] = 0;
  lnnz = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    for(j = L->loffset[i]; j < L->loffset[i+1]; j++)
      A2->lnonzero[lnnz++] = L->lnonzero[j];
    for(j = U->loffset[i]; j < U->loffset[i+1]; j++){
      if(U->lnonzero[j] != global_row) // don't copy the diagonal element twice!
        A2->lnonzero[lnnz++] = U->lnonzero[j];
    }
    A2->loffset[i+1] = lnnz;
  }
  lgp_barrier();
  clear_matrix(U);
  clear_matrix(L);
  free(U); free(L);
  return(A2);
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
sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, int64_t seed){

  
  
}
/*! \brief Generates the lower half of the adjacency matrix (non-local) for an Erdos-Renyi random
 * graph. This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
 * this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the 
 * adjancecy matrix using a geometric random variable.
 *
 * We parallelized this algorithm by noting that the geometric random variable is memory-less. This means, we 
 * can start each PE independently (still need to think about the details of exactly how we did this below).
 *
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param edge_type See edge_type enum. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops enum. Does all or no vertices have self loops.
 * \param seed A random seed.
 * \return A distributed sparsemat_t
 */
sparsemat_t * erdos_renyi_random_graph(int n, double p, edge_type edge_type, self_loops loops, int64_t seed){
  
  int64_t row, col, i;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnnz = 0, lnnz_orig=0;
  int64_t P = p*RAND_MAX;
  double lM = log(RAND_MAX);
  double D = log(1 - p);
  int64_t r;
  int64_t end = n;
  
  /* count lnnz so we can allocate A correctly */
  srand(seed);
  row = MYTHREAD;
  do { r = rand(); } while(r == RAND_MAX);     
  col = 1 + floor((log(RAND_MAX - r) - lM)/D);
  while(row < n){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) ndiag--; // we just hit a diagonal entry (we don't have to generate this one later)
      lnnz_orig++;
      do { r = rand(); } while(r == RAND_MAX);     
      col += 1 + floor((log(RAND_MAX - r) - lM)/D);
    }
    row += THREADS;
    col -= end;
  }
  if(loops == LOOPS) lnnz_orig += ndiag;
  
  lgp_barrier();

  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, lnnz_orig, weighted);
  if(!A){T0_printf("ERROR: erdos_renyi_random_graph: init_matrix failed!\n"); return(NULL);}

  /* reset the seed so we get the same sequence of coin flips */
  srand(seed);
  A->loffset[0] = 0;
  row = MYTHREAD;
  do { r = rand(); } while(r == RAND_MAX);     
  col = 1 + floor((log(RAND_MAX - r) - lM)/D);
  while(row < n){
    int need_diag = (loops == LOOPS);
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) need_diag = 0;
      A->lnonzero[lnnz++] = col;
      do { r = rand(); } while(r == RAND_MAX);     
      col += 1 + floor((log(RAND_MAX - r) - lM)/D);
    }
    if(need_diag) A->nonzero[lnnz++] = row;
    row+=THREADS;
    col -= end;
    A->loffset[row] = lnnz;
  }
  
  if(lnnz != lnnz_orig){
    printf("ERROR: lnnz (%"PRId64") != lnnz_orig (%"PRId64")\n", lnnz, lnnz_orig);
    return(NULL);
  }

  // fill in the weights
  if(weighted){
    for(i = 0; i < lnnz; i++){
      A->lvalue[i] = (double)rand()/RAND_MAX;
    }
  }

  if(loops == LOOPS && (edge_type == DIRECTED || edge_type == DIRECTED_WEIGHTED))
    sort_nonzeros(A); // to get the diagonal entry sorted correctly
  
  return(A);
}

/*! \brief Generates the upper or lower half of the adjacency matrix (non-local) for an Erdos-Renyi random
 * graph. This is the naive O(n^2) algorithm. It flips a coin for each possible edge.
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param edge_type See edge_type enum. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops enum. Do all or no vertices have self loops.
 * \param seed A random seed.
 * \return A distributed sparsemat_t (the lower half of the adjacency matrix).
 */
sparsemat_t * erdos_renyi_random_graph_naive(int n, double p, edge_type edge_type, self_loops loops, int64_t seed){
  
  int64_t row, col, i, j;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnnz, lnnz_orig;
  int64_t P = p*RAND_MAX;
  int64_t end = n;
  
  /* count lnnz so we can allocate A correctly */
  srand(seed);
  lnnz_orig = 0;
  for(row = MYTHREAD; row < n; row += THREADS){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    for(j = 0; j < end; j++){
      if(col == row) continue;
      if(rand() < P)
	lnnz_orig++;
    }    
  }
  if(loops == LOOPS) lnnz_orig += n;
  
  lgp_barrier();

  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, lnnz_orig, weighted);
  if(!A){T0_printf("ERROR: erdos_renyi_random_graph_naive: init_matrix failed!\n"); return(NULL);}

  /* reset the seed so we get the same sequence of coin flips */
  srand(seed);
  A->loffset[0] = 0;  
  lnnz = 0;
  for(row = MYTHREAD; row < n; row += THREADS){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row + (loops == LOOPS);
    for(col = 0; col < row; col++){
      if(col == row && loops == LOOPS){
	A->nonzero[lnnz++] = row;
	continue;
      }
      if(rand() < P){
	A->lnonzero[lnnz++] = col;
      }      
    }
    A->loffset[row/THREADS + 1] = lnnz;
  }

  if(lnnz != lnnz_orig){
    printf("ERROR: erdos_renyi_random_graph_naive: lnnz (%"PRId64") != lnnz_orig (%"PRId64")\n", lnnz, lnnz_orig);
    return(NULL);
  }

  //fill in the weights
  if(weighted){
    int64_t i;
    for(i = 0; i < lnnz; i++){
      A->lvalue[i] = (double)rand()/RAND_MAX;
    }
  }
  
#if 0
  /* if the user wants the whole symmetric matrix, 
     we have made the upper so transpose and add */
  if(mode == 0){ 
    sparsemat_t * L = transpose_matrix(A);
    if(!L){T0_printf("ERROR: gen_er_graph_dist: L is NULL!\n"); return(NULL);}
    
    lnnz = L->lnnz + A->lnnz;
    sparsemat_t * A2 = init_matrix(n, n, lnnz);
    if(!A2){T0_printf("ERROR: gen_er_graph_dist: A2 is NULL!\n"); return(NULL);}

    A2->loffset[0] = 0;
    lnnz = 0;
    for(i = 0; i < L->lnumrows; i++){
      int64_t global_row = i*THREADS + MYTHREAD;
      for(j = L->loffset[i]; j < L->loffset[i+1]; j++)
        A2->lnonzero[lnnz++] = L->lnonzero[j];
      for(j = A->loffset[i]; j < A->loffset[i+1]; j++){
        if(A->lnonzero[j] != global_row) // don't copy the diagonal element twice!
          A2->lnonzero[lnnz++] = A->lnonzero[j];
      }
      A2->loffset[i+1] = lnnz;
    }
    lgp_barrier();
    clear_matrix(A);
    clear_matrix(L);
    free(A);
    free(L);
    A = A2;
  }
  #endif
  return(A);
  
}

/*! \brief Generate kron(B,C) in a distributed matrix (B and C are local matrices and all PEs have the same B and C).
 * \param B A local sparsemat_t (perhaps generated with kron_prod) that is the adjacency matrix of a graph.
 * \param C Another local sparsemat_t that is the adjacency matrix of a graph.
 * \param lower Flag to say whether you want the adjacency matrix to be lower triangular (1) or symmetric (0).
 * \return The adjacency matrix (or lower half of) the kronecker product of B and C.
 * \ingroup spmatgrp
 */
sparsemat_t * kron_prod_dist(sparsemat_t * B, sparsemat_t * C, int64_t lower) {

  int64_t lrow, row, rowB, rowC, col, pos, i, j, k;
  int64_t n = B->lnumrows*C->lnumrows;
  int64_t lnnz = 0;

  /* we first calculate the number of nonzeros on each thread */
  row = MYTHREAD;
  rowB = row / C->lnumrows;
  rowC = row % C->lnumrows;
  while(row < n){
    lnnz += (B->loffset[rowB + 1] - B->loffset[rowB])*(C->loffset[rowC + 1] - C->loffset[rowC]);
    row += THREADS;
    rowB = row / C->lnumrows;
    rowC = row % C->lnumrows;
  }
  
  sparsemat_t * A = init_matrix(n, n, lnnz);
  if(A == NULL){
    printf("ERROR: kron_prod_dist: init_matrix returned NULL with inputs %"PRId64" %"PRId64" %"PRId64"\n", n, n, lnnz);
    return(NULL);
  }
  row = MYTHREAD;
  lrow = 0;
  rowB = row / C->lnumrows;
  rowC = row % C->lnumrows;
  pos = 0;
  A->loffset[0] = 0;
  
  while(row < n){
    
    for(j = B->loffset[rowB]; j < B->loffset[rowB + 1]; j++){
      int64_t Bcol = B->lnonzero[j];
      for(k = C->loffset[rowC]; k < C->loffset[rowC+1]; k++){
        col = Bcol*C->numcols + C->lnonzero[k];
        if(!lower || col < row)
          A->lnonzero[pos++] = col;
        else
          break;
      }      
    }
    A->loffset[lrow+1] = pos;
    
    row += THREADS;
    lrow++;
    rowB = row / C->lnumrows;
    rowC = row % C->lnumrows;
  }
  //if(pos != lnnz){
  //printf("ERROR: kron_prod_dist: pos (%"PRId64") != lnnz (%"PRId64")\n", pos, lnnz);
  //return(NULL);
  //}
  lgp_barrier();
  A->lnnz = pos;
  A->nnz = lgp_reduce_add_l(pos);

#if 0
  SHARED int64_t * sh_data = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t * lsh_data = lgp_local_part(int64_t, sh_data);
  lsh_data[0] = 0;
  lgp_barrier();
  
  /* remove the self loop */
  if(!MYTHREAD){
    fprintf(stderr,"%"PRId64" %"PRId64" %"PRId64"\n", A->lnonzero[A->loffset[0]], A->loffset[0], A->loffset[1]);
    if(A->lnonzero[A->loffset[0]] == 0 && A->loffset[1] != A->loffset[0]){
      fprintf(stderr,"Removing A[0,0]!\n");
      /* set A(0,0) = 0*/
      for(i = 1; i < A->lnumrows+1; i++)
        A->loffset[i]--;
      for(i = 0; i < A->lnnz - 1; i++){
        A->lnonzero[i] = A->lnonzero[i+1];
      }
      A->lnnz--;      
      for(i = 0; i < THREADS; i++)
        lgp_atomic_add(sh_data, i, 1);
    }
  }
  if(MYTHREAD == ((A->numrows-1) % THREADS)){
    if(A->lnonzero[A->lnnz - 1] == (A->numrows-1)){
      fprintf(stderr,"Removing A[m,m]!\n");
      /* set A(m1*m2 - 1, m1*m2 - 1) = 0 */
      A->loffset[A->lnumrows]--;
      A->lnnz--;
      for(i = 0; i < THREADS; i++)
        lgp_atomic_add(sh_data, i, 1);
    }
  }
  lgp_barrier();
  if(lsh_data[0] == 1){
    A->nnz--;
  }
  if(lsh_data[0] > 1){
    T0_fprintf(stderr,"ERROR: kron_prod_dist: two diagonal elements!\n");
    return(NULL);
  }
#endif
  
  return(A);
}


/*! \brief Get the Kronecker product of two local matrices and store the result in a local matrix.
 * \param B A local sparse matrix (which is the adjacency matrix of a graph).
 * \param C A local sparse matrix (which is the adjacency matrix of a graph).
 * \return The adjacency matrix of the Kronecker product of B and C.
 * \ingroup spmatgrp
*/
sparsemat_t * kron_prod(sparsemat_t * B, sparsemat_t * C) {
  int64_t row, col, i, j, k;
  int64_t n = B->lnumrows*C->lnumrows;
  int64_t nnz = B->lnnz*C->lnnz;

  sparsemat_t * A = init_local_matrix(n, n, nnz);

  /* get the number of nonzeros in each row */
  int64_t * tmp_offset = calloc(n + 1, sizeof(int64_t));
  for(row = 0; row < B->lnumrows; row++){
    int64_t d1 = B->loffset[row + 1] - B->loffset[row];
    for(i = 0; i < C->lnumrows; i++){
      int64_t d2 = C->loffset[i + 1] - C->loffset[i]; 
      tmp_offset[1 + row*C->lnumrows + i] = d1*d2;
    }
  }

  /* get the cumulative sum */
  for(i = 0; i < n; i++){
    tmp_offset[i + 1] = tmp_offset[i + 1] + tmp_offset[i];
    A->loffset[i+1] = tmp_offset[i+1];
  }
  
  for(row = 0; row < B->lnumrows; row++){
    for(j = B->loffset[row]; j < B->loffset[row+1]; j++){
      col = B->lnonzero[j];
      for(i = 0; i < C->lnumrows; i++){
        for(k = C->loffset[i]; k < C->loffset[i+1]; k++){
          A->lnonzero[tmp_offset[row*C->lnumrows + i]++] = col*C->numcols + C->lnonzero[k];
        }
      }
    }
  }
  
  free(tmp_offset);

  return(A);
}

/*! \brief Generate the adjacency matrix for K_{1,m} (with 0 or 1 loop edges) 
 * \param m The parameter m in K_{1,m} (this is the number of non-hub vertices).
 * \param mode
 *  - mode == 0: default
 *  - mode == 1: add a self loop to center vertex of star
 *  - mode == 2: add a self loop to an outer vertex
 * \return The adjacency matrix of the graph.
 * \ingroup spmatgrp
*/
sparsemat_t * gen_star(int64_t m, int mode) {
  sparsemat_t * G = init_local_matrix(m + 1, m + 1, 2*m + (mode > 0 ? 1 : 0));
  
  int64_t i, pos = 0;

  if(mode == 1) /* add a self loop to center vertex */    
    G->lnonzero[pos++] = 0;

  /* Add the top row */
  for(i = 0; i < m; i++){
    G->lnonzero[pos++] = i+1;
  }  
  G->loffset[1] = pos;

  /* rest of the rows */
  for(i = 1; i < m+1; i++){
    G->lnonzero[pos++] = 0;
    G->loffset[i+1] = pos;
  }
  
  if(mode == 2){/* add a self loop to last outer vertex */
    G->lnonzero[pos++] = m;
    G->loffset[m+1] = pos;
  }
  
  return(G);
}

/*! \brief Generate the kroncker product of a collection of star graphs.
 * \param M The number of elements in the list. M must be at least 2 
 * \param m The list of parameters: each element m[i] represents a K_{1,m[i]} graph.
 * \param mode Add 0 or 1 loop edges.
 * - mode == 0: default
 * - mode == 1: add a self loop to center vertex of star
 * - mode == 2: add a self loop to an outer vertex
 * \return The local adjacency matrix for the kronecker product of the list of star graphs.
 * \ingroup spmatgrp
*/
sparsemat_t * gen_local_mat_from_stars(int64_t M, int64_t * m, int mode) {
  int64_t i;

  if(M < 1){
    T0_printf("ERROR: gen_local_mat_from_stars: M < 1!\n");
    return(NULL);
  }

  if(M == 1){
    return(gen_star(m[0], mode));
  }
  
  sparsemat_t ** Aarr = calloc(2*M, sizeof(sparsemat_t *));
  for(i = 0; i < M; i++)
    Aarr[i] = gen_star(m[i], mode);

  Aarr[M] = kron_prod(Aarr[0], Aarr[1]);
  for(i = 0; i < M-2; i++)
    Aarr[M + 1 + i] = kron_prod(Aarr[M+i], Aarr[2+i]);
 
  sparsemat_t * A = Aarr[2*M-2];
  Aarr[2*M - 2] = NULL;

  for(i = 0; i < 2*M-2; i++){
    clear_matrix(Aarr[i]); free(Aarr[i]);
  }
#if 0
  /* remove the self loop */
  if(mode == 1){
    /* set A(0,0) = 0*/
    for(i = 0; i < A->lnumrows+1; i++)
      A->loffset[i]--;
    for(i = 0; i < A->lnnz - 1; i++){
      A->lnonzero[i] = A->lnonzero[i+1];
    }
  }else if(mode == 2){
    /* set A(m1*m2 - 1, m1*m2 - 1) = 0 */
    A->loffset[A->lnumrows] -= 1;
  }
  A->lnnz = A->nnz = A->nnz - 1;
#endif
  return(A);
}

/*! \brief sort the non-zeros in each row of a sparse matrix (make it tidy)
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
int sort_nonzeros( sparsemat_t *mat) {
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
int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat) {
  int i,j;

  if( lmat->numrows != rmat->numrows ){
    if(!MYTHREAD)printf("(lmat->numrows = %"PRId64")  != (rmat->numrows = %"PRId64")", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->lnumrows != rmat->lnumrows ){
    fprintf(stderr,"THREAD %03d: (lmat->lnumrows = %"PRId64")  != (rmat->lnumrows = %"PRId64")", 
            MYTHREAD, lmat->lnumrows, rmat->lnumrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    if(!MYTHREAD)printf("(lmat->numcols = %"PRId64")  != (rmat->numcols = %"PRId64")", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    if(!MYTHREAD)printf("(lmat->nnz = %"PRId64")  != (rmat->nnz = %"PRId64")", lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->lnnz != rmat->lnnz ){
    fprintf(stderr,"THREAD %03d: (lmat->lnnz = %"PRId64")  != (rmat->lnnz = %"PRId64")", 
            MYTHREAD, lmat->lnnz, rmat->lnnz );
    return(1);
  }

  if( lmat->loffset[0] != 0 || rmat->loffset[0] != 0 
    || (lmat->loffset[0] != rmat->loffset[0] ) ){
    if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[0] = %"PRId64")  != (rmat->loffset[0] = %"PRId64")", 
       MYTHREAD, lmat->loffset[0], rmat->loffset[0] );
    return(1);
  }

  
  for(i = 0; i < lmat->lnumrows; i++){
    if( lmat->loffset[i+1] != rmat->loffset[i+1] ){
       if(!MYTHREAD)printf("THREAD %03d: (lmat->loffset[%d] = %"PRId64")  != (rmat->loffset[%d] = %"PRId64")", 
          MYTHREAD, i+1, lmat->loffset[i+1], i+1, rmat->loffset[i+1] );
       return(1);
    }
  }
  
  for(j=0; j< lmat->lnnz; j++) {
    if( lmat->lnonzero[j] != rmat->lnonzero[j] ){
      if(!MYTHREAD)printf("THREAD %03d: (lmat->lnonzero[%d] = %"PRId64")  != (rmat->lnonzero[%d] = %"PRId64")", 
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
sparsemat_t * copy_matrix(sparsemat_t *srcmat) {
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
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread) {
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->local = 0;
  mat->numrows  = numrows;
  mat->lnumrows = (numrows + THREADS - MYTHREAD - 1)/THREADS;
  mat->numcols  = numcols;  
  mat->offset   = lgp_all_alloc(mat->numrows + THREADS, sizeof(int64_t));
  if(mat->offset == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %"PRId64" bytes for offset array\n", mat->numrows*8);
    return(NULL);
  }
  mat->loffset  =  lgp_local_part(int64_t, mat->offset);
  int64_t max = lgp_reduce_max_l(nnz_this_thread);
  int64_t total = lgp_reduce_add_l(nnz_this_thread);
  mat->nonzero = lgp_all_alloc(max*THREADS, sizeof(int64_t));
  if(mat->nonzero == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %"PRId64" bytes for nonzero array (max = %"PRId64")\n", max*THREADS*8, max);
    return(NULL);
  }
  mat->lnonzero = lgp_local_part(int64_t, mat->nonzero);
  mat->nnz = total;
  mat->lnnz = nnz_this_thread;

  return(mat);
}

/*! \brief initialize the local parts of a sparsemat_t
 * \param numrows the number of rows
 * \param numcols the number of cols
 * \param nnz the number of nonzeros
 * \ingroup spmatgrp
 */
sparsemat_t * init_local_matrix(int64_t numrows, int64_t numcols, int64_t nnz) {
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->local = 1;
  mat->numrows  = numrows;
  mat->lnumrows = numrows;
  mat->numcols  = numcols;  
  mat->offset   = NULL;
  mat->loffset  = calloc(mat->lnumrows + 1, sizeof(int64_t));
  if(mat->loffset == NULL) 
    return(NULL);
  
  mat->nonzero = NULL;
  mat->lnonzero = calloc(nnz, sizeof(int64_t));
  if(mat->lnonzero == NULL) 
    return(NULL);

  mat->nnz = nnz;
  mat->lnnz = nnz;

  return(mat);
}

/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
void clear_matrix(sparsemat_t * mat) {
  if(mat->local){
    free(mat->lnonzero);
    free(mat->loffset);
  }else{
    lgp_all_free(mat->nonzero);
    lgp_all_free(mat->offset);
  }
}

/*!
 * \brief returns the number of nonzeros in a row of the localize part of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param l_row row index into local version (loffset) of the offset array
 * \return the number of nonzero in the request row
 */
int64_t rowcount_l( sparsemat_t *mat, int64_t l_row ) {
   return( mat->loffset[l_row+1] - mat->loffset[l_row] );
}

/*!
 * \brief returns the number of nonzeros in a row of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param S_row global row index into the shared offset array
 * \return the number of nonzero in the request row
 */
int64_t rowcount_S( sparsemat_t *mat, int64_t S_row ) {
   return( mat->offset[S_row+THREADS] - mat->offset[S_row] );
}

/*!
 * \brief Initialize an nxnz_t struct 
 * \param *mat pointer to the sparse matrix 
 * \return a pointer to the nxnz_t struct associated with the matrix, rest of the state is set to meaningless values.
 */
nxnz_t * init_nxnz( sparsemat_t *mat ) {
  nxnz_t *ret = calloc(1,sizeof(nxnz_t));
  ret->mat = mat;
  ret->row = -1;
  ret->first = 0x7FFFFFFFFFFFFFFF;
  ret->idx = -1;
  ret->stop = -1;
  ret->col = -1;
  return(ret);
};

/*!
 * \brief Sets the state of an nxnz_t struct for local first touch of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return void
 */
void first_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row < nxz->mat->lnumrows );
  nxz->row   = l_row;
  nxz->first = nxz->mat->loffset[l_row];
  nxz->idx   = nxz->mat->loffset[l_row];
  nxz->stop  = nxz->mat->loffset[l_row+1];
  nxz->col   = nxz->mat->lnonzero[nxz->idx];
}

/*!
 * \brief The condition part of a for loop across the local view of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return true if we are at a valid index, false otherwise
 */
bool has_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row == nxz->row );
  if( nxz->idx < nxz->stop )
    return(true);
  return(false);
}

/*!
 * \brief Attempts to move to the next nonzero in the local view of the given row,
 *   by incrementing the index and setting the col field on the nxnz_t.
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return void
 */
void incr_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row == nxz->row );
  nxz->idx += 1;
  nxz->col   = nxz->mat->lnonzero[nxz->idx];
}

/*!
 * \brief Sets the state of an nxnz_t struct for first touch of the given row,
 *   from the shared view of the matrix
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return void
 */
void first_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row < nxz->mat->numrows );
  nxz->row   = S_row;
  nxz->first = lgp_get_int64(nxz->mat->offset, S_row);
  nxz->idx   = nxz->first;
  nxz->stop  = lgp_get_int64(nxz->mat->offset, S_row+THREADS);
  nxz->col   = lgp_get_int64(nxz->mat->nonzero, (nxz->idx) * THREADS +  (nxz->row)%THREADS );
}

/*!
 * \brief The condition part of a for loop across the shared view of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return true if we are at a valid index, false otherwise
 */
bool has_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row == nxz->row );
  if( nxz->idx < nxz->stop )
    return(true);
  return(false);
}

/*!
 * \brief Attempts to move to the next nonzero in the shared view of the given row,
 *   by incrementing the index and setting the col field of the nxnz_t.
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return void
 */
void incr_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row == nxz->row );
  nxz->idx += 1; 
  nxz->col   = lgp_get_int64(nxz->mat->nonzero, (nxz->idx) * THREADS +  (nxz->row)%THREADS );
}

