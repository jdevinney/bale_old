/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// 
//
//  All rights reserved.
//  
//   This file is a part of Bale.  For license information see the
//   LICENSE file in the top level directory of the distribution.
//  
// 
 *****************************************************************/ 

/*! \file triangle_exstack2.upc
 * \brief Demo application that does a triangle on a permuted upper triangular matrix
 * using the exstack2 buffering model. The implementation of triangle_exstack2_push
 * and its helper routine tri_exstack2_push_process are nearly identical to 
 * triangle_exstack_push and tri_exstack_push_process.
 */
#include "triangle.h"

typedef struct pkg_tri_t {
  int64_t w;    
  int64_t vj;
}pkg_tri_t;

/*!
 * \brief pop routine to handle the pushes
 * \param *c a place to return the number of hits.
 * \param *ex2 the extack2 buffers
 * \param *mat the input sparse matrix 
     NB. The nonzero within a row must be increasing
 * \param done the signal to exstack2_proceed that this thread is done
 * \return the return value from exstack2_proceed
 */
static int64_t 
tri_exstack2_push_process(int64_t *c, exstack2_t *ex2, sparsemat_t *mat, int64_t done) {
  int64_t k, cnt = 0;
  struct pkg_tri_t pkg;
    
  while( exstack2_pop(ex2, &pkg, NULL) ){
    for(k = mat->loffset[pkg.vj]; k < mat->loffset[pkg.vj + 1]; k++){
       if( pkg.w == mat->lnonzero[k]){
         cnt ++;
         break;
       }
       if( pkg.w < mat->lnonzero[k]) // requires that nonzeros are increasing
         break;
    }
  }
  *c += cnt;

  return( exstack2_proceed(ex2, done) );
}


/*!
 * \brief This routine implements the exstack2 variant of triangle counting.
 * NB The column indices of the nonzeros in each row must be in increasing order.
 * 
 * \param *count a place to return the number hits
 * \param *sr a place to return the number of pushes 
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \param buf_cnt the size of the exstack2 buffers
 * \return average run time
 */
double triangle_exstack2_push(int64_t * count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg, int64_t buf_cnt) {

  exstack2_t * ex2 = exstack2_init(buf_cnt, sizeof(pkg_tri_t));
  if( ex2 == NULL ) return(-1.0);

  int64_t cnt = 0;
  int64_t numpushed = 0;
  
  int64_t k, kk, pe;
  int64_t l_i, L_i, L_j;

  pkg_tri_t pkg;

  double t1 = wall_seconds();
  if(alg == 0){
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_i > L_j );
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
          pkg.w = L->lnonzero[kk]; 
          if( pkg.w > L_j) 
            break;
          numpushed++;
          if(!exstack2_push(ex2, &pkg, pe)){
            tri_exstack2_push_process(&cnt, ex2, L, 0); 
            numpushed--;
            kk--;
          }
        }
      }
    }
    while ( tri_exstack2_push_process(&cnt, ex2, L, 1) )
      ;
  }else{
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_i > L_j );
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
          pkg.w = U->lnonzero[kk]; 
          numpushed++;
          if(!exstack2_push(ex2, &pkg, pe)){
            tri_exstack2_push_process(&cnt, ex2, U, 0); 
            numpushed--;
            kk--;
          }
        }
      }
    }
    while ( tri_exstack2_push_process(&cnt, ex2, U, 1) )
      ;
  }

  lgp_barrier();
  *sr = numpushed;
  *count = cnt;
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  return(stat->avg);
}

#if 0
static int64_t tri_exstack2_pull_resp_drain(exstack2_t * ex_resp, sparsemat_t * mat, int64_t done){
  while(exstack_pop(ex_resp, &pkg, &fromth)){
    for(ii = mat->loffset[pkg.w]; ii < mat->loffset[pkg.w+1]; ii++){
      if(mat->lnonzero[ii] == pkg.vj)
        cnt++;
      else if(mat->lnonzero[ii] > pkg.vj)
        break;
    }
  }
  exstack2_proceed(ex_req, done)
}

static int64_t tri_exstack2_pull_req_drain(exstack2_t * ex_req, exstack2_t * ex_resp, sparsemat_t * mat, int64_t done){
  int mid_request = 0;
  int done_popping = 0;
  int64_t row, label, requstr;
    
  while(!done_popping){
    if(!mid_request){
      ret = exstack_pop(ex_req, &pkg, &fromth);
      if(ret == 0){
        done_popping = 1;
      }
      else{
        //PE 'fromth' wants me to send him all columns in row 'pkg.vj' each with label 'pkg.w'
        row = pkg.vj;
        label = pkg.w;
        requstr = fromth;
        kk = mat->loffset[row];
      }
    }
      
    if(!done_popping && kk < mat->loffset[row+1]){
      pkg.w = label;
      pkg.vj = mat->lnonzero[kk++];
      mid_request = (kk == mat->loffset[row + 1]) ? 0 : 1;
      exstack_pushes++;
      if(exstack_push(ex_resp, &pkg, requstr) == 1)
        tri_exstack2_pull_resp_drain(ex_resp, mat, 0);
    }
  }
  
  return(exstack2_proceed(ex_resp, done));
  
}

/*!
 * \brief This routine implements the exstack2 variant of triangle counting.
 * NB The column indices of the nonzeros in each row must be in increasing order.
 * 
 * \param *count a place to return the number hits
 * \param *sr a place to return the number of pushes 
 * \param *mat lower triangle of the input matrix
 * \param buf_cnt the size of the exstack2 buffers
 * \return average run time
 */
double triangle_exstack2_pull(int64_t *count, int64_t *sr, sparsemat_t *mat, int64_t buf_cnt){

  exstack_t * ex_req = exstack2_init(buf_cnt, sizeof(pkg_tri_t));
  exstack_t * ex_resp = exstack2_init(buf_cnt, sizeof(pkg_tri_t));
  if( !ex_req || !ex_resp ){return(-1.0);}

  int64_t exstack_pushes = 0;
  int64_t l_i, L_i, L_j, k, kk, ii, pe, fromth;
  pkg_tri_t pkg;
  int64_t cnt = 0;
  double t1 = wall_seconds();
  int ret;
  
  // foreach nonzero in L
  l_i = 0;
  L_i = MYTHREAD;
  k = 0;
  while(k == mat->loffset[l_i + 1]){
    l_i++;
    L_i += THREADS;
  }
  
  //while(exstack2_proceed(ex_req, (k == mat->lnnz)) ){

  while(k < mat->lnnz){
    L_j = mat->lnonzero[k];       // shared name for col j
    assert( L_i > L_j );
    pkg.w = l_i;
    pkg.vj = L_j/THREADS;      
    pe = L_j % THREADS;
    ret = exstack_push(ex_req, &pkg, pe);
    exstack_pushes++;
      
    k++;
    while(k == mat->loffset[l_i + 1]){
      l_i++;
      L_i += THREADS;
    }

    if(ret == 1)
      break;

    tri_exstack2_pull_req_drain(ex_req, ex_resp, mat, 0);
  }

  while(tri_exstack2_req_drain(ex_req, ex_resp, mat, 1));
  
  

  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  *sr = exstack_pushes;
  *count = cnt;
  return(stat->avg);
  
}
#endif
