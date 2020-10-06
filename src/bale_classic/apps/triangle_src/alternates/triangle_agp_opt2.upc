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
/*! \file triangle_agp_opt2.upc
 * \brief Another intuitive implementation of triangle counting 
 * that uses generic global references. This implementation differs from
 * triangle_agp_opt1 only in that we do some buffering when reading remote rows.
 */

#include "triangle.h"

// this depends on the row nonzeros being sorted
// optimizations:
// 1. fetches row nonzeros in buffers (rather than one word at a time)
// 2. search for intersections in two rows is done more intelligently (since nonzeros are sorted)
// I think one thing that helps this model over conveyors/exstack is that it is more cache friendly
// since it is crawling over the l_i row many times all in one go, before moving on to row l_i + 1.
/*!
 * \brief This routine implements another AGP variant of triangle counting
 * \param *count a place to return the counts from this thread
 * \param *sr a place to return the number of shared references
 * \param *L the lower triangle of the input sparse matrix 
 * \param *U the upper triangle of the input sparse matrix 
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 */
double triangle_agp_opt2(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, ii, k, kk, kt, num2pull, L_i, L_j;

/*! \brief pull remote gets in a buffer of this size */
#define PULL_BUF_SIZE 64
  int64_t w[PULL_BUF_SIZE];

  double t1 = wall_seconds();
  //foreach nonzero (i, j) in L

  for(l_i = 0; l_i < L->lnumrows; l_i++){ 
    for(k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
      L_i = l_i*THREADS + MYTHREAD;
      L_j = L->lnonzero[k];
      assert( L_j < L_i );

      numpulled+=2;  // to get the offsets in their shared form
      int64_t row_start, row_end, inc, num2pull, num2pull_now, start;
      if(alg == 0){
        row_start = lgp_get_int64(L->offset, L_j);
        row_end   = lgp_get_int64(L->offset, L_j + THREADS);
        inc = PULL_BUF_SIZE;
        num2pull = row_end - row_start;
        start = L->loffset[l_i];
      }else{
        row_start = lgp_get_int64(U->offset, L_j + THREADS) - 1;
        row_end   = lgp_get_int64(U->offset, L_j);
        inc = -PULL_BUF_SIZE;
        num2pull = row_start - row_end;
        start = U->loffset[l_i + 1] - 1;
      }
      numpulled += num2pull;
      
      for( kk = row_start; kk >= row_end; kk += inc ){
        //num2pull = ((row_end - kk) <= PULL_BUF_SIZE) ? row_end - kk : PULL_BUF_SIZE;
        num2pull_now = (PULL_BUF_SIZE > num2pull ? num2pull : PULL_BUF_SIZE);
        num2pull -= num2pull_now;
        if(alg == 0){
          lgp_memget(w,
                     L->nonzero,
                     num2pull_now*sizeof(int64_t),
                     kk*THREADS + L_j % THREADS);
          
          for( kt = 0; kt < num2pull_now; kt++){ 
            
            for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
              if( w[kt] ==  L->lnonzero[ii] ){ 
                cnt++;
                start = ii + 1;
                break;
              }
              if( w[kt] < L->lnonzero[ii] ){ // the rest are all bigger too, cause L is tidy 
                start = ii; // since pulled row is sorted, we can change start
                break;
              }
            }
          }
        
        }else{
          lgp_memget(w,
                     U->nonzero,
                     num2pull_now*sizeof(int64_t),
                     L_j%THREADS + (kk - num2pull_now)*THREADS);
          for( kt = num2pull_now - 1; kt >= 0; kt--){
            
            //if( w[kt] < L_i){kk = row_end; break;}
            
            for(ii = start; ii >= U->loffset[l_i]; ii--) {
              if( w[kt] == U->lnonzero[ii] ){ 
                cnt++;
                start = ii - 1;
                break;
              }
              if( w[kt] > U->lnonzero[ii] ){ // the rest are all bigger too, cause U is tidy 
                start = ii; // since pulled row is sorted, we can change start
                break;
              }
            }
          }
        }
      }
    }
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  *sr = numpulled; 
  *count = cnt;
  return(stat->avg);
}
