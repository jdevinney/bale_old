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

/*! \file ig_exstack2_goto.upc
 * \brief An exstack2 implementation of indexgather that uses goto's instead of loops.
 */
#include "ig.h"

/*!
 * \brief This routine implements the exstack2 variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array.
 * \param buf_cnt the exstack buffer size in packages
 * \return average run time
 *
 */
double ig_exstack2_goto(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt) {
  double tm;
  int64_t more_req = 1;
  int64_t ret;
  int64_t i;
  int64_t pe, fromth, fromth2;
  minavgmaxD_t stat[1];


  typedef struct index_val_t{
    int64_t i;    
    int64_t idx;
  }index_val_t;

  index_val_t pkg, pkge, pkgr;

  exstack2_t * ex2r = exstack2_init(buf_cnt , sizeof(index_val_t));
  exstack2_t * ex2e = exstack2_init(buf_cnt , sizeof(index_val_t));
  if( (ex2r == NULL) || (ex2e == NULL) ) return(-1.0);

  uint64_t more2req;
  
  void * jmp_req = &&ig_new_req;
  void * jmp_echo = &&ig_new_echo;

  lgp_barrier();
  tm = wall_seconds();

  i=0;
  more2req = 1;

ig_new_req:
  pkgr.i = i;
  pkgr.idx = pckindx[i] >> 16;
  pe = pckindx[i] & 0xffff;

ig_req_retry:
  if ( exstack2_push(ex2r, &pkgr, pe) ) {
    i++;
    if( i == l_num_req )
       goto ig_req_flush;
    goto ig_new_req;
  } else {
    jmp_req = &&ig_req_retry;
    goto *jmp_echo;
  }

ig_req_flush:
  jmp_req = &&ig_req_flush;
  more2req = exstack2_proceed(ex2r, 1);
  goto *jmp_echo;

ig_new_echo:
  if( exstack2_pop(ex2r, &pkg, &fromth2) ){
    pkge.i = pkg.i;
    pkge.idx = ltable[pkg.idx];
ig_echo_retry:
    if( exstack2_push(ex2e, &pkge, fromth2) ) {
       jmp_echo = &&ig_new_echo;
       goto ig_new_echo;
    } else {
      jmp_echo = &&ig_echo_retry;
      goto ig_apply_echo;
    }
  }

ig_apply_echo:
  while(exstack2_pop(ex2e, &pkg, NULL))
    tgt[pkg.i] = pkg.idx;

  if( more2req )
    goto *jmp_req;

  if( exstack2_proceed(ex2e, 1) )
    goto *jmp_echo;

  tm = wall_seconds() - tm;
  lgp_barrier();

  lgp_min_avg_max_d( stat, tm, THREADS );
  exstack2_clear(ex2r);
  exstack2_clear(ex2e);
  free(ex2r);
  free(ex2e);
  return( stat->avg );
}

