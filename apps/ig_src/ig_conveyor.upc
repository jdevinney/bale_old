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

/*! \file ig_conveyor.upc
 * \brief A conveyor implementation of indexgather.
 */
#include "ig.h"

/*!
 * \brief This routine implements the conveyor variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array.
 * \return average run time
 *
 */
double ig_conveyor(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable) {
  double tm;
  int64_t pe, fromth, fromth2;
  int64_t i = 0, from;
  minavgmaxD_t stat[1];
  bool more;

  typedef struct pkg_t {
    int64_t idx;
    int64_t val;
  } pkg_t;
  pkg_t pkg;
  pkg_t *ptr = calloc(1, sizeof(pkg_t));
  
  convey_t* requests = convey_new(sizeof(pkg_t), SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  assert( requests != NULL );
  convey_t* replies = convey_new(sizeof(pkg_t), SIZE_MAX, 0, NULL, 0);
  assert( replies != NULL );

  convey_begin(requests);
  convey_begin(replies);
  lgp_barrier();
  
  tm = wall_seconds();

  i = 0;
  while (more = convey_advance(requests, (i == l_num_req)),
         more | convey_advance(replies, !more)) {

    for (; i < l_num_req; i++) {
      pkg.idx = i;
      pkg.val = pckindx[i] >> 16;
      pe = pckindx[i] & 0xffff;
      if (! convey_push(requests, &pkg, pe))
        break;
    }

    while (convey_pull(requests, ptr, &from) == convey_OK) {
      pkg.idx = ptr->idx;
      pkg.val = ltable[ptr->val];
      if (! convey_push(replies, &pkg, from)) {
        convey_unpull(requests);
        break;
      }
    }

    while (convey_pull(replies, ptr, NULL) == convey_OK)
      tgt[ptr->idx] = ptr->val;
  }

  tm = wall_seconds() - tm;
  free(ptr);
  lgp_barrier();

  lgp_min_avg_max_d( stat, tm, THREADS );
  convey_free(requests);
  convey_free(replies);
  return( stat->avg );
}
