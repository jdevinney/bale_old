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
/*! \file histo_agi.upc
 * \brief The intuitive implementation of histogram that uses global atomics.
 */
#include "histo.h"

/*!
 * \brief This routine implements straight forward, 
 *         single word atomic updates to implement histogram.
 * \param *index array of indices into the shared array of counts.
 * \param l_num_ups the local length of the index array
 * \param *counts SHARED pointer to the count array.
 * \return average run time
 *
 */

double histo_agi(int64_t *index, int64_t l_num_ups,  SHARED int64_t *counts) {
  double tm;
  int64_t i;
  minavgmaxD_t stat[1];

  lgp_barrier();
  tm = wall_seconds();

  for(i = 0; i < l_num_ups; i++) {
#pragma pgas defer_sync
    //_amo_aadd(&counts[index[i]], 1);
    //counts[index[i]] += 1;
    lgp_atomic_add(counts, index[i], 1);
  }
  lgp_barrier();
  tm = wall_seconds() - tm;
  lgp_min_avg_max_d( stat, tm, THREADS );
  
  return( stat->avg );
}

