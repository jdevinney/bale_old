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
/*! \file histo_exstack.upc
 * \brief The exstack implementation of histogram.
 */
#include "histo.h"

/*!
 * \brief This routine implements the exstack classic variant of histogram.
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param T the length of the pcindx array
 * \param *lcounts localized pointer to the count array.
 * \param buf_cnt the size of the exstack buffers in packages
 * \return average run time
 *
 */
//double histo_exstack(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t buf_cnt) {
double histo_exstack(histo_t * data, int64_t buf_cnt){
  int64_t i;
  double tm;
  int64_t pe, col, *colp;
  minavgmaxD_t stat[1];
  exstack_t * ex = exstack_init(buf_cnt, sizeof(int64_t));
  if( ex == NULL) return(-1.0);
  
  lgp_barrier();  
  tm = wall_seconds();
  i = 0UL;

  while( exstack_proceed(ex, (i==data->l_num_ups)) ){
    int64_t popped = 0;
    for( ; i < data->l_num_ups; i++){
      col = data->pckindx[i] >> 20;
      pe  = data->pckindx[i] & 0xfffff;
      assert(pe < THREADS);
      if( !exstack_push(ex, &col, pe) )
        break;
    }

    exstack_exchange(ex);
    
    while((colp = exstack_pull(ex, NULL))){
      popped++;
      assert(*colp < data->lnum_counts);
      data->lcounts[*colp]++;
    }
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  exstack_clear(ex);
  free(ex);
  return( stat->avg );
}

