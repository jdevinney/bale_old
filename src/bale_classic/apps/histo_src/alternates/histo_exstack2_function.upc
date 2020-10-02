
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
/******************************************************************
 *****************************************************************/ 

/*! \file histo_exstack2_function.upc
 * \brief The exstack2 implementation of histogram with the loops hidden in functions.
 */
#include "histo.h"

/*!
 * \brief histo_exstack2_service routine handles the popping and updating the lcounts array
 * \param *ex2 pointer to the exstack2 struct
 * \param *lcounts local view of the counts array
 * \param done_pushing flag to signal the end game
 * \return void
 */
void histo_exstack2_service( exstack2_t *ex2, int64_t *lcounts, int64_t done_pushing) {
  int64_t *l_idx;
  do {
    while((l_idx = exstack2_pull(ex2, NULL)) != 0)
      lcounts[*l_idx]++;
  } while( done_pushing && exstack2_proceed(ex2, done_pushing) );
}

/*!
 * \brief histo_exstack2_request routine pushes the indices 
 * \param *ex2 pointer to the exstack2 struct
 * \param pckidx the packed index holding the requested pe and local index 
 * \param *lcounts local view of the counts array,  needs to be passed to the service function
 * \return void
 */
void histo_exstack2_request(exstack2_t *ex2, int64_t pckidx,  int64_t *lcounts) {
  int64_t pe_idx, l_idx, ret;
  l_idx  =  pckidx >> 16;
  pe_idx =  pckidx & 0xffff;

  if( exstack2_push(ex2, &l_idx, pe_idx) == 1 ) { // last one to work
     do {
       ret = exstack2_send(ex2, pe_idx, 0);
       histo_exstack2_service(ex2, lcounts, 0);   // steady state pushing and popping
     } while( ret == 0 );
  }
}


/*!
 * \brief histo_exstack2_function calls the request and service functions
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_ups number of updates on this thread
 * \param *lcounts local view of the counts array
 * \param buf_cnt the size of the exstack buffers in packages
 * \return average run time
 *
 */

double histo_exstack2_function(int64_t *pckindx, int64_t l_num_ups,  int64_t *lcounts, int64_t buf_cnt) {
  int ret;
  double tm;
  int64_t i, pe, col, idx, *idxp;
  minavgmaxD_t stat[1];
  exstack2_t * ex2 = exstack2_init(buf_cnt, sizeof(int64_t));
  if( ex2 == NULL ) return(-1.0);

  lgp_barrier();
  tm = wall_seconds();

  for( i=0; i<l_num_ups; i++ ){
    histo_exstack2_request(ex2, pckindx[i], lcounts);
  }
  histo_exstack2_service(ex2, lcounts, 1);

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  exstack2_clear(ex2);
  free(ex2);
  return( stat->avg );
}

