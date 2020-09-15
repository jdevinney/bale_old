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

/*! \file ig_exstack_function.upc
 * \brief The classic exstack implementation of indexgather using functions to hide the loops.
 */
#include "ig.h"
  
typedef struct index_val_t {
  int64_t i;    
  int64_t idx;
} index_val_t;

/*!
 * \brief  routine to service the pops
 * \param *ex the exstack struct
 * \param *tgt the target array
 * \param *ltable localized pointer to the count array
 * \param imdone flag for this thread
 * \return average run time
 *
 */
void ig_exstack_service(exstack_t *ex, int64_t *tgt, int64_t *ltable, int64_t imdone ) {
  int64_t fromth;
  index_val_t pkg;

  do {
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg , &fromth)) {
      pkg.idx  = ltable[pkg.idx];
      exstack_push(ex, &pkg, fromth); 
    }   
    lgp_barrier();
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg, NULL)){ 
      tgt[pkg.i] = pkg.idx;
    }
  } while (imdone && exstack_proceed( ex , imdone ));
}


/*!
 * \brief  routine to generate the requests
 * \param *ex the exstack struct
 * \param *pckindx array of "packed" indices
 * \param i the ith request
 * \param *tgt the target array
 * \param *ltable localized pointer to the count array
 * \return average run time
 *
 */
void ig_exstack_request(exstack_t *ex, int64_t *pckindx, int64_t i,  int64_t *tgt, int64_t *ltable) {
  int64_t l_indx, pe;
  index_val_t pkg;

  pkg.i = i;
  pkg.idx = pckindx[i] >> 16;
  pe = pckindx[i] & 0xffff;
  if( exstack_push(ex, &pkg, pe) > 1 ) 
    return;

  ig_exstack_service( ex, tgt, ltable, 0 );
}


/*!
 * \brief This routine implements the exstack classic variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array
 * \param buf_cnt the exstack buffer size in packages
 * \return average run time
 *
 */
double ig_exstack_function(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt) {
  double tm;
  int imdone;
  int64_t ret;
  int64_t room;
  int64_t l_indx, idx, i, i0, j;
  int64_t pe, fromth;
  minavgmaxD_t stat[1];

  index_val_t pkg;

  exstack_t * ex = exstack_init(buf_cnt, sizeof(index_val_t));
  if( ex == NULL ) return(-1.0);
  
  lgp_barrier();
  tm = wall_seconds();

  for(i=0; i<l_num_req; i++) {
    ig_exstack_request( ex, pckindx, i, tgt, ltable ); 
  }
  ig_exstack_service(ex, tgt, ltable, 1);

  lgp_barrier();
  
  tm = wall_seconds() - tm;
  lgp_barrier();
  lgp_min_avg_max_d( stat, tm, THREADS );

  exstack_clear(ex);
  free(ex);
  return( stat->avg );
}
