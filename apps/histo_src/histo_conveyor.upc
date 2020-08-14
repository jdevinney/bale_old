/******************************************************************
//
//
//  Copyright(C) 2019, Institute for Defense Analyses
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

/*! \file histo_conveyor.upc
 * \brief A conveyor implementation of histogram.
 */
#include "histo.h"

/*!
 * \brief This routine implements histogram using conveyors
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param T the length of the pcindx array
 * \param *lcounts localized pointer to the count array.
 * \return average run time
 *
 */

//double histo_conveyor(int64_t *pckindx, int64_t T,  int64_t *lcounts) {
double histo_conveyor(histo_t * data){
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col;
  int64_t pop_col;
  
  minavgmaxD_t stat[1];

  int status = EXIT_FAILURE;
  convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}
  
  ret = convey_begin(conveyor, sizeof(int64_t));
  if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

  lgp_barrier();  
  tm = wall_seconds();
  i = 0UL;
  while(convey_advance(conveyor, i == data->l_num_ups)) {
    for(; i< data->l_num_ups; i++){
      col = data->pckindx[i] >> 20;
      pe  = data->pckindx[i] & 0xfffff;
      assert(pe < THREADS);
      if( !convey_push(conveyor, &col, pe))
	    break;
    }
    while( convey_pull(conveyor, &pop_col, NULL) == convey_OK){
    //while( (pop_col = convey_pull(conveyor, NULL)) != NULL)
      assert(pop_col < data->lnum_counts);
      data->lcounts[pop_col] += 1;
    }
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  convey_free(conveyor);
  return( stat->avg );
}
