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
/*! \file histo_exstack.upc
 * \brief The exstack implementation of histogram.
 */
#include "histo.h"

/*!
 * \brief This routine implements the exstack classic variant of histogram.
 * \param data the histo_t struct that carries all the parameters for the implementations
 * \param buf_cnt the number of packages in the exstack buffers
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

