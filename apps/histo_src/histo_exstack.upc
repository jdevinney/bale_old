/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
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
 * \param bufsiz the exstack bufsiz
 * \return average run time
 *
 */
double histo_exstack(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz) 
{
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col, *colp;
  minavgmaxD_t stat[1];
  exstack_t * ex = exstack_init(bufsiz, sizeof(int64_t));
  if( ex == NULL) return(-1.0);

  lgp_barrier();  
  tm = wall_seconds();
  i = 0UL;

  while( exstack_proceed(ex, (i==T)) ){
    for( ; i < T; i++){
      col = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if( !exstack_push(ex, &col, pe) )
        break;
    }
    
    exstack_exchange(ex);
    
    while(colp = exstack_pull(ex, NULL))
      lcounts[*colp]++;
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  exstack_clear(ex);
  free(ex);
  return( stat->avg );
}

