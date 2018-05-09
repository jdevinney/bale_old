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
/*! \file histo_exstack2.upc
 * \brief The exstack2 implementation of histogram.
 */
#include "histo.h"

/*!
 * \brief This routine implements the exstack2 variant of histogram.
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param T the length of the pcindx array
 * \param *lcounts localized pointer to the count array.
 * \param bufsiz the size of the buffers in exstack2
 * \return average run time
 *
 */
double histo_exstack2(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz) 
{
  int ret;
  double tm;
  int64_t pe, col, idx, *idxp;
  minavgmaxD_t stat[1];
  exstack2_t * ex2 = exstack2_init(bufsiz, sizeof(int64_t));
  if( ex2 == NULL ) return(-1.0);

  lgp_barrier();
  tm = wall_seconds();
  int64_t i = 0;

  while(exstack2_proceed( ex2, i==T )) {
    for( ; i < T; i++){
      col = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if( !exstack2_push(ex2, &col, pe) )
        break;
    }

    while(idxp = exstack2_pull(ex2, NULL))
      lcounts[*idxp]++;
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  exstack2_clear(ex2);
  free(ex2);
  return( stat->avg );
}
