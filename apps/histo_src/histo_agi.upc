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
/*! \file histo_agi.upc
 * \brief The intuitive implementation of histogram that uses global atomics.
 */
#include "histo.h"

/*!
 * \brief This routine implements straight forward, 
 *         single word atomic updates to implement histogram.
 * \param *index array of indices into the shared array of counts.
 * \param T the length of the index array
 * \param *counts SHARED pointer to the count array.
 * \param v the amount to substract from counts array
 *         assume an atomicadd(-v) is equivalent to atomic increment
 *         we can use this rountine to zero (and hence check) the
 *         the results from running other versions of histo.
 * \return average run time
 *
 */
double histo_atomic(int64_t *index, int64_t T,  SHARED int64_t *counts, int64_t v) 
{
  double tm;
  int64_t i;
  minavgmaxD_t stat[1];

  lgp_barrier();
  tm = wall_seconds();

#pragma pgas defer_sync
  for(i = 0; i < T; i++) {
    lgp_atomic_add(counts, index[i], -v);
  }
  lgp_barrier();
  tm = wall_seconds() - tm;
  lgp_min_avg_max_d( stat, tm, THREADS );
  
  return( stat->avg );
}

