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

/*! \file ig_agi.upc
 * \brief The intuitive implementation of indexgather that uses single word gets to shared addresses.
 */
#include "ig.h"

/*!
 * \brief This routine implements the single word get version indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the global array of counts
 * \param T the length of the index array
 * \param *table shared pointer to the shared table array.
 * \return average run time
 *
 */
double ig_gets(int64_t *tgt, int64_t *index, int64_t T,  SHARED int64_t *table) 
{
  int64_t i;
  double tm;
  minavgmaxD_t stat[1];

  lgp_barrier();
  tm = wall_seconds();

  for(i = 0; i < T; i++){
#pragma pgas defer_sync
    tgt[i] = lgp_get_int64(table, index[i]);
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  return( stat->avg );
}

