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

/*! \file ig_conveyor.upc
 * \brief A conveyor implementation of indexgather.
 */
#include "ig.h"

/*!
 * \brief This routine implements the conveyor variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param T the length of the pcindx array
 * \param *ltable localized pointer to the count array.
 * \return average run time
 *
 */
double ig_conveyor(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable)
{
  double tm;
  int64_t pe, fromth, fromth2;
  int64_t i = 0, from;
  minavgmaxD_t stat[1];
  bool more;

  typedef struct pkg_t {
    int64_t idx;
    int64_t val;
  } pkg_t;
  pkg_t pkg;
  pkg_t *ptr = calloc(1, sizeof(pkg_t));

  size_t n_local = 1;
  char *number = getenv("CONVEY_NLOCAL");
  if(number)
    n_local = atoi(number);
  
  convey_t* requests = convey_new(sizeof(pkg_t), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  assert( requests != NULL );
  convey_t* replies = convey_new(sizeof(pkg_t), SIZE_MAX, n_local, NULL, 0);
  assert( replies != NULL );

  convey_begin(requests);
  convey_begin(replies);
  lgp_barrier();
  
  tm = wall_seconds();

  i = 0;
  while (more = convey_advance(requests, (i == T)),
         more | convey_advance(replies, !more)) {

    for (; i < T; i++) {
      pkg.idx = i;
      pkg.val = pckindx[i] >> 16;
      pe = pckindx[i] & 0xffff;
      if (! convey_push(requests, &pkg, pe))
        break;
    }

    while (convey_pull(requests, ptr, &from) == convey_OK) {
      pkg.idx = ptr->idx;
      pkg.val = ltable[ptr->val];
      if (! convey_push(replies, &pkg, from)) {
        convey_unpull(requests);
        break;
      }
    }

    while (convey_pull(replies, ptr, NULL) == convey_OK)
      tgt[ptr->idx] = ptr->val;
  }

  tm = wall_seconds() - tm;
  free(ptr);
  lgp_barrier();

  lgp_min_avg_max_d( stat, tm, THREADS );
  convey_free(requests);
  convey_free(replies);
  return( stat->avg );
}
