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

double histo_conveyor(int64_t *pckindx, int64_t T,  int64_t *lcounts)
{
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col;
  int64_t pop_col;
  size_t n_local = 1; //set to cores/socket
  char *number = getenv("CONVEY_NLOCAL");
  if(number)
    n_local = atoi(number);
  
  minavgmaxD_t stat[1];

  int status = EXIT_FAILURE;
  convey_t* conveyor = convey_new(sizeof(long), SIZE_MAX, n_local, NULL, convey_opt_SCATTER);
  //convey_t* conveyor = convey_new_tensor(1024, sizeof(long), 2, n_local, 2, NULL, convey_opt_DYNAMIC);
  if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}
  
  ret = convey_begin(conveyor);
  if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

  lgp_barrier();  
  tm = wall_seconds();
  i = 0UL;
  while(convey_advance(conveyor, i == T)) {
    for(; i< T; i++){
      col = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if( !convey_push(conveyor, &col, pe))
	    break;
    }
    while( convey_pull(conveyor, &pop_col, NULL) == convey_OK)
    //while( (pop_col = convey_pull(conveyor, NULL)) != NULL)
      lcounts[pop_col] += 1;
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  convey_free(conveyor);
  return( stat->avg );
}
