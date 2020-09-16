/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
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
/*! \file triangle_agi.upc
 * \brief The intuitive implementation of triangle counting 
 * that uses generic global references
 */

#include "sssp.h"

/*!
 * \brief This routine "relaxes" an edge in the Bellman-Ford algorithm
 * \param *tent shared tentative distance array
 * \param J the index of the head of the edge to be relaxed
 * \param new_tent the new tentative distance 
 *        (the weight of the tentative distance to the tail plus the weight of the edge).
 * This is effectively an atomic min operation of the new and old tentative distances.
 */
static void relax_bellman_agi(d_array_t *tent, int64_t J, double new_tent)
{
  // Allows us to use the bits of new and old, either as doubles for the
  // determining the min, or just as bits for the compare_and_swap.
  union{double x; uint64_t u;} old;
  union{double x; uint64_t u;} new;

  new.x = new_tent;
  while(1){
    old.x = lgp_get_double(tent->entry, J);
    if(new.x > old.x) 
      break;
    if( old.u == lgp_cmp_and_swap(tent->entry, J, old.u, new.u) )
      break;
  }
}

/*!
 * \brief This routine implements the Bellman-Ford algorithm with the agi model
 *
 * \param *dist a place to return the distance to the source
 * \param *mat the input matrix
 * \return average run time
 */
double sssp_bellman_agi(d_array_t *tent, sparsemat_t * mat, int64_t v0)
{
  T0_printf(" Running AGI SSSP Bellman-Ford\n");

  double t1 = wall_seconds();

  if(!mat){ T0_printf("ERROR: sssp_bellman_agi: NULL L!\n"); return(-1); }
  
  double tent_tail, new_tent;
  int64_t i, li, k, loop;
  int64_t lnumrows = mat->lnumrows;

  lgp_barrier();
  if( MYTHREAD == 0 ){
    lgp_put_double(tent->entry, v0, 0.0);
  }
  lgp_barrier();

  //dump_tent("AGI: ", tent);

  for(loop=0; loop<mat->numrows; loop++){
    for(li=0; li<lnumrows; li++){ 
      tent_tail = tent->lentry[li];
      for(k = mat->loffset[li]; k < mat->loffset[li+1]; k++){
        new_tent = tent_tail + mat->lvalue[k];
        relax_bellman_agi(tent, mat->lnonzero[k], new_tent);
      }
    }
    lgp_barrier();

    //dump_tent("AGI: ", tent);
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  return(stat->avg);
}

