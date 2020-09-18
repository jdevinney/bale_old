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
/*! \file triangle_agi.upc
 * \brief The intuitive implementation of triangle counting 
 * that uses generic global references
 */

#include "sssp.h"


// C version at this point.
// If we had an atomic min this would be easy and slow


static void relax_edge( double *tw, double *tv, double weight )
 {
   if( *tw > *tv + weight ){
     *tw = *tv + weight;
   }
 }


/*!
 * \brief This routine implements the Bellman-Ford algorithm
 *
 * \param *dist a place to return the distance to the source
 * \param *dmat the input matrix
 * \return average run time
 */

double sssp_bellman_agi(double *tent, sparsemat_t * dmat, int64_t v0)
{
   
  double t1 = wall_seconds();

  if(!dmat){ T0_printf("ERROR: sssp_bellman_agi: NULL L!\n"); return(-1); }
  
  int64_t i, k, loop;
  int64_t lnumrows = dmat->lnumrows;

  for(i=0; i<lnumrows; i++)
    tent[i] = INFINITY;

  tent[v0] = 0.0;

  for(loop=0; loop<numrows; loop++){
    for(i=0; i<numrows; i++){ 
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        relax_edge( &(tent[ mat->nonzero[k] ]), &(tent[i]), mat->value[k]);
      }
    }
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  return(stat->avg);
}

