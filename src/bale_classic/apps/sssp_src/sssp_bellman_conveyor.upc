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

/*! \file sssp_conveyor.upc
 * \brief Application that implements Bellman-Ford with conveyors
 */

#include "sssp.h"

typedef struct conv_bellman_t {
  int64_t i;  // row[i] is "tail" of the edge, in case we want to set backpointers
  int64_t lj; // the local "head" of the edge
  double tw;  // new tentative weight
}conv_bellman_t ;

/*!
 * \brief pop routine to implement relaxing the edges
 * \param tent pointer to the tentative distances array
 * \param *conv the conveyor
 * \param done the signal to convey_advance that this thread is done
 * \return the return value from convey_advance
 */
static int64_t bellman_convey_relax_process(d_array_t *tent, convey_t *conv, int64_t done) 
{
  struct conv_bellman_t pkg;
    
  while(convey_pull(conv, &pkg, NULL) == convey_OK){
    if( tent->lentry[pkg.lj] > pkg.tw ) {
			printf("Convey: replace %ld weight %lg with %lg\n", pkg.lj*THREADS + MYTHREAD, tent->lentry[pkg.lj], pkg.tw);
      tent->lentry[pkg.lj] = pkg.tw;
		}
  }
  return( convey_advance(conv, done) );
}

/*!
* \brief This routine implements the conveyor variant of Bellman-Ford algorithm
 * \param *ltent locat pointer to the array of the tentative weights to each vertex
 * \param *mat the sparse matrix that holds the graph
 * \param v0 is the starting vertex
 * \return average run time
 */
double sssp_bellman_convey(d_array_t *tent, sparsemat_t *mat, int64_t v0) 
{
	printf(" Running Conveyor SSSP Bellman-Ford\n");

  convey_t * conv = convey_new(SIZE_MAX, 0, NULL, 0);
  if(conv == NULL){return(-1);}

  double t1 = wall_seconds();

  conv_bellman_t  pkg;
  int64_t k, pe, li, J;
  int64_t loop;

  set_d_array(tent, INFINITY);
	lgp_barrier();
	if( MYTHREAD == 0 ){
  	lgp_put_double(tent->entry, v0, 0.0);
  }
	lgp_barrier();

  dump_tent("Convey: ", tent);
	lgp_barrier();

  for(loop=0; loop<mat->numrows; loop++){
    convey_begin(conv, sizeof(conv_bellman_t ));
    for(li=0; li < mat->lnumrows; li++){
//      if( ltent[li]  == INFINITY )
//        continue;
      pkg.i = li * THREADS + MYTHREAD;

      for(k=mat->loffset[li]; k< mat->loffset[li + 1]; k++){
        J = mat->nonzero[k];
        pe  = J % THREADS;
        pkg.lj = J / THREADS;
        pkg.tw = tent->lentry[li] + mat->value[k];
        //printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, pkg.i, J, tent->lentry[li],  pkg.tw); 
        if( convey_push(conv, &pkg, pe) != convey_OK ){
          bellman_convey_relax_process(tent, conv, 0); 
          k--;
        }
      }
    }
    while(bellman_convey_relax_process(tent, conv, 1))// keep popping til all threads are done
      ;
    lgp_barrier();
	  dump_tent("Convey : ", tent);
		convey_reset(conv);
  }

  lgp_barrier();

  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  return(stat->avg);
}
