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

/*!
 * \brief Relax the head of the edges delivered by a conveyor buffer
 * \param tent pointer to the tentative distances array
 * \param *conv the conveyor
 * \param done the signal to convey_advance that this thread is done
 * \return the return value from convey_advance
 */
static int64_t bellman_convey_relax_process(d_array_t *tent, convey_t *conv, int64_t done) 
{
  struct sssp_pkg_t pkg;
    
  while(convey_pull(conv, &pkg, NULL) == convey_OK){
    if( tent->lentry[pkg.lj] > pkg.tw ) {
      if(0){printf("Convey: replace %ld weight %lg with %lg\n", pkg.lj*THREADS + MYTHREAD, tent->lentry[pkg.lj], pkg.tw);}
      tent->lentry[pkg.lj] = pkg.tw;
    }
  }
  return( convey_advance(conv, done) );
}

/*!
 * \brief Push the potentially improved weight to the thread handling the head of the edge
 * \param conv the extack buffers
 * \param tent pointer to the tentative distances array to be passed thru to bellman_convey_relax_process
 * \param J the head of the edge, given by it global name
 * \param tw the new weight
 * \return the value from the push
 */
static int64_t bellman_convey_push(convey_t *conv, d_array_t *tent, int64_t J, double tw)
{
  int64_t ret, pe;
  sssp_pkg_t pkg;
  pe     = J % THREADS;
  pkg.lj = J / THREADS;
  pkg.tw = tw;
  if((ret = convey_push(conv, &pkg, pe)) == 0){
    bellman_convey_relax_process(tent, conv, 0);
  }
  return(ret);
}


/*!
* \brief This routine implements the conveyor variant of Bellman-Ford algorithm
 * \param *ltent locat pointer to the array of the tentative weights to each vertex
 * \param *mat the sparse matrix that holds the graph
 * \param v0 is the starting vertex
 * \return average run time
 */
double sssp_bellman_convey(d_array_t *dist, sparsemat_t *mat, int64_t v0) 
{
  int64_t k, li, J, pe;
  int64_t pe_v0, li_v0;
  int64_t loop;
  int64_t changed;
  sssp_pkg_t  pkg;
  d_array_t *tent0, *tent1, *tent2;
  d_array_t *tent_old, *tent_cur, *tent_new, *tent_temp;

  if(0){printf(" Running Conveyor SSSP Bellman-Ford\n");}

  convey_t * conv = convey_new(SIZE_MAX, 0, NULL, 0);
  if(conv == NULL){return(-1);}

  double t1 = wall_seconds();

  pe_v0 = v0 % THREADS;
  li_v0 = v0 / THREADS;

  tent0 = init_d_array(dist->num);
  set_d_array(tent0, INFINITY);
  lgp_barrier();
  if(pe_v0 == MYTHREAD){
    tent0->lentry[li_v0] = 0.0;
  }
  //dump_tent("\nConvey: 0", tent0);

  tent1 = copy_d_array(tent0);
  if(pe_v0 == MYTHREAD ){
    for(k = mat->loffset[li_v0]; k < mat->loffset[li_v0+1]; k++){
      lgp_put_double(tent1->entry, mat->lnonzero[k], mat->lvalue[k]);
    }
  }

  lgp_barrier();

  //dump_tent("Convey: 1", tent1);

  lgp_barrier();
  tent2 = init_d_array(tent1->num);

  tent_old = tent0;
  tent_cur = tent1;
  tent_new = tent2;

  for(loop=0; loop<mat->numrows; loop++){
    changed = 0;
    convey_begin(conv, sizeof(sssp_pkg_t));
    for(li=0; li < mat->lnumrows; li++)
      tent_new->lentry[li] = tent_cur->lentry[li];

    for(li=0; li < mat->lnumrows; li++){
      if(tent_old->lentry[li] == tent_cur->lentry[li])
        continue;
      changed = 1;
      for(k=mat->loffset[li]; k< mat->loffset[li+1]; k++){
        if(0){printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, li * THREADS + MYTHREAD, J, tent_cur->lentry[li],  pkg.tw);}
        if( bellman_convey_push(conv, tent_new, mat->lnonzero[k],  tent_cur->lentry[li] + mat->lvalue[k]) == 0)
          k--;
      }
    }
    while(bellman_convey_relax_process(tent_new, conv, 1))// keep popping til all threads are done
      ;
    lgp_barrier();
    if( lgp_reduce_add_l(changed) == 0 ){                 // keep looping until nothing changes
      replace_d_array(dist, tent_new);
      break;
    }
    //dump_tent("Convey:  ", tent_new);
    tent_temp = tent_old;
    tent_old = tent_cur;
    tent_cur = tent_new;
    tent_new = tent_temp;

    lgp_barrier();
    convey_reset(conv);
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  return(stat->avg);
}
