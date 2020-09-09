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

/*! \file sssp_bellman_exstack.upc
 * \brief Implementation of sssp using exstack
 */

#include "sssp.h"

typedef struct pkg_bellman_e_t{
  int64_t i;  // row[i] is "tail" of the edge, in case we want to set backpointers
  int64_t lj; // the local "head" of the edge
  double tw;  // new tentative weight
}pkg_bellman_e_t;

/*!
 * \brief 
 */

/*!
 * \brief pop routine to implement relaxing the edges
 * \param tent pointer to the tentative distances array
 * \param *ex the extack buffers
 * \param done the signal to exstack_proceed that this thread is done
 * \return the return value from exstack_proceed
 */
static int64_t bellman_exstack_relax_process(d_array_t *tent, exstack_t *ex, int64_t done) 
{
  int64_t fromth;
  pkg_bellman_e_t pkg;

  exstack_exchange(ex);
  
  while(exstack_pop(ex, &pkg, &fromth)){
    if( tent->lentry[pkg.lj] > pkg.tw ){
      if(0){printf("Ex: replace %ld weight %lg with %lg\n", pkg.lj*THREADS + MYTHREAD, tent->lentry[pkg.lj], pkg.tw);}
      tent->lentry[pkg.lj] = pkg.tw;
    }
  }
  return( exstack_proceed(ex, done) );
}

/*!
 * \brief This routine implements the Bellman-Ford algorithm using exstack
 * \param *tent the SHARED array that holds the tentative distances
 * \param *mat the input matrix
 * \param v0 is the the staring row (vertex)
 * \return average run time
 */
double sssp_bellman_exstack(d_array_t *tent, sparsemat_t * mat, int64_t v0)
{
  double t1 = wall_seconds();

  if(!mat){ T0_printf("ERROR: sssp_bellman_exstack: NULL L!\n"); return(-1); }
  
  int64_t k, li, J, pe;
  int64_t pe_v0, li_v0;
  int64_t loop;
  int64_t changed;
  pkg_bellman_e_t pkg;
  d_array_t *tent0, *tent1, *tent2;
  d_array_t *tent_old, *tent_cur, *tent_new, *tent_temp;

  //TODO: Fix the buffer size 
  exstack_t * ex = exstack_init(1024, sizeof(pkg_bellman_e_t));
  if( ex == NULL) return(-1.0);


  pe_v0 = v0 % THREADS;
  li_v0 = v0 / THREADS;

  tent0 = init_d_array(tent->num);
  set_d_array(tent0, INFINITY);
  lgp_barrier();
  if(pe_v0 == MYTHREAD){
    tent0->lentry[li_v0] = 0.0;
  }
  dump_tent("\nExstack: 0", tent0);

  tent1 = copy_d_array(tent0);
  if(pe_v0 == MYTHREAD ){
    for(k = mat->loffset[li_v0]; k < mat->loffset[li_v0+1]; k++){
      lgp_put_double(tent1->entry, mat->lnonzero[k], mat->lvalue[k]);
    }
  }

  lgp_barrier();

  dump_tent("Exstack: 1", tent1);

  lgp_barrier();
  tent2 = init_d_array(tent->num);

  tent_old = tent0;
  tent_cur = tent1;
  tent_new = tent2;
  for(loop=2; loop<mat->numrows; loop++){
    changed = 0;
    for(li=0; li < mat->lnumrows; li++)
      tent_new->lentry[li] = tent_cur->lentry[li];

    for(li=0; li < mat->lnumrows; li++){
      if(tent_old->lentry[li] == tent_cur->lentry[li])
        continue;
      changed = 1;
      pkg.i = li * THREADS + MYTHREAD;
      for(k=mat->loffset[li]; k< mat->loffset[li+1]; k++){
        J = mat->lnonzero[k];
        pe  = J % THREADS;
        pkg.lj = J / THREADS;
        pkg.tw = tent->lentry[li] + mat->lvalue[k];
        if(0){printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, pkg.i, J, tent->lentry[li],  pkg.tw);}
        if( exstack_push(ex, &pkg, pe) == 0 ) {
            bellman_exstack_relax_process(tent_new, ex, 0);
            k--;
        }
      }
    }
    while(bellman_exstack_relax_process(tent_new, ex, 1))
      ;

    lgp_barrier();
    if( lgp_reduce_add_l(changed) == 0 ){
      replace_d_array(tent, tent_cur);
      break;
    }

    dump_tent("EXSTACK:  ", tent_new);

    tent_temp = tent_old;
    tent_old = tent_cur;
    tent_cur = tent_new;
    tent_new = tent_temp;

    lgp_barrier();
    exstack_reset(ex);
  }

  dump_tent("Bellman EXSTACK: ", tent);
  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  return(stat->avg);
}

