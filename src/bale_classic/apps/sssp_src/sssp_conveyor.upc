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

/*! \file toposort_conveyor.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "sssp.h"

static convey_relax_edge(double *tent, sparsemat_t *dmat, int64_t I, int64_t J, double weight_IJ );

typedef struct pkg_bellman_c_t{
  int64_t i;  // row[i] is "tail" of the edge, in case we want to set backpointers
  int64_t lj; // the local "head" of the edge
  double tw;  // new tentative weight
}pkg_bellman_c_t;

/*!
* \brief This routine implements the conveyor variant of Bellman-Ford algorithm
 * \param *tent the tentative weights to each vertex
 * \param *dmat the sparse matrix that holds the graph
 * \param v0 is the starting vertex
 * \return average run time
 */
double sssp_bellman_convey(double *tent, sparsemat_t *dmat, int64_t v0) 
{
  convey_t * conv = convey_new(SIZE_MAX, 0, NULL, 0);
  if(conv == NULL){return(-1);}
  if(convey_begin( conv, sizeof(pkg_bellman_c_t) ) != convey_OK){return(-1);}

  int64_t cnt = 0;
  int64_t numpushed = 0;
  double t1 = wall_seconds();

  pkg_bellman_c_t pkg;
  int64_t k,kk, pe;
  int64_t l_i, I, J;
  int64_t loop;

  for(loop=0; loop<numrows; loop++){
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      I = l_i * THREADS + MYTHREAD;
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        convey_relax_edge(tent, dmat,  I, mat->nonzero[k], mat->value[k]);
      }
    }
    lgp_barrier();
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  return(stat->avg);
}

convey_relax_edge(tent, dmat, I, J, weight_IJ )
{

}
#if 0
            break;
          numpushed++;
          if(convey_push(conv, &pkg, pe) != convey_OK){
            tri_convey_push_process(&cnt, conv, L, 0); 
            kk--;
            numpushed--;
          }
        }
      }
    }
    while ( tri_convey_push_process(&cnt, conv, L, 1) ) // keep popping til all threads are done
      ;
  }else{
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
          pkg.w = U->lnonzero[kk]; 
          numpushed++;
          if(convey_push(conv, &pkg, pe) != convey_OK){
            tri_convey_push_process(&cnt, conv, U, 0); 
            kk--;
            numpushed--;
          }
        }
      }
    }
    while ( tri_convey_push_process(&cnt, conv, U, 1) ) // keep popping til all threads are done
      ;
  }

}
#endif
