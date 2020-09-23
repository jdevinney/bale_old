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

/*! \file sssp.h
 * \brief Demo application that implements Single Source Shortest Path algorithms.
 * 
 */

#ifndef sssp_INCLUDED
#define sssp_INCLUDED
#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <spmat.h>
#include <math.h>
#include <locale.h>

double        sssp_bellman_agi(d_array_t *tent, sparsemat_t * mat, int64_t v0);
double    sssp_bellman_exstack(d_array_t *tent, sparsemat_t * mat, int64_t v0);
double   sssp_bellman_exstack2(d_array_t *tent, sparsemat_t * mat, int64_t v0);
double     sssp_bellman_convey(d_array_t *tent, sparsemat_t * mat, int64_t v0);
double      sssp_delta_exstack(d_array_t *tent, sparsemat_t * mat, int64_t v0, double opt_delta);
double      sssp_delta_exstack2(d_array_t *tent, sparsemat_t * mat, int64_t v0, double opt_delta);
double       sssp_delta_convey(d_array_t *tent, sparsemat_t * mat, int64_t v0, double opt_delta);

void dump_tent(char *str, d_array_t *tent);

// The exstack and conveyor implementations all use the same package struct:

typedef struct sssp_pkg_t {
  //int64_t i; // We don't build the tree of paths from the vertices back along the shortest path to vertex 0.
               // If that were required, we would have to send the tail of the edge be relaxed.
               // This would not change any patterns, only increase the bandwidth demand.
  int64_t lj;  // the local "head" of the edge
  double tw;   // new tentative weight
}sssp_pkg_t ;

// alternates go here

#endif
