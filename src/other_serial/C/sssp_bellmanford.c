/******************************************************************
//
//
// Copyright(C) 2018, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// This material may be reproduced by or for the US Government
// pursuant to the copyright license under the clauses at DFARS
// 252.227-7013 and 252.227-7014.
// 
//
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//  * Neither the name of the copyright holder nor the
//   names of its contributors may be used to endorse or promote products
//   derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.
// 
*****************************************************************/ 

#include "spmat_utils.h"

/*! \file sssp_bellmanford.c
 * \brief These are the function that implement 
 * and support two implementations of Dijsktra's Alg.
 * 
 * It is assumed that the Algorithm is known.
 * Excellent description of the algorithm are easily found.
 * To fix our terminology we say that the algorithm expands 
 * paths to "unvisited" vertices and assigns them tentative weights
 * and once the shortest path is found the vertex is said to be resolved.
 *   We assume:
 *     -- that the graph is weighted with positive weights
 *     -- if the graph is undirected, then mat contains two "directed edges" for
 *        each undirected edge.
 *
 */


/*!
 * \brief This routine implements the most naive version of Dijkstra's algorithm
 *
 * This first implementation uses a simple array to hold the
 * tentative weights and uses tricks with the weights to flag
 * the transition of vertices from unvisited to tentative to resolved.
 *
 * \param *mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 * \return runtime
 */
double sssp_bellmanford_dp(sparsemat_t * mat, double *dist, int64_t v0)
{
	int64_t i, j, k, loop;
	int64_t numrows = mat->numrows;

  double * d1 =  (double*) malloc( numrows * sizeof(double) );
  double * d2 =  (double*) malloc( numrows * sizeof(double) );
  double *oldd, *newd, *swap;
  oldd = d1;
  newd = d2;
	for(i=0; i<numrows; i++)
		oldd[i] = INFINITY;
  oldd[v0] = 0.0;

  for(loop=0; loop<numrows; loop++){
	  for(i=0; i<numrows; i++)
	  	newd[i] = oldd[i];
    for(i=0; i<numrows; i++){ 
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
			  j = mat->nonzero[k];
        if(newd[j] > oldd[i] + mat->value[k])
          newd[j] = oldd[i] + mat->value[k];
      }
    }
    swap = oldd;
    oldd = newd;
    newd = swap;
  }

	for(i=0; i<numrows; i++)
		dist[i] = oldd[i];


	return(0.0);
}


/*!
 * \brief This routine implements the most naive version of Dijkstra's algorithm
 *
 * This first implementation uses a simple array to hold the
 * tentative weights and uses tricks with the weights to flag
 * the transition of vertices from unvisited to tentative to resolved.
 *
 * \param *mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 * \return runtime
 */
double sssp_bellmanford_one(sparsemat_t * mat, double *dist, int64_t v0)
{
	int64_t i, j, k, loop;
	int64_t numrows = mat->numrows;

	for(i=0; i<numrows; i++)
		dist[i] = INFINITY;
  dist[v0] = 0.0;

  for(loop=0; loop<numrows; loop++){
    for(i=0; i<numrows; i++){ 
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
			  j = mat->nonzero[k];
        dist[j] = (dist[j] > dist[i] + mat->value[k]) ? (dist[i] + mat->value[k]) : dist[j];
      }
    }
  }

	return(0.0);
}


