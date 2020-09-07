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

static int relax_edge( double *tent_head, double *tent_tail, double edge_wt )
{
   if( *tent_head > *tent_tail + edge_wt ){
     //printf("relaxing %lg to %lg\n", *tent_head, (*tent_tail + edge_wt));
     *tent_head = *tent_tail + edge_wt;
     return(1);
   }
   return(0);
}


/*!
 * \brief This routine implements the most naive version of the Bellman-Ford algorithm
 *
 * \param tent array to hold the tentative distances
 * \param mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 * \return runtime
 */
double sssp_bellmanford_simple(d_array_t *tent, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i, k, loop;
  int64_t numrows = mat->numrows;
  
  for(i=0; i<numrows; i++)
    tent->entry[i] = INFINITY;
  tent->entry[v0] = 0.0;

  for(loop=0; loop<numrows; loop++){
    for(i=0; i<numrows; i++){ 
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        relax_edge( &(tent->entry[ mat->nonzero[k] ]), &(tent->entry[i]), mat->value[k]);
      }
    }
  }

  return(wall_seconds() - tm);
}



/*!
 * \brief This routine implements the textbook version of 
 * Bellman-Ford as a Dynamic Programming algorithm.
 *
 * \param dist array to hold the tentative distances
 * \param mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 */
double sssp_bellmanford_dynprog(d_array_t *dist, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i,j,k, loop;
  int64_t numrows = mat->numrows;
  int64_t changed;

  double ** tent =  (double**) malloc( numrows * sizeof(double*) );
  for(i=0; i<numrows; i++)
    tent[i] =  (double*) malloc( numrows * sizeof(double) );

  loop = 0;
	for(i=0; i<numrows; i++){
		tent[loop][i] = INFINITY;
  }
  tent[loop][v0] = 0.0;
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}

  loop = 1;
	for(i=0; i<numrows; i++){
		tent[loop][i] = tent[loop-1][i];
  }
  for(k = mat->offset[v0]; k < mat->offset[v0+1]; k++){
    j = mat->nonzero[k];
    tent[loop][j] = mat->value[k];
  }
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}

  for(loop=2; loop<numrows; loop++){
    changed = 0;
    for(i=0; i<numrows; i++)
      tent[loop][i] = tent[loop-1][i];
    for(i=0; i<numrows; i++){
      if( tent[loop-2][i] == tent[loop-1][i] )
        continue;
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        changed |= relax_edge( &(tent[loop][mat->nonzero[k]]), &(tent[loop-1][i]), mat->value[k] );
      }
    }
    if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}
    if(changed == 0){
      break;
    }
  }
  assert( loop <= numrows );
	for(i=0; i<numrows; i++)
	  dist->entry[i] = tent[loop][i];

  for(i=0; i<numrows; i++)
    free(tent[i]);
  free(tent);

  return(wall_seconds() - tm);
}

/*!
 * \brief This routine implements a version of 
 * Bellman-Ford as a Dynamic Programming algorithm.
 *
 * \param dist array to hold the tentative distances
 * \param mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 */
double sssp_bellmanford(d_array_t *dist, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i,k, loop;
  int64_t numrows = mat->numrows;
  int64_t changed;
  double *tent_old, *tent_cur, *tent_new, *tent_temp;

  double *tent0 =  (double*) malloc( numrows * sizeof(double) );
  double *tent1 =  (double*) malloc( numrows * sizeof(double) );
  double *tent2 =  (double*) malloc( numrows * sizeof(double) );
  if( tent0 == NULL || tent1 == NULL || tent2 == NULL) return(-1.0);

  loop = 0;
	for(i=0; i<numrows; i++){
		tent0[i] = INFINITY;
  }
  tent0[v0] = 0.0;
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent0[i]); printf("\n");}

  loop = 1;
	for(i=0; i<numrows; i++){
		tent1[i] = tent0[i];
  }
  for(k = mat->offset[v0]; k < mat->offset[v0+1]; k++){
    tent1[mat->nonzero[k]] = mat->value[k];
  }
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent1[i]); printf("\n");}

  tent_old = tent0;
  tent_cur = tent1;
  tent_new = tent2;
  for(loop=2; loop<numrows; loop++){
    changed = 0;
    for(i=0; i<numrows; i++)
      tent_new[i] = tent_cur[i];
    for(i=0; i<numrows; i++){
      if( tent_old[i] == tent_cur[i] )
        continue;
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        changed |= relax_edge( &(tent_new[mat->nonzero[k]]), &(tent_cur[i]), mat->value[k] );
      }
    }
    if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent_new[i]); printf("\n");}
    if(changed == 0){
      break;
    }
    tent_temp = tent_old;
    tent_old = tent_cur;
    tent_cur = tent_new;
    tent_new = tent_temp;
  }
  assert( loop <= numrows );
	for(i=0; i<numrows; i++)
	  dist->entry[i] = tent_new[i];

  free(tent0);
  free(tent1);
  free(tent2);

  return(wall_seconds() - tm);
}

