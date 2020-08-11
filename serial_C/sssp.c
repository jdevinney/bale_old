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

/*! \file sssp.c
 * \brief Demo application that implements different 
 * Single Source Shortest Path alogrithms.
 *
 * versions of
 * the Delta-stepping algorithm for 
 */

#include "spmat_utils.h"

sparsemat_t * LUmat_from_L(sparsemat_t * L)
{
	int64_t i, l, k, pos;
	sparsemat_t *U  = transpose_matrix(L);
	sparsemat_t *retmat = init_matrix(L->numrows, L->numrows, 2*(L->nnz), 1);

	pos = 0;
	retmat->offset[0] = 0;
	for(i=0; i<L->numrows; i++){
    for(l=L->offset[i]; l<L->offset[i+1]; l++){
			retmat->nonzero[pos] = L->nonzero[l];
			retmat->value[pos] = L->value[l];
			pos++;
		}
    for(k=U->offset[i]; k<U->offset[i+1]; k++){
			retmat->nonzero[pos] = U->nonzero[k];
			retmat->value[pos] = U->value[k];
			pos++;
		}
		retmat->offset[i+1] = pos;
	}
	clear_matrix(U);
	free(U);
  return(retmat);	
}

/*!
 * \brief This routine implements generic serial version of Dijkstra's algorithm
 * \param *L sparsemat_t that holds the graph as a weighted lower triangular matrix
 * \param *U sparsemat_t that holds the upper triangular part of the 
 *           symmetric adjacency matrix
 * \return run time
 * The algorithm requires that the weights are positive.
 * The alogrithm modifies the array dist[], which hold the best know weight of a path
 * from the given vertex to all other vertices. 
 */
double sssp_generic(sparsemat_t * L, int64_t v0)
{
	int64_t i, k, vn;
	int64_t numrows = L->numrows;
	double * dist = malloc(numrows*sizeof(double));
	double minwt, curwt;
	int64_t minidx;

  sparsemat_t *LU = LUmat_from_L(L);
	dump_matrix(LU,20, "LU.out");

	for(i=0; i<numrows; i++)
		dist[i] = INFINITY;
	for(k = LU->offset[v0];  k < LU->offset[v0+1]; k++)
		dist[ LU->nonzero[k] ] = LU->value[k];
  dist[v0] = -0.0;
	printf("dist[v0] = %lg\n", dist[v0]);

	while( 1 ) {
		// find the smallest tentative distance
	  minwt = INFINITY;
		minidx = numrows;
	  for(i=0; i<numrows; i++){
      if( dist[i] > 0 && dist[i] < minwt ){
			  minwt = dist[i];
			  minidx = i;
		  }
		}
		if(minidx == numrows)    // done: all distances have be computed
			break;
		if(minwt == INFINITY)    // fail?: graph is not connected
			break;

    curwt = dist[minidx];
		
		for(k = LU->offset[minidx];  k < LU->offset[minidx+1]; k++){
			vn = LU->nonzero[k];
      if( dist[vn] > curwt + L->value[k] )
				dist[vn] = curwt + L->value[k];
		}
		printf("dist[%ld] = %lg\n", minidx, curwt);
    dist[minidx] = -dist[minidx];
  }

	for(i=0; i<numrows; i++)
		dist[i] = -dist[i];

	for(i=0; i<numrows; i++)
		printf("%ld %g\n",i, dist[i]);

	return(0.0);
}

int main(int argc, char * argv[]) 
{
  #define NUMROWS 20 
  int64_t numrows=NUMROWS;
  double edge_prob = 0.25;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  enum MODEL {GENERIC_Model=1, ALL_Models=2};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;

  sparsemat_t *mat, *tmat;
 
  int64_t dump_files = 1;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows ); break;
    case 's': sscanf(optarg,"%d" ,&seed ); break;
    case 'e': sscanf(optarg,"%lg", &edge_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: sssp: read graph from %s Failed\n", filename); exit(1);}
  } else {
    //mat = erdos_renyi_tri(numrows, er_prob, ER_TRI_L, seed);
		//graph_model model = FLAT;
    mat = random_graph(numrows, FLAT, UNDIRECTED_WEIGHTED, NOLOOPS, edge_prob, seed);
    if(!mat){printf("ERROR: triangles: erdos_renyi_graph Failed\n"); exit(1);}
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of sssp\n");
    fprintf(stderr,"Number of rows    (-n) %ld\n", numrows);
    fprintf(stderr,"random seed     (-s)= %d\n", seed);
    fprintf(stderr,"erdos_renyi_prob   (-e)= %lg\n", edge_prob);
    fprintf(stderr,"models_mask     (-M)= %d\n", models_mask);
    fprintf(stderr,"readgraph      (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files      (-D)=%ld\n", dump_files);
    fprintf(stderr,"quiet        (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }

	// mat is the lower triangle of the adjacency matrix
  tmat = transpose_matrix(mat); 
	// debug info
  if(dump_files) {
		dump_matrix(mat,20, "L.out");
		dump_matrix(tmat,20, "U.out");
	}

	int64_t * deg = malloc(numrows * sizeof(int64_t));
	int64_t i;
	for(i=0; i<numrows; i++)
		deg[i] =0;
	for(i=0; i<numrows; i++){
		deg[i] += mat->offset[i+1] - mat->offset[i];
		deg[i] += tmat->offset[i+1] - tmat->offset[i];
	}
	// worry about the graph being non-connected
  printf("Deg: ");
	for(i=0; i<numrows; i++)
		printf(" %ld, ", deg[i]);
  printf("\n");

  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic          sssp:\n");
      laptime = sssp_generic(mat, 0);
      break;
    default:
      continue;
    }
  }

	printf("done: %g\n", laptime);
  clear_matrix(mat);
  return(0);
}
