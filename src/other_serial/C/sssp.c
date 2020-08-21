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
 * Dijsktra's Algorithm
 * the Delta-stepping algorithm for 
 */

#include "spmat_utils.h"

double sssp_dijsktra_linear(sparsemat_t * mat, double *dist, int64_t v0);
double sssp_dijsktra_heap(sparsemat_t * mat, double *dist, int64_t r0);
double sssp_bellmanford_dp(sparsemat_t *dmat, double *dist, int64_t r0);
double sssp_bellmanford_one(sparsemat_t *dmat, double *dist, int64_t r0);
double sssp_delta_stepping(sparsemat_t *dmat, double *dist, int64_t r0);

int main(int argc, char * argv[]) 
{
  #define NUMROWS 20 
  int64_t numrows=NUMROWS;
  double edge_prob = 0.25;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  enum MODEL {GENERIC_Model=1, DIJSKTRA_HEAP=2, DELTA_STEPPING=4, BELLMAN=8, BELLMAN_ONE=16, ALL_Models=32};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;
  int64_t i;

  sparsemat_t *mat, *dmat;
 
  int64_t dump_files = 1;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%"PRId64 ,&numrows ); break;
    case 's': sscanf(optarg,"%d" ,&seed ); break;
    case 'e': sscanf(optarg,"%lg", &edge_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
#if 0 // test the heap stuff
    PQ_t * pq = init_pqueue(numrows);
    for(i=1; i<numrows; i++) {
      pq->val[i] = (double)(numrows-i);
      pq->row[i] = i;
      pq->node[i] = i;
    }
    pq->tail = numrows;
    print_queue(pq);

    heapify_pqueue(pq);
    exit(1);
#endif



  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: sssp: read graph from %s Failed\n", filename); exit(1);}
  } else {
    //mat = erdos_renyi_tri(numrows, er_prob, ER_TRI_L, seed);
                //graph_model model = FLAT;
    // JGD: SHOULD WE ALLOW BOTH UNDIRECTED and DIRECTED? 
    mat = random_graph(numrows, FLAT, UNDIRECTED_WEIGHTED, NOLOOPS, edge_prob, seed);
    if(!mat){
      printf("ERROR: triangles: erdos_renyi_graph Failed\n"); 
      exit(1);
    }
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of sssp\n");
    fprintf(stderr,"Number of rows    (-n) %"PRId64"\n", numrows);
    fprintf(stderr,"random seed     (-s)= %d\n", seed);
    fprintf(stderr,"erdos_renyi_prob   (-e)= %lg\n", edge_prob);
    fprintf(stderr,"models_mask     (-M)= %d\n", models_mask);
    fprintf(stderr,"readgraph      (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files      (-D)=%"PRId64"\n", dump_files);
    fprintf(stderr,"quiet        (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }

  // why???
  for(i=0; i< mat->nnz; i++) {
    //  mat->value[i] = floor( 20 * mat->value[i]);
  }
  // mat is the lower triangle of the adjacency matrix
  // we use the full symmetric matrix in Dijsktra's algorithm
  dmat = make_symmetric_from_lower(mat);
  // debug info
  if(dump_files) {
    dump_matrix(mat,20, "L.out");
    dump_matrix(dmat,20, "Full.out");
  }


  double * compdist = NULL;
  double * dist = malloc(numrows*sizeof(double));
  for(i=0;i<numrows; i++) dist[i] = INFINITY;
  
  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic          sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_dijsktra_linear(dmat, dist, 0);
      compdist = calloc(numrows, sizeof(double));
      for(i=0; i<numrows; i++){
        compdist[i] = dist[i];
        printf("%"PRId64" %g\n",i, dist[i]);
      }
      break;
    case DIJSKTRA_HEAP:
      if( !quiet ) printf("Dijsktra Heap    sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_dijsktra_heap(dmat, dist, 0);
      if(compdist == NULL){
        compdist = calloc(numrows, sizeof(double));
        for(i=0; i<numrows; i++) compdist[i] = dist[i];
      }else{
        for(i=0; i<numrows; i++) assert(compdist[i] == dist[i]);
      }
      printf("Dijkstra Heap: Success!\n");
      break;
    case DELTA_STEPPING:
      if( !quiet ) printf("Delta Stepping   sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_delta_stepping(dmat, dist, 0);
      if(compdist == NULL){
        compdist = calloc(numrows, sizeof(double));
        for(i=0; i<numrows; i++) compdist[i] = dist[i];
      }else{
        for(i=0; i<numrows; i++) assert(compdist[i] == dist[i]);
      }
      printf("Delta Stepping: Success!\n");
      break;
    case BELLMAN:
      if( !quiet ) printf("Bellman Ford dp  sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_bellmanford_dp(dmat, dist, 5);
      for(i=0; i<numrows; i++)
        printf("%ld %g\n",i, dist[i]);
      break;
    case BELLMAN_ONE:
      if( !quiet ) printf("Bellman Ford one sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_bellmanford_one(dmat, dist, 5);
      for(i=0; i<numrows; i++)
        printf("%ld %g\n",i, dist[i]);
      break;
    default:
      continue;
    }
  }
  
  free(dist);
  free(compdist);
  
  printf("done: %g\n", laptime);
  clear_matrix(mat);
  clear_matrix(dmat);
  return(0);
}
