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

double sssp_dijsktra_linear(d_array_t * tent, sparsemat_t * mat, int64_t v0);
double sssp_dijsktra_heap(d_array_t * tent, sparsemat_t * mat, int64_t r0);
double sssp_bellmanford_simple(d_array_t * tent, sparsemat_t *dmat, int64_t r0);
double sssp_bellmanford_dynprog(d_array_t * tent, sparsemat_t *dmat, int64_t r0);
double sssp_bellmanford(d_array_t * tent, sparsemat_t *dmat, int64_t r0);
double sssp_delta_stepping(d_array_t * tent, sparsemat_t *dmat, int64_t r0);
double sssp_answer_diff(d_array_t *A, d_array_t *B);



double sssp_answer_diff(d_array_t *A, d_array_t *B)
{
  int64_t i;
  double diff = 0.0;

  for(i=0; i<A->num; i++) {
    if( A->entry[i] == INFINITY && B->entry[i] == INFINITY )
      continue;
    diff += (A->entry[i] - B->entry[i]) * (A->entry[i] - B->entry[i]);
  }
  return(sqrt(diff));
}

int main(int argc, char * argv[]) 
{
  #define NUMROWS 20 
  int64_t numrows=NUMROWS;
  double edge_prob = 0.25;
  uint32_t seed = 123456789;
  graph_model model = FLAT;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  enum MODEL {GENERIC_Model=1, DIJSKTRA_HEAP=2, DELTA_STEPPING=4, BELLMAN_SIMPLE=8, BELLMAN=16, ALL_Models=32};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;

  sparsemat_t *dmat;
 
  int64_t dump_files = 1;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:g:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%"PRId64 ,&numrows ); break;
    case 's': sscanf(optarg,"%d" ,&seed ); break;
    case 'e': model = FLAT; sscanf(optarg,"%lg", &edge_prob); break;
    case 'g': model = GEOMETRIC; sscanf(optarg,"%lg", &edge_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
#if 0 // test the heap stuff
    int64_t i;
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
    dmat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: sssp: read graph from %s Failed\n", filename); exit(1);}
  } else {
    dmat = random_graph(numrows, model, DIRECTED_WEIGHTED, NOLOOPS, edge_prob, seed);
    if(!dmat){
      printf("ERROR: sssp: erdos_renyi_graph Failed\n"); 
      exit(1);
    }
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of sssp\n");
    fprintf(stderr,"Number of rows       (-n)= %"PRId64"\n", numrows);
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"Flat edge prob       (-e)= %lg\n", edge_prob);
    fprintf(stderr,"Geometric edge prob  (-g)= %lg\n", edge_prob);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"readgraph            (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files           (-D)= %"PRId64"\n", dump_files);
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }

  // debug info
  if(dump_files) {
    dump_matrix(dmat,20, "Full.out");
  }

#define WRITE_MAT 1
  if(WRITE_MAT){write_matrix_mm(dmat, "ssspout.mm");}

  d_array_t *tent, *comp_tent=NULL;
  tent = init_d_array(numrows);
  set_d_array(tent, INFINITY);

  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic          sssp: ");
      laptime = sssp_dijsktra_linear(tent, dmat, 0);
      comp_tent = init_d_array(numrows);
      copy_d_array(comp_tent, tent);
      printf("n-squared Dijkstra Heap run on default!\n");
      break;
    case DIJSKTRA_HEAP:
      if( !quiet ) printf("Dijsktra Heap    sssp: ");
      laptime = sssp_dijsktra_heap(tent, dmat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Dijkstra Heap: nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Dijkstra Heap: compares successfully!\n");
      }
      break;
    case DELTA_STEPPING:
      if( !quiet ) printf("Delta Stepping   sssp: ");
      laptime = sssp_delta_stepping(tent, dmat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Delta Stepping: nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Delta Stepping: compares successfully!\n");
      }
      
      break;

    case BELLMAN_SIMPLE:
      if( !quiet ) printf("Bellman Ford     sssp: ");
      laptime = sssp_bellmanford_simple(tent, dmat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Bellman-Ford: nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Bellman-Ford (dynamic program): compares successfully!\n");
      }
      break;

    case BELLMAN:
      if( !quiet ) printf("Bellman Ford dp sssp: ");
      laptime = sssp_bellmanford(tent, dmat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Bellman-Ford  DP  : nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Bellman-Ford  DP  : compares successfully!\n");
      }
      break;

    case BELLMAN:
      if( !quiet ) printf("Bellman Ford dp  sssp: ");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_bellmanford_dp(dmat, dist, 0);
      if(compdist == NULL){
        compdist = calloc(numrows, sizeof(double));
        for(i=0; i<numrows; i++) compdist[i] = dist[i];
      }else{
        for(i=0; i<numrows; i++) assert(compdist[i] == dist[i]);
      }
      //printf("Bellman DynProg: Success!\n");
      break;
    default:
      continue;
    }
    if(!quiet) printf("%8.3lf seconds\n", laptime);
  }

  if(WRITE_MAT){write_d_array(comp_tent, "sssp.wts");}
  
  clear_matrix(dmat); free(dmat);
  clear_d_array(tent); free(tent);
  clear_d_array(comp_tent); free(comp_tent);
  return(0);
}
