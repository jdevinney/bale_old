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
 * \brief Demo application that calls different 
 * Single Source Shortest Path alogrithms.
 *
 * versions of
 * Dijsktra's Algorithm
 * the Delta-stepping algorithm for 
 */

#include "spmat_utils.h"
#include "std_options.h"

double sssp_dijsktra_linear(d_array_t * tent, sparsemat_t * mat, int64_t v0);
double sssp_dijsktra_heap(d_array_t * tent, sparsemat_t * mat, int64_t r0);
double sssp_bellmanford_simple(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_bellmanford_dynprog(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_bellmanford(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_delta_stepping(d_array_t * tent, sparsemat_t *mat, int64_t r0, double del);
double sssp_answer_diff(d_array_t *A, d_array_t *B);



/*!
 * \brief Compare two arrays
 * \param *A one array (vector)
 * \param *B the other
 * \return the l_2 norm of the given arrays
 */
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

typedef struct args_t{
  double deltaStep; 
  int64_t V0;
  int64_t alg;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'a': args->alg = atoi(arg); break;
    case 'S': args->deltaStep = atof(arg); break;     
    case 'V': args->V0 = atoi(arg); break;
    case ARGP_KEY_INIT:
      args->deltaStep = 0.0;
      args->V0 = 0;
      args->alg = 3;
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
{
    {"alg", 'a', "flag", 0, "alg: 1==bellman | 2==delta"},
    {"deltaStep", 'S', "STEPSIZE", 0, "user supplied delta step size"},  
    {"V0", 'V', "NUM", 0, "initial vertex"},  
    {0}
};

static struct argp_child children_parsers[] =
{    
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
};


int main(int argc, char * argv[]) 
{
  args_t args ={0};  
  struct argp argp = {options, parse_opt, 0, "C version of Single Source Shortest Path algorithms.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  
  double nz_per_row = args.gstd.nz_per_row;
  double edge_prob = args.gstd.edge_prob;
  int64_t numrows = args.gstd.numrows;
  edge_type edge_type = UNDIRECTED;
  self_loops loops = LOOPS;
  int quiet = args.std.quiet;

  graph_model model = args.gstd.model;
  
  if(args.gstd.readfile == 0){
    resolve_edge_prob_and_nz_per_row(&edge_prob, &nz_per_row, numrows, edge_type, loops);
  }
  
#if 0
  if(!quiet ) {
    fprintf(stderr,"Running C versions of SSSP\n");
    if(args.gstd.readfile == 1)
      fprintf(stderr,"Reading a matrix from file (-f [%s])\n", args.gstd.filename);
    else{
      if(args.gstd.model == FLAT)
        fprintf(stderr,"flat model           (-F)\n");
      else        
        fprintf(stderr,"geometric model      (-G)\n");
      fprintf(stderr,"Number of rows       (-n) %"PRId64"\n", numrows);
      fprintf(stderr,"edge_density         (-e)= %lg\n", edge_prob);
      fprintf(stderr,"nz_per_row           (-z)= %lg\n", nz_per_row);
      fprintf(stderr,"random seed          (-s)= %ld\n",  args.std.seed);
    }
    fprintf(stderr,"models_mask          (-M)= %d\n", args.std.models_mask);
    fprintf(stderr,"dump_files           (-D)= %d\n", args.std.dump_files);
    fprintf(stderr,"---------------------------------------\n");
  }
#endif

  sparsemat_t *mat;
  if(args.gstd.readfile) {
    mat = read_matrix_mm(args.gstd.filename);
    if(!mat){printf("ERROR: sssp: read graph from %s Failed\n", args.gstd.filename); exit(1);}
  } else {
    mat = random_graph(numrows, model, DIRECTED_WEIGHTED, NOLOOPS, edge_prob, args.std.seed);
    if(!mat){ printf("ERROR: sssp: erdos_renyi_graph Failed\n"); exit(1); }
  }

  if(!quiet){
    printf("Input matrix stats:\n");
    spmat_stats(mat);
    fprintf(stderr,"---------------------------------------\n");
  }

  if(args.std.dump_files){
    dump_matrix(mat, 20, "mat.out");
    write_matrix_mm(mat, "ssspout.mm");
  }

  uint32_t use_model;
  enum MODEL {GENERIC_Model=1, DIJSKTRA_HEAP=2, DELTA_STEPPING=4, BELLMAN_SIMPLE=8, BELLMAN=16, ALL_Models=32};
  int64_t models_mask = args.std.models_mask;
  models_mask=ALL_Models - 1;
 
  double laptime = 0.0;
  char model_str[32];

  d_array_t *comp_tent=NULL;
  d_array_t *tent = init_d_array(numrows);

  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    model_str[0] = '\0';
    switch( use_model & models_mask ){
    case GENERIC_Model:
      sprintf(model_str, "Dijsktra Linear");
      set_d_array(tent, INFINITY);
      laptime = sssp_dijsktra_linear(tent, mat, 0);
      break;

    case DIJSKTRA_HEAP:
      sprintf(model_str, "Dijsktra Heap");
      set_d_array(tent, INFINITY);
      laptime = sssp_dijsktra_heap(tent, mat, 0);
      break;

    case DELTA_STEPPING:
      sprintf(model_str, "Delta Stepping");
      set_d_array(tent, INFINITY);
      laptime = sssp_delta_stepping(tent, mat, 0, 0.0);
      break;

    case BELLMAN_SIMPLE:
      sprintf(model_str, "Bellman Ford Simple");
      set_d_array(tent, INFINITY);
      laptime = sssp_bellmanford_simple(tent, mat, 0);
      break;

    case BELLMAN:
      sprintf(model_str, "Bellman Ford DP");
      set_d_array(tent, INFINITY);
      laptime = sssp_bellmanford(tent, mat, 0);
      break;

    default:
      continue;
    }
    if(model_str[0]) {
      if(comp_tent == NULL){
        comp_tent = copy_d_array(tent);
        sprintf(model_str, "%s ()", model_str);
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          sprintf(model_str, "%s (compares)", model_str);
      }
      fprintf(stderr, "%30s: %8.3lf\n", model_str, laptime);
    }
  }

  if(args.std.dump_files){
    write_d_array(comp_tent, "ssspout.wts");
  }
  
  clear_matrix(mat); free(mat);
  clear_d_array(tent); free(tent);
  clear_d_array(comp_tent); free(comp_tent);
  return(0);
}
