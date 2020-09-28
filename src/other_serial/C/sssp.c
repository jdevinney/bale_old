/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file sssp.c
 * \brief The application that calls different Single Source Shortest Path alogrithms.
 *
 * It calls several versions of Dijsktra's algorithm, the Bellman Ford algorithm, and
 * the Delta-stepping algorithm.
 */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"

double sssp_dijsktra_linear(d_array_t * tent, sparsemat_t * mat, int64_t v0);
double sssp_dijsktra_heap(d_array_t * tent, sparsemat_t * mat, int64_t r0);
double sssp_bellmanford_simple(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_bellmanford_dynprog(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_bellmanford(d_array_t * tent, sparsemat_t *mat, int64_t r0);
double sssp_delta_stepping_ptr(d_array_t * tent, sparsemat_t *mat, int64_t r0, double del);
double sssp_delta_stepping_arr(d_array_t * tent, sparsemat_t *mat, int64_t r0, double del);
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

/********************************  argp setup  ************************************/
typedef struct args_t{
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] = {{0}};

static struct argp_child children_parsers[] =
{
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
};


int main(int argc, char * argv[]) 
{
  enum MODEL {GENERIC_Model=1, DIJSKTRA_HEAP=2, DELTA_STEPPING_PTR=4, DELTA_STEPPING_ARR=8, BELLMAN_SIMPLE=16, BELLMAN=32, ALL_Models=64};
  args_t args;  
  struct argp argp = {options, parse_opt, 0, "SSSP for a weighted graph.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  args.gstd.numrows = SSSP_NUM_ROWS;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  //override command line 
  //(note:these will lead to matrices with not quite the right number of nonzeros 
  // if the user also used the -z flag.)
  if (args.gstd.loops == 1) {
    fprintf(stderr,"WARNING: toposort requires 1s on the diagonal.\n");
    args.gstd.loops = 0;
  }
  if (args.gstd.directed == 0) {
    fprintf(stderr,"WARNING: toposort starts with an upper triangalur matrix.\n");
    args.gstd.directed = 1;
  }
  if (args.gstd.weighted == 0) {
    fprintf(stderr,"WARNING: toposort starts with an upper triangalur matrix.\n");
    args.gstd.weighted = 1;
  }

  write_std_graph_options(&args.std, &args.gstd);
  write_std_options(&args.std);
  
  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: SSSP: mat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(L, "sssp_inmat");

#if 0
  
  double nz_per_row = args.gstd.nz_per_row;
  double edge_prob = args.gstd.edge_prob;
  int64_t numrows = args.gstd.numrows;
  edge_type edge_type = UNDIRECTED;
  self_loops loops = LOOPS;
  int quiet = args.std.quiet;

  graph_model model = args.gstd.model;
  int64_t models_mask = args.std.models_mask;
  models_mask=ALL_Models - 1;
  
  if(args.gstd.readfile == 0){
    resolve_edge_prob_and_nz_per_row(&edge_prob, &nz_per_row, numrows, edge_type, loops);
  }
#endif
  
#if 1                 // TODO
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
#endif

#if 0
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
#endif
  double laptime = 0.0;
  uint32_t use_model;
 
  d_array_t *tent, *comp_tent=NULL;
  tent = init_d_array(numrows);
  set_d_array(tent, INFINITY);

  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & args.std.models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic          sssp: ");
      laptime = sssp_dijsktra_linear(tent, mat, 0);
      comp_tent = init_d_array(numrows);
      copy_d_array(comp_tent, tent);
      printf("n-squared Dijkstra Heap run on default!\n");
      break;

    case DIJSKTRA_HEAP:
      if( !quiet ) printf("Dijsktra Heap    sssp: ");
      laptime = sssp_dijsktra_heap(tent, mat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Dijkstra Heap: nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Dijkstra Heap: compares successfully!\n");
      }
      break;

    case DELTA_STEPPING_PTR:
      if( !quiet ) printf("Delta Stepping ptr   : ");
      laptime = sssp_delta_stepping_ptr(tent, mat, 0, 0.0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Delta Stepping: nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Delta Stepping: compares successfully!\n");
      }
      
      break;

    case DELTA_STEPPING_ARR:
      if( !quiet ) printf("Delta Stepping arr   : ");
      laptime = sssp_delta_stepping_arr(tent, mat, 0, 0.0);
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
      laptime = sssp_bellmanford_simple(tent, mat, 0);
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
      laptime = sssp_bellmanford(tent, mat, 0);
      if(comp_tent == NULL){
        comp_tent = init_d_array(numrows);
        copy_d_array( comp_tent, tent);
        printf("Bellman-Ford  DP  : nothing to compare to!\n");
      }else{
        if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
          printf("Bellman-Ford  DP  : compares successfully!\n");
      }
      break;

    default:
      continue;
    }
    if(!quiet) printf("%8.3lf seconds\n", laptime);
  }

  if(args.std.dump_files){
    write_d_array(comp_tent, "ssspout.wts");
  }
  
  clear_matrix(mat); free(mat);
  clear_d_array(tent); free(tent);
  clear_d_array(comp_tent); free(comp_tent);
  return(0);
}
