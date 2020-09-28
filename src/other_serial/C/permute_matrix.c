/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file permute_matrix.c
 * \brief Program that runs the permute_matrix routine from the spmat library
 *
 * Run permute_matrix --help or --usage for insructions on running.
 */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"

/*! \page permute_matrix_page Permute a sparse matrix */

/*!
 * \brief A timing wrapper around the permute_matrix call in spmat_utils
 * \param *A  sparsemat_t holding the given matrix
 * \param *rp the row permutation
 * \param *cp the column permutation
 * \return run time
 * NB: This does charge for the time to initialize the permuted matrix.
 */
double permute_matrix_generic(sparsemat_t *A, int64_t *rp, int64_t *cp)
{
  double tm;
  tm = wall_seconds();
  sparsemat_t *Ap = permute_matrix(A, rp, cp);
  tm = wall_seconds() - tm;
  clear_matrix(Ap);
  return(tm);
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

static struct argp_option options[] =
  {
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
  args_t args = {{0}};
  enum FLAVOR {GENERIC=1, ALL_Models=2};
  args.std.models_mask = ALL_Models-1;
  args.gstd.numrows = PERMUTE_NUM_ROWS;
  struct argp argp = {options, parse_opt, 0, "Permute rows and columns of a sparse matrix", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_graph_options(&args.std, &args.gstd);
  write_std_options(&args.std);
  
  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: triangle: L is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(mat, "permute_inmat");

#if 0
  double nz_per_row = args.gstd.nz_per_row;
  double edge_prob = args.gstd.edge_prob;
  int64_t numrows = args.gstd.numrows;
  edge_type edge_type = DIRECTED;
  self_loops loops = NOLOOPS;
  
  if(!args.gstd.readfile){
    resolve_edge_prob_and_nz_per_row(&edge_prob, &nz_per_row, numrows, edge_type, loops);
  }

  if(!args.std.quiet ) {
    fprintf(stderr,"Running C version of permute matrix\n");
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
  }
  
  sparsemat_t *mat;
  if( args.gstd.readfile ) {
    mat = read_matrix_mm(args.gstd.filename);
    if(!mat){printf("ERROR: Read graph from %s Failed\n", args.gstd.filename); exit(1);}
  } else {
    //mat = erdos_renyi_graph(numrows, er_prob, ER_GRAPH_DIRECT, seed);
    mat = random_graph(numrows, args.gstd.model, edge_type, loops, edge_prob, args.std.seed);
    if(!mat){printf("ERROR: erdos_renyi_graph failed!\n"); exit(1);}
  }
  if(args.std.dump_files) dump_matrix(mat, 20, "mat.out");
  

  if(!args.std.quiet) spmat_stats(mat);
  
#endif
  int64_t *rperm = rand_perm(mat->numrows, args.std.seed + 1); assert(rperm != NULL);
  int64_t *cperm = rand_perm(mat->numrows, args.std.seed + 2); assert(cperm != NULL);

  uint32_t use_model;
  double laptime = 0.0;
  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & args.std.models_mask ) {
    case GENERIC:
      printf("generic permute matrix: ");     //TODO model_str
      laptime = permute_matrix_generic(mat, rperm, cperm);
    break;
    default:
       continue;
    }
    printf("  %8.3lf seconds \n", laptime);
  }
// TODO check the result

  return(0);
}

