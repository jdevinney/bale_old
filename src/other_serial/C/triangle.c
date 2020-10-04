/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file triangle.c
 * \brief Program that counts the number of triangles in a graph 
 * given by its adjacency matrix
 */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"


/* \page triangle_page triangle
Count the triangles in a graph simple graph (a lower triangular matrix)
*/

/*!
 * \brief This routine counts the number of triangles in a graph
 *   given the lower triangular piece of the adjacency matrix
 * \param *triangles a place to write the number of triangles found
 * \param *mat the sparse matrix that holds the graph
 * \return run time
 */
double triangles_matrix(int64_t *triangles, sparsemat_t *mat) 
{
  int64_t j, k, l, numtriangles;
  int64_t u,v;
  numtriangles = 0;

  double t1 = wall_seconds();
  
  // for each non-zero (i,j) in L accumulate the size 
  // of the intersection of row_i and row_j.
  // Note: Because the matrix is tidy,
  // (the columns in each row appear in increasing order)
  // we can find the intersection in a single pass over both rows.

  for(u = 0; u < mat->numrows; u++){ 
    for(j = mat->offset[u]; j < mat->offset[u+1]; j++){
      v = mat->nonzero[j];
      for( l = mat->offset[u], k = mat->offset[v];  k < mat->offset[v+1] && l < mat->offset[u+1];  ){
        if( mat->nonzero[k] == mat->nonzero[l] ){
          numtriangles++;
          k++;
          l++;
        }else if( mat->nonzero[k] > mat->nonzero[l] ){
          l++;
        }else{ // ( mat->nonzero[u] > mat->nonzero[W] )
          k++;
        }
      }
    }
  }

  t1 = wall_seconds() - t1;
 
  *triangles = numtriangles;
  return(t1);
}


/********************************  argp setup  ************************************/
typedef struct args_t{
  int alg;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
  case 'a': args->alg = atoi(arg); break;     
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    state->child_inputs[1] = &args->gstd;
    break;
  }
  return(0);
}

static struct argp_option options[] =
{
  {"triangle_alg", 'a', "ALG", 0, "Algorithm: 0 means L&L*U, 1 means L&U*L"},  
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
  args_t args = {0};
  struct argp argp = {options, parse_opt, 0, "Triangle counting", children_parsers};  
  args.gstd.numrows = 500;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  //override command line (these will lead to matrices with not quite the right number of nonzeros
  // if the user also used the -z flag.
  if ( (args.gstd.loops == 1) || (args.gstd.directed == 1) ) {
    fprintf(stderr,"WARNING: triangles counting requires undirected no-loops graph\n");
    args.gstd.loops = 0;
    args.gstd.directed = 0;
  }

  write_std_graph_options(&args.std, &args.gstd);
  write_std_options(&args.std);
  
  // read in a matrix or generate a random graph
  sparsemat_t * L = get_input_graph(&args.std, &args.gstd);
  if(!L){fprintf(stderr, "ERROR: triangle: L is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(L, "triangle_inmat");

  /* make sure the matrix is in legal form */
  if (args.gstd.readfile) {
    ret = tril(L, -1);
    if (ret)
      fprintf(stderr,"WARNING: input graph was not lower-triangular with zero diagonal. Removing illegal nonzeros.\n");
  } else if (!is_lower_triangular(L, 0)) {
    fprintf(stderr,"ERROR: L is not lower triangular!\n");
    exit(1);  
  }  


  // if KRONECKER, calculate the number of triangles 
  int64_t correct_answer = -1;
  int wrote_num_triangles = 0;
  if (args.gstd.model == KRONECKER) {
    correct_answer = tri_count_kron_graph(args.gstd.kron_mode, args.gstd.kron_spec, args.gstd.kron_num);
    bale_app_write_int(&args.std, "num_triangles", correct_answer);
    wrote_num_triangles = 1;
  }

    
  enum FLAVOR {GENERIC=1, ALL_Models=2};
  uint32_t use_model;
  double laptime = 0.0;
  int64_t tri_cnt;
  char model_str[1024];
  int models_mask = (args.std.models_mask) ? args.std.models_mask : 3;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    tri_cnt = 0;
    switch( use_model & models_mask ){
    case GENERIC:
      sprintf(model_str, "L .& (L * U)");
      laptime = triangles_matrix(&tri_cnt, L);
      break;
    default:
      continue;
    }
    printf(" %12"PRId64" triangles,  %8.3lf seconds\n", tri_cnt, laptime);
  }

  if(!wrote_num_triangles){
    bale_app_write_int(&args.std, "triangles", tri_cnt);
    wrote_num_triangles = 1;
  }
    
  bale_app_write_time(&args.std, model_str, laptime);    

  clear_matrix(L);

  return(0);
}

