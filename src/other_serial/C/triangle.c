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

/*! \file triangle.c
 * \brief Demo application that counts the number of triangles in the graph
 * given by it's adjacency matrix.
 */

#include "spmat_utils.h"
#include "std_options.h"


/*! \page triangle_page 
  Count the triangles in a graph presented by a lower triangular matrix 
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
  int64_t U,V;
  numtriangles = 0;

  double t1 = wall_seconds();
  
  // for each non-zero (i,j) in L accumulate the size of the intersection
  // of row_i and row_j.

  for(U = 0; U < mat->numrows; U++){ 
    for(j = mat->offset[U]; j < mat->offset[U+1]; j++){
      V = mat->nonzero[j];
      for( l = mat->offset[U], k = mat->offset[V];  k < mat->offset[V+1] && l < mat->offset[U+1];  ){  // This requires that the matrix be tidy
        if( mat->nonzero[k] == mat->nonzero[l] ){
          numtriangles++;
          k++;
          l++;
        }else if( mat->nonzero[k] > mat->nonzero[l] ){
          l++;
        }else{ // ( mat->nonzero[U] > mat->nonzero[W] ) {
          k++;
        }
      }
    }
  }

  t1 = wall_seconds() - t1;
 
  *triangles = numtriangles;
  return(t1);
}


typedef struct args_t{
  int kronecker;
  char * kstring;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'K':
      args->kronecker = 1; break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"kronecker", 'K', 0, 0, "Specify the input to be a Kronecker Product graph"},
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

  /* process command line */
  args_t args;
  args.kronecker = 0;
  struct argp argp = {options, parse_opt, 0, "Count the number of triangles in a graph.",
                      children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);

  double nz_per_row = args.gstd.nz_per_row;
  double edge_prob = args.gstd.edge_prob;
  int64_t numrows = args.gstd.numrows;
  edge_type edge_type = UNDIRECTED;
  self_loops loops = LOOPS;
  int quiet = args.std.quiet;
  
  if(args.gstd.readfile == 0){
    resolve_edge_prob_and_nz_per_row(&edge_prob, &nz_per_row, numrows, edge_type, loops);
  }
  
  if(!quiet ) {
    fprintf(stderr,"Running C version of toposort\n");
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
      fprintf(stderr,"random seed          (-s)= %"PRId64"\n",  args.std.seed);
    }
    fprintf(stderr,"models_mask          (-M)= %d\n", args.std.models_mask);
    fprintf(stderr,"dump_files           (-D)= %d\n", args.std.dump_files);
  }

  
  enum FLAVOR {GENERIC=1, ALL=2};
  char * kron_list_def = {"M: k1 k2 k3 ..."};
  char * kron_list = kron_list_def;
  kron_args_t * kron_args = NULL;
  uint32_t use_model;
  
  sparsemat_t *mat = NULL;
  int64_t triangles;

  if( args.gstd.readfile ) {
    mat = read_matrix_mm(args.gstd.filename);
    if(!mat){printf("ERROR: triangles: read graph from %s Failed\n", args.gstd.filename); exit(1);}
  } else {
    if(args.kronecker){
      kron_args = kron_args_init(kron_list);
      mat = kronecker_product_graph(kron_args);
      if(!mat){printf("ERROR: triangles: Kronecker Product generation Failed\n"); exit(1);}
    } else {
      mat = random_graph(numrows, args.gstd.model, UNDIRECTED, NOLOOPS, edge_prob, args.std.seed);
      if(!mat){printf("ERROR: triangles: erdos_renyi_graph generation Failed\n"); exit(1);}
    }
  }
  
  if(!quiet){
    printf("Input matrix stats:\n");
    spmat_stats(mat);
    if(args.kronecker )
      printf("Kronecker model should have %"PRId64" triangles\n", tri_count_kron_graph(kron_args));
  }
    
  if(args.std.dump_files){
    dump_matrix(mat, 20, "mat.out");
  }

  double laptime = 0.0;
  for(use_model=1; use_model < ALL; use_model *=2 ){
    triangles = 0;
    switch( use_model & args.std.models_mask ){
    case GENERIC:
      if( !quiet ) printf("Generic      Triangle: ");
        laptime = triangles_matrix(&triangles, mat);
      break;
    default:
      continue;
    }
    if( !quiet ) 
      printf(" %12"PRId64" triangles : %8.3lf seconds\n", triangles, laptime);
  }
 
  if(args.kronecker) 
    free(kron_args);
  clear_matrix(mat);
  return(0);
}

