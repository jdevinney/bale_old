/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
*****************************************************************/ 

/*! \file transpose_matrix.c
 * \brief Demo program that runs the transpose_matrix from the 
 *  spmat library to provide a framework to study other implementations.
 * 
 * Run transpose_matrix --help or --usage for usage.
 */

#include "spmat_utils.h"
#include "std_options.h"

/*! \page transpose_matrix_page Transpose a given sparse matrix */
double transpose_generic(sparsemat_t *A, int64_t dump_files)
{
  double tm;
  if(dump_files)
    write_matrix_mm(A, "er_orig.mm");

  tm = wall_seconds();
  sparsemat_t * At = transpose_matrix(A);
  tm = wall_seconds() - tm;

  if(dump_files)
    write_matrix_mm(At, "er_tran.mm");
  clear_matrix(At);
  return(tm);
}

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


int main(int argc, char * argv[]){
  
  sparsemat_t *mat;
  double laptime = 0.0;
  uint32_t use_model;
  enum FLAVOR {GENERIC=1, ALL=2};

  /* process command line */
  args_t args;  
  struct argp argp = {options, parse_opt, 0, "Transpose a sparse matrix.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  
  double nz_per_row = args.gstd.nz_per_row;
  double edge_prob = args.gstd.edge_prob;
  int64_t numrows = args.gstd.numrows;
  edge_type edge_type = DIRECTED;
  self_loops loops = NOLOOPS;
  
  if(args.gstd.readfile == 0){
    resolve_edge_prob_and_nz_per_row(&edge_prob, &nz_per_row, numrows, edge_type, loops);
  }
  
  if(!args.std.quiet ) {
    fprintf(stderr,"Running C version of transpose matrix\n");
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
  
  if( args.gstd.readfile ) {
     mat = read_matrix_mm(args.gstd.filename);     
     if(!mat){printf("ERROR: transpose_matrix: read_matrix (%s) failed\n", args.gstd.filename); exit(1);}
  } else {
    mat = random_graph(numrows, args.gstd.model, edge_type, loops, edge_prob, args.std.seed);
    if(!mat){printf("ERROR: transpose_matrix: erdos_renyi_graph failed\n"); exit(1);}
  }

  if(!args.std.quiet) {
    printf("Input matrix stats:\n");
    spmat_stats(mat);
  }
  if(args.std.dump_files)
    dump_matrix(mat, 20, "dump.out");
  
  for( use_model=1; use_model < 2; use_model *=2 ) {
    switch( use_model & args.std.models_mask ) {
    case GENERIC:
    if(!args.std.quiet) printf("transpose matrix : ");
    laptime = transpose_generic(mat, args.std.dump_files);
    break;
    default:
    continue;
    }
    if(!args.std.quiet) printf("  %8.3lf seconds \n", laptime);
  }

  return(0);
}

