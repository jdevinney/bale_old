/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
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
/*! \file sssp.upc
 * \brief Demo application that implements Single Source Shortest Path algorithms
 */

#include "sssp.h"
#include <std_options.h>

/*!
  \page sssp_page Single Source Shortest Path


We are given a weight adjacency matrix for a (directed) graph
and an initial starting vertex. 
The problem is to find the distance 
(the sum of the weights of the edges on the least heavy path)
to all other vertices in the graph.

 */

void dump_tent(char *str, d_array_t *tent)
{
  int64_t i;
  if( MYTHREAD == 0 ){
    printf("%s ", str);
    for(i=0; i<tent->num; i++){
      printf(" %lg ", lgp_get_double(tent->entry, i) );
    }
    printf("\n");
  }
}

double sssp_answer_diff(d_array_t *A, d_array_t *B)
{
  int64_t i;
  double ldiff = 0.0;

  if(A->num != B->num)
    return(INFINITY);

  for(i=0; i<A->lnum; i++) {
    if( A->lentry[i] == INFINITY && B->lentry[i] == INFINITY )
      continue;
    ldiff += (A->lentry[i] - B->lentry[i]) * (A->lentry[i] - B->lentry[i]);
  }
  return(sqrt(lgp_reduce_add_d(ldiff)));
}



typedef struct args_t{
  double deltaStep; 
  int64_t V0;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'S': args->deltaStep = atof(arg); break;     
    case 'V': args->V0 = atoi(arg); break;
    case ARGP_KEY_INIT:
      args->deltaStep = 0.0;
      args->V0 = 0;
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
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
  double t1;
  int64_t i, j;

  /* process command line */
  args_t args = {0};  // initialize args struct to all zero
  struct argp argp = {options, parse_opt, 0,
                      "Parallel Single Source Shortest Path (SSSP).", children_parsers};  
  args.gstd.l_numrows = 50;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  //override command line 
  // SSSP only applies to weighted directed graphs 
  //args.gstd.l_numrows = 100;
  if(args.gstd.loops == 1 || args.gstd.directed == 0 || args.gstd.weighted == 0){
    T0_fprintf(stderr,"WARNING: assume the input graph is directed, weighted [0,1) graph with no loops.\n");
    T0_fprintf(stderr,"Overwriting the options.\n");
    args.gstd.loops = 0;
    args.gstd.directed = 1;
    args.gstd.weighted = 1;
  }


  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }
  
  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){T0_fprintf(stderr, "ERROR: sssp: mat is NULL!\n");return(-1);}

  if(args.V0 < 0 || args.V0 >= mat->numrows){
    T0_fprintf(stderr,"Setting V0 to 0\n");
    args.V0 = 0;
  }
  
  if(args.std.dump_files) write_matrix_mm(mat, "sssp_inmat");

  lgp_barrier();
  
  d_array_t * tent = init_d_array(mat->numrows);
  d_array_t * comp_tent = NULL;

  uint64_t use_model, use_alg;
  double laptime = 0.0;
  char model_str[32];

  T0_printf("delta step = %lf\n", args.deltaStep);

#define USE_BELLMAN (1L<<16)
#define USE_DELTA   (1L<<17)
  
  for( use_alg=(1L<<16); use_alg<(1L<<18); use_alg *=2 ){
    for( use_model=1L; use_model < 32; use_model *=2 ) {
      model_str[0] = '\0';
      switch( (use_model & args.std.models_mask) | use_alg ) {
      case (AGI_Model | USE_BELLMAN):
        sprintf(model_str, "Bellman-Ford AGP");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_agi(tent, mat, 0); 
        break;

      case (EXSTACK_Model | USE_BELLMAN):
        sprintf(model_str, "Bellman-Ford Exstack");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_exstack(tent, mat, args.V0);
        break;

      case (EXSTACK_Model | USE_DELTA):
        sprintf(model_str, "Delta Exstack");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_exstack(tent, mat, args.V0, args.deltaStep);
        break;

      case (EXSTACK2_Model | USE_BELLMAN):
        sprintf(model_str, "Bellman-Ford Exstack2");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_exstack2(tent, mat, 0);
        break;

      case (EXSTACK2_Model | USE_DELTA):
        sprintf(model_str, "Delta Exstack");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_exstack2(tent, mat, args.V0, args.deltaStep);
        break;

      case (CONVEYOR_Model | USE_BELLMAN):
        sprintf(model_str, "Bellman-Ford Conveyor");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_convey(tent, mat, args.V0);
        break;

      case (CONVEYOR_Model | USE_DELTA):
        sprintf(model_str, "Delta Conveyor");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_convey(tent, mat, args.V0, args.deltaStep);
        break;
      }
      if(model_str[0]) {
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          sprintf(model_str, "%s ()", model_str);
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            sprintf(model_str, "%s (compares)", model_str);
        }
        lgp_barrier();
        bale_app_write_time(&args.std, model_str, laptime);
      }
    }
  }
  
  lgp_barrier();

  clear_d_array(tent); free(tent);
  clear_d_array(comp_tent); free(comp_tent);

  bale_app_finish(&args.std);

  return(0);
}
