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

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
sssp [-h][-a 0,1][-e prob][-K str][-f filename]\n\
- -e = p: specify the Edge probability p\n\
- -h print this message\n\
- -M mask is the or of 1,2,4,8,16 for the models: \n\
- -N = n: Specify the number of rows_per_thread in the matrix (if using the random_graph generator)\n\
- -f filename : Specify a filename containing a matrix in MatrixMarket format to read as input\n\
- -b = count: Specify the number of packages in an exstack(2) stack\n\
\n\
\n");
  lgp_finalize();
  lgp_global_exit(0);
}

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
  int alg;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
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
  double t1;
  int64_t i, j;

  /* process command line */
  args_t args = {0};  // initialize args struct to all zero
  struct argp argp = {options, parse_opt, 0,
                      "Parallel Single Source Shortest Path (SSSP).", children_parsers};  
  args.gstd.l_numrows = 100;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  //override command line 
  // SSSP only applies to weighted directed graphs 
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
  
  if(args.std.dump_files) write_matrix_mm(mat, "sssp_inmat");

  lgp_barrier();
  
  T0_fprintf(stderr,"mat has %"PRId64" rows/cols and %"PRId64" nonzeros.\n", mat->numrows, mat->nnz);
  
  T0_fprintf(stderr,"Run sssp ...\n");

  d_array_t * tent = init_d_array(mat->numrows);
  d_array_t * comp_tent = NULL;

  uint64_t use_model, use_alg;
  double laptime = 0.0;


#define USE_BELLMAN (1L<<16)
#define USE_DELTA   (1L<<17)
  
  for( use_alg=(1L<<16); use_alg<(1L<<18); use_alg *=2 ){
    for( use_model=1L; use_model < 32; use_model *=2 ) {

      switch( (use_model & args.std.models_mask) | use_alg ) {
      case (AGI_Model | USE_BELLMAN):
        T0_fprintf(stderr,"    Bellman-Ford  AGI: ");
        //dump_tent("AGI          :",tent);
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_agi(tent, mat, 0); 
        //dump_tent("AGI:",tent);
        comp_tent = copy_d_array(tent);
        T0_fprintf(stderr,"Bellman AGI nothing to compare\n");
        break;

      case (EXSTACK_Model | USE_BELLMAN):
        T0_fprintf(stderr,"  Bellman-Ford Exstack: ");
        set_d_array(tent, INFINITY);
        //dump_tent("Ex Bell      :",tent);
        laptime = sssp_bellman_exstack(tent, mat, 0);
        //dump_tent("Ex Bell      :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Bellman Exstack nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Exstack compares success!\n");
        }
        break;

      case (EXSTACK_Model | USE_DELTA):
        T0_fprintf(stderr,"  Delta Exstack: ");
        set_d_array(tent, INFINITY);
        //dump_tent("Ex Delta     :",tent);
        laptime = sssp_delta_exstack(tent, mat, 0);
        //dump_tent("Ex Delta     :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Delta Exstack nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Delta Exstack compares success!\n");
        }
        break;

      case (EXSTACK2_Model | USE_BELLMAN):
        T0_fprintf(stderr,"  Bellman-Ford: Exstack2: ");
        set_d_array(tent, INFINITY);
        //dump_tent("Ex2 Bell     :",tent);
        laptime = sssp_bellman_exstack2(tent, mat, 0);
        //dump_tent("Ex2 Bell     :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Bellman Exstack2 nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Bellman Exstack2 compares success!\n");
        }
        break;

      case (EXSTACK2_Model | USE_DELTA):
        T0_fprintf(stderr,"  Delta Exstack2: ");
        set_d_array(tent, INFINITY);
        //dump_tent("Ex2 Delta     :",tent);
        laptime = sssp_delta_exstack2(tent, mat, 0);
        //dump_tent("Ex2 Delta     :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Delta Exstack2 nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Delta Exstack2 compares success!\n");
        }
        break;

      case (CONVEYOR_Model | USE_BELLMAN):
      T0_fprintf(stderr,"  Bellman-Ford Convey: ");
        set_d_array(tent, INFINITY);
        //dump_tent("C  Bell      :",tent);
        laptime = sssp_bellman_convey(tent, mat, 0);
        //dump_tent("C  Bell      :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Bellman Conveyor nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Bellman Conveyor compares success!\n");
        }
        break;

      case (CONVEYOR_Model | USE_DELTA):
      T0_fprintf(stderr,"  Delta Convey: ");
        set_d_array(tent, INFINITY);
        //dump_tent("C  Delta     :",tent);
        laptime = sssp_delta_convey(tent, mat, 0);
        //dump_tent("C  Delta     :",tent);
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          T0_fprintf(stderr,"Delta Conveyor nothing to compare\n");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            T0_fprintf(stderr, "Delta Conveyor compares success!\n");
        }
        break;
      }
      
      lgp_barrier();
      T0_fprintf(stderr,"  %8.3lf seconds.\n", laptime);
      // TODO: Check result
    }
  }
  
  lgp_barrier();

  clear_d_array(tent);
  clear_d_array(comp_tent);
  free(tent);
  free(comp_tent);
  lgp_finalize();
  return(0);
}
