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
#include <getopt.h>
#include <libgetput.h>
#include <spmat.h>
#include <spmat_opts.h>
#include <std_options.h>

//#include "alternates/transpose_matrix_alternates.h"

/*! \file transpose_matrix.upc
 * \brief Demo program that runs the variants of transpose_matrix kernel.
 */

/*! 
 * \page transpose_matrix_page Transpose Matrix
 *
 * Demo program that runs the variants of transpose matrix kernel. This program
 * generates a random square asymmetrix (via the Erdos-Renyi model) and then transposes
 * it in parallel.
 *
 * See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Run with the --help, -?, or --usage flags for run details.
 */

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
  lgp_init(argc, argv);
  
  int64_t i;
  int64_t check = 1;
  sparsemat_t * outmat;

  /* process command line */
  int ret = 0;
  args_t args;
  struct argp argp = {options, parse_opt, 0,
                      "Parallel sparse matrix transpose.", children_parsers};
  if(MYTHREAD == 0){
    ret = argp_parse(&argp, argc, argv, ARGP_NO_EXIT, 0, &args);
  }
  ret = distribute_cmd_line(argc, argv, &args, sizeof(args_t), ret);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(!MYTHREAD && !args.std.quiet){
    write_std_graph_options(&args.gstd);
    write_std_options(&args.std);
  }

  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_fprintf(stderr, "ERROR: transpose: inmat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(inmat, "inmat");
    
  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;  
  int64_t use_model;
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {
    case AGI_Model:
      outmat = transpose_matrix_agi(inmat);
      T0_fprintf(stderr, "AGI:     ");
      break;
    case EXSTACK_Model:
      outmat = transpose_matrix_exstack(inmat, args.std.buffer_size);
      T0_fprintf(stderr, "Exstack: ");
      break;
    case EXSTACK2_Model:
      outmat = transpose_matrix_exstack2(inmat, args.std.buffer_size);
      T0_fprintf(stderr, "Exstack2:");
      break;
    case CONVEYOR_Model:
      outmat = transpose_matrix_conveyor(inmat);
      T0_fprintf(stderr, "Conveyor:");
      break;    
    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      break;
    case 0:
      continue;
    }
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_fprintf(stderr, "%8.3lf\n", stat->avg);

    /* correctness check */
    if(0 && check){      
      sparsemat_t * outmatT = transpose_matrix(outmat);
      if(compare_matrix(outmatT, inmat)){
        T0_fprintf(stderr,"ERROR: transpose of transpose does not match!\n");
        if(args.std.dump_files){
          write_matrix_mm(outmat, "outmat");
          write_matrix_mm(outmatT, "outmatT");
        }
      }
      clear_matrix(outmatT);
    }
    clear_matrix(outmat);
  }
  
  clear_matrix(inmat);
  lgp_barrier();
  lgp_finalize();
  return(error);
}


