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

#include <libgetput.h>
#include <spmat.h>
#include "randperm_alternates.h"
#include <std_options.h>

/*! \file randperm.upc
 * \brief Demo program that runs the variants of randperm kernel. This program
 * generates a random permutation in parallel.
 */

/*! 
\page randperm_page Random Permutation

Demo program that runs the variants of randperm kernel. This program
generates a random permutation in parallel. The algorithm used is the 
"dart throwing algorithm" found in 
 P.B.Gibbon, Y.Matias, and V.L.Ramachandran. Efficient low-contention Parallel algorithms.
 J. of Computer and System Sciences, 53:417-442, Dec 1992.

Interestingly, we discovered what looks to be a superior and simpler algorithm. This is implemented
in alternates/randperm_agi_opt.upc.

See files spmat_agi.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
for the source for the kernels.

Run with the --help, -?, or --usage flags for run details.
 */

#include <unistd.h>

typedef struct args_t{
  int64_t l_num_rows;
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'n':
      args->l_num_rows = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"l_perm_size",'n', "NUM", 0, "Per PE length of permutation"},
    {0}
  };

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {0}
  };


int main(int argc, char * argv[]) {
  
  int64_t i;

  /* process command line */
  int ret = 0;
  args_t args;
  args.l_num_rows = 1000000;
  struct argp argp = {options, parse_opt, 0,
                      "Create a random permutation in parallel.", children_parsers};

  ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  
  if(!MYTHREAD){
    bale_app_write_int(&args.std, "num_rows_per_pe", args.l_num_rows);
    write_std_options(&args.std);
  }

  int64_t numrows = args.l_num_rows * THREADS;

  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  SHARED int64_t * out;
  int64_t seed = args.std.seed + MYTHREAD;
  int64_t use_model;
  char model_str[32];
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {

    case AGI_Model:
      sprintf(model_str, "AGI");
      out = rand_permp_agi(numrows, seed);
      break;

    case EXSTACK_Model:
      sprintf(model_str, "Exstack");
      out = rand_permp_exstack(numrows, seed, args.std.buffer_size);
      break;

    case EXSTACK2_Model:
      sprintf(model_str, "Exstack2");
      out = rand_permp_exstack2(numrows, seed, args.std.buffer_size);
      break;

    case CONVEYOR_Model:
      sprintf(model_str, "Conveyor");
      out = rand_permp_conveyor(numrows, seed);
      break;

    case ALTERNATE_Model:
      //T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      sprintf(model_str, "rand_permp_agi_opt");
      out = rand_permp_agi_opt(numrows, seed);
      break;

    case 0:
      continue;
    }
    
    t1 = wall_seconds() - t1;    
    lgp_min_avg_max_d( stat, t1, THREADS );
    bale_app_write_time(&args.std, model_str, stat->avg);
    
    
    if(!is_perm(out, numrows)){
      error++;
      T0_printf("\nERROR: rand_permp_%"PRId64" failed!\n\n", use_model & args.std.models_mask);
    }
    lgp_all_free(out);
  }
  
  if( error ) {
    T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
  }
  lgp_finalize();
  return(error);
}
