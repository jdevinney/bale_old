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

#include <libgetput.h>
#include <spmat.h>
#include "randperm_alternates.h"

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

Usage:
randperm [-h][-n num][-M mask][-s seed]
- -h Print usage banner only
- -b count is the number of packages in an exstack(2) buffer
- -n=num Set the permutation entries per PE to n (default = 1000).
- -M=mask Set the models mask (1,2,4,8,16,32 for AGI,exstack,exstack2,conveyor,alternate)
- -s=seed Set a seed for the random number generation.
 */

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
randperm [-h][-n num][-M mask][-s seed]\n\
 -h Print usage banner only\n\
 -b count is the number of packages in an exstack(2) buffer\n\
 -n num Set the permutation entries per PE to n (default = 1000).\n\
 -M mask Set the models mask (1,2,4,8,16,32 for AGI, exstack, exstack2,conveyor,alternate)\n\
 -s seed Set a seed for the random number generation.\n\
\n");
  lgp_finalize();
  lgp_global_exit(0);
}

#include <unistd.h>
#include <getopt.h>

int main(int argc, char * argv[]) {
  lgp_init(argc, argv);
  
  int64_t i;
  int64_t models_mask = 0xF;
  int printhelp = 0;
  int64_t l_numrows = 1000000;
  int64_t numrows;
  int64_t buf_cnt = 1024;
  int64_t seed = 101892+MYTHREAD;
  int64_t cores_per_node = 1;

  int opt; 
  while( (opt = getopt(argc, argv, "b:c:hn:M:s:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%"PRId64"",&buf_cnt);  break;
    case 'c': sscanf(optarg,"%"PRId64"" ,&cores_per_node); break;
    case 'n': sscanf(optarg,"%"PRId64"",&l_numrows);   break;
    case 'M': sscanf(optarg,"%"PRId64"",&models_mask);  break;
    case 's': sscanf(optarg,"%"PRId64"", &seed); break;      
    default:  break;
    }
  }
  if(printhelp) usage();
  T0_fprintf(stderr,"Running randperm on %d threads\n", THREADS);
  T0_fprintf(stderr,"This is a demo program that runs various implementations of the randperm kernel.\n");
  T0_fprintf(stderr,"Usage:\n");
  T0_fprintf(stderr,"Permutation size per thread (-n) = %"PRId64"\n", l_numrows);
  T0_fprintf(stderr,"models_mask (-M)                 = %"PRId64" or of 1,2,4,8 for atomics,classic,exstack2,conveyor\n", models_mask);
  T0_fprintf(stderr,"buf_cnt (-b)                     = %"PRId64"\n", buf_cnt);
  T0_fprintf(stderr,"seed (-s)                        = %"PRId64"\n", seed);


  numrows = l_numrows * THREADS;

  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  SHARED int64_t * out;

  int64_t use_model;
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & models_mask ) {

    case AGI_Model:
      out = rand_permp_agi(numrows, seed);
      T0_fprintf(stderr,"rand_permp_AGI:           ");
      break;

    case EXSTACK_Model:
      out = rand_permp_exstack(numrows, seed, buf_cnt);
      T0_fprintf(stderr,"rand_permp_EXSTACK:       ");
      break;

    case EXSTACK2_Model:
      out = rand_permp_exstack2(numrows, seed, buf_cnt);
      T0_fprintf(stderr,"rand_permp_EXSTACK2:      ");
      break;

    case CONVEYOR_Model:
      out = rand_permp_conveyor(numrows, seed);
      T0_fprintf(stderr,"rand_permp_CONVEYOR:      ");
      break;

    case ALTERNATE_Model:
      //T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      out = rand_permp_agi_opt(numrows, seed);
      T0_fprintf(stderr,"rand_permp_agi_opt :      ");
      break;

    case 0:
      continue;
    }
    
    t1 = wall_seconds() - t1;    
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_fprintf(stderr,"%8.3lf\n", stat->avg);
    if(!is_perm(out, numrows)){
      error++;
      T0_printf("\nERROR: rand_permp_%"PRId64" failed!\n\n", use_model & models_mask);
    }
    lgp_all_free(out);
  }
  
  if( error ) {
    T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
  }
  lgp_finalize();
  return(error);
}
