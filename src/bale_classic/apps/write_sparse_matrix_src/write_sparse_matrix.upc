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
#include <getopt.h>
#include <libgetput.h>
#include <spmat.h>
#include <std_options.h>

/*! \file write_sparse_matrix.upc
 * \brief Demo program that runs the variants of write_sparse_matrix kernel. This program
 * writes a sparsemat_t to disk in a binary format.
 */

/*! 
\page write_sparse_matrix_page Write Sparse Matrix

Demo program that runs the variants of write_sparse_matrix kernel. It first generates 
a random matrix in FLAT mode and then it writes this matrix to disk
in a directory called 'write_sparse_test'.

We define a sparse matrix dataset to be the following:
 - It lives in a directory of its own
 - It contains one ASCII file called 'metadata' which contains the number of rows, columns and nonzeros in the matrix.
 - It contains N binary 'rowcnt' files. These files contain the number of nonzeros in each row for rows 0..A->numrows
 - It contains N binary 'nonzero' files. These files contain the nonzeros in each row and are ordered by row.

This application is interesting because, as implemented, it requires
us to shuffle the rows of the matrix from cyclic to block. That is,
the nonzeros for row i is stored on PE (i % THREADS) in the
sparsemat_t data structure.  However, we wish PE 0 to write out the
first block of (approx) A->numrows/THREADS rows. One way to do this
would be to call permute_sparse_matrix to get a copy of the matrix
whose rows are distributed in the block layout. However, that would
require 2x the storage space of the matrix. We don't want the
write_sparse_matrix routine to have this requirement. That is where
things get interesting. We have a fixed buffer for writing on each PE
and each PE collects or is sent nonzero data to write in their current
buffer. In AGI (where the PEs just get the data) or synchronous
exstack, this is easy. We don't have a nice way of doing this with
exstack2 or asynchronous conveyors yet. The reason this is a challenge
for asynchronous methods is that PEs can get into a deadlock waiting
for the records they need to complete a write buffer.

See files spmat_agi.upc, spmat_exstack.upc
for the source for the kernels.

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

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };

int main(int argc, char * argv[])
{

  /* process command line */
  args_t args = {0};
  struct argp argp = {NULL, parse_opt, 0,
                      "Parallel sparse matrix transpose.", children_parsers};

  args.gstd.l_numrows = 1000000;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }

  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_fprintf(stderr, "ERROR: transpose: inmat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(inmat, "inmat");
  
  double t1;
  minavgmaxD_t stat[1];
  char * datadir = calloc(64, sizeof(char));
  char model_str[32];
  int64_t use_model;
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {
    case AGI_Model:
      sprintf(datadir,"%s","write_matrix_test_agi");
      write_sparse_matrix_agi(datadir, inmat);
      sprintf(model_str, "AGI");
      break;
    case EXSTACK_Model:
      sprintf(datadir,"%s","write_matrix_test_exstack");
      write_sparse_matrix_exstack(datadir, inmat, args.std.buffer_size);
      //read_sparse_matrix_agi(datadir);
      sprintf(model_str, "Exstack");
      break;
    case EXSTACK2_Model:
      continue;
      //sprintf(model_str, "Exstack2");
      //break;
    case CONVEYOR_Model:
      continue;
      sprintf(model_str, "Conveyor");
      //break;    
    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      break;
    case 0:
      continue;
    }
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    bale_app_write_time(&args.std, model_str, stat->avg);
  }
  
  free(datadir);
  clear_matrix(inmat);
  lgp_barrier();
  bale_app_finish(&args.std);
  return(0);
}
