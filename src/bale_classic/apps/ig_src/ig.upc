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

/*! \file ig.upc
 * \brief Demo program that runs index gather with the various models
 */

#include "ig.h"
#include <std_options.h>

/*!
\page indexgather_page Indexgather

The indexgather application is random gather across the whole memory 
based on local arrays of indices.
Each processor has a local array of read indices 
into a shared array,\f$\tt table\tt\f$, of size \f$\tt T\tt\f$. 
The PE's store the random gets from table into the corresponding
entry in a local target array.

Since we don't overwrite the target array, this is completely order and latency tolerant.

As a PGAS implementation this is a collection of random gets.
As a buffered communication pattern it is both a request and response.

The simplest form of indexgather uses \f$\tt int64\_t\tt\f$'s for all the arrays.
This allows one to reuse communication buffers.
It might be more interesting to consider the case where the table entries were multiple words.

Run with the --help, -?, or --usage flags for run details.
*/


/*! \brief check the target array the local parts of two arrays agree
 * \param use_model the model being used
 * \param *tgt pointer to the target array
 * \param *index array of indices
 * \param l_num_req local number of requests (length of index array)
 * \returns the number of errors.
 */
int64_t ig_check_and_zero(int64_t use_model, int64_t *tgt, int64_t *index, int64_t l_num_req) {
  int64_t errors=0;
  int64_t i;
  lgp_barrier();
  for(i=0; i<l_num_req; i++){
    if(tgt[i] != (-1)*(index[i] + 1)){
      errors++;
      if(errors < 5)  // print first five errors, report all the errors
        fprintf(stderr,"ERROR: model %"PRId64": Thread %d: tgt[%"PRId64"] = %"PRId64" != %"PRId64")\n",
                use_model,  MYTHREAD, i, tgt[i], (-1)*(index[i] + 1));
               //use_model,  MYTHREAD, i, tgt[i],(-1)*(i*THREADS+MYTHREAD + 1) );
    }
    tgt[i] = 0;
  }
  if( errors > 0 )
    fprintf(stderr,"ERROR: %"PRId64": %"PRId64" total errors on thread %d\n", use_model, errors, MYTHREAD);
  lgp_barrier();
  return(errors);
}


typedef struct args_t{
  int64_t l_num_req;
  int64_t l_tbl_size;
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'n':
      args->l_num_req = atol(arg); break;
    case 'T':
      args->l_tbl_size = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"num_requests",'n', "NUM", 0, "Number of reads per PE from the table"},
    {"table_size", 'T', "SIZE", 0, "Number of entries per PE in the table"},
    {0}
  };

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {0}
  };


int main(int argc, char * argv[]) {

  lgp_init(argc, argv);
  
  int64_t i;
  int64_t num_errors = 0L, total_errors = 0L;
  int64_t printhelp = 0;

  /* process command line */
  int ret = 0;
  args_t args;
  if(MYTHREAD == 0){
    args.l_tbl_size = 1000;
    args.l_num_req = 100000;
    struct argp argp = {options, parse_opt, 0,
                        "Many remote reads from a distributed table.", children_parsers};
    ret = argp_parse(&argp, argc, argv, ARGP_NO_EXIT, 0, &args);
  }
  
  ret = distribute_cmd_line(argc, argv, &args, sizeof(args_t), ret);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(!MYTHREAD && !args.std.quiet){
    T0_fprintf(stderr,"Number of Request / PE   (-n): %"PRId64"\n", args.l_num_req );
    T0_fprintf(stderr,"Table size / PE          (-T): %"PRId64"\n\n", args.l_tbl_size);
    write_std_options(&args.std);
  }
  
  int64_t bytes_read_per_request_per_node = 8*2*args.std.cores_per_node;
  
  // Allocate and populate the shared table array 
  int64_t tab_siz = args.l_tbl_size*THREADS;
  SHARED int64_t * table   = lgp_all_alloc(tab_siz, sizeof(int64_t)); assert(table != NULL);
  int64_t *ltable  = lgp_local_part(int64_t, table);
  // fill the table with the negative of its shared index
  // so that checking is easy
  for(i=0; i<args.l_tbl_size; i++)
    ltable[i] = (-1)*(i*THREADS + MYTHREAD + 1);
  
  // As in the histo example, index is used by the _agi version.
  // pckindx is used my the buffered versions
  int64_t *index   = calloc(args.l_num_req, sizeof(int64_t)); assert(index != NULL);
  int64_t *pckindx = calloc(args.l_num_req, sizeof(int64_t)); assert(pckindx != NULL);
  int64_t indx, lindx, pe;
  srand(MYTHREAD+ args.std.seed);
  for(i = 0; i < args.l_num_req; i++){
    indx = rand() % tab_siz;
    index[i] = indx;
    lindx = indx / THREADS;      // the distributed version of indx
    pe  = indx % THREADS;      
    pckindx[i] = (lindx << 16) | (pe & 0xffff); // same thing stored as (local index, thread) "shmem style"
  }

  int64_t *tgt  = calloc(args.l_num_req, sizeof(int64_t)); assert(tgt != NULL);

  lgp_barrier();

  int64_t use_model;
  double laptime = 0.0;
  double volume_per_node = (2*8*args.l_num_req*args.std.cores_per_node)*(1.0E-9);
  double injection_bw = 0.0;

  for( use_model=1L; use_model < 32; use_model *=2 ){
     switch( use_model & args.std.models_mask ){
     case AGI_Model:
        T0_fprintf(stderr,"      AGI: ");
        laptime = ig_agi(tgt, index, args.l_num_req, table);
     break;
     
     case EXSTACK_Model:
        T0_fprintf(stderr,"  Exstack: ");
        laptime =ig_exstack(tgt, pckindx, args.l_num_req,  ltable,  args.std.buffer_size); 
     break;
   
     case EXSTACK2_Model:
        T0_fprintf(stderr," Exstack2: ");
        laptime = ig_exstack2(tgt, pckindx, args.l_num_req,  ltable, args.std.buffer_size);
     break;
   
     case CONVEYOR_Model:
        T0_fprintf(stderr,"Conveyors: ");
        laptime = ig_conveyor(tgt, pckindx, args.l_num_req,  ltable);
     break;
   
     case ALTERNATE_Model:
       T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
       T0_fprintf(stderr,"Alternate: ");
       //laptime = ig_exstack2_cyclic(tgt, pckindx, l_num_req,  ltable, buf_cnt);
       //laptime = ig_exstack2_goto(tgt, pckindx, l_num_req,  ltable, buf_cnt);
       //laptime = ig_exstack_function(tgt, pckindx, l_num_req,  ltable, buf_cnt);
       //laptime = ig_exstack_pkg(tgt, pckindx, l_num_req,  ltable, buf_cnt);
       break;
   
     default:
       continue;
       
     }
     injection_bw = volume_per_node / laptime;
     T0_fprintf(stderr,"  %8.3lf seconds  %8.3lf GB/s/Node\n", laptime, injection_bw);
   
     num_errors += ig_check_and_zero(use_model, tgt, index, args.l_num_req);
  }

  total_errors = lgp_reduce_add_l(num_errors);
  if( total_errors ) {
    T0_fprintf(stderr,"YOU FAILED!!!!\n"); 
  } 

  lgp_barrier();
  lgp_all_free(table);
  free(index);
  free(pckindx);
  free(tgt);
  lgp_finalize();
  return(0);
}

