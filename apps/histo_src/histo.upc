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

/*! \file histo.upc
 * \brief Demo program that computes a histogram of uint64_t's
 *  The number of histogram bins should be large enough that they need 
 *  to be spread across the whole machine
 */

#include "histo.h"

/*!
\page histogram_page Histogram
  
The histogram application computes the histogram of a large number
of \f$\tt int64\_t\tt \f$'s into a large number, \f$\tt M \tt \f$, of buckets.

Each processor has an array of updates which are random indicies 
in the range \f$\tt [0:M] \tt \f$.  The job of histogram is to count the number 
of occurances of each index across all PEs. 
We can accomplish this task by creating a distributed array, counts,  of length \f$\tt M \tt \f$ 
and initializing it to zero. Then each PE can run through its list and 
increments the appropriate entry in the count array.
As described, this technique could lead to a race condition 
where two or more PE's simultaneously update the same entry in count. 
Since incrementing an entry is commutative it is independent of the 
order in which the updates are done.
So, doing the updates atomically is sufficient.

Histogram is extremely order and latency tolerant.
We believe it to be representative of a number of PGAS loops are dominated by random puts.
As a buffered communication pattern, we think of it as one-sided pushes.

Usage:
histo [-h][-b count][-M mask][-n num][-T tabsize][-c num]
- -h prints this help message
- -b count is the number of packages in an exstack(2) buffer
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate
- -n num is the number of updates per thread
- -T tabsize is the table size or number of buckets per thread
- -c num number of cores/node is a scaling factor to adjust the update rate to handle multi-core nodes
*/

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
histo [-h][-b count][-M mask][-n num][-T tabsize][-c num]\n\
- -h prints this help message\n\
- -b count is the number of packages in an exstack(2) buffer\n\
- -M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate\n\
- -T tabsize is the table size or number of buckets per thread\n\
- -c num number of cores/node is a scaling factor to adjust the update rate to handle multi-core nodes\n\
\n");
  lgp_finalize();
  lgp_global_exit(0);
}


int main(int argc, char * argv[]) {
  lgp_init(argc, argv);
  
  int64_t buf_cnt = 1024;
  int64_t models_mask = ALL_Models; // run all the programing models
  int64_t l_num_ups  = 1000000;     // per thread number of requests (updates)
  int64_t lnum_counts = 1000;       // per thread size of the table
  int64_t cores_per_node = 0;       // Default to 0 so it won't give misleading bandwidth numbers

  int64_t i;

  int printhelp = 0;
  int opt; 
  while( (opt = getopt(argc, argv, "hb:M:n:c:T:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%"PRId64"" ,&buf_cnt);  break;
    case 'M': sscanf(optarg,"%"PRId64"" ,&models_mask);  break;
    case 'n': sscanf(optarg,"%"PRId64"" ,&l_num_ups);  break;
    case 'T': sscanf(optarg,"%"PRId64"" ,&lnum_counts);  break;
    case 'c': sscanf(optarg,"%"PRId64"" ,&cores_per_node); break;
    default:  break;
    }
  }
  if( printhelp ) usage(); 

  T0_fprintf(stderr,"Running histo on %d threads\n", THREADS);
  T0_fprintf(stderr,"buf_cnt (number of buffer pkgs)      (-b)= %"PRId64"\n", buf_cnt);
  T0_fprintf(stderr,"Number updates / thread              (-n)= %"PRId64"\n", l_num_ups);
  T0_fprintf(stderr,"Table size / thread                  (-T)= %"PRId64"\n", lnum_counts);
  T0_fprintf(stderr,"models_mask                          (-M)= %"PRId64"\n", models_mask);
  fflush(stderr);


  // Allocate and zero out the counts array
  int64_t num_counts = lnum_counts*THREADS;
  SHARED volatile int64_t * counts = lgp_all_alloc(num_counts, sizeof(int64_t)); assert(counts != NULL);  
  volatile int64_t *lcounts = lgp_local_part(int64_t, counts);
  for(i = 0; i < lnum_counts; i++)
    lcounts[i] = 0L;
  
  // index is a local array of indices into the shared counts array.
  // This is used by the _agi version. 
  // To avoid paying the UPC tax of computing index[i]/THREADS and index[i]%THREADS
  // when using the exstack and conveyor models
  // we also store a packed version that holds the pe (= index%THREADS) and lindx (=index/THREADS)
  int64_t *index   = calloc(l_num_ups, sizeof(int64_t)); assert(index != NULL);
  int64_t *pckindx = calloc(l_num_ups, sizeof(int64_t)); assert(pckindx != NULL);
  int64_t indx, lindx, pe;
  
  srand(MYTHREAD + 120348);
  for(i = 0; i < l_num_ups; i++) {
    //indx = i % num_counts;          //might want to do this for debugging
    indx = rand() % num_counts;
    index[i] = indx;                 
    lindx = indx / THREADS;
    pe  = indx % THREADS;
    pckindx[i]  =  (lindx << 16L) | (pe & 0xffff);
  }
  double volume_per_node = (8*l_num_ups*cores_per_node)*(1.0E-9);
  
  lgp_barrier();

  int64_t use_model;
  double laptime = 0.0;
  double injection_bw = 0.0;
  int64_t num_models = 0L;               // number of models that are executed

  for( use_model=1L; use_model < 32; use_model *=2 ) {

    switch( use_model & models_mask ) {
    case AGI_Model:
      T0_fprintf(stderr,"      AGI: ");
      laptime = histo_agi(index, l_num_ups, (SHARED int64_t *)counts);
      num_models++;
      break;
    
    case EXSTACK_Model:
      T0_fprintf(stderr,"  Exstack: ");
      laptime = histo_exstack(pckindx, l_num_ups, (int64_t *)lcounts, buf_cnt);
      num_models++;
      break;

    case EXSTACK2_Model:
      T0_fprintf(stderr," Exstack2: ");
      laptime = histo_exstack2(pckindx, l_num_ups, (int64_t *)lcounts, buf_cnt);
      num_models++;
      break;

    case CONVEYOR_Model:
      T0_fprintf(stderr,"Conveyors: ");
      laptime = histo_conveyor(pckindx, l_num_ups, (int64_t *)lcounts);
      num_models++;
      break;

    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;      
      //T0_fprintf(stderr,"Alternate: ");
      //laptime = histo_exstack2_cyclic(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack2_goto(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack_function(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack2_function(pckindx, l_num_ups, lcounts, buf_cnt);
      //num_models++;
      break;

    default:
      continue;

    }
    injection_bw = volume_per_node / laptime;
    T0_fprintf(stderr,"  %8.3lf seconds  %8.3lf GB/s injection bandwidth\n", laptime, injection_bw);
  }

  lgp_barrier();

  // Check the results
  // Assume that the atomic add version will correctly zero out the counts array
  for(i = 0; i < l_num_ups; i++) {
#pragma pgas defer_sync
    lgp_atomic_add((SHARED int64_t *)counts, index[i], -num_models);
  }
  lgp_barrier();

  int64_t num_errors = 0, totalerrors = 0;
  for(i = 0; i < lnum_counts; i++) {
    if(lcounts[i] != 0L) {
      num_errors++;
      if(num_errors < 5)  // print first five errors, report number of errors below
        fprintf(stderr,"ERROR: Thread %d error at %"PRId64" (= %"PRId64")\n", MYTHREAD, i, lcounts[i]);
    }
  }
  totalerrors = lgp_reduce_add_l(num_errors);
  if(totalerrors) {
     T0_fprintf(stderr,"FAILED!!!! total errors = %"PRId64"\n", totalerrors);   
  }
  
  lgp_all_free((SHARED int64_t *)counts);
  free(index);
  free(pckindx);
  lgp_finalize();

  return(totalerrors);
}

