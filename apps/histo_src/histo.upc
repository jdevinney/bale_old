/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
 *****************************************************************/ 
/*! \file histo.upc
 * \brief Demo program that computes a histogram of uint64_t's
 *  The number of histogram bins should be large enough that they need 
 *  to be spread across the whole machine
 */

#include "histo.h"

int main(int argc, char * argv[])
{

  int ret;
  int64_t i, j, fromth;
  int64_t pe, row, lindx, *idx, * index, *pckindx, * lcounts;
  int64_t t1, t2, t3;
  int64_t mt1, mt2, mt3;
  int64_t lrequests  = 2000000;            // per thread number of requests (updates)
  int64_t lnum_counts = 1000;              // per thread size of the table
  int64_t correctcount = 0L;
  int64_t buf_cnt = 1024;
  int64_t task_mask = 0xF;
  double NV = 0.0; //Node Volume
  
  lgp_init();
  //setlocale(LC_NUMERIC,"");

  int printhelp = 0;
  int opt; 
  while( (opt = getopt(argc, argv, "hT:U:S:M:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'T': sscanf(optarg,"%ld" ,&lnum_counts);  break;
    case 'U': sscanf(optarg,"%ld" ,&lrequests);  break;
    case 'S': sscanf(optarg,"%ld" ,&buf_cnt);  break;
    case 'M': sscanf(optarg,"%ld" ,&task_mask);  break;
    default:  break;
    }
  }
  T0_fprintf(stderr,"Running histo on %ld threads\n", THREADS);
  T0_fprintf(stderr,"Table size per thread (-T)= %ld\n", lnum_counts);
  T0_fprintf(stderr,"Updates per thread    (-U)= %ld\n", lrequests);
  T0_fprintf(stderr,"buf_cnt (stack size)  (-S)= %ld\n", buf_cnt);
  T0_fprintf(stderr,"task_mask (-M)= %ld or of 1,2,4,8 for atomic,exstack,exstack2,conveyor\n", task_mask);
  fflush(stderr);

  if( printhelp ) 
    return(0);
  
  int64_t num_counts = lnum_counts*THREADS;
  double t_classic, t_exstack2, t_conveyor, t_atomic, t_other;

  // index is a local array of indices into the shared counts array.
  // to avoid dealing with the UPC tax of computing i/THREADS and i%THREADS
  // we will store the index array and the packed version that
  // holds the pe (= index%THREADS) and lindx (=index/THREADS)
  index = calloc(lrequests, sizeof(int64_t)); assert(index != NULL);
  pckindx  = calloc(lrequests, sizeof(int64_t)); assert(pckindx != NULL);
  SHARED int64_t * counts = lgp_all_alloc(num_counts, sizeof(int64_t)); assert(counts != NULL);  
  lcounts = lgp_local_part(int64_t, counts);
  for(i = 0; i < lnum_counts; i++)
    lcounts[i] = 0L;
  
  srand48( MYTHREAD*MYTHREAD*10000 + 5 );
  int64_t indx;
  for(i = 0; i < lrequests; i++) {
    //indx   = lrand48() % num_counts; // a global or shared index into the counts array
    indx = i % num_counts;
    index[i] = indx;
    lindx = indx / THREADS;            // the distributed version of indx
    pe  = indx % THREADS;      
    pckindx[i]  =  (lindx << 16) | (pe & 0xffff);
  }
  
  lgp_barrier();

  if(task_mask & 2L){
    t_classic = histo_exstack(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    //T0_printf("Exstack Classic : %8.3lf    %5.3lf GB/th/s\n", t_classic, (lrequests*8)/(t_classic*(1L<<30)));
    T0_printf("Exstack Classic : %8.3lf seconds\n", t_classic);
    fflush(NULL);
  }

  if(task_mask & 4L){
    t_exstack2 = histo_exstack2(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    //T0_printf("Exstack2        : %8.3lf    %5.3lf GB/th/s\n", t_exstack2, (lrequests*8)/(t_exstack2*(1L<<30)));
    T0_printf("Exstack2        : %8.3lf seconds\n", t_exstack2);
    fflush(NULL);
  }

  if(task_mask & 8L){
    t_conveyor = histo_conveyor(pckindx, lrequests, lcounts);
    correctcount++;
    //T0_printf("Conveyors       : %8.3lf    %5.3lf GB/th/s\n", t_conveyor, (lrequests*8)/(t_conveyor*(1L<<30)));
    T0_printf("Conveyors       : %8.3lf seconds\n", t_conveyor);
    fflush(NULL);
  }

#if 0
  if(task_mask & 16L){
    t_other = histo_exstack2_cyclic(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    T0_printf("Exstack2_cyclic : %8.3lf    %5.3lf GB/th/s\n", t_other, (lrequests*8)/(t_other*(1L<<30)));
    fflush(NULL);
  }

  if(task_mask & 32L){
    t_other = histo_exstack2_goto(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    T0_printf("Exstack2_goto   : %8.3lf    %5.3lf GB/th/s\n", t_other, (lrequests*8)/(t_other*(1L<<30)));
    fflush(NULL);
  }

  if(task_mask & 64L){
    t_other = histo_exstack_function(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    T0_printf("Exstack_function: %8.3lf    %5.3lf GB/th/s\n", t_other, (lrequests*8)/(t_other*(1L<<30)));
    fflush(NULL);
  }

  if(task_mask & 128L){
    t_other = histo_exstack2_function(pckindx, lrequests, lcounts, buf_cnt);
    correctcount++;
    T0_printf("Exstack2_func   : %8.3lf    %5.3lf GB/th/s\n", t_other, (lrequests*8)/(t_other*(1L<<30)));
    fflush(NULL);
  }
#endif

  if(task_mask & 1L){
    t_atomic = histo_atomic(index, lrequests,  counts, correctcount);
    //T0_printf("UPC Atomics     : %8.3lf    %5.3lf GB/th/s\n", t_atomic, (lrequests*8)/(t_atomic*(1L<<30)));
    T0_printf("UPC Atomics     : %8.3lf seconds\n", t_atomic);
    fflush(NULL);
  }
  //fflush(stderr);
  
  lgp_barrier();

  int64_t error = 0, totalerrors = 0;
  if(task_mask & 1L){  // Assume the atomic add version will clean everything up
    for(i = 0; i < lnum_counts; i++) {
      if(lcounts[i] != 0L){
        error++;
        if(error < 5)  // print first five errors, report all the errors below
          fprintf(stderr,"ERROR: Thread %d error at %ld (= %ld)\n", MYTHREAD, i, lcounts[i]);
        break;
      }
    }
    totalerrors = lgp_reduce_add_l(error);
    if(totalerrors)
      T0_fprintf(stderr,"FAILED!!!! total errors = %ld\n", totalerrors);   
  }else{
    T0_fprintf(stderr,"No check...\n");
  }

  lgp_all_free(counts);
  return(totalerrors);
}

