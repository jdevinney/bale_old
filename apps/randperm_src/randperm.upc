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

#include <libgetput.h>
#include <spmat.h>

/*! \file randperm.upc
 * \brief Demo program that runs the variants of randperm kernel.
 */


int main(int argc, char * argv[])
{
  int64_t i, buf_cnt;
  int64_t task_mask = 0xF;
  int printhelp = 0;
  int64_t N, lN = 1000;
  int64_t seed = 101892+MYTHREAD;

  lgp_init();

  int opt; 
  while( (opt = getopt(argc, argv, "hN:M:s:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'N': sscanf(optarg,"%ld" ,&lN    );   break;
      //    case 'S': sscanf(optarg,"%ld" ,&buf_cnt  );   break;
    case 'M': sscanf(optarg,"%ld" ,&task_mask);  break;
    case 's': sscanf(optarg,"%ld", &seed); break;      
    default:  break;
    }
  }

  T0_fprintf(stderr,"Permutation size per thread (-N)= %ld\n", lN);
  //T0_fprintf(stderr,"buf_cnt (stack size)        (-S)= %ld\n", buf_cnt);
  T0_fprintf(stderr,"task_mask (-M)= %ld or of 1,2,4,8 for atomics,classic,exstack2,conveyor\n", task_mask);
  T0_fprintf(stderr,"seed (-s) = %ld\n", seed);

  if(printhelp) 
    lgp_global_exit(0);

  N = lN*THREADS;

  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;

  if(task_mask & 1L) {
    t1 = wall_seconds();
    SHARED int64_t * perm = rand_permp(N, seed, AGI_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("rand_permp_AGI:      %8.3lf\n", stat->avg);
    
    if(!is_perm(perm, N)){
      error++;
      T0_printf("\nERROR: rand_permp_AGI failed!\n\n");
    }
    lgp_all_free(perm);
  }
  
  if(task_mask & 2L) {
    t1 = wall_seconds();
    SHARED int64_t * perm = rand_permp(N, seed, EXSTACK_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("rand_permp_EXSTACK:  %8.3lf\n", stat->avg);
    
    if(!is_perm(perm, N)){
      error++;
      T0_printf("\nERROR: rand_permp_EXSTACK failed!\n\n");
    }
    lgp_all_free(perm);
  }
  
  if(task_mask & 4L) {
    t1 = wall_seconds();
    SHARED int64_t * perm = rand_permp(N, seed, EXSTACK2_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("rand_permp_EXSTACK2: %8.3lf\n", stat->avg);
    
    if(!is_perm(perm, N)){
      error++;
      T0_printf("\nERROR: rand_permp_EXSTACK2 failed!\n\n");
    }
    lgp_all_free(perm);
  }
  if(task_mask & 8L) {
    t1 = wall_seconds();
    SHARED int64_t * perm = rand_permp(N, seed, CONVEYOR_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("rand_permp_CONVEYOR: %8.3lf\n", stat->avg);
    
    if(!is_perm(perm, N)){
      error++;
      T0_printf("\nERROR: rand_permp_CONVEYOR failed!\n\n");
    }
    lgp_all_free(perm);
  }
  
  if( error ) {
    T0_fprintf(stderr,"YOU FAILED!!!!\n"); 
  }
  
  return(error);
}

