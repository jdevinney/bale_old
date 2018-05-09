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

/*! \file permute_matrix.upc
 * \brief Demo program that runs the variants of permute_matrix kernel.
 */


int main(int argc, char * argv[])
{
  int64_t i, buf_cnt;
  int64_t task_mask = 0xF;
  int printhelp = 0;
  int64_t N, lN = 1000;
  int64_t seed = 101892+MYTHREAD;
  sparsemat_t * inmat, * outmat;

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

  T0_fprintf(stderr,"rows per thread (-N)= %ld\n", lN);
  //T0_fprintf(stderr,"buf_cnt (stack size)        (-S)= %ld\n", buf_cnt);
  T0_fprintf(stderr,"task_mask (-M)= %ld or of 1,2,4,8,16,32 for gets,classic,exstack2,conveyor,ex2cyclic,ex2goto\n", task_mask);
  T0_fprintf(stderr,"seed (-s) = %ld\n", seed);

  if(printhelp) 
    lgp_global_exit(0);

  N = lN*THREADS;

  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  
  SHARED int64_t * rp = rand_permp(N, seed, EXSTACK_Model);
  SHARED int64_t * cp = rand_permp(N, seed, EXSTACK_Model);  
  
  inmat = gen_uniform_sparse(N, 30);
  if(inmat == NULL){
    T0_printf("ERROR: inmat is null!\n");
    return(-1);
  }

  if(task_mask & 1L) {
    t1 = wall_seconds();
    outmat = permute_matrix(inmat, rp, cp, AGI_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("permute_matrix_AGI:          %8.3lf\n", stat->avg);    
    clear_matrix(outmat);
  }
  
  if(task_mask & 2L) {
    t1 = wall_seconds();
    outmat = permute_matrix(inmat, rp, cp, EXSTACK_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("permute_matrix_EXSTACK:      %8.3lf\n", stat->avg);    
    clear_matrix(outmat);
  }
  
  if(task_mask & 4L) {
    t1 = wall_seconds();
    outmat = permute_matrix(inmat, rp, cp, EXSTACK2_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("permute_matrix_EXSTACK2:     %8.3lf\n", stat->avg);    
    clear_matrix(outmat);
  }
  if(task_mask & 8L) {
    t1 = wall_seconds();
    outmat = permute_matrix(inmat, rp, cp, CONVEYOR_Model);
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_printf("permute_matrix_CONVEYOR:     %8.3lf\n", stat->avg);    
    clear_matrix(outmat);
  }
  
  clear_matrix(inmat);
  lgp_all_free(rp);
  lgp_all_free(cp);

  return(error);
}


