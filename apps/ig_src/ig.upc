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

/*! \file ig.upc
 * \brief Demo program that runs index gather with the various models
 */
#include "ig.h"

#include <exstack.h>
#include <locale.h>


/*! \brief check that the local parts of two arrays agree
 * \param t first array
 * \param s second array
 * \param len length
 * \returns the number of disagreements.
 */

int64_t ig_check(int64_t *t, int64_t *s, int64_t len)
{
  int64_t errors=0;
  int64_t i;
  for(i=0; i<len; i++){
     errors += (t[i] != s[i]);
  }
  return(errors);
}

int main(int argc, char * argv[])
{
  int64_t i;
  SHARED int64_t * table;
  int64_t * ltable, *tgt, *tgt2;
  int64_t ltab_siz = 10000;
  int64_t tab_siz = ltab_siz*THREADS;
  int64_t lrequests  = 100000;
  int64_t *index, *pckindx, indx, lindx, pe;
  int64_t buf_cnt = 1024;
  int64_t num_errors = 0L, total_errors = 0L;
  int64_t task_mask = 0xF;
  int printhelp = 0;

  lgp_init();

  int opt; 
  while( (opt = getopt(argc, argv, "hT:U:S:M:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'T': sscanf(optarg,"%ld" ,&ltab_siz    );   break;
    case 'U': sscanf(optarg,"%ld" ,&lrequests   );   break;
    case 'S': sscanf(optarg,"%ld" ,&buf_cnt  );   break;
    case 'M': sscanf(optarg,"%ld" ,&task_mask);  break;
    default:  break;
    }
  }

  T0_fprintf(stderr,"Running ig on %ld threads\n", THREADS);
  T0_fprintf(stderr,"Table size per thread       (-T)= %ld\n", ltab_siz);
  T0_fprintf(stderr,"Updates (number of requests (-U)= %ld\n", lrequests );
  T0_fprintf(stderr,"buf_cnt (stack size)        (-S)= %ld\n", buf_cnt);
  T0_fprintf(stderr,"task_mask (-M)= %ld or of 1,2,4,8 for gets,exstack,exstack2,conveyor\n", task_mask);

  if(printhelp) 
    lgp_global_exit(0);


  tab_siz = ltab_siz*THREADS;
  table   = lgp_all_alloc(tab_siz, sizeof(int64_t));
  ltable  = lgp_local_part(int64_t, table);
  tgt     = calloc(lrequests, sizeof(int64_t));
  tgt2    = calloc(lrequests, sizeof(int64_t));
  index   = calloc(lrequests, sizeof(int64_t));
  pckindx = calloc(lrequests, sizeof(int64_t));

  //populate table array and the index array
  srand48( MYTHREAD*MYTHREAD*MYTHREAD*10000 + 5 ); // I don't know why
  for(i=0; i<ltab_siz; i++)
    ltable[i] = lrand48() % tab_siz;
  
  for(i = 0; i < lrequests; i++){
    indx = lrand48() % tab_siz; 
    index[i] = indx;
    lindx = indx / THREADS;      // the distributed version of indx
    pe  = indx % THREADS;      
    pckindx[i] = (lindx << 16) | (pe & 0xffff); // same thing stored as (local index, thread) "shmem style"
  }
  double t_gets, t_exstack, t_exstack2, t_conveyor, t_other;

  if(task_mask & 1L) {
    t_gets = ig_gets(tgt, index, lrequests, table);
    //T0_printf("Single word gets: %8.3lf    %5.3lf GB/th/s\n", t_gets, (2*lrequests*8)/(t_gets*(1L<<30)));
    T0_printf("Single word gets: %8.3lf seconds\n", t_gets);
  }

  if(task_mask & 2L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_exstack = ig_exstack(tgt2, pckindx, lrequests,  ltable,  buf_cnt);
    //T0_printf("Exstack Classic : %8.3lf    %5.3lf GB/th/s\n", t_exstack, (2*lrequests*8)/(t_exstack*(1L<<30)));
    T0_printf("Exstack Classic : %8.3lf    %5.3lf seconds\n", t_exstack);
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }

  if(task_mask & 4L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_exstack2 = ig_exstack2(tgt2, pckindx, lrequests,  ltable,  buf_cnt);
    T0_printf("Exstack2        : %8.3lf seconds\n", t_exstack2);
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }

  if(task_mask & 8L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_conveyor = ig_conveyor(tgt2, pckindx, lrequests,  ltable);
    //T0_printf("Conveyors       : %8.3lf    %5.3lf GB/th/s\n", t_conveyor, (2*lrequests*8)/(t_conveyor*(1L<<30)));
    T0_printf("Conveyors       : %8.3lf seconds\n", t_conveyor);
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }
  
#if 0
  if(task_mask & 16L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_other = ig_exstack2_cyclic(tgt2, pckindx, lrequests,  ltable, buf_cnt);
    T0_printf("Cyclic          : %8.3lf    %5.3lf GB/th/s\n", t_other, (2*lrequests*8)/(t_other*(1L<<30)));
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }

  if(task_mask & 32L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_other = ig_exstack2_goto(tgt2, pckindx, lrequests,  ltable, buf_cnt);
    T0_printf("GOTO            : %8.3lf    %5.3lf GB/th/s\n", t_other, (2*lrequests*8)/(t_other*(1L<<30)));
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }

  if(task_mask & 64L) {
    for(i=0; i<lrequests; i++) tgt2[i] = 0L;
    t_other = ig_exstack_function(tgt2, pckindx, lrequests,  ltable, buf_cnt);
    T0_printf("Function        : %8.3lf    %5.3lf GB/th/s\n", t_other, (2*lrequests*8)/(t_other*(1L<<30)));
    if( task_mask & 1L)
      num_errors += ig_check( tgt, tgt2, lrequests);
  }
#endif
  
  total_errors = lgp_reduce_add_l(num_errors);
  lgp_all_free(table);

  if( (task_mask & 1L) ){
    if( total_errors ) {
      T0_fprintf(stderr,"YOU FAILED!!!!\n"); 
    } 
  } else
    T0_fprintf(stderr,"Not checking\n");

  return(total_errors);
}

