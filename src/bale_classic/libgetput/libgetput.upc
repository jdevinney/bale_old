/******************************************************************
//
//
//  Copyright(C) 2019, Institute for Defense Analyses
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
 ********************************************************g483*********/ 
/*! \file libgetput.upc
 * \brief some standard parallel programming support functions
 */
#include "libgetput.h"

#if __UPC__ && __UPC_ATOMIC__ && !( __cray__ || _CRAYC )
// this is relevant for BUPC or GUPC
#include <upc_atomic.h>
upc_atomicdomain_t * lgp_atomic_domain;
#endif


/*!
 * \brief Wrapper for atomic add to help with ifdef noise
 * \ingroup libgetputgrp
 */
void lgp_atomic_add(SHARED int64_t * ptr, int64_t index, int64_t value) {
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //shmem_int64_atomic_add(&ptr[lindex], value, pe);
  shmem_atomic_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
  _amo_aadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, NULL, UPC_ADD, &ptr[index], &value, NULL);
#endif
}

/*!
 * \brief Wrapper for non-blocking atomic add to help with ifdef noise
 * \ingroup libgetputgrp
 */
void lgp_atomic_add_async(SHARED int64_t * ptr, int64_t index, int64_t value){
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //shmem_int64_atomic_add(&ptr[lindex], value, pe);
  shmem_atomic_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
#pragma pgas defer_sync
  _amo_aadd_upc(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#endif
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_inc(SHARED int64_t * ptr, int64_t index) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_fetch_inc(&ptr[lindex], pe);
  ret = shmem_atomic_fetch_inc(&ptr[lindex], (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_afadd(&ptr[index], 1L);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], 1L);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_INC, &ptr[index], NULL, NULL);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_add(SHARED int64_t * ptr, int64_t index, int64_t value) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_fetch_add(&ptr[lindex], value, pe);
  ret = shmem_atomic_fetch_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_afadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_ADD, &ptr[index], &value, NULL);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic compare and swap to help with ifdef noise
 * \return the old value
 * \ingroup libgetputgrp
 */
int64_t lgp_cmp_and_swap(SHARED int64_t * ptr, int64_t index, int64_t cmp_val, int64_t swap_val) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_compare_swap(&ptr[lindex], cmp_val, swap_val, pe);
  ret = shmem_atomic_compare_swap(&ptr[lindex], cmp_val, swap_val, (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_acswap_upc(&ptr[index], cmp_val, swap_val);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_cswap_relaxed(&ptr[index], cmp_val, swap_val);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_CSWAP, &ptr[index], &cmp_val, &swap_val);
#endif
  return(ret);
}


/******************************************************************************************/
/* COLLECTIVES */
/******************************************************************************************/
#if __UPC__

/*!
 * \brief function to print a banner and do initialization if needed.
 * \param argc  from main
 * \param argv  from main
 * \return the old value
 * \ingroup libgetputgrp
 */
void lgp_init(int argc, char *argv[]) {
  
  setlocale(LC_NUMERIC,"");

#if __UPC_ATOMIC__ && !( __cray__ || _CRAYC )
  lgp_atomic_domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_ADD | UPC_INC | UPC_MAX | UPC_MIN | UPC_CSWAP, 0);
#endif
}


/*!
 * \brief function to shutdown a model if needed
 * \ingroup libgetputgrp
 */
void lgp_finalize(){
  return;
}

#define Define_Reducer( NAME, XTYPE, STYPE, RED_FUNC, UPC_FUNC)         \
  XTYPE NAME (XTYPE myval) {                                            \
    static shared STYPE *dst=NULL, * src;                            \
    if (dst==NULL) {                                                 \
      dst = upc_all_alloc(THREADS, sizeof(STYPE));                   \
      src = upc_all_alloc(THREADS, sizeof(STYPE));                      \
    }                                                                   \
    src[MYTHREAD] = myval;                                              \
    upc_barrier;                                                        \
    RED_FUNC(dst,src,UPC_FUNC, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC); \
    upc_barrier;                                                        \
    return dst[0]; }

Define_Reducer(lgp_reduce_add_l, int64_t, int64_t, upc_all_reduceL, UPC_ADD)
Define_Reducer(lgp_reduce_min_l, int64_t, int64_t, upc_all_reduceL, UPC_MIN)
Define_Reducer(lgp_reduce_max_l, int64_t, int64_t, upc_all_reduceL, UPC_MAX)

Define_Reducer(lgp_reduce_or_int, int, int, upc_all_reduceL, UPC_OR)

Define_Reducer(lgp_reduce_add_d, double, double, upc_all_reduceD, UPC_ADD)
Define_Reducer(lgp_reduce_min_d, double, double, upc_all_reduceD, UPC_MIN)
Define_Reducer(lgp_reduce_max_d, double, double, upc_all_reduceD, UPC_MAX)




/******************************************************************************************/
/******************************************************************************************/
#elif USE_SHMEM

void lgp_init(int argc, char *argv[]) {
  shmem_init();
  setlocale(LC_NUMERIC,"");
}

void lgp_finalize(){
  shmem_finalize();
}
static void *setup_shmem_reduce_workdata(long **psync, size_t xsize) {
  int *work;
  int i;
  
  work=shmem_malloc(_SHMEM_REDUCE_MIN_WRKDATA_SIZE*xsize);
  *psync=shmem_malloc(_SHMEM_REDUCE_SYNC_SIZE*sizeof(long));
  for(i=0;i<_SHMEM_REDUCE_SYNC_SIZE;++i) {
    (*psync)[i]=_SHMEM_SYNC_VALUE;
  }
  shmem_barrier_all();
  return work;
}


/* Macro to define wrapper around shmem reduce functions */
#define Define_Reducer( NAME, XTYPE, STYPE, SHMEM_FUNC)             \
  XTYPE NAME (XTYPE myval) {                                        \
    static STYPE *buff=NULL, *work; static long *sync;              \
    if (buff==NULL) {                                               \
      buff=shmem_malloc(2*sizeof(STYPE));                           \
      work=setup_shmem_reduce_workdata(&sync,sizeof(STYPE));        \
    }                                                               \
    buff[0]=myval;                                                  \
    SHMEM_FUNC(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);        \
    shmem_barrier_all();                                            \
    return buff[1]; }

/*
long lgp_reduce_add_l(long myval){
  static long *buff=NULL, *work; static long *sync;
  if (buff==NULL) {
    buff=shmem_malloc(2*sizeof(long));
    work=setup_shmem_reduce_workdata(&sync,sizeof(long));
  }
  buff[0]=myval;
  shmem_long_sum_to_all(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);
  shmem_barrier_all();
  return buff[1];
}
*/
Define_Reducer(lgp_reduce_add_l, int64_t, int64_t, shmem_longlong_sum_to_all)
Define_Reducer(lgp_reduce_min_l, int64_t, int64_t, shmem_longlong_min_to_all)
Define_Reducer(lgp_reduce_max_l, int64_t, int64_t, shmem_longlong_max_to_all)

Define_Reducer(lgp_reduce_add_d, double, double, shmem_double_sum_to_all)
Define_Reducer(lgp_reduce_min_d, double, double, shmem_double_min_to_all)
Define_Reducer(lgp_reduce_max_d, double, double, shmem_double_max_to_all)

Define_Reducer(lgp_reduce_or_int, int, int, shmem_int_or_to_all)
/*!
* \ingroup libgetputgrp
*/
void lgp_shmem_write_upc_array_int64(SHARED int64_t *addr, size_t index, size_t blocksize, int64_t val) {
  int pe;
  size_t local_index;
  int64_t *local_ptr;

  /* asupc_init tests that (long long) == (int64_t) */

  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(int64_t*)(( (char*)addr ) + local_index);

  shmem_int64_p ( local_ptr, val, pe );
}

/*!
* \ingroup libgetputgrp
*/
int64_t lgp_shmem_read_upc_array_int64(const SHARED int64_t *addr, size_t index, size_t blocksize) {
  int pe;
  size_t local_index;
  int64_t *local_ptr;

  /* asupc_init tests that (long long) == (int64_t) */

  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(int64_t*)(( (char*)addr ) + local_index);

  return shmem_int64_g ( local_ptr, pe );
}

#endif

/******************************************************************************************/
/******************************************************************************************/

/*!
  \brief Compute partial sums across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<=m} x_i\f$ 
* \ingroup libgetputgrp
*/
int64_t lgp_partial_add_l(int64_t x) {

  SHARED int64_t * tmp = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t out = 0;
  
  lgp_put_int64(tmp, MYTHREAD, x);

  lgp_barrier();

  for (int i = 0; i <= MYTHREAD; i++) {
    out += lgp_get_int64(tmp, i);
  }

  lgp_barrier();

  return out;
}

/*! 
  \brief Compute prior partial sums (not including this value) across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<m} x_i\f$ 
  * \ingroup libgetputgrp
  */
int64_t lgp_prior_add_l(int64_t x) {
  return lgp_partial_add_l(x) - x;
}


/*!
 * \brief This routine finds the min average and max of a collection 
 *  of myval's (int64_t's) across all threads
 *
 * This routine is collective and fills in a struct that holds the min average max
 * of the values given by each thread.
 *
 * \param s struct hold the computed reductions 
 * \param myval The value to contributed by MYTHREAD
 * \param dem   the denominator for the average.
 * \return to the minavgmax struct the computed reductions
 * \ingroup libgetputgrp
 *
 */
int64_t lgp_min_avg_max_l(minavgmaxL_t *s, int64_t myval, int64_t dem) {
    long retval=0;

    s->min = lgp_reduce_min_l(myval);
    s->max = lgp_reduce_max_l(myval);
    s->avg = lgp_reduce_add_l(myval) / dem;
    return( retval );
}

/*!
 * \brief This routine finds the min average and max of a collection 
 *  of myval's (int64_t's) across all threads
 *
 * This routine is collective and fills in a struct that holds the min average max
 * of the values given by each thread.
 * \param s struct hold the computed reductions 
 * \param myval The value to contributed by MYTHREAD
 * \param dem   the denominator for the average.
 * \return to the minavgmax struct the computed reductions
 * \ingroup libgetputgrp
 *
 */
int64_t lgp_min_avg_max_d(minavgmaxD_t * s, double myval, int64_t dem){
  long retval = 0;
  s->min = lgp_reduce_min_d(myval);
  s->max = lgp_reduce_max_d(myval);
  s->avg = lgp_reduce_add_d(myval)/dem;  
  return( retval );
}


#if 0
/* utility to print some basic stats about a run */
void dump_header(int argc, char *argv[]) {
  int i;
  char datestr[64];
  time_t dattmp;

  if(MYTHREAD==0) {
    dattmp = time(NULL);
    strcpy(datestr , asctime( localtime( &dattmp ) ) );
        for(i=0; i<strlen(datestr); i++){
           if(datestr[i] == 0xa) 
                   datestr[i] = (char)0;
        }
        fprintf(stderr,"%s: %s\n", argv[0], datestr); 

    fprintf(stderr,"Number of Threads %d\n",THREADS);
    for(i=0; i<argc; i++)
       fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
  }
}
#endif



/*! 
 * \brief This routine uses gettimeofday routine to give access 
 *  to the wall clock timer on most UNIX-like systems.
 * \ingroup libgetputgrp
*/
double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

