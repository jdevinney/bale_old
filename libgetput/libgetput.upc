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
/*! \file libgetput.upc
 * \brief some standard parallel programming support functions
 */

#include "libgetput.h"

/*!
 * \brief Wrapper for atomic add to help with ifdef noise
 * \ingroup libgetputgrp
 */
void lgp_atomic_add(SHARED int64_t * ptr, int64_t index, int64_t value)
{
  long ret;
#if USE_SHMEM
  long lindex = index/shmem_n_pes();
  long pe = index % shmem_n_pes();
  shmem_long_add(&ptr[lindex], value, pe);
  //printf("atomic_add  %ld, to %ld %ld %ld\n", MYTHREAD, pe,  lindex, value);
#elif _CRAYC 
  _amo_aadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#endif

}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_inc(SHARED int64_t * ptr, int64_t index)
{
  long ret;
#if USE_SHMEM
  long lindex = index/shmem_n_pes();
  long pe = index % shmem_n_pes();
  ret = shmem_long_finc(&ptr[lindex], pe);
#elif _CRAYC
  ret = _amo_afadd(&ptr[index], 1L);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], 1L);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_add(SHARED int64_t * ptr, int64_t index, int64_t value)
{
  long ret;
#if USE_SHMEM
  long lindex = index/shmem_n_pes();
  long pe = index % shmem_n_pes();
  ret = shmem_long_fadd(&ptr[lindex], value, pe);
  //printf("atomic_add  %ld, to %ld %ld %ld was %ld\n", MYTHREAD, pe,  lindex, value, ret);
#elif _CRAYC
  ret = _amo_afadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
 int64_t lgp_cmp_and_swap(SHARED int64_t * ptr, int64_t index, int64_t cmp_val, int64_t swap_val)
{
  long ret;
#if USE_SHMEM
  long lindex = index/shmem_n_pes();
  long pe = index % shmem_n_pes();
  ret = shmem_long_cswap(&ptr[lindex], cmp_val, swap_val, pe);
#elif _CRAYC
  ret = _amo_acswap_upc(&ptr[index], cmp_val, swap_val);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_cswap_relaxed(&ptr[index], cmp_val, swap_val);
#endif
  return(ret);
}


/******************************************************************************************/
/* COLLECTIVES */
/******************************************************************************************/
#if __UPC__

void lgp_init(){
  setlocale(LC_NUMERIC,"");
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

Define_Reducer(lgp_reduce_add_l, long, long, upc_all_reduceL, UPC_ADD)
Define_Reducer(lgp_reduce_min_l, long, long, upc_all_reduceL, UPC_MIN)
Define_Reducer(lgp_reduce_max_l, long, long, upc_all_reduceL, UPC_MAX)

Define_Reducer(lgp_reduce_or_int, int, int, upc_all_reduceL, UPC_OR)

Define_Reducer(lgp_reduce_add_d, double, double, upc_all_reduceD, UPC_ADD)
Define_Reducer(lgp_reduce_min_d, double, double, upc_all_reduceD, UPC_MIN)
Define_Reducer(lgp_reduce_max_d, double, double, upc_all_reduceD, UPC_MAX)




/******************************************************************************************/
/******************************************************************************************/
#elif USE_SHMEM

void lgp_init(){
  shmem_init();
  //THREADS = shmem_n_pes();
  //MYTHREAD = shmem_my_pe();
  setlocale(LC_NUMERIC,"");
}

static void *setup_shmem_reduce_workdata(long **psync, size_t xsize)
{
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
Define_Reducer(lgp_reduce_add_l, long, long, shmem_long_sum_to_all)
Define_Reducer(lgp_reduce_min_l, long, long, shmem_long_min_to_all)
Define_Reducer(lgp_reduce_max_l, long, long, shmem_long_max_to_all)

Define_Reducer(lgp_reduce_add_d, double, double, shmem_double_sum_to_all)
Define_Reducer(lgp_reduce_min_d, double, double, shmem_double_min_to_all)
Define_Reducer(lgp_reduce_max_d, double, double, shmem_double_max_to_all)

Define_Reducer(lgp_reduce_or_int, int, int, shmem_int_or_to_all)
/*!
* \ingroup libgetputgrp
*/
void lgp_shmem_write_upc_array_int64(SHARED int64_t *addr, size_t index, size_t blocksize, int64_t val)
{
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
int64_t lgp_shmem_read_upc_array_int64(const SHARED int64_t *addr, size_t index, size_t blocksize)
{
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
* \ingroup libgetputgrp
*/
int64_t lgp_partial_add_l(int64_t x){

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
  \brief Compute partial sums across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<m} x_i\f$ 
  * \ingroup libgetputgrp
  */
int64_t lgp_prior_add_l(int64_t x)
{
  return lgp_partial_add_l(x) - x;
}


#if 0
/*! \brief This routine gets the sum of myval (int) across all PEs
 *
 * This routine is collective and its return is single valued.
 *
 * \param myval The local value to be or'ed into the collection.
 * \return The or of myval across all PEs.
 *
 * \ingroup libgetputgrp
 */

int lgp_reduce_or_int(int myval )
{
  int ret;
  shared int * result = upc_all_alloc(THREADS, sizeof(int));
  shared int * tmp = upc_all_alloc(THREADS, sizeof(int));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceL(result, tmp, UPC_OR, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}

/*! \brief This routine gets the sum of myval across all PEs
 *
 * This routine is collective and its return is single valued.
 *
 * \param myval The local value to add to the collective sum
 * \return The sum of myval across all PEs
 *
 */
long lgp_reduce_add_l(long myval){
  long ret;
  shared long * result = upc_all_alloc(THREADS, sizeof(long));
  shared long * tmp = upc_all_alloc(THREADS, sizeof(long));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceL(result, tmp, UPC_ADD, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}
long lgp_reduce_max_l(long myval){
  long ret;
  shared long * result = upc_all_alloc(THREADS, sizeof(long));
  shared long * tmp = upc_all_alloc(THREADS, sizeof(long));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceL(result, tmp, UPC_MAX, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}
long lgp_reduce_min_l(long myval){
  long ret;
  shared long * result = upc_all_alloc(THREADS, sizeof(long));
  shared long * tmp = upc_all_alloc(THREADS, sizeof(long));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceL(result, tmp, UPC_MIN, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}

double lgp_reduce_add_d(double myval){
  double ret;
  shared double * result = upc_all_alloc(THREADS, sizeof(double));
  shared double * tmp = upc_all_alloc(THREADS, sizeof(double));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceD(result, tmp, UPC_ADD, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}
double lgp_reduce_max_d(double myval){
  double ret;
  shared double * result = upc_all_alloc(THREADS, sizeof(double));
  shared double * tmp = upc_all_alloc(THREADS, sizeof(double));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceD(result, tmp, UPC_MAX, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}
double lgp_reduce_min_d(double myval){
  double ret;
  shared double * result = upc_all_alloc(THREADS, sizeof(double));
  shared double * tmp = upc_all_alloc(THREADS, sizeof(double));  
  tmp[MYTHREAD] = myval;
  upc_barrier;
  upc_all_reduceD(result, tmp, UPC_MIN, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC);
  ret = result[0];
  upc_barrier;
  upc_all_free(result);
  upc_all_free(tmp);
  return(ret);
}
#endif
#if 0
long upc_all_rdc_add_l(long myval )
{
  int t;
  static shared long urdr64[THREADS];
  long retval=0;
  
  upc_barrier(__LINE__);
  urdr64[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  
  for (t = MYTHREAD; t < MYTHREAD + THREADS; t++)
    retval += urdr64[t % THREADS];
  upc_barrier(__LINE__);
  
  return( retval );
}
/*!
 * \brief This routine finds the min of myval (long) across all threads
 *
 * This routine is collective and its return is single valued.
 *
 * \param myval The value to compare
 * \return The minimum value of myval across all threads
 *
 */
int64_t lgp_reduce_min_l(int64_t myval )
{
    int t;
    static shared long urdr64[THREADS];
    long retval=0;

    urdr64[MYTHREAD] = myval;
    upc_barrier(__LINE__);
    
    retval = myval;
    for (t = MYTHREAD; t < MYTHREAD + THREADS; t++){
        if( retval > urdr64[t % THREADS] )
            retval = urdr64[t % THREADS];
    }
    upc_barrier(__LINE__);

    return( retval );
}

/*!
 * \brief This routine finds the max of myval (long) across all threads
 *
 * This routine is collective and its return is single valued.
 *
 * \param myval The value to compare
 * \return The maximum value of myval across all threads
 *
 */
int64_t upc_all_rdc_max_l(int64_t myval ){
  int t;
  static shared long urdr64[THREADS];
  long retval=0;
  
  urdr64[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  
  retval = myval;
  for (t = MYTHREAD; t < MYTHREAD + THREADS; t++){
    if( retval < urdr64[t % THREADS] )
      retval = urdr64[t % THREADS];
  }
  upc_barrier(__LINE__);
  
  return( retval );
}
#endif

#if __UPC__
/*! \brief Compute partial sums across threads.   (sometimes called SCAN)
    In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function is a collective on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i=0}^m x_i\f$ 
  * \ingroup libgetputgrp
*/
long upc_partial_add_l(long x)
{
  static shared long in[THREADS];
  static shared long out[THREADS];

  in[MYTHREAD] = x;

  upc_barrier(__LINE__);
  if (MYTHREAD == 0) {
    long t = 0;

    for (int i = 0; i < THREADS; i++) {
      out[i] = t += in[i];
    }
  }
  upc_barrier(__LINE__);

  return out[MYTHREAD];
}


/*! 
  \brief Compute partial sums across threads.  
  In the formulas below, \a m represents <tt>MYTHREAD</tt>.
  \note This function must be called on all threads.
  \param x input value \f$x_m\f$
  \return \f$\sum_{i<m} x_i\f$ 
  * \ingroup libgetputgrp
*/

long upc_prior_add_l(long x)
{
  return upc_partial_add_l(x) - x;
}

#endif

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
int64_t lgp_min_avg_max_l(minavgmaxL_t *s, int64_t myval, int64_t dem)
{
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
void dump_header(int argc, char *argv[])
{
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

/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */
double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

