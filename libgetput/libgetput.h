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

/*! \file libgetput.h
 * \brief The header file for libgetput library.
 * Mostly reduction operations, wrappers for atomics and basic stats collection.
 *
 * * \ingroup libgetputgrp
*/



#ifndef libgetput_INCLUDED

#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <locale.h>
#include <sys/time.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if _CRAYC
#include <intrinsics.h>
#else
#define _rtc() 0  /*!< the realtime clock */
#endif

# if !HAVE_SHMEM_FREE
// Fall back to old-style names
# define shmem_free shfree
# define shmem_malloc shmalloc
# define shmem_align shmemalign
# define shmem_global_exit exit
# endif

/* defines to support the different models of global and buffered references */
#define AGI_Model 1L      /*!< the Atomic or Generic Interface (straight UPC/SHMEM) */
#define EXSTACK_Model 2L  /*!< the exstack bulk synchronous buffering model */
#define EXSTACK2_Model 4L /*!< the exstack2 asynchronous buffering model */
#define CONVEYOR_Model 8L /*!< the conveyor buffering model */

/*! \struct minavgmaxL_t 
 * \brief A structure to return the global stats computed by lgp_min_avg_max_l
 */
typedef struct minavgmaxL_t{  
  int64_t min;    /*!< int64_t to hold the min. */ 
  int64_t avg;    /*!< int64_t to hold the average. */ 
  int64_t max;    /*!< int64_t to hold the max. */ 
}minavgmaxL_t;

/*! \struct minavgmaxD_t 
 * \brief A structure to return the global stats computed by lgp_min_avg_max_d
 */
typedef struct minavgmaxD_t{
  double min;    /*!< the min. */ 
  double avg;    /*!< the average. */ 
  double max;    /*!< the max. */ 
}minavgmaxD_t;


#if __UPC__

/////////////////////////////////////////////////////////////////
///////                   UPC SECTION                   /////////
/////////////////////////////////////////////////////////////////

#include <upc.h>
#include <upc_collective.h>

#define SHARED shared  /*!< a wrapper so that it can be empty for SHMEM */
/*! \ingroup libgetputgrp */
#define lgp_barrier() upc_barrier /*!< wrapper for global barrier */
/*! \ingroup libgetputgrp */
#define lgp_local_part(t,p) (( t * )((p)+MYTHREAD)) /*!< wrapper for localizing a global pointer */
/*! \ingroup libgetputgrp */
#define lgp_fence() upc_fence /*!< wrapper for memory fence */
/*! \ingroup libgetputgrp */
#define lgp_all_alloc(nblock_on_all_threads,blocksize) (upc_all_alloc((nblock_on_all_threads),(blocksize))) /*!< wrapper for global allocation */
/*! \ingroup libgetputgrp */
#define lgp_all_free(p) (upc_all_free(p)) /*!< wrapper for a free of global memory */
/*! \ingroup libgetputgrp */
#define lgp_memput(dst,src,n,index) (upc_memput( (dst) + (index), (src), (n) ))  /*!< wrapper for upc_memput */
/*! \ingroup libgetputgrp */
#define lgp_put_int64(array,index,val) (((shared int64_t*)(array))[index]=(val)) /*!< wraps a shared write to be compatible with shmem */
/*! \ingroup libgetputgrp */
#define lgp_get_int64(array,index)(((const shared int64_t*)(array))[index]) /*!< wraps a shared read to be compatible with shmem */
/*! \ingroup libgetputgrp */
#define lgp_global_exit(i) upc_global_exit(i) /*! wrapper for global exit */ 
/*! \ingroup libgetputgrp */
#define lgp_poll() upc_poll() /*! wrapper for poll which makes sure the network doesn't live lock */



#else

/////////////////////////////////////////////////////////////////
///////                   SHMEM SECTION                 /////////
/////////////////////////////////////////////////////////////////

#include <shmem.h>
#include <mpp/shmem.h>

#define SHARED    /*!< Empty; Should always refer to a symmetrically-allocated array */

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

#define lgp_barrier() shmem_barrier_all() /*!< wrapper for global barrier */
#define lgp_local_part(t,p) (( t * )(p)) /*!< wrapper for localizing a global pointer */
#define lgp_fence() shmem_quiet() /*!< wrapper for memory fence */
#define lgp_all_alloc(nblock_on_all_threads,blocksize) shmem_malloc( ( ( (size_t)(nblock_on_all_threads)+(shmem_n_pes())-1)/(shmem_n_pes()) ) *(blocksize) ) \
  /*!< wrapper for global allocation */
#define lgp_all_free(p) shmem_free(p) /*!< wrapper for a free of global memory */
#define lgp_global_exit(i) shmem_global_exit(i) /*!< wrapper for shutdown and exit */
#define lgp_poll()  /*!< wrapper for a function needed on some UPC implementations */

/*!< macro magic to make lgp_memput and lgp_memget duplicate the array indexing arithmetic of UPC shared arrays. */
#define lgp_memput_shmemwrap(dst,src,n,pe)                             \
  ( (n)==0 ? (void)(0) :                                                \
    (n)%16==0 ? shmem_put128((dst),(src),(n)/16,(pe)) :                 \
    (n)%8==0 ? shmem_put64((dst),(src),(n)/8,(pe)) :                    \
    (n)%4==0 ? shmem_put32((dst),(src),(n)/4,(pe)) :                    \
    shmem_putmem((dst),(src),(n),(pe)) )                                      /*!< macro magic */
#define lgp_memput_bytes_by_pe(dst,src,n,local_offset,pe)              \
  (lgp_memput_shmemwrap( ((char*)(dst))+(local_offset), (src), (n), (pe) ) ) /*!< macro magic */
#define lgp_memput(dst,src,n,index)                                    \
  lgp_memput_bytes_by_pe( (dst), (src), (n),  sizeof(*(dst))*(((size_t)(index))/shmem_n_pes()), (index)%shmem_n_pes() )    /*!< macro magic */


void    lgp_shmem_write_upc_array_int64(SHARED int64_t *addr, size_t index, size_t blocksize, int64_t val); /*!< macro magic */
int64_t lgp_shmem_read_upc_array_int64(const SHARED int64_t *addr, size_t index, size_t blocksize); /*!< macro magic */
#define shmem_int64_p(addr,val,pe) shmem_longlong_p( (long long*)(addr), (val), (pe) ) /*!< macro magic */
#define shmem_int64_g(addr,pe) shmem_longlong_g( (const long long*)(addr), (pe) ) /*!< macro magic */

#define lgp_put_int64(array,index,val) (lgp_shmem_write_upc_array_int64((array),(index),sizeof(int64_t),(val))) /*!< user callable global single word put */
#define lgp_get_int64(array,index) (lgp_shmem_read_upc_array_int64((array),(index),sizeof(int64_t))) /*!< user callable global single word get */

//extern long THREADS; /*!< number of shmem threads */
//extern long MYTHREAD; /*!< shmem thread id */

#endif

/////////////////////////////////////////////////////////////////
///////           END OF UPC/SHMEM SECTION              /////////
/////////////////////////////////////////////////////////////////


#define T0_printf if(MYTHREAD==0) printf /*!< only execute the printf from thread 0 */
#define T0_fprintf if(MYTHREAD==0) fprintf /*!< only execute the fprintf from thread 0 */

void lgp_init(); /*!< function to do initialization if the programming model needs it */

//////////////////////////
// COLLECTIVES
//////////////////////////

int64_t lgp_reduce_add_l(int64_t myval); /*!< collective reduction add on int64_t's */
int64_t lgp_reduce_min_l(int64_t myval); /*!< collective reduction min on int64_t's */
int64_t lgp_reduce_max_l(int64_t myval); /*!< collective reduction max on int64_t's */

double lgp_reduce_add_d(double myval); /*!< collective reduction add on double's */
double lgp_reduce_min_d(double myval); /*!< collective reduction min on double's */
double lgp_reduce_max_d(double myval); /*!< collective reduction max on double's */

int64_t lgp_partial_add_l(int64_t myval); /*!< collective (inclusive) parallel prefix (scan) operation on int64_t's */ 
int64_t lgp_prior_add_l(int64_t myval); /*!< collective (exclusive) parallel prefix (scan) operation on int64_t's */ 

int64_t lgp_min_avg_max_l(minavgmaxL_t *s, int64_t myval, int64_t dem);  /*!< collective reductions for min, avg, max of int64_t's */
int64_t lgp_min_avg_max_d(minavgmaxD_t * s, double myval, int64_t dem);  /*!< collective reductions for min, avg, max of double's */

///////////////////////////
// ATOMICS
///////////////////////////

void       lgp_atomic_add(SHARED int64_t * ptr, int64_t index, int64_t value); /*!< wrapper of atomic add */
int64_t lgp_fetch_and_inc(SHARED int64_t * ptr, int64_t index); /*!< wrapper of atomic fetch and inc */
int64_t lgp_fetch_and_add(SHARED int64_t * ptr, int64_t index, int64_t value); /*!< wrapper of atomic fetch and add */
int64_t  lgp_cmp_and_swap(SHARED int64_t * ptr, int64_t index, int64_t cmp_val, int64_t swap_val); /*!< wrapper of atomic compare and swap */

//Xvoid dump_header();
//Xvoid dump_defines();

double wall_seconds(); /*!< timer */

#define libgetput_INCLUDED  /*!< std trick */
#endif

