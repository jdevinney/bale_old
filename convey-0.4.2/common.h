// Copyright (c) 2018, Institute for Defense Analyses,
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
//
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.

// SUMMARY: Include configuration settings and enable actimer.

#ifndef CONVEY_COMMON_H
#define CONVEY_COMMON_H

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if MPP_USE_UPC
# include "mpp2upc.h"
#else
# if !HAVE_MPP_UTIL
#  include "mpp2shmem.h"
#  define MPP_USE_SHMEM 1
# else
#  include "mpp_utilV4.h"
# endif
# define PARALLEL(TYPE,FIELD) TYPE FIELD
# define PARALLEL_ALLOC(OBJECT,FIELD,ALLOC,SIZE,TYPE) \
  (OBJECT)->FIELD = (TYPE*) (ALLOC)->grab((ALLOC)->alc8r, \
    (SIZE) * sizeof(TYPE), __FILE__, __LINE__)
# define PARALLEL_DEALLOC(OBJECT,FIELD,ALLOC) \
  (ALLOC)->free((ALLOC)->alc8r, (OBJECT)->FIELD)
#endif

#if ENABLE_PROFILING
// this is redundant because mpp_utilV4.h includes it:
# include "mpp_utilV4_profile.h"
# define CONVEY_PROF_DECL(x) mpp_profile_t x
# define CONVEY_PROF_START mpp_profile_start
# define CONVEY_PROF_STOP mpp_profile_stop
# define CONVEY_SEND_0 PROF_OP_USER(0)
# define CONVEY_SEND_1 PROF_OP_USER_P2P
# define CONVEY_SEND_2 PROF_OP_USER(16)
#else
# define CONVEY_PROF_DECL(x)
# define CONVEY_PROF_START(sample)
# define CONVEY_PROF_STOP(sample,opcode,other_pe,data)
# define CONVEY_SEND_0 (0)
# define CONVEY_SEND_1 (1)
# define CONVEY_SEND_2 (2)
#endif

#if HAVE_MPP_UTIL
# define ACTIMER_MODULE_NAME conveyors
# define ACTIMER_SHARED
# include "actimer.h"
// Use detail level 3 for now, to tidy up mpp_util reporting
# define ACT_START(timer) actimer_start(timer, 3)
# define ACT_STOP(timer) actimer_stop(timer, 3, 1)
#else
# define ACT_START(timer)
# define ACT_STOP(timer)
#endif

#endif
