// Copyright (c) 2018, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the copyright holder nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.


// SUMMARY: Include configuration settings and enable actimer.

#ifndef CONVEY_COMMON_H
#define CONVEY_COMMON_H

#if HAVE_CONFIG_H
#include "config.h"
#endif
#if HAVE_STDATOMIC_H
# include <stdatomic.h>
#endif

#define MATRIX_REMOTE_HOP 1
#if (MPP_USE_UPC || HAVE_SHMEM_PTR) && HAVE_STDATOMIC_H && \
    HAVE__ATOMIC_UINT64_T && (ATOMIC_LLONG_LOCK_FREE >= 2)
# define CONVEY_INTRANODE 1
#endif


#if MPP_USE_UPC
# include "mpp2upc.h"
#else
# if HAVE_MPP_UTIL
#  include "mpp_utilV4.h"
# elif MPP_RAW_MPI
#  include "mpp2mpi.h"
#  define MPP_USE_MPI 1
# else
#  include "mpp2shmem.h"
#  define MPP_USE_SHMEM 1
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

#if MPP_USE_MPI
# if HAVE_MPP_UTIL
#  define mpp_comm_mpi(C) ((C).internal->mpi_comm)
# else
#  define mpp_comm_mpi(C) C
# endif
#endif

#endif
