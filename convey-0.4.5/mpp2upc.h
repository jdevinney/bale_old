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


// SUMMARY: Emulate a limited set of mpp_utilV4 and shmem functions within
// UPC.  Requires that config.h already be included (if HAVE_CONFIG_H).

#ifndef CONVEY_MPP2UPC_H
#define CONVEY_MPP2UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <upc.h>
#if HAVE_UPC_CASTABLE_H
# include <upc_castable.h>
#endif

// A structure member that is a PARALLEL pointer provides access to a
// shared array.  The member itself is a ordinary pointer to the local
// slice of the shared array, but the corresponding shared pointer is also
// present as a member whose name is prefixed by all_.
#define PARALLEL(TYPE,FIELD) TYPE FIELD; shared TYPE all_ ## FIELD
#define PARALLEL_ALLOC(OBJECT,FIELD,ALLOC,SIZE,TYPE) \
  (OBJECT)->all_ ## FIELD = upc_all_alloc(THREADS * (SIZE), sizeof(TYPE)); \
  (OBJECT)->FIELD = ((OBJECT)->all_ ## FIELD == NULL) ? NULL : \
    (TYPE*) & (OBJECT)->all_ ## FIELD [MYTHREAD]
#define PARALLEL_DEALLOC(OBJECT,FIELD,ALLOC) \
  upcx_all_free((OBJECT)->all_ ## FIELD); \
  (OBJECT)->all_ ## FIELD = NULL

#define mpp_accum_long(L) upcx_accum_long(L)
#define mpp_and_long(L) upcx_and_long(L)
#define mpp_or_long(L) upcx_or_long(L)
#define mpp_max_long(L) upcx_max_long(L)
#define mpp_alloc(S) malloc(S)
#define mpp_alloc_align(S,A,B) malloc(S)
#define mpp_free(P) free(P)
#define mpp_broadcast64(P,N,R) upcx_broadcast64(P,N,R)
#define mpp_comm_is_equal(X,Y) (1)
#define mpp_comm_is_world(X) (1)
#define mpp_barrier(X) upc_barrier
#define mpp_exit(N) upc_global_exit(N)
#define mpp_fence() upc_fence
#define mpp_rel_to_abs_proc(C,I) (I)
#define mpp_util_end() exit(0)
#define mpp_util_init(ARGC, ARGV, X) (ARGC)
#define mprint(WHO, LEVEL, ...) mfprint(stdout, WHO, 1, __VA_ARGS__)
#define mprint_np(WHO, LEVEL, ...) mfprint(stdout, WHO, 0, __VA_ARGS__)
#define mfprint(FILE, WHO, PREFIX, ...) \
  do { if ((WHO) == MYTHREAD) upcx_mfprint(FILE, PREFIX, __func__, __VA_ARGS__); } while (0)

#define MY_PROC MYTHREAD
#define PROCS THREADS
#define MPP_COMM_CURR (0)

#define shmem_fence() upc_fence
#define shmem_malloc(S) malloc(S)
#define shmem_free(P) free(P)
#define shmem_ptr(P) upc_cast(P)

typedef void mpp_alltoall_t;
typedef int mpp_comm_t;

// These reduce-to-all functions usually have names like upc_all_rdc_add_l
long upcx_accum_long(long myval);
long upcx_and_long(long myval);
long upcx_or_long(long myval);
long upcx_max_long(long myval);

// This collective deallocator is usually called upc_all_free()
void upcx_all_free(shared void* ptr);

// Need minimal support for printing errors
void upcx_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...);

// Broadcast function...
void upcx_broadcast64(void* data, size_t count, int root);

// upcx_alltoallv: For each pair (i,j), thread i sends a block of data from
// its local slice of send[], starting at byte offset offsets[i], to the
// part of recv[] local to thread j starting at byte offset offsets[j].
// (The offsets[] array is local but should be identical across threads.)
// The number of bytes sent is given by send_bytes[THREADS*j + i], and
// recv_bytes[THREADS*i + j] is set to this value.

int upcx_alltoallv(shared char* recv, shared size_t* recv_bytes,
                   shared char* send, shared size_t* send_bytes,
                   size_t* offsets, int* perm);

#endif
