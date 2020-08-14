// Copyright (c) 2020, Institute for Defense Analyses
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


// SUMMARY: Translate a limited set of mpp_utilV4 functions into trivial
// operations for the single-process case.  Requires that config.h has
// been included (if HAVE_CONFIG_H).

#ifndef CONVEY_MPP2NIL_H
#define CONVEY_MPP2NIL_H

#include <stdio.h>
#include <stdlib.h>

#define mpp_accum_long(L) (L)
#define mpp_and_long(L) (L)
#define mpp_or_long(L) (L)
#define mpp_max_long(L) (L)
#define mpp_min_long(L) (L)

#define mpp_alloc(S) malloc(S)
#define mpp_free(P) free(P)

#if _ISOC11_SOURCE
# define mpp_alloc_align(S,A,B) aligned_alloc(A,S)
#else
# define mpp_alloc_align(S,A,B) xnil_alloc_align(A,S)
#endif

#define mpp_exit(N) exit(N)

#define mpp_broadcast64(P,N,R)
#define mpp_comm_is_equal(X,Y) (1)
#define mpp_comm_is_world(X) (1)
#define mpp_barrier(X)
#define mpp_rel_to_abs_proc(C,I) (I)
#define mpp_util_end() (exit(0))
#define mpp_util_init(ARGC, ARGV, X) (ARGC)
#define mprint(WHO, LEVEL, ...) mfprint(stdout, WHO, 1, __VA_ARGS__)
#define mprint_np(WHO, LEVEL, ...) mfprint(stdout, WHO, 0, __VA_ARGS__)
#define mfprint(FILE, WHO, PREFIX, ...) \
  xnil_mfprint(FILE, PREFIX, __func__, __VA_ARGS__)

#define MY_PROC (0L)
#define PROCS (1L)

typedef void mpp_alltoall_t;
typedef  int mpp_comm_t;
#define MPP_COMM_CURR (0)

// Function prototypes

void* xnil_alloc_align(size_t align, size_t size);
void xnil_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...);

#endif
