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


// SUMMARY: Translate a limited set of mpp_utilV4 functions directly
// into MPI.  Requires that config.h has already been included (if
// HAVE_CONFIG_H).

#ifndef CONVEY_MPP2MPI_H
#define CONVEY_MPP2MPI_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define mpp_accum_long(L) xmpi_accum_long(L)
#define mpp_and_long(L) xmpi_and_long(L)
#define mpp_or_long(L) xmpi_or_long(L)
#define mpp_max_long(L) xmpi_max_long(L)
#define mpp_min_long(L) xmpi_min_long(L)

#define mpp_alloc(S) malloc(S)
#define mpp_free(P) free(P)
#if _ISOC11_SOURCE
# define mpp_alloc_align(S,A,B) aligned_alloc(A,S)
#else
# define mpp_alloc_align(S,A,B) xmpi_alloc_align(A,S)
#endif

#define mpp_exit(N) MPI_Abort(MPI_COMM_WORLD,N)

#define mpp_broadcast64(P,N,R) xmpi_broadcast64(P,N,R)
#define mpp_comm_is_equal(X,Y) (1)
#define mpp_comm_is_world(X) (1)
#define mpp_barrier(X) MPI_Barrier(MPI_COMM_WORLD)
#define mpp_rel_to_abs_proc(C,I) (I)
#define mpp_util_end() xmpi_exit()
#define mpp_util_init(ARGC, ARGV, X) xmpi_init(ARGC,ARGV)
#define mprint(WHO, LEVEL, ...) mfprint(stdout, WHO, 1, __VA_ARGS__)
#define mprint_np(WHO, LEVEL, ...) mfprint(stdout, WHO, 0, __VA_ARGS__)
#define mfprint(FILE, WHO, PREFIX, ...) \
  do { if ((WHO) == MY_PROC) xmpi_mfprint(FILE, PREFIX, __func__, __VA_ARGS__); } while (0)

typedef void mpp_alltoall_t;
extern long xmpi_my_proc;
extern long xmpi_n_procs;
#define MY_PROC xmpi_my_proc
#define PROCS xmpi_n_procs
#define mpp_comm_t MPI_Comm
#define MPP_COMM_CURR MPI_COMM_WORLD

// Function prototypes
int xmpi_init(int argc, char* argv[]);
void xmpi_exit(void);
void* xmpi_alloc_align(size_t align, size_t size);
long xmpi_accum_long(long myval);
long xmpi_and_long(long myval);
long xmpi_or_long(long myval);
long xmpi_max_long(long myval);
long xmpi_min_long(long myval);
void xmpi_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...);
void xmpi_broadcast64(void* data, size_t count, int root);

#endif
