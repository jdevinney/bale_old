// Copyright (c) 2019, Institute for Defense Analyses
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


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpp2shmem.h"


const mpp_comm_t xshmem_comm_world = { .start = 0, .stride = 1 };


static long source;
static long dest;
// use .+1 instead of max(.,1) for simplicity
static long work_array[SHMEM_REDUCE_MIN_WRKDATA_SIZE + 1];
static long sync_array[SHMEM_REDUCE_SYNC_SIZE];

// Note: we could reduce the number of barriers to 1 by alternating
// between two sets of static areas

#define REDUCER(Name,Call) \
long \
xshmem_ ## Name ## _long(long myval) \
{ \
  for (long i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++) \
    sync_array[i] = SHMEM_SYNC_VALUE; \
  source = myval; \
  shmem_barrier_all(); \
  Call(&dest, &source, 1, 0, 0, shmem_n_pes(), work_array, sync_array); \
  long answer = dest; \
  shmem_barrier_all(); \
  return answer; \
}

REDUCER(accum, shmem_long_sum_to_all)
REDUCER(and, shmem_long_and_to_all)
REDUCER(or, shmem_long_or_to_all)
REDUCER(max, shmem_long_max_to_all)
REDUCER(min, shmem_long_min_to_all)
#undef REDUCER


int
xshmem_init(int argc)
{
  shmem_init();
  return argc;
}

void
xshmem_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  if (prefix)
    fprintf(stream, "PE %d: %s> ", (int) shmem_my_pe(), func);
  vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
}

void
xshmem_broadcast64(void* data, size_t count, int root)
{
  static long psync[SHMEM_BCAST_SYNC_SIZE];
  for (int i = 0; i < SHMEM_BCAST_SYNC_SIZE; i++)
    psync[i] = SHMEM_SYNC_VALUE;
  shmem_barrier_all();
  shmem_broadcast64(data, data, count, root, 0, 0, shmem_n_pes(), psync);
  shmem_barrier_all();
}

int
xshmem_alltoallv(char* recv, size_t recv_bytes[],
                 char* send, size_t send_bytes[],
		 const size_t offsets[], const int perm[])
{
  int me = shmem_my_pe();
  int npes = shmem_n_pes();
  shmem_barrier_all();

  for (int i = 0; i < npes; i++) {
    int pe = perm[i];
    if (send_bytes[pe] > 0)
      shmem_putmem(&recv[offsets[me]], &send[offsets[pe]], send_bytes[pe], pe);
    shmem_putmem(&recv_bytes[me], &send_bytes[pe], sizeof(size_t), pe);
  }

  shmem_barrier_all();
  return 0;
}

int
xshmem_alltoall(char* recv, char* send, size_t n_bytes, const int perm[])
{
  int me = shmem_my_pe();
  int npes = shmem_n_pes();
  shmem_barrier_all();
  for (int i = 0; i < npes; i++) {
    int pe = perm[i];
    shmem_putmem(&recv[me * n_bytes], &send[pe * n_bytes], n_bytes, pe);
  }
  shmem_barrier_all();
  return 0;
}

void
xshmem_final(void)
{
  shmem_finalize();
  exit(0);
}

void
xshmem_exit(int status)
{
  if (status == 0)
    exit(0);
  else
    abort();
}
