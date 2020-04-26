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


#include <stdarg.h>
#include <stdio.h>
#if _CRAYC && __UPC__
# include <upc_cray.h>
# include <intrinsics.h>
#endif
#include "mpp2upc.h"


long
upcx_and_long(long myval)
{
  static shared long array[THREADS];
  long result = ~0L;

  upc_barrier(__LINE__);
  array[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  for (int t = MYTHREAD; t < MYTHREAD + THREADS; t++)
    result &= array[(t < THREADS) ? t : t - THREADS];
  upc_barrier(__LINE__);

  return result;
}

long
upcx_or_long(long myval)
{
  static shared long array[THREADS];
  long result = 0;

  upc_barrier(__LINE__);
  array[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  for (int t = MYTHREAD; t < MYTHREAD + THREADS; t++)
    result |= array[(t < THREADS) ? t : t - THREADS];
  upc_barrier(__LINE__);

  return result;
}

long
upcx_accum_long(long myval)
{
  static shared long array[THREADS];
  long sum = 0;

  upc_barrier(__LINE__);
  array[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  for (int t = MYTHREAD; t < MYTHREAD + THREADS; t++)
    sum += array[(t < THREADS) ? t : t - THREADS];
  upc_barrier(__LINE__);

  return sum;
}

long
upcx_max_long(long myval)
{
  static shared long array[THREADS];
  long result = myval;

  upc_barrier(__LINE__);
  array[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  for (int t = MYTHREAD; t < MYTHREAD + THREADS; t++) {
    long value = array[(t < THREADS) ? t : t - THREADS];
    if (result < value)
      result = value;
  }
  upc_barrier(__LINE__);

  return result;
}

long
upcx_min_long(long myval)
{
  static shared long array[THREADS];
  long result = myval;

  upc_barrier(__LINE__);
  array[MYTHREAD] = myval;
  upc_barrier(__LINE__);
  for (int t = MYTHREAD; t < MYTHREAD + THREADS; t++) {
    long value = array[(t < THREADS) ? t : t - THREADS];
    if (result > value)
      result = value;
  }
  upc_barrier(__LINE__);

  return result;
}

void
upcx_all_free(shared void* ptr)
{
#if _CRAYC && __UPC__
  upc_all_free(ptr);
#else
  // This barrier ensures that all threads are finished with the memory.
  upc_barrier(__LINE__);
  if (MYTHREAD == 0 && ptr != NULL)
    upc_free(ptr);
  // This barrier prevents other threads from racing ahead and doing things
  // that conflict with upc_free().  Without it, memory can be corrupted.
  upc_barrier(__LINE__);
#endif
}

void
upcx_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  if (prefix)
    fprintf(stream, "PE %d: %s> ", (int) MYTHREAD, func);
  vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
}

void
upcx_broadcast64(void* data, size_t count, int root)
{
  uint64_t* array = data;
  // FIXME: this implementation is extremely bogus.
  // It only works when count == 32.
  static shared [] uint64_t source[32];
  if (MYTHREAD == root)
    for (int i = 0; i < 32; i++)
      source[i] = array[i];
  upc_barrier(__LINE__);
  for (int i = 0; i < 32; i++)
    array[i] = source[i];
}

int
upcx_alltoallv(shared char* all_recv, shared size_t* all_recv_bytes,
               shared char* all_send, shared size_t* all_send_bytes,
               size_t* offsets, int* perm)
{
  size_t* send_bytes = (size_t*) &all_send_bytes[MYTHREAD];
  size_t* recv_bytes = (size_t*) &all_recv_bytes[MYTHREAD];
  char* send = (char*) &all_send[MYTHREAD];

  upc_barrier(__LINE__);

  for (int i = 0; i < THREADS; i++) {
    size_t thread = perm[i];
    if (send_bytes[thread] > 0)
      upc_memput(&all_recv[THREADS * offsets[MYTHREAD] + thread],
                 &send[offsets[thread]], send_bytes[thread]);
    recv_bytes[thread] = all_send_bytes[THREADS * MYTHREAD + thread];
  }

  upc_barrier(__LINE__);
  return 0;
}
