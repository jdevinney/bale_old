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
#include "mpp2mpi.h"

long xmpi_my_proc = 0;
long xmpi_n_procs = 0;

int
xmpi_init(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  xmpi_my_proc = rank;
  xmpi_n_procs = size;
  MPI_Barrier(MPI_COMM_WORLD);

  return argc;
}

void
xmpi_exit(void)
{
  MPI_Finalize();
  exit(0);
}

void*
xmpi_alloc_align(size_t align, size_t size)
{
  void* ptr = NULL;
  posix_memalign(&ptr, align, size);
  return ptr;
}

#define REDUCER(Name,Op) \
long \
xmpi_ ## Name ## _long(long myval) \
{ \
  long result; \
  MPI_Allreduce(&myval, &result, 1, MPI_LONG, Op, MPI_COMM_WORLD); \
  return result; \
}

REDUCER(accum, MPI_SUM)
REDUCER(and, MPI_BAND)
REDUCER(or, MPI_BOR)
REDUCER(max, MPI_MAX)
REDUCER(min, MPI_MIN)
#undef REDUCER

void
xmpi_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  if (prefix)
    fprintf(stream, "PE %d: %s> ", (int) xmpi_my_proc, func);
  vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
}

void
xmpi_broadcast64(void* data, size_t count, int root)
{
  MPI_Bcast(data, count, MPI_UINT64_T, root, MPI_COMM_WORLD);
}
