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


#ifndef CONVEY_TENSOR_H
#define CONVEY_TENSOR_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

typedef struct tensor tensor_t;
typedef bool (pivot_f)(tensor_t*, buffer_t*);
typedef struct { uint32_t tag; int64_t next; } route_t;
typedef route_t (route_f)(tensor_t*, int64_t);

struct tensor {
  convey_t convey;
  int n_complete;
  int order;  // 1=vector, 2=matrix, 3=tensor
  // Fixed data for dealing with routing tags
  uint32_t my_coord[3];
  uint32_t n_local;
  uint32_t squared;
  divbymul32_t div_local;
  divbymul32_t div_square;
  // Other stuff
  size_t packet_quads;
  void* aligned_item;
  buffer_t* buffer;
  porter_t* porters[3];
  route_f* router;
  pivot_f* pivots[2];
  mpp_comm_t comm;
  int64_t stats[convey_imp_N_STATS];
};

static inline int64_t
origin_from_tag(tensor_t* tensor, int order, uint32_t tag)
{
  if (order == 1)
    return tag;
  else if (order == 2)
    return (tag >> 12) * tensor->n_local + (tag & 0xFFF);
  else
    return (tag >> 16) * tensor->squared
      + (tag >> 8 & 0xFF) * tensor->n_local + (tag & 0xFF);
}

static inline void*
pull_pointer(uint32_t* item, void* temp, size_t item_size)
{
  return (temp && item) ? memcpy(temp, item, item_size) : (void*)item;
}

// Functions that etensor conveyors need to call

int tensor_advance(convey_t* self, bool done);
int tensor_begin(convey_t* self);
int tensor_reset(convey_t* self);
int tensor_free(convey_t* self);
int64_t tensor_statistic(convey_t* self, int which);

#endif
