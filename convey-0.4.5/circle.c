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


// SUMMARY: A simple sorter implementation that keeps a small circular
// buffer of recently pushed items in order to prefetch the buckets they
// belong in.


#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "sorter.h"
#include "bolite.h"


struct sorter {
  area_t* areas;
  char* circular;     // items in transit
  int* targets;       // destination PEs
  size_t item_bytes;
  uint64_t mask;
  // Mutable fields:
  uint64_t head;      // number of items pushed to us
  uint64_t tail;      // number of items we have written out
};

sorter_t*
sorter_new(int n, area_t areas[n], size_t item_bytes, size_t capacity,
           const convey_alc8r_t* alloc, bool dynamic)
{
  if (n == 0 || item_bytes == 0 || capacity == 0 || (capacity & capacity-1))
    return NULL;

  sorter_t* self = malloc(sizeof(sorter_t));
  if (self == 0)
    return NULL;
  *self = (sorter_t) { .areas = areas, .item_bytes = item_bytes, .mask = capacity - 1 };
  self->circular = malloc(capacity * item_bytes);
  self->targets = malloc(capacity * sizeof(int));
  if (!self->circular || !self->targets) {
    sorter_free(self);
    self = NULL;
  }

  return self;
}

bool
sorter_setup(sorter_t* self)
{
  self->head = 0;
  self->tail = 0;
  return true;
}

bool
sorter_push(sorter_t* self, const void* item, int dest)
{
  _prefetch_x(self->areas[dest].next);

  const uint64_t mask = self->mask;
  const size_t size = self->item_bytes;
  uint64_t head = self->head;
  uint64_t tail = self->tail;
  uint64_t index = head & mask;
  char* position = self->circular + index * size;
  bool space = true;

  if (head > tail + mask) {
    int target = self->targets[index];
    area_t* area = &self->areas[target];
    memcpy(area->next, position, size);
    area->next += size;
    space = (area->next < area->limit);
    self->tail = tail + 1;
  }

  self->targets[index] = dest;
  memcpy(position, item, size);
  self->head = head + 1;
  return space;
}

bool
sorter_flush(sorter_t* self)
{
  const size_t size = self->item_bytes;

  while (self->tail < self->head) {
    uint64_t index = self->tail & self->mask;
    int target = self->targets[index];
    area_t* area = &self->areas[target];
    memcpy(area->next, self->circular + index * size, size);
    area->next += size;
    self->tail++;
    if (area->next >= area->limit)
      break;
  }

  return (self->tail == self->head);
}

void
sorter_reset(sorter_t* self)
{
  // nothing to do
}

void
sorter_free(sorter_t* self)
{
  if (self == NULL)
    return;
  free(self->targets);
  free(self->circular);
  free(self);
}
