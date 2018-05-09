// Copyright (c) 2018, Institute for Defense Analyses,
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
//
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.

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
  return temp ? memcpy(temp, item, item_size) : (void*)item;
}

// Functions that etensor conveyors need to call

int tensor_advance(convey_t* self, bool done);
int tensor_begin(convey_t* self);
int tensor_reset(convey_t* self);
int tensor_free(convey_t* self);
int64_t tensor_statistic(convey_t* self, int which);

#endif
