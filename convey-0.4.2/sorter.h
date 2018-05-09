// Copyright (c) 2018, Institute for Defense Analyses,
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
//
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.

#ifndef CONVEY_SORTER_H
#define CONVEY_SORTER_H

#include <stdbool.h>
#include <stddef.h>
#include "convey_alc8r.h"


// A sorter is an object that can efficiently distribute a stream of
// fixed-sized items into an array of buffers.  ...

typedef struct sorter sorter_t;

typedef struct area {
  char* next;   // place to write or read next item
  char* limit;  // upper limit of the area
} area_t;


sorter_t* sorter_new(int n, area_t areas[n], size_t item_bytes, size_t capacity,
                     const convey_alc8r_t* alloc, bool dynamic);

// Returns false if something goes wrong (a memory allocation error).
bool sorter_setup(sorter_t* self);

// Returns false if one of the outgoing buffers is full.  The item is
// always successfully pushed, provided that the caller responds to a
// return value of 'false' by emptying the full buffer(s).
bool sorter_push(sorter_t* self, const void* item, int dest);

// Returns true if all pushed items have been distributed into the client's
// buffers.  After a flush, a reset must precede the next push.
bool sorter_flush(sorter_t* self);

void sorter_reset(sorter_t* self);
void sorter_free(sorter_t* self);


#endif
