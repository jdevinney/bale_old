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

// Returns 0 if there is nothing to flush, -1 if an outgoing buffer
// filled up, and 1 otherwise.
int sorter_flush(sorter_t* self);

void sorter_reset(sorter_t* self);
void sorter_free(sorter_t* self);


#endif
