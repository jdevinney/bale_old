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


#ifndef CONVEY_PORTER_H
#define CONVEY_PORTER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "convey_alc8r.h"


typedef struct porter porter_t;
typedef struct buffer {
  size_t start;     // first valid item starts at data + start
  size_t limit;     // data + limit is just beyond last valid item
  uint32_t data[];  // each item is preceded by a tag (uint32_t) and
                    // is padded to a multiple of 4 bytes
} buffer_t;

// A variable-length item starts with a 4-byte descriptor which equals
// 2 * (item length in bytes) + ticket_flag.  A 'ticket' ends its buffer
// and must be transmitted ASAP.

#define PORTER_DESCR(bytes, flag) (2*(bytes) | (flag))
#define PORTER_TICKET(descr) ((descr) & 1)
#define PORTER_BYTES(descr) ((descr) >> 1)
#define PORTER_QUADS(descr) (((descr) + 22) >> 3)
// that's an optimized version of (2 + ((descr)/2 + 3)/4)

// Creation of a porter is a collective operation.  Locally, the porter is
// able to communicate with the n given PEs (whose indices are relative to
// the current communicator), and my index among their friends is my_rank.
// (For 0 <= i < n, if PE friends[i] sent me a copy of its 'friends' array,
// then copy[my_rank] would be MY_PROC.)  The 'friends' array is absorbed
// by the porter.
//
// The 'multiplicity' is the number of buffers per friend; it must be
// a power of two.  The 'local' flag indicates whether communication is
// thought to be within a node.  The 'opcode' is used for mpp profiling
// of buffer sends.

porter_t* porter_new(int n, int32_t friends[n], int my_rank,
                     size_t item_size, size_t n_items, size_t multiplicity,
                     const convey_alc8r_t* alloc, bool dynamic, bool local,
                     bool steady, int opcode);

bool      porter_setup(porter_t* self);

// Push a fixed-length item of item_size bytes.  Calls to porter_push and
// porter_epush must not be mixed!  A porter may be used for one kind of
// push or the other, but not both.
bool      porter_push(porter_t* self, uint32_t tag, const void* item, int dest);

// Push a variable-length item, preceded by its tag and descriptor.
bool      porter_epush(porter_t* self, uint32_t tag, uint32_t descr,
                       const void* item, int dest);

// Attempt to retrieve a nonempty buffer of items.
buffer_t* porter_borrow(porter_t* self);

// Give back the most recently taken buffer, having adjusted buffer->start.
void      porter_return(porter_t* self);

bool      porter_advance(porter_t* self, bool done);

// Obtain some statistics (since the last setup).  Currently which=0 means
// buffer sends and which=1 means synchronization operations.
int64_t   porter_get_stats(porter_t* self, int which);

void      porter_reset(porter_t* self);

// Tear down a non-NULL porter without waiting for other PEs
// and without freeing the 'relative' array.
void      porter_demolish(porter_t* self);

void      porter_free(porter_t* self);


#endif
