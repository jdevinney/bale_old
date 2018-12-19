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


#ifndef PORTER_IMPL_H
#define PORTER_IMPL_H

#define PORTER_DEBUG 0

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "common.h"
#include "porter.h"

#if PORTER_DEBUG
# define DEBUG_PRINT(...) mprint(MY_PROC, 0, __VA_ARGS__)
#else
# define DEBUG_PRINT(...)
#endif


// FIXME: Describe buffer states here...

// For a "steady" conveyor, which cycles through the ranks to make sure
// all data eventually gets delivered, this is how many cycles to wait
// before moving things along.  It must fit in a uint8_t.
#define PATIENCE 2

typedef struct area {
  uint32_t* next;       // place to write next item
  uint32_t* limit;      // upper limit of the area
} area_t;

// Data kept locally about each outgoing communication channel.
// By definition, produced >= emitted >= delivered.
typedef struct channel {
  uint64_t produced;    // counts outgoing buffers filled and closed
  uint64_t emitted;     // counts outgoing buffers for which sending has begun
  uint64_t delivered;   // counts outgoing buffers actually delivered
  uint64_t urgent;      // if positive, 1 + maximum index of an urgent buffer
  // Must achieve delivered >= urgent in order to send all tickets.
} channel_t;

typedef struct porter_methods porter_methods_t;

struct porter_methods {
  // const porter_methods_t* _super_;
  // Public methods.  Note that setup(), reset(), and demolish() are
  // responsible for allocation and deallocation of send_buffers.
  bool (*setup)(porter_t* self);  
  buffer_t* (*borrow)(porter_t* self);
  void (*turnin)(porter_t* self);
  void (*reset)(porter_t* self);
  void (*demolish)(porter_t* self);

  // Private methods...
  // 0. Check whether the buffer with sequence number n for this
  //    destination (the next one to be sent) is ready to be received.
  bool (*ready)(porter_t* self, int dest, uint64_t n);

  // 1. Send a buffer; return true if the buffer has been delivered.
  bool (*send)(porter_t*, int, uint64_t, size_t, buffer_t*, uint64_t);

  // 2. Ensure progress; send necessary signals; update channels and n_urgent.
  //    If arg >= 0, make progress sending to that destination (at least).
  //    Else try to deliver of all outstanding sends to all destinations.
  //    Return true if delivery was successful.
  bool (*progress)(porter_t*, int);

  // 3. Free any extra data.
  void (*release)(void*);
};

struct porter {
  const porter_methods_t* _class_;
  size_t item_bytes;
  size_t packet_quads;
  size_t buffer_quads;
  size_t buffer_stride;
  int n_ranks;
  int my_rank;          // we write to this position in remote arrays
  int abundance;        // log2(# of buffers per communicant)
  int opcode;           // used for profiling
  // Local memory
  int32_t* relative;    // relative PE numbers
  area_t* send_areas;
  bool* all_sent;       // FIXME: perhaps combine this array with waiting[]
  channel_t* channels;
  uint8_t* waiting;     // tracks how long each link has waited with undelivered data
  // Symmetric memory
  PARALLEL(uint32_t*, send_buffers);
  // State and statistics
  bool endgame;         // no more pushes?
  bool flushed;         // all buffers sent?
  bool drained;         // all buffers received and taken?
  bool dynamic;         // do we alloc on setup and dealloc on reset?
  int n_urgent;         // number of links with undelivered urgent buffers
  int phase;            // (# of advances + my_rank) modulo n_ranks
  int64_t send_count;   // counts sends of buffers
  int64_t sync_count;   // counts calls to quiet and fence
  // Embedded allocator
  convey_alc8r_t alloc;
};

void
porter_record_delivery(porter_t* self, int dest, uint64_t delivered);

#endif
