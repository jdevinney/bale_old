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


#include <stdlib.h>
#include <string.h>

#include "porter_impl.h"
#include "private.h"  // for shmem_free, etc?

#if PORTER_DEBUG
static uint32_t
first_uint32(size_t length, const void* item)
{
  uint32_t first = 0;
  memcpy(&first, item, (length < 4) ? length : 4);
  return first;
}
#endif


/*** Private Methods ***/

static inline buffer_t*
porter_outbuf(porter_t* self, int pe, uint64_t level)
{
  uint32_t* base = self->send_buffers;
  uint64_t index = (pe << self->abundance) + level;
  return (buffer_t*) (base + index * self->buffer_stride);
}

// This function encodes the region delimited by area->next in the buffer
// header, and it updates area->limit to mark the buffer as closed.
// As an invariant, channel->produced counts the number of closed buffers.
static void
porter_close_buffer(porter_t* self, int dest, area_t* area)
{
  area->limit = area->next;
  uint64_t produced = self->channels[dest].produced;
  self->channels[dest].produced = produced + 1;
  uint64_t level = produced & ((UINT64_C(1) << self->abundance) - 1);
  buffer_t* buffer = porter_outbuf(self, dest, level);
  buffer->start = 0;
  buffer->limit = area->next - buffer->data;
}

// This function initiates the sending of a buffer.  A return value of true
// means that the buffer has been delivered, and the caller must record
// this fact.  False means the send is nonblocking; the progress() method
// will eventually force delivery and record the consequences.
static bool
porter_send_buffer(porter_t* self, int dest, uint64_t count, bool last)
{
  // Figure out which buffer we are sending
  const int shift = self->abundance;
  uint64_t level = count & ((UINT64_C(1) << shift) - 1);
  buffer_t* buffer = porter_outbuf(self, dest, level);

  // Work out the buffer size and the signal
  size_t n_bytes = 0;
  uint64_t signal = 2 * count + last;
  if (buffer->limit > 0) {
    n_bytes = sizeof(buffer_t) + buffer->limit * sizeof(uint32_t);
    signal += 2;
  }

  CONVEY_PROF_DECL(_sample);
  CONVEY_PROF_START(&_sample);

  bool arrived = self->_class_->send(self, dest, level, n_bytes, buffer, signal);

  CONVEY_PROF_STOP(&_sample, self->opcode, self->relative[dest], n_bytes + 8);

  return arrived;
}

// The function porter_try_send() tries to send all buffers that have been
// produced, which means they have been closed by porter_close_buffer().
// It then opens a new buffer if possible (the current buffer is closed and
// the next buffer is available) and desirable (either the endgame has not
// begun or the current buffer was sent before the endgame, hence was not
// marked final).
//
// This function may be called with the current buffer open if there are
// urgent buffers to be sent.  For each destination, it must not be called
// again after returning true.

static bool
porter_try_send(porter_t* self, int dest)
{
  self->_class_->progress(self, dest);

  channel_t* channel = &self->channels[dest];
  uint64_t produced = channel->produced;
  uint64_t emitted = channel->emitted;
  uint64_t delivered = channel->delivered;

  // Send as many buffers as we can
  bool final = false;
  while (emitted < produced && self->_class_->ready(self, dest, emitted)) {
    final = self->endgame && (emitted == produced - 1);
    bool arrived = porter_send_buffer(self, dest, emitted, final);
    delivered += arrived;
    if (arrived)
      porter_record_delivery(self, dest, delivered);
    emitted++;  // @emitted
  }
  channel->emitted = emitted;

  // In the endgame, we advance to an empty buffer if we have already
  // emitted all the buffers but did not signal completion.  In this case,
  // because we return false, later calls to porter_try_flush() will close
  // this empty buffer and keep trying to send it (marked as final).

  const uint64_t mask = (UINT64_C(1) << self->abundance) - 1;
  bool drive = !self->endgame || (emitted == produced && !final);
  area_t* area = &self->send_areas[dest];
  if (drive && produced <= delivered + mask && area->next == area->limit) {
    // This code "opens" a new buffer
    buffer_t* buffer = porter_outbuf(self, dest, produced & mask);
    area->next = buffer->data;
    area->limit = (uint32_t*)buffer + self->buffer_quads;
  }

  // Return true if there is nothing left to send
  return final;
}

// A steady porter calls this function periodically for each 'dest' to
// ensure that partially full send buffers don't sit around indefinitely.
// This function is not called during the endgame.
static void
porter_ensure_progress(porter_t* self, int dest)
{
  uint8_t wait = self->waiting[dest];
  if (wait > 0 && wait < PATIENCE) {
    self->waiting[dest] = wait + 1;
    return;
  }

  // Determine whether we have any undelivered data
  channel_t* channel = &self->channels[dest];
  area_t* area = &self->send_areas[dest];
  bool full = (channel->produced > channel->delivered);
  bool partial = false;
  if (area->next < area->limit) {
    uint64_t mask = (UINT64_C(1) << self->abundance) - 1;
    buffer_t* buffer = porter_outbuf(self, dest, channel->produced & mask);
    partial = (area->next != buffer->data);
  }

  if (wait) {
    if (!full && partial)
      porter_close_buffer(self, dest, area);
    porter_try_send(self, dest);
  }
  else if (full || partial)
    self->waiting[dest] = 1;
}

// This function is only called after the endgame has begun.  Then it is
// called repeatedly until it reports, by returning true, that all buffers
// and signals have been delivered.
static bool
porter_try_flush(porter_t* self)
{
  bool done = true;
  int n_ranks = self->n_ranks;
  for (int i = 0; i < n_ranks; i++)
    if (!self->all_sent[i]) {
      area_t* area = &self->send_areas[i];
      if (area->next < area->limit)
        porter_close_buffer(self, i, area);
      bool ok = porter_try_send(self, i);
      DEBUG_PRINT("tried to send to %d, ok = %d\n", i, ok);
      done &= ok;
      self->all_sent[i] = ok;
    }

  // Try to ensure that all emitted buffers are fully delivered
  done &= self->_class_->progress(self, -1);
  return done;
}

// This function is called by porter_advance() when the endgame has not
// begun and there are undelivered urgent buffers.  It does its best to
// deliver them, and in the process updates self->n_urgent.
static void
porter_try_urgent(porter_t* self)
{
  int n_ranks = self->n_ranks;
  for (int i = 0; i < n_ranks; i++)
    if (self->channels[i].urgent > self->channels[i].emitted)
      porter_try_send(self, i);
  self->_class_->progress(self, -1);
}


/*** Methods for Subclasses ***/

void
porter_record_delivery(porter_t* self, int dest, uint64_t delivered)
{
  channel_t* channel = &self->channels[dest];
  if (channel->delivered < channel->urgent &&
      delivered >= channel->urgent)
    self->n_urgent--;
  if (self->waiting)
    self->waiting[dest] = 0;
  channel->delivered = delivered;
}


/*** External Methods ***/

bool
porter_setup(porter_t* self)
{
  long ok = self->_class_->setup(self);
  ok = mpp_and_long(ok);

  if (ok) {
    size_t n = self->n_ranks;
    for (int i = 0; i < n; i++) {
      buffer_t* buffer = porter_outbuf(self, i, 0);
      area_t* area = &self->send_areas[i];
      area->next = buffer->data;
      area->limit = (uint32_t*)buffer + self->buffer_quads;
    }
    for (int i = 0; i < n; i++)
      self->all_sent[i] = (self->relative[i] < 0);
    memset(self->channels, 0, n * sizeof(channel_t));
    if (self->waiting)
      memset(self->waiting, 0, n * sizeof(uint8_t));

    // Current state
    self->endgame = false;
    self->flushed = false;
    self->drained = false;
    self->n_urgent = 0;
    self->phase = self->my_rank;
    self->send_count = 0;
    self->sync_count = 0;
  }

  return ok;
}

// This function is allowed to assume that every item pushed has size
// item_size, and that each buffer holds an integral number of such
// items with no space left over.
bool
porter_push(porter_t* self, uint32_t tag, const void* item, int dest)
{
  area_t* area = &self->send_areas[dest];
  bool room = (area->next < area->limit);
  if (room) {
    DEBUG_PRINT("push  %08x to %u\n", *(uint32_t*)item, dest);
    _prefetch_x(area->next + 24);
    area->next[0] = tag;
    memcpy(area->next + 1, item, self->item_bytes);
    area->next += self->packet_quads;
    if (area->next >= area->limit) {
      porter_close_buffer(self, dest, area);
      porter_try_send(self, dest);
    }
  }
  else
    porter_try_send(self, dest);
  return room;
}

bool
porter_epush(porter_t* self, uint32_t tag, uint32_t descr,
             const void* item, int dest)
{
  area_t* area = &self->send_areas[dest];
  size_t n_quads = PORTER_QUADS(descr);
  bool room = (area->next + n_quads <= area->limit);

  if (room) {
    size_t length = PORTER_BYTES(descr);
    DEBUG_PRINT("epush %zu:%08x to %u\n", length,
                first_uint32(length, item), dest);
    _prefetch_x(area->next + 24);
    area->next[0] = tag;
    area->next[1] = descr;
    memcpy(area->next + 2, item, length);
    area->next += n_quads;
    // Decide whether the buffer is full
    bool ticket = PORTER_TICKET(descr);
    if (ticket) {
      channel_t* channel = &self->channels[dest];
      self->n_urgent += (channel->urgent <= channel->delivered);
      channel->urgent = channel->produced + 1;
    }
    if (ticket || area->next + self->packet_quads >= area->limit) {
      porter_close_buffer(self, dest, area);
      porter_try_send(self, dest);
    }
  }
  else {
    if (area->next < area->limit)
      porter_close_buffer(self, dest, area);
    porter_try_send(self, dest);
  }

  return room;
}

buffer_t*
porter_borrow(porter_t* self)
{
  return self->_class_->borrow(self);
}

void
porter_return(porter_t* self)
{
  self->_class_->turnin(self);
}

bool
porter_advance(porter_t* self, bool done)
{
  if (!self->endgame) {
    if (self->waiting) {
      int phase = self->phase;
      self->phase = (phase+1 < self->n_ranks) ? phase+1 : 0;
      porter_ensure_progress(self, phase);
    }
    if (self->n_urgent > 0)
      porter_try_urgent(self);
    self->endgame = done;
  }

  if (self->endgame && !self->flushed)
    self->flushed = porter_try_flush(self);
  return !self->flushed || !self->drained;
}

void
porter_reset(porter_t* self)
{
  self->_class_->reset(self);
}

int64_t
porter_get_stats(porter_t* self, int which)
{
  if (which == 0)
    return self->send_count;
  else if (which == 1)
    return self->sync_count;
  else
    return 0;
}

// Tear down a non-NULL porter without waiting for other PEs
// and without freeing the 'relative' array.
void
porter_demolish(porter_t* self)
{
  self->_class_->demolish(self);
  free(self->waiting);
  free(self->channels);
  free(self->all_sent);
  free(self->send_areas);
  free(self);
}

void
porter_free(porter_t* self)
{
  if (self == NULL)
    return;
  mpp_barrier(1);  // to be safe
  void* relative = self->relative;
  porter_demolish(self);
  free(relative);  // we absorbed this array
}
