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
#include "private.h"

#if HAVE_STDATOMIC_H && HAVE_USABLE_ATOMICS && (__STDC_NO_ATOMICS__ != 1)
# include <stdatomic.h>
#else
# define atomic_load(p) (*(p))
#endif


// For SHMEM puts, we must use long long to transmit 64-bit integers.
typedef struct put_porter {
  porter_t porter;
  // Local memory: indexed by rank (0 <= rank < n_ranks)
  int32_t* friends;     // absolute PE numbers
  long long* disposed;  // counts incoming buffers received and disposed
  void* extra;          // extra data needed by a specialized porter
  // Symmetric memory: indexed by rank and buffer level, 0 <= level < (1<<abundance)
  PARALLEL(uint32_t*, recv_buffers);
  PARALLEL(uint64_t*, received);   // counts incoming buffers (indexed by PE only)
  PARALLEL(long long*, consumed);  // counts times each outgoing buffer has been consumed
  // Current state
  int n_pending;        // count of buffers waiting to be taken
  int i_pending;        // which of these to consider next
  int32_t* pending;     // sources of pending buffers (first n_pending valid)
  int i_taken;          // source of most recently taken buffer, if any
  buffer_t* taken;      // buffer borrowed most recently
} put_porter_t;


static inline buffer_t*
porter_inbuf(porter_t* self, int pe, uint64_t level)
{
  put_porter_t* putp = (put_porter_t*) self;
  uint32_t* base = putp->recv_buffers;
  uint64_t index = (pe << self->abundance) + level;
  return (buffer_t*) (base + index * self->buffer_stride);
}


/*** Memory Functions ***/

static void*
symmetric_alloc(void* alc8r, size_t size, const char* tag, uint64_t value)
{
  return shmem_malloc(size);
}

static void
symmetric_free(void* alc8r, void* ptr)
{
  if (ptr)
    shmem_free(ptr);
}

static const convey_alc8r_t porter_alc8r = {
  .alc8r = NULL, .grab = &symmetric_alloc, .free = &symmetric_free,
};

static bool
porter_grab_buffers(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  convey_alc8r_t* alloc = &self->alloc;
  size_t size = self->buffer_stride * self->n_ranks << self->abundance;
  PARALLEL_ALLOC(self, send_buffers, alloc, size, uint32_t);
  PARALLEL_ALLOC(putp, recv_buffers, alloc, size, uint32_t);
  return self->send_buffers && putp->recv_buffers;
}

static void
porter_free_buffers(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  convey_alc8r_t* alloc = &self->alloc;
  PARALLEL_DEALLOC(putp, recv_buffers, alloc);
  PARALLEL_DEALLOC(self, send_buffers, alloc);
  putp->recv_buffers = NULL;
  self->send_buffers = NULL;
}


/*** Functions for Blocking Porters ***/


/*** Private Functions ***/

static void
putp_scan_receipts(put_porter_t* putp)
{
  porter_t* self = &putp->porter;
  const int n_ranks = self->n_ranks;
  bool stuck[n_ranks];
  memset(stuck, 0, n_ranks);

  // Compress the current pending list and create set
  int n = putp->n_pending, k = 0;
  for (int i = 0; i < n; i++) {
    int32_t source = putp->pending[i];
    if (source >= 0) {
      putp->pending[k++] = source;
      stuck[source] = true;
    }
  }
  int n_stuck = k;

  // Traverse the incoming buffers (perhaps in a better order?)
  int n_drained = 0;
  for (int source = 0; source < n_ranks; source++) {
    if (stuck[source])
      continue;
    long long disposed = putp->disposed[source];
    uint64_t received = atomic_load(putp->received + source);
    if ((received >> 1) > disposed)
      putp->pending[k++] = source;
    else
      n_drained += (received & 1);
  }

  self->drained = (n_drained == self->n_ranks);
  putp->n_pending = k;
  // Work on new buffers if possible, otherwise start over?
  // self->i_pending = (n_stuck < k) ? n_stuck : 0;
  putp->i_pending = 0;
}


/*** External Methods ***/

static bool
putp_setup(porter_t* self)
{
  bool ok = !self->dynamic || porter_grab_buffers(self);

  if (ok) {
    put_porter_t* putp = (put_porter_t*) self;
    size_t n = self->n_ranks;
    memset(putp->disposed, 0, n * sizeof(long long));
    // Mark receipts from invalid PEs as drained
    for (int i = 0; i < n; i++)
      putp->received[i] = (putp->friends[i] < 0);
    const int shift = self->abundance;
    memset(putp->consumed, 0, (n << shift) * sizeof(long long));

    // Current state
    putp->n_pending = 0;
    putp->i_pending = 0;
    putp->taken = NULL;
  }

  // No barrier needed, because caller will reduce the result
  return ok;
}

static buffer_t*
putp_borrow(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  if (putp->i_pending == putp->n_pending)
    putp_scan_receipts(putp);
  int i = putp->i_pending;
  if (i == putp->n_pending)
    return NULL;

  int source = putp->pending[i];
  putp->i_taken = source;
  uint64_t mask = (1L << self->abundance) - 1;
  uint64_t level = putp->disposed[source] & mask;
  putp->taken = porter_inbuf(self, source, level);
  return putp->taken;
}

static void
putp_return(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  int i = putp->i_pending++;
  buffer_t* taken = putp->taken;

  if (taken->start >= taken->limit) {
    putp->pending[i] = -1;
    int source = putp->i_taken;
    long long disposed = putp->disposed[source]++;  // @disposed
    const int shift = self->abundance;
    uint64_t level = disposed & ((1L << shift) - 1);

    CONVEY_PROF_DECL(_sample);
    CONVEY_PROF_START(&_sample);
    uint64_t slot = (self->my_rank << shift) + level;
    long long value = (disposed >> shift) + 1;
    int32_t pe = putp->friends[source];
#if MPP_USE_SHMEM
    shmem_longlong_p(&putp->consumed[slot], value, pe);  // @consumed
#else
    putp->all_consumed[THREADS * slot + pe] = value;
#endif
    CONVEY_PROF_STOP(&_sample, PROF_OP_PUT, self->relative[source], sizeof(long long));
  }
}

static void
putp_reset(porter_t* self)
{
  // Acknowledgements could still be in flight!
  mpp_barrier(1);
  if (self->dynamic)
    porter_free_buffers(self);
}

static void
putp_demolish(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  if (putp->extra)
    (*self->_class_->release)(putp->extra);
  porter_free_buffers(self);
  convey_alc8r_t* alloc = &self->alloc;
  PARALLEL_DEALLOC(putp, consumed, alloc);
  PARALLEL_DEALLOC(putp, received, alloc);
  free(putp->pending);
  free(putp->disposed);
  free(putp->friends);
}


/*** Methods for Standard Porters ***/

static bool
standard_ready(porter_t* self, int dest, uint64_t count)
{
  put_porter_t* putp = (put_porter_t*) self;
  const int shift = self->abundance;
  const uint64_t mask = (UINT64_C(1) << shift) - 1;
  long long* consumed = &putp->consumed[dest << shift];
  return (consumed[count & mask] >= (count >> shift));
}

static bool
standard_send(porter_t* self, int dest, uint64_t level, size_t n_bytes,
              buffer_t* buffer, uint64_t signal)
{
  put_porter_t* putp = (put_porter_t*) self;
  const int rank = self->my_rank;
  const int pe = putp->friends[dest];

  if (n_bytes > 0) {
#if MPP_USE_SHMEM
    buffer_t* remote = porter_inbuf(self, rank, level);
    shmem_putmem(remote, buffer, n_bytes, pe);
#else
    uint64_t offset = ((rank << self->abundance) + level) * self->buffer_stride;
    upc_memput(&putp->all_recv_buffers[THREADS * offset + pe], buffer, n_bytes);
#endif
    self->send_count++;
  }

  // Ensure previous buffers *and signals* have been delivered
  self->sync_count++;
  shmem_fence();
#if MPP_USE_SHMEM
  shmem_put64(&putp->received[rank], &signal, 1, pe);  // @received
#else
  putp->all_received[THREADS * rank + pe] = signal;
#endif

  return true;
}

static bool
standard_progress(porter_t* self, int dest)
{
  return true;
}

static const porter_methods_t standard_methods = {
  .setup = &putp_setup,
  .borrow = &putp_borrow,
  .turnin = &putp_return,
  .reset = &putp_reset,
  .demolish = &putp_demolish,
  .ready = &standard_ready,
  .send = &standard_send,
  .progress = &standard_progress,
  // release() won't be called, as extra == NULL
};


/*** Methods for Local Porters ***/

// Try to ensure that _Atomic uint64_t not only exists but is lock-free.
#if (MPP_USE_UPC || HAVE_SHMEM_PTR) && HAVE_STDATOMIC_H && \
    HAVE__ATOMIC_UINT64_T && (ATOMIC_LLONG_LOCK_FREE >= 2)

typedef struct nbrhood {
  // The following arrays of pointers have length n_ranks and are indexed by rank
  uint32_t** buffer_ptrs;  // remote 'recv_buffers' translated by shmem_ptr()
  uint64_t** signal_ptrs;  // remote 'received' translated by shmem_ptr()
} nbrhood_t;

static bool
nbrhood_init(put_porter_t* putp, nbrhood_t* nbrhood)
{
  int n = putp->porter.n_ranks;
  bool ok = true;
  for (int i = 0; ok && i < n; i++) {
    int pe = putp->friends[i];
    if (pe >= 0) {
      uint64_t* p = shmem_ptr(putp->received, pe);
      nbrhood->signal_ptrs[i] = p;
      uint32_t* q = shmem_ptr(putp->recv_buffers, pe);
      nbrhood->buffer_ptrs[i] = q;
      ok &= (p != NULL) & (q != NULL);
    }
  }
  return ok;
}

static bool
local_setup(porter_t* self)
{
  put_porter_t* putp = (put_porter_t*) self;
  return putp_setup(self) &&
    (!self->dynamic || nbrhood_init(putp, putp->extra));
}

static bool
local_send(porter_t* self, int dest, uint64_t level, size_t n_bytes,
           buffer_t* buffer, uint64_t signal)
{
  const int rank = self->my_rank;
  const nbrhood_t* nbrhood = ((put_porter_t*)self)->extra;

  // Need local address of remote receive buffers
  if (n_bytes > 0) {
    uint32_t* remote = nbrhood->buffer_ptrs[dest];
    uint64_t index = (rank << self->abundance) + level;
    remote += index * self->buffer_stride;
    memcpy(remote, buffer, n_bytes);
    self->send_count++;
  }

  // Need local address of remote 'received' array
  uint64_t* notify = nbrhood->signal_ptrs[dest] + rank;
  atomic_store(notify, signal);
  return true;
}

static void
local_release(void* extra)
{
  if (extra == NULL)
    return;
  nbrhood_t* nbrhood = extra;
  free(nbrhood->signal_ptrs);
  free(nbrhood->buffer_ptrs);
  free(nbrhood);
}

static const porter_methods_t local_methods = {
  .setup = &local_setup,
  .borrow = &putp_borrow,
  .turnin = &putp_return,
  .reset = &putp_reset,
  .demolish = &putp_demolish,
  .ready = &standard_ready,
  .send = &local_send,
  .progress = &standard_progress,
  .release = &local_release,
};


static const porter_methods_t*
local_prepare(porter_t* self)
{
  // Make sure we have buffers to check access
  if (self->dynamic)
    porter_grab_buffers(self);

  put_porter_t* putp = (put_porter_t*) self;
  int n = self->n_ranks;
  nbrhood_t* nbrhood = malloc(sizeof(nbrhood_t));
  bool ok = (nbrhood != NULL);
  if (ok) {
    nbrhood->buffer_ptrs = malloc(n * sizeof(uint32_t*));
    nbrhood->signal_ptrs = malloc(n * sizeof(uint64_t*));
    ok = nbrhood->buffer_ptrs && nbrhood->signal_ptrs;
    if (ok)
      ok = nbrhood_init(putp, nbrhood);
  }
  ok = mpp_and_long(ok);
  
  if (ok)
    putp->extra = nbrhood;
  else {
    local_release(nbrhood);
    if (putp->friends[self->my_rank] == 0)
      mprint(MY_PROC, 0, "WARNING: failed to set up local porter(%d)\n", n);
  }

  if (self->dynamic)
    porter_free_buffers(self);
  return (ok ? &local_methods : NULL);
}

#else

static const porter_methods_t*
local_prepare(porter_t* self) { return NULL; }

#endif


/*** Functions for Nonblocking Porters ***/

#if HAVE_SHMEM_PUTMEM_NBI

// Maintain (in putp->extra) a vector of the signals that need to be sent.

static bool
nonblock_send(porter_t* self, int dest, uint64_t level, size_t n_bytes,
                 buffer_t* buffer, uint64_t signal)
{
  put_porter_t* putp = (put_porter_t*) self;
  if (n_bytes > 0) {
    const int rank = self->my_rank;
    const int pe = putp->friends[dest];
    buffer_t* remote = porter_inbuf(self, rank, level);
    DEBUG_PRINT("%zu quads to %d, signal = %lu\n", buffer->limit - buffer->start, pe, signal);
    shmem_putmem_nbi(remote, buffer, n_bytes, pe);
    self->send_count++;
  }
  else {
    DEBUG_PRINT("0 quads to %d, signal = %lu\n", putp->friends[dest], signal);
  }

  uint64_t* inflight = putp->extra;
  inflight[dest] = signal;
  return false;
}

static bool
nonblock_progress(porter_t* self, int dest)
{
  // Decide whether it's time to force delivery and send signals.
  // Do this when we have emitted half of our buffers on this channel.
  // [1 buffer => diff > 0; 2 buffers => diff > 0; 4 buffers => diff > 1]
  if (dest >= 0) {
    uint64_t limit = ((UINT64_C(1) << self->abundance) - 1) >> 1;
    channel_t* channel = &self->channels[dest];
    if (channel->emitted <= channel->delivered + limit &&
        channel->urgent <= channel->delivered &&
        (self->waiting == NULL || self->waiting[dest] < PATIENCE))
      return false;
  }

  // Force delivery of all puts
  mpp_quiet();
  self->sync_count++;

  const int n = self->n_ranks;
  const int rank = self->my_rank;
  put_porter_t* putp = (put_porter_t*) self;
  uint64_t* inflight = putp->extra;

  // Update delivery information and send signals
  for (dest = 0; dest < n; dest++) {
    uint64_t signal = inflight[dest];
    if (signal) {
      channel_t* channel = &self->channels[dest];
      porter_record_delivery(self, dest, channel->emitted);
      int pe = putp->friends[dest];
      shmem_put64(&putp->received[rank], &signal, 1, pe);
      DEBUG_PRINT("sent signal %lu to %d\n", signal, pe);
      inflight[dest] = 0;
    }
  }
  return true;
}

static const porter_methods_t nonblock_methods = {
  .setup = &putp_setup,
  .borrow = &putp_borrow,
  .turnin = &putp_return,
  .reset = &putp_reset,
  .demolish = &putp_demolish,
  .ready = &standard_ready,
  .send = &nonblock_send,
  .progress = &nonblock_progress,
  .release = &free,
};

static const porter_methods_t*
nonblock_prepare(porter_t* self)
{
  int n = self->n_ranks;
  ((put_porter_t*)self)->extra = calloc(n, sizeof(uint64_t));
  return &nonblock_methods;
}

#else

static const porter_methods_t*
nonblock_prepare(porter_t* self) { return NULL; }

#endif


/*** Constructor ***/

porter_t*
porter_new(int n, int32_t relative[n], int my_rank,
           size_t item_size, size_t n_items, size_t multiplicity,
           const convey_alc8r_t* alloc, bool dynamic, bool local,
           bool steady, int opcode)
{
  if (n <= 0 || multiplicity == 0 || (multiplicity & multiplicity-1))
    return NULL;

  put_porter_t* putp = malloc(sizeof(put_porter_t));
  if (putp == NULL)
    return NULL;
  porter_t* porter = &putp->porter;

  if (alloc == NULL)
    alloc = &porter_alc8r;

  const size_t packet_quads = 1 + ((item_size + 3) >> 2);
  const size_t buffer_quads = (sizeof(buffer_t) >> 2) + n_items * packet_quads;
  const size_t buffer_stride = buffer_quads + (buffer_quads & 1);
  *putp = (put_porter_t) { .porter = {
      .item_bytes = item_size, .packet_quads = packet_quads,
      .buffer_quads = buffer_quads, .buffer_stride = buffer_stride,
      .n_ranks = n, .my_rank = my_rank, .abundance = _trailz(multiplicity),
      .opcode = opcode, .relative = relative, .dynamic = dynamic, .alloc = *alloc,
    } };

  const size_t m = multiplicity;
  // Local allocations
  porter->send_areas = malloc(n * sizeof(area_t));
  porter->all_sent = malloc(n * sizeof(bool));
  porter->channels = malloc(n * sizeof(channel_t));
  if (steady)
    porter->waiting = malloc(n * sizeof(uint8_t));
  putp->friends = malloc(n * sizeof(int32_t));
  putp->disposed = malloc(n * sizeof(long long));
  putp->pending = malloc(n * sizeof(int32_t));
  // Symmetric allocations
  PARALLEL_ALLOC(putp, received, alloc, n, uint64_t);
  PARALLEL_ALLOC(putp, consumed, alloc, n * m, long long);
  bool ok = (porter->send_areas && porter->all_sent && porter->channels &&
             (!steady || porter->waiting) && putp->friends && putp->disposed &&
             putp->pending && putp->received && putp->consumed);
  if (ok && !dynamic)
    ok = porter_grab_buffers(porter);

  ok = mpp_and_long(ok);
  if (!ok) {
    porter_demolish(porter);
    return NULL;
  }
  // porter_setup() will erase everything

  int64_t n_procs = PROCS;
  for (int i = 0; i < n; i++)
    if (relative[i] < 0 || relative[i] >= n_procs)
      putp->friends[i] = relative[i] = -1;
    else
      putp->friends[i] = mpp_rel_to_abs_proc(MPP_COMM_CURR, relative[i]);

  // Decide which subclass to use
  const porter_methods_t* methods =
    local ? local_prepare(porter) : nonblock_prepare(porter);
  porter->_class_ = (methods != NULL) ? methods : &standard_methods;

  return porter;
}
