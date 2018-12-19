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


#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "convey_impl.h"
#include "private.h"


#if MPP_NO_MIMD

convey_t*
convey_new_tensor(size_t capacity, size_t item_size, int order, size_t n_local,
                  size_t n_buffers, const convey_alc8r_t* alloc, uint64_t options)
{
  return convey_new_simple(capacity, item_size, alloc, NULL, options);
}

#else

#include "porter.h"
#include "tensor.h"

/*** INTERNAL METHODS ***/

// The next three functions convert the PE number into a tag and
// the destination rank for pushing into porters[0].

static route_t
vector_route(tensor_t* vector, int64_t pe)
{
  return (route_t) { .tag = vector->my_coord[0], .next = pe };
}

static route_t
matrix_route(tensor_t* matrix, int64_t pe)
{
  // dest is (x',y'); we are (x,y); tag is (x',y)
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, matrix->div_local);
  uint32_t lower = dest - matrix->n_local * upper;
  uint32_t tag = (upper << 12) | matrix->my_coord[0];
  return (route_t) { .tag = tag, .next = lower };
}

static route_t
tensor_route(tensor_t* tensor, int64_t pe)
{
  // dest is (x',y',z'); we are (x,y,z); tag is (x',z',z)
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, tensor->div_square);
  uint32_t middle = _divbymul32(dest, tensor->div_local);
  uint32_t lower = dest - tensor->n_local * middle;
  middle -= tensor->n_local * upper;
  uint32_t tag = (upper << 16) | (lower << 8) | tensor->my_coord[0];
  return (route_t) { .tag = tag, .next = middle };
}

// The pivot methods have no side effects except on porters and actimers.

static bool
pivot_mid(tensor_t* matrix, buffer_t* buffer)
{
  // tag is (x',y); we are (x,y'); tag becomes (x,y)
  ACT_START(matrix_pivot);
  const uint32_t my_top = matrix->my_coord[1] << 12;
  const size_t quads = matrix->packet_quads;
  uint32_t* packet = &buffer->data[buffer->start];
  uint32_t* limit = &buffer->data[buffer->limit];
  for (; packet < limit; packet += quads) {
    uint32_t tag = packet[0];
    int dest = tag >> 12;
    tag = my_top | (tag & 0xFFF);
    bool ok = porter_push(matrix->porters[1], tag, &packet[1], dest);
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(matrix_pivot);
      return false;
    }
  }
  buffer->start = buffer->limit;
  ACT_STOP(matrix_pivot);
  return true;
}

static bool
pivot_early(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (x',z',z); we are (x,y,y'); tag becomes (x,z',z)
  ACT_START(tensor_early);
  const uint32_t my_top = tensor->my_coord[2] << 16;
  const size_t quads = tensor->packet_quads;
  uint32_t* packet = &buffer->data[buffer->start];
  uint32_t* limit = &buffer->data[buffer->limit];
  for (; packet < limit; packet += quads) {
    uint32_t tag = packet[0];
    int dest = tag >> 16;
    tag = my_top | (tag & 0xFFFF);
    bool ok = porter_push(tensor->porters[1], tag, &packet[1], dest);
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(tensor_early);
      return false;
    }
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_early);
  return true;
}

static bool
pivot_late(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (x,z',z); we are (x',y',y); tag becomes (x,y,z)
  ACT_START(tensor_late);
  const uint32_t my_mid = tensor->my_coord[0] << 8;
  const size_t quads = tensor->packet_quads;
  uint32_t* packet = &buffer->data[buffer->start];
  uint32_t* limit = &buffer->data[buffer->limit];
  for (; packet < limit; packet += quads) {
    uint32_t tag = packet[0];
    int dest = (tag >> 8) & 0xFF;
    tag = (tag & ~UINT32_C(0xFF00)) | my_mid;
    bool ok = porter_push(tensor->porters[2], tag, &packet[1], dest);
    if (!ok) {
      buffer->start = packet - buffer->data;
      ACT_STOP(tensor_late);
      return false;
    }
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_late);
  return true;
}


/*** TENSOR METHODS ***/

static int
tensor_push(convey_t* self, const void* item, int64_t pe)
{
  tensor_t* tensor = (tensor_t*) self;
  route_t _route = tensor->router(tensor, pe);
  bool ok = porter_push(tensor->porters[0], _route.tag, item, _route.next);
  tensor->stats[convey_PUSHES] += ok;
  return ok ? convey_OK : convey_FAIL;
}

static void*
tensor_upull(convey_t* self, int64_t* from)
{
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  const int order = tensor->order;

  if (buffer && buffer->start == buffer->limit) {
    porter_return(tensor->porters[order - 1]);
    buffer = NULL;
  }
  if (!buffer) {
    buffer = porter_borrow(tensor->porters[order - 1]);
    tensor->buffer = buffer;
    if (!buffer)
      return NULL;
  }

  uint32_t* packet = &buffer->data[buffer->start];
  buffer->start += tensor->packet_quads;
  if (from)
    *from = origin_from_tag(tensor, order, packet[0]);
  tensor->stats[convey_PULLS]++;
  return &packet[1];
}

static int
tensor_pull(convey_t* self, void* item, int64_t* from)
{
  void* source = tensor_upull(self, from);
  if (source == NULL)
    return convey_FAIL;
  memcpy(item, source, self->item_size);
  return convey_OK;
}

static void*
tensor_apull(convey_t* self, int64_t* from)
{
  tensor_t* tensor = (tensor_t*) self;
  void* source = tensor_upull(self, from);
  return pull_pointer(source, tensor->aligned_item, self->item_size);
}

static int
tensor_unpull(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  if (buffer && buffer->start > 0) {
    size_t prev = buffer->start - tensor->packet_quads;
    buffer->start = prev;
    tensor->stats[convey_UNPULLS]++;
    return convey_OK;
  }
  return convey_FAIL;
}

int
tensor_advance(convey_t* self, bool done)
{
  tensor_t* tensor = (tensor_t*) self;
  tensor->stats[convey_ADVANCES]++;
  // We won't be called if state is already COMPLETE

  const int order = tensor->order;
  // No shortcuts: must advance every porter to ensure progress
  for (int k = tensor->n_complete; k < order; k++) {
    bool go = porter_advance(tensor->porters[k], done);
    if (!go) {
      tensor->n_complete++;
      continue;  // done is true
    }
    done = false;
    if (k == order - 1)
      break;

    // Try to interleave work in a reasonable way
    for (int loop = 0; loop < 5; loop++) {
      buffer_t* buffer = porter_borrow(tensor->porters[k]);
      if (!buffer)
        break;
      go = (tensor->pivots[k])(tensor, buffer);
      porter_return(tensor->porters[k]);
      if (!go)
        break;
    }
  }

  if (done)
    for (int stat = 0; stat <= 1; stat++) {
      int64_t sum = 0;
      for (int i = 0; i < order; i++)
        sum += porter_get_stats(tensor->porters[i], stat);
      tensor->stats[stat ? convey_SYNCS : convey_COMMS] = sum;
    }
  return done ? convey_DONE : convey_OK;
}

int
tensor_begin(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  int order = tensor->order;
  bool ok = true;
  for (int i = 0; i < order; i++)
    ok &= porter_setup(tensor->porters[i]);
  tensor->n_complete = 0;
  tensor->stats[convey_BEGINS]++;
  return ok ? convey_OK : convey_error_ALLOC;
}

int
tensor_reset(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  for (int i = tensor->order - 1; i >= 0; i--)
    porter_reset(tensor->porters[i]);
  convey_imp_update_stats(tensor->stats);
  return convey_OK;
}

int
tensor_free(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  for (int i = tensor->order - 1; i >= 0; i--)
    porter_free(tensor->porters[i]);
  free(tensor->aligned_item);
  free(self);
  return convey_OK;
}

int64_t
tensor_statistic(convey_t* self, int which)
{
  tensor_t* tensor = (tensor_t*) self;
  return convey_imp_statistic(tensor->stats, which);
}


/*** METHOD TABLES ***/

static const convey_methods_t tensor_fast_methods = {
  .push = &tensor_push,
  .pull = &tensor_pull,
  .unpull = &tensor_unpull,
  .advance = &tensor_advance,
  .begin = &tensor_begin,
  .reset = &tensor_reset,
  .free = &tensor_free,
  .statistic = &tensor_statistic,
  .apull = &tensor_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};

static const convey_methods_t tensor_debug_methods = {
  ._super_ = &tensor_fast_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &tensor_advance,
  .begin = &tensor_begin,
  .reset = &tensor_reset,
  .free = &tensor_free,
  .statistic = &tensor_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};


/*** PRIVATE CONSTRUCTORS ***/

typedef tensor_t* (make_f)(convey_t*, size_t, size_t, size_t, size_t, size_t,
                           const convey_alc8r_t*, size_t, bool, bool, bool);

tensor_t*
vector_new(convey_t* base, size_t capacity, size_t item_size, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           size_t align, bool dynamic, bool quiet, bool steady)
{
  tensor_t* vector = malloc(sizeof(tensor_t));
  int32_t* friends = malloc(n_procs * sizeof(uint32_t));
  size_t my_proc = MY_PROC;
  bool ok = vector && friends;

  if (ok) {
    *vector = (tensor_t) {
      .convey = *base, .order = 1,
      .my_coord = { my_proc },
      .router = &vector_route,
    };
    for (int i = 0; i < n_procs; i++)
      friends[i] = i;
    vector->porters[0] = porter_new(n_procs, friends, my_proc, item_size, capacity,
                                    n_buffers, alloc, dynamic, false, steady, CONVEY_SEND_1);
  }

  if (!ok || !vector->porters[0]) {
    free(friends);
    free(vector);
    CONVEY_REJECT(quiet, ok ? "construction failed" : "a small malloc() failed");
  }
  return vector;
}

tensor_t*
matrix_new(convey_t* base, size_t capacity, size_t item_size, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           size_t align, bool dynamic, bool quiet, bool steady)
{
  size_t n_rows = n_procs / n_local;
  if (n_local >> 12 || n_rows >> 20)
    CONVEY_REJECT(quiet, "PE counts are too large");

  tensor_t* matrix = malloc(sizeof(tensor_t));
  int32_t* friends[2];
  friends[0] = malloc(n_local * sizeof(uint32_t));
  friends[1] = malloc(n_rows * sizeof(uint32_t));
  if (! (matrix && friends[0] && friends[1])) {
    free(friends[1]);
    free(friends[0]);
    free(matrix);
    CONVEY_REJECT(quiet, "a small malloc() failed");
  }

  size_t my_proc = MY_PROC;
  uint32_t me[2] = { my_proc % n_local, my_proc / n_local };
  *matrix = (tensor_t) {
    .convey = *base, .order = 2,
    .my_coord = { me[0], me[1] }, .n_local = n_local,
    .router = &matrix_route, .pivots[0] = &pivot_mid,
  };
  assert(n_local > 1);
  matrix->div_local = _divbymul32_prep(n_local);

  for (int i = 0; i < n_local; i++)
    friends[0][i] = my_proc + (i - me[0]);
  for (int j = 0; j < n_rows; j++)
    friends[1][j] = my_proc + (j - me[1]) * n_local;
  matrix->porters[0] = porter_new(n_local, friends[0], me[0], item_size, capacity,
                                  n_buffers, alloc, dynamic, true, steady, CONVEY_SEND_0);
  matrix->porters[1] = porter_new(n_rows, friends[1], me[1], item_size, capacity,
                                  n_buffers, alloc, dynamic, false, steady, CONVEY_SEND_1);

  if (! (matrix->porters[0] && matrix->porters[1])) {
    for (int k = 0; k < 2; k++)
      if (!matrix->porters[k])
        free(friends[k]);
    tensor_free(&matrix->convey);
    CONVEY_REJECT(quiet, "construction failed");
  }

  return matrix;
}

tensor_t*
tensor_new(convey_t* base, size_t capacity, size_t item_size, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           size_t align, bool dynamic, bool quiet, bool steady)
{
  size_t squared = n_local * n_local;
  size_t n_middle = 1 + (n_procs - 1) / squared;
  if (n_local >> 8 || n_middle >> 16)
    CONVEY_REJECT(quiet, "PE counts are too large");

  tensor_t* tensor = malloc(sizeof(tensor_t));
  int32_t* friends[3];
  friends[0] = malloc(n_local * sizeof(uint32_t));
  friends[1] = malloc(n_middle * sizeof(uint32_t));
  friends[2] = malloc(n_local * sizeof(uint32_t));
  if (! (tensor && friends[0] && friends[1] && friends[2])) {
    for (int i = 0; i < 3; i++)
      free(friends[i]);
    free(tensor);
    CONVEY_REJECT(quiet, "a small malloc() failed");
  }

  // Compute my coordinates
  size_t my_proc = MY_PROC;
  uint32_t me[3];
  me[0] = my_proc % n_local;
  me[1] = (my_proc / n_local) % n_local;
  me[2] = my_proc / squared;

  *tensor = (tensor_t) {
    .convey = *base, .order = 3,
    .my_coord = { me[0], me[1], me[2] }, .n_local = n_local, .squared = squared,
    .router = &tensor_route, .pivots[0] = &pivot_early, .pivots[1] = &pivot_late,
  };
  assert(n_local > 1);
  tensor->div_local = _divbymul32_prep(n_local);
  tensor->div_square = _divbymul32_prep(squared);

  for (int i = 0; i < n_local; i++)
    friends[0][i] = friends[2][i] = my_proc - me[0] + i;
  for (int j = 0; j < n_middle; j++)
    friends[1][j] = j * squared + n_local * me[0] + me[1];
  tensor->porters[0] = porter_new(n_local, friends[0], me[0], item_size, capacity,
                                  n_buffers, alloc, dynamic, true, steady, CONVEY_SEND_0);
  tensor->porters[1] = porter_new(n_middle, friends[1], me[2], item_size, capacity,
                                  n_buffers, alloc, dynamic, false, steady, CONVEY_SEND_1);
  tensor->porters[2] = porter_new(n_local, friends[2], me[0], item_size, capacity,
                                  n_buffers, alloc, dynamic, true, steady, CONVEY_SEND_2);

  if (! (tensor->porters[0] && tensor->porters[1] && tensor->porters[2])) {
    for (int k = 0; k < 3; k++)
      if (!tensor->porters[k])
        free(friends[k]);
    tensor_free(&tensor->convey);
    CONVEY_REJECT(quiet, "construction failed");
  }

  return tensor;
}


/*** PUBLIC CONSTRUCTOR ***/

convey_t*
convey_new_tensor(size_t capacity, size_t item_size, int order, size_t n_local,
                  size_t n_buffers, const convey_alc8r_t* alloc, uint64_t options)
{
  const bool dynamic = (options & convey_opt_DYNAMIC);
  const bool quiet = (options & convey_opt_QUIET);
  const bool reckless = (options & convey_opt_RECKLESS);
  const bool steady = (options & convey_opt_PROGRESS);

  if (capacity == 0 || item_size == 0 || order < 1 || order > 3 ||
      n_buffers == 0 || (n_buffers & n_buffers-1))
    CONVEY_REJECT(quiet, "invalid arguments");
  if (alloc && (!alloc->grab || !alloc->free))
    CONVEY_REJECT(quiet, "alloc is missing one or both methods");
  if (options &~ (convey_opt_RECKLESS | convey_opt_DYNAMIC | convey_opt_QUIET |
                  convey_opt_PROGRESS | 0xFF * convey_opt_NOALIGN))
    CONVEY_REJECT(quiet, "unrecognized option(s)");

  size_t n_procs = PROCS;
  if (order > 1) {
    if (n_local == 0 || n_procs % n_local)
      CONVEY_REJECT(quiet, "n_local must divide PROCS if order > 1");
    if (order == 3 && n_procs <= n_local * n_local)
      order = 2;
    if (order == 2 && n_procs <= n_local)
      order = 1;
    if (n_local == 1)
      order = 1;
  }

  size_t align = (options / convey_opt_NOALIGN) & 0xFF;
  void* slot = NULL;
  if (! convey_prep_aligned(&slot, item_size, 4, align))
    CONVEY_REJECT(quiet, "a small posix_memalign() failed");

  convey_t _base = {
    ._class_ =  reckless ? &tensor_fast_methods : &tensor_debug_methods,
    .features = steady * convey_STEADY, .suppress = UINT64_C(0) - quiet,
    .item_size = item_size, .n_procs = n_procs, .state = convey_DORMANT,
  };
  make_f* maker = ((make_f* []) { &vector_new, &matrix_new, &tensor_new }) [order - 1];
  tensor_t* tensor = (*maker)(&_base, capacity, item_size, n_procs, n_local,
                              n_buffers, alloc, align, dynamic, quiet, steady);
  if (tensor == NULL) {
    free(slot);
    return NULL;
  }

  tensor->packet_quads = 1 + ((item_size + 3) >> 2);
  tensor->aligned_item = slot;
  tensor->comm = MPP_COMM_CURR;
#if defined(__CRAYXT_COMPUTE_LINUX_TARGET) && !MPP_USE_MPI
  if (!quiet && (!getenv("XT_SYMMETRIC_HEAP_SIZE") &&
                 !getenv("SMA_SYMMETRIC_SIZE") &&
                 !getenv("SMA_SYMMETRIC_PARTITION1")))
    mprint(0, 0, "symmetric heap size not set; expect poor performance\n");
#endif

  return &tensor->convey;
}

#endif
