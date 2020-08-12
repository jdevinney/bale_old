#include <assert.h>
#include <string.h>

#include "biconvey_impl.h"
#include "private.h"
#include "simple.h"

typedef struct bisimple bisimple_t;

struct bisimple {
  biconvey_t biconvey;
  convey_t* convey;
  PARALLEL(uint32_t*, replay);
  size_t n_pushed;
  size_t n_pulled;
  size_t n_answers;
};

static int
bisimple_push(biconvey_t* self, const void* query, int64_t to)
{
  bisimple_t* b = (bisimple_t*) self;
  // Cannot push if the replay array is in use by pulls
  if (b->n_pulled < b->n_answers)
    return convey_FAIL;
  int ok = convey_push(b->convey, query, to);
  if (ok == convey_OK) {
    b->replay[b->n_pushed] = (uint32_t) to;
    b->n_pushed += 1;
  }
  return ok;
}

static int
bisimple_pull(biconvey_t* self, void* reply)
{
  bisimple_t* b = (bisimple_t*) self;
  if (b->n_pulled >= b->n_answers)
    return convey_FAIL;

  int64_t from = b->replay[b->n_pulled];
  b->n_pulled += 1;
  simple_t* simple = (simple_t*) b->convey;
  area_t* area = simple->recv + from;
  size_t item_size = b->biconvey.reply_bytes;
  assert(area->next + item_size <= area->limit);
  memcpy(reply, area->next, item_size);
  area->next += item_size;
  return convey_OK;
}

static int
bisimple_advance(biconvey_t* self, bool done)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;

  // If we are done pulling, tell the simple conveyor
  if (b->n_pulled >= b->n_answers)
    simple->pull_from = simple->convey.n_procs;

  // We need to know whether an exchange happened.  We do this in a very
  // hacky way, by watching the communication count.
  int64_t old_stat = simple->stats[convey_COMMS];
  int result = convey_advance(b->convey, done);
  int64_t new_stat = simple->stats[convey_COMMS];
  if (result < 0 || new_stat == old_stat)
    return result;

  // Now we handle all the queries.  The send buffers have been reset.
  const size_t n_procs = simple->convey.n_procs;
  const size_t query_bytes = b->biconvey.query_bytes;
  const size_t reply_bytes = b->biconvey.reply_bytes;
  for (size_t i = 0; i < n_procs; i++) {
    area_t* area = simple->recv + i;
    char* r = simple->send[i].next;
    for (char* q = area->next; q < area->limit; q += query_bytes) {
      b->biconvey.answer(q, r, b->biconvey.context);
      r += reply_bytes;
    }
    simple->send[i].next = r;
    // send[i].limit doesn't matter!
  }

  // Now exchange back, which sets up the recv areas
  int err = simple_alltoallv(simple);
  if (err != convey_OK)
    return err;
  simple_reset_send_buffers(simple);
  b->n_answers = b->n_pushed;
#if 0
  {
    size_t a = b->n_answers;
    fprintf(stderr, "pe %ld says n_answers = %zu\n", MY_PROC, a);
    if (a >= 4)
      fprintf(stderr, "pe %ld replay: %u %u ... %u %u\n", MY_PROC,
	      b->replay[0], b->replay[1], b->replay[a-2], b->replay[a-1]);
  }
#endif
  b->n_pushed = 0;
  b->n_pulled = 0;

  // The result of the first exchange is the correct return value
  return result;
}

static int
bisimple_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;
  int result = convey_begin(b->convey, query_bytes);
  if (result < 0)
    return result;

  if (reply_bytes > query_bytes) {
    simple->capacity = simple->buffer_bytes / reply_bytes;
    if (simple->capacity == 0)
      return convey_error_OFLO;
    simple->buffer_limit = query_bytes * simple->capacity;
    simple_reset_send_buffers(simple);
  }

  size_t max_pushes = simple->convey.n_procs * simple->capacity;
  if (max_pushes > UINT32_MAX)
    return convey_error_TOOBIG;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_ALLOC(b, replay, alloc, max_pushes, uint32_t);
  if (b->replay == NULL)
    return convey_error_ALLOC;
  return convey_OK;
}

static int
bisimple_reset(biconvey_t* self)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_DEALLOC(b, replay, alloc);
  b->replay = NULL;
  return convey_reset(b->convey);
}

static int
bisimple_free(biconvey_t* self)
{
  bisimple_t* b = (bisimple_t*) self;
  if (b->convey->state != convey_DORMANT) {
    int result = bisimple_reset(self);
    if (result < 0)
      return result;
  }

  int result = convey_free(b->convey);
  b->convey = NULL;
  free(b);
  return result;
}


/*** Constructor ***/

static const biconvey_methods_t bisimple_methods = {
  .push = &bisimple_push,
  .pull = &bisimple_pull,
  .advance = &bisimple_advance,
  .begin = &bisimple_begin,
  .reset = &bisimple_reset,
  .free = &bisimple_free,
};

biconvey_t*
biconvey_new_simple(size_t capacity, const convey_alc8r_t* alloc,
		    const convey_mpp_a2a_t* a2a, uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  bisimple_t* b = malloc(sizeof(bisimple_t));
  if (b == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  *b = (bisimple_t) { .biconvey = { ._class_ = &bisimple_methods } };
  options |= convey_opt_PROGRESS;
  b->convey = convey_new_simple(capacity, alloc, a2a, options);
  if (b->convey == NULL) {
    free(b);
    CONVEY_REJECT(quiet, "construction of simple conveyor failed");
  }

  return &b->biconvey;
}
