#include <stdint.h>
#include <string.h>

#include "biconvey_impl.h"
#include "private.h"

typedef struct bitensor bitensor_t;
typedef struct packet packet_t;

struct packet {
  uint32_t token;
  char item[];
};

struct bitensor {
  biconvey_t biconvey;
  size_t reorder_bytes;
  PARALLEL(char*, reorder);
  uint64_t* present;
  packet_t* query;
  packet_t* reply;  // reply->item should be strongly aligned
  uint32_t await;
  uint32_t inflight;
  uint32_t limit;
  bool dynamic;
  convey_alc8r_t alloc;
};


/*** Working Methods ***/

static inline uint32_t
easy_mod(uint32_t x, uint32_t m)
{
  return (x >= m) ? x - m : x;
}

static int
bitensor_push(biconvey_t* self, const void* query, int64_t to)
{
  bitensor_t* b = (bitensor_t*) self;
  uint32_t inflight = b->inflight;
  if (inflight + 1 == b->limit)
    return convey_FAIL;

  convey_t* q = b->biconvey.queries;
  packet_t* packet = b->query;
  packet->token = easy_mod(b->await + inflight, b->limit);
  memcpy(packet->item, query, q->item_size);

  int result = convey_push(q, packet, to);
  if (result == convey_OK)
    b->inflight = inflight + 1;
  return result;
}

static int
bitensor_pull(biconvey_t* self, void* item)
{
  bitensor_t* b = (bitensor_t*) self;
  uint32_t await = b->await;
  uint64_t word = b->present[await >> 6];
  uint64_t bit = UINT64_C(1) << (await & 63);
  if (! (word & bit))
    return convey_FAIL;

  b->present[await >> 6] = word & ~bit;
  b->await = easy_mod(await + 1, b->limit);
  b->inflight -= 1;
  size_t item_size = b->biconvey.replies->item_size;
  memcpy(item, b->reorder + await * item_size, item_size);
  return convey_OK;
}

static int
bitensor_advance(biconvey_t* self, bool done)
{
  bitensor_t* b = (bitensor_t*) self;
  convey_t* q = b->biconvey.queries;
  convey_t* r = b->biconvey.replies;

  int result = convey_advance(q, done);
  if (result < 0)
    return result;
  done = (result == convey_DONE);
  result = convey_advance(r, done);
  if (result < 0)
    return result;

  if (!done) {
    int64_t from;
    packet_t* query;
    packet_t* reply = b->reply;
    while ((query = convey_apull(q, &from)) != NULL) {
      reply->token = query->token;
      b->biconvey.answer(query->item, reply->item, b->biconvey.context);
      // FIXME: can we ensure this push will succeed?
      if (convey_push(r, reply, from) != convey_OK) {
	convey_unpull(q);
	break;
      }
    }
  }

  packet_t* reply;
  while ((reply = convey_apull(r, NULL)) != NULL) {
    uint32_t token = reply->token;
    b->present[token >> 6] |= UINT64_C(1) << (token & 63);
    size_t item_size = r->item_size;
    memcpy(b->reorder + token * item_size, reply->item, item_size);
  }

  if (result == convey_DONE && b->inflight)
    result = convey_NEAR;
  return result;
}


/*** Setup and Teardown ***/

static bool
alloc_reorder(bitensor_t* b)
{
  convey_alc8r_t* alloc = &b->alloc;
  PARALLEL_ALLOC(b, reorder, alloc, b->reorder_bytes, char);
  return (b->reorder != NULL);
}

static void
dealloc_reorder(bitensor_t* b)
{
  convey_alc8r_t* alloc = &b->alloc;
  PARALLEL_DEALLOC(b, reorder, alloc);
  b->reorder = NULL;
}

static int
bitensor_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes)
{
  bitensor_t* b = (bitensor_t*) self;
  // the caller has checked that query_bytes and reply_bytes are positive
  size_t capacity = b->reorder_bytes / reply_bytes;
  if (capacity == 0)
    return convey_error_OFLO;
  if (b->dynamic && !alloc_reorder(b))
    return convey_error_ALLOC;

  b->present = calloc((capacity + 63) / 64, sizeof(uint64_t));
  b->query = malloc(sizeof(packet_t) + query_bytes);
  // FIXME: align this pointer more strongly
  b->reply = malloc(sizeof(packet_t) + reply_bytes);
  if (! (b->present && b->query && b->reply))
    return convey_error_ALLOC;

  b->limit = capacity;
  b->inflight = 0;
  b->await = 0;

  int result = convey_begin(b->biconvey.queries, query_bytes);
  if (result < 0)
    return result;
  return convey_begin(b->biconvey.replies, reply_bytes);
}

static int
bitensor_reset(biconvey_t* self)
{
  bitensor_t* b = (bitensor_t*) self;
  free(b->reply);
  free(b->query);
  free(b->present);
  b->reply = NULL;
  b->query = NULL;
  b->present = NULL;
  if (b->dynamic)
    dealloc_reorder(b);
  return convey_OK;
}

static int
bitensor_free(biconvey_t* self)
{
  bitensor_t* b = (bitensor_t*) self;
  bitensor_reset(&b->biconvey);
  if (!b->dynamic)
    dealloc_reorder(b);
  free(b);
  return convey_OK;
}


/*** Constructor ***/

static const biconvey_methods_t bitensor_methods = {
  .push = &bitensor_push,
  .pull = &bitensor_pull,
  .advance = &bitensor_advance,
  .begin = &bitensor_begin,
  .reset = &bitensor_reset,
  .free = &bitensor_free,
};

biconvey_t*
biconvey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
		    const convey_alc8r_t* alloc, uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  size_t reorder_bytes =
    convey_memory_usage(capacity, false, order, PROCS, n_local, n_buffers);

  bitensor_t* b = malloc(sizeof(bitensor_t));
  if (b == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  *b = (bitensor_t) {
    .biconvey = { ._class_ = &bitensor_methods },
    .reorder_bytes = reorder_bytes,
    .dynamic = (options & convey_opt_DYNAMIC),
    .alloc = *alloc,
  };

  biconvey_t* self = &b->biconvey;
  options |= convey_opt_PROGRESS;
  self->queries = convey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
  self->replies = convey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
  bool ok = (self->queries && self->replies);
  if (ok && !b->dynamic)
    ok = alloc_reorder(b);
  if (!ok) {
    convey_free(self->replies);
    convey_free(self->queries);
    free(b);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  return self;
}
