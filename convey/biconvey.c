#include "biconvey_impl.h"
#include "private.h"
#define PANIC(ERR) convey_panic(self->queries, __func__, ERR)


int
biconvey_push(biconvey_t* self, const void* query, int64_t to)
{
  return self->_class_->push(self, query, to);
}

int
biconvey_pull(biconvey_t* self, void* reply)
{
  return self->_class_->pull(self, reply);
}

int
biconvey_advance(biconvey_t* self, bool done)
{
  if (self == NULL)
    return convey_error_NULL;
  // FIXME: do more error checking...
  int result = self->_class_->advance(self, done);
  return (result < 0) ? PANIC(result) : result;
}

int
biconvey_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes,
	       void (*answer)(const void* query, void* reply, void* context),
	       void* context)
{
  if (self == NULL)
    return convey_error_NULL;
  if (answer == NULL)
    return PANIC(convey_error_NOFUNC);
  if (query_bytes == 0 || reply_bytes == 0)
    return PANIC(convey_error_ZERO);
  if (self->queries->state != convey_DORMANT ||
      self->replies->state != convey_DORMANT)
    return convey_error_STATE;

  self->answer = answer;
  self->context = context;
  int result = self->_class_->begin(self, query_bytes, reply_bytes);
  return (result < 0) ? PANIC(result) : result;
}

int
biconvey_reset(biconvey_t* self)
{
  if (self == NULL)
    return convey_error_NULL;
  int err = convey_reset(self->queries);
  if (err != convey_OK)
    return err;
  err = convey_reset(self->replies);
  if (err != convey_OK)
    return err;
  err = self->_class_->reset(self);
  return (err < 0) ? PANIC(err) : err;
}

int
biconvey_free(biconvey_t* self)
{
  if (self == NULL)
    return convey_OK;

  int err = convey_free(self->replies);
  if (err != convey_OK)
    return err;
  self->replies = NULL;

  err = convey_free(self->queries);
  if (err != convey_OK)
    return err;
  self->queries = NULL;

  self->answer = NULL;
  err = self->_class_->free(self);
  return (err < 0) ? PANIC(err) : err;
}


/*** Generic Constructors ***/

biconvey_t*
biconvey_new(size_t max_bytes, size_t n_local,
	     const convey_alc8r_t* alloc, uint64_t options)
{
  if (max_bytes < SIZE_MAX)
    max_bytes /= 3;
  if (n_local == 0)
    n_local = convey_procs_per_node();

  size_t capacity, n_buffers;
  int sync, order;
  convey_parameters(max_bytes, n_local, &capacity, &n_buffers, &sync, &order);

  if (sync)
    return biconvey_new_simple(capacity, alloc, NULL, options);
  else
    return biconvey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
}
