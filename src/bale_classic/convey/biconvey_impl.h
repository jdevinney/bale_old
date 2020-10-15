#ifndef BICONVEY_IMPL_H
#define BICONVEY_IMPL_H

#include "biconvey.h"
#include "convey_impl.h"

typedef struct biconveyor_methods {
  int (*push)(biconvey_t* self, const void* query, int64_t to);
  int (*pull)(biconvey_t* self, void* reply);
  int (*advance)(biconvey_t* self, bool done);
  int (*begin)(biconvey_t* self, size_t query_bytes, size_t reply_bytes);
  int (*reset)(biconvey_t* self);
  int (*free)(biconvey_t* self);
} biconvey_methods_t;

struct biconveyor {
  const biconvey_methods_t* _class_;
  size_t query_bytes;
  size_t reply_bytes;
  void (*answer)(const void*, void*, void*);
  void* context;
  uint64_t suppress;
};

#endif
