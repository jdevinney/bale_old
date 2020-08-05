/** \file biconvey.h
 * The public API for biconveyors.
 */

#ifndef BICONVEY_H
#define BICONVEY_H

#include "convey.h"

typedef struct biconveyor biconvey_t;


/*** CONSTRUCTORS ***/

biconvey_t*
biconvey_new(size_t max_bytes, size_t n_local,
	     const convey_alc8r_t* alloc, uint64_t options);

biconvey_t*
biconvey_new_simple(size_t capacity, const convey_alc8r_t* alloc,
		    const convey_mpp_a2a_t* a2a, uint64_t options);

biconvey_t*
biconvey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
		    const convey_alc8r_t* alloc, uint64_t options);


/*** BICONVEYOR METHODS ***/

int biconvey_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes,
		   void (*answer)(const void* query, void* reply, void* context),
		   void* context);

int biconvey_push(biconvey_t* self, const void* query, int64_t to);

int biconvey_pull(biconvey_t* self, void* reply);

int biconvey_advance(biconvey_t* self, bool done);

int biconvey_reset(biconvey_t* self);

int biconvey_free(biconvey_t* self);


#endif
