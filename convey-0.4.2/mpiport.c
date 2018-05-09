// Copyright (c) 2018, Institute for Defense Analyses,
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
//
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.

#include <stdlib.h>
#include "porter_impl.h"
#include "private.h"

// Keep track of porter IDs in order to have unique MPI tags
static uint64_t used_ids = 0;

// For simplicity, we assume 2 buffers per link for now.
typedef struct mpi_porter {
  porter_t porter;
  // Local memory
  uint32_t* recv_buffer;
  // MPI-specific stuff
  MPI_Comm comm;
  MPI_Request recv_req;
  MPI_Request* send_req;
  int porter_id;        // tag for message matching
  int mpi_error;        // keep track of first serious error
  // State and statistics
  int n_active;         // counts open incoming channels
  bool gravid;          // does recv_buffer contain data?
} mpi_porter_t;


/*** Memory Functions ***/

static void*
symmetric_alloc(void* alc8r, size_t size, const char* tag, uint64_t value)
{
  return malloc(size);
}

static void
symmetric_free(void* alc8r, void* ptr)
{
  free(ptr);
}

static const convey_alc8r_t porter_alc8r = {
  .alc8r = NULL, .grab = &symmetric_alloc, .free = &symmetric_free,
};

static bool
porter_grab_buffers(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  convey_alc8r_t* alloc = &self->alloc;
  size_t size = self->buffer_stride * self->n_ranks << self->abundance;
  PARALLEL_ALLOC(self, send_buffers, alloc, size, uint32_t);
  PARALLEL_ALLOC(mpip, recv_buffer, alloc, self->buffer_quads, uint32_t);
  return self->send_buffers && mpip->recv_buffer;
}

static void
porter_free_buffers(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  convey_alc8r_t* alloc = &self->alloc;
  PARALLEL_DEALLOC(mpip, recv_buffer, alloc);
  PARALLEL_DEALLOC(self, send_buffers, alloc);
  mpip->recv_buffer = NULL;
  self->send_buffers = NULL;
}


/*** Internal Methods ***/

static void
mpiport_test_code(mpi_porter_t* mpip, int rc)
{
  if (rc != MPI_SUCCESS && mpip->mpi_error == MPI_SUCCESS) {
    mpip->mpi_error = rc;
    // Report the first error
    char string[MPI_MAX_ERROR_STRING+1];
    int length;
    MPI_Error_string(rc, string, &length);
    string[length] = 0;
    mprint(MY_PROC, 0, "MPI ERROR: %s\n", string);
  }
}

static void
mpiport_try_recv(mpi_porter_t* mpip)
{
  int rc = MPI_Irecv(mpip->recv_buffer, mpip->porter.buffer_quads, MPI_UINT32_T,
                     MPI_ANY_SOURCE, mpip->porter_id, mpip->comm, &mpip->recv_req);
  mpiport_test_code(mpip, rc);
}

static bool
mpiport_ready(porter_t* self, int dest, uint64_t emitted)
{
  return self->channels[dest].delivered >= emitted;
}

static bool
mpiport_send(porter_t* self, int dest, uint64_t level,
             size_t n_bytes, buffer_t* buffer, uint64_t signal)
{
  // Use buffer->start as a completion flag
  buffer->start = signal & 1;
  if (n_bytes == 0)
    n_bytes = sizeof(buffer_t);

  DEBUG_PRINT("sending %zu bytes%s to %d\n",
              n_bytes, (signal & 1) ? "!" : "", dest);
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  int rc = MPI_Isend(buffer, n_bytes >> 2, MPI_UINT32_T, self->relative[dest],
                     mpip->porter_id, mpip->comm, &mpip->send_req[dest]);
  mpiport_test_code(mpip, rc);

  self->send_count++;
  return false;  // delivery is asynchronous
}

static bool
mpiport_progress(porter_t* self, int dest)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;

  if (dest >= 0) {
    // Check whether an outstanding send has completed
    MPI_Request* req = &mpip->send_req[dest];
    if (*req != MPI_REQUEST_NULL) {
      MPI_Status _status;
      int ready = 0;
      int rc = MPI_Test(req, &ready, &_status);
      mpiport_test_code(mpip, rc);
      if (!ready)
        return false;
      porter_record_delivery(self, dest, self->channels[dest].delivered + 1);
    }
    return true;
  }

  // Try to complete all outstanding sends
  int n_ranks = self->n_ranks;
  MPI_Status status[n_ranks];
  int index[n_ranks];
  int outcount;
  int rc = MPI_Testsome(n_ranks, mpip->send_req, &outcount, index, status);
  mpiport_test_code(mpip, rc);
  if (outcount != MPI_UNDEFINED)
    for (int i = 0; i < outcount; i++)
      porter_record_delivery(self, index[i],
                             self->channels[index[i]].delivered + 1);

  // Check whether all sends are complete
  bool done = true;
  for (int i = 0; done && i < n_ranks; i++) {
    channel_t* channel = &self->channels[i];
    done = (channel->delivered == channel->emitted);
  }
  return done;
}


/*** Public Methods ***/

static bool
mpiport_setup(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  bool ok = !self->dynamic || porter_grab_buffers(self);

  if (ok) {
    size_t n = self->n_ranks;
    int n_active = 0;
    for (int i = 0; i < n; i++) {
      n_active += (self->relative[i] >= 0);
      mpip->send_req[i] = MPI_REQUEST_NULL;
    }
    mpip->recv_req = MPI_REQUEST_NULL;

    // Current state
    mpip->mpi_error = MPI_SUCCESS;
    mpip->n_active = n_active;
    mpip->gravid = false;
  }

  mpp_barrier(1);
  // Initiate a nonblocking receive operation
  if (ok && mpip->n_active)
    mpiport_try_recv(mpip);
  return ok;
}

static buffer_t*
mpiport_borrow(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  buffer_t* buffer = (buffer_t*) mpip->recv_buffer;
  if (mpip->gravid)
    return buffer; 
  if (mpip->n_active == 0) {
    self->drained = true;
    return NULL;
  }

  while (1) {
    MPI_Status _status;
    int ready = 0;
    int rc = MPI_Test(&mpip->recv_req, &ready, &_status);
    mpiport_test_code(mpip, rc);
    if (!ready)
      return NULL;

    if (buffer->start == 0)
      return buffer;
    mpip->n_active--;
    buffer->start = 0;
    if (buffer->limit > 0)
      return buffer;
    // this is an empty buffer, so try to get another
    if (mpip->n_active == 0) {
      self->drained = true;
      return NULL;
    }
    mpiport_try_recv(mpip);
  }
}

static void
mpiport_return(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  buffer_t* taken = (buffer_t*) mpip->recv_buffer;
  mpip->gravid = (taken->start < taken->limit);
  if (!mpip->gravid && mpip->n_active)
    mpiport_try_recv(mpip);
}

static void
mpiport_reset(porter_t* self)
{
  if (self->dynamic)
    porter_free_buffers(self);
}

static void
mpiport_demolish(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  if (mpip->porter_id >= 0)
    used_ids &= ~(UINT64_C(1) << mpip->porter_id);
  porter_free_buffers(self);
  free(mpip->send_req);
}

static const porter_methods_t mpi_porter_methods = {
  .setup = &mpiport_setup,
  .borrow = &mpiport_borrow,
  .turnin = &mpiport_return,
  .reset = &mpiport_reset,
  .demolish = &mpiport_demolish,
  .ready = &mpiport_ready,
  .send = &mpiport_send,
  .progress = &mpiport_progress,
};


/*** Constructor and Destructor ***/

porter_t*
porter_new(int n, int32_t relative[n], int my_rank,
           size_t item_size, size_t n_items, size_t multiplicity,
           const convey_alc8r_t* alloc, bool dynamic, bool local,
           bool steady, int opcode)
{
  if (n <= 0 || multiplicity == 0 || (multiplicity & multiplicity-1))
    return NULL;

  mpi_porter_t* mpip = malloc(sizeof(mpi_porter_t));
  if (mpip == NULL)
    return NULL;
  porter_t* porter = &mpip->porter;

  if (alloc == NULL)
    alloc = &porter_alc8r;

  const size_t packet_quads = 1 + ((item_size + 3) >> 2);
  const size_t buffer_quads = (sizeof(buffer_t) >> 2) + n_items * packet_quads;
  const size_t buffer_stride = (buffer_quads + 3) & ~(size_t)3;
  *mpip = (mpi_porter_t) { .porter = {
      ._class_ = &mpi_porter_methods,
      .item_bytes = item_size, .packet_quads = packet_quads,
      .buffer_quads = buffer_quads, .buffer_stride = buffer_stride,
      .n_ranks = n, .my_rank = my_rank, .abundance = 1,  // FIXME
      .relative = relative, .alloc = *alloc, .opcode = opcode,
      .dynamic = dynamic, },
    .porter_id = -1,
  };
  mpp_comm_t comm = MPP_COMM_CURR;
  mpip->comm = comm.internal->mpi_comm;

  // Local allocations
  porter->send_areas = malloc(n * sizeof(area_t));
  porter->all_sent = malloc(n * sizeof(bool));
  porter->channels = malloc(n * sizeof(channel_t));
  if (steady)
    porter->waiting = malloc(n * sizeof(uint8_t));
  mpip->send_req = malloc(n * sizeof(MPI_Request));
  bool ok = (porter->send_areas && porter->all_sent && porter->channels &&
             (!steady || porter->waiting) && mpip->send_req);
  if (!dynamic)
    ok &= porter_grab_buffers(porter);

  ok = mpp_and_long(ok);
  if (!ok) {
    // could report an out-of-memory error here
    porter_demolish(porter);
    return NULL;
  }

  // Agree on an unused ID for this porter
  uint64_t avail = mpp_and_long(~used_ids);
  if (avail == 0) {
    // could report an out-of-ids error here
    porter_demolish(porter);
    return NULL;
  }
  int id = _trailz(avail);
  used_ids |= UINT64_C(1) << id;
  mpip->porter_id = id;

  // Standardize the relative[] array
  int32_t n_procs = PROCS;
  for (int i = 0; i < n; i++)
    if (relative[i] < 0 || relative[i] >= n_procs)
      relative[i] = -1;

  // porter_setup() will erase everything
  return porter;
}
