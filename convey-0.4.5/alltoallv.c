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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "common.h"
#include "convey.h"
#include "bolite.h"


static void
conveyor_bug(convey_t* c, const char* call, int error)
{
  mfprint(stderr, MY_PROC, 0, "%s failed: %s (%d)\n", call,
          convey_error_string(c, error), error);
  mpp_exit(-error);
}


/*** One Conveyor, Fixed Size ***/

typedef struct checksum {
  uint64_t sent;
  uint64_t rcvd;
} checksum_t;

static checksum_t
communicate(convey_t* conveyor, brand_t* prng, double load, size_t size,
            convey_t* tally, double reject)
{
  checksum_t sums = { 0, 0 };
  size_t n_words = (size + 7) >> 3;
  uint64_t send[n_words], recv[n_words], valid[n_words];
  memset(send, 0, 8 * n_words);
  memset(recv, 0, 8 * n_words);
  memset(valid, 0, 8 * n_words);
  memset(valid, 0xFF, size);

  const size_t n_procs = PROCS;
  const size_t mask = _maskr(64 - _leadz(n_procs - 1));
  const size_t bulk = ceil(load * n_procs);
  const uint64_t self = MY_PROC;
  // make sure data generation is deterministic
  size_t i, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_begin", rv);

  if (tally == NULL)
    while (1) {
      rv = convey_advance(conveyor, sent == bulk);
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      for (; sent < bulk; sent++) {
        if (!pending) {
          for (i = 0; i < n_words; i++)
            send[i] = brand(prng);
          do pe = brand(prng) & mask;
          while (pe >= n_procs);
        }

        rv = convey_push(conveyor, send, pe);
        pending = (rv == convey_FAIL);
        if (pending)
          break;
        if (rv != convey_OK)
          conveyor_bug(conveyor, "convey_push", rv);

        for (i = 0; i < n_words; i++)
          sums.sent += (2*i + 2*self + 1) * (valid[i] & send[i]);
      }

      int64_t from;
      while ((rv = convey_pull(conveyor, recv, &from)) == convey_OK)
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * (valid[i] & recv[i]);
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_pull", rv);
    }

  else {
    // To test a steady conveyor, we can yoke it with a "tally" conveyor to
    // tell each PE how many incoming messages to expect.  We wait for all
    // messages to arrive before starting the main conveyor's endgame.
    size_t rcvd = 0, goal = 0;
    bool tally_done = false;
    bool stuck = false;  // tally conveyor has a pending push?
    uint32_t one = 1;

    while (1) {
      // Advance the main conveyor
      rv = convey_advance(conveyor, tally_done && (rcvd == goal));
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      // Advance the tally conveyor unless it is already done
      if (!tally_done) {
        rv = convey_advance(tally, (sent == bulk) && !stuck);
        tally_done = (rv == convey_DONE);
        if (rv < 0)
          conveyor_bug(tally, "convey_advance", rv);
      }

      // Push loop.  Stop if we fail to push into either conveyor.
      for (; sent < bulk; sent++) {
        if (!stuck) {
          if (!pending) {
            for (i = 0; i < n_words; i++)
              send[i] = brand(prng);
            do pe = brand(prng) & mask;
            while (pe >= n_procs);
          }

          rv = convey_push(conveyor, send, pe);
          pending = (rv == convey_FAIL);
          if (pending)
            break;
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_push", rv);

          for (i = 0; i < n_words; i++)
            sums.sent += (2*i + 2*self + 1) * (valid[i] & send[i]);
        }

        rv = convey_push(tally, &one, pe);
        if (rv < 0)
          conveyor_bug(tally, "convey_push", rv);
        stuck = (rv == convey_FAIL);
        if (stuck)
          break;
      }

      // Pull from the main conveyor and count what we get
      int64_t from;
      while ((rv = convey_pull(conveyor, recv, &from)) == convey_OK) {
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * (valid[i] & recv[i]);
        rcvd++;
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_pull", rv);

      // Pull the tallies also
      if (!tally_done)
        while (convey_apull(tally, NULL))
          goal++;
    }
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
           sent, bulk);
    mpp_exit(1);
  }

  return sums;
}


/*** One Conveyor, Varying Sizes ***/

static checksum_t
elasticomm(convey_t* conveyor, brand_t* prng, double load, size_t max_size,
           convey_t* tally, double reject)
{
  checksum_t sums = { 0, 0 };
  size_t max_words = (max_size + 7) >> 3;
  uint64_t* send = malloc(max_words * sizeof(uint64_t));
  uint64_t* recv = malloc(max_words * sizeof(uint64_t));

  uint64_t valid[8];
  memset(valid, 0, 8 * sizeof(uint64_t));
  for (int i = 0; i < 8; i++)
    memset(&valid[i], 0xFF, (i ? i : 8));

  const size_t n_procs = PROCS;
  const size_t mask = _maskr(64 - _leadz(n_procs - 1));
  const size_t bulk = ceil(load * n_procs);
  const uint64_t self = MY_PROC;

  size_t i, n_bytes = 0, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_begin", rv);

  if (tally == NULL)
    while (1) {
      rv = convey_advance(conveyor, sent == bulk);
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      for (; sent < bulk; sent++) {
        if (!pending) {
          // For now, use a uniform distribution
          n_bytes = 1 + (size_t)(max_size * dbrand(prng));
          size_t n_words = (n_bytes + 7) >> 3;
          for (i = 0; i < n_words; i++)
            send[i] = brand(prng);
          do pe = brand(prng) & mask;
          while (pe >= n_procs);
        }

        rv = convey_epush(conveyor, n_bytes, send, pe);
        pending = (rv == convey_FAIL);
        if (pending)
          break;
        if (rv != convey_OK)
          conveyor_bug(conveyor, "convey_epush", rv);

        size_t n_words = (n_bytes + 7) >> 3;
        for (i = 0; i < n_words - 1; i++)
          sums.sent += (2*i + 2*self + 1) * send[i];
        sums.sent += (2*i + 2*self + 1) * (send[i] & valid[n_bytes & 7]);
      }

      convey_item_t _item;
      while ((rv = convey_epull(conveyor, &_item)) == convey_OK) {
        if (reject > 0.0 && dbrand(prng) < reject) {
          rv = convey_unpull(conveyor);
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_unpull", rv);
          rv = convey_FAIL;
          break;
        }
        size_t n_words = (_item.bytes + 7) >> 3;
        recv[n_words - 1] = 0;
        memcpy(recv, _item.data, _item.bytes);
        uint64_t from = _item.from;
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * recv[i];
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_epull", rv);
    }

  else {
    // See communicate() for an explanation of this code.
    size_t rcvd = 0, goal = 0;
    bool tally_done = false;
    bool stuck = false;
    uint32_t one = 1;

    while (1) {
      rv = convey_advance(conveyor, tally_done && (rcvd == goal));
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      if (!tally_done) {
        rv = convey_advance(tally, (sent == bulk) && !stuck);
        tally_done = (rv == convey_DONE);
        if (rv < 0)
          conveyor_bug(tally, "convey_advance", rv);
      }

      for (; sent < bulk; sent++) {
        if (!stuck) {
          if (!pending) {
            n_bytes = 1 + (size_t)(max_size * dbrand(prng));
            size_t n_words = (n_bytes + 7) >> 3;
            for (i = 0; i < n_words; i++)
              send[i] = brand(prng);
            do pe = brand(prng) & mask;
            while (pe >= n_procs);
          }

          rv = convey_epush(conveyor, n_bytes, send, pe);
          pending = (rv == convey_FAIL);
          if (pending)
            break;
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_epush", rv);

          size_t n_words = (n_bytes + 7) >> 3;
          for (i = 0; i < n_words - 1; i++)
            sums.sent += (2*i + 2*self + 1) * send[i];
          sums.sent += (2*i + 2*self + 1) * (send[i] & valid[n_bytes & 7]);
        }

        rv = convey_push(tally, &one, pe);
        if (rv < 0)
          conveyor_bug(tally, "convey_push", rv);
        stuck = (rv == convey_FAIL);
        if (stuck)
          break;
      }

      convey_item_t _item;
      while ((rv = convey_epull(conveyor, &_item)) == convey_OK) {
        if (reject > 0.0 && dbrand(prng) < reject) {
          rv = convey_unpull(conveyor);
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_unpull", rv);
          rv = convey_FAIL;
          break;
        }
        size_t n_words = (_item.bytes + 7) >> 3;
        recv[n_words - 1] = 0;
        memcpy(recv, _item.data, _item.bytes);
        uint64_t from = _item.from;
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * recv[i];
        rcvd++;
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_epull", rv);

      if (!tally_done)
        while (convey_apull(tally, NULL))
          goal++;
    }
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
           sent, bulk);
    mpp_exit(1);
  }

  free(recv);
  free(send);
  return sums;
}


/*** Two Conveyors ***/

// Idea: table[pos] = (c * pos) % p
// We sum up (i+1) * table[pos[i]] modulo p, and should get
// c * (sum of (i+1) * pos[i]) modulo p.

#define PRIME UINT64_C(3690948119)
#define MULTIPLIER UINT64_C(2646891621)
#define ENTRIES 100000

// Can't cast from local pointer to shared, so save the shared pointer :(
#if MPP_USE_UPC
static shared uint64_t* all_table;
#endif

static uint64_t*
global_table_init(void)
{
  const size_t n = ENTRIES;
#if MPP_USE_UPC
  all_table = upc_all_alloc(n * THREADS, sizeof(uint64_t));
  uint64_t* table = (uint64_t*) &all_table[MYTHREAD];
#else
  uint64_t* table = mpp_alloc(n * sizeof(uint64_t));
#endif

  for (size_t j = 0; j < n; j++) {
    size_t i = (MY_PROC + j * PROCS) % PRIME;
    table[j] = (i * MULTIPLIER) % PRIME;
  }

  mpp_barrier(1);
  return table;
}

static void
global_table_free(uint64_t* table)
{
#if MPP_USE_UPC
  upcx_all_free(all_table);
#else
  mpp_free(table);
#endif
}

static checksum_t
indexgather(convey_t* request, convey_t* reply, brand_t* prng, double load,
            const uint64_t* source)
{
  typedef struct { int64_t slot; int64_t value; } packet_t;
  assert(sizeof(packet_t) == 16);
  int rv;

  rv = convey_begin(request);
  if (rv != convey_OK)
    conveyor_bug(request, "convey_begin", rv);
  rv = convey_begin(reply);
  if (rv != convey_OK)
    conveyor_bug(reply, "convey_begin", rv);

  const size_t n_procs = PROCS;
  checksum_t sums = { 0, 0 };
  bool pending = false;
  uint64_t index = 0;

  size_t sent = 0;
  size_t bulk = ceil(load * n_procs);
  bool more;
  while (1) {
    rv = convey_advance(request, sent == bulk);
    if (rv < 0)
      conveyor_bug(request, "convey_advance", rv);
    more = (rv != convey_DONE);
    rv = convey_advance(reply, !more);
    if (rv < 0)
      conveyor_bug(reply, "convey_advance", rv);
    if (!more && rv == convey_DONE)
      break;

    for (; sent < bulk; sent++) {
      if (!pending)
        index = dbrand(prng) * (ENTRIES * n_procs);
      uint64_t local = index / n_procs;
      packet_t packet = { .slot = sent, .value = local };
      int64_t pe = index - local * n_procs;

      rv = convey_push(request, &packet, pe);
      pending = (rv == convey_FAIL);
      if (pending)
        break;
      if (rv != convey_OK)
        conveyor_bug(request, "convey_push", rv);

      sums.sent += (sent + 1) * index;
      if (sums.sent >> 63)
        sums.sent %= PRIME;
    }

    packet_t packet;
    int64_t from;
    while ((rv = convey_pull(request, &packet, &from)) == convey_OK) {
      packet.value = source[packet.value];
      rv = convey_push(reply, &packet, from);
      if (rv == convey_FAIL) {
        int rv1 = convey_unpull(request);
        if (rv1 != convey_OK)
          conveyor_bug(request, "convey_unpull", rv1);
        break;
      }
      if (rv != convey_OK)
        conveyor_bug(reply, "convey_push", rv);
    }
    if (rv != convey_FAIL)
      conveyor_bug(request, "convey_pull", rv);

    while ((rv = convey_pull(reply, &packet, NULL)) == convey_OK) {
      sums.rcvd += (packet.slot + 1) * packet.value;
      if (sums.rcvd >> 63)
        sums.rcvd %= PRIME;
    }
    if (rv != convey_FAIL)
      conveyor_bug(reply, "convey_pull", rv);
  }

  rv = convey_reset(reply);
  if (rv != convey_OK)
    conveyor_bug(reply, "convey_reset", rv);
  rv = convey_reset(request);
  if (rv != convey_OK)
    conveyor_bug(request, "convey_reset", rv);

  // Convert the accumulated indices into the correct value
  sums.sent %= PRIME;
  sums.sent = (sums.sent * MULTIPLIER) % PRIME;
  sums.rcvd %= PRIME;
  return sums;
}


/*** Wrapper Function ***/

static checksum_t
basictest(convey_t* conveyor, brand_t* prng, double load, size_t size,
          convey_t* tally, bool elastic, double reject)
{
  if (tally) {
    int rv = convey_begin(tally);
    if (rv != convey_OK)
      conveyor_bug(tally, "convey_begin", rv);
  }

  checksum_t sums = (elastic ? &elasticomm : &communicate)
    (conveyor, prng, load, size, tally, reject);

  if (tally) {
    int rv = convey_reset(tally);
    if (rv != convey_OK)
      conveyor_bug(tally, "convey_reset", rv);
  }

  return sums;
}


/*** Main Program ***/

static double now(void)
{
  struct timeval _t;
  gettimeofday(&_t, NULL);
  return _t.tv_sec + 1.0e-6 * _t.tv_usec;
}

static void
usage(const char* prog)
{
  mprint(0, 0,
         "usage: %s [options] flavor load [item_size]\n"
         "flavor is from { simple, twohop, vector, matrix, tensor }\n"
         "options are:\n"
         "  -a #  autotune mpp_util all-to-all for # seconds { simple }\n"
         "  -b #  use # buffers per link { vector, matrix, tensor }\n"
         "  -c #  set buffer capacity to # items (or bytes, if -e)\n"
         "  -d    allocate buffers at each iteration\n"
         "  -e    send items of varying size { vector, matrix, tensor }\n"
         "  -f #  return # if the constructor fails\n"
         "  -g    use 'indexgather' benchmark\n"
         "  -h    report host name for each process\n"
         "  -n #  do # timed iterations, default 1\n"
         "  -p    use conveyor that guarantees steady progress\n"
         "  -r #  reject (unpull) a fraction # of pulled items\n"
         "  -q    suppress mpp_util timings\n"
         "  -s #  set pseudorandom seed to #\n"
         "  -t #  use # procs per local group { twohop, matrix, tensor }\n"
         "  -u    use mpp_util all-to-all functions { simple }\n"
         "  -v    increase verbosity\n"
         "  -w #  do # warmup iterations, default 1\n"
         "  -x    turn on prefetching for scatters\n"
         "  -y    test conveyor for steadiness\n",
         prog);
  mpp_exit(2);
}

void
print_args(int argc, char* argv[])
{
  mprint(0, 0, "args = ");
  for (int i = 1; i < argc; i++)
    mprint_np(0, 0, "%s%c", argv[i], ((i == argc-1) ? '\n' : ' '));
}

int
alltoallv(int argc, char* argv[])
{
  int autotune = 0;
  bool vendor = true;
  bool dynamic = false;
  bool elastic = false;
  bool gather = false;
  bool hosts = false;
  bool steady = false;
  bool quiet = false;
  bool sorter = false;
  bool yoke = false;
  int failcode = 1;
  size_t capacity = 0;
  size_t n_buffers = 1;
  int64_t n_timed = 1;
  int64_t n_warmup = 1;
  double reject = 0.0;
  uint64_t seed = 0;
  size_t row_procs = 0;
  int verbosity = 0;
  const char* flavor;
  double load;
  size_t item_size;

  print_args(argc, argv);

  for (int opt; (opt = getopt(argc, argv, "a:b:c:def:ghn:pr:qs:t:uvw:xy")) != -1; ) {
    // mprint(0, 0, "option character is '%c'\n", opt);
    switch (opt) {
    case 'a': autotune = atoi(optarg); break;
    case 'b': n_buffers = strtoul(optarg, NULL, 10); break;
    case 'c': capacity = strtoul(optarg, NULL, 10); break;
    case 'd': dynamic = true; break;
    case 'e': elastic = true; break;
    case 'f': failcode = atoi(optarg); break;
    case 'g': gather = true; break;
    case 'h': hosts = true; break;
    case 'n': n_timed = strtoul(optarg, NULL, 10); break;
    case 'p': steady = true; break;
    case 'r': reject = atof(optarg); break;
    case 'q': quiet = true; break;
    case 's': seed = strtoul(optarg, NULL, 10); break;
    case 't': row_procs = strtol(optarg, NULL, 10); break;
    case 'u': vendor = false; break;
    case 'v': verbosity++; break;
    case 'w': n_warmup = strtoul(optarg, NULL, 10); break;
    case 'x': sorter = true; break;
    case 'y': yoke = true; break;
    default: usage(argv[0]);
    }
  }
  if (!seed)
    seed = time(NULL);
  if (! (argc == optind + 2 || (!gather && argc == optind + 3)))
    usage(argv[0]);
  if (yoke && gather) {
    mprint(0, 0, "%s: -y is not compatible with -g\n", argv[0]);
    mpp_exit(2);
  }

  flavor = argv[optind];
  load = atof(argv[optind+1]);
  char* slash = strchr(argv[optind+1], '/');
  if (slash)
    load /= atof(slash + 1);
  item_size = (argc == optind+2) ? 16 : strtoul(argv[optind+2], NULL, 10);
  if (!finite(load) || load <= 0.0 || item_size == 0) {
    mprint(0, 0, "%s: bogus parameters %f, %zu\n", argv[0], load, item_size);
    mpp_exit(2);
  }

  if (row_procs == 0) {
    extern size_t convey_procs_per_node(void);
    row_procs = convey_procs_per_node();
    mprint(0, 0, "procs per group = %zu\n", row_procs);
  }
  if (hosts) {
    char name[64];
    gethostname(name, 64);
    name[63] = 0;
    mprint(MY_PROC, 0, "host = %s\n", name);
  }

  size_t n_procs = PROCS;
  convey_t* conveyor = NULL;
  convey_t* echo = NULL;
  mpp_alltoall_t* a2a = NULL;
  const char* constructor = NULL;
  const uint64_t sortopt = sorter * convey_opt_SCATTER;
  const uint64_t options = convey_opt_RECKLESS | dynamic * convey_opt_DYNAMIC |
    steady * convey_opt_PROGRESS | gather * CONVEY_OPT_ALIGN(8);
  const size_t echo_size = yoke ? sizeof(uint32_t) : item_size;
  const uint64_t echo_opts = options - steady * convey_opt_PROGRESS;

  int order = 0;
  if (strcmp(flavor, "vector") == 0)
    order = 1;
  else if (strcmp(flavor, "matrix") == 0)
    order = 2;
  else if (strcmp(flavor, "tensor") == 0)
    order = 3;

  if (elastic && order) {
    if (!capacity)
      capacity = 9999;  // minimum bytes per buffer
    size_t atom_size = (gather ? 16 : 1);

    constructor = ((const char* []) { "convey_new_etensor(1)",
          "convey_new_etensor(2)", "convey_new_etensor(3)" }) [order - 1];
    conveyor = convey_new_etensor(atom_size, capacity, item_size, order,
                                  row_procs, n_buffers, NULL, options);
    if (gather || yoke)
      echo = convey_new_etensor(atom_size, capacity, echo_size, order,
                                row_procs, n_buffers, NULL, echo_opts);
  }

  else if (order) {
    if (!capacity && order == 1)
      capacity = ceil(load + 4.0 * sqrt(load + 1.0));
    else if (!capacity)
      capacity = 1 + 9999 / (item_size + 4);
    constructor = ((const char* []) { "convey_new_tensor(1)",
          "convey_new_tensor(2)", "convey_new_tensor(3)" }) [order - 1];
    conveyor = convey_new_tensor(capacity, item_size, order,
                                 row_procs, n_buffers, NULL, options);
    if (gather || yoke)
      echo = convey_new_tensor(capacity, echo_size, order,
                               row_procs, n_buffers, NULL, echo_opts);
  }

  else if (strcmp(flavor, "twohop") == 0) {
    if (!capacity && row_procs) {
      size_t col_procs = n_procs / row_procs;
      size_t n_min = (row_procs < col_procs) ? row_procs : col_procs;
      capacity = ceil(load * n_min + 5.0 * sqrt(load * n_min + 1.0));
    }
    constructor = "convey_new_twohop";
    conveyor = convey_new_twohop(capacity, item_size, row_procs, NULL, options);
    if (gather || yoke)
      echo = convey_new_twohop(capacity, echo_size, row_procs, NULL, echo_opts);
  }

  else if (strcmp(flavor, "simple") == 0) {
    if (!capacity)
      capacity = ceil(load + 4.0 * sqrt(load + 1.0));

#if HAVE_MPP_UTIL
    if (!vendor || autotune) {
      a2a = mpp_alltoall_create(1);
      if (!autotune)
        mpp_alltoall_init(a2a, 0, 0, 0, 0, 0);
      else {
        size_t bin_bytes = capacity * item_size;
        size_t all_bytes = n_procs * bin_bytes;
        void* src = mpp_alloc(all_bytes);
        void* dst = mpp_alloc(all_bytes);
        if (!src || !dst) {
          mprint(0, 0, "mpp_alloc failed (2 x %zu bytes)\n", all_bytes);
          mpp_exit(1);
        }
        mpp_alltoall_tune(a2a, dst, src, bin_bytes, 0, 0, autotune);
        mpp_free(dst);
        mpp_free(src);
      }
    }
#endif

    constructor = "convey_new_simple";
    conveyor = convey_new_simple(capacity, item_size, NULL, a2a, options | sortopt);
    // Don't use a sorter for the reply conveyor
    if (gather || yoke)
      echo = convey_new_simple(capacity, echo_size, NULL, a2a, echo_opts);
  }

  else {
    mprint(0, 0, "%s: unknown conveyor flavor '%s'\n", argv[0], flavor);
    mpp_exit(2);
  }

  mprint(0, 0, "%s buffer capacity: %zu\n", flavor, capacity);

  // Check the conveyor
  long built = conveyor && (!(gather | yoke) || echo);
  built = mpp_and_long(built);
  if (!built) {
    mprint(0, 0, "%s failed\n", constructor);
    exit(failcode);
  }

  brand_t _prng;
  brand_init(&_prng, seed ^ (uint64_t)MY_PROC << 32);
  uint64_t* table = gather ? global_table_init() : NULL;

  // Test during warmup
  for (int64_t loop = 0; loop < n_warmup; loop++) {
    checksum_t sums;
    if (gather) {
      sums = indexgather(conveyor, echo, &_prng, load, table);
      if (sums.sent != sums.rcvd)
        mprint(MY_PROC, 1, "local failure: expected %08" PRIx64 " got %08" PRIx64 "\n",
               sums.sent, sums.rcvd);
    }
    else
      sums = basictest(conveyor, &_prng, load, item_size, echo, elastic, reject);
    sums.sent = mpp_accum_long(sums.sent);
    sums.rcvd = mpp_accum_long(sums.rcvd);
    if (sums.sent != sums.rcvd)
      mprint(0, 0, "checksum failure: sent %016" PRIx64 " rcvd %016" PRIx64 "\n",
             sums.sent, sums.rcvd);
  }

  if (n_warmup > 0)
    mprint(0, 0, "finished %" PRId64 " warm-up iterations\n", n_warmup);
  mpp_barrier(1);

  // Now do benchmarking
  double bw_sum = 0.0, bw_sqr = 0.0, bw_max = 0.0, seconds = 0.0;
  double start = now();
  for (int64_t loop = 0; loop < n_timed; loop++) {
    if (gather)
      indexgather(conveyor, echo, &_prng, load, table);
    else
      basictest(conveyor, &_prng, load, item_size, echo, elastic, reject);
    mpp_barrier(1);
    // Compute some timing data
    double finish = now();
    double bw = (n_procs * load) / (finish - start);
    bw *= elastic ? 0.5 * (item_size + 1) : 1.0 * item_size;
    bw_sum += bw;
    bw_sqr += bw * bw;
    bw_max = fmax(bw_max, bw);
    seconds += finish - start;
    start = finish;
  }

  if (n_timed > 0) {
    mprint(0, 0, "elapsed time: %g seconds\n", seconds);
    if (n_timed == 1)
      mprint(0, 0, "bandwidth per PE: %g bytes/sec\n", bw_sum / n_timed);
    else {
      double mean = bw_sum / n_timed;
      double stdv = sqrt((bw_sqr - mean * mean * n_timed) / (n_timed - 1));
      mprint(0, 0, "bandwidth per PE: %g +/- %.1f%% bytes/sec (max %g)\n",
             mean, 100.0 * stdv / mean, bw_max);
    }
  }

  long n_comms = convey_statistic(conveyor, convey_COMMS + convey_CUMULATIVE);
  long n_syncs = convey_statistic(conveyor, convey_SYNCS + convey_CUMULATIVE);
  if (gather) {
    n_comms += convey_statistic(echo, convey_COMMS + convey_CUMULATIVE);
    n_syncs += convey_statistic(echo, convey_SYNCS + convey_CUMULATIVE);
  }
  n_comms = mpp_accum_long(n_comms);
  n_syncs = mpp_accum_long(n_syncs);
  double terms = n_procs * (n_warmup + n_timed);
  mprint(0, 0, "average of %g sends and %g syncs per PE\n",
         n_comms / terms, n_syncs / terms);
  if (reject > 0.0) {
    long n_unpulls = convey_statistic(conveyor, convey_UNPULLS + convey_CUMULATIVE);
    n_unpulls = mpp_accum_long(n_unpulls);
    mprint(0, 0, "average of %g unpulls per PE\n", n_unpulls / terms);
  }

  if (gather)
    global_table_free(table);
  if (echo)
    convey_free(echo);
  convey_free(conveyor);

  return quiet;
}

int
main(int argc, char* argv[])
{
  // Initialize mpp_util without threading support
  argc = mpp_util_init(argc, argv, NULL);

  bool quiet = true;
  if (argc == 2 && strcmp(argv[1], "--") == 0) {
    mprint(0, 0, "total PEs = %ld\n", PROCS);

    static char line[256];
    char* tokens[130];
    tokens[0] = argv[0];

    // Read, broadcast, and parse command lines from standard input,
    // until reaching an empty line or end of input.
    while (1) {
      if (MY_PROC == 0) {
        char* temp = fgets(line, 256, stdin);
        if (!temp)
          line[0] = 0;
      }
      mpp_broadcast64(line, 32, 0);

      int n = 1;
      char* saveptr = NULL;
      while (tokens[n] = strtok_r((n == 1) ? line : NULL, " \t\n", &saveptr))
        n++;
      if (n == 1)
        break;
      // optind = 1;
      optind = 0;  // correct value for glibc
      alltoallv(n, tokens);
      fflush(stdout);
    }

    mprint(0, 0, "job finished\n");
  }
  else
    quiet = alltoallv(argc, argv);

#if HAVE_MPP_UTIL
  if (quiet)
    mprint_set_filter(NULL);
#endif
  mpp_util_end();
}
