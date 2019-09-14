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


#include "example.h"

int
main(int argc, char* argv[])
{
  example_start();

  // Parse command line and environment
  long bins = 10000;
  long load = 100000L;
  uint64_t seed = 1;
  if (argc > 1)
    bins = strtol(argv[1], NULL, 0);
  if (argc > 2)
    load = strtol(argv[2], NULL, 0);
  if (argc > 3)
    seed = strtoull(argv[3], NULL, 0);
  if (MY_PROC == 0)
    printf("command: %s %ld %ld %" PRIu64 "\n"
           "(parameters are: bins, load, seed)\n",
           argv[0], bins, load, seed);

  // Initialize local data
  brand_t _prng;
  brand_init(&_prng, (seed << 32) + MY_PROC);
  long* counts = calloc(bins, sizeof(long));
  long area = PROCS * bins;

  // Build, use, and release the conveyor
  int status = EXIT_FAILURE;
  convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL,
                                  convey_opt_SCATTER | convey_opt_ALERT);
  if (conveyor && counts) {
    convey_begin(conveyor, sizeof(long));

    /*** START OF CONVEYOR LOOP ***/
    long n = 0;
    long index = brand(&_prng) % area;
    while (convey_advance(conveyor, n == load)) {
      for (; n < load; n++) {
        long payload = index / PROCS;
        long pe = index % PROCS;
        if (! convey_push(conveyor, &payload, pe))
          break;
        index = brand(&_prng) % area;
      }

      long* local;
      while (local = convey_apull(conveyor, NULL))
        counts[*local] += 1;
    }
    /*** END OF CONVEYOR LOOP ***/

    convey_reset(conveyor);
    status = EXIT_SUCCESS;

    // Produce a modest amount of output without further communication
    long peak = 0, where = 0;;
    for (long i = 0; i < bins; i++)
      if (counts[i] > peak) {
        peak = counts[i];
        where = i;
      }
    double lambda = load * 1.0 / bins;
    double tail = 0.0;
    for (long i = 0; i < 10; i++)
      tail += exp(-lambda + (peak + i) * log(lambda) - lgamma(peak + i + 1));
    if (tail * area < 4.0) {
      printf("RESULT: %ld[%ld] = %ld\n", (long)(MY_PROC), where, peak);
      fflush(stdout);
    }
  }
  convey_free(conveyor);
  free(counts);

  example_end();
  exit(status);
}
