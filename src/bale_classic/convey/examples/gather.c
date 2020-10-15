// Copyright (c) 2020, Institute for Defense Analyses
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

typedef struct {
  long slot;
  long value;
} packet_t;

int
main(int argc, char* argv[])
{
  example_start();

  // Parse command line and environment
  long dim = 10000;
  if (argc > 1)
    dim = strtol(argv[1], NULL, 0);
  if (MY_PROC == 0)
    printf("command: %s %ld\n(parameter is: local_array_length)\n",
           argv[0], dim);

  // Prepare local data
  int status = EXIT_FAILURE;
  long* index = malloc(dim * sizeof(long));
  long* xedni = malloc(dim * sizeof(long));
  long* source = malloc(dim * sizeof(long));
  long* target = malloc(dim * sizeof(long));
  long* original = malloc(dim * sizeof(long));
  long n_procs = PROCS;
  long my_proc = MY_PROC;
  if (index && source && target) {
    brand_t _prng;
    brand_init(&_prng, 1 + my_proc);
    for (long i = 0; i < dim; i++) {
      original[i] = brand(&_prng);
      source[i] = original[i];
      index[i] = dim * my_proc + i;
      long k = n_procs * i + my_proc;
      xedni[i] = n_procs * (k % dim) + (k / dim);
    }
  }

  convey_t* request = convey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT | convey_opt_SCATTER);
  convey_t* reply = convey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);

  if (request && reply && index && source && target) {
    // Transpose, then perform inverse transpose
    for (int loop = 0; loop < 2; loop++) {
      convey_begin(request, sizeof(packet_t), alignof(packet_t));
      convey_begin(reply, sizeof(packet_t), alignof(packet_t));

      /*** START OF CONVEYOR LOOP ***/
      long n = 0;
      bool more;
      while (more = convey_advance(request, n == dim),
             more | convey_advance(reply, !more)) {\
        for (; n < dim; n++) {
          packet_t packet = { .slot = n, .value = index[n] / n_procs };
          long pe = index[n] % n_procs;
          if (! convey_push(request, &packet, pe))
            break;
        }

        packet_t* p;
        int64_t from;
        while ((p = convey_apull(request, &from)) != NULL) {
          packet_t packet = { .slot = p->slot, .value = source[p->value] };
          if (! convey_push(reply, &packet, from)) {
            convey_unpull(request);
            break;
          }
        }

        while ((p = convey_apull(reply, NULL)) != NULL)
          target[p->slot] = p->value;
      }
      /*** END OF CONVEYOR LOOP ***/

      convey_reset(reply);
      convey_reset(request);

      long* temp;
      temp = source, source = target, target = temp;
      temp = index, index = xedni, xedni = temp;
    }

    // Check that the values are back where they started
    bool ok = true;
    for (long i = 0; ok && i < dim; i++)
      if (source[i] != original[i]) {
        printf("ERROR: %ld[%ld] is wrong\n", my_proc, i);
        ok = false;
      }
    if (ok && my_proc == 0)
      printf("no errors on PE 0\n");
    fflush(stdout);

    if (ok)
      status = EXIT_SUCCESS;
  }

  convey_free(reply);
  convey_free(request);
  free(original);
  free(target);
  free(source);
  free(xedni);
  free(index);

  example_end();
  exit(status);
}
