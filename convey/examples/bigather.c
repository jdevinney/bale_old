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

#include <string.h>
#include "example.h"
#include "biconvey.h"

static void
check_bike(int err, const char* method)
{
  if (err < 0) {
    fprintf(stderr, "biconveyor error in %s: %s\n", method,
	    convey_error_string(NULL, err));
    exit(EXIT_FAILURE);
  }
}

static void
lookup(const void* query, void* reply, void* context)
{
  long* source = context;
  long local;
  memcpy(&local, query, sizeof(long));
  memcpy(reply, &source[local], sizeof(long));
}

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

  biconvey_t* bike = biconvey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);

  if (bike && index && source && target) {
    // Transpose, then perform inverse transpose
    for (int loop = 0; loop < 2; loop++) {
      int err = biconvey_begin(bike, sizeof(long), sizeof(long), &lookup, source);
      check_bike(err, "begin");

      /*** START OF BICONVEYOR LOOP ***/
      long i = 0, j = 0;
      while (biconvey_advance(bike, i == dim)) {
        for (; i < dim; i++) {
	  long query = index[i] / n_procs;
          long pe = index[i] % n_procs;
          if (! biconvey_push(bike, &query, pe))
            break;
	  // fprintf(stderr, "pushed query[%ld] = %ld\n", i, query);
        }
	while (biconvey_pull(bike, &target[j])) {
	  // fprintf(stderr, "got reply[%ld] = %016lx\n", j, target[j]);
	  j++;
	}
      }
      /*** END OF BICONVEYOR LOOP ***/

      err = biconvey_reset(bike);
      check_bike(err, "reset");

      long* temp;
      temp = source, source = target, target = temp;
      temp = index, index = xedni, xedni = temp;
    }

    // Check that the values are back where they started
    bool ok = true;
    for (long i = 0; ok && i < dim; i++)
      if (source[i] != original[i]) {
        printf("ERROR: %ld[%ld] is wrong (%lx should be %lx)\n",
	       my_proc, i, source[i], original[i]);
        ok = false;
      }
    if (ok && my_proc == 0)
      printf("no errors on PE 0\n");
    fflush(stdout);

    if (ok)
      status = EXIT_SUCCESS;
  }

  biconvey_free(bike);
  free(original);
  free(target);
  free(source);
  free(xedni);
  free(index);

  example_end();
  exit(status);
}
