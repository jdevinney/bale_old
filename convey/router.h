// Copyright (c) 2019, Institute for Defense Analyses
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


// The next three functions convert the PE number into a tag and the
// destination rank for pushing into porters[0].  A matrix routing tag is
// 16 bits (local).  A tensor routing tag is either 24 bits (remote) +
// 8 bits (local) or 8 bits + 8 bits.

static ROUTER_HINT route_t
vector_route(tensor_t* vector, int64_t pe)
{
  return (route_t) { .tag = 0, .next = pe };
}

static ROUTER_HINT route_t
matrix_route(tensor_t* matrix, int64_t pe)
{
  // dest is (x',y'), we are (x,y); hop to (x',y), tag is (y')
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, matrix->div_local);
  uint32_t lower = dest - matrix->n_local * upper;
#if MATRIX_REMOTE_HOP == 0
  return (route_t) { .tag = lower, .next = upper };
#else
  return (route_t) { .tag = upper, .next = lower };
#endif
}

static ROUTER_HINT route_t
tensor_route(tensor_t* tensor, int64_t pe)
{
  // dest is (x',y',z'), we are (x,y,z)
  // hop to (x,y,y'), tag is (x',z') [24 bits, 8 bits]
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, tensor->div_square);
  uint32_t middle = _divbymul32(dest, tensor->div_local);
  uint32_t lower = dest - tensor->n_local * middle;
  middle -= tensor->n_local * upper;
  uint32_t tag = (upper << 8) | lower;
  return (route_t) { .tag = tag, .next = middle };
}
