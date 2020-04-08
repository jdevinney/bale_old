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


/** \file convey_alc8r.h
 * The auxiliary API for symmetric memory allocation.
 */


#ifndef CONVEY_ALC8R_H
#define CONVEY_ALC8R_H

#include <stddef.h>
#include <stdint.h>

/** Memory allocation objects that can be passed to conveyor constructors.
 * 
 * Many conveyor constructors can be told how to allocate and deallocate
 * symmetric memory by passing them a pointer to a structure of this type.
 * The structure does not have to stay in scope; the constructor copies it.
 * If the pointer is \c NULL, the conveyor uses standard memory management
 * functions instead.
 */
typedef struct convey_alc8r {
  /// Passed as the first argument of the following functions.
  void* alc8r;
  /// The function for allocating memory: the \a tag and \a value will be the file name and line number of the call.
  void* (*grab)(void* alc8r, size_t size, const char* tag, uint64_t value);
  /// The function for releasing previously allocated memory.
  void  (*free)(void* alc8r, void* ptr);
} convey_alc8r_t;

#endif
