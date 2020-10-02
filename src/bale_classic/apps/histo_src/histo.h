/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 *****************************************************************/ 

/*! \file histo.h
 * \brief header file for the histogram app.
 */
#ifndef HISTO_H
#define HISTO_H
#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <locale.h>

/*! \struct histo_t
 * \brief A structure to carry all the histogram arrays, counts to different implementations,
 *   aids in error checking
 */ 
typedef struct histo_t {
  SHARED int64_t * counts;  /*!< the shared array that holds the histogram counts */
  int64_t * lcounts;        /*!< the local pointer to the per thread parts of counts */
  int64_t num_counts;       /*!< the global size of the counts array */
  int64_t lnum_counts;      /*!< the local size of the counts array */
  int64_t * index;          /*!< the local index array */
  int64_t * pckindx;        /*!< the packed index with the divmod calculation already done */
  int64_t l_num_ups;        /*!< the local number of update to do */
} histo_t;

double histo_agp(histo_t * data);                       /*!< The AGP implementation */
double histo_exstack(histo_t * data, int64_t buf_cnt);  /*!< The EXSTACK implementation */
double histo_exstack2(histo_t * data, int64_t buf_cnt); /*!< The EXSTACK2 implementation */
double histo_conveyor(histo_t * data);                  /*!< The CONVEYOR implementation */

#include "alternates/histo_alternates.h"

#endif
