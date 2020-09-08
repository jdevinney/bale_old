/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
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

/*! \file triangle.h
 * \brief Demo application that counts the number of triangles
 *  in a graph. The graph is stored as a lower triangular sparse matrix.
 */
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <spmat.h>
#include <locale.h>

double   triangle_agi(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);
double   triangle_exstack_push(int64_t *count, int64_t *sr, sparsemat_t *L, sparsemat_t * U, int64_t alg, int64_t bufsiz);
//double   triangle_exstack_pull(int64_t *count, int64_t *sr, sparsemat_t *L, int64_t alg, int64_t bufsiz);
double   triangle_exstack2_push(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg, int64_t bufsiz);
double   triangle_convey_push(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);
//double   triangle_convey_pull(int64_t *count, int64_t *sr, sparsemat_t *mat);

// alternates go here
double   triangle_agi_opt1(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);
double   triangle_agi_opt2(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);
double   triangle_agi_oo(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);
double   triangle_agi_iter(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg);

#endif
