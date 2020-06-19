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

/*! \file demo_spmat.c
 * \brief Program that demonstrates (checks) some of the rountines in
 * spmat_util.c
 */

#include "spmat_utils.h"

int main(int argc, char * argv[]) 
{
  int i;

  int64_t *p = rand_perm(10, 0);
  printf("Identity perm = ");
  for(i=0; i<10; i++)
    printf(" %"PRId64,p[i]);
  printf("\n");
  free(p);

  int64_t *q = rand_perm(12, 1);
  printf("Random perm = ");
  for(i=0; i<12; i++)
    printf(" %"PRId64,q[i]);
  printf("\n");

  printf("Is q a perm?  %s\n", (is_perm(q,12))?"yes":"no");
  q[0] = 12;
  printf("Is q a perm?  %s\n", (is_perm(q,12))?"yes":"no");
  free(q);



  sparsemat_t *graph;
  graph = read_matrix_mm("demo.mm");
  spmat_stats(graph);

  dump_matrix(graph, 4, "dump_4.out");
  dump_matrix(graph, 0, "dump_0.out");

  write_matrix_mm(graph, "tidy_demo.mm.out");
  clear_matrix(graph);

  printf("Generate a Kronecker Product of Stars\n");
  kron_args_t * kron_args = kron_args_init("2: 2 2");
  printf("-- input %s\n", kron_args->str);
  printf("-- mode %"PRId64"\n", kron_args->mode);
  printf("-- num_stars %"PRId64"\n", kron_args->num_stars);
  for(int i=0; i<kron_args->num_stars; i++)
    printf("%"PRId64" ", kron_args->star_size[i]); 
  printf("\n-- numrows %"PRId64"\n", kron_args->numrows);
  printf("Known number of triangles = %"PRId64"\n", tri_count_kron_graph(kron_args));

  sparsemat_t *Kron = kronecker_product_graph(kron_args);
  spmat_stats(Kron);
  dump_matrix(Kron, 0, "kron.out");

  clear_matrix(Kron);
  free(kron_args);

}
