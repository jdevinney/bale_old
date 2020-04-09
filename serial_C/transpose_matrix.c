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

/*! \file transpose_matrix.c
 * \brief Demo program that runs the transpose_matrix from the 
 *  spmat library to provide a framework to study other implementations.
 */

#include "spmat_utils.h"

/*! \page transpose_matrix_page Transpose a given sparse matrix */


double transpose_generic(sparsemat_t *A, int64_t dump_files)
{
  double tm;
  if(dump_files)
    write_matrix_mm(A, "er_orig.mm");

  tm = wall_seconds();
  sparsemat_t *At = transpose_matrix(A);
  tm = wall_seconds() - tm;

  if(dump_files)
    write_matrix_mm(At, "er_tran.mm");
  clear_matrix(At);
  return(tm);
}

int main(int argc, char * argv[])
{
  #define NUMROWS 10000
  int64_t numrows=NUMROWS;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  sparsemat_t *mat;
  int64_t dump_files=0;
  double laptime = 0.0;
  double er_prob = 0.1;
  
  enum MODEL {GENERIC_Model=1, ALL_Models=2};
  uint32_t use_model;
  uint32_t models_mask = ALL_Models - 1;
  int printhelp=0;
  int quiet=0;
  
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:DM:e:f:q")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows );  break;
    case 's': sscanf(optarg,"%d" ,&seed );  break;
    case 'e': er_prob = sscanf(optarg,"%lg", &er_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask);  break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }
  
  if( readgraph ) {
     mat = read_matrix_mm(filename);
     if(!mat){printf("ERROR: transpose_matrix: read_matrix (%s) failed\n", filename); exit(1);}
  } else {
     mat = erdos_renyi_graph(numrows, er_prob, ER_GRAPH_DIRECT, seed);
     if(!mat){printf("ERROR: transpose_matrix: erdos_renyi_graph failed\n"); exit(1);}
  }
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of transpose matrix\n");
    fprintf(stderr,"Number of rows       (-n) %ld\n", numrows);
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"erdos_renyi_prob     (-e)= %lg\n", er_prob);
    fprintf(stderr,"readgraph            (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files           (-D)\n");
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
  
    if(printhelp)
      return(0);
  }
  if(!quiet) {
    printf("Input matrix stats:\n");
    spmat_stats(mat);
  }
  if(dump_files)
    dump_matrix(mat, 20, "dump.out");
  
  for( use_model=1; use_model < 2; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC_Model:
    if(!quiet) printf("transpose matrix : ");
    laptime = transpose_generic(mat, dump_files);
    break;
    default:
    continue;
    }
    if(!quiet) printf("  %8.3lf seconds \n", laptime);
  }
  return(0);
}

