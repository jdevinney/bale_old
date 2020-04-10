/******************************************************************
//
//
// Copyright(C) 2018, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// This material may be reproduced by or for the US Government
// pursuant to the copyright license under the clauses at DFARS
// 252.227-7013 and 252.227-7014.
// 
//
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//  * Neither the name of the copyright holder nor the
//   names of its contributors may be used to endorse or promote products
//   derived from this software without specific prior written permission.
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
// 
*****************************************************************/ 

/*! \file sssp.c
 * \brief Demo application that implements the Delta-stepping algorithm for Single
 * Source Shortest Path.
 */

#include <spmat_utils.h>

int main(){
  int main(int argc, char * argv[]) 
{
  #define NUMROWS 10000 
  int64_t numrows=NUMROWS;
  double er_prob = 0.01;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  enum MODEL {GENERIC_Model=1, ALL_Models=2};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;

  sparsemat_t *mat, *tmat;
 
  int64_t dump_files = 0;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows ); break;
    case 's': sscanf(optarg,"%d" ,&seed ); break;
    case 'e': sscanf(optarg,"%lg", &er_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: triangles: read graph from %s Failed\n", filename); exit(1);}
  } else {
    mat = erdos_renyi_tri(numrows, er_prob, ER_TRI_L, seed);
    if(!mat){printf("ERROR: triangles: erdos_renyi_graph Failed\n"); exit(1);}
  }
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of transpose matrix\n");
    fprintf(stderr,"Number of rows    (-n) %ld\n", numrows);
    fprintf(stderr,"random seed     (-s)= %d\n", seed);
    fprintf(stderr,"erdos_renyi_prob   (-e)= %lg\n", er_prob);
    fprintf(stderr,"models_mask     (-M)= %d\n", models_mask);
    fprintf(stderr,"readgraph      (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files      (-D)\n");
    fprintf(stderr,"quiet        (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }

  tmat = transpose_matrix(mat); 
  if(dump_files){
    dump_matrix(mat, 20, "mat.out");
    dump_matrix(tmat, 20, "tmat.out");
  }

  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic      Triangle: ");
      laptime = sssp(mat);
      break;
    default:
      continue;
    }
  }
 
  clear_matrix(mat);
  return(0);
}



}
