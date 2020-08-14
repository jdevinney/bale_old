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

/*! \file triangle.c
 * \brief Demo application that counts the number of triangles in the graph
 * given by it's adjacency matrix.
 */

#include "spmat_utils.h"

/*! \page triangle_page 
  Count the triangles in a graph presented by a lower triangular matrix 
*/

/*!
 * \brief This routine counts the number of triangles in a graph
 *   given the lower triangular piece of the adjacency matrix
 * \param *triangles a place to write the number of triangles found
 * \param *mat the sparse matrix that holds the graph
 * \return run time
 */
double triangles_matrix(int64_t *triangles, sparsemat_t *mat) 
{
  int64_t j, k, l, numtriangles;
  int64_t U,V;
  numtriangles = 0;

  double t1 = wall_seconds();
  
  // for each non-zero (i,j) in L accumulate the size of the intersection
  // of row_i and row_j.

  for(U = 0; U < mat->numrows; U++){ 
    for(j = mat->offset[U]; j < mat->offset[U+1]; j++){
      V = mat->nonzero[j];
      for( l = mat->offset[U], k = mat->offset[V];  k < mat->offset[V+1] && l < mat->offset[U+1];  ){  // This requires that the matrix be tidy
        if( mat->nonzero[k] == mat->nonzero[l] ){
          numtriangles++;
          k++;
          l++;
        }else if( mat->nonzero[k] > mat->nonzero[l] ){
          l++;
        }else{ // ( mat->nonzero[U] > mat->nonzero[W] ) {
          k++;
        }
      }
    }
  }

  t1 = wall_seconds() - t1;
 
  *triangles = numtriangles;
  return(t1);
}


int main(int argc, char * argv[]) 
{
  #define NUMROWS 10000 
  int64_t numrows=NUMROWS;
  double edge_prob = 0.01;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256] = {"filename"};
  enum FLAVOR {GENERIC=1, ALL=2};
  graph_model model = FLAT;
  char * kron_list_def = {"M: k1 k2 k3 ..."};
  char * kron_list = kron_list_def;
  kron_args_t * kron_args = NULL;
  uint32_t use_model;
  uint32_t models_mask=ALL - 1;
  int printhelp = 0;
  int quiet = 0;

  sparsemat_t *mat;
  int64_t triangles;
 
  int64_t dump_files = 0;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:g:k:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%"SCNd64, &numrows ); break;
    case 's': sscanf(optarg,"%d" , &seed ); break;
    case 'e': model = FLAT;      sscanf(optarg,"%lg", &edge_prob); break;
    case 'g': model = GEOMETRIC; sscanf(optarg,"%lg", &edge_prob); break;
    case 'k': model = KRONECKER; kron_list = optarg; break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg,"%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: triangles: read graph from %s Failed\n", filename); exit(1);}
  } else {
    if(model == KRONECKER){
      kron_args = kron_args_init(kron_list);
      mat = kronecker_product_graph(kron_args);
      if(!mat){printf("ERROR: triangles: Kronecker Product generation Failed\n"); exit(1);}
    } else {
      mat = random_graph(numrows, model, UNDIRECTED, NOLOOPS, edge_prob, seed);
      if(!mat){printf("ERROR: triangles: erdos_renyi_graph generation Failed\n"); exit(1);}
    }
  }
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of triangle counting\n");
    fprintf(stderr,"Number of rows            (-n) %"PRId64"\n", numrows);
    fprintf(stderr,"random seed               (-s)= %d\n", seed);
    fprintf(stderr,"Graph model:\n");
    fprintf(stderr,"  from file               (-f [%s])\n", filename); 
    fprintf(stderr,"  flat model (dens)       (-e [%g]\n", edge_prob);
    fprintf(stderr,"  geometric model (dens)  (-g [%g]\n", edge_prob);
    fprintf(stderr,"  kronecker prod (string) (-k [%s]\n", kron_list);
    fprintf(stderr,"models_mask               (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files                (-D)\n");
    fprintf(stderr,"quiet                     (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }
  if(!quiet){
    printf("Input matrix stats:\n");
    spmat_stats(mat);
    if( model == KRONECKER )
      printf("Kronecker model should have %"PRId64" triangles\n", tri_count_kron_graph(kron_args));
  }
    
  if(dump_files){
    dump_matrix(mat, 20, "mat.out");
  }

  double laptime = 0.0;
  for(use_model=1; use_model < ALL; use_model *=2 ){
    triangles = 0;
    switch( use_model & models_mask ){
    case GENERIC:
      if( !quiet ) printf("Generic      Triangle: ");
        laptime = triangles_matrix(&triangles, mat);
      break;
    default:
      continue;
    }
    if( !quiet ) 
      printf(" %12"PRId64" triangles : %8.3lf seconds\n", triangles, laptime);
  }
 
  if(model == KRONECKER) 
    free(kron_args);
  clear_matrix(mat);
  return(0);
}

