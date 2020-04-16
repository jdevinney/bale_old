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

/*! \file permute_matrix.c
 * \brief Demo program that runs the permute_matrix routine in the spmat library
 */

#include "spmat_utils.h"

/*! \page permute_matrix_page Permute a sparse matrix */

/*!
 * \brief A timing wrapper around the permute_matrix call in spmat_utils
 * \param *A  sparsemat_t holding the given matrix
 * \param *rp the row permutation
 * \param *cp the column permutation
 * \return run time
 * NB: This does charge for the time to initialize the permuted matrix.
 */
double permute_matrix_generic(sparsemat_t *A, int64_t *rp, int64_t *cp)
{
  double tm;
  tm = wall_seconds();
  sparsemat_t *Ap = permute_matrix(A, rp, cp);
  tm = wall_seconds() - tm;
  clear_matrix(Ap);
  return(tm);
}

int main(int argc, char * argv[]) 
{
  int64_t numrows = 10000;
  int64_t numcols = -1;
  double edge_prob = 0.05;
  int32_t dump_files = 0;
  int32_t readgraph = 0;
  char filename[256]={"filename"};
  
  uint32_t seed = 101892;
  double laptime = 0.0;

  enum FLAVOR {GENERIC=1, ALL=2};
  uint32_t use_model;
  uint32_t models_mask = ALL - 1;
  int printhelp=0;
  int quiet=0;

  int opt; 
  while( (opt = getopt(argc, argv, "hn:m:s:M:e:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows );  break;
    case 'm': sscanf(optarg,"%ld" ,&numcols );  break;
    case 's': sscanf(optarg,"%d" ,&seed );  break;
    case 'M': sscanf(optarg,"%d" , &models_mask);  break;
    case 'e': edge_prob = sscanf(optarg,"%lg", &edge_prob); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1;  break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }
  if(numcols == -1)
    numcols = numrows;

  sparsemat_t *mat;
  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: Read graph from %s Failed\n", filename); exit(1);}
  } else {
    //mat = erdos_renyi_graph(numrows, er_prob, ER_GRAPH_DIRECT, seed);
    mat = random_sparse_matrix(numrows, numcols, edge_prob, 0, seed);
    if(!mat){printf("ERROR: erdos_renyi_graph failed!\n"); exit(1);}
  }
  if(dump_files) dump_matrix(mat, 20, "mat.out");

  
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of permute_matrix\n");
    fprintf(stderr,"help                 (-h)\n");
    fprintf(stderr,"numrows              (-n)= %ld\n", mat->numrows);
    fprintf(stderr,"numcols              (-n)= %ld\n", mat->numcols);
    fprintf(stderr,"edge_prob            (-e)= %lg\n", edge_prob);
    fprintf(stderr,"readfile             (-f [%s])\n", filename); 
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files           (-D)\n");
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
    if(printhelp)
      return(0);
  }

  int64_t *rperm = rand_perm(numrows, seed); assert(rperm != NULL);
  int64_t *cperm = rand_perm(numrows, seed); assert(cperm != NULL);

  for( use_model=1; use_model < ALL; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC:
      if(!quiet) printf("generic permute matrix: ");
      laptime = permute_matrix_generic(mat, rperm, cperm);
    break;
    default:
       continue;
    }
    if(!quiet) printf("  %8.3lf seconds \n", laptime);
  }
  return(0);
}

