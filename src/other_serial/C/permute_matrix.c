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
  double edge_density = 0.0;
  int32_t dump_files = 0;
  int32_t readmat = 0;
  char filename[256]={"filename"};
  
  uint32_t seed = 101892;
  double laptime = 0.0;
  graph_model model = FLAT;
  enum FLAVOR {GENERIC=1, ALL=2};
  uint32_t use_model;
  uint32_t models_mask = ALL - 1;
  int printhelp=0;
  int quiet=0;
  int64_t nz_per_row = 10;

  int opt; 
  while( (opt = getopt(argc, argv, "hn:m:s:FGM:e:f:Dqz:")) != -1 ) {
    switch(opt) {
    case 'D': dump_files = 1; break;
    case 'e': sscanf(optarg,"%lg", &edge_density); break;
    case 'f': readmat = 1; sscanf(optarg, "%s", filename); break;
    case 'F': model = FLAT; break;
    case 'G': model = GEOMETRIC; break;
    case 'h': printhelp = 1; break;
    case 'M': sscanf(optarg,"%d" , &models_mask);  break;
    case 'n': sscanf(optarg,"%"SCNd64 ,&numrows );  break;
    case 's': sscanf(optarg,"%d" ,&seed );  break;
    case 'q': quiet = 1; break;
    case 'z': sscanf(optarg, "%"SCNd64, &nz_per_row); break;
    default:  break;
    }
  }
  
  if(!readmat){
    resolve_edge_prob_and_nz_per_row(&edge_density, &nz_per_row, numrows, NOLOOPS);
  }
  
  sparsemat_t *mat;
  if( readmat ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: Read graph from %s Failed\n", filename); exit(1);}
  } else {
    //mat = erdos_renyi_graph(numrows, er_prob, ER_GRAPH_DIRECT, seed);
    mat = random_graph(numrows, model, DIRECTED, NOLOOPS, edge_density, seed);
    if(!mat){printf("ERROR: erdos_renyi_graph failed!\n"); exit(1);}
  }
  if(dump_files) dump_matrix(mat, 20, "mat.out");
  
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of permute_matrix\n");
    fprintf(stderr,"help                 (-h)\n");
    fprintf(stderr,"numrows              (-n)= %"PRId64"\n", mat->numrows);
    fprintf(stderr,"numcols              (-n)= %"PRId64"\n", mat->numcols);
    fprintf(stderr,"edge_density         (-e)= %lg\n", edge_density);
    fprintf(stderr,"nz_per_row           (-z)= %"PRId64"\n", nz_per_row);
    if(model == FLAT)
      fprintf(stderr,"flat model           (-F)\n");
    else
      fprintf(stderr,"geometric model      (-G)\n");
    fprintf(stderr,"readfile             (-f [%s])\n", filename); 
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files           (-D)\n");
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
    if(printhelp)
      return(0);
  }
  if(!quiet)spmat_stats(mat);
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

