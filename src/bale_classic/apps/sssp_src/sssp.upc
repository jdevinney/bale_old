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
/*! \file sssp.upc
 * \brief Demo application that implements Single Source Shortest Path algorithms
 */

#include "sssp.h"

/*!
  \page sssp_page Single Source Shortest Path


We are given a weight adjacency matrix for a (directed) graph
and an initial starting vertex. 
The problem is to find the distance 
(the sum of the weights of the edges on the least heavy path)
to all other vertices in the graph.

 */

static void usage(void) {
  T0_fprintf(stderr,"\
Usage:\n\
sssp [-h][-a 0,1][-e prob][-K str][-f filename]\n\
- -e = p: specify the Edge probability p\n\
- -h print this message\n\
- -M mask is the or of 1,2,4,8,16 for the models: \n\
- -N = n: Specify the number of rows_per_thread in the matrix (if using the random_graph generator)\n\
- -f filename : Specify a filename containing a matrix in MatrixMarket format to read as input\n\
- -b = count: Specify the number of packages in an exstack(2) stack\n\
\n\
\n");
  lgp_finalize();
  lgp_global_exit(0);
}

int main(int argc, char * argv[])
{

  lgp_init(argc, argv);

  int64_t buf_cnt = 1024;
  int64_t models_mask = ALL_Models;  // default is running all models
  int64_t l_numrows = 10000;         // number of a rows per thread
  int64_t read_graph = 0L;           // read graph from a file
  char filename[64];
  int64_t cores_per_node = 0;
  
  double t1;
  int64_t i, j;
  int64_t alg = 0;
  double erdos_renyi_prob = 0.1;
  int64_t nz_per_row = -1;
  graph_model model = FLAT;
  int64_t seed = 1231;
  int printhelp = 0;
  int opt; 
  while( (opt = getopt(argc, argv, "hb:c:M:n:f:FGa:e:K:s:Z:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'a': sscanf(optarg,"%"PRId64"", &alg); break;
    case 'b': sscanf(optarg,"%"PRId64"", &buf_cnt);  break;
    case 'c': sscanf(optarg,"%"PRId64"" ,&cores_per_node); break;
    case 'e': sscanf(optarg,"%lg", &erdos_renyi_prob); break;
    case 'f': read_graph = 1; sscanf(optarg,"%s", filename); break;      
    case 'F': model = FLAT; break;
    case 'G': model = GEOMETRIC; break;  
    case 'M': sscanf(optarg,"%"PRId64"", &models_mask);  break;
    case 'n': sscanf(optarg,"%"PRId64"", &l_numrows); break;
    case 's': sscanf(optarg,"%"PRId64"", &seed); break;
    default:
      T0_fprintf(stderr, "ERROR: Illegal usage\n");
      return(1);
    }
  }
  if( printhelp ) usage(); 

  int64_t numrows = l_numrows * THREADS;
  nz_per_row = erdos_renyi_prob * numrows;
  
  T0_fprintf(stderr,"Running triangle on %d threads\n", THREADS);
  if(!read_graph && !gen_kron_graph){
    T0_fprintf(stderr,"Number of rows per thread   (-N)  %"PRId64"\n", l_numrows);
    T0_fprintf(stderr,"Edge prob                   (-e)  %g\n", erdos_renyi_prob);
    T0_fprintf(stderr,"Graph Model           (-F or -G)  %s\n", (model == FLAT ? "FLAT" : "GEOMETRIC"));
    T0_fprintf(stderr,"Seed                        (-s)  %"PRId64"\n", seed); 
  }
  T0_fprintf(stderr,"Model mask (M) = %"PRId64" (should be 1,2,4, TODO: \n", models_mask);  
  
  lgp_barrier();
  
  if( printhelp )
    return(0);

  double correct_answer = -1;
 
  // TODO: Don't know what we want to require for the input. 
  sparsemat_t *A, *L, *U;
  if(read_graph){
    A = read_matrix_mm_to_dist(filename);
    if(!A)
      lgp_global_exit(1);
    T0_fprintf(stderr,"Reading file %s...\n", filename);
    T0_fprintf(stderr, "A has %"PRId64" rows/cols and %"PRId64" nonzeros.\n", A->numrows, A->nnz);
  }else{
    L = random_graph(numrows, model, DIRECTED_WEIGHTED, 0, erdos_renyi_prob, seed);
  }
  // we should check that A has the required format

  lgp_barrier();
  
  T0_fprintf(stderr,"L has %"PRId64" rows/cols and %"PRId64" nonzeros.\n", L->numrows, L->nnz);
  
  T0_fprintf(stderr,"Run sssp ...\n");

  int64_t use_model;
  double laptime = 0.0;

  // allocate dist array
  
  for( use_model=1L; use_model < 4; use_model *=2 ) {

    switch( use_model & models_mask ) {
    case AGI_BELLMAN:
      T0_fprintf(stderr,"    Bellman-Ford  AGI: ");
      laptime = sssp_bellman_agi(dist, dmat, 0); 
      break;
    
    case EXSTACK_BELLMAN:
      T0_fprintf(stderr,"  Exstack: ");
      laptime = sssp_bellman_exstack(dist, dmat, 0);
      break;
    }
    
    lgp_barrier();
    T0_fprintf(stderr,"  %8.3lf seconds.\n", laptime);
    // TODO: Check result
  }
  
  lgp_barrier();
  lgp_finalize();
  return(0);
}
