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
  sparsemat_t * At = transpose_matrix(A);
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
  int64_t readmat = 0;
  char filename[256]={"filename"};
  sparsemat_t *mat;
  double laptime = 0.0;
  double edge_density = 0.0;
  int64_t nz_per_row = 10;
  graph_model model = FLAT;
  enum FLAVOR {GENERIC=1, ALL=2};
  uint32_t use_model;
  uint32_t models_mask = ALL - 1;
  int printhelp=0;
  int quiet=0;
  int dump_files=0;
  
  int opt; 
  while( (opt = getopt(argc, argv, "hm:n:s:FGDM:e:f:qz:")) != -1 ) {
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
  if(printhelp){
    fprintf(stderr,"Inputs\n");
    fprintf(stderr,"-f <file> : Read a matrix from a MM file\n");
    fprintf(stderr,"Or... the adjacency matrix for a random graph.\n");    
    fprintf(stderr,"-n <numrows>: Number of rows\n");
    fprintf(stderr,"-F or -G: Random graph model (Flat or Geometric)\n");
    fprintf(stderr,"-e <prob> or -z <avg_nz_per_row>: Specify nz density by edge probability or average nonzeros per row\n");
    fprintf(stderr,"-s <seed>: RNG seed\n");
    fprintf(stderr,"-M <mask>: which flavors to run\n");
    fprintf(stderr,"-D : write out matrix (to dump.out)\n");
    fprintf(stderr,"-q : quiet mode\n");
    return(0);
  }
  
  if(!quiet ) {
    fprintf(stderr,"Running C version of transpose matrix\n");
    if(readmat == 1)
      fprintf(stderr,"Reading a matrix from file (-f [%s])\n", filename);
    else{
      if(model == FLAT)
        fprintf(stderr,"flat model           (-F)\n");
      else        
        fprintf(stderr,"geometric model      (-G)\n");
      fprintf(stderr,"Number of rows       (-n) %"PRId64"\n", numrows);
      fprintf(stderr,"edge_density         (-e)= %lg\n", edge_density);
      fprintf(stderr,"nz_per_row           (-z)= %"PRId64"\n", nz_per_row);
      fprintf(stderr,"random seed          (-s)= %d\n", seed);
    }
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files           (-D)= %d\n", dump_files);
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
  }
  
  if( readmat ) {
     mat = read_matrix_mm(filename);     
     if(!mat){printf("ERROR: transpose_matrix: read_matrix (%s) failed\n", filename); exit(1);}
  } else {
    mat = random_graph(numrows, model, DIRECTED, NOLOOPS, edge_density, seed);
    if(!mat){printf("ERROR: transpose_matrix: erdos_renyi_graph failed\n"); exit(1);}
  }

  if(!quiet) {
    printf("Input matrix stats:\n");
    spmat_stats(mat);
  }
  if(dump_files)
    dump_matrix(mat, 20, "dump.out");
  
  for( use_model=1; use_model < 2; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC:
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

