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
#include <getopt.h>
#include <libgetput.h>
#include <spmat.h>

/*! \file write_sparse_matrix.upc
 * \brief Demo program that runs the variants of write_sparse_matrix kernel. This program
 * writes a sparsemat_t to disk in a binary format.
 */

/*! 
\page write_sparse_matrix_page Write Sparse Matrix

Demo program that runs the variants of write_sparse_matrix kernel. It first generates 
a random matrix in FLAT mode and then it writes this matrix to disk
in a directory called 'write_sparse_test'.

We define a sparse matrix dataset to be the following:
 - It lives in a directory of its own
 - It contains one ASCII file called 'metadata' which contains the number of rows, columns and nonzeros in the matrix.
 - It contains N binary 'rowcnt' files. These files contain the number of nonzeros in each row for rows 0..A->numrows
 - It contains N binary 'nonzero' files. These files contain the nonzeros in each row and are ordered by row.

This application is interesting because, as implemented, it requires
us to shuffle the rows of the matrix from cyclic to block. That is,
the nonzeros for row i is stored on PE (i % THREADS) in the
sparsemat_t data structure.  However, we wish PE 0 to write out the
first block of (approx) A->numrows/THREADS rows. One way to do this
would be to call permute_sparse_matrix to get a copy of the matrix
whose rows are distributed in the block layout. However, that would
require 2x the storage space of the matrix. We don't want the
write_sparse_matrix routine to have this requirement. That is where
things get interesting. We have a fixed buffer for writing on each PE
and each PE collects or is sent nonzero data to write in their current
buffer. In AGI (where the PEs just get the data) or synchronous
exstack, this is easy. We don't have a nice way of doing this with
exstack2 or asynchronous conveyors yet. The reason this is a challenge
for asynchronous methods is that PEs can get into a deadlock waiting
for the records they need to complete a write buffer.

See files spmat_agi.upc, spmat_exstack.upc
for the source for the kernels.

 * Usage:
 * write_sparse_matrix [-h][-b count][-M mask][-n num][-T tabsize][-c num]
 * - -h Print usage banner only
 * - -b count a hint of the number of packages in an exstack buffer (may need to be larger)
 * - -e=p Set the Erdos-Renyi probability to p.
 * - -M=m Set the models mask (1,2,4,8,16,32 for gets,classic,exstack2,conveyor,alternate)
 * - -n=n Set the number of rows per PE to n (default = 1000).
 * - -s=s Set a seed for the random number generation.
 * - -Z=z Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).
 */

static void usage(void) {
  T0_fprintf(stderr,"\
This is a demo program that runs various implementations of the write_sparse_matrix kernel.\n\
Usage:\n\
histo [-h][-b count][-M mask][-n num][-T tabsize][-c num]\n\
-h prints this help message\n\
-b count is the number of packages in an exstack(2) buffer\n\
-M mask is the or of 1,2,4,8,16 for the models: agi,exstack,exstack2,conveyor,alternate\n\
-e=p Set the Erdos-Renyi probability to p.\n\
-M=m Set the models mask (1,2,4,8,16,32 for gets,classic,exstack2,conveyor,alternate)\n\
-n=n Set the number of rows per PE to n (default = 1000).\n\
-s=s Set a seed for the random number generation.\n\
-Z=z Set the avg number of nonzeros per row to z (default = 10, overrides Erdos-Renyi p).\n\
\n");
  lgp_global_exit(0);
}

int main(int argc, char * argv[])
{
  lgp_init(argc, argv);
  
  int64_t i;
  int64_t models_mask = 0xF;
  int printhelp = 0;
  double erdos_renyi_prob = 0.0;
  int64_t nz_per_row = -1;
  int64_t buf_cnt = 1024;
  int64_t l_numrows = 10000;
  int64_t numrows;
  int64_t seed = 101892+MYTHREAD;
  sparsemat_t * inmat;
  int64_t cores_per_node = 1; 

  int opt; 
  while( (opt = getopt(argc, argv, "hb:c:e:n:M:s:Z:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%"PRId64"", &buf_cnt);  break;
    case 'c': sscanf(optarg,"%"PRId64"" ,&cores_per_node); break;
    case 'e': sscanf(optarg,"%lf", &erdos_renyi_prob);  break;
    case 'n': sscanf(optarg,"%"PRId64"", &l_numrows);   break;
    case 'M': sscanf(optarg,"%"PRId64"", &models_mask);  break;
    case 's': sscanf(optarg,"%"PRId64"", &seed); break;
    case 'Z': sscanf(optarg,"%"PRId64"", &nz_per_row);  break;
    default:  break;
    }
  }
  if(printhelp) usage();
  
  numrows = l_numrows * THREADS;

  /* set erdos_renyi_prob and nz_per_row to be consistent */
  if(nz_per_row == -1 && erdos_renyi_prob == 0){
    nz_per_row = 10;
  }else if(nz_per_row == -1){
    nz_per_row = erdos_renyi_prob*numrows;
  }
  erdos_renyi_prob = (2.0*(nz_per_row - 1))/numrows;
  if(erdos_renyi_prob > 1.0)
    erdos_renyi_prob = 1.0;

  T0_fprintf(stderr,"Running write_sparse_matrix on %d threads\n", THREADS);
  T0_fprintf(stderr,"buf_cnt (stack size)        (-b) = %"PRId64"\n", buf_cnt);
  T0_fprintf(stderr,"Erdos-Renyi edge probability(-e) = %lf\n", erdos_renyi_prob);
  T0_fprintf(stderr,"rows per PE (-n)                 = %"PRId64"\n", l_numrows);
  T0_fprintf(stderr,"models_mask (-M)                 = %"PRId64" or one of 1,2,4,8,16 for gets,classic,exstack2,conveyor,alternate\n", models_mask);
  T0_fprintf(stderr,"seed (-s)                        = %"PRId64"\n", seed);
  T0_fprintf(stderr,"Avg # of nonzeros per row   (-Z) = %"PRId64"\n", nz_per_row);


  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  
  inmat = random_graph(numrows, FLAT, DIRECTED, 0, erdos_renyi_prob, seed + 2);
  if(inmat == NULL){
    T0_printf("ERROR: inmat is null!\n");
    return(-1);
  }

  char * datadir = calloc(64, sizeof(char));

  int64_t use_model;
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & models_mask ) {
    case AGI_Model:
      sprintf(datadir,"%s","write_matrix_test_agi");
      write_sparse_matrix_agi(datadir, inmat);
      T0_fprintf(stderr, "AGI:     ");
      break;
    case EXSTACK_Model:
      sprintf(datadir,"%s","write_matrix_test_exstack");
      write_sparse_matrix_exstack(datadir, inmat, buf_cnt);
      //read_sparse_matrix_agi(datadir);
      T0_fprintf(stderr, "Exstack: ");
      break;
    case EXSTACK2_Model:
      continue;
      //T0_fprintf(stderr, "Exstack2:");
      //break;
    case CONVEYOR_Model:
      continue;
      //T0_fprintf(stderr, "Conveyor:");
      //break;    
    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      break;
    case 0:
      continue;
    }
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    T0_fprintf(stderr, "%8.3lf\n", stat->avg);  
  }
  
  free(datadir);
  clear_matrix(inmat);
  lgp_barrier();
  lgp_finalize();
  return(error);
}
