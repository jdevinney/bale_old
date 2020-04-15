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
/*! \file toposort_agi.upc
 * \brief The intuitive implementation of toposort that does uses 
 * atomics and generic global references
 */

#include "toposort.h"

/*!
 * \brief This routine implements a cool variant version of the agi version of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_cooler(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) {
  //T0_printf("Running Toposort with UPC ...");
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  int64_t col;

  SHARED int64_t * queue  = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t * lqueue  = lgp_local_part(int64_t, queue);
  SHARED int64_t * rowstat = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowstat = lgp_local_part(int64_t, rowstat);
  
  SHARED int64_t * pivots = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t *lpivots = lgp_local_part(int64_t, pivots);
  lpivots[0] = 0L;
  int64_t i, j;
   
  /* initialize rowstat array */
  for(i = 0; i < mat->lnumrows; i++){
    lrowstat[i] = 0L;
    int64_t cnt = mat->loffset[i+1] - mat->loffset[i];
    for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++)
      lrowstat[i] += (nc + mat->lnonzero[j]);
  }

  lgp_barrier();

  
  // cooler trick: rowstat[i] is the sum (nc + i) for all i in each row,

  double t1 = wall_seconds();

  int64_t rows_left = mat->lnumrows;
  int64_t d1bound = nc+nc;
  int64_t done = nc*100;
  int64_t old_val, pos, S_indx, colcount, S_row;
  while(rows_left) {
    int64_t start, end;
    start = end = 0;
    
    // count the number of d1 rows you see now
    int64_t d1 = 0;
    for(i = 0; i < mat->lnumrows; i++){      
      if(lrowstat[i] < d1bound){
        d1++;
        lqueue[end++] = i;
      }
    }

    // reserve d1 spots in the perm positions
    pos = lgp_fetch_and_add(pivots, 0, d1);
    rows_left -= d1;
    
    while(start < end){
      i = lqueue[start++];
      col = lrowstat[i] - nc;
      old_val = lrowstat[i];
      lrowstat[i] = done; //mark it so we won't see it again

      // add row and col to perms
      lgp_put_int64(rperm, i*THREADS + MYTHREAD, nr - 1 - pos);
      lgp_put_int64(cperm, col,                  nc - 1 - pos);
      pos++;
      
      S_indx = lgp_get_int64(tmat->offset, col) * THREADS + col % THREADS;
      colcount = lgp_get_int64(tmat->offset, col+THREADS) - lgp_get_int64(tmat->offset, col);
      for(j=0; j < colcount; j++) {
        S_row = lgp_get_int64(tmat->nonzero, S_indx + j*THREADS );
        assert((S_row) < mat->numrows);
        lgp_fetch_and_add(rowstat, S_row, -old_val);
      }
    }
    //lgp_barrier();
  }
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  
  if(lgp_get_int64(pivots,0) != nr){
    printf("ERROR! toposort_matrix_upc_orig: found %"PRId64" pivots but expected %"PRId64"!\n", pivots[0], nr);
    exit(1);
  }
  lgp_all_free(queue);
  lgp_all_free(rowstat);
  lgp_all_free(pivots);
  //T0_printf("done\n");
  return(stat->avg);
}
