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
 * \brief This routine implements the agi variant of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_agi(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) {
  //T0_printf("Running Toposort with UPC ...");
  int64_t nr = mat->numrows;

  assert(mat->numrows == mat->numcols);
  int64_t col2;

  SHARED int64_t * queue  = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t * lqueue  = lgp_local_part(int64_t, queue);
  SHARED int64_t * rowsum = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowsum = lgp_local_part(int64_t, rowsum);
  SHARED int64_t * rowcnt = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowcnt = lgp_local_part(int64_t, rowcnt);
  
  SHARED int64_t * S_end = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t *lS_end = lgp_local_part(int64_t, S_end);
  int64_t start, end;

  SHARED int64_t * pivots = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t *lpivots = lgp_local_part(int64_t, pivots);
  lpivots[0] = 0L;
  int64_t i, j;
   
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  start = end = 0;
  for(i = 0; i < mat->lnumrows; i++){
    lrowsum[i] = 0L;
    lrowcnt[i] = mat->loffset[i+1] - mat->loffset[i];
    if(lrowcnt[i] == 1)
      lqueue[end++] = i;    
    for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++)
      lrowsum[i] += mat->lnonzero[j];
  }
  //S_end[MYTHREAD] = end;
  lS_end[0] = end;
  lgp_barrier();

  // we a pick a row with a single nonzero = col.
  // setting rperm[pos] = row and cprem[pos] = col
  // moves that nonzero to the diagonal.
  // Now, cross out that row and col by decrementing 
  //  the rowcnt for any row that contains that column
  // repeat


  // Do it in levels:
  // Per level -- everyone picks their degree one rows
  //              and claims their spots in the permutations,
  //              then they all go back and update the 
  //              appropriate rowcnts for that level
  //                 
  // cool trick: rowsum[i] is the sum of all the column indices in row i,
  //             so rowsum[i] is the column indice we want when rowcnt = 1;

  double t1 = wall_seconds();

  int work_to_do = lgp_reduce_add_l(end - start);
  int64_t pos, l_row, S_col, S_row, S_indx, colcount, level = 0;
  int64_t old_row_sum, old_row_cnt, l_pos;

  while(work_to_do) {
    level++;
    while( start < end ){
      l_row = lqueue[start++];
      S_col = lrowsum[l_row];  // see cool trick

      // claim our spot on the diag     
      pos = lgp_fetch_and_inc(pivots, 0);
      lgp_put_int64(rperm, l_row*THREADS + MYTHREAD, nr - 1 - pos);
      lgp_put_int64(cperm, S_col, nr - 1 - pos);

      // use the global version of tmat to look at this column (tmat's row)
      // tmat->offset[S_col] is the offset local to S_col%THREADS
      S_indx = lgp_get_int64(tmat->offset, S_col) * THREADS + S_col % THREADS;
      colcount = lgp_get_int64(tmat->offset, S_col+THREADS) - lgp_get_int64(tmat->offset,S_col);
      for(j=0; j < colcount; j++) {
        S_row = lgp_get_int64(tmat->nonzero, S_indx + j*THREADS );
        assert((S_row) < mat->numrows);
        old_row_cnt = lgp_fetch_and_add(rowcnt, S_row, -1L);
        old_row_sum = lgp_fetch_and_add(rowsum, S_row, (-1L)*S_col);
        if( old_row_cnt == 2L ) {
          l_pos = lgp_fetch_and_inc(S_end, S_row % THREADS);
          lgp_put_int64(queue, l_pos*THREADS + S_row%THREADS , S_row / THREADS);fflush(0);
        }
      }
    }
    lgp_barrier();
    assert( lS_end[0] >= end );
    end = lS_end[0];
    work_to_do = lgp_reduce_add_l(end - start);
  }

  level = lgp_reduce_max_l(level);
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  
  if(lgp_get_int64(pivots,0) != nr){
    printf("ERROR! toposort_matrix_upc_orig: found %"PRId64" pivots but expected %"PRId64"!\n", lgp_get_int64(pivots, 0), nr);
    exit(1);
  }
  T0_fprintf(stderr, "num levels = %ld ", level+1);
  lgp_all_free(queue);
  lgp_all_free(rowsum);
  lgp_all_free(rowcnt);
  lgp_all_free(S_end);
  lgp_all_free(pivots);
  //T0_printf("done\n");
  return(stat->avg);
}
