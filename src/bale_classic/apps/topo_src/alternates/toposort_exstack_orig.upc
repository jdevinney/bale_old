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

/*! \file toposort_exstack_orig.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"

typedef struct pkg_rowcol_t{
  int64_t row;    
  int64_t col;
}pkg_rowcol_t;

/*!
 * \brief This routine implements the exstack variant of toposort. It uses atomic ops
 * to construct the permutations.
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_exstack_orig(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat){
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  int64_t lnr = (nr + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnc = (nc + THREADS - MYTHREAD - 1)/THREADS;
  int64_t i,j,row,col,curr_col,pe,fromth,ret, pos;

  int64_t items_processed = 0;
  uint64_t type_mask = 0x8000000000000000;
  int64_t start, end;
  int64_t * lrowqueue  = calloc(lnr, sizeof(int64_t));
  int64_t * lcolqueue  = calloc(lnc, sizeof(int64_t));
  int64_t * lrowsum    = calloc(lnr, sizeof(int64_t));
  int64_t * lrowcnt    = calloc(lnr, sizeof(int64_t));

  pkg_rowcol_t pkg;

  exstack_t * ex = exstack_init(1024, sizeof(pkg_rowcol_t));
  if( ex == NULL ){return(-1.0);}
  
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  int64_t rownext, rowlast;
  int64_t colnext, collast;
  int64_t colstart,colend;
  rownext = rowlast = colnext = collast = colstart = colend = 0;

  for(i = 0; i < mat->lnumrows; i++){
    lrowsum[i] = 0L;
    lrowcnt[i] = mat->loffset[i+1] - mat->loffset[i];
    if(lrowcnt[i] == 1)
      lrowqueue[rowlast++] = i;
    for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++)
      lrowsum[i] += mat->lnonzero[j];
  }
  
  lgp_barrier();    
  
  SHARED int64_t * pivots = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t * lpivots = lgp_local_part(int64_t, pivots);
  lpivots[0] = 0L;
  lgp_put_int64(pivots, 0, 0);

  lgp_barrier();

  double t1 = wall_seconds();

  int64_t level = 0;
  while(exstack_proceed(ex, (items_processed == (lnr+lnc)))){
    
    int full = 0;
    
    /* push the column for each degree one row */
    while(rownext < rowlast){
      row = pkg.row = lrowqueue[rownext++];
      pkg.row |= type_mask;
      pkg.col = lrowsum[row];
      pe = pkg.col % THREADS;

      // claim our spot on the diag
      pos = lgp_fetch_and_inc(pivots,0);
      lgp_put_int64(rperm, pkg.row*THREADS + MYTHREAD, nr - 1 - pos);
      lgp_put_int64(cperm, pkg.col, nc - 1 - pos);

      items_processed++;
      if(exstack_push(ex, &pkg, pe) == 1){
        full = 1;
        break;
      }
    }
    
    if(!full){
      while(colnext <= collast){
        if(colstart == colend){
          if(colnext == collast) break;
          items_processed++;
          curr_col = lcolqueue[colnext++];
          colstart = tmat->loffset[curr_col];
          colend = tmat->loffset[curr_col + 1];
        }
        row = tmat->lnonzero[colstart];
        pkg.row = row/THREADS;
        pkg.col = curr_col*THREADS + MYTHREAD;
        pe = row % THREADS;
        if(exstack_push(ex, &pkg, pe) == 0)
          break;
        colstart++;
      
      }
    }
    
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      if(pkg.row & type_mask){
        lcolqueue[collast++] = pkg.col/THREADS;
      }else{
        lrowsum[pkg.row] -= pkg.col;
        lrowcnt[pkg.row]--;
        if(lrowcnt[pkg.row] == 1){
          lrowqueue[rowlast++] = pkg.row;
        }
      }
    }
  }
  
  t1 = wall_seconds() - t1;
  minavgmaxD_t stat[1];
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  free(lrowcnt);
  free(lrowsum);
  free(lrowqueue);
  free(lcolqueue);
  
  return(stat->avg);
}

