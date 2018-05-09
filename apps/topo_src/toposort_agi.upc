/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
 *****************************************************************/ 
/*! \file toposort_agi.upc
 * \brief The intuitive implementation of toposort that does uses 
 * atomics and generic global references
 */

#include "toposort.h"

#define toposort_matrix_upc_oo toposort_matrix_upc
//#define toposort_matrix_upc_orig toposort_matrix_upc

/* To play with some higher level abstractions for queues and 
 * moving across the rows of a sparse matrix, we have written 
   toposort_matrix_upc_oo which uses the following support functions.
 * the toposort_matrix_upc_orig  uses more the explicit code.
 * the main toposort functions just calls toposort_matrix_upc,
   hence the macros above.
 */

/*! \struct swlrQ_t 
 * \brief A structure to hold a set of queues that can be filled by other threads
    and emptied only by MYTHREAD
    We used shared arrays for the items and both shared and a local variable to hold the 
    indices for the heads of the queues. Since one dequeues item locally, we only need
    a local variable for the tails.
    The race condition to add items to a queue is handled with a fetch_and_add to the shared head.
    The read after write race condition to dequeue items up to head is side stepped by using a 
    snap shot of the queues and only updating the snap shots between barriers. 
  */  
typedef struct swlrQ {
   SHARED int64_t * S_queue;  //!< space for all the queues
   int64_t * l_queue;         //!< localized pointers to the individual queues
   SHARED int64_t *S_head;    //!< shared version of the heads of each queue (inherently suffers races conditions)
   int64_t l_head;            //!< snap shot version of the head of the local queue
                              //  this is only set by qrab_swlrQ.
   int64_t l_tail;            //!< tail of the local queue
} swlrQ_t;

/*!
 * \brief initialize the shared write local read queue
 * \param nitems the number of items on each local queue
 */
swlrQ_t * 
init_swlrQ(int64_t nitems)
{
  swlrQ_t * sq = calloc(1, sizeof(swlrQ_t));
  sq->S_queue  = lgp_all_alloc(((nitems/THREADS+1)*THREADS), sizeof(int64_t));
  sq->l_queue  = lgp_local_part(int64_t, sq->S_queue);

  sq->S_head = lgp_all_alloc(THREADS, sizeof(int64_t));
  lgp_put_int64(sq->S_head, MYTHREAD, 0);
  sq->l_head = 0;
  sq->l_tail = 0;

  lgp_barrier();
  return(sq);
}

/*!
 * \brief clears the shared write local read queue
 * \param *sq pointer the queue
 */
void *
clear_swlrQ( swlrQ_t * sq )
{
  if( !sq ) return(NULL);
  lgp_all_free(sq->S_queue);
  lgp_all_free(sq->S_head);
  return(NULL);
}

/*!
 * \brief grabs a local snap shot version of the queue by setting the "local head" of the queue
 *      Allows one to work with local queue while other threads are use S_head to add more 
 *      items to the queue.
 * \param *nitems place to hold the number of items on the local queue
 * \return the total number of items in all the threads queues
 */
int64_t
grab_swlrQ(swlrQ_t * sq, int64_t *nitems)
{
   lgp_barrier();  // finish inflight pushes
   sq->l_head = lgp_get_int64(sq->S_head, MYTHREAD);
   *nitems = sq->l_head - sq->l_tail;
   return( lgp_reduce_add_l(*nitems) ); // implicit barrier here is necessary
}

/*!
 * \brief pushs an item (int64_t) on to another threads queue
 * \param owner the target thread for the item
 * \param item the int64_t we are pushing
 * \return true (might return false if we did some error checking)
 */
bool
en_swlrQ(swlrQ_t * sq, int64_t owner, int64_t item)
{
   int64_t l_pos;
   
   // the race for the spot on the another threads queue
   // is handled with an atomic add to it head
   l_pos = lgp_fetch_and_inc(sq->S_head, owner);
   lgp_put_int64(sq->S_queue, l_pos*THREADS + owner, item);
   //printf("%d: >> %d to %d into %d\n",MYTHREAD, item, owner, l_pos);
   return true;
}

/*!
 * \brief pull the tail element from a thread's queue.
 * \param *sq pointer to the shared write local read queue struct
 * \param *ret_item a place to put the int64_t item from the queue
 * \return false if the queue was empty, else true
 */
bool
de_swlrQ(swlrQ_t * sq, int64_t *ret_item )
{
  if( sq->l_head > sq->l_tail ) {
     *ret_item = sq->l_queue[sq->l_tail];
     sq->l_tail++;
     //printf("%d: << %d\n",MYTHREAD, *ret_item);
     return true;
  } else {
     //printf("%d: << empty\n",MYTHREAD);
     return false;
  }
}

/*!
 * \brief returns the number of nonzeros in a row of the localize part of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param l_row row index into local version (loffset) of the offset array
 * \return the number of nonzero in the request row
 */
int64_t 
rowcount_l( sparsemat_t *mat, int64_t l_row )
{
   return( mat->loffset[l_row+1] - mat->loffset[l_row] );
}

/*!
 * \brief returns the number of nonzeros in a row of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param S_row global row index into the shared offset array
 * \return the number of nonzero in the request row
 */
int64_t 
rowcount_S( sparsemat_t *mat, int64_t S_row )
{
   return( mat->offset[S_row+THREADS] - mat->offset[S_row] );
}

/*!
 * \brief enumerates the nonzeros in a row of the localize part of a sparse matrix
 * \param *ret int64_t place to put the return value
 * \param *mat pointer to the sparse matrix 
 * \param l_row row index into local version of the sparse matrix
 * \return true if returning a valid nonzero, false at the end of the row
 * Hidden State: uses l_enum_row, l_enum_idx, l_enum_nstop as state variables
 * to keep track of which row is be referenced and how far along the row each
 * successive call has moved
 */
bool
row_enum_l( int64_t *ret, sparsemat_t *mat, int64_t l_row )
{
  assert( l_row < mat->lnumrows );
  if( mat->l_enum_row == l_row ){ // keep listing them
    mat->l_enum_idx += 1;
  } else {                                // first call to this row
    mat->l_enum_row = l_row;
    mat->l_enum_idx    = mat->loffset[l_row];
    mat->l_enum_nstop  = mat->loffset[l_row+1];
  }
  if( mat->l_enum_idx < mat->l_enum_nstop ){
    *ret = mat->lnonzero[mat->l_enum_idx];
    return true;
  } else {
    *ret = -1;
    mat->l_enum_row = -1;
    mat->l_enum_idx = -1;
    return false;
  }
}

/*!
 * \brief enumerates the nonzeros in a row of the of a sparse matrix as a global data structure
 * \param *ret int64_t place to put the return value
 * \param *mat pointer to the sparse matrix 
 * \param S_row global row index into rows of the sparse matrix
 * \return true if returning a valid nonzero, false at the end of the row
 * Hidden State: uses S_enum_row, S_enum_idx, S_enum_nstop as state variables
 * to keep track of which row is be referenced and how far along the row each
 * successive call has moved
 */
bool
row_enum_S( int64_t *ret, sparsemat_t *mat, int64_t S_row )
{
  if( mat->S_enum_row == S_row ){ // keep enumerating them
    mat->S_enum_idx += 1;
  } else {   // first call to this row
    mat->S_enum_row = S_row;
    mat->S_enum_idx = lgp_get_int64(mat->offset, S_row);
    mat->S_enum_nstop  = lgp_get_int64(mat->offset, S_row+THREADS);
  }
  if( mat->S_enum_idx < mat->S_enum_nstop ){
    *ret = lgp_get_int64(mat->nonzero, mat->S_enum_idx * THREADS +  S_row%THREADS );
    return true;
  } else {
    *ret = -1;
    mat->S_enum_row = -1;
    mat->S_enum_idx = -1;
    return false;
  }
}

/*!
 * \brief This routine implements an agi variant of toposort with code to encapsulate working 
 *   with the queue of degree one rows and moving across the rows a transpose matrix.
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time

 */
double toposort_matrix_upc_oo(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat)
{
  //T0_fprintf(stderr,"Running Toposort with UPC ....");
  int64_t l_row, S_col, S_row;
  int64_t old_row_sum, old_row_cnt;
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;

  swlrQ_t * sq = init_swlrQ(nr);

  int64_t *l_rperm = lgp_local_part(int64_t, rperm);
  SHARED int64_t * rowsum = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowsum = lgp_local_part(int64_t, rowsum);
  SHARED int64_t * rowcnt = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowcnt = lgp_local_part(int64_t, rowcnt);
  
  double t1 = wall_seconds();
   
  /* initialize rowsum, rowcnt, and queue holding the degree one rows */
  for(int64_t i = 0; i < mat->lnumrows; i++){
    lrowcnt[i] = rowcount_l(mat, i);
    if(lrowcnt[i] == 1)
      en_swlrQ(sq, MYTHREAD, i);
    lrowsum[i] = 0L;
    while( row_enum_l(&S_col, mat, i) )   // foreach S_col in l_row i of mat
      lrowsum[i] += S_col;
  }

  int64_t mypivs, n_pivs; // number of pivots on my queue, all queues (done when n_pivs == 0)
  int64_t pivs_found = 0; // total pivots found so far
  int64_t pos;
  while( n_pivs = grab_swlrQ(sq, &mypivs) ) {  // getting a local copy of the queue defines a level

    pos = pivs_found + lgp_prior_add_l(mypivs); // pivots in this level and thread start here
    while( de_swlrQ(sq, &l_row) ) {
      S_col = lrowsum[l_row];  // see cool trick

      l_rperm[l_row]  = nr - 1 - pos;
      lgp_put_int64(cperm, S_col, nc - 1 - pos);
      pos++;

      while( row_enum_S(&S_row, tmat, S_col ) ) {       // foreach S_row in S_col of tmat
        old_row_cnt = lgp_fetch_and_add(rowcnt, S_row, -1L);
        old_row_sum = lgp_fetch_and_add(rowsum, S_row, (-1L)*S_col);
        if( old_row_cnt == 2L ) {
          en_swlrQ(sq, S_row%THREADS, S_row/THREADS);
        }
      }
    }
    pivs_found += n_pivs;
  }
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  if(pivs_found != nr){
    printf("ERROR! toposort_matrix_upc: found %ld pivots but expected %ld!\n", pivs_found, nr);
    exit(1);
  }

  clear_swlrQ(sq); free(sq);
  lgp_all_free(rowsum);
  lgp_all_free(rowcnt);
  return(stat->avg);
}


/*!
 * \brief This routine implements the agi variant of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */
double toposort_matrix_upc_orig(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat)
{
  //T0_printf("Running Toposort with UPC ...");
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
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
    //printf("rowsum[%ld] = %ld\n", i*THREADS + MYTHREAD, lrowsum[i]);fflush(0);
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
      lgp_put_int64(cperm, S_col, nc - 1 - pos);

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
          lgp_put_int64(queue, l_pos*THREADS + S_row%THREADS , S_row / THREADS);
        }
      }
    }
    lgp_barrier();
    assert( lS_end[0] >= end );
    end = lS_end[0];
    work_to_do = lgp_reduce_add_l(end - start);
  }
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  
  if(lgp_get_int64(pivots,0) != nr){
    printf("ERROR! toposort_matrix_upc: found %ld pivots but expected %ld!\n", pivots[0], nr);
    exit(1);
  }
  lgp_all_free(queue);
  lgp_all_free(rowsum);
  lgp_all_free(rowcnt);
  lgp_all_free(S_end);
  lgp_all_free(pivots);
  //T0_printf("done\n");
  return(stat->avg);
}
