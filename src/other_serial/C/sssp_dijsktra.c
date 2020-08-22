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

#include "spmat_utils.h"

/*! \file sssp_dijsktra.c
 * \brief These are the function that implement and support two implementations 
 * of Dijsktra's Alg.
 * 
 * It is assumed that the Algorithm is known.
 * Excellent description of the algorithm are easily found.
 * To fix our terminology we say that the algorithm expands 
 * paths to "unvisited" vertices and assigns them tentative weights
 * and once the shortest path is found the vertex is said to be resolved.
 *   We assume:
 *     -- that the graph is weighted with positive weights
 *     -- if the graph is undirected, then mat contains two "directed edges" for
 *        each undirected edge.
 *
 */


/*!
 * \brief This routine implements the most naive version of Dijkstra's algorithm
 *
 * This first implementation uses a simple array to hold the
 * tentative weights and uses tricks with the weights to flag
 * the transition of vertices from unvisited to tentative to resolved.
 *
 * \param *mat sparsemat_t that holds the graph. 
 * \param v0 is the starting vertex
 * \return runtime
 */
double sssp_dijsktra_linear(sparsemat_t * mat, double *dist, int64_t v0)
{
  int64_t i, k, rn;
  int64_t numrows = mat->numrows;
  double minwt, curwt;
  int64_t minidx;

  for(i=0; i<numrows; i++)
    dist[i] = INFINITY;
  for(k = mat->offset[v0];  k < mat->offset[v0+1]; k++)
    dist[ mat->nonzero[k] ] = mat->value[k];
  dist[v0] = -0.0;
  //printf("dist[v0] = %lg\n", dist[v0]);

  while( 1 ) {
    // find the smallest tentative distance
    minwt = INFINITY;
    minidx = numrows;
    for(i=0; i<numrows; i++){
      if(dist[i] > 0 && dist[i] < minwt){
        minwt = dist[i];
        minidx = i;
      }
    }
    // done if all connected vertices have been resolved
    if( (minidx == numrows) || (minwt == INFINITY) ) 
      break;

    // update tentative distance from the current vertex
    curwt = dist[minidx];
    for(k = mat->offset[minidx];  k < mat->offset[minidx+1]; k++){
      rn = mat->nonzero[k];
      if(dist[rn] > curwt + mat->value[k])
				dist[rn] = curwt + mat->value[k];
		}
		//printf("dist[%"PRId64"] = %lg\n", minidx, curwt);
    dist[minidx] = -dist[minidx];
  }

  for(i=0; i<numrows; i++)
    dist[i] = -dist[i];

  return(0.0);
}

/*! Brief The min-heap implementation.
 * 
 * In the heap implementation the tentative weights are stored in
 * a min-heap data structure that is coupled with the rows of the matrix.
 * We will refer to the entries in the heap as nodes and the vertices 
 * in the graph as rows (of the matrix).  
 * The coupling is maintained by the arrays row[] and node[]. 
 * row[] and node[] are inverses of each other. 
 * row[n] gives a way look up the row when the heap node n reaches the root of the heap
 * node[r] gives a way look up the heap node that currently handles row r's tentative weight
 */
typedef struct PQ_t {
  int64_t numrows; // the number of vertices in the graph, rows in the matrix
  int64_t tail;    // the first node not actively in the queue;
  double  *wt;     // the corresponding weight that is being prioritized
  int64_t *row;    // the index of the row being represented by the given node in the heap
  int64_t *node;   // the index into the heap where the weight for the given row is stored
}PQ_t;

PQ_t * init_pqueue(int64_t numrows);
void bubble_up(PQ_t *pq, int64_t k);
int delete_root(PQ_t * pq);

void  heapify_pqueue(PQ_t * pq);
void print_queue(PQ_t * pq);
int64_t check_pqueue(PQ_t * pq);

/*! Brief initialize the heap
 * Note the indexing into the queue will be "1-up" so the parent and children are easy to calculate.
 * A given parent, p, has kids 2*p and 2*p+1, and the parent of kid, k, is \floor(k/2) 
 * As a convenience, we will use the zero node in the heap for the initial row (starting vertex).
 * We allocate an extra node so that numrows is a safe index.
 * That allows us to use node[] to flag rows: 
 *   node[row] == numrows means the row is unvisited
 *   node[row] == 0 means the row is resolved
 *   otherwise node[row] is the index of the entry in the heap that holds the tentative weight of row.
 */  
PQ_t * init_pqueue(int64_t numrows)
{
  PQ_t * pq = calloc(1, sizeof(PQ_t));
  pq->numrows = numrows;
  pq->tail = 1;
  pq->wt  = (double  *)calloc(numrows+1, sizeof(double));
  pq->row  = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  pq->node = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  return(pq);
}

/*! Brief bubble up a given node to return the heap to a legal state
 * Since a tentative weight only changes if it gets smaller
 * and we change only one node at a time, it is enough 
 * to bubble up the changed node until it is not smaller than its parent.
 * \param pq the priority queue
 * \param nd the index of the node that has changed
 */  
void bubble_up(PQ_t *pq, int64_t nd)
{
   double w;
   int64_t kid_row, par_row; 
   while( nd > 1 ){
     if( pq->wt[ nd/2 ] <= pq->wt[nd] )
       return; 
     //printf("swap kid %"PRId64" and parent %"PRId64"\n", nd/2, nd);
     w = pq->wt[nd];
     pq->wt[nd] = pq->wt[nd/2];
     pq->wt[nd/2] = w;
     
     kid_row = pq->row[nd];
     par_row = pq->row[nd/2];
     assert( kid_row != pq->numrows); 
     assert( par_row != pq->numrows); 
     pq->row[nd] = par_row;
     pq->row[nd/2] = kid_row;
     
     pq->node[kid_row] = nd/2;
     pq->node[par_row] = nd;
     
     nd = nd/2;
  }    
  return;
}

/*! Brief remove the root node of the heap 
 *  Replace it with the last node in the heap (making the heap one node shorter).
 *  Then restore the heap property by bubbling the root node down until
 *  it is not bigger than either of its children.
 *  \param pq the priority queue
 *  \return 0 if the heap is empty, 1 otherwise
 */  
int delete_root(PQ_t * pq)
{
  double w;
  int64_t kid_row, par_row; 
  int64_t par_nd, kid_nd;

  pq->tail--;                 // now the index of the last active node in the heap
  if( pq->tail == 1)          // the heap is now empty
    return(0);
  pq->node[pq->row[1]] = 0;           // mark this row as resolved
  pq->row[1] = pq->row[ pq->tail];
  pq->wt[1] = pq->wt[ pq->tail];
  pq->wt[pq->tail] = 99.0;            // pq->tail is now available
  pq->row[pq->tail] = pq->numrows;
  //printf("delete root:\n");
  //print_queue(pq);

  // recursively, if the parent node is not less than both it's children, then swap parent with smaller child
  // note: if right child is not in the heap (ie. pq->tail is the right child), 
  // then that nodes exists (even though it is not active) and its wt is INFINITY.  
  // We don't have to determine and make a special case if there is only a left child.
  par_nd = 1;
  //printf("bubble down: %lg : %lg, %lg\n", pq->wt[par_nd], pq->wt[2*par_nd], pq->wt[2*par_nd + 1]);
  while( (2*par_nd + 1 <= pq->tail) && 
         ((pq->wt[ 2*par_nd ] <  pq->wt[par_nd]) || (pq->wt[ 2*par_nd + 1 ] <  pq->wt[par_nd]) ) ) {
    if( pq->wt[2*par_nd] < pq->wt[2*par_nd + 1] ){
      kid_nd = 2*par_nd;
    } else {
      kid_nd = 2*par_nd + 1;
    }
    //printf("bubble down: swap parent  %"PRId64" with kid %"PRId64"\n", par_nd, kid_nd);
    w = pq->wt[par_nd];
    pq->wt[par_nd] = pq->wt[kid_nd];
    pq->wt[kid_nd] = w;
    
    par_row = pq->row[par_nd];
    kid_row = pq->row[kid_nd];
    pq->row[par_nd] = kid_row;
    pq->row[kid_nd] = par_row;

    pq->node[par_row] = kid_nd;
    pq->node[kid_row] = par_nd;
   
    par_nd = kid_nd;
  }
  //print_queue(pq);

  return(1);
}

/*!
 * \brief The implementation of Dijkstra's algorithm that uses a heap to prioritize the tentative vertices.
 * \param *mat sparsemat_t that holds the graph. 
 * \param r0 is the starting row (vertex)
 * \return run time
 */
double sssp_dijsktra_heap(sparsemat_t * mat, double *dist, int64_t r0)
{
  int64_t i, k;
  int64_t numrows = mat->numrows;

  // initialize the heap 
  PQ_t * pq = init_pqueue(numrows);
  for(i=0; i<numrows+1; i++){
    pq->wt[i] = 99.0;
    pq->row[i] = numrows;
    pq->node[i] = numrows;
  }
  pq->wt[0] = 0.0;
  pq->row[0] = r0;
  pq->node[r0] = 0;

  dist[r0] = 0.0;
  //printf(">>>>>>>>>>>>> dist[%"PRId64"] = %lg\n", r0, dist[r0]);

  int64_t rn, row, nd;
  double e_wt, vn_wt;
  // start the heap with weights to vertices adjacent to r0
  pq->tail = 1;
  for(k = mat->offset[r0];  k < mat->offset[r0+1]; k++){
    row = mat->nonzero[k];
    e_wt = mat->value[k];
    nd = pq->node[row] = pq->tail;
    pq->wt[nd] = e_wt;
    pq->row[nd] = row;
    pq->tail++;
    //printf("explore (%2"PRId64",%2"PRId64"): new row %"PRId64" with wt %lg at node %"PRId64"\n", r0, row, row, pq->wt[nd], nd);
    bubble_up(pq, nd);
    //print_queue(pq);
  }
  //print_queue(pq);
  //printf("\n");

  while(1){
    rn     = pq->row[1];
    vn_wt  = pq->wt[1];
    dist[rn] = vn_wt;
    //printf(">>>>>>>>>>>>> dist[%"PRId64"] = %lg\n", rn, vn_wt);
    for(k = mat->offset[rn];  k < mat->offset[rn+1]; k++){
      row = mat->nonzero[k];
      e_wt  = mat->value[k];
      nd   = pq->node[row];
      if(nd == 0){                  // row is done
        //printf("explore (%2"PRId64",%2"PRId64"): done with row (%"PRId64")\n", rn, row, row);
        continue;
      }
      if( nd == numrows ) {                      // row is new
        nd = pq->tail;
        pq->node[row] = nd;
        pq->wt[nd] = vn_wt + e_wt;
        pq->row[nd] = row;
        pq->tail++;
        //printf("explore (%2"PRId64",%2"PRId64"): new row %"PRId64" with wt %lg at node %"PRId64"\n", rn, row, row, pq->wt[nd], nd);
        bubble_up(pq, nd);
        //print_queue(pq);
      } else if(pq->wt[nd] > (vn_wt + e_wt)) {          // improved the weight
        pq->wt[nd] = vn_wt + e_wt;
        //printf("explore (%2"PRId64",%2"PRId64"): updated row %"PRId64" with wt %lg at node %"PRId64"\n", rn, row, row, pq->wt[nd], nd);
        bubble_up(pq, nd);
        //print_queue(pq);
      } else {
        //printf("explore (%2"PRId64",%2"PRId64"): no change\n", rn, row);
      }
    }
    if(delete_root(pq) == 0)                           // nothing else in the heap
      break;
  };

  return(0.0);
}


/*
 * These functions were useful during the development of the heap based implementation.
 * We left them here for debugging and exploring the code.
 */

// start at the top and bubble down
// going to assume that a bunch of the node are  out of order
void heapify_pqueue(PQ_t * pq)
{
  int64_t k;

  //printf("heapify:\n");
  for(k=2;  k < pq->tail; k++) {
    bubble_up(pq, k);
  }
  //print_queue(pq);
}


// check that the wt of a parent is never greater than either child
int64_t check_pqueue(PQ_t * pq)
{
  int64_t p, ret_ok=1;
  for(p=1; p <= (pq->tail)/2; p++){
    if(  ((2*p) < pq->tail && pq->wt[p] > pq->wt[(2*p)]) 
       ||((2*p+1) < pq->tail &&  pq->wt[p] > pq->wt[(2*p+1)] ) ) {
      ret_ok = 0;
      break;
    }
  }
  return(ret_ok);
}

// just dump the entries in the heap it three consecutive lines.
// numrows needs to be small
void print_queue(PQ_t * pq)
{
  int i;
  printf("wt:   "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3lg ", pq->wt[i]);
  }
  printf("\nrow:  "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3"PRId64" ", pq->row[i]);
  }
  printf("\nnode: ");
  for(i=0; i<pq->numrows+1; i++){
    printf("%3"PRId64" ", pq->node[i]);
  }
  printf("\n\n");
}
  
