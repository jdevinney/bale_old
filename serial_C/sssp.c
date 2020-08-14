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

/*! \file sssp.c
 * \brief Demo application that implements different 
 * Single Source Shortest Path alogrithms.
 *
 * versions of
 * Dijsktra's Algorithm
 * the Delta-stepping algorithm for 
 */

#include "spmat_utils.h"

double sssp_dijsktra_linear(sparsemat_t * mat, double *dist, int64_t v0);
double sssp_dijsktra_heap(sparsemat_t * mat, double *dist, int64_t r0);
#if 0
/*!
 * \brief This routine implements the most naive version of Dijkstra's algorithm
 * \param *mat sparsemat_t that holds the graph. 
 *   We assume:
 *     -- that the graph is weighted with positive weights
 *     -- if the graph is undirected, then mat contains two "directed edges" for
 *        each undirected edge.
 * \param v0 is the starting vertex
 * \return run time
 */
double sssp_dijsktra_linear(sparsemat_t * mat, double *dist, int64_t v0)
{
	int64_t i, k, vn;
	int64_t numrows = mat->numrows;
	double minwt, curwt;
	int64_t minidx;

	for(i=0; i<numrows; i++)
		dist[i] = INFINITY;
	for(k = mat->offset[v0];  k < mat->offset[v0+1]; k++)
		dist[ mat->nonzero[k] ] = mat->value[k];
  dist[v0] = -0.0;
	printf("dist[v0] = %lg\n", dist[v0]);

	while( 1 ) {
		// find the smallest tentative distance
	  minwt = INFINITY;
		minidx = numrows;
	  for(i=0; i<numrows; i++){
      if( dist[i] > 0 && dist[i] < minwt ){
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
			vn = mat->nonzero[k];
      if( dist[vn] > curwt + mat->value[k] )
				dist[vn] = curwt + mat->value[k];
		}
		//printf("dist[%ld] = %lg\n", minidx, curwt);
    dist[minidx] = -dist[minidx];
  }

	for(i=0; i<numrows; i++)
		dist[i] = -dist[i];

	return(0.0);
}


typedef struct PQ_t {
  int64_t numrows; // the first node not in the queue;
  int64_t tail;    // the first node not in the queue;
  double  *val;    // the corresponding value that is being prioritized
  int64_t *row;    // the index of the row being represented by a node in the queue
  int64_t *node;   // the index into the queue where the given row is stored
}PQ_t;

void  heapify_pqueue(PQ_t * pq);
void print_queue(PQ_t * pq);
int64_t check_pqueue(PQ_t * pq);

PQ_t * init_pqueue(int64_t numrows);
void bubble_up(PQ_t *pq, int64_t k);
int delete_root(PQ_t * pq);

// Note the queue will be "1-up" so the parent and children are easy to calculate.
// We never put the root vertex into the queue so there will be lot of space.
PQ_t * init_pqueue(int64_t numrows)
{
  PQ_t * pq = calloc(1, sizeof(PQ_t));
  pq->numrows = numrows;
  pq->tail = 1;
  pq->val  = (double  *)calloc(numrows+1, sizeof(double));
  pq->row  = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  pq->node = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  return(pq);
}

// bubble up the kth node
void bubble_up(PQ_t *pq, int64_t k)
{
   double v;
   int64_t kid_row, par_row; 
   while( k > 1 ){
     if( pq->val[ k/2 ] <= pq->val[k] )
       return; 
     printf("swap %ld and %ld\n", k, k/2);
     //swap the kid node with the parent node
     v = pq->val[k];
     pq->val[k] = pq->val[k/2];
     pq->val[k/2] = v;
     
     kid_row = pq->row[k];
     par_row = pq->row[k/2];
     assert( kid_row != pq->numrows); 
     assert( par_row != pq->numrows); 
     pq->row[k] = par_row;
     pq->row[k/2] = kid_row;
     
     pq->node[kid_row] = k/2;
     pq->node[par_row] = k;
     
     k = k/2;
  }    
  return;
}

// delete the root node of the queue 
int delete_root(PQ_t * pq)
{
  double v;
  int64_t kid_row, par_row; 
  int64_t par_nd, kid_nd;

  pq->tail--;
  if( pq->tail == 1)                   // the heap is now empty
    return(0);
  pq->node[pq->row[1]] = 0;            // mark this row as done
  pq->row[1] = pq->row[ pq->tail];
  pq->val[1] = pq->val[ pq->tail];
  pq->val[pq->tail] = 99.0;
  pq->row[pq->tail] = pq->numrows;
  printf("delete root:\n");
  print_queue(pq);
  // should be bubble down from the top
  //heapify_pqueue(pq);
  //bubble_down(pq);

  par_nd = 1;
  // recursively, if par node is not less than both it's children, then swap parent with smaller child
  // note: if right child is not in the heap (ie. pq->tail is the right child), then that nodes exists 
  // and its val is INFINITY.  We don't have worry if there is only a left child.
  printf("bubble down: %lg : %lg, %lg\n", pq->val[par_nd], pq->val[2*par_nd], pq->val[2*par_nd + 1]);
  while( (2*par_nd + 1 <= pq->tail) && 
         ((pq->val[ 2*par_nd ] <  pq->val[par_nd]) || (pq->val[ 2*par_nd + 1 ] <  pq->val[par_nd]) ) ) {
    if( pq->val[2*par_nd] < pq->val[2*par_nd + 1] ){
      kid_nd = 2*par_nd;
    } else {
      kid_nd = 2*par_nd + 1;
    }
    //swap the kid node with the parent node
    printf("bubble down: swap %ld with %ld\n", par_nd, kid_nd);
    v = pq->val[par_nd];
    pq->val[par_nd] = pq->val[kid_nd];
    pq->val[kid_nd] = v;
    
    par_row = pq->row[par_nd];
    kid_row = pq->row[kid_nd];
    pq->row[par_nd] = kid_row;
    pq->row[kid_nd] = par_row;

    pq->node[par_row] = kid_nd;
    pq->node[kid_row] = par_nd;
   
    par_nd = kid_nd;
  }
  print_queue(pq);

  return(1);
}

// start at the top and bubble down
// going to assume that a bunch of sort is out of order
void heapify_pqueue(PQ_t * pq)
{
  int64_t k;

  printf("heapify:\n");
  for(k=2;  k < pq->tail; k++) {
    bubble_up(pq, k);
  }
  print_queue(pq);
}


// check that the val of a parent is always never greater than either child
int64_t check_pqueue(PQ_t * pq)
{
  int64_t p, ret_ok=1;
  for(p=1; p <= (pq->tail)/2; p++){
    if(  ((2*p) < pq->tail && pq->val[p] > pq->val[(2*p)]) 
       ||((2*p+1) < pq->tail &&  pq->val[p] > pq->val[(2*p+1)] ) ) {
      ret_ok = 0;
      break;
    }
  }
  return(ret_ok);
}

void print_queue(PQ_t * pq)
{
  int i;
  printf("val:  "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3lg ", pq->val[i]);
  }
  printf("\nrow:  "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3ld ", pq->row[i]);
  }
  printf("\nnode: ");
  for(i=0; i<pq->numrows+1; i++){
    printf("%3ld ", pq->node[i]);
  }
  printf("\n\n");
}
  
/*!
 * \brief This routine implements generic serial version of Dijkstra's algorithm
 * \param *mat sparsemat_t that holds the graph. 
 *   We assume:
 *     -- that the graph is weighted 
 *     -- if the graph is undirected, then mat contains two "directed edges" for
 *        each undirected edge.
 * \param r0 is the starting row (vertex)
 * \return run time
 * The algorithm requires that the weights are positive.
 */
double sssp_dijsktra_heap(sparsemat_t * mat, double *dist, int64_t r0)
{
	int64_t i, k;
	int64_t numrows = mat->numrows;

  PQ_t * pq = init_pqueue(numrows);

  for(i=1; i<numrows; i++) 
    pq->val[i] = 99.0;

  int64_t pos=1;
  for(i=0; i<numrows+1; i++){
    if( i == r0){
      pq->val[0] = 0.0;
      pq->row[0] = r0;
      pq->node[r0] = 0;
      continue;
    } 
    pq->val[pos] = 99.0;
    pq->row[pos] = numrows;
    pq->node[i] = numrows;
    pos++;
  }
  print_queue(pq);
  printf("\n");

  int64_t vn, vert, nd;
  double val, vnval;
  int64_t root = 0;
  pq->tail = 1;

  while(1){
    vn     = pq->row[root];
    vnval  = pq->val[root];
    dist[vn] = vnval;
    printf(">>>>>>>>>>>>> dist[%ld] = %lg\n", vn, vnval);
    for(k = mat->offset[vn];  k < mat->offset[vn+1]; k++){
      vert = mat->nonzero[k];
      val  = mat->value[k];
      nd   = pq->node[vert];
      if(nd == 0){                  // row is done
        printf("explore (%2ld,%2ld): done with vert (%ld)\n", vn, vert, vert);
        continue;
      }
      if( nd == numrows ) {        // row is new
        pq->node[vert] = pq->tail;
        nd = pq->node[vert];
        pq->val[nd] = vnval + val;
        pq->row[nd] = vert;
        pq->tail++;
        printf("explore (%2ld,%2ld): new vert %ld with val %lg at node %ld\n", vn, vert, vert, pq->val[nd], nd);
        bubble_up(pq, nd);
        print_queue(pq);
      } else if( pq->val[nd] > (vnval + val)) {
        pq->val[nd] = vnval + val;
        pq->row[nd] = vert;
        printf("explore (%2ld,%2ld): updated vert %ld with val %lg at node %ld\n", vn, vert, vert, pq->val[nd], nd);
        bubble_up(pq, nd);
        print_queue(pq);
      } else {
        printf("explore (%2ld,%2ld): no change\n", vn, vert);
      }
    }

    if( root == 0 ) {
      root = 1;
      continue;
    }
    if(delete_root(pq) == 0 )
      break;
  };

	return(0.0);
}
#endif

int main(int argc, char * argv[]) 
{
  #define NUMROWS 20 
  int64_t numrows=NUMROWS;
  double edge_prob = 0.25;
  uint32_t seed = 123456789;
  int64_t readgraph = 0;
  char filename[256]={"filename"};
  enum MODEL {GENERIC_Model=1, DIJSKTRA_HEAP=2, ALL_Models=4};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;
  int64_t i;

  sparsemat_t *mat, *dmat;
 
  int64_t dump_files = 1;
 
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:e:M:f:Dq")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows ); break;
    case 's': sscanf(optarg,"%d" ,&seed ); break;
    case 'e': sscanf(optarg,"%lg", &edge_prob); break;
    case 'M': sscanf(optarg,"%d" , &models_mask); break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'D': dump_files = 1; break;
    case 'q': quiet = 1; break;
    default: break;
    }
  }
 
#if 0 // test the heap stuff
    PQ_t * pq = init_pqueue(numrows);
    for(i=1; i<numrows; i++) {
      pq->val[i] = (double)(numrows-i);
      pq->row[i] = i;
      pq->node[i] = i;
    }
    pq->tail = numrows;
    print_queue(pq);

    heapify_pqueue(pq);
    exit(1);
#endif



  if( readgraph ) {
    mat = read_matrix_mm(filename);
    if(!mat){printf("ERROR: sssp: read graph from %s Failed\n", filename); exit(1);}
  } else {
    //mat = erdos_renyi_tri(numrows, er_prob, ER_TRI_L, seed);
		//graph_model model = FLAT;
    mat = random_graph(numrows, FLAT, UNDIRECTED_WEIGHTED, NOLOOPS, edge_prob, seed);
    if(!mat){printf("ERROR: triangles: erdos_renyi_graph Failed\n"); exit(1);}
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of sssp\n");
    fprintf(stderr,"Number of rows    (-n) %ld\n", numrows);
    fprintf(stderr,"random seed     (-s)= %d\n", seed);
    fprintf(stderr,"erdos_renyi_prob   (-e)= %lg\n", edge_prob);
    fprintf(stderr,"models_mask     (-M)= %d\n", models_mask);
    fprintf(stderr,"readgraph      (-f [%s])\n", filename); 
    fprintf(stderr,"dump_files      (-D)=%ld\n", dump_files);
    fprintf(stderr,"quiet        (-q)= %d\n", quiet);
 
    if(printhelp)
      return(0);
  }

  for(i=0; i< mat->nnz; i++) {
    mat->value[i] = floor( 20 * mat->value[i]);
  }
	// mat is the lower triangle of the adjacency matrix
  // we use the full symmetric matrix in Dijsktra's algorithm
  dmat = make_symmetric_from_lower(mat);
	// debug info
  if(dump_files) {
		dump_matrix(mat,20, "L.out");
		dump_matrix(dmat,20, "Full.out");
	}

	double * dist = malloc(numrows*sizeof(double));
  for(i=0;i<numrows; i++) dist[i] = INFINITY;

  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic          sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_dijsktra_linear(dmat, dist, 5);
	    for(i=0; i<numrows; i++)
		    printf("%ld %g\n",i, dist[i]);
      break;
    case DIJSKTRA_HEAP:
      if( !quiet ) printf("Dijsktra Heap    sssp:\n");
      for(i=0;i<numrows; i++) dist[i] = INFINITY;
      laptime = sssp_dijsktra_heap(dmat, dist, 5);
	    for(i=0; i<numrows; i++)
		    printf("%ld %g\n",i, dist[i]);
      break;
    default:
      continue;
    }
  }

	printf("done: %g\n", laptime);
  clear_matrix(mat);
  clear_matrix(dmat);
  return(0);
}
