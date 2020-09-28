/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

#include "spmat_utils.h"

typedef struct llnode_t{
  int64_t index;
  struct llnode_t * next;
  struct llnode_t * prev;
}llnode_t;

typedef struct buckets_t{
  llnode_t ** B;
  llnode_t * nodes;
  int64_t * in_bucket;
  int64_t num_buckets;
  double delta;
}buckets_t;


// Remove a specific node from bucket i.
void remove_node_from_bucket_ptr(llnode_t * v, int64_t i, buckets_t * buckets)
{
  i = i % buckets->num_buckets;
  assert(buckets->in_bucket[v->index] == i);
  buckets->in_bucket[v->index] = -1;
  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v->index, i);
  if(v->prev){
    v->prev->next = v->next;
  }else
    buckets->B[i] = v->next;
  if(v->next != NULL)
    v->next->prev = v->prev;
}

// Prepend a node into a bucket
void insert_node_in_bucket_ptr(llnode_t * w, int64_t i, buckets_t * buckets)
{
  int64_t actual_i = (i % buckets->num_buckets);

  Dprintf("Adding %"PRId64" to bucket %"PRId64"\n", w->index, actual_i);
  
  if(buckets->in_bucket[w->index] == actual_i)  // this node is already in this bucket
    return; 

  assert(buckets->in_bucket[w->index] == -1);

  w->next = NULL;
  w->prev = NULL;
  if(buckets->B[actual_i]){
    w->next = buckets->B[actual_i];
    buckets->B[actual_i]->prev = w;
  }
  buckets->B[actual_i]= w;
  buckets->in_bucket[w->index] = actual_i;
}


// Remove a node from bucket i and put it into bucket j
void move_node_from_bucket_i_to_j_ptr(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets)
{
  llnode_t * node;
  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
          w->index, buckets->in_bucket[w->index],i, j);
  if(i >= 0){
    i = i % buckets->num_buckets;
    if(buckets->in_bucket[w->index] == i){
      node = buckets->B[i]; 
      while(node != NULL){ // find w in B[i] if it is there
	      if(node == w){
	        remove_node_from_bucket_ptr(w, i, buckets);
	        break;
	      }
	      node = node->next;
      }
    }
  }
  insert_node_in_bucket_ptr(w, j, buckets);   // insert w into the front of B[j] 
}


// relax an edge from a node to a node w (the current tenative distance to w is tent[windex]).
// the candidate distance to w is cand_dist.
void relax_ptr(int64_t windex, double cand_dist, double * tent, buckets_t * buckets)
{
  if ( cand_dist < tent[windex] ){
    Dprintf("relax w=%"PRId64" cand_dist = %lf < %lf?\n", windex, cand_dist, tent[windex]);
    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
    int64_t j;
    if(tent[windex] == INFINITY) 
      j = -1;
    else 
      j = (int64_t)floor(tent[windex]/buckets->delta);
    move_node_from_bucket_i_to_j_ptr(&(buckets->nodes[windex]), j, (int64_t)floor(cand_dist/buckets->delta), buckets);
    tent[windex] = cand_dist;
  }
}

// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_stepping_ptr(d_array_t *dist, sparsemat_t * mat, int64_t r0, double del){
  int64_t i, j;

  double tm = wall_seconds();
  
  assert(r0 < mat->numrows);
  
  char * deleted = calloc(mat->numrows, sizeof(char));
  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  double * tent = dist->entry;
  for(i = 0; i < mat->numrows; i++){
    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
      max_degree = (mat->offset[i+1] - mat->offset[i]);
    tent[i] = INFINITY;
  }
  assert(max_degree > 0);
  delta = 1.0/max_degree;
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->nnz; i++)
    if(max_edge_weight < mat->value[i])
      max_edge_weight = mat->value[i];
  Dprintf("max edge weight = %lf\n", max_edge_weight);
  
  /* set up buckets as an array of linked lists */
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  buckets_t * buckets = calloc(1, sizeof(buckets_t));
  buckets->num_buckets = num_buckets;
  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
  buckets->nodes = calloc(mat->numrows, sizeof(llnode_t));
  buckets->in_bucket = calloc(mat->numrows, sizeof(int64_t));
  buckets->delta = delta;
  
  Dprintf("num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", num_buckets, delta, max_degree);
  for(i = 0; i < num_buckets; i++){
    buckets->B[i] = NULL;
  }
  for(i = 0; i < mat->numrows; i++){
    buckets->nodes[i].index = i;
    buckets->nodes[i].next = NULL;
    buckets->nodes[i].prev = NULL;
    buckets->in_bucket[i] = -1;
  }
  
  /* set the source distance to 0 */
  insert_node_in_bucket_ptr(&(buckets->nodes[r0]), 0, buckets);
  tent[r0] = 0.0;
  
  /* main loop */
  int64_t current = 0;
  while(1){

    /* find the minimum indexed non-empty bucket */
    for(i = 0; i < buckets->num_buckets; i++)
      if(buckets->B[(current + i) % buckets->num_buckets] != NULL)
        break;
    
    if(i == num_buckets)
      break;
    current = (current + i) % num_buckets;

    Dprintf("Starting inner loop: working on bucket %"PRId64"\n", current);
    int64_t start = 0;
    int64_t end = 0;
    llnode_t * v;
    // inner loop
    while(buckets->B[current] != NULL){
      v = buckets->B[current];
      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v->index, current);

      /* relax light edges from v */
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] <= delta){	  
          relax_ptr(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
        }
      } 
      
      /* insert v into R if it is not already there */
      if(deleted[v->index] == 0){
        deleted[v->index] = 1;
        R[end++] = v->index;
        Dprintf("deleted %"PRId64"s\n", v->index);
      }
      
      remove_node_from_bucket_ptr(v, current, buckets);
    }

    /* relax heavy requests edges for everything in R */
    while(start < end){
      v = &(buckets->nodes[R[start++]]);
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] > delta){
          relax_ptr(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
        }
      }      
    }
    current++;
  }// end main loop
  free(buckets->B);
  free(buckets->nodes);
  free(deleted);
  free(R);
  
  return(wall_seconds() - tm);
}

typedef struct ds_t{
  int64_t *next;      // next in linked list
  int64_t *prev;      // prev in linked list
  int64_t *in_bucket; // which bucket it is in, else -1
  int64_t *deleted;   // deleted means resolved?
  double  *tent;      // the tentative weight for the vertex
  int64_t *B;         // array of Buckets: B[i] is the index of the first node on the list, or -1 if empty
  int64_t num_buckets;
  double  delta;      // delta
}ds_t;

// debugging function
void dump_bucket_arr(ds_t *ds, int64_t i_m)
{
  int64_t v, w ;
  printf("Bucket[%ld] =", i_m);
  if( ds->B[i_m] == -1 ){
    printf("empty \n");
    return;
  }
  v = w = ds->B[i_m];
  do {
    printf("%ld ", w);
  } while( (w = ds->next[w]) != v );
  printf("\n");
  return;
}

// Prepend node v into a bucket i
void insert_node_in_bucket_arr(ds_t *ds, int64_t v, int64_t i_m)
{
  int64_t w;                   // node on list, given by ds->B[i_m]
  Dprintf("Adding %"PRId64" to bucket %"PRId64" of %"PRId64"\n", v, i_m, ds->num_buckets);
  
  assert(i_m >= -1 && i_m < ds->num_buckets);
  
  if(ds->in_bucket[v] == i_m){      // it is ok if this node is already in this bucket
    return; 
  }
  assert(ds->in_bucket[v] == -1);   // better not be in a different bucket

  ds->next[v] = v;
  ds->prev[v] = v;
  if(ds->B[i_m] != -1){    
    w = ds->B[i_m];             // w is "first" on the list, insert v before w              
    Dprintf("non-empty: w=%ld, prev=%ld\n", w, ds->prev[w]);
    ds->prev[v] = ds->prev[w];            
    ds->next[ds->prev[w]] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    Dprintf("   v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );
  }
  ds->B[i_m] = v;                       // set v to be the "first" on the list
  ds->in_bucket[v] = i_m;
}


// Remove node v from its bucket (need not be in any bucket)
static void remove_node_from_bucket_arr(ds_t *ds, int64_t v)
{
  int64_t i_m, w;

  i_m = ds->in_bucket[v];
  assert(i_m >= -1 && i_m < ds->num_buckets);
  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v, i_m);
  if(DEBUG && !( ((ds->next[v] == v) && (ds->prev[v] == v)) || (ds->next[v] != ds->prev[v]) )){
    printf("v, next, prev = %ld %ld %ld\n", v, ds->next[v], ds->prev[v] );
  }

  if(i_m == -1)     // v wasn't in a bucket
    return;

  ds->in_bucket[v] = -1;
  if((ds->next[v] == v) && (ds->prev[v] == v)){  // the only thing on the list
    ds->B[i_m] = -1;
    return;
  }

  w = ds->next[v];
  ds->prev[w] = ds->prev[v];
  ds->next[ds->prev[v]] = w;
  ds->B[i_m] = w;
  return;
}

// relax an edge to the head vertex, given the new tentative distance
// (= the tentative distance to the tail plus the weight of the edge).
// the candidate distance to w is cand_dist.
void relax_arr(ds_t *ds, int64_t w, double cand_dist)
{
  int64_t iold, inew;
  Dprintf("relax head %"PRId64" cand_dist = %lf < %lf?\n", w, cand_dist, ds->tent[w]);
  if ( cand_dist < ds->tent[w] ){
    iold = ds->in_bucket[w];
    inew = ((int64_t)floor(cand_dist/ds->delta)) % (ds->num_buckets);

    assert(iold >= -1 && iold < ds->num_buckets);
    assert(inew >= -1 && inew < ds->num_buckets);
    Dprintf("winner: %"PRId64"  move from bucket %"PRId64" to %"PRId64"\n", w, iold, inew);
    if( iold != inew ){
      if(iold >= 0)
        remove_node_from_bucket_arr(ds, w);
      insert_node_in_bucket_arr(ds, w, inew);
    }
    ds->tent[w] = cand_dist;
  }
}

// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_stepping_arr(d_array_t *dist, sparsemat_t * mat, int64_t r0, double del)
{
  int64_t i, i_m, k;
  int64_t v;
  int64_t rbi;     // the real bucket index

  double tm = wall_seconds();

  assert((r0 >= 0) && (r0<mat->numrows));
  
  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  for(i = 0; i < mat->numrows; i++){
    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
      max_degree = (mat->offset[i+1] - mat->offset[i]);
  }
  assert(max_degree > 0);
  delta = 1.0/max_degree;
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->nnz; i++)
    if(max_edge_weight < mat->value[i])
      max_edge_weight = mat->value[i];
  Dprintf("max edge weight = %lf\n", max_edge_weight);
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  
  // allocate and initialize the all the arrays
  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  ds->next = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->deleted = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  ds->tent = (double *)malloc(mat->numrows * sizeof(double)); assert(ds->tent != NULL);
  for(i = 0; i < mat->numrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->tent[i]      = INFINITY;
  }
  Dprintf("Allocate buckets\n");
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t)); assert(ds->B != NULL);
  // bucket indices are mod num_buckets.
  for(i_m = 0; i_m < num_buckets; i_m++){
    ds->B[i_m] = -1;
  }
  ds->delta = delta;

  // set the distance to r0 equal to 0.0
  Dprintf("Set source node %ld\n", r0);
  ds->tent[r0] = 0.0;
  insert_node_in_bucket_arr(ds, r0, 0);
  Dprintf("putting source r0=%ld in bucket %ld\n", r0, 0L);
  
  /* main loop */
  rbi = 0;
  while(1){
    Dprintf("Outer loop rbi = %"PRId64"\n", rbi);
    // find the minimum indexed non-empty bucket 
    for(i = 0; i < ds->num_buckets; i++){
      if(ds->B[(rbi + i) % ds->num_buckets] != -1)
        break;
    }
    if(i == num_buckets)  //All buckets are empty, we are done
      break;
    rbi = rbi + i;
    i_m = rbi % ds->num_buckets;
    Dprintf("Starting inner loop: working on bucket %"PRId64"\n", i_m);
    if(DEBUG) dump_bucket_arr(ds, i_m);

    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    for( v=ds->B[i_m]; ds->in_bucket[v] == i_m; v=ds->next[v]){
      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v, i_m);

      remove_node_from_bucket_arr(ds, v);
      
      /* relax light edges from v */
      for(k = mat->offset[v]; k < mat->offset[v + 1]; k++){
        if(mat->value[k] <= delta){	  
          relax_arr(ds, mat->nonzero[k], ds->tent[v] + mat->value[k]);
        }
      } 
      
      /* insert v into R if it is not already there */
      if(ds->deleted[v] == 0){
        ds->deleted[v] = 1;
        R[end++] = v;
        Dprintf("deleted %"PRId64"\n", v);
      }
    }

    /* relax heavy requests edges for everything in R */
    for(start=0; start<end; start++){
      v = R[start];
      for(k = mat->offset[v]; k < mat->offset[v + 1]; k++){
        if(mat->value[k] > delta){
          relax_arr(ds, mat->nonzero[k], ds->tent[v] + mat->value[k]);
        }
      }      
    }
  }

  // Copy the answers to dist array
  for(v = 0; v < mat->numrows; v++){
    dist->entry[v] = ds->tent[v];
  }

  free(ds->next);
  free(ds->prev);
  free(ds->in_bucket);
  free(ds->deleted);
  free(ds->tent);
  free(ds);
  free(R);
  
  return(wall_seconds() - tm);
}
