#include "spmat_utils.h"


typedef struct llnode_t{
  int64_t index;
  struct llnode_t * next;
  struct llnode_t * prev;
}llnode_t;

typedef struct buckets_t{
  llnode_t ** B;
  llnode_t * nodes;
  int64_t * level;
  char * empty;
  char * in_bucket;
  int64_t num_buckets;
}buckets_t;

void remove_node_from_bucket(llnode_t * v, int64_t i, buckets_t * buckets){
  llnode_t * next = v->next;
  assert(buckets->empty[i] == 0);
  //printf("Removing %"PRId64" from bucket %"PRId64"\n", v->index, i);
  if(v->prev)
    v->prev->next = next;
  else
    buckets->B[i] = v->next;
  if(next != NULL)
    next->prev = v->prev;
  if(buckets->B[i] == NULL)
    buckets->empty[i] = 1;
  buckets->in_bucket[v->index] = 0;
}

void insert_node_in_bucket(llnode_t * w, int64_t i, buckets_t * buckets){
  int64_t actual_i = (i % buckets->num_buckets);
  //printf("Adding %"PRId64" to bucket %"PRId64"\n", w->index, actual_i);
  if(buckets->empty[actual_i] == 1){
    assert(buckets->B[actual_i] == NULL);
    buckets->empty[actual_i] = 0;
    buckets->level[actual_i] = (actual_i % buckets->num_buckets);
  }else{
    //printf("%lld %lld\n", buckets->level[actual_i], (actual_i % buckets->num_buckets));
    assert(buckets->level[actual_i] == (actual_i % buckets->num_buckets));
  }
  w->next = NULL;
  w->prev = NULL;
  if(buckets->B[actual_i]){
    w->next = buckets->B[actual_i];
    buckets->B[actual_i]->prev = w;
  }
  buckets->B[actual_i]= w;
  buckets->in_bucket[w->index] = 1;
}

void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets){


  if(i >= 0 && buckets->in_bucket[w->index]){    
    //printf("Move from %lld to %lld\n",i, j);
    i = i % buckets->num_buckets;
    llnode_t * node = buckets->B[i];
    
    /* remove w from B[i] if it is there */
    while(node != NULL){
      if(node == w){
	remove_node_from_bucket(w, i, buckets);
	break;
      }
      node = node->next;
    }
  }
  
  /* insert w into the front of B[j] */
  insert_node_in_bucket(w, j, buckets);

}


void relax(int64_t windex, double x, double delta, double * tent, buckets_t * buckets){
  if ( x < tent[windex] ){
    //printf("relax w=%"PRId64" x = %lf < %lf?\n", windex, x, tent[windex]);
    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
    int64_t j;
    if(tent[windex] == INFINITY) j = -1;
    else j = (int64_t)floor(tent[windex]/delta);
    move_node_from_bucket_i_to_j(&(buckets->nodes[windex]), j, (int64_t)floor(x/delta), buckets);
    tent[windex] = x;
  }
}


double sssp_delta_stepping(sparsemat_t * mat, double * dist, int64_t r0){
  int64_t i, j;

  assert(r0 < mat->numrows);
  
  char * deleted = calloc(mat->numrows, sizeof(char));
  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  double * tent = dist;
  for(i = 0; i < mat->numrows; i++){
    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
      max_degree = (mat->offset[i+1] - mat->offset[i]);
    tent[i] = INFINITY;
  }
  
  double max_edge_weight = 0;
  for(i = 0; i < mat->nnz; i++)
    if(max_edge_weight < mat->value[i])
      max_edge_weight = mat->value[i];
  printf("max edge weight = %lf\n", max_edge_weight);
  assert(max_degree > 0.0);
  delta = 1.0/max_degree;
  
  /* set up buckets as an array of linked lists */
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  buckets_t * buckets = calloc(1, sizeof(buckets_t));
  buckets->num_buckets = num_buckets;
  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
  buckets->nodes = calloc(mat->numrows, sizeof(llnode_t));
  buckets->in_bucket = calloc(mat->numrows, sizeof(char));
  buckets->level = calloc(num_buckets, sizeof(int64_t));
  buckets->empty = calloc(num_buckets, sizeof(char));

  // TODO: need to be able to tell if a node is in a bucket or not!
  
  printf("num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", num_buckets, delta, max_degree);
  for(i = 0; i < num_buckets; i++){
    buckets->B[i] = NULL;
    buckets->empty[i] = 1;
    buckets->level[i] = -1;
  }
  for(i = 0; i < mat->numrows; i++){
    buckets->nodes[i].index = i;
    buckets->nodes[i].next = NULL;
    buckets->nodes[i].prev = NULL;
  }
  
  /* set the source distance to 0 */
  insert_node_in_bucket(&(buckets->nodes[r0]), 0, buckets);
  tent[r0] = 0;
  
  /* main loop */
  int64_t min_bucket = 0;
  int64_t num_deleted = 0;
  while(1){

    /* find the minimum indexed non-empty bucket */
    for(min_bucket = 0; min_bucket < buckets->num_buckets; min_bucket++)
      if(buckets->empty[min_bucket] == 0)
        break;
    //printf("Starting inner loop: working on bucket %"PRId64"\n", min_bucket);
    if(min_bucket == num_buckets)
      break;
    
    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    llnode_t * v = buckets->B[min_bucket];
    
    while(v != NULL){
      //printf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v->index, min_bucket);

      /* take v out of B[min_bucket]??? */
      remove_node_from_bucket(v, min_bucket, buckets);
      
      /* find light requests for v */
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] <= delta){
          relax(mat->nonzero[j], tent[v->index] + mat->value[j], delta, tent, buckets);
        }
      }      
      
      /* insert v into R if it is not already there */
      if(deleted[v->index] == 0){
        deleted[v->index] = 1;
        R[end++] = v->index;
        num_deleted++;
        //printf("deleted %"PRId64"s\n", v->index);
      }
      
      v = v->next;
    }// end inner loop

    /* find heavy requests for everything in R */
    while(start < end){
      v = &(buckets->nodes[R[start++]]);
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] > delta){
          relax(mat->nonzero[j], tent[v->index] + mat->value[j], delta, tent, buckets);
        }
      }      
    }
    
  }// end main loop
  free(buckets->B);
  free(buckets->nodes);
  free(buckets->empty);
  free(buckets->level);
  free(deleted);
  free(R);
  
  return(0.0);
}
