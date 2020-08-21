#include "spmat_utils.h"


typedef struct llnode_t{
  int64_t index;
  struct llnode_t * next;
  struct llnode_t * prev;
}llnode_t;

typedef struct request_t{
  int64_t w;
  double x;
}request_t;

#if 0
request_t find_requests(int64_t bin, int light, llnode_t ** bins, sparsemat_t * mat){
  int64_t j;
  llnode_t * node;
  while(node != NULL){
    int64_t v = node->index;
    for(j = mat->offset[v]; j < mat->offset[v+1]; j++){
      if
    }
    node = node->next;
  }
  
}
#endif

double relax_edge(sparsemat_t * mat, int64_t u, int64_t v, double * tent){

  /* find weight of edge(u, v) */
  double c_uv = -1;
	int64_t i;
  for(i = mat->offset[u]; i < mat->offset[u+1]; i++)
    if(mat->nonzero[i] == v)
      c_uv = mat->value[i];

  if(c_uv == -1){
    printf("ERROR: relax edge: illegal u,v pair (%"PRId64" %"PRId64")\n", u, v);
    return (1);
  }

  if( tent[v] < (tent[u] + c_uv)){
    tent[v] = tent[u] + c_uv;
  }
  return(0);
}


void remove_node_from_bucket(llnode_t * v, int64_t i, llnode_t ** B){
  llnode_t * next = v->next;
  printf("Removing %"PRId64" from bucket %"PRId64"\n", v->index, i);
  if(v->prev)
    v->prev->next = next;
  else
    B[i] = v->next;
  if(next != NULL)
    next->prev = v->prev;
}



void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, llnode_t * nodes, llnode_t ** B){
  if(i == j) return;
  llnode_t * node = B[i];
  
  /* remove w from B[i] if it is there */  
  while(node != NULL){
    if(node == w){      
      remove_node_from_bucket(w, i, B);
      break;
    }    
    node = node->next;
  }
  
  /* insert w into the front of B[j] */
  printf("Adding %"PRId64" to bucket %"PRId64"\n", w->index, j);
  w->next = NULL;
  w->prev = NULL;
  if(B[j]){
    w->next = B[j];
    B[j]->prev = w;
  }
  B[j]= w;
}


void relax(int64_t windex, double x, double delta, double * tent, llnode_t * nodes, llnode_t ** B){
  printf("relax w=%"PRId64" x = %lf < %lf?\n", windex, x, tent[windex]);
  if ( x < tent[windex] ){
    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
    move_node_from_bucket_i_to_j(&nodes[windex], (int64_t)floor(tent[windex]/delta), (int64_t)floor(x/delta), nodes, B);
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
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta);
  printf("num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", num_buckets, delta, max_degree);
  llnode_t ** B = calloc(num_buckets, sizeof(llnode_t *));
  for(i = 0; i < num_buckets; i++)
    B[i] = NULL;

  llnode_t * nodes = calloc(mat->numrows, sizeof(llnode_t));
  for(i = 0; i < mat->numrows; i++){
    nodes[i].index = i;
    nodes[i].next = NULL;
    nodes[i].prev = NULL;
  }
  
  /* set the source distance to 0 */
  tent[r0] = 0;
  B[0] = &nodes[r0];
    
  /* main loop */
  int64_t min_bucket = 0;
  int64_t num_deleted = 0;
  while(1){

    /* find the minimum indexed non-empty bucket */
    for(min_bucket = 0; min_bucket < num_buckets; min_bucket++)
      if(B[min_bucket])
        break;
    printf("Starting inner loop: working on bucket %"PRId64"\n", min_bucket);
    if(min_bucket == num_buckets)
      break;
    
    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    llnode_t * v = B[min_bucket];
    
    while(v != NULL){
      printf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v->index, min_bucket);

      /* take v out of B[min_bucket]??? */
      remove_node_from_bucket(v, min_bucket, B);
      
      /* find light requests for v */
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] <= delta){
          relax(mat->nonzero[j], tent[v->index] + mat->value[j], delta, tent, nodes, B);
        }
      }      
      
      /* insert v into R if it is not already there */
      if(deleted[v->index] == 0){
        deleted[v->index] = 1;
        R[end++] = v->index;
        num_deleted++;
        printf("deleted %"PRId64"s\n", v->index);
      }
      
      v = v->next;
    }// end inner loop

    /* find heavy requests for everything in R */
    while(start < end){
      v = &nodes[R[start++]];
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] > delta){
          relax(mat->nonzero[j], tent[v->index] + mat->value[j], delta, tent, nodes, B);
        }
      }      
    }
    
  }// end main loop

  
  
  return(0.0);
}
