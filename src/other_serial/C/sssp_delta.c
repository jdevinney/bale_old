#include "spmat_utils.h"
#if 0
/**/typedef struct llnode_t{
/**/  int64_t index;
/**/  struct llnode_t * next;
/**/  struct llnode_t * prev;
/**/}llnode_t;
/**/
/**/typedef struct buckets_t{
/**/  llnode_t ** B;
/**/  llnode_t * nodes;
/**/  int64_t * in_bucket;
/**/  int64_t num_buckets;
/**/  double delta;
/**/}buckets_t;
/**/
/**/
/**/// Remove a specific node from bucket i.
/**/void remove_node_from_bucket(llnode_t * v, int64_t i, buckets_t * buckets)
/**/{
/**/  i = i % buckets->num_buckets;
/**/  assert(buckets->in_bucket[v->index] == i);
/**/  buckets->in_bucket[v->index] = -1;
/**/  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v->index, i);
/**/  if(v->prev){
/**/    v->prev->next = v->next;
/**/  }else
/**/    buckets->B[i] = v->next;
/**/  if(v->next != NULL)
/**/    v->next->prev = v->prev;
/**/}
/**/
/**/// Prepend a node into a bucket
/**/void insert_node_in_bucket(llnode_t * w, int64_t i, buckets_t * buckets)
/**/{
/**/  int64_t actual_i = (i % buckets->num_buckets);
/**/
/**/  Dprintf("Adding %"PRId64" to bucket %"PRId64"\n", w->index, actual_i);
/**/  
/**/  if(buckets->in_bucket[w->index] == actual_i)  // this node is already in this bucket
/**/    return; 
/**/
/**/  assert(buckets->in_bucket[w->index] == -1);
/**/
/**/  w->next = NULL;
/**/  w->prev = NULL;
/**/  if(buckets->B[actual_i]){
/**/    w->next = buckets->B[actual_i];
/**/    buckets->B[actual_i]->prev = w;
/**/  }
/**/  buckets->B[actual_i]= w;
/**/  buckets->in_bucket[w->index] = actual_i;
/**/}
/**/
/**/
/**/// Remove a node from bucket i and put it into bucket j
/**/void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets)
/**/{
/**/  llnode_t * node;
/**/  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
/**/          w->index, buckets->in_bucket[w->index],i, j);
/**/  if(i >= 0){
/**/    i = i % buckets->num_buckets;
/**/    if(buckets->in_bucket[w->index] == i){
/**/      node = buckets->B[i]; 
/**/      while(node != NULL){ // find w in B[i] if it is there
/**/	      if(node == w){
/**/	        remove_node_from_bucket(w, i, buckets);
/**/	        break;
/**/	      }
/**/	      node = node->next;
/**/      }
/**/    }
/**/  }
/**/  insert_node_in_bucket(w, j, buckets);   // insert w into the front of B[j] 
/**/}
/**/
/**/
/**/// relax an edge from a node to a node w (the current tenative distance to w is tent[windex]).
/**/// the candidate distance to w is cand_dist.
/**/void relax(int64_t windex, double cand_dist, double * tent, buckets_t * buckets)
/**/{
/**/  if ( cand_dist < tent[windex] ){
/**/    Dprintf("relax w=%"PRId64" cand_dist = %lf < %lf?\n", windex, cand_dist, tent[windex]);
/**/    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
/**/    int64_t j;
/**/    if(tent[windex] == INFINITY) 
/**/      j = -1;
/**/    else 
/**/      j = (int64_t)floor(tent[windex]/buckets->delta);
/**/    move_node_from_bucket_i_to_j(&(buckets->nodes[windex]), j, (int64_t)floor(cand_dist/buckets->delta), buckets);
/**/    tent[windex] = cand_dist;
/**/  }
/**/}
/**/
/**/// This is the delta stepping algorithm as it appears in
/**/// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
/**/// U. Meyer and P. Sanders.
/**/
/**/double sssp_delta_stepping(d_array_t *dist, sparsemat_t * mat, int64_t r0){
/**/  int64_t i, j;
/**/
/**/  double tm = wall_seconds();
/**/  
/**/  assert(r0 < mat->numrows);
/**/  
/**/  char * deleted = calloc(mat->numrows, sizeof(char));
/**/  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
/**/  
/**/  /* calculate delta and set tentative distances to infinity */  
/**/  double delta = 0.0;
/**/  int64_t max_degree = 0;
/**/  double * tent = dist->entry;
/**/  for(i = 0; i < mat->numrows; i++){
/**/    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
/**/      max_degree = (mat->offset[i+1] - mat->offset[i]);
/**/    tent[i] = INFINITY;
/**/  }
/**/  assert(max_degree > 0);
/**/  delta = 1.0/max_degree;
/**/  
/**/  double max_edge_weight = 0.0;
/**/  for(i = 0; i < mat->nnz; i++)
/**/    if(max_edge_weight < mat->value[i])
/**/      max_edge_weight = mat->value[i];
/**/  Dprintf("max edge weight = %lf\n", max_edge_weight);
/**/  
/**/  /* set up buckets as an array of linked lists */
/**/  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
/**/  buckets_t * buckets = calloc(1, sizeof(buckets_t));
/**/  buckets->num_buckets = num_buckets;
/**/  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
/**/  buckets->nodes = calloc(mat->numrows, sizeof(llnode_t));
/**/  buckets->in_bucket = calloc(mat->numrows, sizeof(int64_t));
/**/  buckets->delta = delta;
/**/  
/**/  Dprintf("num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", num_buckets, delta, max_degree);
/**/  for(i = 0; i < num_buckets; i++){
/**/    buckets->B[i] = NULL;
/**/  }
/**/  for(i = 0; i < mat->numrows; i++){
/**/    buckets->nodes[i].index = i;
/**/    buckets->nodes[i].next = NULL;
/**/    buckets->nodes[i].prev = NULL;
/**/    buckets->in_bucket[i] = -1;
/**/  }
/**/  
/**/  /* set the source distance to 0 */
/**/  insert_node_in_bucket(&(buckets->nodes[r0]), 0, buckets);
/**/  tent[r0] = 0.0;
/**/  
/**/  /* main loop */
/**/  int64_t current = 0;
/**/  while(1){
/**/
/**/    /* find the minimum indexed non-empty bucket */
/**/    for(i = 0; i < buckets->num_buckets; i++)
/**/      if(buckets->B[(current + i) % buckets->num_buckets] != NULL)
/**/        break;
/**/    
/**/    if(i == num_buckets)
/**/      break;
/**/    current = (current + i) % num_buckets;
/**/
/**/    Dprintf("Starting inner loop: working on bucket %"PRId64"\n", current);
/**/    int64_t start = 0;
/**/    int64_t end = 0;
/**/    llnode_t * v;
/**/    // inner loop
/**/    while(buckets->B[current] != NULL){
/**/      v = buckets->B[current];
/**/      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v->index, current);
/**/
/**/      /* relax light edges from v */
/**/      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
/**/        if(mat->value[j] <= delta){	  
/**/          relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
/**/        }
/**/      } 
/**/      
/**/      /* insert v into R if it is not already there */
/**/      if(deleted[v->index] == 0){
/**/        deleted[v->index] = 1;
/**/        R[end++] = v->index;
/**/        Dprintf("deleted %"PRId64"s\n", v->index);
/**/      }
/**/      
/**/      remove_node_from_bucket(v, current, buckets);
/**/    }
/**/
/**/    /* relax heavy requests edges for everything in R */
/**/    while(start < end){
/**/      v = &(buckets->nodes[R[start++]]);
/**/      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
/**/        if(mat->value[j] > delta){
/**/          relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
/**/        }
/**/      }      
/**/    }
/**/    current++;
/**/  }// end main loop
/**/  free(buckets->B);
/**/  free(buckets->nodes);
/**/  free(deleted);
/**/  free(R);
/**/  
/**/  return(wall_seconds() - tm);
/**/}

#else

#define DEBUG 1
typedef struct node_t{
  int64_t next;      // next in linked list
  int64_t prev;      // prev in linked list
  int64_t in_bucket; // which bucket it is in, else -1
  int64_t deleted;   // deleted means resolved?
  double  tent;      // the tentative weight for the vertex
}node_t;

typedef struct ds_t{
  int64_t num_buckets;
  int64_t * B;         // array of Buckets: B[i] is the index of the first node on the list, or -1 if empty
  node_t  * V;         // array of nodes to hold vertex information
  
  double  delta;       // delta
}ds_t;


// Remove node v from bucket i.
static void remove_node_from_bucket(ds_t *ds, int64_t v, int64_t i)
{
  int64_t i_m, first, last;
  i_m = i % ds->num_buckets;
  assert(ds->V[v].in_bucket == i_m);
  assert(ds->B[i_m] != -1);
  ds->V[v].in_bucket = -1;

  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v, i_m);

  if(ds->V[v].prev == -1){               // v is the first node of this list
    first = ds->V[v].next;                // new first, or -1 if list is now empty
    ds->B[i_m] = first;                   
    if(first >= 0) 
      ds->V[first].prev = -1;             // mark the new first as the first on the list
    return;
  }
  if(ds->V[v].next == -1){               // v is the last node on the list
    last = ds->V[v].prev;                // new last,  != -1 by previous case
    ds->V[last].next = -1;               // mark the new last as the end of the list
    return;
  } 
  //(ds->V[v].prev != -1) && (ds->V[v].next != -1)  // in the middle 
  ds->V[ds->V[v].prev].next = ds->V[v].next; 
  ds->V[ds->V[v].next].prev = ds->V[v].prev; 
  return;
}

// Prepend node v into a bucket i
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i)
{
  int64_t first;         // first node on list, if the list is non-empty
  int64_t i_m = (i % ds->num_buckets);

  Dprintf("Adding %"PRId64" to bucket %"PRId64"\n", v, i_m);
  
  if(ds->V[v].in_bucket == i_m){      // this node is already in this bucket
    return; 
  }
  assert(ds->V[v].in_bucket == -1);   // better not be in a different bucket

  ds->V[v].next = -1;
  ds->V[v].prev = -1;
  if(ds->B[i_m] != -1){               // something in the list
    first = ds->B[i_m];               // was the first on the list
    ds->V[v].next = first;            
    ds->V[first].prev = v;            
  }
  ds->B[i_m] = v;
  ds->V[v].in_bucket = i_m;
}


// Remove a node from bucket iold to bucket inew
// if the node is not in a bucket, just insert it into bucket inew
// if v is not in the old bucket, this is just an insert into the new bucket
void move_node_between_buckets(ds_t * ds, int64_t v, int64_t iold, int64_t inew)
{
  int64_t i_m_old, ll;
  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
          v, ds->V[v].in_bucket, iold, inew);
  if( ds->V[v].in_bucket >= 0 ){
    i_m_old = iold % ds->num_buckets;
    assert(i_m_old == ds->V[v].in_bucket);
    ll = ds->B[i_m_old];
    do {
	     if( v == ll ){
	       remove_node_from_bucket(ds, v, iold);
	       break;
	     }
    }while((ll = ds->V[ll].next) >= 0);
  }
  assert(ll != -1);
  insert_node_in_bucket(ds, v, inew); 
}

// relax an edge to the head vertex, given the new tentative distance
// (= the tentative distance to the tail plus the weight of the edge).
// the candidate distance to w is cand_dist.
void relax(ds_t *ds, int64_t j, double cand_dist)
{
  int64_t b;
  if ( cand_dist < ds->V[j].tent ){
    Dprintf("relax head %"PRId64" cand_dist = %lf < %lf?\n", j, cand_dist, ds->V[j].tent);
    /* if j is in B[floor((ds->V[j].tent)/delta)], remove it from that bucket */
    if(ds->V[j].tent == INFINITY) 
      b = -1;
    else 
      b = (int64_t)floor((ds->V[j].tent)/ds->delta);
    move_node_between_buckets(ds, j, b, (int64_t)floor(cand_dist/ds->delta));
    ds->V[j].tent = cand_dist;
  }
}

// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_stepping(d_array_t *dist, sparsemat_t * mat, int64_t r0)
{
  int64_t i, i_m, j, k;
  int64_t v;
  int64_t cur; // current bucket

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
  
  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t *));
  assert(ds != NULL);
  ds->V = (node_t *)calloc(mat->numrows, sizeof(node_t));
  assert(ds->V != NULL);
  for(i = 0; i < mat->numrows; i++){
    ds->V[i].next      = -1;
    ds->V[i].prev      = -1;
    ds->V[i].in_bucket = -1;
    ds->V[i].deleted   =  0;
    ds->V[i].tent      = INFINITY;
  }
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t));
  assert(ds->B != NULL);
  for(i_m = 0; i_m < num_buckets; i_m++){
    ds->B[i_m] = -1;
  }
  ds->delta = delta;

  /* set the source distance to r0 */
  ds->V[r0].tent = 0.0;
  cur = 0;
  insert_node_in_bucket(ds, r0, cur);
  
  /* main loop */
  for(cur=0;  ; cur++){
    // find the minimum indexed non-empty bucket 
    for(i_m = 0; i_m < ds->num_buckets; i_m++){
      if(ds->B[(cur + i_m) % ds->num_buckets] != -1)
        break;
    }
    if(i_m == num_buckets)  //No non-empty buckets, we are done
      break;
    i = cur + i_m;
    Dprintf("Starting inner loop: working on i = %"PRId64" in bucket %"PRId64"\n", i, i_m);

    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    for( j=ds->B[i_m]; ds->V[j].in_bucket == i_m; j=ds->V[j].next){
      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", j, i_m);

      remove_node_from_bucket(ds, j, i_m);
      
      /* relax light edges from v */
      for(k = mat->offset[j]; k < mat->offset[j + 1]; k++){
        if(mat->value[k] <= delta){	  
          relax(ds, mat->nonzero[k], ds->V[j].tent + mat->value[k]);
        }
      } 
      
      /* insert v into R if it is not already there */
      if(ds->V[j].deleted == 0){
        ds->V[j].deleted = 1;
        R[end++] = j;
        Dprintf("deleted %"PRId64"s\n", j);
      }
    }

    /* relax heavy requests edges for everything in R */
    for(start=0; start<end; start++){
      j = R[start];
      for(k = mat->offset[j]; k < mat->offset[j + 1]; k++){
        if(mat->value[k] > delta){
          relax(ds, mat->nonzero[k], ds->V[j].tent + mat->value[k]);
        }
      }      
    }
  }// end main loop

  Dprintf("copy tentative to distance array\n");
  // Copy the answers to dist array
  for(v = 0; v < mat->numrows; v++){
    dist->entry[v] = ds->V[v].tent;
  }

  printf("Delta   : ");
  for(i=0; i<mat->numrows; i++)
    printf(" %lf", dist->entry[i]);
  printf("\n");

  Dprintf("free stuff\n");

  free(ds->B);
//  free(ds->V);
  free(ds);
  free(R);
  
  return(wall_seconds() - tm);
}
#endif
