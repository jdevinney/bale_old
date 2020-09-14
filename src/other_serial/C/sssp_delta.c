#include "spmat_utils.h"

#define POINTER 0
#define STRUCT 0
#define ARRAY 1

#if POINTER
/*P*/typedef struct llnode_t{
/*P*/  int64_t index;
/*P*/  struct llnode_t * next;
/*P*/  struct llnode_t * prev;
/*P*/}llnode_t;
/*P*/
/*P*/typedef struct buckets_t{
/*P*/  llnode_t ** B;
/*P*/  llnode_t * nodes;
/*P*/  int64_t * in_bucket;
/*P*/  int64_t num_buckets;
/*P*/  double delta;
/*P*/}buckets_t;
/*P*/
/*P*/
/*P*/// Remove a specific node from bucket i.
/*P*/void remove_node_from_bucket(llnode_t * v, int64_t i, buckets_t * buckets)
/*P*/{
/*P*/  i = i % buckets->num_buckets;
/*P*/  assert(buckets->in_bucket[v->index] == i);
/*P*/  buckets->in_bucket[v->index] = -1;
/*P*/  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v->index, i);
/*P*/  if(v->prev){
/*P*/    v->prev->next = v->next;
/*P*/  }else
/*P*/    buckets->B[i] = v->next;
/*P*/  if(v->next != NULL)
/*P*/    v->next->prev = v->prev;
/*P*/}
/*P*/
/*P*/// Prepend a node into a bucket
/*P*/void insert_node_in_bucket(llnode_t * w, int64_t i, buckets_t * buckets)
/*P*/{
/*P*/  int64_t actual_i = (i % buckets->num_buckets);
/*P*/
/*P*/  Dprintf("Adding %"PRId64" to bucket %"PRId64"\n", w->index, actual_i);
/*P*/  
/*P*/  if(buckets->in_bucket[w->index] == actual_i)  // this node is already in this bucket
/*P*/    return; 
/*P*/
/*P*/  assert(buckets->in_bucket[w->index] == -1);
/*P*/
/*P*/  w->next = NULL;
/*P*/  w->prev = NULL;
/*P*/  if(buckets->B[actual_i]){
/*P*/    w->next = buckets->B[actual_i];
/*P*/    buckets->B[actual_i]->prev = w;
/*P*/  }
/*P*/  buckets->B[actual_i]= w;
/*P*/  buckets->in_bucket[w->index] = actual_i;
/*P*/}
/*P*/
/*P*/
/*P*/// Remove a node from bucket i and put it into bucket j
/*P*/void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets)
/*P*/{
/*P*/  llnode_t * node;
/*P*/  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
/*P*/          w->index, buckets->in_bucket[w->index],i, j);
/*P*/  if(i >= 0){
/*P*/    i = i % buckets->num_buckets;
/*P*/    if(buckets->in_bucket[w->index] == i){
/*P*/      node = buckets->B[i]; 
/*P*/      while(node != NULL){ // find w in B[i] if it is there
/*P*/	      if(node == w){
/*P*/	        remove_node_from_bucket(w, i, buckets);
/*P*/	        break;
/*P*/	      }
/*P*/	      node = node->next;
/*P*/      }
/*P*/    }
/*P*/  }
/*P*/  insert_node_in_bucket(w, j, buckets);   // insert w into the front of B[j] 
/*P*/}
/*P*/
/*P*/
/*P*/// relax an edge from a node to a node w (the current tenative distance to w is tent[windex]).
/*P*/// the candidate distance to w is cand_dist.
/*P*/void relax(int64_t windex, double cand_dist, double * tent, buckets_t * buckets)
/*P*/{
/*P*/  if ( cand_dist < tent[windex] ){
/*P*/    Dprintf("relax w=%"PRId64" cand_dist = %lf < %lf?\n", windex, cand_dist, tent[windex]);
/*P*/    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
/*P*/    int64_t j;
/*P*/    if(tent[windex] == INFINITY) 
/*P*/      j = -1;
/*P*/    else 
/*P*/      j = (int64_t)floor(tent[windex]/buckets->delta);
/*P*/    move_node_from_bucket_i_to_j(&(buckets->nodes[windex]), j, (int64_t)floor(cand_dist/buckets->delta), buckets);
/*P*/    tent[windex] = cand_dist;
/*P*/  }
/*P*/}
/*P*/
/*P*/// This is the delta stepping algorithm as it appears in
/*P*/// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
/*P*/// U. Meyer and P. Sanders.
/*P*/
/*P*/double sssp_delta_stepping(d_array_t *dist, sparsemat_t * mat, int64_t r0){
/*P*/  int64_t i, j;
/*P*/
/*P*/  double tm = wall_seconds();
/*P*/  
/*P*/  assert(r0 < mat->numrows);
/*P*/  
/*P*/  char * deleted = calloc(mat->numrows, sizeof(char));
/*P*/  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
/*P*/  
/*P*/  /* calculate delta and set tentative distances to infinity */  
/*P*/  double delta = 0.0;
/*P*/  int64_t max_degree = 0;
/*P*/  double * tent = dist->entry;
/*P*/  for(i = 0; i < mat->numrows; i++){
/*P*/    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
/*P*/      max_degree = (mat->offset[i+1] - mat->offset[i]);
/*P*/    tent[i] = INFINITY;
/*P*/  }
/*P*/  assert(max_degree > 0);
/*P*/  delta = 1.0/max_degree;
/*P*/  
/*P*/  double max_edge_weight = 0.0;
/*P*/  for(i = 0; i < mat->nnz; i++)
/*P*/    if(max_edge_weight < mat->value[i])
/*P*/      max_edge_weight = mat->value[i];
/*P*/  Dprintf("max edge weight = %lf\n", max_edge_weight);
/*P*/  
/*P*/  /* set up buckets as an array of linked lists */
/*P*/  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
/*P*/  buckets_t * buckets = calloc(1, sizeof(buckets_t));
/*P*/  buckets->num_buckets = num_buckets;
/*P*/  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
/*P*/  buckets->nodes = calloc(mat->numrows, sizeof(llnode_t));
/*P*/  buckets->in_bucket = calloc(mat->numrows, sizeof(int64_t));
/*P*/  buckets->delta = delta;
/*P*/  
/*P*/  Dprintf("num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", num_buckets, delta, max_degree);
/*P*/  for(i = 0; i < num_buckets; i++){
/*P*/    buckets->B[i] = NULL;
/*P*/  }
/*P*/  for(i = 0; i < mat->numrows; i++){
/*P*/    buckets->nodes[i].index = i;
/*P*/    buckets->nodes[i].next = NULL;
/*P*/    buckets->nodes[i].prev = NULL;
/*P*/    buckets->in_bucket[i] = -1;
/*P*/  }
/*P*/  
/*P*/  /* set the source distance to 0 */
/*P*/  insert_node_in_bucket(&(buckets->nodes[r0]), 0, buckets);
/*P*/  tent[r0] = 0.0;
/*P*/  
/*P*/  /* main loop */
/*P*/  int64_t current = 0;
/*P*/  while(1){
/*P*/
/*P*/    /* find the minimum indexed non-empty bucket */
/*P*/    for(i = 0; i < buckets->num_buckets; i++)
/*P*/      if(buckets->B[(current + i) % buckets->num_buckets] != NULL)
/*P*/        break;
/*P*/    
/*P*/    if(i == num_buckets)
/*P*/      break;
/*P*/    current = (current + i) % num_buckets;
/*P*/
/*P*/    Dprintf("Starting inner loop: working on bucket %"PRId64"\n", current);
/*P*/    int64_t start = 0;
/*P*/    int64_t end = 0;
/*P*/    llnode_t * v;
/*P*/    // inner loop
/*P*/    while(buckets->B[current] != NULL){
/*P*/      v = buckets->B[current];
/*P*/      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v->index, current);
/*P*/
/*P*/      /* relax light edges from v */
/*P*/      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
/*P*/        if(mat->value[j] <= delta){	  
/*P*/          relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
/*P*/        }
/*P*/      } 
/*P*/      
/*P*/      /* insert v into R if it is not already there */
/*P*/      if(deleted[v->index] == 0){
/*P*/        deleted[v->index] = 1;
/*P*/        R[end++] = v->index;
/*P*/        Dprintf("deleted %"PRId64"s\n", v->index);
/*P*/      }
/*P*/      
/*P*/      remove_node_from_bucket(v, current, buckets);
/*P*/    }
/*P*/
/*P*/    /* relax heavy requests edges for everything in R */
/*P*/    while(start < end){
/*P*/      v = &(buckets->nodes[R[start++]]);
/*P*/      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
/*P*/        if(mat->value[j] > delta){
/*P*/          relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
/*P*/        }
/*P*/      }      
/*P*/    }
/*P*/    current++;
/*P*/  }// end main loop
/*P*/  free(buckets->B);
/*P*/  free(buckets->nodes);
/*P*/  free(deleted);
/*P*/  free(R);
/*P*/  
/*P*/  return(wall_seconds() - tm);
/*P*/}

#endif

#if STRUCT 

/*S*/#define DEBUG 1
/*S*/typedef struct node_t{
/*S*/  int64_t next;      // next in linked list
/*S*/  int64_t prev;      // prev in linked list
/*S*/  int64_t in_bucket; // which bucket it is in, else -1
/*S*/  int64_t deleted;   // deleted means resolved?
/*S*/  double  tent;      // the tentative weight for the vertex
/*S*/}node_t;
/*S*/
/*S*/typedef struct ds_t{
/*S*/  int64_t num_buckets;
/*S*/  int64_t * B;         // array of Buckets: B[i] is the index of the first node on the list, or -1 if empty
/*S*/  node_t  * V;         // array of nodes to hold vertex information
/*S*/  
/*S*/  double  delta;       // delta
/*S*/}ds_t;
/*S*/
/*S*/
/*S*/// Remove node v from bucket i.
/*S*/static void remove_node_from_bucket(ds_t *ds, int64_t v, int64_t i)
/*S*/{
/*S*/  int64_t i_m, first, last;
/*S*/  i_m = i % ds->num_buckets;
/*S*/  assert(ds->V[v].in_bucket == i_m);
/*S*/  assert(ds->B[i_m] != -1);
/*S*/  ds->V[v].in_bucket = -1;
/*S*/
/*S*/  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v, i_m);
/*S*/
/*S*/  if(ds->V[v].prev == -1){               // v is the first node of this list
/*S*/    first = ds->V[v].next;                // new first, or -1 if list is now empty
/*S*/    ds->B[i_m] = first;                   
/*S*/    if(first >= 0) 
/*S*/      ds->V[first].prev = -1;             // mark the new first as the first on the list
/*S*/    return;
/*S*/  }
/*S*/  if(ds->V[v].next == -1){               // v is the last node on the list
/*S*/    last = ds->V[v].prev;                // new last,  != -1 by previous case
/*S*/    ds->V[last].next = -1;               // mark the new last as the end of the list
/*S*/    return;
/*S*/  } 
/*S*/  //(ds->V[v].prev != -1) && (ds->V[v].next != -1)  // in the middle 
/*S*/  ds->V[ds->V[v].prev].next = ds->V[v].next; 
/*S*/  ds->V[ds->V[v].next].prev = ds->V[v].prev; 
/*S*/  return;
/*S*/}
/*S*/
/*S*/// Prepend node v into a bucket i
/*S*/void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i)
/*S*/{
/*S*/  int64_t first;         // first node on list, if the list is non-empty
/*S*/  int64_t i_m = (i % ds->num_buckets);
/*S*/
/*S*/  Dprintf("Adding %"PRId64" to bucket %"PRId64"\n", v, i_m);
/*S*/  
/*S*/  if(ds->V[v].in_bucket == i_m){      // this node is already in this bucket
/*S*/    return; 
/*S*/  }
/*S*/  assert(ds->V[v].in_bucket == -1);   // better not be in a different bucket
/*S*/
/*S*/  ds->V[v].next = -1;
/*S*/  ds->V[v].prev = -1;
/*S*/  if(ds->B[i_m] != -1){               // something in the list
/*S*/    first = ds->B[i_m];               // was the first on the list
/*S*/    ds->V[v].next = first;            
/*S*/    ds->V[first].prev = v;            
/*S*/  }
/*S*/  ds->B[i_m] = v;
/*S*/  ds->V[v].in_bucket = i_m;
/*S*/}
/*S*/
/*S*/
/*S*/// Remove a node from bucket iold to bucket inew
/*S*/// if the node is not in a bucket, just insert it into bucket inew
/*S*/// if v is not in the old bucket, this is just an insert into the new bucket
/*S*/void move_node_between_buckets(ds_t * ds, int64_t v, int64_t iold, int64_t inew)
/*S*/{
/*S*/  int64_t i_m_old, ll;
/*S*/  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
/*S*/          v, ds->V[v].in_bucket, iold, inew);
/*S*/  if( ds->V[v].in_bucket >= 0 ){
/*S*/    i_m_old = iold % ds->num_buckets;
/*S*/    assert(i_m_old == ds->V[v].in_bucket);
/*S*/    ll = ds->B[i_m_old];
/*S*/    do {
/*S*/	     if( v == ll ){
/*S*/	       remove_node_from_bucket(ds, v, iold);
/*S*/	       break;
/*S*/	     }
/*S*/    }while((ll = ds->V[ll].next) >= 0);
/*S*/  }
/*S*/  assert(ll != -1);
/*S*/  insert_node_in_bucket(ds, v, inew); 
/*S*/}
/*S*/
/*S*/// relax an edge to the head vertex, given the new tentative distance
/*S*/// (= the tentative distance to the tail plus the weight of the edge).
/*S*/// the candidate distance to w is cand_dist.
/*S*/void relax(ds_t *ds, int64_t j, double cand_dist)
/*S*/{
/*S*/  int64_t b;
/*S*/  if ( cand_dist < ds->V[j].tent ){
/*S*/    Dprintf("relax head %"PRId64" cand_dist = %lf < %lf?\n", j, cand_dist, ds->V[j].tent);
/*S*/    /* if j is in B[floor((ds->V[j].tent)/delta)], remove it from that bucket */
/*S*/    if(ds->V[j].tent == INFINITY) 
/*S*/      b = -1;
/*S*/    else 
/*S*/      b = (int64_t)floor((ds->V[j].tent)/ds->delta);
/*S*/    move_node_between_buckets(ds, j, b, (int64_t)floor(cand_dist/ds->delta));
/*S*/    ds->V[j].tent = cand_dist;
/*S*/  }
/*S*/}
/*S*/
/*S*/// This is the delta stepping algorithm as it appears in
/*S*/// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
/*S*/// U. Meyer and P. Sanders.
/*S*/
/*S*/double sssp_delta_stepping(d_array_t *dist, sparsemat_t * mat, int64_t r0)
/*S*/{
/*S*/  int64_t i, i_m, j, k;
/*S*/  int64_t v;
/*S*/  int64_t cur; // current bucket
/*S*/
/*S*/  double tm = wall_seconds();
/*S*/
/*S*/  assert((r0 >= 0) && (r0<mat->numrows));
/*S*/  
/*S*/  int64_t * R = calloc(mat->numrows, sizeof(int64_t));
/*S*/  
/*S*/  /* calculate delta and set tentative distances to infinity */  
/*S*/  double delta = 0.0;
/*S*/  int64_t max_degree = 0;
/*S*/  for(i = 0; i < mat->numrows; i++){
/*S*/    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
/*S*/      max_degree = (mat->offset[i+1] - mat->offset[i]);
/*S*/  }
/*S*/  assert(max_degree > 0);
/*S*/  delta = 1.0/max_degree;
/*S*/  
/*S*/  double max_edge_weight = 0.0;
/*S*/  for(i = 0; i < mat->nnz; i++)
/*S*/    if(max_edge_weight < mat->value[i])
/*S*/      max_edge_weight = mat->value[i];
/*S*/  Dprintf("max edge weight = %lf\n", max_edge_weight);
/*S*/  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
/*S*/  
/*S*/  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t));
/*S*/  assert(ds != NULL);
/*S*/  ds->V = (node_t *)calloc(mat->numrows, sizeof(node_t));
/*S*/  assert(ds->V != NULL);
/*S*/  for(i = 0; i < mat->numrows; i++){
/*S*/    ds->V[i].next      = -1;
/*S*/    ds->V[i].prev      = -1;
/*S*/    ds->V[i].in_bucket = -1;
/*S*/    ds->V[i].deleted   =  0;
/*S*/    ds->V[i].tent      = INFINITY;
/*S*/  }
/*S*/  ds->num_buckets = num_buckets;
/*S*/  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t));
/*S*/  assert(ds->B != NULL);
/*S*/  for(i_m = 0; i_m < num_buckets; i_m++){
/*S*/    ds->B[i_m] = -1;
/*S*/  }
/*S*/  ds->delta = delta;
/*S*/
/*S*/  /* set the source distance to r0 */
/*S*/  ds->V[r0].tent = 0.0;
/*S*/  cur = 0;
/*S*/  insert_node_in_bucket(ds, r0, cur);
/*S*/  
/*S*/  /* main loop */
/*S*/  for(cur=0;  ; cur++){
/*S*/    // find the minimum indexed non-empty bucket 
/*S*/    for(i_m = 0; i_m < ds->num_buckets; i_m++){
/*S*/      if(ds->B[(cur + i_m) % ds->num_buckets] != -1)
/*S*/        break;
/*S*/    }
/*S*/    if(i_m == num_buckets)  //No non-empty buckets, we are done
/*S*/      break;
/*S*/    i = cur + i_m;
/*S*/    Dprintf("Starting inner loop: working on i = %"PRId64" in bucket %"PRId64"\n", i, i_m);
/*S*/
/*S*/    // inner loop
/*S*/    int64_t start = 0;
/*S*/    int64_t end = 0;
/*S*/    for( j=ds->B[i_m]; ds->V[j].in_bucket == i_m; j=ds->V[j].next){
/*S*/      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", j, i_m);
/*S*/
/*S*/      remove_node_from_bucket(ds, j, i_m);
/*S*/      
/*S*/      /* relax light edges from v */
/*S*/      for(k = mat->offset[j]; k < mat->offset[j + 1]; k++){
/*S*/        if(mat->value[k] <= delta){	  
/*S*/          relax(ds, mat->nonzero[k], ds->V[j].tent + mat->value[k]);
/*S*/        }
/*S*/      } 
/*S*/      
/*S*/      /* insert v into R if it is not already there */
/*S*/      if(ds->V[j].deleted == 0){
/*S*/        ds->V[j].deleted = 1;
/*S*/        R[end++] = j;
/*S*/        Dprintf("deleted %"PRId64"s\n", j);
/*S*/      }
/*S*/    }
/*S*/
/*S*/    /* relax heavy requests edges for everything in R */
/*S*/    for(start=0; start<end; start++){
/*S*/      j = R[start];
/*S*/      for(k = mat->offset[j]; k < mat->offset[j + 1]; k++){
/*S*/        if(mat->value[k] > delta){
/*S*/          relax(ds, mat->nonzero[k], ds->V[j].tent + mat->value[k]);
/*S*/        }
/*S*/      }      
/*S*/    }
/*S*/  }// end main loop
/*S*/
/*S*/  Dprintf("copy tentative to distance array\n");
/*S*/  // Copy the answers to dist array
/*S*/  for(v = 0; v < mat->numrows; v++){
/*S*/    dist->entry[v] = ds->V[v].tent;
/*S*/  }
/*S*/
/*S*/  assert(ds->V != NULL);
/*S*/  //free(ds->V);
/*S*/  free(ds->B);
/*S*/  free(ds);
/*S*/  free(R);
/*S*/  
/*S*/  return(wall_seconds() - tm);
/*S*/}
#endif

#if ARRAY 

#define DEBUG 1
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


// Prepend node v into a bucket i
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m)
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
    int64_t pv = ds->prev[w];
    Dprintf("non-empty: w=%ld, prev=%ld\n", w, pv);
    ds->next[pv] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    ds->prev[v] = pv;            
    Dprintf("   v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );
  }
  ds->B[i_m] = v;                       // set v to be the "first" on the list
  ds->in_bucket[v] = i_m;
}


// Remove node v from its bucket (need not be in any bucket)
static void remove_node_from_bucket(ds_t *ds, int64_t v)
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

#if 0
// Remove a node from bucket iold to bucket inew
// if the node is not in a bucket, just insert it into bucket inew
// if v is not in the old bucket, this is just an insert into the new bucket
void move_node_between_buckets(ds_t * ds, int64_t v, int64_t iold, int64_t inew)
{
  int64_t i_m_old, w;
  Dprintf("Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", 
          v, ds->in_bucket[v], iold, inew);
  if( ds->in_bucket[v] >= 0 ){
    i_m_old = iold % ds->num_buckets;
    assert(i_m_old == ds->in_bucket[v]);
    // linear search the bucket
    w = ds->B[i_m_old]; 

    ll >= 0; ll = ds->next[ll]){

    for(ll = ds->B[i_m_old]; ll >= 0; ll = ds->next[ll]){
	    if( v == ll ){
	      remove_node_from_bucket(ds, v, iold);
	      break;
	    }
    }
  }
  assert(ll != -1);    // should have found v in the bucket
  insert_node_in_bucket(ds, v, inew); 
}

#endif

int64_t home_bucket(ds_t *ds, double dist);
// relax an edge to the head vertex, given the new tentative distance
// (= the tentative distance to the tail plus the weight of the edge).
// the candidate distance to w is cand_dist.
void relax(ds_t *ds, int64_t w, double cand_dist)
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
        remove_node_from_bucket(ds, w);
      insert_node_in_bucket(ds, w, inew);
    }
    ds->tent[w] = cand_dist;
  }
}

void dump_bucket(ds_t *ds, int64_t i_m)
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
// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_stepping(d_array_t *dist, sparsemat_t * mat, int64_t r0)
{
  int64_t i, i_m, k;
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
  insert_node_in_bucket(ds, r0, 0);
  Dprintf("putting source r0=%ld in bucket %ld\n", r0, 0L);
  
  /* main loop */
  for(cur=0;  ; cur++){
    Dprintf("Outer loop cur = %"PRId64"\n", cur);
    // find the minimum indexed non-empty bucket 
    for(i = 0; i < ds->num_buckets; i++){
      if(ds->B[(cur + i) % ds->num_buckets] != -1)
        break;
    }
    if(i == num_buckets)  //All buckets are empty, we are done
      break;
    i_m = (cur + i) % ds->num_buckets;
    Dprintf("Starting inner loop: working on bucket %"PRId64"\n", i_m);
    if(DEBUG) dump_bucket(ds, i_m);

    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    for( v=ds->B[i_m]; ds->in_bucket[v] == i_m; v=ds->next[v]){
      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v, i_m);

      remove_node_from_bucket(ds, v);
      
      /* relax light edges from v */
      for(k = mat->offset[v]; k < mat->offset[v + 1]; k++){
        if(mat->value[k] <= delta){	  
          relax(ds, mat->nonzero[k], ds->tent[v] + mat->value[k]);
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
          relax(ds, mat->nonzero[k], ds->tent[v] + mat->value[k]);
        }
      }      
    }
  }// end main loop

  Dprintf("copy tentative to distance array\n");
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
#endif
