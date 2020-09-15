/*! \file sssp_delta_exstack.upc
 * \brief Implementation of delta stepping using exstack
 */

#include "sssp.h"

#define DPRT 1

#if 0
/**/typedef struct pkg_delta_e_t{
/**/  int64_t i;  // row[i] is "tail" of the edge, in case we want to set backpointers
/**/  int64_t lj; // the local "head" of the edge
/**/  double tw;  // new tentative weight
/**/}pkg_delta_e_t;
/**/
/**/
/**/typedef struct llnode_t{
/**/  int64_t index;
/**/  struct llnode_t * next;
/**/  struct llnode_t * prev;
/**/}llnode_t;
/**/
/**/
/**/typedef struct buckets_t{
/**/  llnode_t ** B;
/**/  llnode_t * nodes;
/**/  int64_t * level;
/**/  char * empty;
/**/  int64_t * in_bucket;
/**/  int64_t num_buckets;
/**/  double delta;
/**/}buckets_t;
/**/
/**/typedef struct all_buckets_t {
/**/  SHARED buckets_t * Sbuckets;
/**/  buckets_t * lbuckets;
/**/}all_buckets_t;
/**/
/**/// Remove a specific node from bucket i.
/**/void remove_node_from_bucket(llnode_t * v, int64_t i, buckets_t * buckets){
/**/  llnode_t * next = v->next;
/**/  i = i % buckets->num_buckets;
/**/  assert(buckets->in_bucket[v->index] == i);
/**/  assert(buckets->empty[i] == 0);
/**/  buckets->in_bucket[v->index] = -1;
/**/  buckets->level[i] = -1;
/**/  if(DPRT){printf("%02d: Removing %"PRId64" from bucket %"PRId64"\n",MYTHREAD, v->index, i);}
/**/  if(v->prev)
/**/    v->prev->next = next;
/**/  else
/**/    buckets->B[i] = v->next;
/**/  if(next != NULL)
/**/    next->prev = v->prev;
/**/  if(buckets->B[i] == NULL)
/**/    buckets->empty[i] = 1;
/**/}
/**/
/**/// Prepend a node into a bucket
/**/void insert_node_in_bucket(llnode_t * w, int64_t i, buckets_t * buckets){
/**/  int64_t actual_i = (i % buckets->num_buckets);
/**/
/**/  if(DPRT){printf("%02d:Adding %"PRId64" to bucket %"PRId64"\n", MYTHREAD, w->index, actual_i);}
/**/  
/**/  if(buckets->in_bucket[w->index] == actual_i){
/**/    //this node is already in this bucket
/**/    if(buckets->level[actual_i] != (actual_i % buckets->num_buckets))
/**/      if(DPRT){printf("%02d:%"PRId64" %"PRId64"\n", MYTHREAD, buckets->level[actual_i], (actual_i % buckets->num_buckets));}
/**/    assert(buckets->level[actual_i] == (actual_i % buckets->num_buckets));
/**/    return; 
/**/  }
/**/
/**/  if(buckets->in_bucket[w->index] != -1){
/**/    printf("AAAhhh: %"PRId64"\n", buckets->in_bucket[w->index]);
/**/  }
/**/  assert(buckets->in_bucket[w->index] == -1);
/**/  
/**/
/**/  if(buckets->empty[actual_i] == 1){
/**/    // this bucket is empty right now
/**/    assert(buckets->B[actual_i] == NULL);
/**/    buckets->empty[actual_i] = 0;
/**/  }else{
/**/    if(DPRT){printf("%02d:%"PRId64" %"PRId64"\n", MYTHREAD, buckets->level[actual_i], (actual_i % buckets->num_buckets));}
/**/    assert(buckets->level[actual_i] == -1 || buckets->level[actual_i] == (actual_i % buckets->num_buckets));
/**/  }
/**/  
/**/  buckets->level[actual_i] = (actual_i % buckets->num_buckets);
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
/**/void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets){
/**/
/**/  if(DPRT){printf("%02d:Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", MYTHREAD, w->index, buckets->in_bucket[w->index],i, j);}
/**/  if(i >= 0){
/**/    i = i % buckets->num_buckets;
/**/    if(buckets->in_bucket[w->index] == i){
/**/      
/**/      llnode_t * node = buckets->B[i];
/**/      
/**/      /* remove w from B[i] if it is there */
/**/      while(node != NULL){
/**/        if(node == w){
/**/          remove_node_from_bucket(w, i, buckets);
/**/          break;
/**/        }
/**/        node = node->next;
/**/      }
/**/    }
/**/  }
/**/  
/**/  /* insert w into the front of B[j] */
/**/  insert_node_in_bucket(w, j, buckets);
/**/}
/**/
/**/// CLAIM:  while Thread foo is foolin with the Buckets,  it can't be reading the buckets for pushes.
/**///       So there is no race on the list.
/**///      But when the Thread starts pushing again, the list could be different.
/**///      Different, but not Busted.
/**/
/**/
/**/// relax an edge from a node to a node w (the current tenative distance to w is tent[windex]).
/**/// the candidate distance to w is cand_dist.
/**/static void local_relax(int64_t windex, double cand_dist, d_array_t * tent, buckets_t * buckets){
/**/  
/**/  if ( cand_dist < tent->lentry[windex] ){
/**/    if(DPRT){printf("%02d:relax w=%"PRId64" cand_dist = %lf < %lf?\n", MYTHREAD, windex, cand_dist, tent->lentry[windex]);}
/**/    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
/**/    int64_t j;
/**/    if(tent->lentry[windex] == INFINITY) 
/**/      j = -1;
/**/    else 
/**/      j = (int64_t)floor(tent->lentry[windex]/buckets->delta);
/**/    move_node_from_bucket_i_to_j(&(buckets->nodes[windex]), j, (int64_t)floor(cand_dist/buckets->delta), buckets);
/**/    tent->lentry[windex] = cand_dist;
/**/  }
/**/}
/**/
/**//*!
/**/ * \brief pop routine to implement relaxing the edges
/**/ * \param tent pointer to the tentative distances array
/**/ * \param *ex the extack buffers
/**/ * \param done the signal to exstack_proceed that this thread is done
/**/ * \return the return value from exstack_proceed
/**/ */
/**/static int64_t delta_exstack_relax_process(d_array_t *tent, exstack_t *ex, all_buckets_t *Buckets,  int64_t done) 
/**/{
/**/  buckets_t * buckets = Buckets->lbuckets;
/**/  int64_t fromth;
/**/  pkg_delta_e_t pkg;
/**/
/**/  exstack_exchange(ex);
/**/  
/**/  while(exstack_pop(ex, &pkg, &fromth)){
/**/    local_relax(pkg.lj, pkg.tw, tent, buckets);
/**/  }
/**/  return( exstack_proceed(ex, done) );
/**/}
/**/
/**/
/**/
/**/// This is the delta stepping algorithm as it appears in
/**/// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
/**/// U. Meyer and P. Sanders.
/**/
/**/double sssp_delta_exstack(d_array_t *tent, sparsemat_t * mat, int64_t r0)
/**/{
/**/  int64_t i, li, j;
/**/  int64_t J, pe;
/**/  pkg_delta_e_t pkg;
/**/  int64_t lr0 = r0 / THREADS;   // the local row (vertex) name for r0
/**/  int64_t rbi;
/**/
/**/  //TODO: Fix the buffer size 
/**/  exstack_t * ex = exstack_init(2, sizeof(pkg_delta_e_t));
/**/  if( ex == NULL) return(-1.0);
/**/
/**/  double tm = wall_seconds();
/**/  
/**/  assert(r0 < mat->numrows);
/**/  
/**/  char * deleted = calloc(mat->lnumrows, sizeof(char));
/**/  int64_t * R = calloc(mat->lnumrows, sizeof(int64_t));
/**/  
/**/  /* calculate delta and set tentative distances to infinity */  
/**/  double delta = 0.0;
/**/  int64_t max_degree = 0;
/**/  for(li = 0; li < mat->lnumrows; li++){
/**/    if(max_degree < (mat->loffset[i+1] - mat->loffset[i]))
/**/      max_degree = (mat->loffset[i+1] - mat->loffset[i]);
/**/  }
/**/  max_degree = lgp_reduce_max_l(max_degree);
/**/  assert(max_degree > 0);
/**/  delta = 1.0/max_degree;
/**/  if(DPRT){printf("%02d:delta = %lf\n", MYTHREAD, delta);}
/**/  
/**/  double max_edge_weight = 0.0;
/**/  for(li = 0; li < mat->lnnz; li++)
/**/    if(max_edge_weight < mat->lvalue[li])
/**/      max_edge_weight = mat->lvalue[li];
/**/  max_edge_weight = lgp_reduce_max_d(max_edge_weight);
/**/
/**/  if(DPRT){printf("%02d:max edge weight = %lf\n", MYTHREAD, max_edge_weight);}
/**/  
/**/  if(DPRT){printf("%02d:Init Buckets\n", MYTHREAD);}
/**/
/**/
/**/  /* set up buckets as an array of linked lists */
/**/
/**/  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
/**/
/**/#if 0
/**/  /* set up buckets as an array of linked lists */
/**/  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
/**/  buckets_t * buckets = calloc(1, sizeof(buckets_t));
/**/#endif
/**/  all_buckets_t * Buckets= calloc(1, sizeof(all_buckets_t));
/**/  Buckets->Sbuckets = lgp_all_alloc(THREADS, sizeof(buckets_t));
/**/  Buckets->lbuckets = lgp_local_part(buckets_t, Buckets->Sbuckets);
/**/
/**/  buckets_t * buckets = Buckets->lbuckets;
/**/
/**/  buckets->num_buckets = num_buckets;
/**/  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
/**/  buckets->nodes = calloc(mat->lnumrows, sizeof(llnode_t));
/**/  buckets->in_bucket = calloc(mat->lnumrows, sizeof(int64_t));
/**/  buckets->level = calloc(num_buckets, sizeof(int64_t));
/**/  buckets->empty = calloc(num_buckets, sizeof(char));
/**/  buckets->delta = delta;
/**/  
/**/  if(DPRT){printf("%02d:num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", MYTHREAD, num_buckets, delta, max_degree);}
/**/  for(i = 0; i < num_buckets; i++){
/**/    buckets->B[i] = NULL;
/**/    buckets->empty[i] = 1;
/**/    buckets->level[i] = -1;
/**/  }
/**/  for(i = 0; i < mat->lnumrows; i++){
/**/    buckets->nodes[i].index = i;
/**/    buckets->nodes[i].next = NULL;
/**/    buckets->nodes[i].prev = NULL;
/**/    buckets->in_bucket[i] = -1;
/**/  }
/**/
/**/  lgp_barrier();
/**/  for(li = 0; li < mat->lnumrows; li++)
/**/    tent->lentry[li] = INFINITY;
/**/
/**/  
/**/  /* set the source distance to 0 */
/**/  if( (r0 % THREADS) == MYTHREAD) {
/**/    insert_node_in_bucket(&(buckets->nodes[lr0]), 0, buckets);
/**/    tent->lentry[lr0] = 0.0;
/**/  }
/**/
/**/  lgp_barrier();
/**/  dump_tent(">>Delta Exstack:", tent);
/**/  lgp_barrier();
/**/
/**/  /* main loop */
/**/  //int64_t min_bucket = 0;
/**/  int64_t all_cur=0, current_bucket = 0;
/**/  int64_t num_deleted = 0;
/**/  rbi = 0;
/**/
/**/  while(1){
/**/
/**/    /* find the minimum indexed non-empty bucket */
/**/    for(i = 0; i < num_buckets; i++)
/**/      if(buckets->empty[(current_bucket + i) % num_buckets] == 0)
/**/        break;
/**/    
/**/    if(i == num_buckets)
/**/      break;
/**/    current_bucket = (current_bucket + i) % num_buckets;
/**/    if(DPRT){printf("%02d: Starting inner loop: working on bucket %"PRId64"\n", MYTHREAD, current_bucket);}
/**/    //min_bucket = i + 1;
/**/    // inner loop
/**/    int64_t start = 0;
/**/    int64_t end = 0;
/**/    llnode_t * v = buckets->B[current_bucket];
/**/    
/**/    while(v != NULL){
/**/      if(DPRT){printf("%02d: Processing Node %"PRId64" in Bucket %"PRId64"\n", MYTHREAD, v->index, current_bucket);}
/**/
/**/      /* take v out of B[current_bucket]??? */
/**/      remove_node_from_bucket(v, current_bucket, buckets);
/**/      
/**/      /* relax light edges from v */
/**/      for(j = mat->loffset[v->index]; j < mat->loffset[v->index + 1]; j++){
/**/        if(mat->lvalue[j] <= delta){      
/**/          J = mat->lnonzero[j];
/**/          pe  = J % THREADS;
/**/          pkg.lj = J / THREADS;
/**/          pkg.tw = tent->lentry[v->index] + mat->lvalue[j];
/**/          //printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, pkg.i, J, ltent[li],  pkg.tw); 
/**/          if( exstack_push(ex, &pkg, pe) == 0 ) {
/**/            delta_exstack_relax_process(tent, ex, Buckets, 0);
/**/            j--;
/**/          }
/**/          //relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
/**/        }
/**/      }
/**/      
/**/      /* insert v into R if it is not already there */
/**/      if(deleted[v->index] == 0){
/**/        deleted[v->index] = 1;
/**/        R[end++] = v->index;
/**/        num_deleted++;
/**/        if(DPRT){printf("%02d: deleted %"PRId64"s\n", MYTHREAD, v->index);}
/**/      }
/**/      
/**/      //v = v->next;
/**/      v = buckets->B[current_bucket];
/**/    }// end inner loop
/**/    delta_exstack_relax_process(tent, ex, Buckets, 1);
/**/
/**/    /* relax heavy requests edges for everything in R */
/**/    while(start < end){
/**/      v = &(buckets->nodes[R[start++]]);
/**/      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
/**/        if(mat->value[j] > delta){
/**/          local_relax(mat->nonzero[j], tent->lentry[v->index] + mat->lvalue[j], tent, buckets);
/**/        }
/**/      }      
/**/    }
/**/    current_bucket++;
/**/    lgp_barrier();
/**/    exstack_reset(ex);
/**/  }// end main loop
/**/
/**/  dump_tent("  Delta Exstack:", tent);
/**/
/**/  free(buckets->B);
/**/  free(buckets->nodes);
/**/  free(buckets->empty);
/**/  free(buckets->level);
/**/  free(deleted);
/**/  free(R);
/**/  
/**/  return(wall_seconds() - tm);
/**/}
/**/
/**/
#endif

#include "sssp.h"

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

typedef struct all_ds_t {
  SHARED ds_t * Sds;
  ds_t * lds;
}all_ds_t;


typedef struct pkg_delta_e_t{
  int64_t i;  // tail of the edge (in case one wanted to make back pointers
  int64_t lj; // local "head" of the edge on remote pe
  double tw;  // candidate tentative weight
}pkg_delta_e_t;


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
  if(DPRT){printf("Adding %"PRId64" to bucket %"PRId64" of %"PRId64"\n", v, i_m, ds->num_buckets);}
  
  assert(i_m >= -1 && i_m < ds->num_buckets);
  
  if(ds->in_bucket[v] == i_m){      // it is ok if this node is already in this bucket
    return; 
  }
  assert(ds->in_bucket[v] == -1);   // better not be in a different bucket

  ds->next[v] = v;
  ds->prev[v] = v;
  if(ds->B[i_m] != -1){    
    w = ds->B[i_m];             // w is "first" on the list, insert v before w              
    if(DPRT){printf("non-empty: w=%ld, prev=%ld\n", w, ds->prev[w]);}
    ds->prev[v] = ds->prev[w];            
    ds->next[ds->prev[w]] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    if(DPRT){printf("   v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );}
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
  if(DPRT){printf("Removing %"PRId64" from bucket %"PRId64"\n", v, i_m);}
  if(DPRT && !( ((ds->next[v] == v) && (ds->prev[v] == v)) || (ds->next[v] != ds->prev[v]) )){
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
  if(DPRT){printf("relax head %"PRId64" cand_dist = %lf < %lf?\n", w, cand_dist, ds->tent[w]);}
  if ( cand_dist < ds->tent[w] ){
    iold = ds->in_bucket[w];
    inew = ((int64_t)floor(cand_dist/ds->delta)) % (ds->num_buckets);

    assert(iold >= -1 && iold < ds->num_buckets);
    assert(inew >= -1 && inew < ds->num_buckets);
    if(DPRT){printf("winner: %"PRId64"  move from bucket %"PRId64" to %"PRId64"\n", w, iold, inew);}
    if( iold != inew ){
      if(iold >= 0)
        remove_node_from_bucket_arr(ds, w);
      insert_node_in_bucket_arr(ds, w, inew);
    }
    ds->tent[w] = cand_dist;
  }
}

static int64_t delta_exstack_relax_process(ds_t *ds, exstack_t *ex, int64_t done) 
{
  int64_t fromth;
  pkg_delta_e_t pkg;

  exstack_exchange(ex);
  
  while(exstack_pop(ex, &pkg, &fromth)){
    relax_arr(ds, pkg.lj, pkg.tw);
  }
  return( exstack_proceed(ex, done) );
}


// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_exstack(d_array_t *dist, sparsemat_t * mat, int64_t r0)
{
  int64_t i, i_m, k;
  int64_t v;
  int64_t all_rbi, rbi;     // the real bucket index
  int64_t J, pe;
  pkg_delta_e_t pkg;


  //TODO: Fix the buffer size 
  exstack_t * ex = exstack_init(1024, sizeof(pkg_delta_e_t));
  if( ex == NULL) return(-1.0);
  double tm = wall_seconds();

  assert((r0 >= 0) && (r0<mat->numrows));
  
  int64_t * R = calloc(mat->lnumrows, sizeof(int64_t));
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  for(i = 0; i < mat->lnumrows; i++){
    if(max_degree < (mat->loffset[i+1] - mat->loffset[i]))
      max_degree = (mat->loffset[i+1] - mat->loffset[i]);
  }
  max_degree = lgp_reduce_max_l(max_degree);
  assert(max_degree > 0);
  delta = 1.0/max_degree;
  if(DPRT){printf("%02d:delta = %lf\n", MYTHREAD, delta);}
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->lnnz; i++)
    if(max_edge_weight < mat->lvalue[i])
      max_edge_weight = mat->lvalue[i];
  max_edge_weight = lgp_reduce_max_d(max_edge_weight);
  if(DPRT){printf("%02d:max edge weight = %lf\n", MYTHREAD, max_edge_weight);}

  if(DPRT){printf("%02d:Init Buckets\n", MYTHREAD);}
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  
  // Set up all the arrays to be local.
  // We will make a global index from local and pe if needed
  // allocate and initialize the all the arrays
  all_ds_t *DS = calloc(1,sizeof(all_ds_t));
  DS->Sds = lgp_all_alloc(THREADS, sizeof(ds_t));
  DS->lds = lgp_local_part(ds_t, DS->Sds);

  ds_t * ds = DS->lds;
  //ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  ds->next = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->deleted = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  //ds->tent = (double *)malloc(mat->lnumrows * sizeof(double)); assert(ds->tent != NULL);
  ds->tent = dist->lentry;
  for(i = 0; i < mat->lnumrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->tent[i]      = INFINITY;
  }
  if(DPRT){printf("Allocate buckets\n");}
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t)); assert(ds->B != NULL);
  // bucket indices are mod num_buckets.
  for(i_m = 0; i_m < num_buckets; i_m++){
    ds->B[i_m] = -1;
  }
  ds->delta = delta;

  // set the distance to r0 (as a global index) equal to 0.0
  if( (r0 % THREADS) == MYTHREAD) {
    r0 = r0/THREADS;
    ds->tent[r0] = 0.0;
    insert_node_in_bucket_arr(ds, r0, 0);
    if(DPRT){printf("%02d: Set source node %ld\n", MYTHREAD, r0);}
  }

  lgp_barrier();
  dump_tent("Delta Exstack Init:", dist);
  lgp_barrier();

  int64_t imdone, alldone;
  /* main loop */
  rbi = 0;   // my real bucket i
  all_rbi = 0;  // global real bucket i
  while(1){
    if(DPRT){printf("Outer loop rbi = %"PRId64"\n", rbi);}
    // find the minimum indexed non-empty bucket 
    for(i = 0; i < ds->num_buckets; i++){
      if(ds->B[(rbi + i) % ds->num_buckets] != -1)
        break;
    }
    if(i == num_buckets){  //All buckets are empty, we are done
      imdone = 1;
      alldone = lgp_reduce_min_l(imdone);
      if(alldone)
        break;
    }
    rbi = rbi + i;
    all_rbi = lgp_reduce_min_l(rbi);

    if( all_rbi == rbi ) {
      i_m = rbi % ds->num_buckets;
      if(DPRT){printf("%02d:Starting inner loop: working on bucket %"PRId64"\n",MYTHREAD, i_m);}
      if(DPRT) dump_bucket_arr(ds, i_m);

      // inner loop
      int64_t start = 0;
      int64_t end = 0;
      for( v=ds->B[i_m]; ds->in_bucket[v] == i_m; v=ds->next[v]){
        if(DPRT){printf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v, i_m);}

        remove_node_from_bucket_arr(ds, v);
        
        /* relax light edges from v */
        for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
          if(mat->lvalue[k] <= delta){	  
            J = mat->lvalue[k];
            pe  = J % THREADS;
            pkg.lj = J / THREADS;
            pkg.tw = ds->tent[v] + mat->lvalue[k];
            if(DPRT){printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", rbi, MYTHREAD, pkg.i, J, ds->tent[v],  pkg.tw);}
            if( exstack_push(ex, &pkg, pe) == 0 ) {
              delta_exstack_relax_process(ds, ex, 0);
              k--;
            }
          }
        } 
        
        /* insert v into R if it is not already there */
        if(ds->deleted[v] == 0){
          ds->deleted[v] = 1;
          R[end++] = v;
          if(DPRT){printf("%02d: deleted %"PRId64"\n", MYTHREAD, v);}
        }
      }

      /* relax heavy requests edges for everything in R */
      for(start=0; start<end; start++){
        v = R[start];
        for(k = mat->offset[v]; k < mat->offset[v + 1]; k++){
          if(mat->lvalue[k] > delta){	  
            J = mat->lvalue[k];
            pe  = J % THREADS;
            pkg.lj = J / THREADS;
            pkg.tw = ds->tent[v] + mat->lvalue[k];
            if(DPRT){printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", rbi, MYTHREAD, pkg.i, J, ds->tent[v],  pkg.tw);}
            if( exstack_push(ex, &pkg, pe) == 0 ) {
              delta_exstack_relax_process(ds, ex, 0);
              k--;
            }
          }
          //if(mat->value[k] > delta){
            //relax_arr(ds, mat->nonzero[k], ds->tent[v] + mat->value[k]);
          //}
        }      
      }
    }
    delta_exstack_relax_process(ds, ex, 1);
    lgp_barrier();
    exstack_reset(ex);
  }

  lgp_barrier();
  // Copy the answers to dist array
  //for(v = 0; v < mat->lnumrows; v++){
  //  dist->entry[v] = ds->tent[v];
  //}

  //free(ds->next);
  //free(ds->prev);
  //free(ds->in_bucket);
  //free(ds->deleted);
  //free(ds->tent);
  //free(ds);
  free(R);
  
  return(wall_seconds() - tm);
}
