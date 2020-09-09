/*! \file sssp_delta_exstack.upc
 * \brief Implementation of delta stepping using exstack
 */

#include "sssp.h"

#define DPRT 1

typedef struct pkg_delta_e_t{
  int64_t i;  // row[i] is "tail" of the edge, in case we want to set backpointers
  int64_t lj; // the local "head" of the edge
  double tw;  // new tentative weight
}pkg_delta_e_t;


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
  int64_t * in_bucket;
  int64_t num_buckets;
  double delta;
}buckets_t;

typedef struct all_buckets_t {
  SHARED buckets_t * Sbuckets;
  buckets_t * lbuckets;
}all_buckets_t;

// Remove a specific node from bucket i.
void remove_node_from_bucket(llnode_t * v, int64_t i, buckets_t * buckets){
  llnode_t * next = v->next;
  i = i % buckets->num_buckets;
  assert(buckets->in_bucket[v->index] == i);
  assert(buckets->empty[i] == 0);
  buckets->in_bucket[v->index] = -1;
  buckets->level[i] = -1;
  if(DPRT){printf("%02d: Removing %"PRId64" from bucket %"PRId64"\n",MYTHREAD, v->index, i);}
  if(v->prev)
    v->prev->next = next;
  else
    buckets->B[i] = v->next;
  if(next != NULL)
    next->prev = v->prev;
  if(buckets->B[i] == NULL)
    buckets->empty[i] = 1;
}

// Prepend a node into a bucket
void insert_node_in_bucket(llnode_t * w, int64_t i, buckets_t * buckets){
  int64_t actual_i = (i % buckets->num_buckets);

  if(DPRT){printf("%02d:Adding %"PRId64" to bucket %"PRId64"\n", MYTHREAD, w->index, actual_i);}
  
  if(buckets->in_bucket[w->index] == actual_i){
    //this node is already in this bucket
    if(buckets->level[actual_i] != (actual_i % buckets->num_buckets))
      if(DPRT){printf("%02d:%"PRId64" %"PRId64"\n", MYTHREAD, buckets->level[actual_i], (actual_i % buckets->num_buckets));}
    assert(buckets->level[actual_i] == (actual_i % buckets->num_buckets));
    return; 
  }

  if(buckets->in_bucket[w->index] != -1){
    printf("AAAhhh: %"PRId64"\n", buckets->in_bucket[w->index]);
  }
  assert(buckets->in_bucket[w->index] == -1);
  

  if(buckets->empty[actual_i] == 1){
    // this bucket is empty right now
    assert(buckets->B[actual_i] == NULL);
    buckets->empty[actual_i] = 0;
  }else{
    if(DPRT){printf("%02d:%"PRId64" %"PRId64"\n", MYTHREAD, buckets->level[actual_i], (actual_i % buckets->num_buckets));}
    assert(buckets->level[actual_i] == -1 || buckets->level[actual_i] == (actual_i % buckets->num_buckets));
  }
  
  buckets->level[actual_i] = (actual_i % buckets->num_buckets);
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
void move_node_from_bucket_i_to_j(llnode_t * w, int64_t i, int64_t j, buckets_t * buckets){

  if(DPRT){printf("%02d:Move (%"PRId64" which is in bucket %"PRId64") from %"PRId64" to %"PRId64"\n", MYTHREAD, w->index, buckets->in_bucket[w->index],i, j);}
  if(i >= 0){
    i = i % buckets->num_buckets;
    if(buckets->in_bucket[w->index] == i){
      
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
  }
  
  /* insert w into the front of B[j] */
  insert_node_in_bucket(w, j, buckets);
}

// relax an edge from a node to a node w (the current tenative distance to w is tent[windex]).
// the candidate distance to w is cand_dist.
static void local_relax(int64_t windex, double cand_dist, d_array_t * tent, buckets_t * buckets){
  
  if ( cand_dist < tent->lentry[windex] ){
    if(DPRT){printf("%02d:relax w=%"PRId64" cand_dist = %lf < %lf?\n", MYTHREAD, windex, cand_dist, tent->lentry[windex]);}
    /* if w is in B[floor(tent[windex]/delta)], remove it from that bucket */
    int64_t j;
    if(tent->lentry[windex] == INFINITY) 
      j = -1;
    else 
      j = (int64_t)floor(tent->lentry[windex]/buckets->delta);
    move_node_from_bucket_i_to_j(&(buckets->nodes[windex]), j, (int64_t)floor(cand_dist/buckets->delta), buckets);
    tent->lentry[windex] = cand_dist;
  }
}

/*!
 * \brief pop routine to implement relaxing the edges
 * \param tent pointer to the tentative distances array
 * \param *ex the extack buffers
 * \param done the signal to exstack_proceed that this thread is done
 * \return the return value from exstack_proceed
 */
static int64_t delta_exstack_relax_process(d_array_t *tent, exstack_t *ex, all_buckets_t *Buckets,  int64_t done) 
{
  buckets_t * buckets = Buckets->lbuckets;
  int64_t fromth;
  pkg_delta_e_t pkg;

  exstack_exchange(ex);
  
  while(exstack_pop(ex, &pkg, &fromth)){
    local_relax(pkg.lj, pkg.tw, tent, buckets);
  }
  return( exstack_proceed(ex, done) );
}



// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_exstack(d_array_t *tent, sparsemat_t * mat, int64_t r0)
{
  int64_t i, li, j;
  int64_t J, pe;
  pkg_delta_e_t pkg;
  int64_t lr0 = r0 / THREADS;   // the local row (vertex) name for r0

  //TODO: Fix the buffer size 
  exstack_t * ex = exstack_init(1024, sizeof(pkg_delta_e_t));
  if( ex == NULL) return(-1.0);

  double tm = wall_seconds();
  
  assert(r0 < mat->numrows);
  
  char * deleted = calloc(mat->lnumrows, sizeof(char));
  int64_t * R = calloc(mat->lnumrows, sizeof(int64_t));
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  for(li = 0; li < mat->lnumrows; li++)
    tent->lentry[i] = INFINITY;

  for(li = 0; li < mat->lnumrows; li++){
    if(max_degree < (mat->loffset[i+1] - mat->loffset[i]))
      max_degree = (mat->loffset[i+1] - mat->loffset[i]);
  }
  max_degree = lgp_reduce_max_l(max_degree);
  assert(max_degree > 0);
  delta = 1.0/max_degree;
  if(DPRT){printf("%02d:delta = %lf\n", MYTHREAD, delta);}
  
  double max_edge_weight = 0.0;
  for(li = 0; li < mat->lnnz; li++)
    if(max_edge_weight < mat->lvalue[li])
      max_edge_weight = mat->lvalue[li];
  max_edge_weight = lgp_reduce_max_d(max_edge_weight);

  if(DPRT){printf("%02d:max edge weight = %lf\n", MYTHREAD, max_edge_weight);}
  
  if(DPRT){printf("%02d:Init Buckets\n", MYTHREAD);}


  /* set up buckets as an array of linked lists */

  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;

#if 0
  /* set up buckets as an array of linked lists */
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  buckets_t * buckets = calloc(1, sizeof(buckets_t));
#endif
  all_buckets_t * Buckets= calloc(1, sizeof(all_buckets_t));
  Buckets->Sbuckets = lgp_all_alloc(THREADS, sizeof(buckets_t));
  Buckets->lbuckets = lgp_local_part(buckets_t, Buckets->Sbuckets);

  buckets_t * buckets = Buckets->lbuckets;

  buckets->num_buckets = num_buckets;
  buckets->B = calloc(num_buckets, sizeof(llnode_t *));
  buckets->nodes = calloc(mat->lnumrows, sizeof(llnode_t));
  buckets->in_bucket = calloc(mat->lnumrows, sizeof(int64_t));
  buckets->level = calloc(num_buckets, sizeof(int64_t));
  buckets->empty = calloc(num_buckets, sizeof(char));
  buckets->delta = delta;
  
  if(DPRT){printf("%02d:num_buckets = %"PRId64" delta = %lf max_degree = %"PRId64"\n", MYTHREAD, num_buckets, delta, max_degree);}
  for(i = 0; i < num_buckets; i++){
    buckets->B[i] = NULL;
    buckets->empty[i] = 1;
    buckets->level[i] = -1;
  }
  for(i = 0; i < mat->lnumrows; i++){
    buckets->nodes[i].index = i;
    buckets->nodes[i].next = NULL;
    buckets->nodes[i].prev = NULL;
    buckets->in_bucket[i] = -1;
  }

  
  /* set the source distance to 0 */
  if( (r0 % THREADS) == MYTHREAD) {
    insert_node_in_bucket(&(buckets->nodes[lr0]), 0, buckets);
    tent->lentry[lr0] = 0.0;
  }

  lgp_barrier();

  /* main loop */
  //int64_t min_bucket = 0;
  int64_t current_bucket = 0;
  int64_t num_deleted = 0;
  while(1){

    /* find the minimum indexed non-empty bucket */
    for(i = 0; i < num_buckets; i++)
      if(buckets->empty[(current_bucket + i) % num_buckets] == 0)
        break;
    
    if(i == num_buckets)
      break;
    current_bucket = (current_bucket + i) % num_buckets;
    if(DPRT){printf("%02d: Starting inner loop: working on bucket %"PRId64"\n", MYTHREAD, current_bucket);}
    //min_bucket = i + 1;
    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    llnode_t * v = buckets->B[current_bucket];
    
    while(v != NULL){
      if(DPRT){printf("%02d: Processing Node %"PRId64" in Bucket %"PRId64"\n", MYTHREAD, v->index, current_bucket);}

      /* take v out of B[current_bucket]??? */
      remove_node_from_bucket(v, current_bucket, buckets);
      
      /* relax light edges from v */
      for(j = mat->loffset[v->index]; j < mat->loffset[v->index + 1]; j++){
        if(mat->lvalue[j] <= delta){      
          J = mat->lnonzero[j];
          pe  = J % THREADS;
          pkg.lj = J / THREADS;
          pkg.tw = tent->lentry[v->index] + mat->lvalue[j];
          //printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, pkg.i, J, ltent[li],  pkg.tw); 
          if( exstack_push(ex, &pkg, pe) == 0 ) {
            delta_exstack_relax_process(tent, ex, Buckets, 0);
            j--;
          }
          //relax(mat->nonzero[j], tent[v->index] + mat->value[j], tent, buckets);
        }
      }
      
      /* insert v into R if it is not already there */
      if(deleted[v->index] == 0){
        deleted[v->index] = 1;
        R[end++] = v->index;
        num_deleted++;
        if(DPRT){printf("%02d: deleted %"PRId64"s\n", MYTHREAD, v->index);}
      }
      
      //v = v->next;
      v = buckets->B[current_bucket];
    }// end inner loop
    delta_exstack_relax_process(tent, ex, Buckets, 1);

    /* relax heavy requests edges for everything in R */
    while(start < end){
      v = &(buckets->nodes[R[start++]]);
      for(j = mat->offset[v->index]; j < mat->offset[v->index + 1]; j++){
        if(mat->value[j] > delta){
          local_relax(mat->nonzero[j], tent->lentry[v->index] + mat->lvalue[j], tent, buckets);
        }
      }      
    }
    current_bucket++;
    lgp_barrier();
    exstack_reset(ex);
  }// end main loop
  free(buckets->B);
  free(buckets->nodes);
  free(buckets->empty);
  free(buckets->level);
  free(deleted);
  free(R);
  
  return(wall_seconds() - tm);
}
