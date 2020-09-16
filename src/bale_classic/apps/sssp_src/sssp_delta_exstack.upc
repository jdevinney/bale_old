/*! \file sssp_delta_exstack.upc
 * \brief Implementation of delta stepping using exstack
 */

#include "sssp.h"

#define DPRT 0
#define D0PRT (!MYTHREAD && DPRT)

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

typedef struct pkg_delta_e_t{
  int64_t i;  // tail of the edge (in case one wanted to make back pointers
  int64_t lj; // local "head" of the edge on remote pe
  double tw;  // candidate tentative weight
}pkg_delta_e_t;


// debugging function
void dump_bucket_arr(ds_t *ds, int64_t i_m)
{
  int64_t v, w ;
  char lineout[2048];

  lineout[0] = '\0';
  sprintf(lineout+strlen(lineout), "%02d: Bucket[%ld] =", MYTHREAD, i_m);
  if( ds->B[i_m] == -1 ){
    printf("%s empty \n",lineout);
    return;
  }
  v = w = ds->B[i_m];
  do {
    sprintf(lineout+strlen(lineout), "%ld ", w);
  } while( (w = ds->next[w]) != v );
  printf("%s\n",lineout);
  return;
}

// Prepend node v into a bucket i
void insert_node_in_bucket_arr(ds_t *ds, int64_t v, int64_t i_m)
{
  int64_t w;                   // node on list, given by ds->B[i_m]
  //if(DPRT){printf("%02d: Adding %"PRId64" to bucket %"PRId64" of %"PRId64"\n", MYTHREAD, v, i_m, ds->num_buckets);}
  
  assert(i_m >= -1 && i_m < ds->num_buckets);
  
  if(ds->in_bucket[v] == i_m){      // it is ok if this node is already in this bucket
    return; 
  }
  assert(ds->in_bucket[v] == -1);   // better not be in a different bucket

  ds->next[v] = v;
  ds->prev[v] = v;
  if(ds->B[i_m] != -1){    
    w = ds->B[i_m];             // w is "first" on the list, insert v before w              
    //if(DPRT){printf("%02d: non-empty: w=%ld, prev=%ld\n", MYTHREAD, w, ds->prev[w]);}
    ds->prev[v] = ds->prev[w];            
    ds->next[ds->prev[w]] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    //if(DPRT){printf("%02d:    v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", MYTHREAD, v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );}
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
  //if(DPRT){printf("%02d: Removing %"PRId64" from bucket %"PRId64"\n", MYTHREAD, v, i_m);}
  if(00 && DPRT && !( ((ds->next[v] == v) && (ds->prev[v] == v)) || (ds->next[v] != ds->prev[v]) )){
    printf("%02d: ERROR:  v, next, prev = %ld %ld %ld\n", MYTHREAD, v, ds->next[v], ds->prev[v] );
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
  //if(DPRT){printf("%02d: relax head %"PRId64" cand_dist = %lf < %lf?\n", MYTHREAD, w, cand_dist, ds->tent[w]);}
  if ( cand_dist < ds->tent[w] ){
    iold = ds->in_bucket[w];
    inew = ((int64_t)floor(cand_dist/ds->delta)) % (ds->num_buckets);

    assert(iold >= -1 && iold < ds->num_buckets);
    assert(inew >= -1 && inew < ds->num_buckets);
    //if(DPRT){printf("%02d: winner: %"PRId64"  move from bucket %"PRId64" to %"PRId64"\n", MYTHREAD, w, iold, inew);}
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
  //if(DPRT){printf("%02d: exstack relax process\n",MYTHREAD);}

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
  int64_t rbi;     // the real bucket index
  int64_t J, pe;
  pkg_delta_e_t pkg;



  //TODO: Fix the buffer size 
  exstack_t * ex = exstack_init(32, sizeof(pkg_delta_e_t));
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
  if(D0PRT){printf("%02d:delta = %lf\n", MYTHREAD, delta);}
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->lnnz; i++)
    if(max_edge_weight < mat->lvalue[i])
      max_edge_weight = mat->lvalue[i];
  max_edge_weight = lgp_reduce_max_d(max_edge_weight);
  if(D0PRT){printf("%02d:max edge weight = %lf\n", MYTHREAD, max_edge_weight);}

  if(D0PRT){printf("%02d:Init Buckets\n", MYTHREAD);}
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  
  // Set up all the ds arrays to be local, ie. shared nothing.
  // The only way to affect change is local read/writes and exstack pushes/pops.
  // We will make a global for vertices to match the rows in the matrix 
  // from the local index and the pe if needed.

  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  ds->next = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->deleted = (int64_t *)malloc(mat->lnumrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  ds->tent = (double *)malloc(mat->lnumrows * sizeof(double)); assert(ds->tent != NULL);
  for(i = 0; i < mat->lnumrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->tent[i]      = INFINITY;
  }
  if(D0PRT){printf("Allocate buckets\n");}
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
  //for(v = 0; v < mat->lnumrows; v++){ dist->lentry[v] = ds->tent[v]; }
  //lgp_barrier();
  //dump_tent("Delta Exstack Init:", dist);
  //lgp_barrier();

  /* main loop */
  rbi = 0;   // collectively, the real min bucket 
  int64_t all_i = 0;
  while(1){
    if(DPRT){printf("%02d: Outer loop rbi = %"PRId64"\n", MYTHREAD, rbi);}

    // find the minimum indexed non-empty bucket 
    for(i = 0; i < ds->num_buckets; i++){
      if(ds->B[(rbi + i) % ds->num_buckets] != -1)
        break;
    }
    if(DPRT){printf("%02d: my min bucket %ld\n", MYTHREAD, rbi+i);}

    all_i = lgp_reduce_min_l(i);

    //if(all_i == 0){ // B[rbi] is empty on all threads
    //  rbi = rbi + 1;
    //  continue;
    //}

    if(all_i == num_buckets){  //All buckets are empty, we are done
        break;
    }
    rbi = rbi + all_i;

    i_m = rbi % ds->num_buckets;
    if(DPRT){printf("%02d:Starting inner loop: working on bucket %ld, %ld\n",MYTHREAD, rbi, i_m);}
    if(DPRT) dump_bucket_arr(ds, i_m);

    // inner loop
    int64_t start = 0;
    int64_t end = 0;
    for( v=ds->B[i_m]; v>=0 && ds->in_bucket[v] == i_m; v=ds->next[v]){
      if(DPRT){printf("%02d: Processing Node %"PRId64" in Bucket %"PRId64"\n",MYTHREAD, v, i_m);}

      remove_node_from_bucket_arr(ds, v);
      
      /* relax light edges from v */
      pkg.i = v*THREADS + MYTHREAD; 
      if(0&&DPRT){printf("%02d: v=%ld has degree %ld\n", MYTHREAD, pkg.i, mat->loffset[v+1]-mat->loffset[v]);}
      for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
        if(mat->lvalue[k] <= delta){	  
          J = mat->lnonzero[k];
          pe  = J % THREADS;
          pkg.lj = J / THREADS;
          pkg.tw = ds->tent[v] + mat->lvalue[k];
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
        //if(DPRT){printf("%02d: deleted %"PRId64"\n", MYTHREAD, v);}
      }
    }

    /* relax heavy requests edges for everything in R */
    for(start=0; start<end; start++){
      v = R[start];
      pkg.i = v*THREADS + MYTHREAD; 
      for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
        if(mat->lvalue[k] > delta){	  
          J = mat->lnonzero[k];
          pe  = J % THREADS;
          pkg.lj = J / THREADS;
          pkg.tw = ds->tent[v] + mat->lvalue[k];
          if( exstack_push(ex, &pkg, pe) == 0 ) {
            delta_exstack_relax_process(ds, ex, 0);
            k--;
          }
        }
      }
    }
    if(DPRT){printf("%02d: Finishing inner loop: rbi %ld\n",MYTHREAD, rbi);}
    while( delta_exstack_relax_process(ds, ex, 1) )
      ;
    lgp_barrier();
    exstack_reset(ex);
  }

  lgp_barrier();
  //Copy the answers to dist array
  for(v = 0; v < mat->lnumrows; v++){
    dist->lentry[v] = ds->tent[v];
  }
  //dump_tent("Delta Exstack Done:", dist);

  free(ds->next);
  free(ds->prev);
  free(ds->in_bucket);
  free(ds->deleted);
  free(ds->tent);
  free(ds);
  free(R);
  
  return(wall_seconds() - tm);
}
