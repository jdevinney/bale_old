/*! \file sssp_delta_exstack.upc
 * \brief Implementation of delta stepping using exstack
 */

#include "sssp.h"
#include "sssp_delta_common.h"


// debugging function
void dump_bucket(ds_t *ds, int64_t i_m)
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
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m)
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
void remove_node_from_bucket(ds_t *ds, int64_t v)
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
  ds->B[i_m] = w;         //move the B[i] pointer to w, 
  return;
}

// relax an edge to the head vertex, given the new tentative distance
// (= the tentative distance to the tail plus the weight of the edge).
// the candidate distance to w is cand_dist.
void local_relax(ds_t *ds, int64_t w, double cand_dist)
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
        remove_node_from_bucket(ds, w);
      insert_node_in_bucket(ds, w, inew);
    }
    ds->tent[w] = cand_dist;
  }
}

void allocate_and_initialize_delta_stepping_struct(ds_t *ds, int64_t lnumrows, int64_t num_buckets, double delta)
{
  int64_t i;
  ds->next = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->deleted = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  ds->R = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->R != NULL);
  ds->tent = (double *)malloc(lnumrows * sizeof(double)); assert(ds->tent != NULL);
  for(i = 0; i < lnumrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->tent[i]      = INFINITY;
  }
  if(D0PRT){printf("Allocate buckets\n");}
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t)); assert(ds->B != NULL);

  for(i = 0; i < num_buckets; i++){
    ds->B[i] = -1;
  }
  ds->delta = delta;
}

void clear_ds_struct(ds_t *ds)
{
  free(ds->next);
  free(ds->prev);
  free(ds->in_bucket);
  free(ds->deleted);
  free(ds->R);
  free(ds->tent);
}

