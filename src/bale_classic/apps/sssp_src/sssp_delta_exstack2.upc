/*! \file sssp_delta_exstack2.upc
 * \brief Implementation of delta stepping using exstack2
 */

#include "sssp.h"
#include "sssp_delta_common.h"

static int64_t delta_exstack2_relax_process(ds_t *ds, exstack2_t *ex2, int64_t done) 
{
  int64_t fromth;
  sssp_pkg_t pkg;
  //if(DPRT){printf("%02d: exstack2 relax process\n",MYTHREAD);}

  while(exstack2_pop(ex2, &pkg, &fromth)){
    local_relax(ds, pkg.lj, pkg.tw);
  }
  return( exstack2_proceed(ex2, done) );
}


// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_exstack2(d_array_t *dist, sparsemat_t * mat, int64_t r0)
{
  int64_t i, i_m, k;
  int64_t v;
  int64_t rbi;     // the real bucket index
  int64_t J, pe;
  sssp_pkg_t pkg;



  //TODO: Fix the buffer size 
  exstack2_t * ex2 = exstack2_init(64, sizeof(sssp_pkg_t));
  if( ex2 == NULL) return(-1.0);

  double tm = wall_seconds();

  assert((r0 >= 0) && (r0<mat->numrows));
  
  double delta;
  int64_t num_buckets;
  calculate_delta_and_num_buckets(&delta, &num_buckets, mat);

  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  allocate_and_initialize_delta_stepping_struct(ds, mat->lnumrows, num_buckets, delta);

  // set the distance to r0 (as a global index) equal to 0.0
  if( (r0 % THREADS) == MYTHREAD) {
    r0 = r0/THREADS;
    ds->tent[r0] = 0.0;
    insert_node_in_bucket(ds, r0, 0);
    if(DPRT){printf("%02d: Set source node %ld\n", MYTHREAD, r0);}
  }

  lgp_barrier();

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

    if(all_i == num_buckets){  //All buckets are empty, we are done
        break;
    }
    rbi = rbi + all_i;

    i_m = rbi % ds->num_buckets;
    if(DPRT){printf("%02d:Starting inner loop: working on bucket %ld, %ld\n",MYTHREAD, rbi, i_m);}
    if(DPRT) dump_bucket(ds, i_m);

    // inner loop
    int64_t start = 0;
    int64_t end = 0;

    while( lgp_reduce_max_l(ds->B[i_m]) > -1){
      while( ds->B[i_m] >= 0 ) {
        v = ds->B[i_m]; 
        if(DPRT){printf("%02d: Processing Node %"PRId64" in Bucket %"PRId64"\n",MYTHREAD, v, i_m);}

        remove_node_from_bucket(ds, v);
        /* relax light edges from v */
        if(0&&DPRT){printf("%02d: v=%ld has degree %ld\n", MYTHREAD, v*THREADS + MYTHREAD, mat->loffset[v+1]-mat->loffset[v]);}
        for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
          if(mat->lvalue[k] <= delta){	  
            J = mat->lnonzero[k];
            pe  = J % THREADS;
            pkg.lj = J / THREADS;
            pkg.tw = ds->tent[v] + mat->lvalue[k];
            if( exstack2_push(ex2, &pkg, pe) == 0 ) {
              delta_exstack2_relax_process(ds, ex2, 0);
              k--;
            }
          }
        } 
        if(ds->deleted[v] == 0){  // insert v into R if it is not already there
          ds->deleted[v] = 1;
          ds->R[end++] = v;
          //if(DPRT){printf("%02d: deleted %"PRId64"\n", MYTHREAD, v);}
        }
      }
      while( delta_exstack2_relax_process(ds, ex2, 1)) 
        ;
      lgp_barrier();
      exstack2_reset(ex2);
    } 
    lgp_barrier();

    /* relax heavy requests edges for everything in R */
    for(start=0; start<end; start++){
      v = ds->R[start];
      for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
        if(mat->lvalue[k] > delta){	  
          J = mat->lnonzero[k];
          pe  = J % THREADS;
          pkg.lj = J / THREADS;
          pkg.tw = ds->tent[v] + mat->lvalue[k];
          if( exstack2_push(ex2, &pkg, pe) == 0 ) {
            delta_exstack2_relax_process(ds, ex2, 0);
            k--;
          }
        }
      }
    }
    if(DPRT){printf("%02d: Finishing inner loop: rbi %ld\n",MYTHREAD, rbi);}
    while( delta_exstack2_relax_process(ds, ex2, 1) )
      ;
    lgp_barrier();
    exstack2_reset(ex2);
    lgp_barrier();
  }

  lgp_barrier();
  //Copy the answers to dist array
  for(v = 0; v < mat->lnumrows; v++){
    dist->lentry[v] = ds->tent[v];
  }
  //dump_tent("Delta exstack2 Done:", dist);

  clear_ds_struct(ds); free(ds);
  
  return(wall_seconds() - tm);
}
