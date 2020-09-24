/*! \file sssp_delta_exstack.upc
 * \brief Implementation of delta stepping using exstack
 */

#include "sssp.h"
#include "sssp_delta_common.h"


static int64_t delta_exstack_relax_process(ds_t *ds, exstack_t *ex, int64_t done) 
{
  int64_t fromth;
  sssp_pkg_t pkg;

  exstack_exchange(ex);
  while(exstack_pop(ex, &pkg, &fromth)){
    local_relax(ds, pkg.lj, pkg.tw);
  }
  return( exstack_proceed(ex, done) );
}

int64_t delta_exstack_push(exstack_t *ex, int64_t gidx, double tent_wt)
{
  struct sssp_pkg_t pkg;
  int64_t pe; 
  pe = gidx % THREADS;
  pkg.lj  = gidx / THREADS;
  pkg.tw = tent_wt;
  return(exstack_push(ex, &pkg, pe));
}

// This is the delta stepping algorithm as it appears in
// the paper "Delta-stepping: a parallelizable shortest path algorithm" by
// U. Meyer and P. Sanders.

double sssp_delta_exstack(d_array_t *dist, sparsemat_t * mat, int64_t r0, double opt_delta)
{
  int64_t i, i_m, k;
  int64_t v;
  int64_t J, pe;
  struct sssp_pkg_t pkg;


  //TODO: Fix the buffer size 
  exstack_t * ex = exstack_init(32, sizeof(sssp_pkg_t));
  if( ex == NULL) return(-1.0);
  double tm = wall_seconds();

  assert((r0 >= 0) && (r0<mat->numrows));
  
  double delta;
  int64_t num_buckets;
  calculate_delta_and_num_buckets(&delta, &num_buckets, mat, opt_delta);
  
  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  allocate_and_initialize_delta_stepping_struct(ds, mat->lnumrows, num_buckets, delta);

  sparsemat_t *light = get_light_edges(mat, delta);
  sparsemat_t *heavy = get_heavy_edges(mat, delta);
  write_matrix_mm(mat, "mat_mat");
  write_matrix_mm(light, "mat_light");
  write_matrix_mm(heavy, "mat_heavy");
  lgp_barrier();

  if( (r0 % THREADS) == MYTHREAD) {    // set the distance to r0 (as a global index) equal to 0.0
    r0 = r0/THREADS;
    ds->tent[r0] = 0.0;
    insert_node_in_bucket(ds, r0, 0);
  }
  lgp_barrier();

  int64_t grmb = 0;    // the global real min bucket (collectively the min non-empty bucket)
                       // we use this to handle the problem of min bucket (mod num_buckets)
  int64_t all_i = 0;
  int64_t start = 0;
  int64_t end = 0;

  // MAIN LOOP: exists when all buckets are empty
  while(1){      
    for(i = 0; i < ds->num_buckets; i++){          // find the local min (done mod num_buckets) non-empty bucket 
      if(ds->B[(grmb + i) % ds->num_buckets] != -1)
        break;
    }
    all_i = lgp_reduce_min_l(i);                   // global min index
    if(all_i == num_buckets)                       // All buckets are empty, we are done
        break;
    grmb = grmb + all_i;
    i_m = grmb % ds->num_buckets;

    start = 0;
    end = 0;
    while( lgp_reduce_max_l(ds->B[i_m]) > -1){    // repeatly check that B[i_m] on all threads are empty

      while( ds->B[i_m] >= 0 ) {                  // removes vertices from this bucket, but the relaxation process
                                                  // may add vertices to this bucket on any thread
        v = ds->B[i_m]; 
        remove_node_from_bucket(ds, v);

#if 0
        for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){        // relax light edges from v 
          if(mat->lvalue[k] <= delta){	  
            if( delta_exstack_push(ex, mat->lnonzero[k],  ds->tent[v] + mat->lvalue[k]) == 0 ) {
              delta_exstack_relax_process(ds, ex, 0);
              k--;
            }
          }
        } 
#else 
        for(k = light->loffset[v]; k < light->loffset[v + 1]; k++){        // relax light edges from v 
          if( delta_exstack_push(ex, light->lnonzero[k],  ds->tent[v] + light->lvalue[k]) == 0 ) {
            delta_exstack_relax_process(ds, ex, 0);
            k--;
          }
        } 
#endif
        if(ds->deleted[v] == 0){  // insert v into R if it is not already there
          ds->deleted[v] = 1;
          ds->R[end++] = v;
        }
      }
      while( delta_exstack_relax_process(ds, ex, 1)) 
        ;
      lgp_barrier();
      exstack_reset(ex);
    } 
    lgp_barrier();

    for(start=0; start<end; start++){           // relax heavy requests edges for everything in R 
      v = ds->R[start];
#if 0
      for(k = mat->loffset[v]; k < mat->loffset[v + 1]; k++){
        if(mat->lvalue[k] > delta){	  
          if( delta_exstack_push(ex, mat->lnonzero[k],  ds->tent[v] + mat->lvalue[k]) == 0 ) {
            delta_exstack_relax_process(ds, ex, 0);
            k--;
          }
        }
      }
#else
      for(k = heavy->loffset[v]; k < heavy->loffset[v + 1]; k++){
        if( delta_exstack_push(ex, heavy->lnonzero[k],  ds->tent[v] + heavy->lvalue[k]) == 0 ) {
          delta_exstack_relax_process(ds, ex, 0);
          k--;
        }
      }
#endif
    }
    while( delta_exstack_relax_process(ds, ex, 1) )
      ;
    exstack_reset(ex);
  }

  lgp_barrier();
  for(v = 0; v < mat->lnumrows; v++){              //Copy the answers to dist array
    dist->lentry[v] = ds->tent[v];
  }

  clear_ds_struct(ds); free(ds);

  return(wall_seconds() - tm);
}
