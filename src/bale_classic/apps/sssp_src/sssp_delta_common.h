/*! \file sssp_delta_convey.upc
 * \brief Implementation of delta stepping using convey
 */

#ifndef sssp_delta_common_INCLUDED
#define sssp_delta_common_INCLUDED

#define DPRT 0
#define D0PRT (!MYTHREAD && DPRT)

typedef struct ds_t{
  int64_t *next;        // next in linked list
  int64_t *prev;        // prev in linked list
  int64_t *in_bucket;   // which bucket it is in, else -1
  int64_t *B;           // array of Buckets: B[i] is the index of the first node on the list, or -1 if empty
  int64_t num_buckets;  //
  double  *tent;        // the tentative weight for the vertex
  int64_t *deleted;     // deleted means resolved?
  int64_t *R;           // queue to hold tail vertices that need to relax their heavy edges
  double  delta;        // 
}ds_t;

void dump_bucket(ds_t *ds, int64_t i_m);
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m);
void remove_node_from_bucket(ds_t *ds, int64_t v);
void local_relax(ds_t *ds, int64_t w, double cand_dist);


void calculate_delta_and_num_buckets(double *delta, int64_t *num_buckets, sparsemat_t *mat);
void allocate_and_initialize_delta_stepping_struct(ds_t *ds, int64_t lnumrows, int64_t num_buckets, double delta);
void clear_ds_struct(ds_t *ds);

#endif  // sssp_delta_common_INCLUDED
