/******************************************************************
//
//
//  Copyright(C) 2018, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 *****************************************************************/
/*! \file triangle.upc
 * \brief Demo application that counts triangles in a graph.
 */

#include "triangle.h"
#include <spmat_opts.h>
#include <std_options.h>
/*!
  \page triangles_page Triangles

This uses matrix algebra approach to counting triangles in a graph.

The adjacency matrix, <b>A</b>, for the graph is a {0,1} matrix
where the rows and cols correspond to the vertices
and \f$a_{ij} \f$ = <b>A[i][j]</b> is 1 exactly when there is a edge between 
vertices <i>v_i</i> and <i>v_j</i>.

The triangle with vertices <i>{v_i, v_j, v_k}</i> has associated 
edges <i>{v_i, v_j}</i>, <i>{v_j, v_k}</i> and <i>{v_k, v_i}</i> 
which correspond to non-zero entries 
\f$a_{ij}\f$,
\f$a_{jk}\f$, and
\f$a_{ki}\f$
in the adjacency matrix.  
Hence the sum
\f$ \sum_{i,j,k} a_{ij}a_{jk}a_{ki} \f$ counts the triangles in the graph.
However, it counts each triangle 6 times according to the 6 symmetries of a triangle
and the 6 symmetric ways to choose the three nonzeros in <b>A</b>.
To count each triangle once, we compute the sum
\f[ \sum_{i=1}^{n}\sum_{j=1}^{i-1}\sum_{k=1}^{j-1} a_{ij}a_{jk}a_{ik} = 
    \sum_{i=1}^{n}\sum_{j=1}^{i-1} a_{ij} \sum_{k=1}^{j-1} a_{jk}a_{ik}. \f]

This picks out a unique labelling from the 6 possible and it means that 
all the information we need about edges is contained in the lower triangular 
part of symmetric adjacency matrix.  We call this matrix <b>L</b>.

The mathematical expression: 
for each nonzero \f$ a_{ij} \f$ compute the dot product 
of row \f$ i\f$ and row \f$ j \f$ becomes
\verbatim
  For each non-zero L[i][j] 
     compute the size of the intersection of the nonzeros in row i and row j
\endverbatim

Kronecker product graphs were implemented specifically for this app.
See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al.
for more information on Kronecker product graphs. 
 */


typedef struct args_t{
  int alg;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'a': args->alg = atoi(arg); break;     
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"triangle_alg", 'a', "ALG", 0, "Algorithm: 0 means L&L*U, 1 means L&U*L"},  
    {0}
  };

static struct argp_child children_parsers[] =
  {    
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };



int main(int argc, char * argv[]) {

  lgp_init(argc, argv);

  double t1;
  int64_t i, j;

  /* process command line */
  int ret = 0;
  args_t args;
  struct argp argp = {options, parse_opt, 0,
                      "Parallel sparse matrix transpose.", children_parsers};
  args.alg = 0;
  if(MYTHREAD == 0){
    ret = argp_parse(&argp, argc, argv, ARGP_NO_EXIT, 0, &args);
  }
  ret = distribute_cmd_line(argc, argv, &args, sizeof(args_t), ret);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(!MYTHREAD && !args.std.quiet){
    T0_fprintf(stderr,"Running on %d PEs\n", THREADS);
    write_std_graph_options(&args.gstd);
    write_std_options(&args.std);
  }
  
  // read in a matrix or generate a random graph
  sparsemat_t * L = get_input_graph(&args.std, &args.gstd);
  if(!L){T0_fprintf(stderr, "ERROR: transpose: L is NULL!\n");return(-1);}

  if(!is_lower_triangular(L, 0)){
    if(args.gstd.readfile){
      tril(L, -1);
    }else{
      T0_fprintf(stderr,"ERROR: L is not lower triangular!\n");
      lgp_global_exit(1);
    }
  }  
  
  lgp_barrier();
  
  /* calculate the number of triangles */
  double correct_answer = -1;
  if(args.gstd.model == KRONECKER){
    correct_answer = calculate_num_triangles(args.gstd.kron_mode, args.gstd.kron_spec, args.gstd.kron_num);
    //T0_fprintf(stderr, "Pre-calculated answer = %"PRId64"\n", (int64_t)correct_answer);
  }
  
  sparsemat_t * U;
  if(args.alg == 1)
    U = transpose_matrix(L);

  lgp_barrier();

  int64_t tri_cnt;           // partial count of triangles on this thread
  int64_t total_tri_cnt;     // the total number of triangles on all threads
  int64_t sh_refs;         // number of shared reference or pushes
  int64_t total_sh_refs;
  SHARED int64_t * cc = lgp_all_alloc(L->numrows, sizeof(int64_t));
  int64_t * l_cc = lgp_local_part(int64_t, cc);
  for(i = 0; i < L->lnumrows; i++)
    l_cc[i] = 0;
  lgp_barrier();
  
  /* calculate col sums */
  for(i = 0; i < L->lnnz; i++){
    lgp_fetch_and_inc(cc, L->lnonzero[i]);
  }
  
  lgp_barrier();
  
  int64_t rtimesc_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];        
    rtimesc_calc += deg*l_cc[i];
  }

  /* calculate sum (r_i choose 2) */
  int64_t rchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];
    rchoose2_calc += deg*(deg-1)/2;
  }
  
  /* calculate sum (c_i choose 2) */
  int64_t cchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = l_cc[i];
    cchoose2_calc += deg*(deg-1)/2;
  }
  int64_t pulls_calc = 0;
  int64_t pushes_calc = 0;
  if(args.alg == 0){
    pulls_calc = lgp_reduce_add_l(rtimesc_calc);
    pushes_calc = lgp_reduce_add_l(rchoose2_calc);
  }else{
    pushes_calc = lgp_reduce_add_l(rtimesc_calc);
    pulls_calc = lgp_reduce_add_l(cchoose2_calc);
  }

  lgp_all_free(cc);
  
  T0_fprintf(stderr,"Calculated: Pulls = %"PRId64"\n            Pushes = %"PRId64"\n\n",pulls_calc, pushes_calc);
  
  int64_t use_model;
  double laptime = 0.0;
  
  for( use_model=1L; use_model < 32; use_model *=2 ) {

    tri_cnt = 0;
    total_tri_cnt = 0;
    sh_refs = 0;
    total_sh_refs = 0;

    switch( use_model & args.std.models_mask ) {
    case AGI_Model:
      T0_fprintf(stderr,"      AGI: ");
      laptime = triangle_agi(&tri_cnt, &sh_refs, L, U, args.alg); 
      break;
    
    case EXSTACK_Model:
      T0_fprintf(stderr,"  Exstack: ");
      laptime = triangle_exstack_push(&tri_cnt, &sh_refs, L, U, args.alg, args.std.buffer_size);
      break;

    case EXSTACK2_Model:
      T0_fprintf(stderr," Exstack2: ");
      laptime = triangle_exstack2_push(&tri_cnt, &sh_refs, L, U, args.alg, args.std.buffer_size);
      break;

    case CONVEYOR_Model:
      T0_fprintf(stderr," Conveyor: ");
      laptime = triangle_convey_push(&tri_cnt, &sh_refs, L, U, args.alg);
      break;

    case ALTERNATE_Model:
      T0_fprintf(stderr,"ALTERNATE: ");      
      laptime = triangle_agi_iter(&tri_cnt, &sh_refs, L, U, args.alg);
      break;
    case 0:
      continue;
    }
    
    lgp_barrier();
    total_tri_cnt = lgp_reduce_add_l(tri_cnt);
    total_sh_refs = lgp_reduce_add_l(sh_refs);
    T0_fprintf(stderr,"  %8.3lf seconds: %16"PRId64" triangles", laptime, total_tri_cnt);
    T0_fprintf(stderr,"%16"PRId64" shared refs\n", total_sh_refs);
    if((correct_answer >= 0) && (total_tri_cnt != (int64_t)correct_answer)){
      T0_fprintf(stderr, "ERROR: Wrong answer!\n");
    }
    
    if(correct_answer == -1)
      correct_answer = total_tri_cnt;
    
  }
  
  lgp_barrier();
  lgp_finalize();
  return(0);
}
