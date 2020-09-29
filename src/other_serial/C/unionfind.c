/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/


/*! \file unionfind.c
 * \brief Program that uses a unionfind data structure to find connected components in a graph
 *
 * Run unionfind --help or --usage for insructions on running.
 */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"

/*! \page unionfind_page 
 Demo that uses the unionfind data structure
*/

/*! \brief Generate a matrix for the lower triangular section of the adjacency matrix of a Erdös–Renyi random graph.
 * \param numrows the number of rows (and columns) in the produced matrix
 * \param er_prob the probability that a particular edge exists in the graph
 * \param seed the seed for random number generator that determines the original matrix and the permutations
 * \param dump_files is a debugging flag
 * \return the lower triangular matrix
 */
sparsemat_t * generate_concomp_input(int64_t numrows, double edge_prob,  uint32_t seed, int dump_files) 
{
  //sparsemat_t * mat = erdos_renyi_tri(numrows, edge_prob, ER_TRI_L, seed);
  sparsemat_t * mat = random_graph(numrows, FLAT, UNDIRECTED, NOLOOPS, edge_prob, seed);
  if(!mat){
    fprintf(stderr, "ERROR: generate_concomp_iinput: failed!\n");
    return(NULL);
  }

  // What if we permute it here
  spmat_stats(mat);

  if(dump_files) dump_matrix(mat, 20, "orig_uf.out");

  return( mat );
}

/*! \struct comp_tree_t 
 * \brief the nodes of the disjoint union trees
 */
typedef struct comp_tree_t {
  int64_t parent;  //!< pointer to the nodes parent in the tree
  int64_t rank;    //!< if the node is the root, size of the component
} comp_tree_t;

/*! \brief Find the name (the root of the parent tree) for the component containing a given node.
 * \param cc the component tree
 * \param x a given node
 * \return the index of the root of the tree
 */
static int64_t comp_find( comp_tree_t *cc,  int64_t x)
{
  int64_t p;

  p = cc[x].parent;
  if( p == x )
    return(x);

  cc[x].parent = comp_find( cc, p );
  return( cc[x].parent );
};

/*! \brief Merge or take the union of two component trees, based on the rank of the components.
 * \param cc the component tree
 * \param r the root of one tree
 * \param s the root of the other
 * \param verbose print flag
 * Discussion about why this is good.
 */
static void comp_rank_union( comp_tree_t *cc, int64_t r, int64_t s, int verbose)
{
   if( r == s )
     return;
   if( cc[r].rank < cc[s].rank ) {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64"\n", r, s); 
     cc[r].parent = s;
   } else if( cc[s].rank < cc[r].rank ) {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64"\n", s, r); 
     cc[s].parent = r;
   } else {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64" rank++\n", r, s); 
     cc[r].parent = s;
     cc[s].rank  += 1;
   }
}

/*! \brief Merge or take the union of two component trees, based on the rank of the components.
 * \param cc the component tree
 * \param r the root of one tree
 * \param s the root of the other
 * \param e the node that caused us to realize the components were connected,
 *          essentially a random node in one of the trees
 * \param verbose print flag
 * Discussion about why this is way sub-optimal.
 */
static void comp_bad_union( comp_tree_t *cc, int64_t r, int64_t s, int64_t e, int verbose)
{
   if( r == s )
     return;
   if(verbose) printf("set root %"PRId64" to limb  %"PRId64"\n", r, e); 
   cc[r].parent = e;
}

/*!
 * \brief prints the comp_tree
 * \param *prefix a string one can use to keep the output straight
 * \param *cc the comp_tree_t data structure 
 * \param numverts the number of vertices 
 * \param verbose print flag
 */
void dump_comp_tree( char *prefix, comp_tree_t *cc, int64_t numverts, int verbose)
{
  int64_t i;

  if(verbose) printf("tree: %s ",prefix);
  for(i=0; i<numverts; i++) {
    if(verbose) printf(" %2"PRId64"",cc[i].parent);
  }
  if(verbose) printf("\n");
}

/*!
 * \brief This routine implements finds the connected components 
 * \param *numcomps address to return the number of components
 * \param *cc the comp_tree_t data structure 
 * \param *graph the matrix that specifies the graph
 * \return run time
 */
double concomp(int64_t *numcomps, comp_tree_t * cc, sparsemat_t *graph, int verbose, int which_union)
{
  int64_t i, j, k, r, s;
  
  double t1 = wall_seconds();
   
  dump_comp_tree(" ",cc, graph->numrows, verbose);
  for(i = 0; i < graph->numrows; i++){ 
    for(k=graph->offset[i]; k < graph->offset[i+1]; k++) {
      j = graph->nonzero[k];
  
      if( verbose > 0 ) printf("looking at (%"PRId64",%"PRId64"): ", i, j);
  
      r = comp_find(cc,i);
      s = comp_find(cc,j);
  
      if( verbose > 0 ) printf("parentof %"PRId64" is %"PRId64" , parentof %"PRId64" is %"PRId64"\n", i,r,j,s);
  
      if( verbose > 1 ) dump_comp_tree(">", cc, graph->numrows, verbose);
  
      if( which_union == 0 ) {
        comp_rank_union(cc, r, s, verbose);
      } else {
        comp_bad_union(cc, r, s, j, verbose);
      }
      if( verbose > 1 ) dump_comp_tree("<",cc, graph->numrows, verbose);
      if( verbose > 0 ) printf("parentof %"PRId64" is %"PRId64" , parentof %"PRId64" is %"PRId64"\n", r,comp_find(cc,r),s,comp_find(cc,s));
    }
  }
  for(i = 0; i < graph->numrows; i++){ 
    cc[i].parent = comp_find(cc,i);
  }
  dump_comp_tree("-",cc, graph->numrows, verbose);
  
  int64_t *comp_size = (int64_t *) calloc(graph->numrows, sizeof(int64_t));
  for(i = 0; i < graph->numrows; i++){ 
    comp_size[cc[i].parent] += 1;
  }
  
  for(i = 0; i < graph->numrows; i++){ 
    if( comp_size[i] > 0 ){
      *numcomps += 1;
      if(verbose) printf("component %"PRId64" has size %"PRId64"\n", i, comp_size[i]);
    }
  }
  
  t1 = wall_seconds() - t1;
  return(t1);
}


/********************************  argp setup  ************************************/
typedef struct args_t{
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
{
  {0}
};

static struct argp_child children_parsers[] =
{
  {&std_options_argp, 0, "Standard Options", -2},
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
  {0}
};

int main(int argc, char * argv[])
{  
  args_t args = {{0}};
  enum FLAVOR {GENERIC=1, ALL_Models=2};
  args.std.models_mask = ALL_Models-1;
  args.gstd.numrows = UNIONFIND_NUM_ROWS;
  struct argp argp = {options, parse_opt, 0, "Union Find for connected components", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_graph_options(&args.std, &args.gstd);
  write_std_options(&args.std);
  
  //override command line (these will lead to matrices with not quite the right number of nonzeros
  // if the user also used the -z flag.
  if ( (args.gstd.loops == 1) || (args.gstd.directed == 1) ) {
    fprintf(stderr,"WARNING: triangles counting requires undirected no-loops graph\n");
    args.gstd.loops = 0;
    args.gstd.directed = 0;
  }


  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: triangle: mat is NULL!\n");return(-1);}

  comp_tree_t * cc = calloc(mat->numrows, sizeof(comp_tree_t));
  int64_t i;
  for( i = 0; i< mat->numrows; i++){
    cc[i].parent = i;
    cc[i].rank   = 1;
  }
  
  if(args.std.dump_files) write_matrix_mm(mat, "unionfind_inmat");
  int verbose=0;

  int64_t num_components = 0;
  double laptime;
  uint32_t use_model;  
  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & args.std.models_mask ) {
    case GENERIC:
      printf("generic unionfind: ");
      laptime = concomp(&num_components, cc, mat, verbose, 1);
      break;
    default:
      continue;
    }
    printf(" number of components: %8"PRId64" in  %8.3lf seconds \n", num_components, laptime);
  }
  
  clear_matrix(mat);
  return(0);
}

