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

/*! \file unionfind.c
 * \brief Demo that uses a unionfind data structure to find connected components in a graph
 */

#include "spmat_utils.h"

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
    fprintf(stderr, "ERROR: generate_concomp_input: failed!\n");
    return(NULL);
  }

  // What if we permute it here
  spmat_stats(mat);

  if(dump_files) dump_matrix(mat, 20, "orig_uf.out");

  return( mat );
}

/*! \struct comp_tree_t 
 * \brief  guts of the whole thing
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
     if(verbose) printf("set parent of %ld to %ld\n", r, s); 
     cc[r].parent = s;
   } else if( cc[s].rank < cc[r].rank ) {
     if(verbose) printf("set parent of %ld to %ld\n", s, r); 
     cc[s].parent = r;
   } else {
     if(verbose) printf("set parent of %ld to %ld rank++\n", r, s); 
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
   if(verbose) printf("set root %ld to limb  %ld\n", r, e); 
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
    if(verbose) printf(" %2ld",cc[i].parent);
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
  
      if( verbose > 0 ) printf("looking at (%ld,%ld): ", i, j);
  
      r = comp_find(cc,i);
      s = comp_find(cc,j);
  
      if( verbose > 0 ) printf("parentof %ld is %ld , parentof %ld is %ld\n", i,r,j,s);
  
      if( verbose > 1 ) dump_comp_tree(">", cc, graph->numrows, verbose);
  
      if( which_union == 0 ) {
        comp_rank_union(cc, r, s, verbose);
      } else {
        comp_bad_union(cc, r, s, j, verbose);
      }
      if( verbose > 1 ) dump_comp_tree("<",cc, graph->numrows, verbose);
      if( verbose > 0 ) printf("parentof %ld is %ld , parentof %ld is %ld\n", r,comp_find(cc,r),s,comp_find(cc,s));
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
      if(verbose) printf("component %ld has size %ld\n", i, comp_size[i]);
    }
  }
  
  t1 = wall_seconds() - t1;
  return(t1);
}

int main(int argc, char * argv[])
{
  
  double er_prob = 0.01;
  uint32_t seed =  123456789;
  #define NUMROWS 20 
  int64_t numrows=NUMROWS;
  int64_t num_components = 0;
  int readgraph = 0;
  char filename[256]={"filename"};

  enum MODEL {GENERIC_Model=1, ALL_Models=2};
  uint32_t use_model;
  uint32_t models_mask = ALL_Models - 1;
  int printhelp=0;
  int dump_files = 0;
  int quiet=0;

  sparsemat_t *graph;
  double laptime;
  
  
  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:M:e:Df:q")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%ld" ,&numrows );  break;
    case 's': sscanf(optarg,"%d" ,&seed );  break;
    case 'M': sscanf(optarg,"%d" , &models_mask);  break;
    case 'e': sscanf(optarg,"%lg", &er_prob); break;
    case 'D': dump_files = 1; break;
    case 'f': readgraph = 1; sscanf(optarg, "%s", filename); break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }
  
  if( readgraph ) {
    graph = read_matrix_mm(filename);
    if(!graph){printf("ERROR: Read graph from %s Failed\n", filename); exit(1);}
  }  else {
    graph = generate_concomp_input(numrows, er_prob, seed, dump_files);
    if(!graph){printf("ERROR: graph is NULL!\n"); exit(1);}
  }
  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of union_find\n");
    fprintf(stderr,"number of nodes      (-n) %ld\n", graph->numrows);
    fprintf(stderr,"erdos_renyi_prob     (-e)= %lg\n", er_prob);
    fprintf(stderr,"readfile             (-f [%s])\n", filename); 
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"dump_files           (-D)\n");
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
    if(printhelp)
      return(0);
  } 
  comp_tree_t * cc = calloc(graph->numrows, sizeof(comp_tree_t));
  int64_t i;
  for( i = 0; i< graph->numrows; i++){
    cc[i].parent = i;
    cc[i].rank   = 1;
  }
  
  if(dump_files) 
    dump_matrix(graph, 20, "dump.out");

  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC_Model:
      if(!quiet) printf("generic unionfind: ");
      laptime = concomp(&num_components, cc, graph, quiet, 1);
      break;
    default:
      continue;
    }
    if(!quiet) printf(" number of components: %8ld in  %8.3lf seconds \n", num_components, laptime);
  }
  
  clear_matrix(graph);
  return(0);
}

