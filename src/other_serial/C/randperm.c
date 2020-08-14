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

/*! \file randperm.c
 * \brief Demo program that generates a random permutation.
 * The well known serial algorithm is from Knuth (Vol FIXME).
 * The "dart throwing" and the "sorting algorithm" can be parallelized.
 */

#include "spmat_utils.h"

/*! \page randperm_page 
 * Generate a uniform random permutation.  We consider three different algorithm:
  - the standard array method known at least by the names
     Fisher-Yates or Knuth shuffle.
  - the dart board algorithm 
  - random sort algorithm
*/

/*! \brief A timing wrapper around the rand_perm routine 
 * that is in spmat_utils.
 * \param l length of the permutation array
 * \param s seed for the random number generator.
 */
double randperm_generic(int64_t l, uint32_t s) 
{
  double tm;
  int64_t *p;
  tm = wall_seconds();
 
  p = rand_perm(l, s);          // from the sparsemat library
 
  tm = wall_seconds() - tm;
  if(!is_perm(p, l)){
    fprintf(stderr, "\nERROR: randperm_generic failed\n\n");
    exit(1);
  }
  free(p);
  return(tm);
}

/*!
 * \brief The "dart board" algorithm
 * \param l length of the permutation array
 * \param s seed for the random number generator.
 * We pick a dart board (an array) that is bigger than the 
 * permutation needed. Then randomly throw darts at the 
 * dart board, re-throwing any dart that hits a entry 
 * that is already occupied. Then we squeeze out the holes.
 * We picked the dartboard to be twice the size of the array 
 * so that even the last dart has a 50/50 of hitting an open entry.
 */
double randperm_dart(int64_t l, uint32_t s) 
{
  double tm;
  int64_t i, j, d;
  tm = wall_seconds();
  
  int64_t * dartboard = calloc(2*l, sizeof(int64_t));
 
  // fill the dartboard with empty holes
  for(i=0; i<2*l; i++)
    dartboard[i] = -1;
 
  // randomly throw darts (entries 0 through l-1) 
  // at the dartboard until they all stick (land in 
  // an empty slot)
  srand(s);
  for (i=0; i < l;  ) {
    d = rand() % (2*l);
    if( dartboard[d] == -1 ){
      dartboard[d] = i;
      i++;
    }
  }
  // squeeze out the holes
  int64_t * p = calloc(l, sizeof(int64_t));
  for(i=0, j=0; i<l; i++){
    while( dartboard[j] == -1 ) 
      j++;
    p[i] = dartboard[j];
    j++;
    if( j > 2*l ) {
      fprintf(stderr, "ERROR: randperm_dart ran out of darts\n");
      exit(2);
    }
  }
  tm = wall_seconds() - tm;
  free(dartboard);
 
  if(!is_perm(p, l)){
     fprintf(stderr, "ERROR: randperm_dart not a permutation\n");
     exit(1);
  }
  free(p);
  return(tm);
}

/*! \struct idxkey_t 
 * \brief structure used by the randperm_sort routine
 * NB. Repeated key are bad, but tolerated. They would
 *  be ok if ties were broken randomly. 
 */
typedef struct idxkey_t {
  int64_t idx; //!< sequential index 0,...,l-1 
  int64_t key; //!< random key
} idxkey_t;


/*! \brief the compare function for qsort called in the 
 * randperm_sort routine
 */
static int rp_comp( const void *a, const void *b) {
  idxkey_t * ak = (idxkey_t *)a;
  idxkey_t * bk = (idxkey_t *)b;
  return( (ak->key) - (bk->key) );
}

/*!
 * \brief The sorting algorithm to produce a random perm
 * \param l length of the permutation array
 * \param s seed for the random number generator.
 * We form an (index, key) pair. Then randomly fill the key
 * and then sort on the key. Then read out the permuted index.
 */
double randperm_sort(int64_t l, uint32_t seed) 
{
  double tm;
  int64_t i;
  
  tm = wall_seconds();
  
  idxkey_t * idxkey = calloc(l, sizeof(idxkey_t));
  
  srand(seed);
  for(i=0; i<l; i++) {
    idxkey[i].idx = i;
    idxkey[i].key = rand();
  }
  qsort( idxkey, l, sizeof(idxkey_t), rp_comp );
  
  int64_t * p = calloc(l, sizeof(int64_t));
  for(i=0; i<l; i++)
    p[i] = idxkey[i].idx;
  
  tm = wall_seconds() - tm;
  if(!is_perm(p, l)){
     fprintf(stderr, "ERROR: sorting not a permutation\n");
     exit(1);
  }
  free(p);
  free(idxkey);
  return(tm);
}

int main(int argc, char * argv[])
{
  int64_t length = (1L<<23);
  uint32_t seed = 101892;
  double laptime = 0.0;

  enum MODEL {GENERIC_Model=1, DART_Model=2, SORT_Model=4, ALL_Models=8};
  uint32_t use_model;
  uint32_t models_mask = ALL_Models - 1;
  int printhelp=0;
  int quiet=0;

  int opt; 
  while( (opt = getopt(argc, argv, "hn:s:M:q")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'n': sscanf(optarg,"%"SCNd64 ,&length    );   break;
    case 's': sscanf(optarg,"%d", &seed); break;      
    case 'M': sscanf(optarg,"%d" , &models_mask);  break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of rand_perm\n");
    fprintf(stderr,"help                 (-h)\n");
    fprintf(stderr,"permutation size     (-n)= %"PRId64"\n", length);
    fprintf(stderr,"random seed          (-s)= %d\n", seed);
    fprintf(stderr,"models_mask          (-M)= %d\n", models_mask);
    fprintf(stderr,"quiet                (-q)= %d\n", quiet);
    if(printhelp)
      return(0);
  }

  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC_Model:
      if(!quiet) printf("Generic  perm: ");
      laptime = randperm_generic(length, seed);
      break;
    case DART_Model:
      if(!quiet) printf("Dart     perm: ");
      laptime = randperm_dart(length, seed);
      break;
    case SORT_Model:
      if(!quiet) printf("Sort     perm: ");
      laptime = randperm_sort(length, seed);
      break;
    default:
      continue;
    }
    if(!quiet) printf("  %8.3lf seconds \n", laptime);
  }
	return(0);
}

