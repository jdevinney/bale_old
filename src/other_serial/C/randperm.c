/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file randperm.c
\brief Program that generates a random permutation.

Run randperm --help or --usage for insructions on running.
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"

/*! \brief A timing wrapper around the rand_perm routine that is in `spmat_utils.c`
\param len length of the permutation array
\param s seed for the random number generator.

The library uses the standard serial algorithm of Knuth-Yates.
*/
double randperm_generic(int64_t len, uint32_t s) 
{
  double tm;
  int64_t *p;
  tm = wall_seconds();
 
  p = rand_perm(len, s);
 
  tm = wall_seconds() - tm;
  if(!is_perm(p, len)){
    fprintf(stderr, "\nERROR: randperm_generic failed\n\n");
    exit(1);
  }
  free(p);
  return(tm);
}

/*!  \brief The "dart board" algorithm
\param len length of the permutation array
\param s seed for the random number generator.

We pick a dart board (an array) that is bigger than the permutation needed,
say twice as big. Then randomly throw darts at the  dart board, 
re-throwing any dart that hits a entry that is already occupied. 
Then we squeeze out the holes.
We picked the dartboard to be twice the size of the array 
so that even the last dart has a 50/50 of hitting an open entry.
*/
double randperm_dart(int64_t len, uint32_t s) 
{
  double tm;
  int64_t i, j, d;
  
  int64_t * dartboard = calloc(2*len, sizeof(int64_t));
  int64_t * p = calloc(len, sizeof(int64_t));
  tm = wall_seconds();
 
  // fill the dartboard with empty holes
  for(i=0; i<2*len; i++)
    dartboard[i] = -1;
 
  // randomly throw darts (entries 0 through len-1) 
  // at the dartboard until they all stick (land in an empty slot)
  rand_seed(s);
  for (i=0; i < len;  ) {
    d = rand_int64(2*len);
    if( dartboard[d] == -1 ){
      dartboard[d] = i;
      i++;
    }
  }
  // squeeze out the holes
  for(i=0, j=0; i<len; i++){
    while( dartboard[j] == -1 ) 
      j++;
    p[i] = dartboard[j];
    j++;
    if( j > 2*len ) {
      fprintf(stderr, "ERROR: randperm_dart ran out of darts\n");
      exit(2);
    }
  }
  tm = wall_seconds() - tm;
  free(dartboard);
 
  if(!is_perm(p, len)){
     fprintf(stderr, "ERROR: randperm_dart not a permutation\n");
     exit(1);
  }
  free(p);
  return(tm);
}

/*! 
\brief structure used by the randperm_sort routine

NB. Repeated key are bad, but tolerated. 
They would be ok if ties were broken randomly or if doubles were real numbers. 
*/
typedef struct idxkey_t {
  int64_t idx; //!< sequential index 0,...,len-1 
  double  key; //!< random key 
} idxkey_t;


/*! \brief the compare function for qsort called in the randperm_sort routine */
static int rp_comp( const void *a, const void *b) {
  idxkey_t * ak = (idxkey_t *)a;
  idxkey_t * bk = (idxkey_t *)b;
  return( (ak->key) - (bk->key) );
}

/*!  \brief The sorting algorithm to produce a random perm
\param len length of the permutation array
\param seed seed for the random number generator.

We form an array of (index, key) pairs. Then we randomly fill the keys
and sort on the keys. Then we read the permutation from the indices
*/
double randperm_sort(int64_t len, uint32_t seed) 
{
  double tm;
  int64_t i;
  
  
  idxkey_t *idxkey = calloc(len, sizeof(idxkey_t));
  int64_t *p = calloc(len, sizeof(int64_t));
  tm = wall_seconds();
  
  rand_seed(seed);
  for(i=0; i<len; i++) {
    idxkey[i].idx = i;
    idxkey[i].key = rand_double();
  }
  qsort( idxkey, len, sizeof(idxkey_t), rp_comp );
  
  for(i=0; i<len; i++)
    p[i] = idxkey[i].idx;
  
  tm = wall_seconds() - tm;
  if(!is_perm(p, len)){
     fprintf(stderr, "ERROR: sorting not a permutation\n");
     exit(1);
  }
  free(idxkey);
  free(p);
  return(tm);
}

/********************************  argp setup  ************************************/
typedef struct args_t{
  int64_t perm_size;
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'n':
      args->perm_size = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"perm_size", 'n', "NUM",  0, "Size of permutation to create."},
    {0}
  };

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {0}
  };


int main(int argc, char * argv[])
{
  enum MODEL {GENERIC_Model=1, DART_Model=2, SORT_Model=4, ALL_Models=8};
  args_t args={0};
  args.perm_size = RANDPERM_SIZE;
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Create a random permutation.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  printf("Randperm Serial C\n");                            //TODO 
  printf("permuatation size: %ld\n", args.perm_size);
  printf("----------------------\n");

  uint32_t use_model;
  double laptime = 0.0;
  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & args.std.models_mask ) {
    case GENERIC_Model:
      printf("Generic  perm: ");                 // TODO model_str
      laptime = randperm_generic(args.perm_size, args.std.seed);
      break;
    case DART_Model:
      printf("Dart     perm: ");                 // TODO model_str
      laptime = randperm_dart(args.perm_size, args.std.seed);
      break;
    case SORT_Model:
      printf("Sort     perm: ");                 // TODO model_str
      laptime = randperm_sort(args.perm_size, args.std.seed);
      break;
    default:
      continue;
    }
    printf("  %8.3lf seconds \n", laptime);
  }
	return(0);
}

