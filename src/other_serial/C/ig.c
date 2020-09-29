/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file ig.c
 * \brief Program that runs C version of index gather
 *
 * Run ig --help or --usage for insructions on running.
 */
/*! \page ig_page Indexgather */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_sizes.h"


/*! \brief check that the indexgather worked. 
 * THIS REQUIRES that the source array, table, is set to just be minus the index
 * \param tgt the target array
 * \param index the indices that drive the gather
 * \param len length
 * \returns the number of disagreements.
 */

int64_t ig_check_and_zero(int64_t *tgt, int64_t *index, int64_t len)
{
  int64_t errors=0;
  int64_t i;
  for(i=0; i<len; i++){
    if( tgt[i] != -index[i] ) {
      errors ++;
      if( errors < 5 )  // print the first 5 errors and count the rest
        printf("  error tgt[%"PRId64"] = %"PRId64" != %"PRId64"\n", i, tgt[i], -index[i] );
    }
    tgt[i] = 0;
  }
  if( errors ) 
    printf(" total of %"PRId64" errors\n", errors);
  return(errors);
}


/*!
 * \brief This is the generic serial version of indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the source array of counts
 * \param num_req the length of the index array (number of updates)
 * \param *table the array from which the values are gathered
 * \return run time
 * This exercises a streaming load of index, then random loads from table       // TODO move to README
 * and a streaming store to tgt.
 */
double ig_generic(int64_t *tgt, int64_t *index, int64_t num_req,  int64_t *table) 
{
  int64_t i;
  double tm;

  tm = wall_seconds();

  for(i = 0; i < num_req; i++)
    tgt[i] = table[ index[i] ];

  tm = wall_seconds() - tm;
  return( tm );
}

/*!
 * \brief This routine implements a buffered version of indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the source array of counts
 * \param num_req the length of the index array (number of updates)
 * \param *table the array from which the values are gathered
 * \param log_tab_size the log of size of the table (number of bits in an index)
 * \return run time
 * The idea is to buffer up the indices based on their high bits.
 * Hopefully there will be a difference between doing full random loads and
 * doing loads that are close to each another. 
 */
double ig_buffered(int64_t *tgt, int64_t *index, int64_t num_req,  int64_t *table, int64_t table_size) 
{
  double tm;
  int64_t i, j;
  int64_t nbits;
  int64_t sort_shift;
#define LOG_NUM_BUFFERS 6  
#define NUM_BUFFERS (1L<<LOG_NUM_BUFFERS)
#define BUFFER_SIZE 128

  int64_t s, cnts[NUM_BUFFERS]; 
  int64_t table_idx[NUM_BUFFERS][BUFFER_SIZE];
  int64_t tgt_idx[NUM_BUFFERS][BUFFER_SIZE];

  assert(table_size > 0);
  nbits = 0;
  while(table_size>>nbits){   // shift table size to find the number of bit in it
    nbits++;
  }
  // We will put indices into buffer according to their top LOG_NUM_BUFFERS 
  sort_shift = nbits - LOG_NUM_BUFFERS;
  sort_shift = (sort_shift > 0) ? sort_shift : 0;

  for(i = 0; i < NUM_BUFFERS; i++)
    cnts[i] = 0L; 

  tm = wall_seconds();

  for(i = 0; i < num_req; i++){
    s = index[i] >> sort_shift;
    assert( (0 <= s) && (s<NUM_BUFFERS));
    assert( (0 <= cnts[s]) && (cnts[s]<BUFFER_SIZE));
    tgt_idx[s][cnts[s]] = i;
    table_idx[s][cnts[s]] = index[i];
    cnts[s]++;

    if( cnts[s] >= BUFFER_SIZE ) { 
      for(j = 0; j < cnts[s]; j++){
        tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
      }
      cnts[s] = 0;
    }
  }
  for(s = 0; s < NUM_BUFFERS; s++){
    for(j = 0; j < cnts[s]; j++){
      tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
    } 
  }

  tm = wall_seconds() - tm;
  return( tm );
}

/********************************  argp setup  ************************************/
typedef struct args_t{
  int64_t num_req;
  int64_t tbl_size;
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'n':
      args->num_req = atol(arg); break;
    case 'T':
      args->tbl_size = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"num_requests", 'n', "NUM",  0, "Number of requests to the table"},
    {"table_size",   'T', "SIZE", 0, "Number of entries in look-up table"},
    {0}
  };

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {0}
  };


int main(int argc, char * argv[])
{
  enum MODEL {GENERIC_Model=1, BUF_Model=2, ALL_Models=4};
  args_t args;
  args.tbl_size = IG_TABLE_SIZE;
  args.num_req = IG_NUM_UPDATES;
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Index gather from a table", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  int64_t *table, *tgt;
  int64_t *index;

  table   = calloc(args.tbl_size, sizeof(int64_t));
  tgt     = calloc(args.num_req, sizeof(int64_t));
  index   = calloc(args.num_req, sizeof(int64_t));

  printf("Index Gather Serial C\n");             // TODO delete
  printf("num_requests: %ld\n", args.num_req);
  printf("table size:  %ld\n", args.tbl_size);
  printf("----------------------\n");

  //populate table array and the index array
  // just fill the table with minus the index, so we check it easily
  int64_t i;
  for(i=0; i<args.tbl_size; i++)
    table[i] = -i;

  rand_seed(args.std.seed);
  for(i = 0; i < args.num_req; i++)
    index[i] = rand_int64(args.tbl_size);

  uint32_t use_model;
  int64_t errors = 0L;
  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & args.std.models_mask ){
    case GENERIC_Model:
      printf("Generic  IG: ");                                // TODO model_str
      laptime = ig_generic(tgt, index, args.num_req, table);
      errors += ig_check_and_zero(tgt, index, args.num_req);
      break;
    case BUF_Model:
      printf("Buffered IG: ");                                // TODO model_str
      laptime =ig_buffered(tgt, index, args.num_req,  table,  args.tbl_size); 
      errors += ig_check_and_zero(tgt, index, args.num_req);
      break;
    default:
      continue;
    }
    printf("  %8.3lf seconds\n", laptime);
  }

  free(table);
  free(tgt);
  free(index);

  if( errors ) {
    fprintf(stderr,"YOU FAILED!!!!\n"); 
  } 

  return(errors);
}


