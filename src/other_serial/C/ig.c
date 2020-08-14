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

/*! \file ig.c
 * \brief Demo program that runs index gather
 */

#include "spmat_utils.h"

/*! \page ig_page Indexgather */

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
 * \brief This routine implements generic serial version of indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the source array of counts
 * \param num_req the length of the index array (number of updates)
 * \param *table the array from which the values are gathered
 * \return run time
 * This exercises a streaming load of index, then random loads from table
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
 * \brief This routine implements a buffered  version of indexgather
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
double ig_buffered(int64_t *tgt, int64_t *index, int64_t num_req,  int64_t *table, int64_t log_tab_size) 
{
  int64_t i, j, s;
  double tm;

  int64_t cnts[64]; 
  int64_t table_idx[64][128];
  int64_t tgt_idx[64][128];

  for(i = 0; i < 64; i++)
    cnts[i] = 0L; 

  tm = wall_seconds();

  for(i = 0; i < num_req; i++){
    s = (index[i] >> (log_tab_size - 6));     
    assert( (0 <= s) && (s<64));
    assert( (0 <= cnts[s]) && (cnts[s]<128));
    tgt_idx[s][cnts[s]] = i;
    table_idx[s][cnts[s]] = index[i];
    cnts[s]++;

    if( cnts[s] >= 128 ) { 
      for(j = 0; j < cnts[s]; j++){
        tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
      }
      cnts[s] = 0;
    }
  }
  for(s = 0; s < 64; s++){
    for(j = 0; j < cnts[s]; j++){
      tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
    } 
  }

  tm = wall_seconds() - tm;
  return( tm );
}


int main(int argc, char * argv[])
{
  int64_t *table, *tgt;
  int64_t log_tab_size = 21;
  int64_t tbl_size = (1<<log_tab_size);
  int64_t num_req  = 100000;
  int64_t *index;
  uint32_t seed = 0x1234567;
  enum MODEL {GENERIC_Model=1, BUF_Model=2, ALL_Models=4};
  uint32_t use_model;
  uint32_t models_mask=ALL_Models - 1;
  int printhelp = 0;
  int quiet = 0;

  int opt; 
  while( (opt = getopt(argc, argv, "hT:n:s:M:q")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'T': sscanf(optarg,"%"SCNd64"" ,&log_tab_size);   break;
    case 'n': sscanf(optarg,"%"SCNd64"" ,&num_req);   break;
    case 's': sscanf(optarg,"%d" ,&seed);  break;
    case 'M': sscanf(optarg,"%d" ,&models_mask);  break;
    case 'q': quiet = 1; break;
    default:  break;
    }
  }

  if( printhelp || !quiet ) {
    fprintf(stderr,"Running C version of ig\n");
    fprintf(stderr,"help                  (-h)\n");
    fprintf(stderr,"Log of the Table size (-T)= %"PRId64"\n", log_tab_size);
    fprintf(stderr,"number of requests    (-n)= %"PRId64"\n", num_req );
    fprintf(stderr,"models_mask           (-M)= %d\n", models_mask);
    fprintf(stderr,"random seed           (-s)= %d\n", seed);
    fprintf(stderr,"quiet                 (-q)= %d\n", quiet);
    if(printhelp) 
      return(0);
  }

  table   = calloc(tbl_size, sizeof(int64_t));
  tgt     = calloc(num_req, sizeof(int64_t));
  index   = calloc(num_req, sizeof(int64_t));

  //populate table array and the index array
  int64_t i;
  for(i=0; i<tbl_size; i++)  // just fill the table with minus the index, so we check it easily
    table[i] = -i;

  srand(seed);
  for(i = 0; i < num_req; i++)
    index[i] = rand() % tbl_size; 

  int64_t errors = 0L;
  double laptime = 0.0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & models_mask ){
    case GENERIC_Model:
      if( !quiet ) printf("Generic  IG: ");
      laptime = ig_generic(tgt, index, num_req, table);
      errors += ig_check_and_zero(tgt, index, num_req);
      break;
    case BUF_Model:
      if( !quiet ) printf("Buffered IG: ");
      laptime =ig_buffered(tgt, index, num_req,  table,  log_tab_size); 
      errors += ig_check_and_zero(tgt, index, num_req);
      break;
    default:
      continue;
    }
    if( !quiet ) printf("  %8.3lf seconds\n", laptime);
  }

  free(table);
  free(tgt);
  free(index);

  if( errors ) {
    fprintf(stderr,"YOU FAILED!!!!\n"); 
  } 

  return(errors);
}


