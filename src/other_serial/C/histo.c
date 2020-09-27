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

/*! \file histo.c\\
 * \brief Demo program that computes a histogram of uint64_t's
 *  The intent is to histogram a large number of words into
 *  a large number of bins.
 * 
 * Run histo --help or --usage for insructions on running.
*/
    
    #include "spmat_utils.h"
    #include "std_options.h"
    
/*! \page histo_page 
 * Statements about why histo is not as lame as it seems.
 */

/*! \brief This routine is a generic histogram of a large number of 
   items (int64_t's) into a large table
  - \param *index array of indices into the array of counts
  - \param num_ups the length of the index array (number of updates)
  - \param *counts pointer to the count array
  - \param v the amount to add to the counts array for each update.
           we can use v = (-1) to check against other implementations.
  - \return run time
 * This is an exercise is a streaming read of index and random writes into counts. 
 */

    double histo_generic(int64_t *index, int64_t num_ups,  int64_t *counts, int64_t v) 
    {
      double tm;
      int64_t i;
    
      tm = wall_seconds();
    
      for(i = 0; i < num_ups; i++)
        counts[index[i]] += v;
    
      tm = wall_seconds() - tm;
      return( tm );
    }
    
/*!
 \brief This routine does a buffered version of histogram.
  - \param *index array of indices into the array of counts
  - \param num_ups the length of the index array (number of updates)
  - \param *counts pointer to the count array
  - \param num_index_bits  FIXME
          we can use v = (-1) to check against other implementations.
  - \return run time

The exercise is a streaming read of indices and random writes into counts. 
 *
 */

    double histo_buffered(int64_t *index, int64_t num_ups,  int64_t *counts, int64_t num_index_bits) 
    {
      double tm;
      int64_t i,j;
      
      int64_t s, cnts[64]; 
      int64_t counts_idx[64][128];
      
      for(i = 0; i < 64; i++)
         cnts[i] = 0L; 
      
      tm = wall_seconds();
       
      for(i = 0; i < num_ups; i++){
        s = (index[i] >> (num_index_bits-6));     
        assert( (0 <= s) && (s<64));
        assert( (0 <= cnts[s]) && (cnts[s]<128));
        counts_idx[s][cnts[s]] = index[i];
        cnts[s]++;
        if( cnts[s] >= 128 ) { 
          for(j = 0; j < cnts[s]; j++)
            counts[ counts_idx[s][j] ]++;
          cnts[s] = 0;
        }
      }
      // empty any remaining elements in the buffers
      for(s = 0; s < 64; s++)
        for(j = 0; j < cnts[s]; j++)
          counts[ counts_idx[s][j] ]++;
      
      tm = wall_seconds() - tm;
      return( tm );
    }


/*! \brief comparison function to support histo_sorted
 */
    static int comp(const void *a, const void *b) 
    {
      return( *(uint64_t *)a - *(uint64_t *)b );
    }
    
/*!
 * \brief This routine does a sorting version of histogram.
 * Given the random index array, we first sort it so that the 
 * updates to the counts array occur in order.
 * NB: we don't charge for the sorting.
 * \param *index array of indices into the array of counts
 * \param num_ups the length of the index array (number of updates)
 * \param *counts pointer to the count array
 * \return run time
 */
    double histo_sorted(int64_t *index, int64_t num_ups,  int64_t *counts) 
    {
      double tm;
      int64_t i;
    
      qsort( index, num_ups, sizeof(int64_t),  comp );
    
      tm = wall_seconds();
      for(i = 0; i < num_ups; i++) 
        counts[index[i]] += 1;
      
      tm = wall_seconds() - tm;
      return( tm );
    }

/*!
 * \brief The next few structs and functions are for the command line parsing.
 */

    typedef struct args_t{
      int64_t num_ups;
      int64_t tbl_size;
      std_args_t std;
    }args_t;
    
    static int parse_opt(int key, char * arg, struct argp_state * state){
      args_t * args = (args_t *)state->input;
      switch(key)
        {
        case 'n':
          args->num_ups = atol(arg); break;
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
        {"num_updates",'n', "NUM", 0, "Number of updates to the histogram table"},
        {"table_size", 'T', "SIZE", 0, "Number of entries in the histogram table"},
        {0}
      };
    
    static struct argp_child children_parsers[] =
      {
        {&std_options_argp, 0, "Standard Options", -2},
        {0}
      };

/*!
 * \brief Driver program that calls the apps and collects the timing info.
 */    

    int main(int argc, char * argv[]) 
    {
      int64_t num_index_bits = 21;
    
      /* process command line */
      args_t args;
      args.tbl_size = 1L<<num_index_bits;
      args.num_ups = 100000;
      struct argp argp = {options, parse_opt, 0, "Accumulate updates into a table.", children_parsers};
      argp_parse(&argp, argc, argv, 0, 0, &args);
      
      enum MODEL {GENERIC_Model=1, BUF_Model=2, SORT_Model=4, ALL_Models=8};
      uint32_t use_model;
      uint32_t models_mask = args.std.models_mask;
      uint32_t quiet = args.std.quiet;
        
      num_index_bits = 0;
      while(args.tbl_size >> num_index_bits)
        num_index_bits++;
      
      // index is an array of indices into the counts array.
      int64_t * index  = calloc(args.num_ups, sizeof(int64_t)); assert(index != NULL);
      int64_t * counts = calloc(args.tbl_size, sizeof(int64_t)); assert(counts != NULL);  
    
      if(!quiet){
        printf("Histogram Serial C\n");
        printf("num_updates: %ld\n", args.num_ups);
        printf("table size:  %ld\n", args.tbl_size);
        printf("----------------------\n");
      }
      srand( args.std.seed );
      int64_t i;
      for(i = 0; i < args.num_ups; i++)
        index[i] = rand() % args.tbl_size; // index into the counts array
      
      double laptime = 0.0;
      int32_t num_runs = 0;
      for(use_model=1; use_model < ALL_Models; use_model *=2 ){
        switch( use_model & models_mask ) {
        case GENERIC_Model:
          if(!quiet) printf("Generic  histo: ");
          laptime = histo_generic(index, args.num_ups, counts, 1);
          num_runs++;
          break;
        case BUF_Model:
          if(!quiet) printf("Buffered histo: ");
          laptime = histo_buffered(index, args.num_ups, counts, num_index_bits);
          num_runs++;
          break;
        case SORT_Model:
          if(!quiet) printf("Sorted   histo: ");
          laptime = histo_sorted(index, args.num_ups, counts);
          num_runs++;
          break;
        default:
          continue;
        }
        if(!quiet) printf("  %8.3lf seconds \n", laptime);
      }
      
      // Check the buffered and sorted against the generic again
      int64_t errors = 0;
      laptime = histo_generic(index, args.num_ups,  counts, -num_runs);
      for(i = 0; i < args.tbl_size; i++) {
        if(counts[i] != 0L){
          errors++;
          if(errors < 5)  // print first five errors, report number of errors below
            fprintf(stderr,"ERROR: first five errors at %"PRId64" (= %"PRId64")\n", i, counts[i]);
        }
      }
      if(errors)
        fprintf(stderr,"FAILED!!!! total errors = %"PRId64"\n", errors);   
      
      free(index);
      free(counts);
      return(errors);
    }

