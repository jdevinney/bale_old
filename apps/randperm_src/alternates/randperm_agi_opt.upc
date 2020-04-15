#include <libgetput.h>
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
#include <exstack.h>

SHARED int64_t * rand_permp_agi_opt(int64_t N, int seed) {  
  int64_t * ltarget, *lperm;
  int64_t r, i, j, t, idx, pe;
  int64_t pos, numdarts, numtargets, lnumtargets;
  int64_t buf_cnt = 1024;
  double t1 = wall_seconds();
  if( seed != 0 ) srand48( seed );  
  
  //T0_printf("Entering rand_permp_atomic...");fflush(0);
  SHARED int64_t * inbox = lgp_all_alloc(THREADS*THREADS, sizeof(int64_t));
  if( inbox == NULL ) return(NULL);
  int64_t * linbox = lgp_local_part(int64_t, inbox);
  
  SHARED int64_t * perm = lgp_all_alloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  lperm = lgp_local_part(int64_t, perm);

  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = 2*N;
  int64_t l_M = (M + THREADS - MYTHREAD - 1)/THREADS;
  
  SHARED int64_t * target = lgp_all_alloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  ltarget = lgp_local_part(int64_t, target);
  
  for(i=0; i<l_M; i++)
    ltarget[i] = -1L;
  lgp_barrier();

  int64_t * PE_hist = calloc(THREADS, sizeof(int64_t));

  // figure out which PE each "dart" will land on and initialize perm
  for(i = 0; i < l_N; i++){
    r = lrand48() % THREADS;
    PE_hist[r]++;
    lperm[i] = r;
  }

  //tell the other PEs how many darts you will be sending them.
  for(i = 0; i < THREADS; i++)
    lgp_put_int64(inbox, MYTHREAD*THREADS + i, PE_hist[i]);
  
  lgp_barrier();

  int64_t total_local_darts = 0;
  for(i = 0; i < THREADS; i++)
    total_local_darts += linbox[i];
  
  int64_t * tmp_perm = calloc(total_local_darts, sizeof(int64_t));

  lgp_barrier();

  //T0_fprintf(stderr,"time phase 1: %lf\n", wall_seconds() - t1);
  //t1 = wall_seconds();

  //send the darts!
  exstack_t * ex = exstack_init( buf_cnt, sizeof(int64_t));
  if( ex == NULL ) return(NULL);

  pos = i = 0;
  while(exstack_proceed(ex, (i == l_N))){
    while(i < l_N){
      idx = i*THREADS + MYTHREAD;
      pe = lperm[i];
      if( !exstack_push(ex, &idx, pe) )
        break;
      i++;
    }

    exstack_exchange(ex);
    
    while(exstack_pop(ex, &idx, NULL)){
      tmp_perm[pos++] = idx;
    }
  }

  exstack_clear(ex);
  
  
  if(pos != total_local_darts){
    fprintf(stderr,"ERROR: randperm_agi_opt: pos != total_local_darts (%"PRId64" %"PRId64")\n", pos, total_local_darts);return(NULL);}
  pos = lgp_reduce_add_l(pos);
  if(pos != N){fprintf(stderr,"ERROR: randperm_agi_opt: total_darts != N (%"PRId64" %"PRId64")\n", pos, N);return(NULL);}

  //T0_fprintf(stderr,"time phase 2: %lf\n", wall_seconds() - t1);
  //t1 = wall_seconds();

  // This next section is the slowest section and could be improved
  
  // shuffle tmp_perm
  for(i = 0; i < total_local_darts; i++){
    j = i + lrand48() % (total_local_darts - i);
    t = tmp_perm[j];
    tmp_perm[j] = tmp_perm[i];
    tmp_perm[i] = t;
  }
    
  pos = lgp_prior_add_l(total_local_darts);    // my first index in the perm array is the number 
                                               // of elements contained on the smaller threads
  for(i = 0; i < total_local_darts; i++){
    lgp_put_int64(perm, pos, tmp_perm[i]);
    pos++;  
  }

  lgp_barrier();

  lgp_all_free(inbox);
  free(tmp_perm);

  lgp_barrier();
  //T0_fprintf(stderr,"time phase 3: %lf\n", wall_seconds() - t1);
  //t1 = wall_seconds();
  return(perm);
}
