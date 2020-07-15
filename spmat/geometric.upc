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

/*! \file geometric.upc
 * \brief Functions to support generating random geometric graphs.
 */
#include <spmat.h>
#include <exstack.h>

typedef struct point_t{
  double x;
  double y;
  int64_t index;
}point_t;

// returns the square of the L2 distance
double dist(point_t a, point_t b){
  return((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

// Input: a global point index and the total number of points
// Output: A local offset and PE number when the PEs are BLOCK distributed.
int64_t get_point_pe_and_offset(int64_t pindex, int64_t n, int64_t * pe){
  int64_t idx;
  int64_t upper_points_per_pe = n / THREADS + ((n % THREADS > 0) ? 1 : 0);
  int64_t rem = n % THREADS;
  if( pindex / upper_points_per_pe < rem ){
    *pe = pindex / upper_points_per_pe;
    idx = pindex % upper_points_per_pe;
  }else{
    *pe = (pindex - rem) / (upper_points_per_pe - 1);
    idx= (pindex - rem) % (upper_points_per_pe - 1);
  }
  return(idx);
}

// This takes edges with global indices and remaps the points to PEs in a BLOCK fashion.
// The edges are relabled to populate the distributed sparse adjacency matrix.
// pe: is output: the destination pe
// n: the number of points total
// inedge: the edge to distribute.
edge_t distribute_edge(int64_t * pe, int64_t n, edge_t * inedge){
  int64_t rowpe, colpe;
  int64_t row_off, col_off;
  edge_t e;

  row_off = get_point_pe_and_offset(inedge->row, n, &rowpe);
  col_off = get_point_pe_and_offset(inedge->col, n, &colpe);

  // maybe write this in SHMEM style
  e.row = row_off*THREADS + rowpe;
  e.col = col_off*THREADS + colpe;
  *pe = rowpe;
  return(e);
}

edge_t append_edges_between_sectors(int64_t this_sec_idx,
                                    int64_t other_sec_idx,
                                    SHARED point_t * points,
                                    SHARED int64_t * sector_offsets,
                                    double r2,
                                    int directed,
                                    edge_list_t * el){
  edge_t e;
  int64_t j, l;
  point_t tmp;
  int64_t * lsec_offsets = lgp_local_part(int64_t, sector_offsets);
  point_t * lpoints = lgp_local_part(point_t, points);
  int64_t other_sec_pe = other_sec_idx % THREADS;
  int64_t first_point_this = lsec_offsets[this_sec_idx];
  int64_t npts_this = lsec_offsets[this_sec_idx + 1] - first_point_this;

  // set lpoints to the beginning of points for this sector.
  lpoints = &lpoints[first_point_this];
  
  if(this_sec_idx == other_sec_idx){
    // We are finding edges between points within the same sector.
    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = (first_point_this + j)*THREADS + MYTHREAD;
      for(l = 0; l < j; l++){
        int64_t point_indexB = (first_point_this + l)*THREADS + MYTHREAD;
        if(dist(lpoints[j], lpoints[l]) < r2){
          if(directed && ((double)rand()/RAND_MAX > 0.5)){
            append_edge(el, point_indexB, point_indexA);
          }else{
            append_edge(el, point_indexA, point_indexB);
          }
        }
      }
    }
  }else{
    
    int64_t first_point_other = lgp_get_int64(sector_offsets, other_sec_idx);
    int64_t first_point_other_next = lgp_get_int64(sector_offsets, other_sec_idx + THREADS);
    int64_t num_pts_other = first_point_other_next - first_point_other;
    point_t tmp;
    first_point_other = first_point_other*THREADS + other_sec_pe;

    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = (first_point_this + j)*THREADS + MYTHREAD;
      point_t this_point = lpoints[j];
      for(l = 0; l < num_pts_other; l++){
        int64_t point_indexB = (first_point_other + l)*THREADS + other_sec_pe;
        lgp_memget(&tmp, points, sizeof(point_t), first_point_other + l*THREADS);
        if(dist(this_point, tmp) < r2){
          if(directed && (rand()/RAND_MAX > 0.5)){
            append_edge(el, point_indexB, point_indexA);
          }else{
            append_edge(el, point_indexA, point_indexB);
          }
        }
      }
    }
  }
}

/*! \brief Generates the adjacency matrix for a random geometric graph. 
  * 
  * See https://en.wikipedia.org/wiki/Random_geometric_graph
  * Each vertex corresponds to a point randomly placed in the unit square. Two vertices
  * are adjancent if their corresponding points are within distance r of each other.
  * 
  * \param n The number of vertices
  * \param r The size of the neighborhoods that determine edges.
  * \param type See edge_type. If undirected, this routine returns a lower triangular matrix.
  * \param loops See self_loops. Are there self loops?
  * \param flattening_mode 0: every sector gets the same number of points (at most off by 1)
  *                        1: points are uniformly distributed over the unit square and then
  *                            the rows are minimally shuffled so that each PE gets equal 
  *                            number of points. This maintains more locality.
  *                        2: points are uniformly distributed over the unit square and then
  *                            the rows are completely shuffled so that each PE gets equal.
  *                            This destroys most of the locality inherent in geometric 
  *                            graphs.
  * \param seed A seed for the RNG. This should be a single across all PEs (it will be modified by each PE individually).
  * \return An adjacency matrix (or lower portion of in the undirected case).
  */
sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed){
  
  // We generate n points (p_i) uniformly at random over the unit square.
  // Each point corresponds to a vertex (v_i).
  // There is an edge between v_A and v_B if the distance between the p_A and p_B is less than r.
  // To do this in parallel, we break up the unit square into chunks that are rxr grid.  
  // Each processor is responsible for some number of chunks.
  // We calculate edges by comparing distances between each point and every other point in
  // its own sector and in neighboring sectors.

  int flattening_mode = 0;
  int64_t i, j, k, l;
  double r2 = r*r;
  int64_t nsectors_across = ceil(1.0/r);
  int64_t nsectors = nsectors_across*nsectors_across;
  T0_printf("GEOMETRIC with r = %lf number of sectors = %ld\n", r, nsectors);

  
  // The PEs will get sectors in a round-robin fashion
  // Step 1. We count how many points land in each sector.
  // If flattening mode = 0, this is just n/nsectors.
  // O.w. each PE generates n/NPES random sector destinations
  // to create a distributed histogram of counts in each sector.
  SHARED int64_t * counts = lgp_all_alloc(nsectors,sizeof(int64_t));
  int64_t * lcounts = lgp_local_part(int64_t, counts);
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  if(flattening_mode == 0){
    for(i = 0; i < lnsectors; i++)
      lcounts[i] = (ln + lnsectors - i - 1)/lnsectors;
  }else{
    srand(seed + MYTHREAD + 1);

    for(i = 0; i < lnsectors; i++)
      lcounts[i] = 0;
    lgp_barrier();    

    for(i = 0; i < ln; i++){
      int64_t row = rand() % nsectors_across;
      int64_t col = rand() % nsectors_across;
      assert(row*nsectors_across + col < n);
      lgp_atomic_add(counts, row*nsectors_across + col, 1L);
    }
  }
  
  lgp_barrier();
  
  // Step 2. Figure out how many points landed in your sectors
  SHARED int64_t * first_pt_this_sector = lgp_all_alloc(nsectors, sizeof(int64_t));
  SHARED int64_t * sector_offsets = lgp_all_alloc(nsectors + THREADS, sizeof(int64_t));
  int64_t * lsector_offsets = lgp_local_part(int64_t, sector_offsets);
  
  int64_t my_total_points = 0;
  lsector_offsets[0] = 0;
  for(i = 0; i < lnsectors; i++){
    my_total_points += lcounts[i];
    lsector_offsets[i+1] = lsector_offsets[i] + lcounts[i];
  }
  
  // allocate space for the points
  int64_t max_lnum_points = lgp_reduce_max_l(my_total_points);
  int64_t sum = lgp_reduce_add_l(my_total_points);
  assert(sum == n);
  
  SHARED point_t * points = lgp_all_alloc(max_lnum_points, sizeof(point_t));  
  point_t * lpoints = lgp_local_part(point_t, points);
  
  // Step 3. Each PE generates the points for their own sectors.
  int64_t pt = 0;
  for(i = 0; i < lnsectors; i++){
    for(j = 0; j < lcounts[i]; j++){
      assert(pt < my_total_points);
      lpoints[pt].x = ((double)rand()/RAND_MAX)*r; //sector width + sector_start
      lpoints[pt].y = ((double)rand()/RAND_MAX)*r; //
      lpoints[pt].index = pt + lsector_offsets[i];
      pt++;
    }
  }
  assert(pt == my_total_points);
  lgp_barrier();

  int weighted = (edge_type == DIRECTED_WEIGHTED || edge_type == UNDIRECTED_WEIGHTED);
  int directed = (edge_type == DIRECTED) || (edge_type == DIRECTED_WEIGHTED);
  
  // Determine Edges
  int64_t space = ceil(1.1*my_total_points*(n*M_PI*r*r)/2.0);
  edge_list_t * el = init_edge_list(n, n, space);
  if(el == NULL){
    printf("ERROR: geometric graph: el is NULL\n");
    return(NULL);
  }
  int64_t point_index = 0;
  int64_t row, col;
  double val;
  for(i = 0; i < lnsectors; i++){
    int64_t global_sector_index = (i*THREADS + MYTHREAD);
    int64_t sector_row = global_sector_index / nsectors;
    int64_t sector_col = global_sector_index % nsectors;    

    // add loops
    if(loops == LOOPS){
      for(i = 0; i < lcounts[i]; i++){
        int64_t point_index = (lsector_offsets[i] + j)*THREADS + MYTHREAD;
        append_edge(el, point_index, point_index);
      }
    }
    
    
    // add edges between points within this sector
    append_edges_between_sectors(global_sector_index, global_sector_index,
                                 points, sector_offsets, r2, directed, el);

    // add edges between this sector and its neighbor to the West
    if(sector_col > 0){
      int64_t nbr_sector_index = global_sector_index - 1;
      append_edges_between_sectors(global_sector_index, nbr_sector_index,
                                   points, sector_offsets, r2, directed, el);
    }

    // add edges between this sector and its neighbor to the NorthWest
    if(sector_row > 0 && sector_col > 0){
      int64_t nbr_sector_index = (sector_row - 1)*nsectors_across + (sector_col - 1);
      append_edges_between_sectors(global_sector_index, nbr_sector_index,
                                   points, sector_offsets, r2, directed, el);
    }
    
    // add edges between this sector and its neighbor to the North
    if(sector_row > 0){
      int64_t nbr_sector_index = (sector_row - 1)*nsectors_across + sector_col;
      append_edges_between_sectors(global_sector_index, nbr_sector_index,
                                   points, sector_offsets, r2, directed, el);
    }

    // add edges between this sector and its neighbor to the NorthEast
    if(sector_row > 0 && sector_col < (nsectors_across - 1)){
      int64_t nbr_sector_index = (sector_row - 1)*nsectors_across + sector_col + 1;
      append_edges_between_sectors(global_sector_index, nbr_sector_index,
                                   points, sector_offsets, r2, directed, el);
    }
    

  }
  
  T0_printf("After comparing inter-point distances: lnnz = %ld\n",el->num);

  // Redistribute points to PEs. We will dole out points block-wise to
  // flatten the number of points per PE.

  // Do a histogram like thing to get row counts
  SHARED int64_t * row_counts = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lrow_counts = lgp_local_part(int64_t, row_counts);
  for(i = 0; i < ln; i++){
    lrow_counts = 0;
  }

  lgp_barrier();
  
  for(i = 0; i < el->num; i++)
    lgp_atomic_add(row_counts, el->edges[i].row, 1L);

  lgp_barrier();

  // now get the row counts for the points you will receive
  int64_t lrow = 0;
  lrow_counts = calloc(ln, sizeof(int64_t));  
  for(i = 0; i < lnsectors; i++){
    for(j = 0; j < lcounts[i]; j++){
      int64_t point_index = (lsector_offsets[i] + j)*THREADS + MYTHREAD;
      lrow_counts[lrow++] = lgp_get_int64(row_counts, point_index);
    }
  }
  lgp_barrier();
  
  lgp_all_free(row_counts);

  sparsemat_t * A = init_matrix(n, n, lnnz, weighted);
  if(!A){
    T0_printf("ERROR: geometric_random_graph: init_matrix failed.\n");
    return(NULL);
  }
  
  A->loffset[0] = 0;
  for(i = 1; i <= ln; i++){
    A->loffset[i] = A->loffset[i - 1] + lrow_counts[i-1];
    lrow_counts[i-1] = 0;
  }

  // now distribute the edges to the proper PEs and populate the matrix struct
  exstack_t * ex = exstack_init(buf_cnt, sizeof(edge_t));
  if(ex == NULL){return(NULL);}
  while(exstack_proceed(ex, (i == el->num))){
    for(i = 0; i < el->num; i++){
      edge_t e = distribute_edge(&pe, n, &el->edges[i]);
      if(exstack_push(ex, &e, pe) == 0L)
        break;
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &edge, NULL)){
      row = edge.row / THREADS;
      pos = A->loffset[row] + lrow_counts[row]++;
      A->nonzero[pos] = edge.col;
      if(weighted) A->value = (double)rand()/RAND_MAX;
    }
  }
  lgp_barrier();

  exstack_clear(&ex);
  free(el);
  free(lrow_counts);
  
  sort_nonzeros(A);

  lgp_barrier();

  return(A);
}


