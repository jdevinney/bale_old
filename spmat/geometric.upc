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
int64_t sector_max;


// returns the square of the L2 distance
double dist(point_t a, point_t b){

  double d = (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
  //printf("comparing (%lf,%lf) and (%lf, %lf)...distance = %lf\n",a.x,a.y,b.x,b.y, d);
  return(d);
}

// 
// Input: an item's global index and the total number of points and the layout.
// The layout is either BLOCK or CYLIC.
//
// Output: A local offset and PE number.
int64_t get_pe_and_offset(int64_t pindex, int64_t n, int64_t * pe, layout_t layout){

  if(layout == CYCLIC){
    *pe = pindex % THREADS;
    return(pindex/THREADS + pindex % THREADS);
  }
  int64_t idx;
  int64_t upper_points_per_pe = n / THREADS + ((n % THREADS > 0) ? 1 : 0);
  int64_t rem = n % THREADS;

  if( (rem == 0) || (pindex / upper_points_per_pe < rem) ){
    *pe = pindex / upper_points_per_pe;
    idx = pindex % upper_points_per_pe;
  }else{
    *pe = (pindex - rem) / (upper_points_per_pe - 1);
    idx= (pindex - rem) % (upper_points_per_pe - 1);
  }
  //printf("Input %ld %ld: rem = %ld  upper = %ld pe = %ld idx = %ld\n", pindex, n, rem, upper_points_per_pe, *pe, idx);
  return(idx);
}

int64_t get_global_index_from_pe_and_offset(int64_t pe, int64_t offset){
  return(offset*THREADS + pe);
}

// This takes edges with global indices and remaps the points to PEs in a BLOCK fashion.
// The edges are relabled to populate the distributed sparse adjacency matrix.
// pe: is output: the destination pe
// n: the number of points total
// inedge: the edge to distribute.

// Or could we just name points by sector and offset?
// I say this because global point index is going to be annoying to compute.

edge_t distribute_edge(int64_t * pe, int64_t n, edge_t * inedge){
  int64_t rowpe, colpe;
  int64_t row_off, col_off;
  edge_t e;


  assert(inedge->row < n);
  row_off = get_pe_and_offset(inedge->row, n, &rowpe, BLOCK);

  // maybe write this in SHMEM style
  e.row = get_global_index_from_pe_and_offset(rowpe, row_off);
  e.col = inedge->col;
  //printf("Distributing edge %ld %ld...to %ld\n", inedge->row, inedge->col, e.row);
  *pe = rowpe;
  //printf("to %ld %ld on pe %ld\n", e.row, e.col, rowpe);
  return(e);
}

void append_edges_between_sectors(uint64_t this_sec_idx,
                                    uint64_t other_sec_idx,
                                    int64_t my_first_sector,
                                    int64_t nsectors,
                                    SHARED point_t * points,
                                    SHARED int64_t * first_point_in_sector,
                                    SHARED int64_t * counts,
                                    double r2,
                                    int directed,
                                    edge_list_t * el){
  assert(this_sec_idx < nsectors);
  assert(other_sec_idx < nsectors);
  edge_t e;
  int64_t j, l;
  point_t tmp;
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  int64_t * lcounts = lgp_local_part(int64_t, counts);
  int64_t * lfirst_point_in_sector = lgp_local_part(int64_t, first_point_in_sector);
  point_t * lpoints = lgp_local_part(point_t, points);
  
  //int64_t other_sec_pe = other_sec_idx % THREADS;
  int64_t this_sec_local_idx = this_sec_idx - my_first_sector;
  point_t * lpoints_this = &lpoints[this_sec_local_idx*sector_max];
  int64_t first_point_this = lfirst_point_in_sector[this_sec_local_idx];
  int64_t npts_this = lcounts[this_sec_local_idx];

  if((other_sec_idx >= my_first_sector) &&
     (other_sec_idx < my_first_sector + lnsectors)){

    // We are finding edges between points within two local sectors.
    int64_t other_sec_local_idx = other_sec_idx - my_first_sector;
    point_t * lpoints_other = &lpoints[other_sec_local_idx*sector_max];
    int64_t first_point_other = lfirst_point_in_sector[other_sec_local_idx];
    int64_t npts_other = lcounts[other_sec_local_idx];
    int64_t stop_pt = npts_other;
    
    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = first_point_this + j;
      //assert(point_indexA < 10);
      // If we are looking within the same sector, only compare
      // points where l < j.
      if(other_sec_idx == this_sec_idx)
        stop_pt = j;
      for(l = 0; l < stop_pt; l++){
        int64_t point_indexB = first_point_other + l;
        // assert(point_indexB < 10);
        if(dist(lpoints_this[j], lpoints_other[l]) < r2){
          if(directed && ((double)rand()/RAND_MAX > 0.5)){
            append_edge(el, point_indexB, point_indexA);
          }else{
            append_edge(el, point_indexA, point_indexB);
          }
        }
      }
    }
  }else{
    /* the other sector is not on our PE */
    int64_t other_sec_pe;
    int64_t other_sec_local_idx = get_pe_and_offset(other_sec_idx, nsectors, &other_sec_pe, BLOCK);
    int64_t first_point_other = lgp_get_int64(first_point_in_sector, other_sec_local_idx*THREADS + other_sec_pe);
    int64_t num_pts_other = lgp_get_int64(counts, get_global_index_from_pe_and_offset(other_sec_pe,other_sec_local_idx));
    if(num_pts_other == 0)
      return;
    /* get all the points from the other sector in one transfer */
    point_t * other_sec_pts = calloc(num_pts_other, sizeof(point_t));
    int64_t sector_start = sector_max*other_sec_local_idx*THREADS + other_sec_pe;
    lgp_memget(other_sec_pts, points, num_pts_other*sizeof(point_t), sector_start);

    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = first_point_this + j;
      for(l = 0; l < num_pts_other; l++){
        int64_t point_indexB = first_point_other + l;
        if(dist(lpoints_this[j], other_sec_pts[l]) < r2){
          if(directed && (rand()/RAND_MAX > 0.5)){
            append_edge(el, point_indexB, point_indexA);
          }else{
            append_edge(el, point_indexA, point_indexB);
          }
        }
      }
    }
    free(other_sec_pts);
  }
}

// The purpose of this function is to produce a histogram with the number
// of points in each sector. The function returns the counts array.
// We do this by having each PE generate n/THREADS random sector indices,
// incrementing the count for each.
SHARED int64_t * calculate_npoints_per_sector(int64_t n, int64_t nsectors){
  int64_t i;
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  SHARED int64_t * counts = lgp_all_alloc(nsectors,sizeof(int64_t));
  int64_t * lcounts = lgp_local_part(int64_t, counts);
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  
  for(i = 0; i < lnsectors; i++) lcounts[i] = 0;

  lgp_barrier();    
  
  for(i = 0; i < ln; i++){
    lgp_atomic_add(counts, rand() % nsectors, 1L);
  }
    
  lgp_barrier();
  
  return(counts);
}

// This is the Hillis and Steele algorithm as presented on wikipedia.
SHARED int64_t * prefix_scan(SHARED int64_t * counts, int64_t n){
  int64_t i, j;
  SHARED int64_t * cumsum = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lcumsum = lgp_local_part(int64_t, cumsum);
  for(i = 0; i < (n + THREADS - MYTHREAD - 1)/THREADS; i++)
    lcumsum[i] = 0;

  lgp_barrier();
  
  int64_t tn = n;
  int64_t lgn = 0;
  while(tn > 0){
    tn = tn >> 1;
    lgn++;
  }
  for(i = 0; i < lgn; i++){
    for(j = 0; j < n; j++){
      int64_t twoj = 1L<<j;
      if((j >= twoj) && ((j % n) == MYTHREAD)){
        lcumsum[j] += lgp_get_int64(cumsum, j-twoj);
      }
    }
    lgp_barrier();
  }
  return(cumsum);
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

  int64_t i, j, k, l;
  double r2 = r*r;
  int64_t nsectors_across = floor(1.0/r);
  double sector_width = 1.0/nsectors_across;
  int64_t nsectors = nsectors_across*nsectors_across;
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  T0_printf("GEOMETRIC GRAPH: r = %lf number of sectors = %ld sector_width = %lf\n",
            r, nsectors, sector_width);
  T0_printf("                 edge_type = %ld loops = %d\n", edge_type, loops);

  
  srand(seed + MYTHREAD + 1);
  // TODO: permute matrix at end to get Zmorton order (which would improve locality)
  //       or round-robin point order (which would destroy locality)

  // Step 1. Figure out how many points land in every sector.
  SHARED int64_t * counts = calculate_npoints_per_sector(n, nsectors);
  int64_t * lcounts = lgp_local_part(int64_t, counts);
  
  // Step 2. Figure out how many points landed in your sectors.
  // Note: We are not attaching sectors to anywhere in particular in
  // the unit square yet.  
  int64_t my_total_points = 0;
  int64_t my_sector_max = 0;
  for(i = 0; i < lnsectors; i++){
    my_total_points += lcounts[i];
    my_sector_max = (my_sector_max < lcounts[i] ? lcounts[i] : my_sector_max);
  }
  int64_t sum = lgp_reduce_add_l(my_total_points);
  assert(sum == n);
  
  // Step 3. Allocate space for global array of points.
  // sector_max is a global variable ...yuck.
  sector_max = lgp_reduce_max_l(my_sector_max);
  int64_t max_lnsectors = (nsectors + THREADS - 1)/THREADS; 
  SHARED point_t * points = lgp_all_alloc(sector_max*max_lnsectors*THREADS, sizeof(point_t));
  point_t * lpoints = lgp_local_part(point_t, points);
  
  // Step 4. Each PE generates the points for their own sectors.
  // The points in each sector are just random points in the rxr square
  // of the sector. We still have not placed the sectors in the unit square.
  int64_t pt = 0;
  for(i = 0; i < lnsectors; i++){
    pt = i*sector_max;
    for(j = 0; j < lcounts[i]; j++){    
      lpoints[pt].x = ((double)rand()/RAND_MAX)*sector_width;
      lpoints[pt].y = ((double)rand()/RAND_MAX)*sector_width; 
      pt++;
    }
    //Sort the points in each sector lexiographically
    qsort(&lpoints[i*sector_max], lcounts[i], sizeof(point_t), point_comp);
  }

  lgp_barrier();

  // Step 5.
  // We assign sectors to PEs and order sectors here. Sectors are assigned to PEs
  // in BLOCK fashion. The sectors will be laid out in the unit square in columns
  // starting with the sector whose lower left index is (0,0), proceeding up the column,
  // and then moving over to the bottom of the next column and repeating.

  // Here we compute the global index of the first point in each sector.
  int64_t my_first_sector = 0;
  for(i = 0; i < MYTHREAD; i++)
    my_first_sector += (nsectors + THREADS - i - 1)/THREADS;

  int64_t points_before_my_block = lgp_prior_add_l(my_total_points);
  SHARED int64_t * first_point_in_sector = lgp_all_alloc(max_lnsectors*THREADS, sizeof(int64_t));
  int64_t * lfirst_point_in_sector = lgp_local_part(int64_t, first_point_in_sector);

  lfirst_point_in_sector[0] = points_before_my_block;
  for(i = 1; i < lnsectors; i++){
    lfirst_point_in_sector[i] = lfirst_point_in_sector[i-1] + lcounts[i-1];
  }
  lgp_barrier();
#if 0
  for(i = 0; i < nsectors; i++)
    T0_printf("%ld ", lgp_get_int64(counts, i));
  T0_printf("\n");
  for(i = 0; i < nsectors; i++)
    T0_printf("%ld ", lgp_get_int64(first_point_in_sector, i));
  T0_printf("\n");
#endif
  
  // Step 5b. Correct points to have global x and y coordinates.
  for(i = 0; i < lnsectors; i++){
    int64_t sector = my_first_sector + i;
    double x_off = (sector / nsectors_across) * sector_width;
    double y_off = (sector % nsectors_across) * sector_width;
    pt = i*sector_max;
    //printf("PE %d: sector %ld (%lf %lf) (width = %lf)\n", MYTHREAD, i+my_first_sector,
    //x_off, y_off, sector_width);
    for(j = 0; j < lcounts[i]; j++){    
      lpoints[pt].x += x_off;
      lpoints[pt].y += y_off;      
      //printf("%lf %lf\n", lpoints[pt].x, lpoints[pt].y);
      pt++;
    }
  }

  int weighted = (edge_type == DIRECTED_WEIGHTED || edge_type == UNDIRECTED_WEIGHTED);
  int directed = (edge_type == DIRECTED) || (edge_type == DIRECTED_WEIGHTED);
  //T0_printf("Calculating edges...\n");fflush(0);
  
  // Step 6. Determine Edges  
  int64_t space = ceil(1.1*my_total_points*(n*M_PI*r*r)/2);
  //printf("PE %d: Initial allocation: %ld\n", MYTHREAD, space);
  edge_list_t * el = init_edge_list(space);
  if(el == NULL){
    printf("ERROR: geometric graph: el is NULL\n");
    return(NULL);
  }

  /* For each of our local sectors we compute all
   * edges 1) between points in that sector and 
   * 2) between points in that sector and sectors to the NW, W, SW, and S.
   * In a lower-triangular matrix, edges always go from p1 -> p2 where p2 < p1.
   * 
   * It is also mostly true that p2.x < p1.x or 
   * p2.x==p1.x and p2.y < p1.y. 
   *
   * To make this always true, we would need to compare to both N and S sectors and 
   * only do some of the comparisons. So maybe we don't need to sort
   * points within a sector?
   */  
  int64_t point_index = 0;
  int64_t sector, row, col;
  double val;
  for(sector = my_first_sector; sector < my_first_sector + lnsectors; sector++){
    int64_t sector_row = sector % nsectors_across;
    int64_t sector_col = sector / nsectors_across;

    // add loops
    if(loops == LOOPS){
      int64_t fpts = lfirst_point_in_sector[sector - my_first_sector];
      for(i = 0; i < lcounts[sector - my_first_sector]; i++){
        int64_t point_index = fpts + i;
        assert(point_index < n);
        append_edge(el, point_index, point_index);
      }
    }
    
    // add edges between points within this sector
    append_edges_between_sectors(sector, sector,
                                 my_first_sector, nsectors,
                                 points, first_point_in_sector, counts, r2, directed, el);
    
    // add edges between this sector and its neighbor to the NorthWest
    if(sector_row < (nsectors_across-1) && sector_col > 0){
      int64_t nbr_sector_index = (sector_col - 1)*nsectors_across + (sector_row + 1);
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, directed, el);
    }

    // add edges between this sector and its neighbor to the West
    if(sector_col > 0){
      int64_t nbr_sector_index = sector - nsectors_across;
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, directed, el);
    }
    
    // add edges between this sector and its neighbor to the SouthWest
    if(sector_row > 0 && sector_col > 0){
      int64_t nbr_sector_index = (sector_col - 1)*nsectors_across + (sector_row - 1);
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, directed, el);
    }

    // add edges between this sector and its neighbor to the South
    if(sector_row > 0){
      int64_t nbr_sector_index = sector - 1;
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, directed, el);
    }
    
  }

  lgp_barrier();
  lgp_all_free(points);
  
  //printf("PE %d: After comparing inter-point distances: lnnz = %ld\n",MYTHREAD, el->num);
  int64_t nnz = lgp_reduce_add_l(el->num);
  //T0_printf("Matrix has %ld nonzeros total\n", nnz);
  
  // We need to redistribute points to PEs to
  // flatten the number of points per PE. We do this by assigning points
  // to PEs in a BLOCK fashion. Note this is not a huge change from before
  // where we assigned sectors to PEs in block fashion. But sectors
  // do not all have the same number of points, so we need to smooth things out.
  
  // Step 7. Do a histogram to get vertex degrees (aka row_counts)
  SHARED int64_t * row_counts = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lrow_counts = lgp_local_part(int64_t, row_counts);
  for(i = 0; i < ln; i++) lrow_counts[i] = 0;
  
  lgp_barrier();
  
  int64_t pe;
  for(i = 0; i < el->num; i++){
    //if(el->edges[i].row >=n)
    //printf("illegal edge: %ld: %ld %ld\n", i, el->edges[i].row, el->edges[i].col);
    assert(el->edges[i].row < n);
    int64_t offset = get_pe_and_offset(el->edges[i].row, n, &pe, BLOCK);    
    assert(offset*THREADS + pe < n);
    lgp_atomic_add(row_counts, get_global_index_from_pe_and_offset(pe, offset), 1L);
    //lgp_atomic_add(row_counts, offset*THREADS + pe, 1L);
  }

  lgp_barrier();
  
  
  // Step 9. Sum the vertex degrees for the points you will receive
  int64_t lnnz = 0;
  for(i = 0; i < ln; i++) lnnz += lrow_counts[i];

  lgp_barrier();
  lgp_all_free(first_point_in_sector);
  
  // Step 10. Initialilze the sparse matrix
  T0_printf("Trying to init matrix with %ld nonzeros on PE 0\n", lnnz);fflush(0);
  sparsemat_t * A = init_matrix(n, n, lnnz, weighted);
  if(!A){
    T0_printf("ERROR: geometric_random_graph: init_matrix failed.\n");
    return(NULL);
  }

  // initialize the offsets for the sparse matrix
  A->loffset[0] = 0;
  for(i = 1; i <= ln; i++){
    A->loffset[i] = A->loffset[i - 1] + lrow_counts[i-1];
    lrow_counts[i-1] = 0;
  }
  assert(A->loffset[ln] == lnnz);

  // Step 11. Distribute the edges to the proper PEs
  // and populate the matrix struct
  exstack_t * ex = exstack_init(128, sizeof(edge_t));
  edge_t edge;
  if(ex == NULL){return(NULL);}
  i = 0;
  while(exstack_proceed(ex, (i == el->num))){
    for(; i < el->num; i++){
      edge_t e = distribute_edge(&pe, n, &el->edges[i]);
      //printf("Pushing %ld %ld to pe %ld\n", e.row, e.col, pe);
      if(exstack_push(ex, &e, pe) == 0L)
        break;
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &edge, NULL)){
      //row = edge.row / THREADS;
      row = edge.row/THREADS;
      assert(row < ln);
      int64_t pos = A->loffset[row] + lrow_counts[row]++;
      assert(pos < A->lnnz);
      assert(edge.col < n);
      A->lnonzero[pos] = edge.col;
      if(weighted) A->lvalue[pos] = (double)rand()/RAND_MAX;
    }
  }
  lgp_barrier();

  exstack_clear(ex);
  free(el);
  lgp_all_free(row_counts);

  sort_nonzeros(A);
  
  // print the matrix
#if 0
  T0_printf("Printing matrix in geometric...");
  for(i = 0; i < A->lnumrows; i++){
    T0_printf("row %ld: ",i);
    for(j = A->loffset[i]; j < A->loffset[i+1]; j++){
      T0_printf("%ld ", A->lnonzero[j]);
    }
    T0_printf("\n");
  }
#endif

  lgp_barrier();

  return(A);
}


