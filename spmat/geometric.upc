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
int64_t get_pe_and_offset(int64_t pindex, int64_t n, int64_t * pe, layout layout){

  if(layout == CYCLIC){
    *pe = pindex % THREADS;
    return(pindex/THREADS);
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
                                  self_loops loops,
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
      for(l = 0; l < npts_other; l++){

        // skip self loops if not specified
        if((other_sec_idx == this_sec_idx)
           && (loops != LOOPS)
           && (l == j))
          continue;
        
        int64_t point_indexB = first_point_other + l;
        // assert(point_indexB < 10);
        if(dist(lpoints_this[j], lpoints_other[l]) < r2){
          //printf("PE %d: L (%ld) appending edge %ld %ld\n", MYTHREAD,this_sec_idx, point_indexA, point_indexB);
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
          //printf("PE %d: appending edge %ld %ld\n", MYTHREAD, point_indexA, point_indexB);
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


// This routine appends to the edge_list (el) all edges for points in a given sector.
// It does this by looking for point pairs that are closer than r, where one point lies
// in "sector" and the other point lies in any of the 8 neighboring sectors or "sector" itself.
// 
void append_edges_for_sector(int64_t sector,
                             int64_t my_first_sector,
                             int64_t nsectors,
                             SHARED point_t * points,
                             SHARED int64_t * first_point_in_sector,
                             SHARED int64_t * counts,
                             double r,
                             int directed,
                             self_loops loops,
                             edge_list_t * el){

  int64_t nsectors_across = floor(1.0/r);
  double r2 = r*r;
  int64_t sector_row = sector % nsectors_across;
  int64_t sector_col = sector / nsectors_across;
  int64_t i, j;

  // Append edges for between this sector and other sectors
  // We only need to look at all neighboring sectors
  // (including diagonals) AND the sector itself.
  for(i = -1; i < 2; i++){
    for(j = -1; j < 2; j++){
      int64_t nbrsector_col = sector_col + i;
      int64_t nbrsector_row = sector_row + j;
      if(nbrsector_col < 0 || nbrsector_col == nsectors_across)
        continue;
      if(nbrsector_row < 0 || nbrsector_row == nsectors_across)
        continue;
      int64_t nbr_sector_index = (nbrsector_col)*nsectors_across + (nbrsector_row);
      //      if(MYTHREAD == 1)
        //        printf("appending edges between %ld and %ld\n", sector, nbr_sector_index);
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, directed, loops, el);      
    }
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
  
  lgp_barrier();

  // Step 4.
  // We assign sectors to PEs and order sectors here. Sectors are assigned to PEs
  // in BLOCK (i.e not CYCLIC) fashion. The sectors will be laid out in the unit square in columns
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


  // Step 4. Each PE generates the points for their own sectors.
  int64_t pt = 0;
  for(i = 0; i < lnsectors; i++){
    int64_t sector = my_first_sector + i;
    double x_off = (sector / nsectors_across) * sector_width;
    double y_off = (sector % nsectors_across) * sector_width;
    pt = i*sector_max;
    for(j = 0; j < lcounts[i]; j++){    
      lpoints[pt].x = ((double)rand()/RAND_MAX)*sector_width + x_off;
      lpoints[pt].y = ((double)rand()/RAND_MAX)*sector_width + y_off; 
      pt++;
    }
    //Sort the points in each sector lexiographically
    //qsort(&lpoints[i*sector_max], lcounts[i], sizeof(point_t), point_comp);
  }


  int weighted = (edge_type == DIRECTED_WEIGHTED || edge_type == UNDIRECTED_WEIGHTED);
  int directed = (edge_type == DIRECTED) || (edge_type == DIRECTED_WEIGHTED);
  //T0_printf("Calculating edges...\n");fflush(0);
  
  // Step 6. Determine Edges
  int64_t space = ceil(1.1*my_total_points*(n*M_PI*r*r));
  printf("PE %d: Initial allocation: %ld\n", MYTHREAD, space);
  edge_list_t * el = init_edge_list(space);
  if(el == NULL){
    printf("ERROR: geometric graph: el is NULL\n");
    return(NULL);
  }

  /* For each of our local sectors we compute all possible edges */  
  for(i = 0; i < lnsectors; i++){
    append_edges_for_sector(my_first_sector + i, my_first_sector,
                            nsectors, points, first_point_in_sector,
                            counts, r, directed, loops, el);
    
  }

  lgp_barrier();
  lgp_all_free(points);
  
  //printf("PE %d: After comparing inter-point distances: lnnz = %ld\n",MYTHREAD, el->num);
  int64_t nnz = lgp_reduce_add_l(el->num);
  T0_printf("Symmetric matrix has %ld nonzeros total\n", nnz);
  

  // Step 7. We need to redistribute points to PEs to flatten
  // the number of points per PE. We do this by assigning points
  // to PEs in a BLOCK fashion. Note this is not a huge change from before
  // where we assigned sectors to PEs in block fashion. But sectors
  // do not all have the same number of points, so we need to smooth things out.
  //
  // At the same time, we are preparing to create the adjacency matrix for our graph.
  // Bale matrices are only stored in CYCLIC row layout, so we need to reinterpret the
  // labels on points. So the points that PE i owns will become points i, i + THREADS, i + 2*THREADS, etc.
  //
  // Once we relabel points, we can get rid of the unnecessary edges (those that represent nonzeros
  // above the diagonal).
  //

  int64_t pe;
  int64_t toss = 0, keep = 0;
  for(i = 0; i < el->num; i++){
    // Calling get_pe_and_offset gives us the PE and offset of a point if distribute the points
    // to PEs evenly and in BLOCK fashion.
    // We then convert that pe and offset into a new global index based on CYCLIC row layout.
    // Can we name these functions something better?
    int64_t roffset = get_pe_and_offset(el->edges[i].row, n, &pe, BLOCK);
    int64_t row_index = get_global_index_from_pe_and_offset(pe, roffset);
    int64_t coffset = get_pe_and_offset(el->edges[i].col, n, &pe, BLOCK);
    int64_t col_index = get_global_index_from_pe_and_offset(pe, coffset);
    //printf("PE %d: %ld,%ld goes to %ld,%ld\n", MYTHREAD,  el->edges[i].row, el->edges[i].col, row_index, col_index);
    if ( col_index > row_index) {
      el->edges[i].row = -1; // mark as NULL since this edge represents a 1 above the diagonal in the adjacency matrix.
      toss++;
    }else{
      keep++;
      el->edges[i].row = row_index;
      el->edges[i].col = col_index;
    }    
  }
  
  // Step 8. Do a histogram to get vertex degrees (aka row_counts)
  SHARED int64_t * row_counts = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lrow_counts = lgp_local_part(int64_t, row_counts);
  for(i = 0; i < ln; i++) lrow_counts[i] = 0;
  
  lgp_barrier();
  
  for(i = 0; i < el->num; i++){
    if(el->edges[i].row < 0) continue;
    lgp_atomic_add(row_counts, el->edges[i].row, 1L);
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
  // and populate the matrix struct.
  exstack_t * ex = exstack_init(128, sizeof(edge_t));
  edge_t edge;
  if(ex == NULL){return(NULL);}
  i = 0;
  while(exstack_proceed(ex, (i == el->num))){
    for(; i < el->num; i++){
      //edge_t e = distribute_edge(&pe, n, &el->edges[i]);
      edge_t e = el->edges[i];
      if(e.row < 0) continue;
      pe = e.row % THREADS;
      if(pe == MYTHREAD){
        int64_t row = e.row/THREADS;
        int64_t pos = A->loffset[row] + lrow_counts[row]++;
        A->lnonzero[pos] = e.col;
        if(weighted) A->lvalue[pos] = (double)rand()/RAND_MAX;
      }else{
        //printf("Pushing %ld %ld to pe %ld\n", e.row, e.col, pe);
        if(exstack_push(ex, &e, pe) == 0L)
          break;
      }
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &edge, NULL)){
      int64_t row = edge.row/THREADS;
      assert(row < ln);
      int64_t pos = A->loffset[row] + lrow_counts[row]++;
      assert(pos < A->lnnz);
      assert(edge.col < n);
      A->lnonzero[pos] = edge.col;
      if(weighted) A->lvalue[pos] = (double)rand()/RAND_MAX;
    }
  }

  lgp_barrier();
  
  for(i = 0; i < A->lnumrows; i++){
    assert(lrow_counts[i] == (A->loffset[i+1] - A->loffset[i]));
  }
  
  exstack_clear(ex);
  free(el);
  lgp_all_free(row_counts);
  
  sort_nonzeros(A);
  
  // print the matrix
#if 0
  T0_printf("Printing matrix in geometric...");
  for(i = 0; i < A->numrows; i++){
    T0_printf("row %ld: ",i);
    for(j = A->offset[i]; j < A->offset[i+THREADS]; j++){
      T0_printf("%ld ", A->nonzero[j*THREADS + i%THREADS]);
    }
    T0_printf("\n");
  }
#endif

  lgp_barrier();

  return(A);
}


