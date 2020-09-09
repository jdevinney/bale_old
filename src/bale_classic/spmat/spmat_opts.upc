#include "spmat_opts.h"

static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  switch(key){
  case 'd': args->directed = 1; break;
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; args->filename = arg; break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'l': args->loops = 1; break;
  case 'n': args->l_numrows = atol(arg); break;
  case 'w': args->weighted = 1; break;
  case 'z': args->nz_per_row = atof(arg); break;
  case ARGP_KEY_INIT:
    args->edge_prob = 0.0;
    args->readfile = 0;
    args->model = FLAT;
    args->l_numrows = 10000;
    args->nz_per_row = 10.0;
    args->directed = 0;
    args->weighted = 0;
    args->loops = 0;
    break;
  }
  return(0);
}

static struct argp_option graph_options[] =
  {
    {0, 0, 0, 0, "Input (as file):", 5},
    {"readfile",   'f', "FILE",  0, "Read input from a file"},
    {0, 0, 0, 0, "Input (as random graph):", 6},
    {"l_numrows",  'n', "NUM",   0, "Number of rows per PE in the matrix"},
    {"directed",   'd', 0,       0, "Specify a directed graph"},
    {"edge_prob",  'e', "EDGEP", 0, "Probability that an edge appears"},
    {"flat",       'F', 0,       0, "Specify flat random graph model"},
    {"geometric",  'G', 0,       0, "Specify geometric random graph model"},
    {"loops",      'l', 0,       0, "Specify you want to force loops into graph"},
    {"weighted",   'w', 0,       0, "Specify you want the edges to be weighted"},
    {"nz_per_row", 'z', "NZPR",  0, "Avg. number of nonzeros per row"},
    {0}
  };

struct argp std_graph_options_argp = {
  graph_options, graph_parse_opt, 0, 0, 0
};


sparsemat_t * get_input_graph(std_args_t * sargs, std_graph_args_t * gargs){

  sparsemat_t * mat;
  
  if(!gargs->readfile){

    int64_t numrows = gargs->l_numrows * THREADS;
    int64_t seed = MYTHREAD + sargs->seed;
    edge_type et;
    if(gargs->directed){
      if(gargs->weighted)
        et = DIRECTED_WEIGHTED;
      else
        et = DIRECTED;
    }else{
      if(gargs->weighted)
        et = UNDIRECTED_WEIGHTED;
      else
        et = UNDIRECTED;
    }
    self_loops loops = NOLOOPS;
    if(gargs->loops)      
      loops = LOOPS;

    resolve_edge_prob_and_nz_per_row(&gargs->edge_prob, &gargs->nz_per_row,
                                     numrows, et, loops);
    mat = random_graph(numrows, gargs->model, et, loops, gargs->edge_prob, seed + 2);
    
  }else{

    mat = read_matrix_mm_to_dist(gargs->filename);
  }
  if(!mat){T0_printf("ERROR: get_input_mat: mat is NULL!\n"); lgp_global_exit(1);}

  T0_printf("Input matrix has %"PRId64" rows and %"PRId64" nonzeros\n", mat->numrows, mat->nnz);

  return(mat);
}

void write_std_graph_options(std_graph_args_t * gargs){
  if(!gargs->readfile){
    T0_fprintf(stderr, "Generating a %s graph (-F or -G).\n",
               (gargs->model == FLAT ? "FLAT" : "GEOMETRIC"));
    T0_fprintf(stderr,"%s, %s, %s\n",
               (gargs->directed ? "Directed": "Undirected"),
               (gargs->weighted ? "Weighted": "Unweighted"),
               (gargs->loops ? "Loops": "No Loops"));
              
               
    T0_fprintf(stderr,"Number of rows per PE    (-n): %"PRId64"\n", gargs->l_numrows);
    T0_fprintf(stderr,"Avg # nnz per row        (-z): %2.2lf\n", gargs->nz_per_row);
    T0_fprintf(stderr,"Edge probability         (-e): %lf\n", gargs->edge_prob);

  }else{
    T0_fprintf(stderr,"Reading input from %s\n", gargs->filename);
  }
}
