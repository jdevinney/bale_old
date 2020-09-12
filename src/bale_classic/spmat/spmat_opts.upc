#include "spmat_opts.h"

static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  switch(key){
  case 'd': args->directed = 1; break;
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; args->filename = arg; break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'K':
    args->model = KRONECKER;
    args->kron_string = arg;
    char * ptr = arg;
    //T0_fprintf(stderr, "%s\n", arg);
    sscanf(ptr, "%d:", &args->kron_mode);
    //printf("mode %d\n", args->kron_mode);
    int i = 0;
    int tmp;
    ptr+=2;
    while(i < 63 && sscanf(ptr, "%d", &args->kron_spec[i])){
      //fprintf(stderr,"%d\n", args->kron_spec[i]);
      i++;
      ptr++; 
      if(*ptr != 'x'){
        //fprintf(stderr,"nope");
        break;
      }
      ptr++; 
    }

    args->kron_num = i;
    break;
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
    args->kron_string = NULL;
    break;
  case ARGP_KEY_END:
    if(args->directed)
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->l_numrows*THREADS,
                                       (args->weighted ? DIRECTED_WEIGHTED : DIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    else{
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->l_numrows*THREADS,
                                       (args->weighted ? UNDIRECTED_WEIGHTED : UNDIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    }
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
    {"kronecker",  'K', "KSTR",  0, "Specify a Kronecker product graph.\nKSTR must be a string of the form MODE:S1xS2x...Sk where MODE is 0, 1, or 2 and the Si are small integers that specify the stars whose product is the kronecker product graph. For instance -K 0:3x4x5 specifies MODE 0 and takes the product of K_{1,3}, K_{1,4}, and K_{1,5}. MODE 0 : No triangles. MODE 1: Many triangles. MODE 2: Few triangles. "},
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
    
    if(gargs->model == KRONECKER){

      // generate a random kronecker graph from a string (mode:#x#x...#)
      mat = generate_kronecker_graph_from_spec(gargs->kron_mode,
                                               gargs->kron_spec,
                                               gargs->kron_num);
      
    }else{
      
      // Generate a random FLAT or GEOMETRIC graph
      int64_t numrows = gargs->l_numrows * THREADS;
      int64_t seed = MYTHREAD + sargs->seed;
      edge_type et;
      if(gargs->directed){
        et = (gargs->weighted ? DIRECTED_WEIGHTED : DIRECTED);
      }else{
        et = (gargs->weighted ? UNDIRECTED_WEIGHTED: UNDIRECTED);
      }
      self_loops loops = (gargs->loops ? LOOPS : NOLOOPS);
      
      mat = random_graph(numrows, gargs->model, et, loops, gargs->edge_prob, seed + 2);
    }
  }else{
    
    // Read a matrix from a file
    mat = read_matrix_mm_to_dist(gargs->filename);
    
  }
  if(!mat){T0_printf("ERROR: get_input_mat: mat is NULL!\n"); lgp_global_exit(1);}

  if(sargs->json_output == NULL){
    T0_fprintf(stderr,"Input matrix:\n");
    T0_fprintf(stderr,"----------------------------------------------------\n");
    T0_fprintf(stderr,"\t%"PRId64" rows\n\t%"PRId64" columns\n\t%"PRId64" nonzeros\n\n",
               mat->numcols, mat->numrows, mat->nnz);
  }else{
    FILE * jp = fopen(sargs->json_output, 'a');    
    T0_fprintf(fp, "matrix_numrows: \"%"PRId64"\",\n", mat->numrows);
    T0_fprintf(fp, "matrix_numcols: \"%"PRId64"\",\n", mat->numcols);
    T0_fprintf(fp, "matrix_nnz: \"%"PRId64"\",\n", mat->nnz);
    fclose(jp);
  }

  return(mat);
}

void write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs){
  T0_fprintf(stderr,"Input Graph/Matrix parameters:\n");
  T0_fprintf(stderr,"----------------------------------------------------\n");
  if(!gargs->readfile){
    char model[32];
    if(gargs->model == FLAT)
      sprintf(model, "FLAT (-F)");
    else if(gargs->model == GEOMETRIC)
      sprintf(model, "GEOMETRIC (-G)");
    else if(gargs->model == KRONECKER)
      sprintf(model, "KRONECKER (-K)");
    
    T0_fprintf(stderr, "Graph model: %s.\n", model);
    T0_fprintf(stderr,"%s, %s, %s\n",
               (gargs->directed ? "Directed": "Undirected"),
               (gargs->weighted ? "Weighted": "Unweighted"),
               (gargs->loops ? "Loops": "No Loops"));
              
               
    T0_fprintf(stderr,"Number of rows per PE    (-n): %"PRId64"\n", gargs->l_numrows);
    T0_fprintf(stderr,"Avg # nnz per row        (-z): %2.2lf\n", gargs->nz_per_row);
    T0_fprintf(stderr,"Edge probability         (-e): %lf\n\n", gargs->edge_prob);

  }else{
    T0_fprintf(stderr,"Reading input from %s\n\n", gargs->filename);
  }
}
