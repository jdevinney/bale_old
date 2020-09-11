#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include "std_options.h"

static int std_parse_opt(int key, char * arg, struct argp_state * state){
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'b': args->buffer_size = atol(arg); break;
  case 'c': args->cores_per_node = atoi(arg); break;
  case 'D': args->dump_files = 1; break;
  case 'M': args->models_mask = atol(arg); break;
  case 's': args->seed = atol(arg); break;
  case 'q': args->quiet = 1; break;    
  case ARGP_KEY_INIT:
    args->buffer_size = 1024;
    args->cores_per_node = 0;
    args->quiet = 0;
    args->seed = 122222;
    args->models_mask = ALL_Models;
    args->dump_files = 0;
    break;
  }
  return(0);
}

static struct argp_option std_options[] =
  {
    {"buffer_size",   'b', "BUF", 0, "Aggregation buffer size"},
    {"cores_per_node",'c', "CPN", 0, "Specify cores per node for network injection rate statistics"},
    {"dump_files",    'D', 0,     0, "Dump files for debugging"},
    {"models_mask",   'M', "MASK",0, "Which flavors to run."},
    {"seed",          's', "SEED",0, "Seed for RNG"},
    {"quiet",         'q', 0,     0, "No output during program execution"},
    {0}
  };

struct argp std_options_argp = {
  std_options, std_parse_opt, 0, 0, 0
};

void write_std_options(std_args_t * sargs){
  fprintf(stderr,"Standard options:\n");
  fprintf(stderr,"----------------------------------------------------\n");
  fprintf(stderr,"buf_cnt (buffer size)    (-b): %"PRId64"\n", sargs->buffer_size);
  fprintf(stderr,"seed                     (-s): %"PRId64"\n", sargs->seed);
  fprintf(stderr,"cores_per_node           (-c): %d\n", sargs->cores_per_node);
  fprintf(stderr,"Models Mask              (-M): %d\n\n", sargs->models_mask);
}

#if 0
static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  switch(key){
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; args->filename = arg; break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'n': args->l_numrows = atol(arg); break;
  case 'z': args->nz_per_row = atof(arg); break;
  case ARGP_KEY_INIT:
    args->edge_prob = 0.0;
    args->readfile = 0;
    args->model = FLAT;
    args->l_numrows = 10000;
    args->nz_per_row = 10.0;
    break;
  }
  return(0);
}

static struct argp_option graph_options[] =
  {
    {0, 0, 0, 0, "Input (as file):", 5},
    {"readfile",   'f', "FILE",  0, "Read input from a file"},
    {0, 0, 0, 0, "Input (as random graph):", 6},
    {"l_numrows",    'n', "NUM",   0, "Number of rows per PE in the matrix"},
    {"edge_prob",  'e', "EDGEP", 0, "Probability that an edge appears"},
    {"flat",       'F', 0,       0, "Specify flat random graph model"},
    {"geometric",  'G', 0,       0, "Specify geometric random graph model"},
    {"nz_per_row", 'z', "NZPR",  0, "Avg. number of nonzeros per row"},
    {0}
  };

struct argp std_graph_options_argp = {
  graph_options, graph_parse_opt, 0, 0, 0
};
#endif



