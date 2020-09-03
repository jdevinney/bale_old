#include "std_options.h"

static int std_parse_opt(int key, char * arg, struct argp_state * state){
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'D': args->dump_files = 1; break;
  case 'M': args->models_mask = atol(arg); break;
  case 's': args->seed = atol(arg); break;
  case 'q': args->quiet = 1; break;
  case ARGP_KEY_INIT:
    args->quiet = 0;
    args->seed = 122222;
    args->models_mask = 15;
    args->dump_files = 0;
    break;
  }
  return(0);
}

static struct argp_option std_options[] =
  {
    {"dump_files", 'D', 0, 0, "Dump files for debugging"},
    {"models_mask", 'M', "MASK", 0, "Which flavors to run."},
    {"seed", 's', "SEED", 0, "Seed for RNG"},
    {"quiet", 'q', 0, 0, "No output during program execution"},
    {0}
  };

struct argp std_options_argp = {
  std_options, std_parse_opt, 0, 0, 0
};


static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  switch(key){
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; args->filename = arg; break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'n': args->numrows = atol(arg); break;
  case 'z': args->nz_per_row = atof(arg); break;
  case ARGP_KEY_INIT:
    args->edge_prob = 0.0;
    args->readfile = 0;
    args->model = FLAT;
    args->numrows = 1000;
    args->nz_per_row = 10.0;
    break;
  }
  return(0);
}

static struct argp_option graph_options[] =
  {
    {0, 0, 0, 0, "Input (as file):", 1},
    {"readfile",   'f', "FILE",  0, "Read input from a file"},
    {0, 0, 0, 0, "Input (as random graph):", 2},
    {"numrows",    'n', "NUM",   0, "Number of rows in a matrix"},
    {"edge_prob",  'e', "EDGEP", 0, "Probability that an edge appears"},
    {"flat",       'F', 0,       0, "Specify flat random graph model"},
    {"geometric",  'G', 0,       0, "Specify geometric random graph model"},
    {"nz_per_row", 'z', "NZPR",  0, "Avg. number of nonzeros per row"},
    {0}
  };

struct argp std_graph_options_argp = {
  graph_options, graph_parse_opt, 0, 0, 0
};
