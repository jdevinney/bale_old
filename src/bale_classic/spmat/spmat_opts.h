#include <std_options.h>
#include "spmat.h"

typedef struct std_graph_args_t{
  int64_t l_numrows;
  int readfile;  
  char * filename;
  graph_model model;
  double edge_prob;
  int directed;
  int loops;
  int weighted;
  double nz_per_row;
} std_graph_args_t;

extern struct argp std_graph_options_argp;

sparsemat_t *       get_input_graph(std_args_t * sargs, std_graph_args_t * gargs);
void write_std_graph_options(std_graph_args_t * gargs);
