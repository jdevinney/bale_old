#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H
#include <argp.h>
#include <stdlib.h>
#include <stdint.h>
#include "spmat_utils.h"

typedef struct std_args_t{
  int dump_files;
  int models_mask;
  int quiet;
  int64_t seed;
}std_args_t;

typedef struct std_graph_args_t{
  int64_t numrows;
  int readfile;  
  char * filename;
  graph_model model;
  double edge_prob;
  double nz_per_row;
} std_graph_args_t;

extern struct argp std_options_argp;
extern struct argp std_graph_options_argp;
#endif
