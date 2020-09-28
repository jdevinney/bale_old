#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H

#include <stdlib.h>
#include <stdint.h>
#include "spmat_utils.h"
#include <argp.h>

#define C_BALE_VERSION 3.0

typedef struct std_args_t{
  int dump_files;
  int json;
  char json_output[128];
  int models_mask;
  int quiet;
  int64_t seed;
}std_args_t;

extern struct argp std_options_argp;
typedef struct std_graph_args_t{
  int64_t numrows;
  int readfile;  
  char filename[128];
  graph_model model;
  double edge_prob;
  int directed;
  int loops;
  int weighted;
  double nz_per_row;
  char kron_string[128];
  int kron_spec[64];
  int kron_num;
  int kron_mode;
} std_graph_args_t;

extern struct argp std_graph_options_argp;

/* defines to support the different models of global and buffered references */

sparsemat_t *  get_input_graph(std_args_t * sargs, std_graph_args_t * gargs);
void           write_std_options(std_args_t * sargs);
void           write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs);

int  bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs);
void bale_app_finish(std_args_t * sargs);
void bale_app_write_int(std_args_t * sargs, char * key, int64_t val);
void bale_app_write_time(std_args_t * sargs, char * model_str, double time);

#endif
