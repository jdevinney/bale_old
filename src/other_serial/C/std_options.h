#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H
#include <argp.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct std_args_t{
  int dump_files;
  int models_mask;
  int quiet;
  int64_t seed;
}std_args_t;

extern struct argp std_options_argp;
#endif
