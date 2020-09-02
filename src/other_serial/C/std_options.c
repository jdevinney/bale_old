#include "std_options.h"

static int parse_opt(int key, char * arg, struct argp_state * state){
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'D': args->dump_files = 1; break;
  case 'M': args->models_mask = atol(arg); break;
  case 's': args->seed = atol(arg); break;
  case 'q': args->quiet = 1; break;
  case ARGP_KEY_INIT:
    args->quiet = 0;
    args->seed = 122222;
    args->models_mask = 1;
    args->dump_files = 0;
    break;
  }
  return(0);
}

static struct argp_option options[] =
  {
    {"dump_files", 'D', 0, 0, "Dump files for debugging"},
    {"models_mask", 'M', "MASK", 0, "Which flavors to run."},
    {"seed", 's', "SEED", 0, "Seed for RNG"},
    {"quiet", 'q', 0, 0, "No output during program execution"},
    {0}
  };

struct argp std_options_argp = {
  options, parse_opt, 0, 0, 0
};

