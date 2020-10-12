/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file opts_demo.c
\brief Play with options
*/

#include "spmat_utils.h"
#include "std_options.h"


/********************************  argp setup  ************************************/
// Starting with defines both set to 0.
#define GRAPH_OPTS 0
#define APP_OPTS 0
// GRAPH_OPTS==0 and APP_OPTS==0 give you just the standard options with default values
//  ./opts_demo
//  ./opts_demo --help
//  ./opts_demo -s 131
//  ./opts_demo -s 131 -M15
//  ./opts_demo -s 131 -M15 -j opts_demo.json
// Then find models_mask and uncomment the line to set your own default 
// Next try (this should fail until we add the GRAPH options)
// ./opts_demo -F
// define GRAPH_OPTS 1 
// ./opts_demo --help
// define APP_OPTS 1  to see how to add an option for an individual app.
// ./opts_demo
// ./opts_demo -W 2
// ./opts_demo -W 2 -Y"You bet ja"
//
// Set GRAPH_OPT to 0 and REUSE_OPT to 1
#define REUSE_OPT 1

typedef struct args_t {
#if APP_OPTS
  int64_t Wacky;  
  char Yes_str[128];
#endif
#if REUSE_OPT
  int64_t num_things;
#endif
  std_args_t std;
#if GRAPH_OPTS
  std_graph_args_t gstd;
#endif
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
#if APP_OPTS
  case 'W': args->Wacky = atoi(arg); break;
  case 'Y': strcpy(args->Yes_str, arg); break;
#endif
#if REUSE_OPT
  case 'N': args->num_things = atoi(arg); break;
#endif
  case ARGP_KEY_INIT:
#if APP_OPTS
    args->Wacky = 0;
    strcpy(args->Yes_str,"Yes sir.");
#endif
    state->child_inputs[0] = &args->std;
#if GRAPH_OPTS
    state->child_inputs[1] = &args->gstd;
#endif
    break;
  }
  return(0);
}

static struct argp_option options[] = {
#if APP_OPTS
  {"Wacky",        'W', "NUM", 0, "Wacky FLAG"},
  {"Yes_str",      'Y', "YES",  0, "Yes string"},
#endif
#if REUSE_OPT
  {"num_things",      'N', "NUM",  0, "numthings, not the graph numrows"},
#endif 
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
#if GRAPH_OPTS
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
#endif
  {0}
};


int main(int argc, char * argv[]) 
{
  args_t args = {0};
  args.std.models_mask = 63;  // default value for model_mask can be overridden here
#if GRAPH_OPTS
  args.gstd.numrows = 100;
#endif
  struct argp argp = {options, parse_opt, 0, "opts_demo", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_options(&args.std);
#if GRAPH_OPTS
  write_std_graph_options(&args.std, &args.gstd);
#endif 
  
  if(args.std.dump_files) printf("dumpfile is on\n");
  
  if(args.std.json) printf("All bale init/write rountines send everything to json file %s\n", args.std.json_output);

#if APP_OPTS
  if(args.Wacky)
    printf("%s! This is Wacky times %ld\n", args.Yes_str, args.Wacky);
#endif

#if REUSE_OPT
    printf("option -N now means num_things = %ld\n", args.num_things);
#endif

  double laptime = 3.141592;
  bale_app_write_time(&args.std, "playing with options", laptime);
  bale_app_finish(&args.std);
  return(0);
}
