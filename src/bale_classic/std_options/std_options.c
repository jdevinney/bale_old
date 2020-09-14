#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include <libgetput.h>
#include "std_options.h"

static int std_parse_opt(int key, char * arg, struct argp_state * state){
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'b': args->buffer_size = atol(arg); break;
  case 'c': args->cores_per_node = atoi(arg); break;
  case 'D': args->dump_files = 1; break;
  case 'j': args->json_output = arg; break;
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
    args->json_output = NULL;
    break;
  }
  return(0);
}

static struct argp_option std_options[] =
  {
    {"buffer_size",   'b', "BUF", 0, "Aggregation buffer size"},
    {"cores_per_node",'c', "CPN", 0, "Specify cores per node for network injection rate statistics"},
    {"dump_files",    'D', 0,     0, "Dump files for debugging"},
    {"json_output",   'j', "FILE",0, "Output results to a json file, rather than to stderr"},
    {"models_mask",   'M', "MASK",0, "Which flavors to run."},
    {"seed",          's', "SEED",0, "Seed for RNG"},
    {"quiet",         'q', 0,     0, "No output during program execution"},
    {0}
  };

struct argp std_options_argp = {
  std_options, std_parse_opt, 0, 0, 0
};

void write_std_options(std_args_t * sargs){
  if(sargs->json_output == NULL){
    fprintf(stderr,"Standard options:\n");
    fprintf(stderr,"----------------------------------------------------\n");
    fprintf(stderr,"buf_cnt (buffer size)    (-b): %"PRId64"\n", sargs->buffer_size);
    fprintf(stderr,"seed                     (-s): %"PRId64"\n", sargs->seed);
    fprintf(stderr,"cores_per_node           (-c): %d\n", sargs->cores_per_node);
    fprintf(stderr,"Models Mask              (-M): %d\n\n", sargs->models_mask);
  }else{
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp, "\"buf_cnt\": \"%"PRId64"\",\n", sargs->buffer_size);
    fclose(jp);
  }
}


int bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs){

  lgp_init(argc, argv);
  int ret = 0;
  if(MYTHREAD == 0){
    ret = argp_parse(argp, argc, argv, ARGP_NO_EXIT, 0, args);
  }
  
  ret = distribute_cmd_line(argc, argv, &args, arg_len, ret);
  if(ret) return(ret);

  time_t now = time(NULL);
  struct tm *date = localtime(&now);
  if(sargs->json_output){
  /* open the json file */ 
    if(MYTHREAD == 0){
      FILE * fp = fopen(sargs->json_output, "a");
      fprintf(fp, "{\n\"bale_version\": \"%4.2f\",\n", BALE_VERSION);
      fprintf(fp, "\"date\": \"%04d-%02d-%02d.%02d:%02d\",\n",
                date->tm_year+1990, date->tm_mon, date->tm_mday,
                date->tm_hour, date->tm_min);
      fprintf(fp, "\"app\": \"%s\",\n", argv[0]);
      fclose(fp);
    }
  }else{
    T0_fprintf(stderr,"\n***************************************************************\n");
#if __UPC__
    T0_fprintf(stderr,"Bale Version %4.2f (UPC %ld): %04d-%02d-%02d.%02d:%02d\n",
               BALE_VERSION,
               __UPC_VERSION__,
               date->tm_year+1990, date->tm_mon, date->tm_mday,
               date->tm_hour, date->tm_min);
#elif USE_SHMEM
    T0_fprintf(stderr,"Bale Version %4.2f (OpenShmem version %d.%d): %04d-%02d-%02d.%02d:%02d\n",
               BALE_VERSION,
               SHMEM_MAJOR_VERSION, SHMEM_MINOR_VERSION,
               date->tm_year+1990, date->tm_mon+1, date->tm_mday,
               date->tm_hour, date->tm_min); 
#endif
    
    int i;
    
    T0_fprintf(stderr,"Running command on %d PEs:", THREADS);
    for(i=0; i<argc;i++){
      T0_fprintf(stderr," %s", argv[i]);
    }
    T0_fprintf(stderr,"\n");
    T0_fprintf(stderr,"***************************************************************\n\n");
  }
  return(0);
}

void bale_app_finish(std_args_t * sargs){
  if(sargs->json_output && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "r+");
    fseek(jp,-2,SEEK_END);
    fprintf(jp,"\n}\n");
    fclose(jp);
  }
  lgp_finalize();
}

void bale_app_write_time(std_args_t * sargs, char * model_str, double time){
  if(sargs->json_output && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%lf\",\n", model_str, time);
    fclose(jp);
  }else{
    T0_fprintf(stderr, "%s:   %8.3lf\n", model_str, time);
  }
}
