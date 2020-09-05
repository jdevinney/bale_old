#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H
#include <argp.h>
#include <stdlib.h>
#include <stdint.h>
#include <libgetput.h>
#include "spmat_enums.h"

typedef struct std_args_t{
  int64_t buffer_size;
  int cores_per_node;
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


/* defines to support the different models of global and buffered references */
#define AGI_Model        1L /*!< the Atomic or Generic Interface (straight UPC/SHMEM) */
#define EXSTACK_Model    2L /*!< the exstack bulk synchronous buffering model */
#define EXSTACK2_Model   4L /*!< the exstack2 asynchronous buffering model */
#define CONVEYOR_Model   8L /*!< the conveyor buffering model */
#define ALTERNATE_Model 16L /*!< an alternate model (meant for user supplied code) */
#define ALL_Models      15L /*!< default for running all models */


void share_args(void * args, size_t n);
int check_for_exit(int argc, char * argv[], int ret);
#endif
