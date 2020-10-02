/******************************************************************
//
//
//  Copyright(C) 2020-2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// 
//
//  All rights reserved.
//  
//   This file is a part of Bale.  For license information see the
//   LICENSE file in the top level directory of the distribution.
//  
// 
 *****************************************************************/ 

/*! \file std_options.h
 * \brief Header file for std_options library. 
 * The std_options library is a support library for all bale apps, mostly for command line parsing.
 */
#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H

#include <stdlib.h>
#include <stdint.h>
#include <spmat.h>
#include <argp.h>

typedef struct std_args_t{
  int64_t buffer_size;
  int cores_per_node;
  int dump_files;
  int json;
  char json_output[128];
  int models_mask;
  int64_t seed;
}std_args_t;

extern struct argp std_options_argp;

typedef struct std_graph_args_t{
  int64_t l_numrows;
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
#define AGP_Model        1L /*!< the Atomic or Generic Interface (straight UPC/SHMEM) */
#define EXSTACK_Model    2L /*!< the exstack bulk synchronous buffering model */
#define EXSTACK2_Model   4L /*!< the exstack2 asynchronous buffering model */
#define CONVEYOR_Model   8L /*!< the conveyor buffering model */
#define ALTERNATE_Model 16L /*!< an alternate model (meant for user supplied code) */
#define ALL_Models      15L /*!< default for running all models */


sparsemat_t *  get_input_graph(std_args_t * sargs, std_graph_args_t * gargs);
void           write_std_options(std_args_t * sargs);
void           write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs);

int  bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs);
void bale_app_finish(std_args_t * sargs);
void bale_app_write_int(std_args_t * sargs, char * key, int64_t val);
void bale_app_write_time(std_args_t * sargs, char * model_str, double time);

#endif


