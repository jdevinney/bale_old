/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
// 
//
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
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


