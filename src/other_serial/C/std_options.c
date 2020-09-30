/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For licence information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

#include "std_options.h"

static int std_parse_opt(int key, char * arg, struct argp_state * state) // TODO:doxygen 
{
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'D': args->dump_files = 1; break;
  case 'j': args->json = 1; strcpy(args->json_output,arg); break;
  case 'M': args->models_mask = atol(arg); break;
  case 's': args->seed = atol(arg); break;
  case 'q': args->quiet = 1; break;
  case ARGP_KEY_INIT:
    args->quiet = 0;
    args->seed = 122222;
    args->dump_files = 0;
    args->json = 0;
    break;
  }
  return(0);
}

static struct argp_option std_options[] =  // TODO:doxygen 
{
  {"dump_files", 'D', 0, 0, "Dump files for debugging"},
  {"models_mask", 'M', "MASK", 0, "Which flavors to run."},
  {"seed", 's', "SEED", 0, "Seed for RNG"},
  {"quiet", 'q', 0, 0, "No output during program execution"},
  {0}
};

struct argp std_options_argp =  // TODO:doxygen 
{
  std_options, std_parse_opt, 0, 0, 0
};

// This function writes some of the standard options and their values to the screen (or json file).
void write_std_options(std_args_t * sargs){ // TODO:doxygen 
  if(sargs->json == 0){
    fprintf(stderr,"Standard options:\n");
    fprintf(stderr,"----------------------------------------------------\n");
    fprintf(stderr,"seed                     (-s): %"PRId64"\n", sargs->seed);
    fprintf(stderr,"Models Mask              (-M): %d\n\n", sargs->models_mask);
  }else{
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp, "\"rng_seed\": \"%"PRId64"\",\n", sargs->seed);
    fclose(jp);
  }
}


// the arg_parse function for std_graph_args_t.
static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  switch(key){
  case 'd': args->directed = 1; break;
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; strcpy(args->filename, arg); break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'K':
    args->model = KRONECKER;
    strcpy(args->kron_string,arg);
    char * ptr = arg;
    //T0_fprintf(stderr, "%s\n", arg);
    sscanf(ptr, "%d:", &args->kron_mode);
    //printf("mode %d\n", args->kron_mode);
    int i = 0;
    ptr+=2;
    while(i < 63 && sscanf(ptr, "%d", &args->kron_spec[i])){
      //fprintf(stderr,"%d\n", args->kron_spec[i]);
      i++;
      ptr++; 
      if(*ptr != 'x'){
        break;
      }
      ptr++; 
    }

    args->kron_num = i;
    break;
  case 'l': args->loops = 1; break;
  case 'n': args->numrows = atol(arg); break;
  case 'w': args->weighted = 1; break;
  case 'z': args->nz_per_row = atof(arg); break;
  case ARGP_KEY_INIT:
    args->edge_prob = 0.0;
    args->readfile = 0;
    args->model = FLAT;
    args->numrows = 100;
    args->nz_per_row = 10.0;
    break;
  case ARGP_KEY_END:
    if(args->directed) {
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->numrows,
                                       (args->weighted ? DIRECTED_WEIGHTED : DIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    } else {
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->numrows,
                                       (args->weighted ? UNDIRECTED_WEIGHTED : UNDIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    }
    break;
  }
  return(0);
}

static struct argp_option graph_options[] =
{
  {0, 0, 0, 0, "Input (as file):", 5},
  {"readfile",   'f', "FILE",  0, "Read input from a file"},
  {0, 0, 0, 0, "Input (as random graph):", 6},
  {"numrows",    'n', "NUM",   0, "Number of rows in a matrix"},
  {"directed",   'd', 0,       0, "Specify a directed graph"},
  {"edge_prob",  'e', "EDGEP", 0, "Probability that an edge appears. Use this or -z option to control the density of the graph."},
  {"flat",       'F', 0,       0, "Specify flat random graph model"},
  {"geometric",  'G', 0,       0, "Specify geometric random graph model"},
  {"kronecker",  'K', "KSTR",  0, "Specify a Kronecker product graph.\nKSTR must be a string of the form MODE:S1xS2x...Sk where MODE is 0, 1, or 2 and the Si are small integers that specify the stars whose product is the kronecker product graph. For instance -K 0:3x4x5 specifies MODE 0 and takes the product of K_{1,3}, K_{1,4}, and K_{1,5}. MODE 0 : No triangles. MODE 1: Many triangles. MODE 2: Few triangles. "},
  {"loops",      'l', 0,       0, "Specify you want to force loops into graph"},
  {"weighted",   'w', 0,       0, "Specify you want the edges to be weighted"},
  {"nz_per_row", 'z', "NZPR",  0, "Average number of nonzeros per row. Specify this or -e option to control the density of the graph. Default = 10.0"},
  {0}
};

struct argp std_graph_options_argp =
{
  graph_options, graph_parse_opt, 0, 0, 0
};

/* this function looks at the std_args_t and std_graph_args_t structs to decide
   whether to read a matrix or generate a random matrix. The function returns
   the read or generated matrix. It also writes some matrix statistics to stderr,
   or the output json file.
*/
sparsemat_t * get_input_graph(std_args_t * sargs, std_graph_args_t * gargs) //TODO: move to spmat_utils ??? //TODO: doxygen
//TODO should this be in spmat_utils  maybe need std_enums
{
  sparsemat_t * mat;
  
  if (!gargs->readfile) {
    if (gargs->model == KRONECKER) {
      // generate a random kronecker graph from a string (mode:#x#x...#)
      mat = generate_kronecker_graph_from_spec(gargs->kron_mode, gargs->kron_spec, gargs->kron_num, gargs->weighted);
    } else {
      // Generate a random FLAT or GEOMETRIC graph
      //int64_t numrows = gargs->numrows;
      edge_type et;
      if (gargs->directed) {
        et = (gargs->weighted ? DIRECTED_WEIGHTED : DIRECTED);
      } else {
        et = (gargs->weighted ? UNDIRECTED_WEIGHTED: UNDIRECTED);
      }
      self_loops loops = (gargs->loops ? LOOPS : NOLOOPS);
      
      mat = random_graph(gargs->numrows, gargs->model, et, loops, gargs->edge_prob, sargs->seed);

    }
  } else {
    // Read a matrix from a file
    mat = read_matrix_mm(gargs->filename);
  }
  if (!mat) {fprintf(stderr, "ERROR: get_input_mat: mat is NULL!\n"); exit(1); }

  if (sargs->json == 0) {
    fprintf(stderr,"Input matrix:\n");
    fprintf(stderr,"----------------------------------------------------\n");
    fprintf(stderr,"\t%"PRId64" rows\n\t%"PRId64" columns\n\t%"PRId64" nonzeros\n\n",
                   mat->numcols, mat->numrows, mat->nnz);
  } else {
      FILE * jp = fopen(sargs->json_output, "a");    
      fprintf(jp, "\"matrix_numrows\": \"%"PRId64"\",\n", mat->numrows);
      fprintf(jp, "\"matrix_numcols\": \"%"PRId64"\",\n", mat->numcols);
      fprintf(jp, "\"matrix_nnz\": \"%"PRId64"\",\n", mat->nnz);
      fclose(jp);
  }
  return(mat);
}


// Writes some parameters for the input matrix before it is generated to stderr or an output json file.
void write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs)
{
  if (!gargs->readfile) {
    char model[32];
    if (gargs->model == FLAT)
      sprintf(model, "FLAT        (-F)");
    else if (gargs->model == GEOMETRIC)
      sprintf(model, "GEOMETRIC   (-G)");
    else if (gargs->model == KRONECKER)
      sprintf(model, "KRONECKER   (-K)");
    
    if (sargs->json) {
      FILE * jp = fopen(sargs->json_output, "a");
      fprintf(jp,"\"graph_model\": \"%.*s\",\n", 9, model);
      fclose(jp);
    } else {
      fprintf(stderr,"Input Graph/Matrix parameters:\n");
      fprintf(stderr,"----------------------------------------------------\n");
      fprintf(stderr, "Graph model: %s.\n",model);
      fprintf(stderr,"%s, %s, %s\n",
                 (gargs->directed ? "Directed": "Undirected"),
                 (gargs->weighted ? "Weighted": "Unweighted"),
                 (gargs->loops ? "Loops": "No Loops"));
      fprintf(stderr,"Number of rows           (-N): %"PRId64"\n", gargs->numrows);
      fprintf(stderr,"Avg # nnz per row        (-z): %2.2lf\n", gargs->nz_per_row);
      fprintf(stderr,"Edge probability         (-e): %lf\n\n", gargs->edge_prob);
    }
  } else {
    if (sargs->json) {
      FILE * jp = fopen(sargs->json_output, "a");
      fprintf(jp,"\"input_file\": \"%s\",\n", gargs->filename);
      fclose(jp);
    }else{
      fprintf(stderr,"Reading input from %s\n\n", gargs->filename);
    }
  }
}


// TODO: look at argv for the strings "--help", "--usage", or "-?".
//   If those strings are found, we know, main needs to exit after parsing the command line.
int check_for_exit(int argc, char * argv[], int ret) //TODO: doxygen 
{
  int i;
  for(i = 0; i < argc; i++){
    if ((strcmp(argv[i], "--help") == 0) || 
        (strcmp(argv[i], "-?") == 0) ||
        (strcmp(argv[i], "--usage") == 0)) {
      return 1;
    }
  }

  return(0);
}


/*! \brief This function initializes the environment, parses the command line, and prints
 * some basic app data to stderr or json.
 * 
 * \param argc Standard argc
 * \param argv Standard argv
 * \param args The arguments struct from the app. This struct will be different in general for each app.
 * \param arg_len The sizeof(the arguments struct)
 * \param argp The argp struct which must be initialized in main.
 * \param sargs A pointer to an std_args_t, which must be a field in args.
 * \return 0 if success, < 0 if error, > 0 if main should exit after this call, but main should return(0)
 */
int bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs)
{
    argp_parse(argp, argc, argv, ARGP_NO_EXIT, 0, args);

  // print header 
  time_t now = time(NULL);
  struct tm *date = localtime(&now);
  if (sargs->json) {
      FILE * fp = fopen(sargs->json_output, "w");
      fprintf(fp, "{\n\"C_bale_version\": \"%4.2f\",\n", C_BALE_VERSION); //}
      fprintf(fp, "\"date\": \"%04d-%02d-%02d.%02d:%02d\",\n",
                date->tm_year+1990, date->tm_mon, date->tm_mday,
                date->tm_hour, date->tm_min);
      fprintf(fp, "\"app\": \"%s\",\n", argv[0]);
      fclose(fp);
  } else {
    fprintf(stderr,"\n***************************************************************\n");
    fprintf(stderr,"C Bale Version %4.2f : %04d-%02d-%02d.%02d:%02d\n", C_BALE_VERSION,
                   date->tm_year+1990, date->tm_mon, date->tm_mday, date->tm_hour, date->tm_min);
    
    int i;
    fprintf(stderr,"Running command: ");
    for(i=0; i<argc;i++){
      fprintf(stderr," %s", argv[i]);
    }
    fprintf(stderr,"\n***************************************************************\n\n");
  }
  return(0);
}


/*! \brief Finish off the json output file if appropriate. Call lgp_finalize()
 */
void bale_app_finish(std_args_t * sargs)
{
  if (sargs->json) {    
    FILE * jp = fopen(sargs->json_output, "r+");
    fseek(jp,-2,SEEK_END); //{
    fprintf(jp,"\n}");
    fclose(jp);
  }
}

/*! \brief Write the time of an app implmentation to stderr, or the json file.
 *
 */
void bale_app_write_time(std_args_t * sargs, char * model_str, double time)
{
  if (sargs->json) {    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%lf\",\n", model_str, time);
    fclose(jp);
  }else{
    fprintf(stderr, "%10s: %8.3lf\n", model_str, time);
  }
}

void bale_app_write_int(std_args_t * sargs, char * key, int64_t val)
{
  if (sargs->json) {    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%"PRId64"\",\n", key, val);
    fclose(jp);
  } else {
    fprintf(stderr, "%10s: %"PRId64"\n", key, val);
  }
}
