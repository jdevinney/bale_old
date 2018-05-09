// Copyright (c) 2018, Institute for Defense Analyses,
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
//
// This material may be reproduced by or for the U.S. Government 
// pursuant to the copyright license under the clauses at DFARS 
// 252.227-7013 and 252.227-7014.

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "convey.h"
#if HAVE_CONFIG_H
# include "config.h"
#endif
#include "bolite.h"

#if MPP_USE_UPC
# define PROCS THREADS
# define MY_PROC MYTHREAD
# define example_start()
# define example_end()
#elif !HAVE_MPP_UTIL
# include "shmem.h"
# define PROCS shmem_n_pes()
# define MY_PROC shmem_my_pe()
# define example_start() shmem_init()
# define example_end() shmem_finalize()
#else
# include "mpp_utilV4.h"
# define example_start() argc = mpp_util_init(argc, argv, NULL)
# define example_end() mpp_util_fin()
#endif
