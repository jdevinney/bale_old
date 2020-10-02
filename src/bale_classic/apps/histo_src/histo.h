/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
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

/*! \file histo.h
 * \brief header file for the histogram app.
 */
#ifndef HISTO_H
#define HISTO_H
#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <locale.h>

// We package all the state up into this struct so we can do some error checking.
typedef struct histo_t{
  SHARED int64_t * counts;
  int64_t * index;
  int64_t * lcounts;
  int64_t * pckindx;
  int64_t num_counts;
  int64_t lnum_counts;
  int64_t l_num_ups;
} histo_t;

double histo_agp(histo_t * data);
double histo_exstack(histo_t * data, int64_t buf_cnt);
double histo_exstack2(histo_t * data, int64_t buf_cnt);
double histo_conveyor(histo_t * data);

#include "alternates/histo_alternates.h"

#endif
