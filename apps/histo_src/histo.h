/******************************************************************
 * Copyright 2014, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 * POC: Bale <bale@super.org>
 * Please contact the POC before disseminating this code.
 *****************************************************************/ 

/*! \file histo.h
 * \brief header file for the histogram app.
 */

#include <exstack.h>
#include <convey.h>
#include <locale.h>

double histo_atomic(int64_t *index, int64_t T,  SHARED int64_t *counts, int64_t v);
double histo_conveyor(int64_t *pckindx, int64_t T,  int64_t *lcounts); 
double histo_exstack(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);
double histo_exstack2(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);
double histo_exstack2_goto(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);
double histo_exstack2_cyclic(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);
double histo_exstack_function(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);
double histo_exstack2_function(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t bufsiz);

