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

/*! \file ig.h
 * \brief Demo program that computes an indexed gather of elements from a 
 *  shared array into local arrays.
 *  The size of the source array should be large enough that the elements 
 *  need to be spread across the whole machine
 */

#include <exstack.h>
#include <convey.h>

double ig_gets(int64_t *tgt, int64_t *index, int64_t T,  SHARED int64_t *table);
double ig_exstack(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable, int64_t buf_cnt); 
double ig_exstack2(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable, int64_t buf_cnt); 
double ig_conveyor(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable);
double ig_exstack2_cyclic(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable, int64_t buf_cnt);
double ig_exstack2_goto(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable, int64_t buf_cnt);
double ig_exstack_function(int64_t *tgt, int64_t *pckindx, int64_t T,  int64_t *ltable, int64_t buf_cnt);

