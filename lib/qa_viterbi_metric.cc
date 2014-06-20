/* -*- c++ -*- */
/*
 * Copyright 2012 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */
 
#include <stdio.h>
#include <cppunit/TestAssert.h>
#include <gnuradio/attributes.h>
#include <gnuradio/random.h>
#include <volk_fec/volk_fec.h>

void qa_viterbi_metric_fi::t1()
{
  volk_fec_environment_init();
  clock_t start, end;
  double total;
  int vlen = 2048;
  int iters = 1000;
  std::vector<int> trace(S*K);
  size_t align = volk_fec_get_alignment();
  float *alpha;
  alpha = (float*)volk_fec_malloc(sizeof(float)*2*S, align);
  size_t align = volk_get_alignment();
  int PS[] = {0,1,2,3,0,1,2,3};
  int PI[] = {0,1,2,3,0,1,2,3};
  int OS[] = {0,1,2,3,0,1,2,3}; 
  float *in;
  rx_samples = (float*) volk_malloc(sizeof(float)*num_symbols, align);
  for(int i = 0; i < num_symbols; i++)
  {
    // Create random samples ranging from -4 to 4
    rx_samples[i] = (float) (((float)::random()/RANDOM_MAX-0.5)*10);
  }

  for(int i = 0; i < num_runs; i++)
  {
    norm = 1e9;
    clock_t t_volk, t;
    t_volk = clock();
    //volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32i_32f_manual(&trace[0], &norm, &alpha[0], &in[0], 2, 0, 1,&PS[0], 
    //                                &OS[0], &PI[0], S, "sse");
    t_volk = clock() - t_volk;
    time_volk += ((float) t_volk);

    // t = clock();
    // volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32i_32f_manual(&trace[0], &norm, &alpha[0], &in[0], 2, 0, 1,&PS[0], 
    //                                 &OS[0], &PI[0], S, "generic");
    // t = clock() - t;
    // time_generic += ((float) t);
  }
}

