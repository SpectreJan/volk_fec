#ifndef INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_u_H
#define INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_u_H


static inline void norm(float *metric, int num_states)
{
  float min = metric[0];
  int minmi = 0;
  int s;
  for(s = 1; s < num_states; s++)
  {
    if (metric[s] < min)
    {
      min = metric[s];
    }
  }
  for(s = 0; s < num_states; s++)
  {
    metric[s] -= min;
  }
}

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC
/*!
  \brief Update viterbi state metric
  \param trace Trace Matrix
  \param alpha State Metric Buffer
  \param gamma Branch Metrics Matrix
  \param length Trellis length
  \param OS Output Matrix
  \param num_points The number of States in the trellis
*/
/* Update state metric and use butterfly structure

a_n[s    ]******a_{n+1}[s]
           *   *
            * *
             *
            * *
           *   *
a_n[s+S/2]******a_{n+1}[s+1]
*/
static inline void volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_generic(unsigned int *trace,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          const int length,                                                        
                                                          const int *OS,
                                                          unsigned int num_points)
{
  int nibits = length;
  int num_states = num_points;
  int decision0, decision1; // Decision that is stored in Trace Pointer
  float m0,m1,m2,m3;
  float *metrics_new, *metrics_old;
  float *tmp;
  int n, s;
  
  //for every Infobit = Trellisstep
  metrics_new = &alpha[num_states];
  metrics_old = &alpha[0];
  //printf("NumStates: %d\nInfobits: %d\n", num_states, nibits);
  for(n = 0; n < nibits; n++)
  {
    // Do a Butterfly for 2 States
    for(s= 0; s < num_states/2; s++)
    {
      m0 = metrics_old[2*s  ] + gamma[OS[4*s]*nibits+n];
      m1 = metrics_old[2*s+1] + gamma[OS[4*s+2]*nibits+n];
      m2 = metrics_old[2*s  ] + gamma[OS[4*s+1]*nibits+n];
      m3 = metrics_old[2*s+1] + gamma[OS[4*s+3]*nibits+n];
      
      //printf("%f %f %f %f\n", m0, m1, m2, m3);
      // Get decision for trace based on minimum euclidean distance
      decision0 = (m0 - m1) > 0;
      decision1 = (m2 - m3) > 0;

      // Save metric for next iteration
      metrics_new[s] = decision0 ? m1 : m0;
      metrics_new[s+num_states/2] = decision1 ? m3 : m2;

      // Save trace by packing 8 decision into one char
      trace[n*num_states+s]   = decision0;
      trace[n*num_states+s+num_states/2] = decision1;
    }

    // Normalize metrics and swap pointers for next trellisstep
    norm(metrics_new, num_states);
    tmp = metrics_old;
    metrics_old = metrics_new;
    metrics_new = tmp;
  }
  alpha = &metrics_old[0];
}


#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32f_x2_32i_x3_viterbi_metric_32i_s32f_u_H */

#ifndef INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32f_a_H
#define INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32f_a_H

//#include <time.h>
struct timespec start_loop, end_loop, start_instr, end_instr, res; 
double dur_loop, dur_instruction;

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"

/*!
  \brief Update viterbi state metric
  \param trace Trace Matrix
  \param alpha State Metric Buffer
  \param gamma Branch Metrics Matrix
  \param length Trellis length
  \param OS Output Matrix
  \param num_points The number of States in the trellis
*/

/* Update state metric and use butterfly structure

a_n[2*s    ]******a_{n+1}[s]
           *   *
            * *
             *
            * *
           *   *
a_n[2*s+1]******a_{n+1}[s+num_states/2]
*/

static inline void volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_a_sse4(unsigned int *trace,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          const int length,
                                                          const int *OS,
                                                          unsigned int num_points)
{
  int nibits = length;
  int quarter_points = num_points/4;
  __m128i decision0[1], decision1[1]; // Decision that is stored in Trace Pointer
  __m128 m[4];
  __m128 metrics_new128[2], metrics_old128[2];
  __m128 gamma128[4];
  void *tmp;
  float *metrics_new, *metrics_old;
  unsigned int *trace_ptr;
  int n, s;

  metrics_new = &alpha[num_points];
  metrics_old = &alpha[0];
  trace_ptr = &trace[0];
  for(n = 0; n < nibits; n++)
  {
    for(s = 0; s < quarter_points/2; s++)
    {
      // Set all registers for branch metrics
      gamma128[0] = _mm_set_ps(gamma[OS[16*s+12]*nibits+n], 
                               gamma[OS[16*s+8 ]*nibits+n],
                               gamma[OS[16*s+4 ]*nibits+n], 
                               gamma[OS[16*s   ]*nibits+n]);
      gamma128[1] = _mm_set_ps(gamma[OS[16*s+14]*nibits+n],
                               gamma[OS[16*s+10]*nibits+n], 
                               gamma[OS[16*s+6 ]*nibits+n], 
                               gamma[OS[16*s+2 ]*nibits+n]);
      gamma128[2] = _mm_set_ps(gamma[OS[16*s+13]*nibits+n], 
                               gamma[OS[16*s+9 ]*nibits+n],
                               gamma[OS[16*s+5 ]*nibits+n], 
                               gamma[OS[16*s+1 ]*nibits+n]);
      gamma128[3] = _mm_set_ps(gamma[OS[16*s+15]*nibits+n], 
                               gamma[OS[16*s+11]*nibits+n],
                               gamma[OS[16*s+7 ]*nibits+n], 
                               gamma[OS[16*s+3 ]*nibits+n]);

      metrics_old128[0] = _mm_set_ps(metrics_old[8*s +6 ], metrics_old[8*s+4], 
                                     metrics_old[8*s+2], metrics_old[8*s]);
      metrics_old128[1] = _mm_set_ps(metrics_old[8*s+7], metrics_old[8*s+5], 
                                     metrics_old[8*s+3], metrics_old[8*s+1]);
      
      // Calculate path metrics
      m[0] = _mm_add_ps(metrics_old128[0], gamma128[0]);
      m[1] = _mm_add_ps(metrics_old128[1], gamma128[1]);
      m[2] = _mm_add_ps(metrics_old128[0], gamma128[2]);
      m[3] = _mm_add_ps(metrics_old128[1], gamma128[3]);

      // Calculate Decisions and cast to __m128i
      *decision0 = (__m128i) _mm_cmpgt_ps(_mm_sub_ps(m[0], m[1]), _mm_setzero_ps());
      *decision1 = (__m128i) _mm_cmpgt_ps(_mm_sub_ps(m[2], m[3]), _mm_setzero_ps());

      // Is m[i] or m[i+1] minimum, calculated by complicated bitshift operatiońs to gain speed (see Karn)
      metrics_new128[0] = _mm_or_ps(_mm_and_ps((__m128)*decision0, m[1]), _mm_andnot_ps((__m128)*decision0,m[0])); 
      metrics_new128[1] = _mm_or_ps(_mm_and_ps((__m128)*decision1, m[3]), _mm_andnot_ps((__m128)*decision1,m[2]));
      *decision0 = _mm_and_si128(*decision0, _mm_set1_epi32(1));
      *decision1 = _mm_and_si128(*decision1, _mm_set1_epi32(1));


      // Store decisions and path metrics
      _mm_store_si128((__m128i*)&trace_ptr[n*num_points+4*s], *decision0);
      _mm_store_si128((__m128i*)&trace_ptr[n*num_points+4*s+num_points/2], *decision1);
        
      _mm_store_ps(&metrics_new[4*s], metrics_new128[0]);
      _mm_store_ps(&metrics_new[4*s+num_points/2], metrics_new128[1]);     
    }

    // Normalize metrics and swap pointers for next trellisstep
    norm(metrics_new, num_points);
    tmp = (void*) metrics_old;
    metrics_old = metrics_new;
    metrics_new = (float*) tmp;  
  }
  alpha = metrics_new;
}

#endif /* LV_HAVE_SSE */
#endif /* INCLUDED_volk_fec_32f_x2_32i_x2_viterbi_metric_32i_s32f_a_H */