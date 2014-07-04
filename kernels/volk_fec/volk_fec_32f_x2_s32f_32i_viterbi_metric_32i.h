#ifndef INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_u_H
#define INCLUDED_volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_u_H


static inline void norm(float *metric, int num_states)
{
  float min = metric[0];
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
  \param norm Norm value to normalize state metrics
  \param alpha State Metric Buffer
  \param gamma Branch Metrics Matrix
  \param in_len Number of input symbols
  \param step Current trellisstep
  \param length Trellis length
  \param PS Previous State Matrix
  \param OS Output Matrix
  \param Previous Input Matrix
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
  void *tmp;
  int n, s;
  //for every Infobit = Trellisstep
  metrics_new = &alpha[num_states];
  metrics_old = alpha;
  for(n = 0; n < nibits; n++)
  {
    // Do a Butterfly for 2 States
    for(s= 0; s < num_states/2; s++)
    {
      m0 = metrics_old[s] + gamma[OS[s<<1]+n];
      m1 = metrics_old[s+num_states/2] + gamma[OS[(s+num_states/2)<<1] + n];
      m2 = metrics_old[s] + gamma[OS[s<<1 + 1]+n];
      m3 = metrics_old[s+num_states/2] + gamma[OS[(s+num_states/2)<<1 + 1] + n];
      
      decision0 = (m0 - m1) > 0;
      decision1 = (m2 - m3) > 0;

      metrics_new[2*s] = decision0 ? m0 : m1;
      metrics_new[2*s+1] = decision0 ? m0 : m1;
    }

    norm(metrics_new, num_states);
    tmp = (void*) metrics_old;
    metrics_old = metrics_new;
    metrics_new = (float*) tmp;
  }

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
  \param norm Norm value to normalize state metrics
  \param alpha State Metric Buffer
  \param gamma Branch Metrics Matrix
  \param in_len Number of input symbols
  \param step Current trellisstep
  \param length Trellis length
  \param PS Previous State Matrix
  \param OS Output Matrix
  \param Previous Input Matrix
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

static inline void volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_a_sse4(unsigned int *trace,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          const int length,
                                                          const int *OS,
                                                          unsigned int num_points)
{

}

#endif /* LV_HAVE_SSE */
#endif /* INCLUDED_volk_fec_32f_x2_32i_x2_viterbi_metric_32i_s32f_a_H */