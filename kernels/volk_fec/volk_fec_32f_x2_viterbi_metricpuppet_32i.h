#ifndef INCLUDED_volk_fec_32f_x2_viterbi_metricpuppet_32i_u_H
#define INCLUDED_volk_fec_32f_x2_viterbi_metricpuppet_32i_u_H

#include <inttypes.h>
#include <stdio.h>

#include <volk_fec/volk_fec.h>
#include <volk_fec/volk_fec_32f_x2_s32f_32i_viterbi_metric_32i.h>

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
static inline void volk_fec_32f_x2_viterbi_metricpuppet_32i_generic(unsigned int *trace,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_manual(trace, alpha, gamma, 
                                                     num_points/32, &OS[0], 16, "generic"); 
}


#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32f_x2_viterbi_metric_32i_u_H */

#ifndef INCLUDED_volk_fec_32f_x2_viterbi_metricpuppet_32i_a_H
#define INCLUDED_volk_fec_32f_x2_viterbi_metricpuppet_32i_a_H

#include <inttypes.h>
#include <stdio.h>

#include <volk_fec/volk_fec.h>
#include <volk_fec/volk_fec_32f_x2_s32f_32i_viterbi_metric_32i.h>

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

static inline void volk_fec_32f_x2_viterbi_metricpuppet_32i_a_sse4(unsigned int *trace,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  volk_fec_32f_x2_s32f_32i_viterbi_metric_32i_manual(trace, alpha, gamma, 
                                                     num_points/32, &OS[0], 16, "a_sse4");
}

#endif /* LV_HAVE_SSE */
#endif /* INCLUDED_volk_fec_32f_x2_viterbi_metric_32i_a_H */