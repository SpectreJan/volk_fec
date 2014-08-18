#ifndef INCLUDED_VOLK_FEC_32f_calc_branch_metricpuppet_32f_u_H
#define INCLUDED_VOLK_FEC_32f_calc_branch_metricpuppet_32f_u_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef LV_HAVE_GENERIC


static inline void volk_fec_32f_calc_branch_metricpuppet_32f_generic(const float *input_llr,
                                                                    float *log_probabilities,
                                                                    const int num_points)
{
  volk_fec_32f_s32f_32i_calc_branch_metric_32f_manual(input_llr, 2, log_probabilities, 
                                                      num_points/4, "generic");
} 

#endif /*LV_HAVE_GENERIC*/
#endif /*INCLUDED_VOLD_FEC_32f_calc_branch_metric_32f_u_H*/

#ifndef INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metricpuppet_32f_a_H
#define INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metricpuppet_32f_a_H

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"

static inline void volk_fec_32f_calc_branch_metricpuppet_32f_a_sse4(const float *input_llr,
                                                                    float *log_probabilities,
                                                                    const int num_points)
{
  volk_fec_32f_s32f_32i_calc_branch_metric_32f_manual(input_llr, 2, log_probabilities, 
                                                      num_points/4, "a_sse4");
}
#endif /* LV_HAVE_SSE4_1 */
#endif /* INCLUDED_VOLD_FEC_32f_calc_branch_metricpuppet_32f_a_H */