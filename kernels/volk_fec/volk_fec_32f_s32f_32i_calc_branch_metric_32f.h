#ifndef INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metric_32f_u_H
#define INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metric_32f_u_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef LV_HAVE_GENERIC


static inline void volk_fec_32f_s32f_32i_calc_branch_metric_32f_generic(const float *input_llr,
                                                                    const int n_cbit,
                                                                    float *log_probabilities,
                                                                    const int num_points)
{
  int o = 0;
  int n = 0;    
  int k = 0;

  int O = 1<<(n_cbit);
  for(k = 0; k < num_points; k++)
  { 
    //for all possible encoder outputs

    for(o = 0; o < O/2; o++)
    {
      log_probabilities[k*(1<<n_cbit)+o] = 0;
      for(n = 0; n < n_cbit; n++)
      {
        log_probabilities[k*O+o] += 0.5*input_llr[k*n_cbit+n]*(2*(o&(1<<(n_cbit-n-1)))-1);
      }
    }
    for(o = 0; o < O/2; o++)
    {
      log_probabilities[k*O+o+O/2] = -log_probabilities[k*O+O/2-o-1];
    }
  }

} 
#endif /*LV_HAVE_GENERIC*/
#endif /*INCLUDED_VOLD_FEC_32f_calc_branch_metric_32f_u_H*/

#ifndef INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H
#define INCLUDED_VOLK_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"

static inline void volk_fec_32f_s32f_32i_calc_branch_metric_32f_a_sse4(const float *input_llr,
                                                                    const int n_cbit,
                                                                    float *log_probabilities,
                                                                    const int num_points)
{
  int o = 0;
  int n = 0;
  int k = 0;
  int output_iter = (1<<n_cbit)>>2;
  float *out_ptr = log_probabilities;

  __m128 log_prob128, input_llr128, codebit128;
  
  // For all input LLRs
  for(k = 0; k < num_points; k++)
  {
    // For all possible encoder outputs/codewords
    for(o = 0; o < output_iter; o++)
    {
      log_prob128 = _mm_setzero_ps();
      // For all codebits in a codeword
      for(n = 0; n < n_cbit; n++)
      {
        // Load input LLR
        input_llr128 = _mm_set1_ps(input_llr[k*n_cbit+n]);
        // Set codewordbits
        codebit128 = _mm_set_ps(2*(((4*o+3)>>(n_cbit-n-1))&1)-1,
                                2*(((4*o+2)>>(n_cbit-n-1))&1)-1,
                                2*(((4*o+1)>>(n_cbit-n-1))&1)-1,
                                2*(((4*o)>>(n_cbit-n-1))&1)-1);

        // Apply Hanzo equation 14.13 (page 564), basically you correlate your input with all possible codewords
        log_prob128 = _mm_add_ps(log_prob128, _mm_mul_ps(_mm_set1_ps(0.5), _mm_mul_ps(codebit128, input_llr128)));
      }
      // Store results and move pointer
      _mm_store_ps(out_ptr, log_prob128);
      out_ptr += 4;  
    }
  }
}
#endif /* LV_HAVE_SSE4_1 */
#endif /* INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H */