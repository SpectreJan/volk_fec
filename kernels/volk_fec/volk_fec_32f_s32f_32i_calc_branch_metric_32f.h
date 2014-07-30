#ifndef INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_u_H
#define INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_u_H

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
  for(; k < num_points; k++)
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