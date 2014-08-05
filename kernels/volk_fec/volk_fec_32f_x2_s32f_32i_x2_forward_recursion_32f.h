#ifndef INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H
#define INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H


#include <stdlib.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC

static inline void volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f_generic(float *alpha,
                                              float *gamma, const int length, 
                                              const int *OS, const int *shuffle,
                                              unsigned int num_points)
{
  int n = 0;
  int s = 0;
  float m0,m1,m2,m3;
  int t[4];
  for(n = 0; n < length; n++)
  {
    for(s = 0; s < num_points/2; s++)
    {
      // Get shuffled index for branch transitions
      t[0] = OS[4*s+shuffle[4*s]];
      t[1] = OS[4*s+2+shuffle[4*s+2]];
      t[2] = OS[4*s+shuffle[4*s+1]];
      t[3] = OS[4*s+2+shuffle[4*s+3]];

      // Calc state metrics
      m0 = alpha[n*num_points+2*s  ] + gamma[n*4+t[0]];
      m1 = alpha[n*num_points+2*s+1] + gamma[n*4+t[1]];
      m2 = alpha[n*num_points+2*s  ] + gamma[n*4+t[2]];
      m3 = alpha[n*num_points+2*s+1] + gamma[n*4+t[3]];

      // Chose the one with maximum correlation/probabiliy
      alpha[(n+1)*num_points+s] = ((m1-m0>0) ? m1 : m0);
      alpha[(n+1)*num_points+s+num_points/2] = ((m3-m2>0) ? m3 : m2); 
    }
  } 
}

#endif /*LV_HAVE_GENERIC*/
#endif /*INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H*/

#ifndef INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H
#define INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"
static inline void volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f_a_sse4(float *alpha,
                                              float *gamma, const int length, 
                                              const int *OS, const int *shuffle,
                                              unsigned int num_points)
{
  __m128 m[4];
  __m128 metrics_new128[2], metrics_old128[2];
  __m128 gamma128[4];
  float *metrics_new, *metrics_old;
  int t[16];
  int n,s,i; // counting variables
  int quarter_points = num_points/4;

  metrics_new = &alpha[num_points];
  metrics_old = &alpha[0];
  for(n = 0; n < length; n++)
  {
    for(s = 0; s < quarter_points/2; s++)
    {
      gamma128[0] = _mm_set_ps(gamma[OS[16*s+12+shuffle[16*s+12]]+n*4], 
                               gamma[OS[16*s+8+shuffle[16*s+8]]+n*4],
                               gamma[OS[16*s+4+shuffle[16*s+4]]+n*4], 
                               gamma[OS[16*s+shuffle[16*s]]+n*4]);
      
      gamma128[1] = _mm_set_ps(gamma[OS[16*s+14+shuffle[16*s+14]]+n*4],
                               gamma[OS[16*s+10+shuffle[16*s+10]]+n*4], 
                               gamma[OS[16*s+6+shuffle[16*s+6]]+n*4], 
                               gamma[OS[16*s+2+shuffle[16*s+2]]+n*4]);
      
      gamma128[2] = _mm_set_ps(gamma[OS[16*s+12+shuffle[16*s+13]]+n*4], 
                               gamma[OS[16*s+8+shuffle[16*s+9]]+n*4],
                               gamma[OS[16*s+4+shuffle[16*s+5]]+n*4], 
                               gamma[OS[16*s+shuffle[16*s+1]]+n*4]);
      
      gamma128[3] = _mm_set_ps(gamma[OS[16*s+14+shuffle[16*s+15]]+n*4], 
                               gamma[OS[16*s+10+shuffle[16*s+11]]+n*4],
                               gamma[OS[16*s+6+shuffle[16*s+7]]+n*4], 
                               gamma[OS[16*s+2+shuffle[16*s+3]]+n*4]);

      metrics_old128[0] = _mm_set_ps(metrics_old[8*s +6 ], metrics_old[8*s+4], 
                                     metrics_old[8*s+2], metrics_old[8*s]);
      metrics_old128[1] = _mm_set_ps(metrics_old[8*s+7], metrics_old[8*s+5], 
                                     metrics_old[8*s+3], metrics_old[8*s+1]);

      m[0] = _mm_add_ps(metrics_old128[0], gamma128[0]);
      m[1] = _mm_add_ps(metrics_old128[1], gamma128[1]);
      m[2] = _mm_add_ps(metrics_old128[0], gamma128[2]);
      m[3] = _mm_add_ps(metrics_old128[1], gamma128[3]);

      metrics_new128[0] = _mm_max_ps(m[1], m[0]); 
      metrics_new128[1] = _mm_max_ps(m[3], m[2]);

      _mm_store_ps(&metrics_new[4*s], metrics_new128[0]);
      _mm_store_ps(&metrics_new[4*s+num_points/2], metrics_new128[1]); 

    }
    metrics_old += num_points;
    metrics_new += num_points;
  }


}
#endif /*LV_HAVE_SSE4_1*/
#endif /* INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H */