#ifndef INCLUDED_volk_fec_32f_x2_s32f_x3_32i_x3_s32f_viterbi_metric_32f_u_H
#define INCLUDED_volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32f_u_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC

static inline void volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32i_32f_generic(unsigned int *trace,
                                                          float *norm,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          const int in_len,
                                                          const int step,
                                                          const int length,                                                        
                                                          const int *PS,
                                                          const int *OS,
                                                          const int *PI,
                                                          unsigned int num_points)
{

  float m;
  float min_tmp;
  int minmi;
  unsigned int s = 0;
  unsigned int i = 0;

  int alphai = step%2;
  for(s = 0; s < num_points; s++)
  { 
    minmi = 0;
    min_tmp = 1e9;
    for(i = 0; i < in_len; i++)
    {
      m = alpha[alphai*num_points+PS[s*in_len+i]] + gamma[OS[PS[s*in_len+i]*in_len+PI[s*in_len+i]]*length+step];
      if(m < min_tmp)
      {
        min_tmp = m;
        minmi = i; 
      } 
    }

    trace[step*num_points + s]= minmi;
    alpha[((alphai+1)%2)*num_points+s] = min_tmp;
    if(min_tmp<norm[0]) norm[0]=min_tmp;
  }
}


#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32f_x2_32i_x3_viterbi_metric_32i_s32f_u_H */

#ifndef INCLUDED_volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32f_a_H
#define INCLUDED_volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32f_a_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"

static inline void volk_fec_32f_x2_s32f_x3_32i_x3_viterbi_metric_32i_32f_a_sse4(unsigned int *trace,
                                                          float *norm,
                                                          float *alpha, 
                                                          const float *gamma,
                                                          const int in_len,
                                                          const int step,
                                                          const int length,                                                        
                                                          const int *PS,
                                                          const int *OS,
                                                          const int *PI,
                                                          unsigned int num_points)
{
  unsigned int s = 0;
  unsigned int i = 0;
  unsigned int quarter_points = num_points/4;
  unsigned int *trace_ptr = (int*) &trace[step*num_points]; 
  float norm_buf[4];
  __m128 minm, minm_tmp, gamma128, alpha128, norm128, input;
  __m128 minmi_tmp;
  __m128i minmi;

  norm128 = _mm_load_ps1(norm);
  int alphai = step%2;
  
  for(s = 0; s < quarter_points; s++)
  {
    //minmi_tmp = _mm_set_ps(0,0,0,0);
    minm = _mm_set_ps(1e9,1e9,1e9,1e9);
    for(i = 0; i < in_len; i++)
    {
      
     gamma128 = _mm_set_ps(gamma[OS[PS[((s*4+3)*in_len)+i]*in_len+PI[((s*4+3)*in_len)+i]]*length+step],
                           gamma[OS[PS[((s*4+2)*in_len)+i]*in_len+PI[((s*4+2)*in_len)+i]]*length+step],
                           gamma[OS[PS[((s*4+1)*in_len)+i]*in_len+PI[((s*4+1)*in_len)+i]]*length+step],
                           gamma[OS[PS[((s)*in_len*4)+i]*in_len+PI[(s*in_len*4)+i]]*length+step]);
      
      alpha128 = _mm_set_ps(alpha[alphai*num_points+PS[((s*4+3)*in_len)+i]],
                            alpha[alphai*num_points+PS[((s*4+2)*in_len)+i]],
                            alpha[alphai*num_points+PS[((s*4+1)*in_len)+i]],
                            alpha[alphai*num_points+PS[(s*in_len*4)+i]]);
      
      input = _mm_set_ps1(i);
      minm_tmp = _mm_add_ps(alpha128, gamma128);
      minmi_tmp = _mm_cmpgt_ps(minm, minm_tmp);
      minmi_tmp = _mm_min_ps(minmi_tmp, input);
      minm = _mm_min_ps(minm, minm_tmp);

            
    }
    
    norm128 = _mm_min_ps(minm, norm128);
    _mm_store_ps(&alpha[((alphai+1)%2)*num_points+s*4], minm);

    minmi = _mm_cvtps_epi32(minmi_tmp);
    _mm_store_si128((__m128i*)trace_ptr, minmi);
    trace_ptr += 4;
//
  }
  _mm_store_ps(norm_buf, norm128);
  for(i = 0; i < 4; i++)
  {
    if(norm_buf[i] < norm[0]) norm[0] = norm_buf[i];
  }

  s = quarter_points*4;
  for(; s < num_points; s++)
  { 
    int minmi_l = 0;
    int min_tmp = 1e9;
    int m;
    for(i = 0; i < in_len; i++)
    {
      m = alpha[alphai*num_points+PS[s*in_len+i]] + gamma[OS[PS[s*in_len+i]*in_len+PI[s*in_len+i]]*length+step];
      if(m < min_tmp)
      {
        min_tmp = m;
        minmi_l = i; 
      } 
    }

    trace[step*num_points + s]= minmi_l;
    alpha[((alphai+1)%2)*num_points+s] = min_tmp;
    if(min_tmp<norm[0]) norm[0]=min_tmp;
  }

}

#endif /* LV_HAVE_SSE */
#endif /* INCLUDED_volk_fec_32f_x2_32i_x3_viterbi_metric_32i_s32f_a_H */