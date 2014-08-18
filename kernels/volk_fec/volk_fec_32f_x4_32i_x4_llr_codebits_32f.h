#ifndef INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_u_H
#define INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_u_H

#define INF 1e9
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef LV_HAVE_GENERIC

static inline void volk_fec_32f_x4_32i_x4_llr_codebits_32f_generic(float *alpha, float *gamma,
                                                                float *beta, float *llr,
                                                                const int n_cbit, const int length,
                                                                const int *OS, const int *shuffle,
                                                                unsigned int num_points)
{
  int O =  1<<n_cbit;
  float m0,m1,m2,m3, sum_tmp0, sum_tmp1, sum_tmp2, sum_tmp3;
  int t0,t1,t2,t3, cbit0, cbit1, cbit2, cbit3;

  float *metrics;
  metrics = (float*) malloc(sizeof(float)*4);

  float *beta_new = &beta[num_points];
  float *beta_old = &beta[0];
  float *beta_tmp = NULL;
  
  int n = 0;
  int s = 0;
  int c = 0;
 
  for(n = length-1; n >= 0; n--)
  {
    for(c = 0; c < n_cbit*2; c++)
    {
      metrics[c] = -INF;
    }
    for(s = 0; s < num_points/2; s++)
    {
      // Get shuffled index for branch transitions
      t0 = OS[4*s+shuffle[4*s]];
      t1 = OS[4*s+shuffle[4*s+1]];
      t2 = OS[4*s+2+shuffle[4*s+2]];
      t3 = OS[4*s+2+shuffle[4*s+3]];



      // Calc state metrics
      m0 = beta_old[s] + gamma[n*O+t0];
      m1 = beta_old[s+num_points/2] + gamma[n*4+t1];
      m2 = beta_old[s] + gamma[n*O+t2];
      m3 = beta_old[s+num_points/2] + gamma[n*4+t3];

      // Chose metric with highest correlation/probability
      beta_new[2*s] = ((m1-m0>0) ? m1 : m0);
      beta_new[2*s+1] = ((m3-m2>0) ? m3 : m2);

      sum_tmp0 = alpha[n*num_points+2*s]+m0;
      sum_tmp1 = alpha[n*num_points+2*s]+m1;
      sum_tmp2 = alpha[n*num_points+2*s+1]+m2;
      sum_tmp3 = alpha[n*num_points+2*s+1]+m3;

      for(c = 0; c < n_cbit; c++)
      {
        cbit0 = (t0>>(n_cbit-c-1))&1;
        cbit1 = (t1>>(n_cbit-c-1))&1;
        cbit2 = (t2>>(n_cbit-c-1))&1;
        cbit3 = (t3>>(n_cbit-c-1))&1;
        
        metrics[2*c+cbit0] = (metrics[2*c+cbit0] - sum_tmp0 > 0) ? metrics[2*c+cbit0] : sum_tmp0;
        metrics[2*c+cbit1] = (metrics[2*c+cbit1] - sum_tmp1 > 0) ? metrics[2*c+cbit1] : sum_tmp1;
        metrics[2*c+cbit2] = (metrics[2*c+cbit2] - sum_tmp2 > 0) ? metrics[2*c+cbit2] : sum_tmp2;
        metrics[2*c+cbit3] = (metrics[2*c+cbit3] - sum_tmp3 > 0) ? metrics[2*c+cbit3] : sum_tmp3;
      } 
    }
      for(c = 0; c < n_cbit; c++)
      {
        llr[n*n_cbit+c] = metrics[2*c+1]-metrics[2*c];
      }
      beta_tmp = beta_old;
      beta_old = beta_new;
      beta_new = beta_tmp;
  }
  free(metrics);
}

#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_u_H */

#ifndef INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_a_H
#define INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_a_H
#define INF 1e9
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <volk_fec/volk_fec.h>

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"
static inline void volk_fec_32f_x4_32i_x4_llr_codebits_32f_a_sse4(float *alpha, float *gamma,
                                                                float *beta, float *llr,
                                                                const int n_cbit, const int length,
                                                                const int *OS, const int *shuffle,
                                                                unsigned int num_points)
{
  __m128 m[4], m01, m23;
  __m128 beta_new128[2], beta_old128[2], gamma128[4];
  __m128 beta_tmp1, beta_tmp2;
  __m128 *alpha128_1;
  __m128 *alpha128_0; 
  __m128 *sum_tmp0;
  __m128 *sum_tmp1;
  __m128 *sum_tmp2;
  __m128 *sum_tmp3;
  float *metrics;
  metrics = (float*) volk_fec_malloc(sizeof(float)*2*n_cbit, volk_fec_get_alignment());
  float sum[16];
  sum_tmp0 = (__m128*) &sum[0];
  sum_tmp1 = (__m128*) &sum[4];
  sum_tmp2 = (__m128*) &sum[8];
  sum_tmp3 = (__m128*) &sum[12];
  float *beta_new = &beta[num_points];
  float *beta_old = &beta[0];
  float *beta_tmp = NULL;

   
  int t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;
  int b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15;
  int n = 0;
  int s = 0;
  int c = 0;
  int i = 0;

  int quarter_points = num_points/4;

  for(n = length-1; n >= 0; n--)
  {
    for(c = 0; c < n_cbit*2; c++)
      {
        metrics[c] = -INF;
      }
    for(s = 0; s < quarter_points/2; s++)
    {
      // Set the right transition and load up the registers
      t0  = OS[16*s+shuffle[16*s]];
      t1  = OS[16*s+4+shuffle[16*s+4]];
      t2  = OS[16*s+8+shuffle[16*s+8]];  
      t3  = OS[16*s+12+shuffle[16*s+12]];

      gamma128[0] = _mm_set_ps(gamma[t3+n*4], 
                               gamma[t2+n*4],
                               gamma[t1+n*4], 
                               gamma[t0+n*4]);

      t4  = OS[16*s+shuffle[16*s+1]];
      t5  = OS[16*s+4+shuffle[16*s+5]];
      t6  = OS[16*s+8+shuffle[16*s+9]];
      t7  = OS[16*s+12+shuffle[16*s+13]];

      gamma128[1] = _mm_set_ps(gamma[t7+n*4], 
                               gamma[t6+n*4],
                               gamma[t5+n*4], 
                               gamma[t4+n*4]);

      t8  =  OS[16*s+2+shuffle[16*s+2]];
      t9  = OS[16*s+6+shuffle[16*s+6]];
      t10 = OS[16*s+10+shuffle[16*s+10]];
      t11 = OS[16*s+14+shuffle[16*s+14]];
      
      gamma128[2] = _mm_set_ps(gamma[t11+n*4],
                               gamma[t10+n*4], 
                               gamma[t9+n*4], 
                               gamma[t8+n*4]);
      
      t12 = OS[16*s+2+shuffle[16*s+3]];
      t13 = OS[16*s+6+shuffle[16*s+7]];
      t14 = OS[16*s+10+shuffle[16*s+11]];
      t15 = OS[16*s+14+shuffle[16*s+15]];

      gamma128[3] = _mm_set_ps(gamma[t15+n*4], 
                               gamma[t14+n*4],
                               gamma[t13+n*4], 
                               gamma[t12+n*4]);

      beta_old128[0] = _mm_load_ps(&beta_old[4*s]);
      beta_old128[1] = _mm_load_ps(&beta_old[4*s+num_points/2]);
      
      
      m[0] = _mm_add_ps(beta_old128[0], gamma128[0]);
      m[1] = _mm_add_ps(beta_old128[1], gamma128[1]);
      m[2] = _mm_add_ps(beta_old128[0], gamma128[2]);
      m[3] = _mm_add_ps(beta_old128[1], gamma128[3]);

      beta_new128[0] = _mm_max_ps(m[0], m[1]); 
      beta_new128[1] = _mm_max_ps(m[2], m[3]); 

      beta_tmp1 = _mm_unpacklo_ps(beta_new128[0], beta_new128[1]);
      beta_tmp2 = _mm_unpackhi_ps(beta_new128[0], beta_new128[1]);

      _mm_store_ps(&beta_new[8*s], beta_tmp1);
      _mm_store_ps(&beta_new[8*s+4], beta_tmp2);      

      alpha128_0 = (__m128*) &alpha[n*num_points+8*s];
      alpha128_1 = (__m128*) &alpha[n*num_points+8*s+4];

      *sum_tmp0 = _mm_unpacklo_ps(m[0], m[2]) + *alpha128_0;
      *sum_tmp1 = _mm_unpackhi_ps(m[0], m[2]) + *alpha128_1;
      *sum_tmp2 = _mm_unpacklo_ps(m[1], m[3]) + *alpha128_0;
      *sum_tmp3 = _mm_unpackhi_ps(m[1], m[3]) + *alpha128_1;

      for(c = 0; c < n_cbit; c++)
      {
        b0 = (t0>>(n_cbit-c-1))&1;
        b1 = (t8>>(n_cbit-c-1))&1;
        b2 = (t1>>(n_cbit-c-1))&1;
        b3 = (t9>>(n_cbit-c-1))&1;
        b4 = (t2>>(n_cbit-c-1))&1;
        b5 = (t10>>(n_cbit-c-1))&1;
        b6 = (t3>>(n_cbit-c-1))&1;
        b7 = (t11>>(n_cbit-c-1))&1;
        b8 = (t4>>(n_cbit-c-1))&1;
        b9 = (t12>>(n_cbit-c-1))&1;
        b10 = (t5>>(n_cbit-c-1))&1;
        b11 = (t13>>(n_cbit-c-1))&1;
        b12 = (t6>>(n_cbit-c-1))&1;
        b13 = (t14>>(n_cbit-c-1))&1;
        b14 = (t7>>(n_cbit-c-1))&1;
        b15 = (t15>>(n_cbit-c-1))&1;
   //   

        metrics[2*c+b0] = (metrics[2*c+b0] - sum[0] > 0) ? metrics[2*c+b0] : sum[0];
        metrics[2*c+b1] = (metrics[2*c+b1] - sum[1] > 0) ? metrics[2*c+b1] : sum[1];
        metrics[2*c+b2] = (metrics[2*c+b2] - sum[2] > 0) ? metrics[2*c+b2] : sum[2];
        metrics[2*c+b3] = (metrics[2*c+b3] - sum[3] > 0) ? metrics[2*c+b3] : sum[3];
        metrics[2*c+b4] = (metrics[2*c+b4] - sum[4] > 0) ? metrics[2*c+b4] : sum[4];
        metrics[2*c+b5] = (metrics[2*c+b5] - sum[5] > 0) ? metrics[2*c+b5] : sum[5];
        metrics[2*c+b6] = (metrics[2*c+b6] - sum[6] > 0) ? metrics[2*c+b6] : sum[6];
        metrics[2*c+b7] = (metrics[2*c+b7] - sum[7] > 0) ? metrics[2*c+b7] : sum[7];
        metrics[2*c+b8] = (metrics[2*c+b8] - sum[8] > 0) ? metrics[2*c+b8] : sum[8];
        metrics[2*c+b9] = (metrics[2*c+b9] - sum[9] > 0) ? metrics[2*c+b9] : sum[9];
        metrics[2*c+b10] = (metrics[2*c+b10] - sum[10] > 0) ? metrics[2*c+b10] : sum[10];
        metrics[2*c+b11] = (metrics[2*c+b11] - sum[11] > 0) ? metrics[2*c+b11] : sum[11];
        metrics[2*c+b12] = (metrics[2*c+b12] - sum[12] > 0) ? metrics[2*c+b12] : sum[12];
        metrics[2*c+b13] = (metrics[2*c+b13] - sum[13] > 0) ? metrics[2*c+b13] : sum[13];
        metrics[2*c+b14] = (metrics[2*c+b14] - sum[14] > 0) ? metrics[2*c+b14] : sum[14];
        metrics[2*c+b15] = (metrics[2*c+b15] - sum[15] > 0) ? metrics[2*c+b15] : sum[15];

      }
    }

    for(c = 0; c < n_cbit; c++)
      llr[n*n_cbit+c] = metrics[2*c+1]-metrics[2*c];      

    beta_tmp = beta_old;
    beta_old = beta_new;
    beta_new = beta_tmp;
  }
  volk_fec_free(metrics);
}

#endif /* LV_HAVE SSE */
#endif /* INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_a_H */