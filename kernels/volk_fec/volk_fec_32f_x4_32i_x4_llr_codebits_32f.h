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

     //printf("beta[%d]: %f\nbeta[%d]: %f\n", 2*s,beta_new[2*s],2*s+1, beta_new[2*s+1]);
//      if(n == 0)
//        printf("Beta_new[%d] %f Beta_new[%d] %f\n",2*s, beta_new[2*s], 2*s+1, beta_new[2*s+1]);

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
        
        // 
        metrics[2*c+cbit0] = (metrics[2*c+cbit0] - sum_tmp0 > 0) ? metrics[2*c+cbit0] : sum_tmp0;
        metrics[2*c+cbit1] = (metrics[2*c+cbit1] - sum_tmp1 > 0) ? metrics[2*c+cbit1] : sum_tmp1;
        metrics[2*c+cbit2] = (metrics[2*c+cbit2] - sum_tmp2 > 0) ? metrics[2*c+cbit2] : sum_tmp2;
        metrics[2*c+cbit3] = (metrics[2*c+cbit3] - sum_tmp3 > 0) ? metrics[2*c+cbit3] : sum_tmp3;
       //printf("metrics[%d] %f metrics[%d] %f\n", 2*c, metrics[2*c], 2*c+1, metrics[2*c+1]);
        //
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
  __m128 m[4];
  __m128 beta_new128[2], beta_old128[2], sum_tmp[4], gamma128[4];
  __m128 beta_tmp1, beta_tmp2;
  __m128 *metrics128;
  float *metrics;
  metrics = (float*) volk_fec_malloc(sizeof(float)*8*(1<<n_cbit), volk_fec_get_alignment());

  metrics128 = (__m128*) metrics;

  float *beta_new = &beta[num_points];
  float *beta_old = &beta[0];
  float *beta_tmp = NULL;

  int n = 0;
  int s = 0;
  int c = 0;
  int i = 0;
  int cb;
  int quarter_points = num_points/4;

  for(n = 0; n < (length/4); n++)
  {
    for(i = 0; i < 4; i++)
    {
      cb = length-n*4-i-1;
/*      for(c = 0; c < n_cbit*2; c++)
      {
        metrics[c] = -INF;
      }*/
      for(s = 0; s < quarter_points/2; s++)
      {
        gamma128[0] = _mm_set_ps(gamma[OS[16*s+12+shuffle[16*s+12]]+cb*4], 
                                 gamma[OS[16*s+8+shuffle[16*s+8]]+cb*4],
                                 gamma[OS[16*s+4+shuffle[16*s+4]]+cb*4], 
                                 gamma[OS[16*s+shuffle[16*s]]+cb*4]);

        gamma128[1] = _mm_set_ps(gamma[OS[16*s+12+shuffle[16*s+13]]+cb*4], 
                                 gamma[OS[16*s+8+shuffle[16*s+9]]+cb*4],
                                 gamma[OS[16*s+4+shuffle[16*s+5]]+cb*4], 
                                 gamma[OS[16*s+shuffle[16*s+1]]+cb*4]);
        
        gamma128[2] = _mm_set_ps(gamma[OS[16*s+14+shuffle[16*s+14]]+cb*4],
                                 gamma[OS[16*s+10+shuffle[16*s+10]]+cb*4], 
                                 gamma[OS[16*s+6+shuffle[16*s+6]]+cb*4], 
                                 gamma[OS[16*s+2+shuffle[16*s+2]]+cb*4]);
        
        gamma128[3] = _mm_set_ps(gamma[OS[16*s+14+shuffle[16*s+15]]+cb*4], 
                                 gamma[OS[16*s+10+shuffle[16*s+11]]+cb*4],
                                 gamma[OS[16*s+6+shuffle[16*s+7]]+cb*4], 
                                 gamma[OS[16*s+2+shuffle[16*s+3]]+cb*4]);

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


/*        if(cb == 0)
        {
          printf("Beta_new[%d] %f Beta_new[%d] %f\nBeta_new[%d] %f Beta_new[%d] %f\n",s*8, beta_new[s*8], s*8+1, beta_new[s*8+1], s*8+2, beta_new[s*8+2], s*8+3, beta_new[s*8+3]);
          printf("\n");
          printf("Beta_new[%d] %f Beta_new[%d] %f\nBeta_new[%d] %f Beta_new[%d] %f\n",s*8+4, beta_new[s*8+4], s*8+5, beta_new[s*8+5], s*8+6, beta_new[s*8+6], s*8+7, beta_new[s*8+7]);
        
          printf("G %f G %f G %f G %f\n", gamma[OS[16*s+12+shuffle[16*s+12]]+cb*4], 
                                 gamma[OS[16*s+8+shuffle[16*s+8]]+cb*4],
                                 gamma[OS[16*s+4+shuffle[16*s+4]]+cb*4], 
                                 gamma[OS[16*s+shuffle[16*s]]+cb*4]);
          printf("G %f G %f G %f G %f\n", gamma[OS[16*s+12+shuffle[16*s+13]]+cb*4], 
                     gamma[OS[16*s+8+shuffle[16*s+9]]+cb*4],
                     gamma[OS[16*s+4+shuffle[16*s+5]]+cb*4], 
                     gamma[OS[16*s+shuffle[16*s+1]]+cb*4]);
          printf("G %f G %f G %f G %f\n", gamma[OS[16*s+14+shuffle[16*s+14]]+cb*4], 
                                 gamma[OS[16*s+10+shuffle[16*s+10]]+cb*4],
                                 gamma[OS[16*s+6+shuffle[16*s+6]]+cb*4], 
                                 gamma[OS[16*s+2+shuffle[16*s+2]]+cb*4]);
          printf("G %f G %f G %f G %f\n", gamma[OS[16*s+14+shuffle[16*s+15]]+cb*4], 
                     gamma[OS[16*s+10+shuffle[16*s+11]]+cb*4],
                     gamma[OS[16*s+6+shuffle[16*s+7]]+cb*4], 
                     gamma[OS[16*s+2+shuffle[16*s+3]]+cb*4]);


        }*/
       // sum_tmp[0] = m[0] + alpha_old[];
       // sum_tmp[1] = m[1] + alpha_old[];
       // sum_tmp[2] = m[2] + alpha_old[];
       // sum_tmp[3] = m[3] + alpha_old[];

/*        for(c = 0; c < n_cbit; c++)
        {
          cbit0 = (t0>>(n_cbit-c-1))&1;
          cbit1 = (t1>>(n_cbit-c-1))&1;
          cbit2 = (t2>>(n_cbit-c-1))&1;
          cbit3 = (t3>>(n_cbit-c-1))&1;
          
          // 
          metrics[2*c+cbit0] = (metrics[2*c+cbit0] - sum_tmp0 > 0) ? metrics[2*c+cbit0] : sum_tmp0;
          metrics[2*c+cbit1] = (metrics[2*c+cbit1] - sum_tmp1 > 0) ? metrics[2*c+cbit1] : sum_tmp1;
          metrics[2*c+cbit2] = (metrics[2*c+cbit2] - sum_tmp2 > 0) ? metrics[2*c+cbit2] : sum_tmp2;
          metrics[2*c+cbit3] = (metrics[2*c+cbit3] - sum_tmp3 > 0) ? metrics[2*c+cbit3] : sum_tmp3;
         //printf("metrics[%d] %f metrics[%d] %f\n", 2*c, metrics[2*c], 2*c+1, metrics[2*c+1]);
          //
        }*/

      }
/*      for(c = 0; c < n_cbit/4; c++)
      {
        llr[n*n_cbit+c] = metrics[2*c+1]-metrics[2*c];
      }
      for(c = c*4; c < n_cbit; c++)
      {
        llr[n*n_cbit+c] = metrics[2*c+1]-metrics[2*c];
      }*/
      beta_tmp = beta_old;
      beta_old = beta_new;
      beta_new = beta_tmp;
    }
  }
  volk_fec_free(metrics);
}

#endif /* LV_HAVE SSE */
#endif /* INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_a_H */