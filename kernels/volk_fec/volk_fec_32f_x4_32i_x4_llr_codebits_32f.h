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
  int n = 0;
  int s = 0;
  int c = 0;

  int O =  1<<n_cbit;
  float m0,m1,m2,m3, sum_tmp0, sum_tmp1, sum_tmp2, sum_tmp3;
  int t0,t1,t2,t3, cbit0, cbit1, cbit2, cbit3;

  float *metrics;
  metrics = (float*) malloc(sizeof(float)*4);

  float *beta_new = &beta[num_points];
  float *beta_old = &beta[0];
  float *beta_tmp = NULL;
 
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
      t1 = OS[4*s+2+shuffle[4*s+1]];
      t2 = OS[4*s+shuffle[4*s+2]];
      t3 = OS[4*s+2+shuffle[4*s+3]];



      // Calc state metrics
      m0 = beta_old[s] + gamma[n*O+t0];
      m1 = beta_old[s+num_points/2] + gamma[n*O+t1];
      m2 = beta_old[s] + gamma[n*O+t2];
      m3 = beta_old[s+num_points/2] + gamma[n*O+t3];

      // Chose metric with highest correlation/probability
      beta_new[2*s] = ((m1-m0>0) ? m1 : m0);
      beta_new[2*s+1] = ((m3-m2>0) ? m3 : m2);

      if (n == 4)
      {
        printf("beta[%d]: %f beta[%d]: %f\n", s,beta_old[s],s+num_points/2, beta_old[s+num_points/2]);
        printf("alpha[%d] %f alpha[%d] %f\n", 2*s, alpha[n*num_points+2*s], 2*s+1, alpha[n*num_points+2*s+1]);
        printf("gamma[%d] %f gamma[%d] %f gamma[%d] %f gamma[%d] %f\n\n", t0, gamma[n*O+t0], t1, gamma[n*O+t1], t2, gamma[n*O+t2], t3, gamma[n*O+t3]);
      }
     //printf("beta[%d]: %f\nbeta[%d]: %f\n", 2*s,beta_new[2*s],2*s+1, beta_new[2*s+1]);

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
        metrics[2*c+cbit1] = (metrics[2*c+cbit1] - sum_tmp1 > 0) ? metrics[2*c+cbit1] : sum_tmp1;;
        metrics[2*c+cbit2] = (metrics[2*c+cbit2] - sum_tmp2 > 0) ? metrics[2*c+cbit2] : sum_tmp2;;
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