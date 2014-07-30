#ifndef INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H
#define INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H


#include <stdlib.h>
#include <stdio.h>
#include "volk_fec_32f_x2_s32f_32i_viterbi_metric_32i.h"

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
      t[0] = OS[4*s+shuffle[4*s]];
      t[1] = OS[4*s+2+shuffle[4*s+2]];
      t[2] = OS[4*s+shuffle[4*s+1]];
      t[3] = OS[4*s+2+shuffle[4*s+3]];

      m0 = alpha[n*num_points+2*s  ] + gamma[n*4+t[0]];
      m1 = alpha[n*num_points+2*s+1] + gamma[n*4+t[1]];
      m2 = alpha[n*num_points+2*s  ] + gamma[n*4+t[2]];
      m3 = alpha[n*num_points+2*s+1] + gamma[n*4+t[3]];
      if(n==0 && s == 0)
      {
        printf("a[0] %f g[] %f t[0] %d \n", alpha[n*num_points+2*s],  gamma[n*4+t[0]], t[0]);
        printf("a[1] %f g[] %f t[1] %d\n", alpha[n*num_points+2*s+1],  gamma[n*4+t[1]],t[1]);
        printf("a[0] %f g[] %f t[2]%d\n", alpha[n*num_points+2*s],  gamma[n*4+t[2]], t[2]);
        printf("a[1] %f g[] %f t[3] %d\n", alpha[n*num_points+2*s+1],  gamma[n*4+t[3]], t[3]);
      }


      alpha[(n+1)*num_points+s] = ((m1-m2>0) ? m1 : m0);
      alpha[(n+1)*num_points+s+num_points/2] = ((m3-m1>0) ? m3 : m2); 
    }
  } 
}

#endif /*LV_HAVE_GENERIC*/
#endif /*INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H*/