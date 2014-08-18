#ifndef INCLUDED_VOLD_FEC_32f_forward_recursionpuppet_32f_u_H
#define INCLUDED_VOLD_FEC_32f_forward_recursionpuppet_32f_u_H


#include <stdlib.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC
#include <volk_fec/volk_fec.h>
#include <volk_fec/volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f.h>

static inline void volk_fec_32f_forward_recursionpuppet_32f_generic(float *alpha,
                                              float *gamma, unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  int shuffle[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};

  volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f_manual(alpha,
                                              gamma, num_points/18, 
                                              &OS[0], &shuffle[0],
                                              16, "generic");

}

#endif /*LV_HAVE_GENERIC*/
#endif /*INCLUDED_VOLD_FEC_32f_s32f_32i_forward_recursion_32f_u_H*/

#ifndef INCLUDED_VOLD_FEC_32f_forward_recursionpuppet_32f_32f_a_H
#define INCLUDED_VOLD_FEC_32f_forward_recursionpuppet_32f_32f_a_H

#include <volk_fec/volk_fec.h>
#include <volk_fec/volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f.h>

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"
static inline void volk_fec_32f_forward_recursionpuppet_32f_a_sse4(float *alpha,
                                              float *gamma, unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  int shuffle[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f_manual(alpha,
                                                    gamma, num_points/18, 
                                                    &OS[0], &shuffle[0],
                                                    16, "a_sse4");
}
#endif /*LV_HAVE_SSE4_1*/
#endif /* INCLUDED_VOLD_FEC_32f_s32f_32i_calc_branch_metric_32f_a_H */