#ifndef INCLUDED_VOLD_FEC_32f_x3_llr_codebitspuppet_32f_u_H
#define INCLUDED_VOLD_FEC_32f_x3_llr_codebitspuppet_32f_u_H

#define INF 1e9
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef LV_HAVE_GENERIC

static inline void volk_fec_32f_x3_llr_codebitspuppet_32f_generic(float *alpha, float *gamma,
                                                                float *beta, float *llr,
                                                                unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  int shuffle[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  volk_fec_32f_x4_32i_x4_llr_codebits_32f_manual(alpha, gamma, beta, llr,
                                                  2, num_points/9, 
                                                  &OS[0], &shuffle[0],
                                                  16, "generic"); 
}

#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_VOLD_FEC_32f_x3_32i_x4_llr_codebitspuppet_32f_u_H */

#ifndef INCLUDED_VOLD_FEC_32f_x3_llr_codebitspuppet_32f_a_H
#define INCLUDED_VOLD_FEC_32f_x3_llr_codebitspuppet_32f_a_H
#define INF 1e9
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <volk_fec/volk_fec.h>

#ifdef LV_HAVE_SSE4_1
#include "smmintrin.h"
static inline void volk_fec_32f_x3_llr_codebitspuppet_32f_a_sse4(float *alpha, float *gamma,
                                                                float *beta, float *llr,
                                                                unsigned int num_points)
{
  int OS[] = {0, 3, 3, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 3, 3, 0, 2, 1, 1, 2, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 1, 2};
  int shuffle[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  volk_fec_32f_x4_32i_x4_llr_codebits_32f_manual(alpha, gamma, beta, llr,
                                                  2, num_points/9, 
                                                  &OS[0], &shuffle[0],
                                                  16, "a_sse4"); 
}

#endif /* LV_HAVE SSE */
#endif /* INCLUDED_VOLD_FEC_32f_x4_32i_x4_llr_codebits_32f_a_H */