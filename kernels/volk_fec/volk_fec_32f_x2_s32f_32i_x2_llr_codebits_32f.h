#ifndef INCLUDED_VOLD_FEC_32f_s32f_32i_llr_codebits_32f_u_H
#define INCLUDED_VOLD_FEC_32f_s32f_32i_llr_codebits_32f_u_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "volk_fec_32f_x2_s32f_32i_viterbi_metric_32i.h"

#ifdef LV_HAVE_GENERIC

static inline void volk_fec_32f_x2_s32f_32i_x2_llr_codebits_32f_generic(float *beta,
                                              float llr, float *gamma, 
                                              const int length, const int *OS, 
                                              const int *shuffle, unsigned int num_points)
{

}

#endif
#endif