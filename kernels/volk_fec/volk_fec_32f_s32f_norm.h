#ifndef INCLUDED_volk_fec_32f_s32f_norm_32f_u_H
#define INCLUDED_volk_fec_32f_s32f_norm_32f_u_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC
/*! 
  \brief normalize vector (so it doesn't explode ;))
  \param vec_buffer Input
  \param norm Normalizationparameter
*/

static inline void volk_fec_32f_s32f_norm_generic(float *vec_buffer, float norm, unsigned int num_points)
{
	unsigned int number = 0;
	for(number = 0; number < num_points; number++)
	{
		vec_buffer[number] -= norm;
	}
}

#endif /* LV_HAVE_GENERIC */

#endif /* INCLUDED_volk_fec_32f_s32f_norm_32f_u_H */

#ifndef INCLUDED_volk_fec_32f_s32f_norm_32f_a_H
#define INCLUDED_volk_fec_32f_s32f_norm_32f_a_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_SSE
#include "xmmintrin.h"

static inline void volk_fec_32f_s32f_norm_a_sse(float *vec_buffer, float norm, unsigned int num_points)
{
	const unsigned int quarter_points = num_points/4;
	unsigned int number = 0;
	float *in_ptr = vec_buffer;
	__m128 buffer, norm_val;

	norm_val = _mm_set_ps1(norm);
	for(number = 0; number < quarter_points; number++)
	{
		buffer = _mm_load_ps(in_ptr);
		buffer = _mm_sub_ps(buffer, norm_val);
		_mm_store_ps(in_ptr, buffer);
		in_ptr += 4;
	}

	for(number = quarter_points*4; number < num_points; number++)
	{
		vec_buffer[number] -= norm;
	}
}

#endif /* LV_HAVE_SSE */
#endif /* INCLUDED_volk_fec_32f_s32f_norm_32f_a_H */

