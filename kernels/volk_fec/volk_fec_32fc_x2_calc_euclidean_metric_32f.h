#ifndef INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_u_H
#define INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_u_H

#include <inttypes.h>
#include <stdio.h>
#include <volk_fec/volk_fec_complex.h>

#ifdef LV_HAVE_GENERIC
/*! 
  \brief calculates euclidean metric for two input vectors and stores their
         their result a the third vector
  \param O FSM Output cardinality
  \param a_vector Input samples
  \param b_vector Symbol table
  \param c_vector Branch metrics are stored in this vector
*/

static inline void volk_fec_32fc_x2_calc_euclidean_metric_32f_generic(float *c_vector,
                                                    const lv_32fc_t *a_vector,
                                                    const lv_32fc_t *b_vector,
                                                    const int O,
                                                    unsigned int num_points)
{

  unsigned int number = 0;
  for(number = 0; number < num_points; number++)
  {
    int o = 0;
    for(o = 0; o < O; o++)
    {
      c_vector[number+o*num_points] = 0.0;

      lv_32fc_t s = a_vector[number] - b_vector[o];
      c_vector[number+o*num_points] = lv_creal(s)*lv_creal(s)+lv_cimag(s)*lv_cimag(s); 
    }
  }
}

#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_u_H */

#ifndef INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_a_H
#define INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_a_H

#include <inttypes.h>
#include <stdio.h>
#include <volk_fec/volk_fec_complex.h>

#ifdef LV_HAVE_SSE4_1

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>

static inline void volk_fec_32fc_x2_calc_euclidean_metric_32f_a_sse4(float *c_vector,
                                                    const lv_32fc_t *a_vector,
                                                    const lv_32fc_t *b_vector,
                                                    const int O,
                                                    unsigned int num_points)
{
  const unsigned int quarter_points = num_points / 4;
  unsigned int number = 0;
  const lv_32fc_t *input_ptr =  a_vector;
  float *metrics_ptr = c_vector;
  __m128 input_val1, input_val2, table_val, distance1, distance2, mul_val1, mul_val2;
  __m128 metrics_val;

  int o = 0;

  for(o = 0; o < O; o++)
  {
    input_ptr = (lv_32fc_t*) a_vector;
    table_val = _mm_set_ps(lv_cimag(b_vector[o]),
                           lv_creal(b_vector[o]),
                           lv_cimag(b_vector[o]),
                           lv_creal(b_vector[o]));

    for(number = 0; number < quarter_points; number++)
    {
      input_val1 = _mm_load_ps((float*)input_ptr); // load a_1r, a_1i, a_2r, a_2i
      //input_val1 = _mm_shuffle_ps(input_val1, input_val1, 0xD8); //re-arrange to a_1r, a_2r, a_1i, a_2i
      input_ptr += 2;     
      input_val2 = _mm_load_ps((float*) input_ptr);
      //input_val2 = _mm_shuffle_ps(input_val2, input_val2, 0xD8);
      input_ptr += 2;

      distance1 = _mm_sub_ps(input_val1, table_val);
      distance2 = _mm_sub_ps(input_val2, table_val);

      mul_val1 = _mm_mul_ps(distance1, distance1);
      mul_val2 = _mm_mul_ps(distance2, distance2);

      metrics_val = _mm_hadd_ps(mul_val1, mul_val2);
      _mm_store_ps(metrics_ptr, metrics_val);
      metrics_ptr += 4;
    }

    // for(number = quarter_points*4; number < num_points; number++)
    // {
    //   c_vector[number+o*num_points] = 0.0;

    //   lv_32fc_t s = a_vector[number] - b_vector[o];
    //   c_vector[number+o*num_points] += lv_creal(s)*lv_creal(s)+lv_cimag(s)*lv_cimag(s); 
    //   metrics_ptr++;
    //   printf("In der Schleife[%d]\n", number);
    //   input_ptr++;       
    // }
  }
} 

#endif /* LV_HAVE_SSE4_1 */

#ifdef LV_HAVE_AVX
#include <immintrin.h>

static inline void volk_fec_32fc_x2_calc_euclidean_metric_32f_a_avx(float *c_vector,
                                                    const lv_32fc_t *a_vector,
                                                    const lv_32fc_t *b_vector,
                                                    const int O,
                                                    unsigned int num_points)
{
  const unsigned int eigth_points = num_points / 8;
  unsigned int number = 0;
  const lv_32fc_t *input_ptr =  a_vector;
  float *metrics_ptr = c_vector;
  __m256 input_val1, input_val2, table_val, distance1, distance2, mul_val1, mul_val2;
  __m256 metrics_val, metrics_tmp;

  int o = 0;

  for(o = 0; o < O; o++)
  {
    input_ptr = a_vector;
    for(number = 0; number < eigth_points; number ++)
    {
      input_val1 = _mm256_load_ps((float*)input_ptr); // load a_1r, a_1i, a_2r, a_2i
      //input_val1 = _mm_shuffle_ps(input_val1, input_val1, 0xD8); //re-arrange to a_1r, a_2r, a_1i, a_2i
      input_ptr += 4;     
      input_val2 = _mm256_load_ps((float*) input_ptr);
      //input_val2 = _mm_shuffle_ps(input_val2, input_val2, 0xD8);
      input_ptr += 4;
      metrics_tmp = _mm256_setzero_ps();

     table_val = _mm256_set_ps(lv_cimag(b_vector[o]),
                               lv_creal(b_vector[o]),
                               lv_cimag(b_vector[o]),
                               lv_creal(b_vector[o]),
                               lv_cimag(b_vector[o]),
                               lv_creal(b_vector[o]),
                               lv_cimag(b_vector[o]),
                               lv_creal(b_vector[o]));

      distance1 = _mm256_sub_ps(input_val1, table_val);
      distance2 = _mm256_sub_ps(input_val2, table_val);

      distance1 = _mm256_mul_ps(distance1, distance1);
      distance2 = _mm256_mul_ps(distance2, distance2);

      mul_val1 = _mm256_permute2f128_ps(distance1, distance2, 0x20);
      mul_val2 = _mm256_permute2f128_ps(distance1, distance2, 0x31), 


      metrics_val = _mm256_hadd_ps(mul_val1, mul_val2);
      //metrics_val = __mm256_shuffle_ps(metrics_val, metrics_val, )
      metrics_val = _mm256_add_ps(metrics_val, metrics_tmp);
      _mm256_store_ps(metrics_ptr, metrics_val);
      metrics_tmp = _mm256_load_ps(metrics_ptr);
      metrics_ptr += 8;
    }
    
    for(number = eigth_points*8; number < num_points; number++)
    {
      c_vector[number+o*num_points] = 0.0;

      lv_32fc_t s = a_vector[number] - b_vector[o];
      c_vector[number+o*num_points] += lv_creal(s)*lv_creal(s)+lv_cimag(s)*lv_cimag(s); 

      metrics_ptr++;
      input_ptr++;  
    } 
  }
}

#endif /* LV_HAVE_AVX */
#endif /* INCLUDED_volk_fec_32fc_x2_calc_euclidean_metric_32f_a_H */