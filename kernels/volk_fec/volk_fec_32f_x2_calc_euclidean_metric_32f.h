#ifndef INCLUDED_volk_fec_32f_x2_calc_euclidean_metric_32f_u_H
#define INCLUDED_volk_fec_32f_x2_calc_euclidean_metric_32f_u_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_GENERIC
/*! 
  \brief calculates euclidean metric for two input vectors and stores their
         their result a the third vector
  \param O FSM Output cardinality
  \param D Dimensionality of input vector
  \param table Vector with valid symbols
  \param input Input vector
  \param metric Metrics are stored in this vector
*/

static inline void volk_fec_32f_x2_calc_euclidean_metric_32f_generic(float *c_vector,
                                                    const float *a_vector,
                                                    const float *b_vector,
                                                    const int O,
                                                    const int D,
                                                    unsigned int num_points)
{

  unsigned int number = 0;
  for(number = 0; number < num_points; number++)
  {
    int o = 0;
    for(o = 0; o < O; o++)
    {
      int d = 0;
      c_vector[number+o*num_points] = 0.0;
      for(d = 0; d < D; d++)
      {
        float s = a_vector[number] - b_vector[o*D+d];
        c_vector[number+o*num_points] += s*s; 
      }
    }
  }
}
#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32f_x2_calc_euclidean_metric_32f_u_H */

#ifndef INCLUDED_volk_fec_32f_x2_calc_euclidean_metric_32f_a_H
#define INCLUDED_volk_fec_32f_x2_calc_euclidean_metric_32f_a_H

#include <inttypes.h>
#include <stdio.h>

#ifdef LV_HAVE_SSE4_1

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
/*! 
  \brief calculates euclidean metric for two input vectors and stores their
         their result a the third vector
  \param c_vector Output metrics vector 
  \param a_vector Input vector
  \param b_vector Symbol table
  \param O Symbol cardinality
  \param D Smybol dimensionality
*/

static inline void volk_fec_32f_x2_calc_euclidean_metric_32f_a_sse4(float *c_vector,
                                                    const float *a_vector,
                                                    const float *b_vector,
                                                    const int O,
                                                    const int D,
                                                    unsigned int num_points)
{
  const unsigned int quarter_points = num_points / 4;
  unsigned int number = 0;
  float *input_ptr = (float*) a_vector;
  float *metrics_ptr = c_vector;
  __m128 input_val, distance, table_val;
  __m128 metrics_val, metrics_tmp; 

  int o = 0;
  int d = 0;

  for (o = 0; o < O; o++)
  {
    input_ptr = (float*) a_vector;
    for(number = 0; number < quarter_points; number++)
    {
      input_val  = _mm_load_ps(input_ptr);
      metrics_tmp = _mm_setzero_ps();
      for(d = 0; d < D; d++)
      {    
        //fill table registers
        table_val = _mm_load1_ps(&b_vector[o*D+d]);

        // Calculate Distance and split it into 2 32 bit values
        distance    = _mm_sub_ps(input_val, table_val);

        metrics_val = _mm_mul_ps( distance, distance); 
        metrics_val = _mm_add_ps(metrics_val, metrics_tmp);
        //store results and move pointer 
        _mm_store_ps(metrics_ptr, metrics_val); 
        metrics_tmp = _mm_load_ps(metrics_ptr);
      }
      metrics_ptr += 4;        
      input_ptr += 4;
    }

    for (number = quarter_points*4; number < num_points; number++)
    {
      c_vector[o*num_points+number] = 0.0;
      for(d = 0; d < D; d++)
      {
        float s = a_vector[number] - b_vector[o*D+d];
        c_vector[o*num_points+number] += s*s;
      }
      metrics_ptr++;
      input_ptr++;
    } 
  }
}


#endif /* LV_HAVE_SSE4_1 */

#ifdef LV_HAVE_AVX
#include <immintrin.h>
/*! 
  \brief calculates euclidean metric for two input vectors and stores their
         their result a the third vector
  \param c_vector Output metrics vector 
  \param a_vector Input vector
  \param b_vector Symbol table
  \param O Symbol cardinality
  \param D Smybol dimensionality
*/
static inline void volk_fec_32f_x2_calc_euclidean_metric_32f_a_avx(float *c_vector,
                                                    const float *a_vector,
                                                    const float *b_vector,
                                                    const int O,
                                                    const int D,
                                                    unsigned int num_points)
{
  const unsigned int eigth_points = num_points / 8;
  unsigned int number = 0;
  float *input_ptr = (float*) a_vector;
  float *metrics_ptr = c_vector;
  __m256 input_val, distance, table_val;
  __m256 metrics_val, metrics_tmp; 

  int o = 0;
  int d = 0;

  for (o = 0; o < O; o++)
  {
    input_ptr = (float*) a_vector;
    for(number = 0; number < eigth_points; number++)
    {
      input_val  = _mm256_load_ps(input_ptr);
      metrics_tmp = _mm256_setzero_ps();
      for(d = 0; d < D; d++)
      {    
        //fill table registers
        table_val = _mm256_set1_ps(b_vector[o*D+d]);

        // Calculate Distance and split it into 2 32 bit values
        distance    = _mm256_sub_ps(input_val, table_val);

        metrics_val = _mm256_mul_ps( distance, distance); 
        metrics_val = _mm256_add_ps(metrics_val, metrics_tmp);
        //store results and move pointer 
        _mm256_store_ps(metrics_ptr, metrics_val);
        _mm256_load_ps(metrics_ptr); 
      }
      metrics_ptr += 8;        
      input_ptr += 8;   
    }

    for (number = eigth_points*8; number < num_points; number++)
    {
      c_vector[o*num_points+number] = 0.0;
      for(d = 0; d < D; d++)
      {
        float s = a_vector[number] - b_vector[o*D+d];
        c_vector[o*num_points+number] += s*s;
      }
      metrics_ptr++;
      input_ptr++;
    } 
  }
}
#endif /* LV_HAVE_AVX */
#endif /* INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_a_H */
