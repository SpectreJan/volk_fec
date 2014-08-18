#ifndef INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_u_H
#define INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_u_H

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

static inline void volk_fec_32i_x2_calc_euclidean_metric_32f_generic(float *c_vector,
                                                    const int *a_vector,
                                                    const int *b_vector,
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
        int s = a_vector[number] - b_vector[o*D+d];
        c_vector[number+o*num_points] += s*s; 
      }
    }
  }
}
#endif /* LV_HAVE_GENERIC */
#endif /* INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_u_H */

#ifndef INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_a_H
#define INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_a_H

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

static inline void volk_fec_32i_x2_calc_euclidean_metric_32f_a_sse4(float *c_vector,
                                                    const int *a_vector,
                                                    const int *b_vector,
                                                    const int O,
                                                    const int D,
                                                    unsigned int num_points)
{
  const unsigned int quarter_points = num_points / 4;
  unsigned int number = 0;
  int *input_ptr = (int*) a_vector;
  float *metrics_ptr = c_vector;
  __m128i input_val, distance, table_val, metrics_int;
  __m128 metrics_val, metrics_tmp; 

  int o = 0;
  int d = 0;

  for (o = 0; o < O; o++)
  {
    input_ptr = (int*) a_vector;
    for(number = 0; number < quarter_points; number++)
    {
      input_val  = _mm_load_si128 ((__m128i*)input_ptr);
      metrics_tmp = _mm_setzero_ps();
      for(d = 0; d < D; d++)
      {    
        //fill table registers
        table_val = _mm_set_epi32(b_vector[o*D+d], b_vector[o*D+d],
                                  b_vector[o*D+d], b_vector[o*D+d]);

        // Calculate Distance and split it into 2 32 bit values
        distance    = _mm_sub_epi32(input_val, table_val);

        metrics_int = _mm_mullo_epi32( distance, distance); 

        // convert int to float
        metrics_val = _mm_cvtepi32_ps(metrics_int);
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
      c_vector[o*num_points+number] = 0;
      for(d = 0; d < D; d++)
      {
        int s = a_vector[number] - b_vector[o*D+d];
        c_vector[o*num_points+number] += s*s;
      }
      metrics_ptr++;
      input_ptr++;
    } 
  }
}


#endif /* LV_HAVE_SSE4_1 */
#endif /* INCLUDED_volk_fec_32i_x2_calc_euclidean_metric_32f_a_H */
