#ifndef INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_u_H
#define INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_u_H

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

static inline void volk_fec_16i_x2_calc_euclidean_metric_16i_generic(short *c_vector,
                                                    const short *a_vector,
                                                    const short *b_vector,
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
#endif /* INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_u_H */

#ifndef INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_a_H
#define INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_a_H

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

static inline void volk_fec_16i_x2_calc_euclidean_metric_16i_a_sse4(short *c_vector,
                                                    const short *a_vector,
                                                    const short *b_vector,
                                                    const int O,
                                                    const int D,
                                                    unsigned int num_points)
{
  const unsigned int eight_points = num_points / 8;
  unsigned int number = 0;
  short *input_ptr = (short*) a_vector;
  short *metrics_ptr = (short*)c_vector;
  __m128i input_val, distance, table_val;
  __m128i metrics_val; 

  int o = 0;
  int d = 0;

  for (o = 0; o < O; o++)
  {
    input_ptr = (short*) a_vector;
    for(number = 0; number < eight_points; number++)
    {
      for(d = 0; d < D; d++)
      {    
        input_val  = _mm_loadu_si128 ((__m128i*)input_ptr);
        //fill table registers
        table_val = _mm_set_epi16(b_vector[o*D+d], b_vector[o*D+d],
                                  b_vector[o*D+d], b_vector[o*D+d],
                                  b_vector[o*D+d], b_vector[o*D+d],
                                  b_vector[o*D+d], b_vector[o*D+d]);

        // Calculate Distance and split it into 2 32 bit values
        distance    = _mm_sub_epi16(input_val, table_val);

        metrics_val = _mm_mullo_epi16( distance, distance);

        //store results and move pointer
        _mm_storeu_si128((__m128i*)metrics_ptr, metrics_val);
        metrics_ptr += 8;        
        input_ptr += 8; 
      }
    }
    // for (number = eight_points*8; number < num_points; number++)
    // {
    //  printf("number: %d   number_points: %d\n", number, num_points);
    //  c_vector[o*num_points+number] = 0;
    //  for(d = 0; d < D; d++)
    //  {  
    //    int s = a_vector[number] - b_vector[o*D+d];
    //    c_vector[o*num_points+number] += s*s;
    //  }  
    //  metrics_ptr++;
    //  input_ptr++;
    // }

  }
}


#endif /* LV_HAVE_SSE4_1 */
#endif /* INCLUDED_volk_fec_16i_x2_calc_euclidean_metric_16i_a_H */
