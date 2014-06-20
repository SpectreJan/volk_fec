#include "qa_utils.h"
#include <volk_fec/volk_fec.h>
#include <boost/test/unit_test.hpp>
VOLK_RUN_TESTS(volk_fec_32f_s32f_normalize, 1e-4, 100, 20462, 1);
VOLK_RUN_TESTS(volk_fec_32f_s32f_norm, 1e-4, 3.0, 20462, 1);

//VOLK_RUN_TESTS(volk_fec_16i_x5_add_quad_16i_x4, 1e-4, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_16i_branch_4_state_8, 1e-4, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_32f_s32f_calc_metricpuppet_32f, 1e-4, 4.0, 20462, 1);
//VOLK_RUN_TESTS(volk_fec_16i_max_star_16i, 0, 0, 20462, 10000);
//VOLK_RUN_TESTS(volk_fec_16i_max_star_horizontal_16i, 0, 0, 20462, 10000);
//VOLK_RUN_TESTS(volk_fec_16i_permute_and_scalar_add, 1e-4, 0, 2046, 1000);
//VOLK_RUN_TESTS(volk_fec_16i_x4_quad_max_star_16i, 1e-4, 0, 2046, 1000);
//VOLK_RUN_TESTS(volk_fec_16i_32fc_dot_prod_32fc, 1e-4, 0, 204602, 1);
//VOLK_RUN_TESTS(volk_fec_32fc_x2_conjugate_dot_prod_32fc, 1e-4, 0, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_32fc_s32f_x2_power_spectral_density_32f, 1e-4, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_32f_s32f_32f_fm_detect_32f, 1e-4, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_32u_popcnt, 0, 0, 2046, 10000);
//VOLK_RUN_TESTS(volk_fec_64u_popcnt, 0, 0, 2046, 10000);

