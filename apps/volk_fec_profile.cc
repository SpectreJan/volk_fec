#include "qa_utils.h"

#include <volk_fec/volk_fec.h>
#include <volk_fec/volk_fec_prefs.h>

#include <ciso646>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
    // Adding program options
    boost::program_options::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")
      ("benchmark,b",
            boost::program_options::value<bool>()->default_value( false )
                                                ->implicit_value( true ),
            "Run all kernels (benchmark mode)")
      ("tests-regex,R",
            boost::program_options::value<std::string>(),
            "Run tests matching regular expression.")
      ;

    // Handle the options that were given
    boost::program_options::variables_map vm;
    bool benchmark_mode;
    std::string kernel_regex;
    bool store_results = true;
    try {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
        benchmark_mode = vm.count("benchmark")?vm["benchmark"].as<bool>():false;
        if ( vm.count("tests-regex" ) ) {
            kernel_regex = vm["tests-regex"].as<std::string>();
            store_results = false;
            std::cout << "Warning: using a regexp will not save results to a config" << std::endl;
        }
        else {
            kernel_regex = ".*";
            store_results = true;
        }
    } catch (boost::program_options::error& error) {
        std::cerr << "Error: " << error.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return 1;
    }
    /** --help option
*/
    if ( vm.count("help") )
    {
      std::cout << "The VOLK profiler." << std::endl
                << desc << std::endl;
      return 0;
    }


    // Run tests
    std::vector<std::string> results;

    //VOLK_PROFILE(volk_fec_16i_x5_add_quad_16i_x4, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_16i_branch_4_state_8, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_16i_max_star_16i, 0, 0, 204602, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_16i_max_star_horizontal_16i, 0, 0, 204602, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_16i_permute_and_scalar_add, 1e-4, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_16i_x4_quad_max_star_16i, 1e-4, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32fc_x2_conjugate_dot_prod_32fc, 1e-4, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32fc_s32f_x2_power_spectral_density_32f, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32f_s32f_32f_fm_detect_32f, 1e-4, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32u_popcnt, 0, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_64u_popcnt, 0, 0, 2046, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32fc_s32fc_multiply_32fc, 1e-4, lv_32fc_t(1.0, 0.5), 204602, 1000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32f_s32f_norm, 1e-4, 100, 204602, 10000, &results, benchmark_mode, kernel_regex);
    //VOLK_PROFILE(volk_fec_32f_s32f_normalize, 1e-4, 100, 204602, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_16i_x2_calc_euclidean_metricpuppet_16i, volk_fec_16i_x2_calc_euclidean_metric_16i, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_16i_x2_calc_euclidean_metricpuppet_32f, volk_fec_16i_x2_calc_euclidean_metric_32f, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32i_x2_calc_euclidean_metricpuppet_32f, volk_fec_32i_x2_calc_euclidean_metric_32f, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32f_x2_calc_euclidean_metricpuppet_32f, volk_fec_32f_x2_calc_euclidean_metric_32f, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32fc_x2_calc_euclidean_metricpuppet_32f, volk_fec_32fc_x2_calc_euclidean_metric_32f, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32f_x2_viterbi_metricpuppet_32i, volk_fec_32f_x2_s32f_32i_viterbi_metric_32i, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32f_calc_branch_metricpuppet_32f, volk_fec_32f_s32f_32i_calc_branch_metric_32f, 1e-4, 0, 204800, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32f_forward_recursionpuppet_32f, volk_fec_32f_x2_s32f_32i_x2_forward_recursion_32f, 1e-4, 0, 204750, 10000, &results, benchmark_mode, kernel_regex);
    VOLK_PUPPET_PROFILE(volk_fec_32f_x3_llr_codebitspuppet_32f, volk_fec_32f_x4_32i_x4_llr_codebitspuppet_32f, 1e-4, 0, 102420, 10000, &results, benchmark_mode, kernel_regex);
    // Until we can update the config on a kernel by kernel basis
    // do not overwrite volk_fec_config when using a regex.
    if(store_results) {
        char path[1024];
        volk_fec_get_config_path(path);

        const fs::path config_path(path);

        if (not fs::exists(config_path.branch_path()))
        {
            std::cout << "Creating " << config_path.branch_path() << "..." << std::endl;
            fs::create_directories(config_path.branch_path());
        }

        std::cout << "Writing " << config_path << "..." << std::endl;
        std::ofstream config(config_path.string().c_str());
        if(!config.is_open()) { //either we don't have write access or we don't have the dir yet
            std::cout << "Error opening file " << config_path << std::endl;
        }

        config << "\
#thi    s file is generated by volk_fec_profile.\n\
#the     function name is followed by the preferred architecture.\n\
";

        BOOST_FOREACH(std::string result, results) {
            config << result << std::endl;
        }
        config.close();
    }
    else {
        std::cout << "Warning: config not generated" << std::endl;
    }
}
