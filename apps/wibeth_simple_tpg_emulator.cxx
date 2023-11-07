/**
 * @file TestTPGAlgorithms.cxx Main file for testing different tpg algorithms 
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "fddetdataformats/WIBEthFrame.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"

#include "iomanager/IOManager.hpp"

#include "readoutlibs/utils/RateLimiter.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp"

#include "fdreadoutlibs/wibeth/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessRSAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaive.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

#include "triggeralgs/TriggerPrimitive.hpp"

// Local
#include "tpg_apps_utilities.hpp"
#include "tpg_apps_issues.hpp"



// system
#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>


#include <cstring>
#include <immintrin.h>
#include <cstdio> 
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>
#include <memory>




// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "Test TPG algorithms" };
    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");

    std::string select_algorithm = "SimpleThreshold";
    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS). Default: SimpleThreshold");
  
    std::string select_implementation = "AVX";
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE). Default: AVX");

    int duration_test = 120; 
    app.add_option("-d,--duration-test", duration_test, "Duration (in seconds) to run the test. Default value is 120.");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default value is 500.");

    int core_number = 0;
    app.add_option("-c,--core", core_number, "Set core number of the executing TPG thread. Default value is 0.");

    bool save_adc_data = false;
    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data (first frame only)");

    bool save_trigprim = false;
    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");


    CLI11_PARSE(app, argc, argv);



    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================

    const int wibeth_frame_size = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::fixed_payload_size;
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    // Read only 10 MB of data
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, wibeth_frame_size); 

    m_source_buffer->read(file_path_input);
    auto& source = m_source_buffer->get(); 
    const int total_num_frames = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of DUNE WIBEth frames in the input file: " << total_num_frames << std::endl;

    
    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    tpg_emulator emulator(save_adc_data, save_trigprim, false, select_algorithm, select_implementation, "") ;
    emulator.set_tpg_threshold(tpg_threshold);
    emulator.set_CPU_affinity(core_number);
    emulator.initialize();


    // Setup the rate limiter for WIBEth frames
    auto limiter = dunedaq::readoutlibs::RateLimiter(31);
    limiter.init();



    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0; 
    uint64_t frame_repeat_index = 0;
    auto start_test = std::chrono::high_resolution_clock::now();  

    // Loop over the DUNEWIB Ethernet frames in the file
    while (wibeth_frame_index < total_num_frames ){      

      // current WIBEth frame
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(source.data() + wibeth_frame_index*wibeth_frame_size);

      emulator.execute_tpg(fp);

      ++wibeth_frame_index;


      // If end of the file is reached, restart the index counter
      if (wibeth_frame_index == total_num_frames) {
        wibeth_frame_index = 0;
	      frame_repeat_index++;
      }

      if (frame_repeat_index % 500  == 0) {
        // Calculate elapsed time in seconds  
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - start_test).count();  
        std::cout << "Elapsed time [s]: " << elapsed_seconds << std::endl;      
	
	      frame_repeat_index = 0;        

        // stop the testing after a time a condition
        if (elapsed_seconds > duration_test) {
          wibeth_frame_index = total_num_frames;
        }
      }

	

      limiter.limit();

    }
    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << emulator.get_total_hit_number() << " hits." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;

}


