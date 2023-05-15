/**
 * @file TestSWTPGAlgorithms.cxx Main file for testing different swtpg algorithms 
 * @author Adam Abed Abud (adam.abed.abud@cern.ch)
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "detdataformats/wib2/WIB2Frame.hpp"

#include "iomanager/IOManager.hpp"

#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib2/WIB2FrameProcessor.hpp"

#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"



// Local
#include "SwtpgTest/SwtpgBase.hpp"
#include "SwtpgTest/SwtpgNaive.hpp"
#include "SwtpgTest/SwtpgAvx.hpp"
#include "SwtpgTest/RSNaive.hpp"
#include "SwtpgTest/RSAvx.hpp"


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




int 
main(int argc, char** argv)
{

    CLI::App app{ "Test SWTPG algorithms" };

    // Set default input frame file
    std::string frame_file_path = "./frames_wib2.bin";
    app.add_option("-f,--frame_file_path", frame_file_path, "Path to the input frame file");

    std::string select_algorithm;
    app.add_option("-a,--algorithm", select_algorithm, "SWTPG Algorithm (SWTPG / RS)");
  
    std::string select_implementation;
    app.add_option("-i,--implementation", select_implementation, "SWTPG implementation (AVX / NAIVE)");

    int num_frames = -1;
    app.add_option("-n,--num_frames", num_frames, "Number of frames to read. Default: select all frames.");

    int swtpg_threshold = 100;
    app.add_option("-t,--swtpg_threshold", swtpg_threshold, "Value of the SWTPG threshold");

    std::string swtpg_channel_map = "VDColdboxChannelMap";
    app.add_option("-m,--swtpg_channel_map", swtpg_channel_map, "Name of the SWTPG offline channel map (or None), default: VDColdboxChannelMap");

    bool save_adc_data{false};
    app.add_option("--save_adc_data", save_adc_data, "Save ADC data (true/false)");

    bool save_hit_data{false};
    app.add_option("--save_hit_data", save_hit_data, "Save hit data (true/false)");


    CLI11_PARSE(app, argc, argv);

    std::unique_ptr<SwtpgBase> algo;

    if (select_algorithm == "SWTPG") {
      if (select_implementation == "NAIVE") {
        algo = std::make_unique<SwtpgNaive>(save_adc_data, save_hit_data);
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        algo = std::make_unique<SwtpgAvx>(save_adc_data, save_hit_data);
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else {
        std::cout << "Select a valid algorithm implementation. Use --help for further details." << std::endl;
        return 1;
      }
    } else if (select_algorithm == "RS") {
      if (select_implementation == "NAIVE") {
        algo = std::make_unique<RSNaive>(save_adc_data, save_hit_data);
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        algo = std::make_unique<RSAvx>(save_adc_data, save_hit_data);
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else {
        std::cout << "Select a valid algorithm implementation. Use --help for further details." << std::endl;
        return 1;
      }
    } else {
      std::cout << "Select at least an algorithm. Use --help for further details." << std::endl;
      return 1;
    }


    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    const std::string input_file = frame_file_path;
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
 
    m_source_buffer->read(input_file);
    auto& source = m_source_buffer->get(); 
    int num_superchunks = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of superchunks in the input file: " << num_superchunks << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (num_superchunks < num_frames) {
      std::cout << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_frames == -1) {
      num_frames = num_superchunks;
    }


    // =================================================================
    //                       Setup the SWTPG
    // =================================================================
    
    bool first_hit = true;




    // =================================================================
    //                       Process the DUNEWIB superchunks
    // =================================================================
    int superchunk_index = 0;   

    // Loop over the DUNEWIB superchunks in the file
    while (superchunk_index<num_frames ){

      // current superchunk
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + superchunk_index*swtpg_wib2::SUPERCHUNK_FRAME_SIZE);

      // Reset the memory buffers
      algo->reset(first_hit, swtpg_threshold, swtpg_channel_map);
      // Find the SWTPG hits
      algo->find_hits(fp, first_hit);
      first_hit = false;        

      
      ++superchunk_index;
      std::cout << "Executing superchunk number " << superchunk_index << " out of " << num_superchunks << std::endl;

    }

    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << algo->m_total_swtpg_hits << " hits." << std::endl;
    
    

    std::cout << "\n\nFinished testing." << std::endl;



}


