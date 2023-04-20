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

#include "readoutlibs/utils/RateLimiter.hpp"
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



// =================================================================
//                       PUBLIC VARIABLES
// =================================================================
int m_capacity_dest_queue = 100;
dunedaq::iomanager::FollyMPMCQueue<uint16_t*> m_dest_queue{"dest_queue", m_capacity_dest_queue};
dunedaq::fdreadoutlibs::WIB2FrameHandler fh(0, m_dest_queue);

unsigned int total_hits = 0;
bool first_hit = true;


// =================================================================
//                       FUNCTIONS
// =================================================================

// Set CPU affinity of the processing thread
void SetAffinityThread(int executorId) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(executorId, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
       std::cerr << "Error calling pthread_setaffinity_np Readout: " << rc << "\n";
    }
}


// =================================================================
//                       TPG FUNCTIONS
// =================================================================

void extract_nhits(uint16_t* primfind_it) {
  SetAffinityThread(1);

  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; 
  unsigned int nhits = 0;

  while (*primfind_it != swtpg_wib2::MAGIC) {
    for (int i = 0; i < 16; ++i) {
      chan[i] = *primfind_it++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_end[i] = *primfind_it++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_charge[i] = *primfind_it++;
    }
    for (int i = 0; i < 16; ++i) {        
      hit_tover[i] = *primfind_it++; 
    }  
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (int i = 0; i < 16; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wib2::MAGIC) 
        //std::cout << "Channel number: " << chan[i] << std::endl;
        //std::cout << "Hit charge: " << hit_charge[i] << std::endl
        ++nhits;
      }
    }    

  //std::cout << "Found " << nhits << " hits " << std::endl;
  //total_hits += nhits;

  m_dest_queue.push(std::move(primfind_it), std::chrono::milliseconds(0));    
  
}


uint16_t* execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp, bool save_adc_data) {

  SetAffinityThread(0);

  // Parse the DUNEWIB frames
  //auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp);
  //uint64_t timestamp = wfptr->get_timestamp();      
  swtpg_wib2::MessageRegisters registers_array;
  swtpg_wib2::expand_wib2_adcs(fp, &registers_array, 0);

  
  if (first_hit) {                     
    fh.m_tpg_processing_info->setState(registers_array);
    first_hit = false;    

    //Save ADC info
    //if (save_adc_data){
    //  save_raw_data(registers_array, timestamp, -1, "SwtpgNaive");
    //}


  }  
       
  
  fh.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = fh.get_primfind_dest();
  *destination_ptr = swtpg_wib2::MAGIC;
  fh.m_tpg_processing_info->output = destination_ptr;
  swtpg_wib2::process_window_avx2(*fh.m_tpg_processing_info, 0);                  

  return destination_ptr;

}


// =================================================================
//                       MAIN
// =================================================================


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

    // Populate the queue of primfind destinations
    for (size_t i=0; i<m_capacity_dest_queue; ++i) {
      uint16_t* dest = new uint16_t[1000000];
      m_dest_queue.push(std::move(dest), std::chrono::milliseconds(0));    
    }

    auto limiter = dunedaq::readoutlibs::RateLimiter(166);
    limiter.init();

    fh.initialize(swtpg_threshold);

    // =================================================================
    //                       Process the DUNEWIB superchunks
    // =================================================================
    int superchunk_index = 0;   
    // Loop over the DUNEWIB superchunks in the file
    //while (superchunk_index < num_frames ){      
    while (true){        
      
      // current superchunk
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + superchunk_index*swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
      uint16_t* destination_ptr = execute_tpg(fp, save_adc_data);
      extract_nhits(destination_ptr);


      ++superchunk_index;
      if (superchunk_index % 1000 == 0) {
        std::cout << "Executing superchunk number " << superchunk_index << " out of " << num_superchunks << std::endl;
      }
      if (superchunk_index == num_frames) {
        superchunk_index = 0;
      }
      limiter.limit();

    }

    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << total_hits << " hits." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;



}


