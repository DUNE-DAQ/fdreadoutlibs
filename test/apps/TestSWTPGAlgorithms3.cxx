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
//                       GENERAL
// =================================================================

struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};


// Set CPU affinity of the processing thread
void SetAffinityThread(int executorId) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(executorId, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
       std::cerr << "Error   calling pthread_setaffinity_np Readout: " << rc << "\n";
    }
}


// =================================================================
//                       TPG emulator class
// =================================================================


class tpg_emulator {
public:
    tpg_emulator(std::string file_path, int num_frames, int tpg_threshold) {
        m_frame_file_path = file_path;
        m_num_frames = num_frames;
        m_tpg_threshold = tpg_threshold;
    }; 

    void setup() {

      //                       READ THE FRAMES.BIN FILE
      const std::string input_file = m_frame_file_path;
      
      m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
   
      m_source_buffer->read(input_file);
      int num_superchunks = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 
  
      std::cout << "Number of superchunks in the input file: " << num_superchunks << std::endl;
  
      // Check if the selected number of frames is <= than the ones available in the input file
      if (num_superchunks < m_num_frames) {
        throw "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file.";        
      } else if (m_num_frames == -1) {
        m_num_frames = num_superchunks;
      }


      //                       SETUP THE SWTPG    
      // Populate the queue of primfind destinations
      for (size_t i=0; i<m_capacity_dest_queue; ++i) {
        uint16_t* dest = new uint16_t[100000000];
        m_dest_queue.push(std::move(dest), std::chrono::milliseconds(0));    
      }

      // Initialize the frame handler will setup the processing info class
      fh->initialize(m_tpg_threshold);


    }; 

    void execute_tpg() {
      
      SetAffinityThread(m_tpg_execution_core);

      // Add the rate limiter and fix it to 166 kHz 
      auto limiter = dunedaq::readoutlibs::RateLimiter(166);
      limiter.init();
  
      auto start_test = std::chrono::high_resolution_clock::now();
      //auto now = start_test;
      //while (std::chrono::duration_cast<std::chrono::seconds>(now - start_test).count() < m_execution_time) {
      while (true) {  

        auto& source = m_source_buffer->get(); 
        // current superchunk
        auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + m_superchunk_index*swtpg_wib2::SUPERCHUNK_FRAME_SIZE);

        // Parse the DUNEWIB frames
        auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp);
        uint64_t timestamp = wfptr->get_timestamp();   

        swtpg_wib2::MessageRegisters registers_array;
        swtpg_wib2::expand_wib2_adcs(fp, &registers_array, 0);
      
        if (m_first_hit) {                     
          fh->m_tpg_processing_info->setState(registers_array);
          m_first_hit = false; 
          m_extract_tp_thread_should_run.store(true);
      
          //Save ADC info
          //if (m_save_adc_data){
          //  save_raw_data(registers_array, timestamp, -1, "SwtpgNaive");
          //}          
        }  

  
        fh->m_tpg_processing_info->input = &registers_array;
        uint16_t* destination_ptr = fh->get_primfind_dest();
        *destination_ptr = swtpg_wib2::MAGIC;
        fh->m_tpg_processing_info->output = destination_ptr;
        swtpg_wib2::process_window_avx2(*fh->m_tpg_processing_info, 0);       

        // Insert output of the AVX processing into the swtpg_output 
        swtpg_output swtpg_processing_result = {destination_ptr, timestamp};

        m_tphandler_queue.try_push(std::move(swtpg_processing_result), std::chrono::milliseconds(0));


        ++m_superchunk_index;
        //std::cout << m_superchunk_index << std::endl;

        if (m_superchunk_index == m_num_frames) {
          m_superchunk_index = 0;
          // Calculate elapsed time in seconds  
          auto now = std::chrono::high_resolution_clock::now();
          auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - start_test).count();  
          std::cout << "Elapsed time: " << elapsed_seconds << std::endl;
          //now = std::chrono::high_resolution_clock::now();

        }

        //m_dest_queue.try_push(std::move(destination_ptr), std::chrono::milliseconds(0));    

 
        limiter.limit();
      }  



    }; 
    void extract_tp() {
      SetAffinityThread(m_extract_tp_core);
    
      swtpg_output result_from_swtpg; 
      unsigned int nhits = 0;
      while (m_extract_tp_thread_should_run.load()) {
      while(m_tphandler_queue.can_pop()) {
        if(m_tphandler_queue.try_pop(result_from_swtpg, std::chrono::milliseconds(0))) {
    
    
        uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; 
        nhits = 0;
      
          while (*result_from_swtpg.output_location != swtpg_wib2::MAGIC) {
            for (int i = 0; i < 16; ++i) {
              chan[i] = *result_from_swtpg.output_location++; 
            }
            for (int i = 0; i < 16; ++i) {
              hit_end[i] = *result_from_swtpg.output_location++; 
            }
            for (int i = 0; i < 16; ++i) {
              hit_charge[i] = *result_from_swtpg.output_location++;
            }
            for (int i = 0; i < 16; ++i) {        
              hit_tover[i] = *result_from_swtpg.output_location++; 
            }  
            
            // Now that we have all the register values in local
            // variables, loop over the register index (ie, channel) and
            // find the channels which actually had a hit, as indicated by
            // nonzero value of hit_charge
            for (int i = 0; i < 16; ++i) {
              if (hit_charge[i] && chan[i] != swtpg_wib2::MAGIC) {
                //std::cout << "Channel number: " << chan[i] << std::endl;
                //std::cout << "Hit charge: " << hit_charge[i] << std::endl
                ++nhits;
              }
            } // loop over 16 registers   
          } // while not magic        

          //std::cout << "Found " << nhits << " hits " << std::endl;
          m_dest_queue.push(std::move(result_from_swtpg.output_location), std::chrono::milliseconds(0));    

        } // try_pop tphandler  
      } // can_pop tphandler  
      } // should run
          
  

    }; 
    
private:
    std::string m_frame_file_path;
    int m_num_frames;
    int m_tpg_threshold;
    int m_execution_time = 10;
    int m_superchunk_index = 0;
    int m_tpg_execution_core = 0;
    int m_extract_tp_core = 1;

    bool m_first_hit = true;
    bool m_save_adc_data = false;
    std::atomic<bool> m_extract_tp_thread_should_run{false};


    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    
    
    int m_capacity_dest_queue = 100;
    dunedaq::iomanager::FollyMPMCQueue<uint16_t*> m_dest_queue{"dest_queue", m_capacity_dest_queue};

    // Declare frame handler class
    std::unique_ptr<dunedaq::fdreadoutlibs::WIB2FrameHandler> fh = std::make_unique<dunedaq::fdreadoutlibs::WIB2FrameHandler>(0, m_dest_queue);

    size_t m_capacity_mpmc_queue = 300000000; 
    dunedaq::iomanager::FollyMPMCQueue<swtpg_output> m_tphandler_queue{"tphandler_queue", m_capacity_mpmc_queue};


    
};



// =================================================================
//                       MAIN
// =================================================================
// Usage: ./TestSWTPGAlgorithms3  --frame_file_path /nfs/sw/work_dirs/aabedabu/dunedaq-v4.0.0.candidate/dev/configurations/frames_wib2.bin --algorithm SWTPG --implementation AVX  --swtpg_threshold 500


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



    class tpg_emulator emulator01(frame_file_path, num_frames, swtpg_threshold);

    emulator01.setup();
    //emulator01.execute_tpg();




    //auto future = std::async(std::launch::async, [&emulator01]() {
    //    emulator01.extract_tp();
    //});



    std::thread execution_thread(&tpg_emulator::execute_tpg, &emulator01);
    execution_thread.detach();    

    std::thread extract_tp_thread(&tpg_emulator::extract_tp, &emulator01);
    extract_tp_thread.join();    



    std::cout << "\n\nFinished testing." << std::endl;


}