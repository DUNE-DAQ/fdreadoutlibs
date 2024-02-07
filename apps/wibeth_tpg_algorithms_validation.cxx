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
//#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

#include "triggeralgs/TriggerPrimitive.hpp"


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

struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};

int WIBEth_FRAME_SIZE = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::fixed_payload_size;

int duration_test = 120; // default value
unsigned int total_hits = 0;
bool first_hit = true;

dunedaq::fdreadoutlibs::WIBEthFrameHandler fh;

std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;

std::string select_algorithm = "";
std::string select_implementation = "";
bool save_adc_data = false;
bool save_trigprim = false;
std::string out_suffix = "";

// =================================================================
//                       FUNCTIONS and UTILITIES
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

// Function save the TP data to a file 
void save_hit_data( triggeralgs::TriggerPrimitive trigprim, std::string source_name, std::string out_suffix=""){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + source_name + "_" + date_time_str + out_suffix + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  /*out_file << trigprim.time_start 
           << " " << trigprim.time_over_threshold 
           << " " << trigprim.time_peak 
           << " " << trigprim.channel 
           << " " << trigprim.adc_integral 
           << " "  << trigprim.adc_peak << "\n";*/
  // DBG temporary for the sake of comparison 
  out_file << trigprim.channel
	   << "," << trigprim.time_start 
           << "," << trigprim.time_over_threshold 
           << "," << trigprim.time_peak 
           << "," << trigprim.adc_integral 
           << ","  << trigprim.adc_peak << "\n";  


  out_file.close();
}

// Function to save raw ADC data to a file (only for debugging) 
void save_raw_data(swtpg_wibeth::MessageRegisters register_array, 
	       uint64_t t0, int channel_number,
           std::string source_name)
{
  std::ofstream out_file;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  
  std::string file_name;
  if (channel_number == -1) {
    file_name = "all_channels_" + source_name + "_data" + date_time_str + ".txt";
  } else {
    file_name = "Channel_" + std::to_string(channel_number) + "_" + source_name + "_data" + date_time_str + ".txt";
  }
  out_file.open(file_name.c_str(), std::ofstream::app);

  uint64_t t_current= t0 ; 
  
  const uint16_t* input16 = register_array.data();
  for (int ichan = 0; ichan < static_cast<int>(swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER); ++ichan) {
    const size_t register_index = ichan / swtpg_wibeth::SAMPLES_PER_REGISTER;
    if (register_index >= swtpg_wibeth::NUM_REGISTERS_PER_FRAME)
       continue;

    // Parse only selected channel number. To select all channels choose -1
    if (ichan == channel_number || channel_number == -1) { 
   
      const size_t register_offset = ichan % swtpg_wibeth::SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * swtpg_wibeth::SAMPLES_PER_REGISTER * swtpg_wibeth::FRAMES_PER_MSG;
  
      for (size_t iframe = 0; iframe<swtpg_wibeth::FRAMES_PER_MSG; ++iframe) {
    
        const size_t msg_index = iframe / swtpg_wibeth::FRAMES_PER_MSG; 
        const size_t msg_time_offset = iframe % swtpg_wibeth::FRAMES_PER_MSG;
        // The index in uint16_t of the start of the message we want // NOLINT 
        const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + swtpg_wibeth::SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;
    
        int16_t adc_value = input16[index];
        std::cout << adc_value << std::endl;
        out_file << ichan << "," <<  adc_value << "," << t_current << std::endl;
        t_current += 32;
      } 

    }
  }
  out_file.close();


}

// =================================================================
//                       TPG FUNCTIONS
// =================================================================
void extract_hits_naive(uint16_t* output_location, uint64_t timestamp, std::string out_suffix) {

    constexpr int clocksPerTPCTick = 32;
    const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
    uint16_t chan, hit_end, hit_peak_adc, hit_charge, hit_tover, hit_peak_time;

    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    //std::cout << "DBG ProcessingInfo: timeWindowNumFrames " << fh.m_tpg_processing_info->timeWindowNumFrames << std::endl;
    //std::cout << "DBG nhits " << fh.m_tpg_processing_info->nhits << std::endl;

    size_t i = 0;
    //while (*output_location != swtpg_wibeth::MAGIC) {
    for (size_t n=0; n<fh.m_tpg_processing_info->nhits; n++) {
      chan            = *output_location++;  // TYPES!!!  int32_t  would be better match to triggeralgs::TriggerPrimitive
      hit_end         = *output_location++;  // TYPES!!!  int16_t  better than  uint16_t
      hit_charge      = *output_location++;  // TYPES!!!  uint32_t  would be best
      hit_tover       = *output_location++;  // TYPES!!!  uint32_t  would be best
      hit_peak_adc    = *output_location++;
      hit_peak_time   = *output_location++;  // TYPES!!!  uint32_t  would be best

      i += 1;
      chan = nreg*(chan/nreg)+indices[chan%nreg];
     
      //std::cout << "DBG chan " << n << ": " << chan << std::endl;
      //std::cout << "DBG hit_end " << n << ": " << hit_end << std::endl;
      //std::cout << "DBG hit_tover " << n << ": " << hit_tover << std::endl;

      uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end - (int64_t)hit_tover);
      uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time;

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = tp_t_peak;

      trigprim.time_over_threshold = uint64_t((hit_tover - 1) * clocksPerTPCTick);

      trigprim.channel = chan;
      trigprim.adc_integral = hit_charge;
      trigprim.adc_peak = hit_peak_adc;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_hit_data(trigprim, "NAIVE", out_suffix);
      }
      ++total_hits;

    }
}


void extract_hits_avx(uint16_t* output_location, uint64_t timestamp, std::string out_suffix) {

  constexpr int clocksPerTPCTick = 32;
  const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
  uint16_t chan[nreg], left[nreg], hit_end[nreg], hit_peak_adc[nreg], hit_charge[nreg], hit_tover[nreg], hit_peak_time[nreg];

  //std::cout << "DBG nhits " << fh.m_tpg_processing_info->nhits << std::endl;

  for (size_t n=0; n<fh.m_tpg_processing_info->nhits; n++) {
    for (std::size_t i = 0; i < nreg; ++i) {
      chan[i] = *output_location++; 
      //std::cout << "DBG chan " << chan[i] << std::endl;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      left[i] = *output_location++;
      //std::cout << "DBG left " << left[i] << std::endl;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_end[i] = *output_location++; 
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_charge[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {        
      hit_tover[i] = *output_location++; 
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_adc[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_time[i] = *output_location++;
    }
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (std::size_t i = 0; i < nreg; ++i) {
      if (hit_charge[i] && left[i] == std::numeric_limits<std::uint16_t>::max() 
          && chan[i] != swtpg_wibeth::MAGIC) {
        //std::cout << "Channel number: " << chan[i] << std::endl;
        //std::cout << "Hit charge: " << hit_charge[i] << std::endl;

        uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end[i] - (int64_t)hit_tover[i]);
        uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time[i];

          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          //
          //TLOG_DEBUG(0) << "Hit: " << tp_t_begin << " " << offline_channel;
          triggeralgs::TriggerPrimitive trigprim;
          trigprim.time_start = tp_t_begin;
          trigprim.time_peak = tp_t_peak;
          trigprim.time_over_threshold = uint64_t((hit_tover[i] - 1) * clocksPerTPCTick);
          trigprim.channel = chan[i];
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_peak_adc[i];
          trigprim.detid = 666;           
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
          trigprim.version = 1;

          if (save_trigprim){
            save_hit_data(trigprim, "AVX", out_suffix);
          }          

        ++total_hits;
      }
    } // loop over 16 registers   

  } // while not magic    
}


void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  // Set CPU affinity of the TPG thread
  SetAffinityThread(0);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();      
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);
  //swtpg_wibeth::parse_wibeth_adcs(&registers_array);

  
  if (first_hit) {                     
    fh.m_tpg_processing_info->setState(registers_array);
    first_hit = false;    

    // Save ADC info
    if (save_adc_data){
      //save_raw_data(registers_array, timestamp, -1, select_algorithm + "_" + select_implementation);
    }


  }  
       
  fh.m_tpg_processing_info->nhits = 0;
  fh.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = fh.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  fh.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*fh.m_tpg_processing_info);
       
  // Insert output of the AVX processing into the swtpg_output 
  //swtpg_output swtpg_processing_result = {destination_ptr, timestamp};
    
  if (select_implementation == "AVX") {
    extract_hits_avx(destination_ptr, timestamp, out_suffix);
  } else if(select_implementation == "NAIVE") {  
    extract_hits_naive(destination_ptr, timestamp, out_suffix);
  }

}


// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "Test TPG algorithms" };

    // Set default input frame file
    std::string frame_file_path = "./wibeth-frames.bin";
    app.add_option("-f,--frame-file-path", frame_file_path, "Path to the input frame file");

    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");
  
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE)");
    
    app.add_option("-d,--duration-test", duration_test, "Duration (in seconds) to run the test");

    int num_frames_to_read = -1;
    app.add_option("-n,--num-frames-to-read", num_frames_to_read, "Number of frames to read. Default: select all frames.");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold");

    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data");

    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    // additional options 
    bool enable_repeat_timer = true;
    app.add_option("-r, --enable_repeat_timer", enable_repeat_timer, "Repeat frame processing while certain time elapsed (true/false). Default: true");
    app.add_option("-s ,--out_suffix", out_suffix, "Append string to output hit file name");
    app.set_config("--config", "app.cfg", "Read a configuration file. Default: app.cfg", false);

    std::cout << "DBG sizeof(word_t): " << sizeof(uint64_t) << std::endl;
    std::cout << "DBG s_bits_per_adc: " << dunedaq::fddetdataformats::WIBEthFrame::s_bits_per_adc << std::endl;
    std::cout << "DBG s_bits_per_word: " << dunedaq::fddetdataformats::WIBEthFrame::s_bits_per_word << std::endl;
    std::cout << "DBG s_time_samples_per_frame: " << dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame << std::endl;
    std::cout << "DBG s_channels_per_half_femb: " << dunedaq::fddetdataformats::WIBEthFrame::s_channels_per_half_femb  << std::endl;
    std::cout << "DBG s_num_channels: " << dunedaq::fddetdataformats::WIBEthFrame::s_num_channels << std::endl;
    std::cout << "DBG s_num_adc_words_per_ts: " << dunedaq::fddetdataformats::WIBEthFrame::s_num_adc_words_per_ts << std::endl;
    std::cout << "DBG s_num_adc_words: " << dunedaq::fddetdataformats::WIBEthFrame::s_num_adc_words << std::endl;


    CLI11_PARSE(app, argc, argv);

    if (select_algorithm == "SimpleThreshold") {
      if (select_implementation == "NAIVE") {
        m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else {
        std::cout << "Select a valid algorithm implementation. Use --help for further details." << std::endl;
        return 1;
      }
    } else if (select_algorithm == "AbsRS") {
      if (select_implementation == "NAIVE") {        
        //m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive_RS<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
        //std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
      } else if (select_implementation == "AVX") {
        m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_rs_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
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
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, WIBEth_FRAME_SIZE); // AAA: hardcoded values!! 
 
    m_source_buffer->read(input_file);
    auto& source = m_source_buffer->get(); 
    int total_num_frames = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of DUNE WIBEth frames in the input file: " << total_num_frames << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (total_num_frames < num_frames_to_read) {
      std::cout << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_frames_to_read == -1) {
      num_frames_to_read = total_num_frames;
    }

    
    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    auto limiter = dunedaq::readoutlibs::RateLimiter(31);
    limiter.init();

    fh.initialize(tpg_threshold);

    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0; 
    uint64_t frame_repeat_index = 0;
    auto start_test = std::chrono::high_resolution_clock::now();  

    // Loop over the DUNEWIB Ethernet frames in the file
    while (wibeth_frame_index < num_frames_to_read ){      

      // current WIBEth frame
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(source.data() + wibeth_frame_index*WIBEth_FRAME_SIZE);

      execute_tpg(fp);

      ++wibeth_frame_index;

      if (!enable_repeat_timer) {
        if (wibeth_frame_index == num_frames_to_read) {
          continue;
        }
      }

      // If end of the file is reached, restart the index counter
      if (wibeth_frame_index == num_frames_to_read) {
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
          wibeth_frame_index = num_frames_to_read;
        }

        limiter.limit();
      }

    }
    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << total_hits << " hits." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;

}


