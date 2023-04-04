/**
 * @file utils.hpp  Utilities for the SWTPG code testing
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TEST_SWTPG_UTILITIES_HPPP_
#define TEST_SWTPG_UTILITIES_HPPP_

#include "triggeralgs/TriggerPrimitive.hpp"

#include <iostream>
#include <sstream>
#include <ctime>
#include <cstring>


typedef swtpg_wib2::RegisterArray<swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::FRAMES_PER_MSG> MessageRegisters;
using timestamp_t = std::uint64_t; 

// Function save the TP data to a file 
void save_hit_data( triggeralgs::TriggerPrimitive trigprim, std::string source_name ){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + source_name + "_" + date_time_str + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  out_file << trigprim.channel << "," << trigprim.time_start << "," << trigprim.time_over_threshold << "," 
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak << "\n";  

  out_file.close();
}


// Function to save raw ADC data to a file (only for debugging) 
void save_raw_data(swtpg_wib2::MessageRegisters register_array, 
	       uint64_t t0, size_t channel_number,
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
  for (size_t ichan = 0; ichan < swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER; ++ichan) {
    const size_t register_index = ichan / swtpg_wib2::SAMPLES_PER_REGISTER;
    if (register_index < 0 || register_index >= swtpg_wib2::NUM_REGISTERS_PER_FRAME)
       continue;

    // Parse only selected channel number. To select all channels choose -1
    if (ichan == channel_number || channel_number == -1) { 
   
      const size_t register_offset = ichan % swtpg_wib2::SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * swtpg_wib2::SAMPLES_PER_REGISTER * swtpg_wib2::FRAMES_PER_MSG;
  
      for (size_t iframe = 0; iframe<swtpg_wib2::FRAMES_PER_MSG; ++iframe) {
    
        const size_t msg_index = iframe / 12;
        const size_t msg_time_offset = iframe % 12;
        // The index in uint16_t of the start of the message we want // NOLINT 
        const size_t msg_start_index = msg_index * swtpg_wib2::ADCS_SIZE / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + swtpg_wib2::SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;
    
        int16_t adc_value = input16[index];
        //std::cout << adc_value << std::endl;
        out_file << ichan << "," <<  adc_value << "," << t_current << std::endl;
        t_current += 32;
      } 

    }
  }
  out_file.close();


}

#endif // TEST_SWTPG_UTILITIES_HPPP_
