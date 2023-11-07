/**
 * @file WIBEthTestBench.cxx  Example of expanding a WIBEth frame using AVX2 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/WIBEthFrame.hpp"
#include "fdreadoutlibs/wibeth/tpg/FrameExpand.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>


#include <cstring>
#include <immintrin.h>
#include <cstdio> // For printf
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>



swtpg_wibeth::RegisterArray<4*64> unpack_wibeth( dunedaq::fddetdataformats::WIBEthFrame& frame)
{

    // Number of time samples (TS) per frame
    int time_samples_per_frame = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;
   
    // Number of ADC words per TS
    int num_adc_words_per_ts = dunedaq::fddetdataformats::WIBEthFrame::s_num_adc_words_per_ts;

    
    swtpg_wibeth::RegisterArray<4*64> ret;

    // Define a pointer to walk the rows of the WIBEth frame which is 2D array.
    //dunedaq::fddetdataformats::WIBEthFrame::word_t (*ptr)[14] = frame.adc_words;
    auto frame_words_ptr = frame.adc_words;

    for (int i = 0; i < time_samples_per_frame; i++) {

      // The register index is used to decide on which of the 
      // 4 registers we want to unpack the ADC messages
      int reg_index = 0;

      for (int j = 0; j < num_adc_words_per_ts; j++) {

          // The words repeat every 7 iterations. 
          // In this way we can use the same unpacking 
          // function (unpack_one_register) as for the DUNE WIBs
          if (j%7 == 0 ) {
            dunedaq::fddetdataformats::WIBEthFrame::word_t * first_half = (*(frame_words_ptr + i) + j);
                               

            // Unpack one register and add it to the register array
            ret.set_ymm(i+reg_index*time_samples_per_frame, swtpg_wibeth::unpack_one_register(first_half));

            reg_index += 1;

            // Critical part: increment the cursor by 224 bits to get the second part of the first time sample
            // 224 corresponds to 16 (U blocks or ADCs) times 14 which are the bits per ADC in the frame. 
            // Check the WIBEth spreadsheet for further details
            char* cursor = (char*) first_half;
            cursor += 224 / 8; // divide by 8 to get the results in bytes
            dunedaq::fddetdataformats::WIBEthFrame::word_t * second_half = (dunedaq::fddetdataformats::WIBEthFrame::word_t*) cursor;
            // Unpack another register and add it to the register array
            ret.set_ymm(i+reg_index*time_samples_per_frame, swtpg_wibeth::unpack_one_register(second_half));

            reg_index += 1;

    
          }       
          
      } // loop over number of words 
    } // loop over time frames
          
        

    return ret;
}




int main()
{
    
    // Test WIBEth unpacking

    dunedaq::fddetdataformats::WIBEthFrame frame;

    // Zero it out first
    std::memset(&frame, 0, sizeof(dunedaq::fddetdataformats::WIBEthFrame));

    // Number of channels per frame
    int num_channels = dunedaq::fddetdataformats::WIBEthFrame::s_channels_per_half_femb;

    // Number of time samples (TS) per frame
    int time_samples_per_frame = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;
   


    // First, loop over time frames
    for (int itime=0; itime<time_samples_per_frame; ++itime) {
      // Loop over all the channels  
      for(int i=0; i<num_channels; ++i){
        frame.set_adc(i, itime, i);           
      }
    }  

    // Check that the ADC values have been properly set
    //for (int itime=0; itime<time_samples_per_frame; ++itime) {
    //  for(int j=0; j<num_channels; ++j){
    //      uint16_t adc_val = frame.get_adc(j, itime);
    //      std::cout << "Index " << j << " time sample " << itime << " : " << adc_val << std::endl;
    //  }
    //}


    
    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    // Unpack the WIB Eth frame
    swtpg_wibeth::RegisterArray<4*64> unpacked = unpack_wibeth(frame);

    // AAA: TODO: the following quantities are hardcode for the moment. 
    int NREGISTERS = swtpg_wibeth::NUM_REGISTERS_PER_FRAME;
    int SAMPLES_PER_REGISTER = swtpg_wibeth::SAMPLES_PER_REGISTER;
    int TIME_FRAMES = time_samples_per_frame;
    
    bool success=true;

	    
    for (size_t j = 0; j < NREGISTERS * SAMPLES_PER_REGISTER; ++j) {
      // Index of the ADC frame 
      int in_index=16*(j/16)+indices[j%16];
      
      const size_t register_offset = j % SAMPLES_PER_REGISTER; 
      const size_t register_index = j / SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * TIME_FRAMES;

      int16_t out_val;
      for (size_t itime = 0; itime < TIME_FRAMES; ++itime) {
        const size_t msg_index = itime / TIME_FRAMES;
        const size_t msg_time_offset = itime % TIME_FRAMES;
        
        const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        // The index in uint16_t of the start of the message we want. 
        const size_t index = msg_start_index + offset_within_msg;

        //const uint16_t* input16 = unpacked.data(); // NOLINT
        //out_val = input16[index]
        out_val = unpacked.uint16(index);
        // Input value into the frame
        uint16_t in_val=frame.get_adc(in_index, itime);

        std::cout << "Time sample: " << itime << " index OUTPUT: " << index << "    value:   " << out_val << std::endl;
        std::cout << "Time sample: " << itime << " index  INPUT: " << in_index << "    value:   " << in_val << std::endl;
        std::cout << "============" << std::endl;

        if(in_val!=out_val){
            success=false;  
        }
      }
    }

    if (success) {
        std::cout << "\n\nALL THE ADC VALUES ARE MATCHING!\n\n" << std::endl;
    }
    

}

