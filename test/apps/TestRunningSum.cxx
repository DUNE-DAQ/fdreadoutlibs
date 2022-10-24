/**
 * @file WIB2TestBench.cxx  Example of expanding a WIB2 frame using AVX2 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "detdataformats/wib2/WIB2Frame.hpp"
#include "utils/swtpg_avx.hpp"
#include "utils/swtpg_running_sum_avx.hpp"

#include "utils/swtpg_naive.hpp"
#include "utils/FrameExpand.hpp"

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
#include <memory>







int main()
{



    // PARAMS of the SWTPG     
    const uint16_t m_tpg_threshold = 5;                    // units of sigma 
    const uint8_t m_tpg_tap_exponent = 6; 
    const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
    std::vector<int16_t> m_tpg_taps;                       // firwin_int(7, 0.1, multiplier);


    uint16_t* m_primfind_dest = nullptr;  
    if (m_primfind_dest == nullptr) {
      m_primfind_dest = new uint16_t[100000]; 
    }
    int16_t* m_tpg_taps_p = nullptr;
    if (m_tpg_taps_p == nullptr) {
      m_tpg_taps_p = new int16_t[m_tpg_taps.size()];
    }
    m_tpg_taps = firwin_int(7, 0.1, m_tpg_multiplier);
    m_tpg_taps.push_back(0);    


    std::unique_ptr<ProcessingInfo<NUM_REGISTERS_PER_FRAME>> tpg_processing_info = std::make_unique<ProcessingInfo<NUM_REGISTERS_PER_FRAME>>(
      nullptr,
      FRAMES_PER_MSG,
      0,
      NUM_REGISTERS_PER_FRAME,
      m_primfind_dest,
      m_tpg_taps_p,
      (uint8_t)m_tpg_taps.size(), 
      m_tpg_tap_exponent,
      m_tpg_threshold,
      0,
      0
    );



    dunedaq::detdataformats::wib2::WIB2Frame frame;

    // Zero it out first
    std::memset(&frame, 0, sizeof(dunedaq::detdataformats::wib2::WIB2Frame));

    for(int i=0; i<256; ++i){
        frame.set_adc(i, 0x3a0+i);
    }

    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    // Frame expansion
    RegisterArray<16> unpacked=unpack(frame);

    tpg_processing_info->setState(unpacked);

    tpg_processing_info->input = &unpacked;
    *tpg_processing_info->output = MAGIC;
    process_window_avx2_running_sum(*tpg_processing_info);
    process_window_avx2(*tpg_processing_info);
    
    process_window_naive(*tpg_processing_info);

    
    printf("\n");





  



    std::cout << "Finished testing." << std::endl;




}


