/**
 * @file WIB2TestBench.cxx  Example of expanding a WIB2 frame using AVX2 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "detdataformats/wib2/WIB2Frame.hpp"

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

#include "iomanager/IOManager.hpp"

#include "fdreadoutlibs/wib2/tpg/FrameExpand512.hpp"
#include "fdreadoutlibs/wib2/WIB2FrameProcessor.hpp"
#include "fdreadoutlibs/wib2/tpg//ProcessAVX512.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

using namespace swtpg_wib2;
using namespace dunedaq::fdreadoutlibs;

typedef RegisterArray512<16> RegisterSet;
# define RAW_WORDS_PER_REG 14
# define NREGISTERS 16



RegisterSet unpack( dunedaq::detdataformats::wib2::WIB2Frame& frame)
{
    RegisterSet ret;
    for(size_t i=0; i<ret.size(); ++i){
        ret.set_ymm(i, unpack_one_register512(frame.adc_words+RAW_WORDS_PER_REG*i));
    }
    return ret;
}


bool in_out_test(const std::array<uint16_t, 256>& vals)
{
    dunedaq::detdataformats::wib2::WIB2Frame frame;
    // Zero it out first
    std::memset(&frame, 0, sizeof(dunedaq::detdataformats::wib2::WIB2Frame));

    for(size_t i=0; i<vals.size(); ++i){
        frame.set_adc(i, vals[i]);
    }

    bool success=true;
    for(size_t i=0; i<vals.size(); ++i){
        uint16_t out_val=frame.get_adc(i);
        uint16_t in_val=vals[i];
        if(out_val!=in_val){
            success=false;
            printf("Failure at index %zu. Input: %03x, output: %03x\n", i, in_val, out_val);
        }
    }
    return success;
}


int main()
{


    // Some preliminary testing when creating an input array 	
    {
        printf("Test 1:\n");
        std::array<uint16_t, 256> vals;
        for(size_t i=0; i<vals.size(); ++i) vals[i]=i;
        in_out_test(vals);
    }

    {
        printf("Test 2:\n");
        std::array<uint16_t, 256> vals;
        for(size_t i=0; i<vals.size(); ++i) vals[i]=0x3a0+i;
        in_out_test(vals);
    }
      

    
    // Test WIB unpacking

    dunedaq::detdataformats::wib2::WIB2Frame frame;

    // Zero it out first
    std::memset(&frame, 0, sizeof(dunedaq::detdataformats::wib2::WIB2Frame));
    std::cout << "Setting adc values to:" << std::hex;
    for(int i=0; i<256; ++i){
        frame.set_adc(i, 0x3a0+i);
	std::cout << " " << 0x3a0+i;
    }
    std::cout << std::dec << std::endl;

    std::cout << "============  Calling unpack =========\n";
    RegisterSet unpacked=unpack(frame);

    std::array<int, 32> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14,
                                16, 17, 18, 19, 20, 21, 22, 23, 31, 24, 25, 26, 27, 28, 29, 30};
    print512_as16(unpacked.ymm(0));
    std::cout << std::endl;
    for (int i=0; i<32; ++i) std::cout << std::hex << unpacked.uint16(i) << " ";
    std::cout << std::endl;

    //print256(unpacked);
    //printf("\n");
    int nindices=indices.size();
    bool success=true;
    for(int i=0; i<256; ++i){
        int in_index=nindices*(i/nindices)+indices[i%nindices];
        uint16_t in_val=frame.get_adc(in_index);
	// std::cout << in_val << std::endl;
        uint16_t out_val=unpacked.uint16(i);
        if(in_val!=out_val){
           printf("%03d %03x %03x\n", in_index, in_val, out_val);
            success=false;
        }
    }
    printf("Success? %d\n", success);


    dunedaq::detdataformats::wib2::WIB2Frame frames[swtpg_wib2::FRAMES_PER_MSG];
    int count=0;
    for (unsigned int f=0; f<swtpg_wib2::FRAMES_PER_MSG; f++) {
       std::cout << "New frame starts with " << std::hex << 0x3a0+count << std::dec << std::endl;
       if (f==0) {
         for (int i=0; i<256; ++i) {
           frames[f].set_adc(i, 666);
        } } else {
          for (int i=0; i<256; ++i){
            frames[f].set_adc(i, i);        
          }
       }
    }

    int threshold_value = 120;
    WIB2FrameHandler fh(0);
    fh.initialize(threshold_value);
    swtpg_wib2::MessageRegisters512 registers_array;
    std::cout << "Calling expand_wib2_adcs()\n";
    expand_wib2_adcs_512(reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(&frames),
                     &registers_array, 0);


    for (size_t j = 0; j < NREGISTERS * swtpg_wib2::SAMPLES_PER_REGISTER; ++j) {
      const size_t register_offset = j % swtpg_wib2::SAMPLES_PER_REGISTER; 
      const size_t register_index = j / swtpg_wib2::SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * swtpg_wib2::SAMPLES_PER_REGISTER * swtpg_wib2::FRAMES_PER_MSG;

      int16_t ped;
      for (size_t itime = 0; itime < swtpg_wib2::FRAMES_PER_MSG; ++itime) {
        const size_t msg_index = itime / 12;
        const size_t msg_time_offset = itime % 12;
        // The index in uint16_t of the start of the message we want. 
        // AAA: 6144 is the sizeof(MessageADCs) which is hardcoded here. TODO: put
        // MessageADCs in a location that can be included and used in the ProcessingInfo 
        const size_t msg_start_index = msg_index * swtpg_wib2::ADCS_SIZE / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + swtpg_wib2::SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;

        const uint16_t* input16 = registers_array.data(); // NOLINT

        ped = input16[index];
        std::cout << ped << " ";
      }
      std::cout << std::endl;;
    }       
     


    /*                 
    fh.m_tpg_processing_info->setState(registers_array);
    fh.m_tpg_processing_info->input = &registers_array;
    auto primPtr=fh.get_primfind_dest();
    *fh.get_primfind_dest() = swtpg_wib2::MAGIC;
    
    std::cout << "Calling process_window_avx512()\n";
    swtpg_wib2::process_window_avx512(*fh.m_tpg_processing_info);

    int nfound=0;
    while (*primPtr != swtpg_wib2::MAGIC) {
      nfound++;
      primPtr++;
    }
    std::cout << "Found " << nfound << " words of primitives\n";

    fh.reset();
    */

    std::cout << "Finished" << std::endl;




}

