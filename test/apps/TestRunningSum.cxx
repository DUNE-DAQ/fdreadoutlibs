/**
 * @file WIB2TestBench.cxx  Example of expanding a WIB2 frame using AVX2 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "detdataformats/wib2/WIB2Frame.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"

// local
#include "utils/swtpg_avx.hpp"
#include "utils/swtpg_running_sum_avx.hpp"
#include "utils/swtpg_naive.hpp"
#include "utils/FrameExpand.hpp"

// system
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


typedef RegisterArray<NUM_REGISTERS_PER_FRAME * FRAMES_PER_MSG> MessageRegisters;


// =================================================================
//                       PARAMS for SWTPG
// =================================================================
const uint16_t m_tpg_threshold = 5;                    // units of sigma 
const uint8_t m_tpg_tap_exponent = 6; 
const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
std::vector<int16_t> tpg_taps;                       // firwin_int(7, 0.1, multiplier)
int16_t* tpg_taps_p = nullptr;
uint16_t* primfind_dest = nullptr;  
std::unique_ptr<ProcessingInfo<NUM_REGISTERS_PER_FRAME>> tpg_processing_info;




void reset_TPG_objects(bool first_hit){
  if (first_hit){
    tpg_taps = firwin_int(7, 0.1, m_tpg_multiplier); //coefficients associated with nf 
    tpg_taps.push_back(0);
    if (tpg_taps_p == nullptr){
    tpg_taps_p = new int16_t[tpg_taps.size()];
    }
    if (primfind_dest == nullptr){
      primfind_dest = new uint16_t[100000];
    }
    for (size_t i = 0; i < tpg_taps.size(); ++i) {
    tpg_taps_p[i] = tpg_taps[i];
    }
        
    //Initialise the TPG processing info struct ready for next WIB Superchunk
    tpg_processing_info = std::make_unique<ProcessingInfo<NUM_REGISTERS_PER_FRAME>>(
										       nullptr, //input                                        
										       FRAMES_PER_MSG, //timeWindowNumFrames = 12
										       0, //first register in frame      
										       NUM_REGISTERS_PER_FRAME,
										       primfind_dest, // output                                
										       tpg_taps_p, // values of (nf) coeffictients /taps
										       (uint8_t)tpg_taps.size(), //n taps
										       m_tpg_tap_exponent, // multiplier(?)
										       m_tpg_threshold, //adc max(?)
										       0, //nhits
										       0 );//absTimeModNTAPS
  }
  tpg_processing_info->nhits = 0;     
}


void find_hits(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp, bool first_hit, bool dump_adc_data, bool dump_hit_data) {

    // Parse the DUNEWIB frames
    auto wf = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp);
    uint64_t timestamp = wf->get_timestamp();       
    
    // Frame expansion
    MessageRegisters registers_array;
    //expand_wib2_adcs(fp, &registers_array, register_selection); 


}  



int main()
{

    // =================================================================
    //                       SINGLE FRAME
    // =================================================================
    /*
    dunedaq::detdataformats::wib2::WIB2Frame frame;
    std::memset(&frame, 0, sizeof(dunedaq::detdataformats::wib2::WIB2Frame));
    for(int i=0; i<256; ++i){
        frame.set_adc(i, 0x3a0+i);
    }
    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};
    // Frame expansion
    RegisterArray<16> unpacked=unpack(frame);

    reset_TPG_objects(true);

    tpg_processing_info->setState(unpacked);
    tpg_processing_info->input = &unpacked;
    *tpg_processing_info->output = MAGIC;
    
    process_window_avx2_running_sum(*tpg_processing_info);
    process_window_avx2(*tpg_processing_info);
    process_window_naive(*tpg_processing_info);
    
    printf("\n");
    */

    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    const std::string input_file = "/nfs/sw/work_dirs/aabedabu/daq_config/frames_wib2.bin";
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, SUPERCHUNK_FRAME_SIZE);
 
    m_source_buffer->read(input_file);
    auto& source = m_source_buffer->get(); 
    int num_superchunks = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of superchunks in the input file: " << num_superchunks << std::endl;

    bool first_hit = true;
    int offset = 0;   
    
    //choose whether to save adc data and hit data to a file 
    bool dump_adc_data = false;  
    bool dump_hit_data = true;  

    // Loop over the DUNEWIB superchunks in the file
    while (offset<num_superchunks){
      // current superchunk
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + offset*SUPERCHUNK_FRAME_SIZE);
      
      //  **Find collection hits for the current fp
      reset_TPG_objects(first_hit);
      //find_collection_hits(fp, first_hit, dump_adc_data, dump_hit_data);
  
      first_hit = false;    
      ++offset;
    }


    std::cout << "Finished testing." << std::endl;




}


