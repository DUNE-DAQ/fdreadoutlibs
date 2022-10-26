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
#include "detchannelmaps/TPCChannelMap.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessNaive.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wib2/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wib2/tpg/DesignFIR.hpp"
#include "fdreadoutlibs/wib2/tpg/RegisterToChannelNumber.hpp"


// system
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


typedef swtpg_wib2::RegisterArray<swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::FRAMES_PER_MSG> MessageRegisters;
using timestamp_t = std::uint64_t; 

// =================================================================
//                       PARAMS for SWTPG
// =================================================================
const uint16_t m_tpg_threshold = 5;                    // units of sigma 
const uint8_t m_tpg_tap_exponent = 6; 
const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
std::vector<int16_t> tpg_taps;                       // firwin_int(7, 0.1, multiplier)
int16_t* tpg_taps_p = nullptr;
uint16_t* primfind_dest = nullptr;  
std::unique_ptr<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>> tpg_processing_info;

// Map from expanded AVX register position to offline channel number
swtpg_wib2::RegisterChannelMap register_channel_map; 
std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> channel_map;
std::string ch_map_name = "HDColdboxChannelMap";
//choose whether to save adc data and hit data to a file 
bool dump_adc_data = false;  
bool dump_hit_data = false;  


int link_no; 
int crate_no;
int slot_no;

void reset_TPG_objects(bool first_hit) {
  if (first_hit){
    tpg_taps = swtpg_wib2::firwin_int(7, 0.1, m_tpg_multiplier); //coefficients associated with nf 
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
    tpg_processing_info = std::make_unique<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>>(
										       nullptr, //input                                        
										       swtpg_wib2::FRAMES_PER_MSG, //timeWindowNumFrames = 12
										       0, //first register in frame      
										       swtpg_wib2::NUM_REGISTERS_PER_FRAME,
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



// Function save the TP data to a file 
void save_hit_data( triggeralgs::TriggerPrimitive trigprim ){
  std::ofstream out_file; 
  out_file.open("TP_dump.txt", std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  out_file << trigprim.channel << "," << trigprim.time_start << "," << trigprim.time_over_threshold << "," 
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak << "\n";  

  out_file.close();
}


// Function to save raw ADC data to a file (only for debugging) 
void save_raw_data(swtpg_wib2::MessageRegisters register_array, 
	       uint64_t t0, int set_index, int column_index)
{
  std::ofstream out_file;
  out_file.open("single_channel_data.txt", std::ofstream::app);

  uint64_t t_current= t0 ; 
  
  for ( size_t j = 0; j < swtpg_wib2::FRAMES_PER_MSG; j++){
    uint16_t adc = register_array.uint16(set_index*swtpg_wib2::FRAMES_PER_MSG + j, column_index);
    out_file << adc << "," << t_current << std::endl;
    t_current += 32;
  }
  out_file.close();
}




unsigned int process_swtpg_hits(uint16_t* primfind_it, timestamp_t timestamp)   {


    constexpr int clocksPerTPCTick = 32;

    uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; // NOLINT(build/unsigned)
    unsigned int nhits = 0;

    while (*primfind_it != swtpg_wib2::MAGIC) {
      // First, get all of the register values (including those with no hit) into local variables
      for (int i = 0; i < 16; ++i) {
        chan[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_end[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_charge[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_tover[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }

      // Now that we have all the register values in local
      // variables, loop over the register index (ie, channel) and
      // find the channels which actually had a hit, as indicated by
      // nonzero value of hit_charge
      for (int i = 0; i < 16; ++i) {
        //std::cout << hit_charge[i] << "  " << chan[i] << std::endl;
        if (hit_charge[i] && chan[i] != swtpg_wib2::MAGIC) {
          // This channel had a hit ending here, so we can create and output the hit here
          //const uint16_t offline_channel = register_channel_map.channel[chan[i]];
          //std::cout << chan[i] << "   "  << register_channel_map.channel[chan[i]] << std::endl;
          const uint16_t offline_channel = chan[i];

          uint64_t tp_t_begin =                                                        // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - hit_tover[i]);       // NOLINT(build/unsigned)
          uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);      // NOLINT(build/unsigned)

          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          // sed -n -e 's/.*Hit: \(.*\) \(.*\).*/\1 \2/p' log.txt  > hits.txt
          //
          //TLOG_DEBUG(0) << "Hit: " << tp_t_begin << " " << offline_channel;

          triggeralgs::TriggerPrimitive trigprim;
          trigprim.time_start = tp_t_begin;
          trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
          trigprim.time_over_threshold = hit_tover[i] * clocksPerTPCTick;
          trigprim.channel = offline_channel;
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_charge[i] / 20;
          trigprim.detid =
            link_no; // TODO: convert crate/slot/link to SourceID Roland Sipos rsipos@cern.ch July-22-2021
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
          trigprim.version = 1;

          if (dump_hit_data){
            save_hit_data(trigprim);
          }          

          //std::cout << "TP makes sense? -> hit_t_begin:" << tp_t_begin << " hit_t_end:" << tp_t_end
          //         << " time_peak:" << (tp_t_begin + tp_t_end) / 2 << std::endl;
          
          ++nhits;
        }
      }
    }
    return nhits;
  }



void find_hits(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp, bool first_hit) {

    // Parse the DUNEWIB frames
    auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp);
    uint64_t timestamp = wfptr->get_timestamp();       
    
    // Frame expansion
    MessageRegisters registers_array;
    expand_wib2_adcs(fp, &registers_array, 0); // register selector is zero because we are using 16 registers

    //Save ADC info? (only for debugging)
    if (dump_adc_data){
      //which register and column (channel) to save
      int i_register = 5;
      int i_column   = 3;
      save_raw_data(registers_array, timestamp, i_register, i_column);
    }


    
    if (first_hit) {

      channel_map = dunedaq::detchannelmaps::make_map(ch_map_name);
      register_channel_map = swtpg_wib2::get_register_to_offline_channel_map_wib2(wfptr, channel_map, 0);
      tpg_processing_info->setState(registers_array);

      // Debugging statements 
      link_no = wfptr->header.link;
      crate_no = wfptr->header.crate;
      slot_no = wfptr->header.slot;
      std::cout << "Got first item, link/crate/slot=" << link_no << "/" << crate_no << "/" << slot_no << std::endl;      

      std::stringstream ss;
      ss << " Channels are:\n";
      for(size_t i=0; i<swtpg_wib2::NUM_REGISTERS_PER_FRAME*swtpg_wib2::SAMPLES_PER_REGISTER; ++i){
        ss << i << "\t" << register_channel_map.channel[i] << "\n";
      }
      //std::cout << ss.str();      
      

    } // end if (first_hit)

    // Execute the SWTPG algorithm
    tpg_processing_info->input = &registers_array;
    *primfind_dest = swtpg_wib2::MAGIC;
    //swtpg_wib2::process_window_avx2(*tpg_processing_info);
    swtpg_wib2::process_window_naive(*tpg_processing_info);
    
    unsigned int nhits = process_swtpg_hits(primfind_dest, timestamp);
    if ( nhits > 0 ) {
      //std::cout << "Non null hits: " << nhits << " for ts: " << timestamp << std::endl;    
    }



    



}  




int main()
{
    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    const std::string input_file = "/nfs/sw/work_dirs/aabedabu/daq_config/frames_wib2.bin";
    std::unique_ptr<dunedaq::readoutlibs::FileSourceBuffer> m_source_buffer;
    m_source_buffer = std::make_unique<dunedaq::readoutlibs::FileSourceBuffer>(10485100, swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
 
    m_source_buffer->read(input_file);
    auto& source = m_source_buffer->get(); 
    int num_superchunks = m_source_buffer->num_elements(); // file_ size/chunk_size = 180 

    std::cout << "Number of superchunks in the input file: " << num_superchunks << std::endl;

    bool first_hit = true;
    int offset = 0;   
    

    // Loop over the DUNEWIB superchunks in the file
    while (offset<num_superchunks){
      // current superchunk
      auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter*>(source.data() + offset*swtpg_wib2::SUPERCHUNK_FRAME_SIZE);
      
      //  **Find collection hits for the current fp
      reset_TPG_objects(first_hit);
      find_hits(fp, first_hit);
  
      first_hit = false;    
      ++offset;
    }


    std::cout << "\n\nFinished testing." << std::endl;




}


