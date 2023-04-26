/**
 * @file SwtpgNaive.hpp  class that performs the naive implementation of SWTPG
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TEST_SWTPG_NAIVE_HPPP_
#define TEST_SWTPG_NAIVE_HPPP_


// DUNE-DAQ
#include "detdataformats/wib2/WIB2Frame.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "detchannelmaps/TPCChannelMap.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessNaive.hpp"

#include "fdreadoutlibs/wib2/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wib2/tpg/DesignFIR.hpp"
#include "fdreadoutlibs/wib2/tpg/RegisterToChannelNumber.hpp"


// Local
#include "SwtpgBase.hpp"
#include "utils.hpp"


// system
#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>

#include <immintrin.h>
#include <cstdio> 
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>
#include <memory>



class SwtpgNaive : public SwtpgBase {


public: 
  explicit SwtpgNaive(bool dump_adc_data, bool dump_hit_data) {
    m_dump_adc_data = dump_adc_data;
    m_dump_hit_data = dump_hit_data;
  }
  SwtpgNaive(const SwtpgNaive&) = delete;
  SwtpgNaive& operator=(const SwtpgNaive&) = delete;

  void reset(bool first_hit, int threshold_value) {
    if (first_hit){
      m_tpg_taps = swtpg_wib2::firwin_int(7, 0.1, m_tpg_multiplier); //coefficients associated with nf 
      m_tpg_taps.push_back(0);
      if (m_tpg_taps_p == nullptr){
        m_tpg_taps_p = new int16_t[m_tpg_taps.size()];
      }
      if (m_primfind_dest == nullptr){
        m_primfind_dest = new uint16_t[100000];
      }
      for (size_t i = 0; i < m_tpg_taps.size(); ++i) {
        m_tpg_taps_p[i] = m_tpg_taps[i];
      }
    
    m_tpg_threshold = threshold_value;  
        
    //Initialise the TPG processing info struct ready for next WIB Superchunk
    m_tpg_processing_info = std::make_unique<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>>(
										       nullptr, //input                                        
										       swtpg_wib2::FRAMES_PER_MSG, //timeWindowNumFrames = 12
										       0, //first register in frame      
										       swtpg_wib2::NUM_REGISTERS_PER_FRAME,
										       m_primfind_dest, // output                                
										       m_tpg_taps_p, // values of (nf) coeffictients /taps
										       (uint8_t)m_tpg_taps.size(), //n taps
										       m_tpg_tap_exponent, // multiplier(?)
										       m_tpg_threshold, //adc max(?)
										       0, //nhits
										       0 );//absTimeModNTAPS
  }

  // If not first hit then reset hits count to zero
  m_tpg_processing_info->nhits = 0;     
  }




unsigned int extract_swtpg_hits_naive(uint16_t* primfind_it, timestamp_t timestamp)   {


    constexpr int clocksPerTPCTick = 32;
    //uint16_t chan[100], hit_end[100], hit_charge[100], hit_tover[100]; 
    uint16_t chan, hit_end, hit_charge, hit_tover, hit_peak_adc, hit_peak_time; 
    unsigned int nhits = 0;

    size_t i = 0;
    while (*primfind_it != swtpg_wib2::MAGIC) {
      chan   = *primfind_it++;
      hit_end    = *primfind_it++;
      hit_charge  = *primfind_it++;
      hit_tover     = *primfind_it++;
      hit_peak_adc  = *primfind_it++;
      hit_peak_time = *primfind_it++;

      i += 1;
      const uint16_t offline_channel = m_register_channel_map.channel[chan ];
      uint64_t tp_t_begin =                                                        
        timestamp + clocksPerTPCTick * (int64_t(hit_end ) - hit_tover );       
      uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end );      
      uint64_t tp_t_peak = 
	timestamp + clocksPerTPCTick * int64_t(hit_peak_time);

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = tp_t_peak;

      trigprim.time_over_threshold = hit_tover  * clocksPerTPCTick;


      trigprim.channel = offline_channel;
      trigprim.adc_integral = hit_charge;
      trigprim.adc_peak = hit_peak_adc;
      trigprim.detid =
        m_link_no; // TODO: convert crate/slot/link to SourceID Roland Sipos rsipos@cern.ch July-22-2021
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (m_dump_hit_data) {
        save_hit_data(trigprim, "SwtpgNaive");
      }
      ++nhits;
      

    }

    return nhits;
  }






void find_hits(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp, bool first_hit) {

    // Parse the DUNEWIB frames
    auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp);
    uint64_t timestamp = wfptr->get_timestamp();       


    // Print ADCs before expansion    
    /*
    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14}; 
    for (size_t iframe = 0; iframe < swtpg_wib2::FRAMES_PER_MSG; ++iframe) {
    const dunedaq::detdataformats::wib2::WIB2Frame* frame =
      reinterpret_cast<const dunedaq::detdataformats::wib2::WIB2Frame*>(fp) + iframe; // NOLINT
      for(int i=0; i<256; ++i){        
          int in_index=16*(i/16)+indices[i%16];
          uint16_t in_val=frame->get_adc(in_index);
          std::cout << "Index " << i << "  " <<  in_val  <<  std::endl;        
      }
    }  
    */


    // Frame expansion
    MessageRegisters registers_array;
    expand_wib2_adcs(fp, &registers_array, 0); // register selector is zero because we are using 16 registers

    //Save ADC info? (only for debugging)
    if (m_dump_adc_data){
      save_raw_data(registers_array, timestamp, -1, "SwtpgNaive");
    }


    
    if (first_hit) {

      // Save only one hit data
      //if (m_dump_adc_data){
      //  save_raw_data(registers_array, timestamp, -1, "SwtpgNaive");
      //}


      m_channel_map = dunedaq::detchannelmaps::make_map(m_ch_map_name);
      m_register_channel_map = swtpg_wib2::get_register_to_offline_channel_map_wib2(wfptr, m_channel_map, 0);
      m_tpg_processing_info->setState(registers_array);

      // Debugging statements 
      m_link_no = wfptr->header.link;
      m_crate_no = wfptr->header.crate;
      m_slot_no = wfptr->header.slot;
      std::cout << "Got first item, link/crate/slot=" << m_link_no << "/" << m_crate_no << "/" << m_slot_no << std::endl;      

      std::stringstream ss;
      ss << " Channels are:\n";
      for(size_t i=0; i<swtpg_wib2::NUM_REGISTERS_PER_FRAME*swtpg_wib2::SAMPLES_PER_REGISTER; ++i){
        ss << i << "\t" << m_register_channel_map.channel[i] << "\n";
      }
      //std::cout << ss.str();      
      

    } // end if (first_hit)

    // Execute the SWTPG algorithm
    m_tpg_processing_info->input = &registers_array;
    *m_primfind_dest = swtpg_wib2::MAGIC;
    swtpg_wib2::process_window_naive(*m_tpg_processing_info);

    std::cout << "Finished processing window " << std::endl;
    
    unsigned int nhits = extract_swtpg_hits_naive(m_primfind_dest, timestamp);

    if ( nhits > 0 ) {
      std::cout << "Non null hits: " << nhits << " for ts: " << timestamp << std::endl;    
    }
    
    m_total_swtpg_hits += nhits; 
    
  
    
}  


private: 

  std::unique_ptr<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>> m_tpg_processing_info;

  bool m_dump_adc_data;
  bool m_dump_hit_data;

};

#endif // TEST_SWTPG_NAIVE_HPPP_
