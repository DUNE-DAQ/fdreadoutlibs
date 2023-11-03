/**
 * @file tpg_apps_utilities.hpp 
 * Utility functions needed for standalone
 * TPG applications
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
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
#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"
#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"
#include "detchannelmaps/TPCChannelMap.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "hdf5libs/HDF5RawDataFile.hpp"
#include "logging/Logging.hpp"


using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;


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
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak <<  ","  << trigprim.type << "\n";  

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

  std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};
  for (auto ichan = 0; ichan < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++ichan) {

    int in_index=16*(ichan/16)+indices[ichan%16];

    const size_t register_index = ichan / swtpg_wibeth::SAMPLES_PER_REGISTER;
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
    
        //int16_t adc_value = input16[index];    
        int16_t adc_value = register_array.uint16(index);
        out_file << " Time " << iframe << " channel " <<  ichan << " ADC_value " <<  adc_value <<  " timestamp " << t_current << std::endl;

        t_current += 32; 
      } 

    }
  }
  out_file.close();


}

void
save_tp(const dunedaq::trgdataformats::TriggerPrimitive& prim, bool save_trigprim)
{
  std::ofstream out_file; 
  std::ostringstream oss;

  if (save_trigprim) {   
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);  
    
    oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
    auto date_time_str = oss.str();
  
    std::string file_name = "TriggerPrimitiveRecord_dump_" + date_time_str + ".txt";
    out_file.open(file_name.c_str(), std::ofstream::app);
  
    //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
    out_file << prim.channel << "," << prim.time_start << "," << prim.time_over_threshold << "," 
	     << prim.time_peak << "," << prim.adc_integral << ","  << prim.adc_peak << "," << prim.type << "\n";  
  
  }
  out_file.close();
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Saved TriggerPrimitives output file." ;
  
}


void print_tps(std::unique_ptr<dunedaq::daqdataformats::Fragment>&& frag, 
               int TP_index, 
               unsigned int& total_tp_hits, 
               bool save_trigprim)
{
  size_t payload_size = frag->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader);
  size_t n_tps = payload_size / sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Trigger Primitive number " << TP_index << " with SourceID[" << frag->get_element_id() << "] has " << n_tps << " TPs";
  total_tp_hits = total_tp_hits + n_tps ; 
  size_t remainder = payload_size % sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  assert(remainder == 0);
  const dunedaq::trgdataformats::TriggerPrimitive* prim = reinterpret_cast<dunedaq::trgdataformats::TriggerPrimitive*>(frag->get_data());
  
  for (size_t i = 0; i < n_tps; ++i) {
    save_tp(*prim, save_trigprim);
    ++prim;
  }
}



// =================================================================
//                       TPG FUNCTIONS
// =================================================================

void extract_hits_avx(uint16_t* output_location, uint64_t timestamp,
                      std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER>& register_channels,
                      unsigned int& total_hits,
                      bool save_trigprim) {

  constexpr int clocksPerTPCTick = 32;
  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; 

  while (*output_location != swtpg_wibeth::MAGIC) {
    for (int i = 0; i < 16; ++i) {
      chan[i] = *output_location++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_end[i] = *output_location++; 
    }
    for (int i = 0; i < 16; ++i) {
      hit_charge[i] = *output_location++;
    }
    for (int i = 0; i < 16; ++i) {        
      hit_tover[i] = *output_location++; 
    }  
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (int i = 0; i < 16; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wibeth::MAGIC) {
        //std::cout << "Channel number: " << chan[i] << std::endl;
        //std::cout << "Hit charge: " << hit_charge[i] << std::endl;


        uint64_t tp_t_begin =                                                        // NOLINT(build/unsigned)
          timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_tover[i]));       // NOLINT(build/unsigned)
        uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);      // NOLINT(build/unsigned)
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
        trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
        trigprim.time_over_threshold = hit_tover[i] * clocksPerTPCTick;      
        trigprim.channel = register_channels[chan[i]]; //offline channel map
        trigprim.adc_integral = hit_charge[i];
        trigprim.adc_peak = hit_charge[i] / 20;
        trigprim.detid = 666;          
        trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
        trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
        trigprim.version = 1;
        if (save_trigprim){
          save_hit_data(trigprim, "AVX");
        }          

        ++total_hits;
      }
    } // loop over 16 registers 
  } // while not magic   

}

