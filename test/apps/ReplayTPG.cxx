/**
 * @file ReplayTPG.cxx Main file for testing different tpg algorithms 
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
#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

#include "triggeralgs/TriggerPrimitive.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"

#include "hdf5libs/HDF5RawDataFile.hpp"

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


using namespace dunedaq::hdf5libs;


// =================================================================
//                       PUBLIC VARIABLES
// =================================================================

struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};

unsigned int total_hits = 0;
unsigned int total_hits_trigger_primitive = 0;
bool first_hit = true;

dunedaq::fdreadoutlibs::WIBEthFrameHandler fh;

std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;

std::string select_algorithm = "";
std::string select_implementation = "";
bool save_adc_data = false;
bool save_trigprim = false;

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
	   << trigprim.time_peak << "," << trigprim.adc_integral << ","  << trigprim.adc_peak << "\n";  

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
  for (auto ichan = 0; ichan < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++ichan) {
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
        out_file << "Time " << iframe << " channel " <<  ichan << " ADC_value " <<  adc_value <<  " timestamp " << t_current << std::endl;
        t_current += 32;
      } 

    }
  }
  out_file.close();


}

void
save_tp(const dunedaq::trgdataformats::TriggerPrimitive& prim, size_t offset = 0)
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
	     << prim.time_peak << "," << prim.adc_integral << ","  << prim.adc_peak << "," << prim.detid << "," << prim.type << "\n";  
  
  }
  out_file.close();
  
}


void print_tps(std::unique_ptr<dunedaq::daqdataformats::Fragment>&& frag, int TP_index,
               size_t offset = 0)
{
  size_t payload_size = frag->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader);
  size_t n_tps = payload_size / sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  std::cout << "Trigger Primitive number " << TP_index << " with SourceID[" << frag->get_element_id() << "] has " << n_tps << " TPs" << std::endl;
  total_hits_trigger_primitive = total_hits_trigger_primitive + n_tps ; 
  size_t remainder = payload_size % sizeof(dunedaq::trgdataformats::TriggerPrimitive);
  assert(remainder == 0);
  const dunedaq::trgdataformats::TriggerPrimitive* prim = reinterpret_cast<dunedaq::trgdataformats::TriggerPrimitive*>(frag->get_data());
  
  for (size_t i = 0; i < n_tps; ++i) {
    save_tp(*prim, offset);
    ++prim;
  }
}



// =================================================================
//                       TPG FUNCTIONS
// =================================================================
void extract_hits_naive(uint16_t* output_location, uint64_t timestamp) {

    constexpr int clocksPerTPCTick = 32;
    uint16_t chan, hit_end, hit_charge, hit_tover; 

    size_t i = 0;
    while (*output_location != swtpg_wibeth::MAGIC) {
      chan   = *output_location++;
      hit_end    = *output_location++;
      hit_charge  = *output_location++;
      hit_tover     = *output_location++;


      //if (hit_charge && chan != swtpg_wibeth::MAGIC) {
      //  std::cout << "Channel number: " << chan << std::endl;
      //  std::cout << "Hit charge: " << hit_charge << std::endl;
      //}

      
      i += 1;
      uint64_t tp_t_begin =                                                        
        timestamp + clocksPerTPCTick * (int64_t(hit_end ) - hit_tover );       
      uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end );      

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;

      trigprim.time_over_threshold = hit_tover  * clocksPerTPCTick;


      trigprim.channel = chan;
      trigprim.adc_integral = hit_charge ;
      trigprim.adc_peak = hit_charge  / 20;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_hit_data(trigprim, "NAIVE");
      }
      ++total_hits;

    }
}

void extract_hits_avx(uint16_t* output_location, uint64_t timestamp) {

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
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - hit_tover[i]);       // NOLINT(build/unsigned)
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
          trigprim.channel = chan[i];
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



void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  std::cout << "EXECUTING TPG" << std::endl;

  // Set CPU affinity of the TPG thread
  SetAffinityThread(0);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();      
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);

  if (first_hit) {                     
    fh.m_tpg_processing_info->setState(registers_array);
    first_hit = false;    
    // Save ADC info
    //if (save_adc_data){
    //  save_raw_data(registers_array, timestamp, -1, select_algorithm + "_" + select_implementation);
    //}
  } 
  save_raw_data(registers_array, timestamp, -1, select_algorithm + "_" + select_implementation); 
       
  fh.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = fh.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  fh.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*fh.m_tpg_processing_info);
       
  // Insert output of the AVX processing into the swtpg_output 
  //swtpg_output swtpg_processing_result = {destination_ptr, timestamp};
    
  if (select_implementation == "AVX") {
    extract_hits_avx(destination_ptr, timestamp);
  } else if(select_implementation == "NAIVE") {  
    extract_hits_naive(destination_ptr, timestamp);
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
    std::string file_path_input = "./test.hdf5";
    app.add_option("-f,--file_path_input", file_path_input, "Path to the input file");


    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");
  
    app.add_option("-i,--implementation", select_implementation, "TPG implementation (AVX / NAIVE)");

    int num_TR_to_read = -1;
    app.add_option("-n,--num_TR_to_read", num_TR_to_read, "Number of Trigger Records to read. Default: select all TRs.");

    int swtpg_threshold = 500;
    app.add_option("-t,--swtpg_threshold", swtpg_threshold, "Value of the TPG threshold");

    app.add_option("--save_adc_data", save_adc_data, "Save ADC data (true/false)");

    app.add_option("--save_trigprim", save_trigprim, "Save trigger primitive data (true/false)");


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
        m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_naive_RS<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
        std::cout << "Created an instance of the " << select_algorithm << " algorithm ( " << select_implementation << " )" << std::endl;
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
    //                       Setup the SWTPG
    // =================================================================

    fh.initialize(swtpg_threshold);


    
    // =================================================================
    //                       READ THE HDF5 FILE
    // =================================================================

    // open our file reading
    const std::string ifile_name = file_path_input;
    HDF5RawDataFile h5_raw_data_file(ifile_name);

    auto recorded_size = h5_raw_data_file.get_attribute<size_t>("recorded_size");
    auto record_type = h5_raw_data_file.get_record_type();

    auto run_number = h5_raw_data_file.get_attribute<unsigned int>("run_number");
    auto file_index = h5_raw_data_file.get_attribute<unsigned int>("file_index");
  
    auto creation_timestamp = h5_raw_data_file.get_attribute<std::string>("creation_timestamp");
    auto app_name = h5_raw_data_file.get_attribute<std::string>("application_name");

    uint32_t n_ch = dunedaq::fddetdataformats::WIBEthFrame::s_num_channels;
    uint32_t n_smpl = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;

    std::cout << "Run number: " << run_number << std::endl;
    std::cout << "Recorded size [bytes]: " << recorded_size << std::endl;
    std::cout << "Recorded type: " << record_type << std::endl;


    auto records = h5_raw_data_file.get_all_record_ids();

    if (records.empty()) {
      std::cout << "\n\nNO TRIGGER RECORDS FOUND" << std::endl;    
      return 0;
    }


    // Check if the selected number of frames is <= than the ones available in the input file
    /*
    num_trigger_records = records.size()

    if (num_trigger_records < num_TR_to_read) {
      std::cout << "\n**ERROR**: Select a valid number of Trigger Records that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_TR_to_read == -1) {
      num_TR_to_read = num_trigger_records;
    }
    */

    int record_idx_TR = 0;
    int record_idx_TP = 0;
    auto start_test = std::chrono::high_resolution_clock::now();  


    for (auto const& rid : records) {
      for(auto const& frag_dataset : h5_raw_data_file.get_fragment_dataset_paths(rid)) {
        auto frag_ptr = h5_raw_data_file.get_frag_ptr(frag_dataset);

        if (frag_ptr->get_fragment_type() == dunedaq::daqdataformats::FragmentType::kWIBEth) {

        //auto id = frag_ptr->get_element_id();
        //auto element_id = id.id;
        int num_frames =
          (frag_ptr->get_size() - sizeof(dunedaq::daqdataformats::FragmentHeader)) / sizeof(dunedaq::fddetdataformats::WIBEthFrame);                

        std::cout << "Trigger Record number " << record_idx_TR << " [ " << frag_ptr->get_element_id() << "] has "   << num_frames << " frames" << std::endl;
        ++record_idx_TR;

        // Calculate the total number of elements in the array
        size_t total_elements = static_cast<size_t>(n_ch) * n_smpl * num_frames;
        
        // Create a vector of unsigned 16-bit integers (uint16_t)
        std::vector<uint16_t> result(total_elements);
        
        // Access the data through a pointer
        uint16_t* ptr_res = result.data();        

        for (int i = 0; i < num_frames; ++i) {

          auto fr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(
            static_cast<char*>(frag_ptr->get_data()) + i * sizeof(dunedaq::fddetdataformats::WIBEthFrame)
          );

          // Unpack the ADC values
          for (size_t j=0; j<n_smpl; ++j){
            for (size_t k=0; k<n_ch; ++k){
              ptr_res[(n_smpl*n_ch) * i + n_ch*j + k] = fr->get_adc(k, j);             
            } // loop over channels
          } // loop over time samples    

          if (record_idx_TR == 25) {
            for (size_t j=0; j<n_smpl; ++j){
              for (size_t k=0; k<n_ch; ++k){
                std::cout << "Channel " << k << " time sample " << j << " ADC " <<  fr->get_adc(k, j)  << std::endl;
                
              } // loop over channels
            } // loop over time samples            
          }
          
          auto fp = reinterpret_cast<dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter*>(result.data());
          if (record_idx_TR == 25) {
            execute_tpg(fp);
          }
          






        } // end loop over number of frames      
        
        } 
        if (frag_ptr->get_fragment_type() == dunedaq::daqdataformats::FragmentType::kTriggerPrimitive) {
          if (record_idx_TP == 25) {
            print_tps(std::move(frag_ptr), record_idx_TP);
          }
          ++record_idx_TP;
        }
        

        
      } // loop over all the fragments in a single trigger record    
    } // loop over all trigger records

    // Calculate elapsed time in seconds  
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_test).count();  
    std::cout << "Elapsed time for reading input file [ms]: " << elapsed_milliseconds << std::endl;      



    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << total_hits << " hits." << std::endl;
    std::cout << "Found in total  (from Trigger Primitive objects) " << total_hits_trigger_primitive << " TPs." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;

}


