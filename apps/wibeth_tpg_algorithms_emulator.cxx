/**
 * @file ReplayTPG.cxx 
 * Replay tpg algorithms on trigger records
 * 
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

// DUNE-DAQ
#include "fddetdataformats/WIBEthFrame.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"
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

// Local
#include "tpg_apps_utilities.hpp"
#include "tpg_apps_issues.hpp"

using namespace dunedaq::hdf5libs;
using namespace tpgtools;

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

// =================================================================
//                       PUBLIC VARIABLES
// =================================================================

unsigned int total_hits = 0;
unsigned int total_hits_trigger_primitive = 0;
bool first_hit = true;

dunedaq::fdreadoutlibs::WIBEthFrameHandler fh;

std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;
std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> channel_map;

// Map from expanded AVX register position to offline channel number
swtpg_wibeth::RegisterChannelMap register_channel_map; 

// Mapping from expanded AVX register position to offline channel number
std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER> m_register_channels = {};

std::string select_algorithm = "";
std::string select_channel_map = "None";
bool save_adc_data = false;
bool save_trigprim = false;



// =================================================================
//                       TPG 
// =================================================================


void execute_tpg(const dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter* fp) {

  // Set CPU affinity of the TPG thread
  SetAffinityThread(0);

  // Parse the WIBEth frames
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp);
  uint64_t timestamp = wfptr->get_timestamp();     

  // Frame expansion
  swtpg_wibeth::MessageRegisters registers_array;
  swtpg_wibeth::expand_wibeth_adcs(fp, &registers_array);

  
  if (first_hit) {                         
    fh.m_tpg_processing_info->setState(registers_array);
    first_hit = false;    
    // Save ADC info
    if (save_adc_data){
      save_raw_data(registers_array, timestamp, -1, select_algorithm );
    }

  }


  fh.m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = fh.get_hits_dest();
  *destination_ptr = swtpg_wibeth::MAGIC;
  fh.m_tpg_processing_info->output = destination_ptr;
  m_assigned_tpg_algorithm_function(*fh.m_tpg_processing_info);
       

  // Parse the output from the TPG    
  extract_hits_avx(fh.m_tpg_processing_info->output, timestamp, m_register_channels, total_hits, save_trigprim);

}


// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "Test TPG algorithms" };

    // Set default input frame file
    std::string file_path_input = "";
    app.add_option("-f,--file-path-input", file_path_input, "Path to the input file");


    app.add_option("-a,--algorithm", select_algorithm, "TPG Algorithm (SimpleThreshold / AbsRS)");
  
    app.add_option("-m,--channel-map", select_channel_map, "Select a valid channel map: None, VDColdboxChannelMap, ProtoDUNESP1ChannelMap, PD2HDChannelMap, HDColdboxChannelMap, FiftyLChannelMap, etc.");


    int num_TR_to_read = -1;
    app.add_option("-n,--num-TR-to-read", num_TR_to_read, "Number of Trigger Records to read. Default: select all TRs. ");

    int tpg_threshold = 500;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold");

    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data");

    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");


    CLI11_PARSE(app, argc, argv);

    // =================================================================
    //                       Setup the SWTPG
    // =================================================================

    tpg_emulator emulator(tpg_threshold, save_adc_data, save_trigprim, select_algorithm, select_channel_map) ;
    emulator.initialize();

    
    // =================================================================
    //                       READ THE HDF5 FILE
    // =================================================================

    // open our file reading
    const std::string ifile_name = file_path_input;
    if (ifile_name.empty()) {
      throw std::runtime_error("Please select a valid input file.");
    } 
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
      TLOG() << "*** NO TRIGGER RECORDS FOUND" ;    
      return 0;
    }

    int record_idx_TR = 0;
    int record_idx_TP = 0;

    uint32_t element_id; // part of the source id

    // Measure the time taken to read the whole trigger record file
    auto start_test = std::chrono::high_resolution_clock::now();  

   
    for (auto const& rid : records) {
      for(auto const& frag_dataset : h5_raw_data_file.get_fragment_dataset_paths(rid)) {
        auto frag_ptr = h5_raw_data_file.get_frag_ptr(frag_dataset);

        if (record_idx_TR <= num_TR_to_read || num_TR_to_read == -1 ) {  

          //process_fragment(std::move(frag_ptr));   
          emulator.process_fragment(std::move(frag_ptr), record_idx_TR);
        
        } // if statement number of trigger records
        

        
      } // loop over all the fragments in a single trigger record    
    } // loop over all trigger records

    // Calculate elapsed time in seconds 
    // AAA: to be remove it if not useufl?
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_test).count();  



    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Elapsed time for reading input file [ms]: " << elapsed_milliseconds;
    std::cout << "Found in total " << emulator.get_total_hit_number() << " hits" << std::endl;
    //std::cout << "Found in total  (from Trigger Primitive objects) " << total_hits_trigger_primitive << " TPs" << std::endl;
    


}


