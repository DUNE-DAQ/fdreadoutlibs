/**
 * @file wibeth_tpg_pattern_generator.cxx Main file for generating tpg patterns 
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
//#include "fdreadoutlibs/wibeth/tpg/ProcessNaiveRS.hpp"


#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"

#include "triggeralgs/TriggerPrimitive.hpp"


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

// =================================================================
//                     PATTERN GENERATION INFO
// =================================================================
struct PattgenInfo
{
  PattgenInfo(dunedaq::fddetdataformats::WIBEthFrame* output_frame_
		, int time_tick_offset_
		, int input_ch_
		, int iframe_
		)
  : output_frame(output_frame_)
  , time_tick_offset(time_tick_offset_)
  , input_ch(input_ch)
  , iframe(iframe_)
  {}
  dunedaq::fddetdataformats::WIBEthFrame* output_frame;  
  int time_tick_offset;
  int input_ch;
  int iframe;
};
class PattgenHandler {
  public:
  PattgenHandler(){};
  ~PattgenHandler(){};
  std::unique_ptr<PattgenInfo> m_pattgen_info;

  void initialize() {
    m_pattgen_info = std::make_unique<PattgenInfo>(nullptr
		    , 0
		    , 0
		    , 0
		    );
  }

};



// =================================================================
//                       PUBLIC VARIABLES
// =================================================================

struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};

int WIBEth_FRAME_SIZE = dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter::fixed_payload_size;

int duration_test = 120; // default value
unsigned int total_hits = 0;
uint64_t first_timestamp = 0;

dunedaq::fdreadoutlibs::WIBEthFrameHandler fh;
PattgenHandler ph;

std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;
std::function<void(PattgenInfo& info)> m_assigned_pattgen_function;


std::string select_pattern = "";
//std::string select_implementation = "";
bool save_adc_data = false;
bool save_trigprim = false;
std::string out_suffix = "";

// =================================================================
//                     PATTERN GENERATION FUNCTIONS
// =================================================================
void
pattgen_function_golden(PattgenInfo& info) {
  std::cout << "********** GENERATED PATTERN: PATT_GOLDEN " << std::endl; 
  int patt_time[9]{0};
  int patt_adc[9]{500, 502, 504, 505, 506, 505, 504, 502, 500};
  patt_time[0] = info.time_tick_offset;
  std::cout << "DBG pattgen info " << info.time_tick_offset << std::endl; 
  int npatt = sizeof(patt_time) / sizeof(*patt_time);
  for (int ipatt=1; ipatt<npatt; ipatt++) {
    patt_time[ipatt] = info.time_tick_offset+ipatt < 64 ? info.time_tick_offset+ipatt : info.time_tick_offset-64+ipatt;
  } 

  for (int itime=0; itime<64; ++itime) {
    for (int ch=0; ch<64; ++ch) {
      info.output_frame->set_adc(ch, itime, 0); // set pedestal value, make it configurable
    }
    //if (i==0 && itime < npatt && std::find(patt_time, patt_time + npatt, itime) == patt_time + npatt ) {
    if (info.iframe==0 && itime < patt_time[0]) {
      std::cout << "Nothing to do for first frame" << std::endl;
    } else {
      for (int ipatt=0; ipatt<sizeof(patt_time)/sizeof(patt_time[0]); ipatt++) {
        if (itime == patt_time[ipatt]) info.output_frame->set_adc(info.input_ch, itime, patt_adc[ipatt]);
      }
    }
    uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
    //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << info.input_ch << " \t\tTimeSample: " << itime <<  std::endl;
  }
}

void
pattgen_function_pulse(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: PULSE " << std::endl;
    for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          info.output_frame->set_adc(ch, itime, 0);
        }
        if (itime == 0 && info.iframe==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {
          info.output_frame->set_adc(info.input_ch, itime, 666);
        }
        uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
        //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}
void
pattgen_function_edge_square(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: EDGE_SQUARE " << std::endl;
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        std::cout << "Nothing to do for first frame" << std::endl;
      } else {	
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 500);      
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 500);      
      }
      uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
      //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}
void
pattgen_function_edge_left(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: EDGE_LEFT " << std::endl;
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        std::cout << "Nothing to do for first frame" << std::endl;
      } else {	
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 500);      
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 501);      
      }
      uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
      //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}
void
pattgen_function_edge_right(PattgenInfo& info) {
    std::cout << "********** GENERATED PATTERN: EDGE_RIGHT " << std::endl;
    for (int itime=0; itime<64; ++itime) {
      for (int ch=0; ch<64; ++ch) {
        info.output_frame->set_adc(ch, itime, 0);
      }
      if (itime >= 0 && itime<=62 && info.iframe==0) {
        std::cout << "Nothing to do for first frame" << std::endl;
      } else {	
        if (itime == 0) info.output_frame->set_adc(info.input_ch, itime, 501);      
        if (itime == 63) info.output_frame->set_adc(info.input_ch, itime, 500);      
      }
      uint16_t adc_val = info.output_frame->get_adc(info.input_ch, itime);
      //std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
    }
}

// =================================================================
//                       FUNCTIONS and UTILITIES
// =================================================================

// Function save the TP data to a file 
void save_hit_data( triggeralgs::TriggerPrimitive trigprim, std::string source_name, std::string out_suffix=""){
  std::ofstream out_file; 

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);  
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M");
  auto date_time_str = oss.str();

  std::string file_name = "TP_dump_" + source_name + "_" + date_time_str + out_suffix + ".txt";
  out_file.open(file_name.c_str(), std::ofstream::app);

  //offline channel, start time, time over threshold [ns], peak_time, ADC sum, amplitude    
  /*out_file << trigprim.time_start 
           << " " << trigprim.time_over_threshold 
           << " " << trigprim.time_peak 
           << " " << trigprim.channel 
           << " " << trigprim.adc_integral 
           << " "  << trigprim.adc_peak << "\n";*/
  // DBG temporary for the sake of comparison 
  out_file << trigprim.channel
	   << "," << trigprim.time_start 
           << "," << trigprim.time_over_threshold 
           << "," << trigprim.time_peak 
           << "," << trigprim.adc_integral 
           << ","  << trigprim.adc_peak << "\n";  


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
  for (int ichan = 0; ichan < static_cast<int>(swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER); ++ichan) {
    const size_t register_index = ichan / swtpg_wibeth::SAMPLES_PER_REGISTER;
    if (register_index >= swtpg_wibeth::NUM_REGISTERS_PER_FRAME)
       continue;

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
    
        int16_t adc_value = input16[index];
        std::cout << adc_value << std::endl;
        out_file << ichan << "," <<  adc_value << "," << t_current << std::endl;
        t_current += 32;
      } 

    }
  }
  out_file.close();
}

class FrameFile
{
public:

    FrameFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[sizeof(dunedaq::fddetdataformats::WIBEthFrame)])
    {
        if(m_file.bad() || m_file.fail() || !m_file.is_open()){
            throw std::runtime_error(std::string("Bad file ")+std::string(filename));
        }
        // Calculate the length of the file
        m_file.seekg(0, m_file.end);
        m_length = m_file.tellg();
        m_file.seekg(0, m_file.beg);
        if(m_length==0){
            throw std::runtime_error("Empty file");
        }
        //if(m_length%sizeof(dunedaq::fddetdataformats::WIBEthFrame)!=0){
        //    throw std::runtime_error("File does not contain an integer number of frames");
        //}
        m_n_frames=m_length/sizeof(dunedaq::fddetdataformats::WIBEthFrame);

	// Reinterpret the frame as WIBEthFrame
	
	m_file.read(m_buffer, m_file.eof());
        m_wibeth_frame = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }

    ~FrameFile()
    {
        m_file.close();
        delete[] m_buffer;
    }

    // Length of the file in bytes
    size_t length() const {return m_length;}
    // Number of frames in the file
    size_t num_frames() const { return m_n_frames; }

    dunedaq::fddetdataformats::WIBEthFrame* get_wibeth_frame() const { return m_wibeth_frame; }

    dunedaq::fddetdataformats::WIBEthFrame* frame(size_t i)
    {
        if(i>=num_frames()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the fragment into the buffer
        m_file.read(m_buffer,sizeof(dunedaq::fddetdataformats::WIBEthFrame));
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(m_buffer);
    }
    
protected:
    std::ifstream m_file;
    char* m_buffer;
    dunedaq::fddetdataformats::WIBEthFrame* m_wibeth_frame = nullptr;
    size_t m_length;
    size_t m_n_frames;
};

void hit_finder(std::vector<uint16_t>& adcs, std::vector<std::vector<int>>& out, const int& channel, 
		const int& threshold, const uint64_t timestamp) {

  //int threshold = (int)m_threshold;

      // index of ADCs above threshold
      //int threshold = (int)m_threshold;
      std::cout << "DBG threshold " << threshold << std::endl;
      std::cout << "DBG ADC sample length " << adcs.size() << std::endl;
      for(long unsigned int m = 0; m < adcs.size(); m++){
        std::cout << "DBG ADC sample, value " << m << ", " << adcs[m] << std::endl;
      }

      int m_tov_min = 2;

      bool dbg = false;
      int sz = adcs.size();
      if (sz == 64) {
        for (auto& t : adcs) std::cout << "DBG pattgen 64 adcs " << t << std::endl;
      }

      std::vector<int> igt;
      std::vector<uint16_t>::iterator it_adcs = adcs.begin();
      while ((it_adcs = find_if(it_adcs, adcs.end(), [&threshold](int x){return x > threshold; })) != adcs.end())
      {
        igt.push_back(distance(adcs.begin(), it_adcs));
        // cout << distance(adcs.begin(), it_adcs) << endl;
        it_adcs++;
      }

      if (dbg && sz == 64) std::cout << "DBG pattgen 64 igt " << igt.size() << std::endl;

      // cout << "igt.size = " << igt.size() << endl;
      // cout << (igt.size() < m_tov_min) << endl;

      if(igt.size() < m_tov_min){
        // cout << "m_tov_min condition not fulfilled for whole packet!" << endl;
        return void();
      }

      std::vector<int> igt_diff;
      adjacent_difference (igt.begin(), igt.end(), back_inserter(igt_diff));
      igt_diff.erase(igt_diff.begin());

      // find start and end of hits
      std::vector<int> istart;
      std::vector<int> iend;
      istart.push_back(0);
      std::vector<int>::iterator it_igt = igt_diff.begin();
      while ((it_igt = find_if(it_igt, igt_diff.end(), [ ](int x){return x != 1; })) != igt_diff.end())
      {
        istart.push_back(distance(igt_diff.begin(), it_igt)+1);
        iend.push_back(distance(igt_diff.begin(), it_igt));
        it_igt++;
      }
      iend.push_back(igt.size()-1);

      std::vector<int> start;
      std::vector<int> end;
      std::vector<int> hitcontinue;
      for(long unsigned int i = 0; i < istart.size(); i++){
        start.push_back(igt[istart[i]]);
        end.push_back(igt[iend[i]]);
        if(end[i] == 63){
          hitcontinue.push_back(1);
        }else{
          hitcontinue.push_back(0);
        }
      }

      // find hit sums
      std::vector<int> sums;
      for(long unsigned int j = 0; j < start.size(); j++){
	std::vector<uint16_t>::iterator it_start = adcs.begin()+start[j];
	std::vector<uint16_t>::iterator it_end = adcs.begin()+end[j]+1;
        int sum = accumulate(it_start, it_end, 0);
        sums.push_back(sum);
      }

      // find peak adcs and times
      std::vector<int> peak_adcs;
      std::vector<int> peak_times;
      for(long unsigned int j = 0; j < start.size(); j++){
	std::vector<uint16_t>::iterator it_start = adcs.begin()+start[j];
	std::vector<uint16_t>::iterator it_end = adcs.begin()+end[j]+1;
	std::vector<uint16_t>::iterator max = max_element(it_start, it_end);
        int peak = *max;
        int time = distance(adcs.begin(), max);
        peak_adcs.push_back(peak);
        peak_times.push_back(time);
      }

      if (dbg && sz == 64) {
	std::cout << "DBG pattgen 64 hit sizes " << channel << " : " << start.size() << " " << end.size() << " " << peak_times.size() << " " << sums.size() << " " << peak_adcs.size() << " " << hitcontinue.size() << std::endl;
      }

      // check output hits fullfil the m_tov_min condition
      for(long unsigned int k = 0; k < start.size(); k++){
        if (peak_adcs[k] > 16384) {
          if (dbg && sz == 64) std::cout << "DBG pattgen 64 peak adc is too big " << peak_adcs[k] << std::endl;
          continue; // 2**14
        }
        if (end[k]-start[k] >= m_tov_min-1){
          //std::vector<int> aux = {start[k], end[k], peak_times[k], peak_adcs[k], sums[k], hitcontinue[k]};
          //std::vector<int> aux = {k, start[k], end[k], peak_times[k], peak_adcs[k], sums[k], hitcontinue[k]};
          // start, end, peak_time, channel, sum_adc, peak_adc, hitcotninue
	  std::vector<int> aux = {start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]};
          out.push_back(aux);

	  std::cout << "# hit " << "[" << adcs.size() << "]"<< k << ", start time " << start[k] << ", end time " << end[k] << ", peak time " << peak_times[k] << ", peak adc " << peak_adcs[k] << ", sum adc " << sums[k] << ", hit continue " << hitcontinue[k] << ", timestamp " << timestamp << std::endl;
          
	  uint64_t ts_tov = (end[k]-start[k])*32;
	  uint64_t ts_start = start[k] * 32 + timestamp;
	  uint64_t ts_peak = peak_times[k] * 32 + timestamp;
	  std::cout << channel << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << sums[k] << "," << peak_adcs[k] << std::endl;
	}
      }
}



// =================================================================
//                       TPG FUNCTIONS
// =================================================================
void extract_hits_naive(uint16_t* output_location, uint64_t timestamp, std::string out_suffix) {

    constexpr int clocksPerTPCTick = 32;
    const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
    uint16_t chan, hit_end, hit_peak_adc, hit_charge, hit_tover, hit_peak_time;

    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    //std::cout << "DBG ProcessingInfo: timeWindowNumFrames " << fh.m_tpg_processing_info->timeWindowNumFrames << std::endl;
    //std::cout << "DBG nhits " << fh.m_tpg_processing_info->nhits << std::endl;

    size_t i = 0;
    //while (*output_location != swtpg_wibeth::MAGIC) {
    for (size_t n=0; n<fh.m_tpg_processing_info->nhits; n++) {
      chan            = *output_location++;     // TYPES!!!  int32_t  would be better match to triggeralgs::TriggerPrimitive
      hit_end         = *output_location++;  // TYPES!!!  int16_t  better than  uint16_t
      hit_charge      = *output_location++;     // TYPES!!!  uint32_t  would be best
      hit_tover       = *output_location++;     // TYPES!!!  uint32_t  would be best
      hit_peak_adc    = *output_location++;
      hit_peak_time   = *output_location++;     // TYPES!!!  uint32_t  would be best

      i += 1;
      chan = nreg*(chan/nreg)+indices[chan%nreg];
     
      //std::cout << "DBG chan " << n << ": " << chan << std::endl;
      //std::cout << "DBG hit_end " << n << ": " << hit_end << std::endl;

      uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end - (int64_t)hit_tover);
      uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time;

      triggeralgs::TriggerPrimitive trigprim;
      trigprim.time_start = tp_t_begin;
      trigprim.time_peak = tp_t_peak;

      trigprim.time_over_threshold = uint64_t((hit_tover - 1) * clocksPerTPCTick);

      trigprim.channel = chan;
      trigprim.adc_integral = hit_charge;
      trigprim.adc_peak = hit_peak_adc;
      trigprim.detid = 666; 
      trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
      trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
      trigprim.version = 1;
    
      if (save_trigprim) {
        save_hit_data(trigprim, "NAIVE", out_suffix);
      }
      ++total_hits;

    }
}


void extract_hits_avx(uint16_t* output_location, uint64_t timestamp, std::string out_suffix) {

  constexpr int clocksPerTPCTick = 32;
  const constexpr std::size_t nreg = swtpg_wibeth::SAMPLES_PER_REGISTER;
  uint16_t chan[nreg], hit_end[nreg], hit_peak_adc[nreg], hit_charge[nreg], hit_tover[nreg], hit_peak_time[nreg];

  //std::cout << "DBG nhits " << fh.m_tpg_processing_info->nhits << std::endl;

  for (size_t n=0; n<fh.m_tpg_processing_info->nhits; n++) {
    for (std::size_t i = 0; i < nreg; ++i) {
      chan[i] = *output_location++; 
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_end[i] = *output_location++; 
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_charge[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {        
      hit_tover[i] = *output_location++; 
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_adc[i] = *output_location++;
    }
    for (std::size_t i = 0; i < nreg; ++i) {
      hit_peak_time[i] = *output_location++;
    }
    
    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (std::size_t i = 0; i < nreg; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wibeth::MAGIC) {
        //std::cout << "Channel number: " << chan[i] << std::endl;
        //std::cout << "Hit charge: " << hit_charge[i] << std::endl;

        uint64_t tp_t_begin = timestamp + clocksPerTPCTick * ((int64_t)hit_end[i] - (int64_t)hit_tover[i]);
        uint64_t tp_t_peak  = tp_t_begin + clocksPerTPCTick * hit_peak_time[i];

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
          trigprim.time_peak = tp_t_peak;
          trigprim.time_over_threshold = uint64_t((hit_tover[i] - 1) * clocksPerTPCTick);
          trigprim.channel = chan[i];
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_peak_adc[i];
          trigprim.detid = 666;           
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
          trigprim.version = 1;

          if (save_trigprim){
            save_hit_data(trigprim, "AVX", out_suffix);
          }          

        ++total_hits;
      }
    } // loop over 16 registers   

  } // while not magic    
}


//void execute_tpg(dunedaq::fddetdataformats::WIBEthFrame* output_frame, const int ifr, const int input_ch) {
void execute_tpg(PattgenInfo& info) {

  // Parse the WIBEth frames
  uint64_t timestamp = info.output_frame->get_timestamp();      
  uint64_t expected_timestamp = first_timestamp + info.iframe*2048;
  if (info.iframe>0 && timestamp != expected_timestamp) {
    int64_t ts_diff = timestamp - expected_timestamp;
    std::cout << " |______ TIMESTAMP WILL BE OVERWRITTEN! OLD: " << timestamp << ", NEW: " << expected_timestamp <<  std::endl;
    info.output_frame->set_timestamp(expected_timestamp);
  }
  
  m_assigned_pattgen_function(*ph.m_pattgen_info);

}


// =================================================================
//                       MAIN
// =================================================================
int 
main(int argc, char** argv)
{

    CLI::App app{ "Generate TPG test patterns" };

    // Set default input frame file
    std::string frame_file_path = "./wibeth-frames.bin";
    app.add_option("-f,--frame-file-path", frame_file_path, "Path to the input frame file");

    int num_frames_to_read = -1;
    app.add_option("-n,--num-frames-to-read", num_frames_to_read, "Number of frames to read. Default: select all frames.");

    size_t input_ch = 0;
    app.add_option("-i,--input_channel", input_ch, "Input channel number for adding fake hit. Default: 0");

    int tpg_threshold = 499;
    app.add_option("-t,--tpg-threshold", tpg_threshold, "Value of the TPG threshold. Default: 499");

    app.add_flag("--save-adc-data", save_adc_data, "Save ADC data");

    app.add_flag("--save-trigprim", save_trigprim, "Save trigger primitive data");

    // additional options 
    app.add_option("-s ,--out_suffix", out_suffix, "Append string to output hit file name");

    int time_tick_offset = 1;
    app.add_option("-o,--time-tick-offset", time_tick_offset, "Time tick of pattern start. Default: 1 (max:63)");

    app.add_option("-p,--select-pattern", select_pattern, "Test Pattern (patt_golden, patt_pulse, patt_edge_square, patt_edge_left, patt_edge_right). Default: patt_golden");
 
    app.set_config("--config", "patt.cfg", "Read a configuration file. Default: patt.cfg", false);

    CLI11_PARSE(app, argc, argv);

    if (select_pattern == "patt_golden") {
      m_assigned_pattgen_function = &pattgen_function_golden;
    } else if (select_pattern == "patt_pulse") {
      m_assigned_pattgen_function = &pattgen_function_pulse;
    } else if (select_pattern == "patt_edge_square") {
      m_assigned_pattgen_function = &pattgen_function_edge_square;
    } else if (select_pattern == "patt_edge_left") {
      m_assigned_pattgen_function = &pattgen_function_edge_left;
    } else if (select_pattern == "patt_edge_right") {
      m_assigned_pattgen_function = &pattgen_function_edge_right;
    }
    std::cout << "Created an instance of the " << select_pattern << " pattern generator " << std::endl;

    
    // =================================================================
    //                       READ THE FRAMES.BIN FILE
    // =================================================================
    FrameFile input_file = FrameFile(frame_file_path.c_str()); 
    int total_num_frames = input_file.num_frames();

    std::cout << "Input file name: " << frame_file_path << std::endl;
    std::cout << "Size of the input file: " << input_file.length() << std::endl;
    std::cout << "Input channel number: " << input_ch << std::endl;
    std::cout << "Number of DUNE WIBEth frames in the input file: " << total_num_frames << std::endl;

    // Check if the selected number of frames is <= than the ones available in the input file
    if (total_num_frames < num_frames_to_read) {
      std::cout << "\n**ERROR**: Select a valid number of frames that is less or equal to the ones available in the input file." << std::endl;
      return 1;
    } else if (num_frames_to_read == -1) {
      num_frames_to_read = total_num_frames;
    }

    //std::unique_ptr<PattgenInfo<> m_pattgen_info;
    //m_pattgen_info = std::make_unique<PattgenInfo<>(500);
    ph.initialize();
    ph.m_pattgen_info->time_tick_offset = time_tick_offset;

    // =================================================================
    //                  Process the DUNE WIBEth frames
    // =================================================================
    int wibeth_frame_index = 0; 
    auto start_test = std::chrono::high_resolution_clock::now();  

    first_timestamp = input_file.frame(0)->get_timestamp();
 
    std::string out_prefix = select_pattern + "_" + std::to_string(time_tick_offset);    

    // Loop over the DUNEWIB Ethernet frames in the file
    std::fstream output_file;
    output_file.open(out_prefix+"_wibeth_output.bin", std::ios::app | std::ios::binary);
    //dunedaq::fddetdataformats::WIBEthFrame* output_frame;
    while (wibeth_frame_index < num_frames_to_read ){      

      // current WIBEth frame
      ph.m_pattgen_info->output_frame = input_file.frame(wibeth_frame_index);
      //ph.m_pattgen_info->output_frame = output_frame; 
      ph.m_pattgen_info->iframe = wibeth_frame_index;
      ph.m_pattgen_info->input_ch = input_ch;
      //execute_tpg(output_frame, wibeth_frame_index, input_ch);
      execute_tpg(*ph.m_pattgen_info);

      ++wibeth_frame_index;

      output_file.write(reinterpret_cast<char*>(ph.m_pattgen_info->output_frame), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );
    }
    output_file.close();

  // =================================================================
  //  Run TPG Hit Finder oon the generated pattern DUNE WIBEth frames
  // =================================================================

  if (save_trigprim) {
    std::cout << "Run pedestal subtraction " << std::endl;

    bool do_pedsub = true;  // pedestal subtraction is mandatory but can be switched off only for debugging
    bool do_filter = false; // no FIR filter applied 

    // pedsub as stream - and reused headers for writing binary file 
    dunedaq::fddetdataformats::WIBEthFrame* output_frame_pedsub; 
    FrameFile input_file_fake = FrameFile((out_prefix+"_wibeth_output.bin").c_str()); 
    std::cout << "Size of the input file wibeth_output.bin: " << input_file_fake.length() << std::endl;
    std::cout << "Number of frames " << num_frames_to_read << std::endl;

    std::fstream output_file_pedsub;
    if (do_filter) {
      out_suffix = "_fir";
    }
    output_file_pedsub.open(out_prefix+"_wibeth_output_pedsub"+out_suffix+".bin", std::ios::app | std::ios::binary);
    // pedestal subtraction
    int16_t median_ssr[64] = { 0 };
    int16_t accum_ssr[64] = { 0 };
    int16_t median = 0;
    int16_t accum = 0;
    // filtering
    /*
    const size_t NTAPS = 8;
    uint16_t absTimeModNTAPS;
    int16_t* prev_samp = new int16_t[64*64];
    for (size_t it = 0; it < 64 * 64; ++it) {
      for (size_t jt = 0; jt < NTAPS; ++jt) {
        prev_samp[it * NTAPS + jt] = 0;
      }
    }
    const uint8_t tap_exponent = 6;
    int16_t multiplier{1 << tap_exponent};
    int16_t adcMax{INT16_MAX / multiplier};
    std::vector<int16_t> tpg_taps = swtpg_wibeth::firwin_int(7, 0.1, multiplier);
    int16_t* taps = new int16_t[tpg_taps.size()];
    if (do_filter) {
      for (size_t it = 0; it < tpg_taps.size(); ++it) {
        std::cout << "DBG fir filter taps " << it << " : " << tpg_taps[it] << std::endl;
        taps[it] = tpg_taps[it];
      }
      std::cout << "DBG fir filter consts: multiplier, adcMax, INT16_MAX: " << multiplier << ", " << adcMax << ", " << INT16_MAX << std::endl;
    }
    */

    for (size_t i=0; i<num_frames_to_read; i++) {
      std::cout << "========== FRAME_NUM " << i <<  std::endl;
      output_frame_pedsub = input_file_fake.frame(i);
      // pedsub
      for (int ch=0; ch<64; ++ch) {
        if (i==0) {
          median_ssr[ch] = output_frame_pedsub->get_adc(ch, 0);
        }	
        median = median_ssr[ch];
        accum = accum_ssr[ch];
        for (int itime=0; itime<64; ++itime) {
	  int16_t sample = output_frame_pedsub->get_adc(ch, itime);
	  //std::cout << "DBG pedsub write to file sample 1 " << sample << ", " << median << ", " << accum << std::endl;
	  if (ch == input_ch) {
	    //std::cout << "DBG frugal before " << median << ", " << accum << ", " << sample << std::endl;
            swtpg_wibeth::frugal_accum_update(median, sample, accum, 10);
	    //std::cout << "DBG frugal after " << median << ", " << accum << ", " << sample << std::endl;
            median_ssr[ch] = median;
            accum_ssr[ch] = accum;
	    sample -= median;
	    //std::cout << "DBG frugal after sample " << sample << std::endl;
	    // filtering
	    /*if (do_filter) {
              int16_t filt = swtpg_wibeth::fir_filter(sample, adcMax, NTAPS, absTimeModNTAPS, taps, prev_samp); 
	      std::cout << "Filtering before / after: " << sample << " / " << filt << " where " << absTimeModNTAPS << std::endl;
	      sample = filt >> tap_exponent;
            }*/
	  }
	  //std::cout << "DBG pedsub write to file sample 2: " << sample << ", " << median << ", " << accum << ", " << adcMax << std::endl;
	  // NB TDAQ ERROR ADC value out of range 
	  if (sample > 0 && sample < INT16_MAX) {
	    output_frame_pedsub->set_adc(ch, itime, sample);
	  } else {
            output_frame_pedsub->set_adc(ch, itime, 0);
          }
          //std::cout << "Pedsub ADC value: " << sample << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime << " \t\tiCh: " << ch <<  std::endl;
        }
      }
      output_file_pedsub.write(reinterpret_cast<char*>(output_frame_pedsub), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );  
    }
    output_file_pedsub.close();

    // pedestal subtraction - ADCs stored in long vector - works for small number of frames -> text file 
    // todo


    std::cout << "Run hit finding " << std::endl;

    // hit finding - ADCs stored in long vector - works for small number of frames -> text file
    dunedaq::fddetdataformats::WIBEthFrame* input_frame_pedsub; 
    std::string input_file_name = out_prefix+"_wibeth_output.bin";

    std::string file_name_hits = out_prefix+"_wibeth_output_hits.txt";
    if (do_pedsub) {
      input_file_name = out_prefix+"_wibeth_output_pedsub"+out_suffix+".bin";
      file_name_hits = out_prefix+"_wibeth_output_pedsub"+out_suffix+"_hits.txt";
    }
    FrameFile input_file_pedsub = FrameFile(input_file_name.c_str()); 
    std::cout << "Size of the input file " << input_file_name << ": " << input_file_pedsub.length() << std::endl;
    std::cout << "Number of frames " << num_frames_to_read << std::endl;

    std::ofstream output_file_pedsub_hits;
    output_file_pedsub_hits.open(file_name_hits.c_str(), std::ofstream::app);

    for (int ch=0; ch<64; ++ch) {
      if (ch != input_ch) continue;
      //std::vector<int16_t> adcs;
      //std::vector<int> adcs;
      std::vector<uint16_t> adcs;
      std::vector<std::vector<int>> tmp_out;
      uint64_t timestamp;
      for (size_t i=0; i<num_frames_to_read; i++) {
        std::cout << "========== FRAME_NUM " << i << " max " << std::vector<int>().max_size() <<  std::endl;
        //std::cout << "========== FRAME_NUM " << i << " max " << std::vector<std::vector<int> >().max_size() <<  std::endl;
        input_frame_pedsub = input_file_pedsub.frame(i);
	if (i==0) {
	  timestamp = input_frame_pedsub->get_timestamp();
	}
        for (int itime=0; itime<64; ++itime) {
	  //int16_t adc = input_frame_pedsub->get_adc(ch, itime);
	  //int adc = (int)input_frame_pedsub->get_adc(ch, itime);
	  uint16_t adc = input_frame_pedsub->get_adc(ch, itime);
	  adcs.push_back(adc);
        }
      }
      hit_finder(adcs, tmp_out, ch, tpg_threshold, timestamp);
      for (auto& it : tmp_out) {
        // start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]
        uint64_t ts_tov = (it[1]-it[0])*32;
        uint64_t ts_start = it[0] * 32 + timestamp;
        uint64_t ts_peak = it[2] * 32 + timestamp;
        output_file_pedsub_hits << it[3] << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << it[4] << "," << it[5] << "\n";
	++total_hits;
      }
    }
    output_file_pedsub_hits.close();
  } // end save_trigprim


    //
    std::cout << "\n\n===============================" << std::endl;
    std::cout << "Found in total " << total_hits << " hits." << std::endl;
    
    std::cout << "\n\nFinished testing." << std::endl;

}


