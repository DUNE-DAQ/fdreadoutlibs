/**
 * @file WIBEthBinaryFrameReader.cxx: binary frame reader
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/WIBEthFrame.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessNaive.hpp"

// system
#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <sstream>

#include <cstring>
#include <immintrin.h>
#include <cstdio> // For printf
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>
#include <cstdlib> // abs

#include "hdf5libs/HDF5RawDataFile.hpp"
#include "hdf5libs/hdf5filelayout/Nljs.hpp"

#include "logging/Logging.hpp"

#include <iostream>
#include <sstream>
#include <string>

#include "fdreadoutlibs/wibeth/tpg/DesignFIR.hpp"
#include "TestTPGAlgorithmsWIBEth.hpp"

using namespace dunedaq::hdf5libs;
using namespace dunedaq::daqdataformats;

std::map<std::string, std::string> cfg;
int threshold = 0;

void
print_usage()
{
  TLOG() << "Usage: WIBEthBinaryFrameReader <input_file_name> <input_channel>";
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


void hit_finder(std::vector<int>& adcs, std::vector<std::vector<int>>& out, const int& channel, 
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
      std::vector<int>::iterator it_adcs = adcs.begin();
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
	std::vector<int>::iterator it_start = adcs.begin()+start[j];
	std::vector<int>::iterator it_end = adcs.begin()+end[j]+1;
        int sum = accumulate(it_start, it_end, 0);
        sums.push_back(sum);
      }

      // find peak adcs and times
      std::vector<int> peak_adcs;
      std::vector<int> peak_times;
      for(long unsigned int j = 0; j < start.size(); j++){
	std::vector<int>::iterator it_start = adcs.begin()+start[j];
	std::vector<int>::iterator it_end = adcs.begin()+end[j]+1;
	std::vector<int>::iterator max = max_element(it_start, it_end);
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



int main(int argc, char** argv)
{
 
  //if (argc != 3) {
  //  print_usage();
  //  return 1;
  //}
  //const std::string ifile_name = std::string(argv[1]);
  //const size_t input_ch = atoi(argv[2]);
  
  CLI::App app{ "Generate TPG test patterns" };

  std::string ifile_name = "./wibeth-frames.bin";
  app.add_option("-f,--frame_file_path", ifile_name, "Path to the input frame file. Default: ./wibeth-frames.bin");

  size_t input_ch = 0;
  app.add_option("-i,--input_channel", input_ch, "Input channel number for adding fake hit. Default: 0");
  
  std::string app_cfg_fn = "patt.cfg";
  app.add_option("-c", app_cfg_fn, "App config file. Default: patt.cfg");  
  

  std::ifstream file(app_cfg_fn);
  AppCfg ac = AppCfg(app_cfg_fn);
  ac.parse(cfg);
  ac.print(cfg);
  
  // Read file
  FrameFile input_file = FrameFile(ifile_name.c_str()); 

  std::cout << "Size of the input file " << input_file.length() << std::endl;
  std::cout << "Number of frames " << input_file.num_frames() << std::endl;
 
  // Number of frames to process
  int num_frames = input_file.num_frames(); 
  if (cfg.count("patt_time_offset") == 1) {
    num_frames = stoi(cfg["num_frames"]);
  }


  // Write output file 
  std::fstream output_file;
  std::string patt_name, name_suffix; 
  int patt_time_offset = 1;
  if (cfg.count("patt_time_offset") == 1) {
    patt_time_offset = stoi(cfg["patt_time_offset"]);
  }
 
  for (auto& it : cfg) {
    if (it.second == "true") {
      if (it.first.substr(0,5) == "patt_") {
        patt_name = it.first.substr(0) + "_" + std::to_string(patt_time_offset);
      }
    }
  }
  std::cout << "Pattern name from configuration: " << patt_name << std::endl;

  output_file.open(patt_name+"_wibeth_output.bin", std::ios::app | std::ios::binary);

  uint64_t timestamp_first = input_file.frame(0)->get_timestamp();

  dunedaq::fddetdataformats::WIBEthFrame* output_frame; 
  for (size_t i=0; i<num_frames; i++) {
    std::cout << "========== FRAME_NUM " << i <<  std::endl;
    output_frame = input_file.frame(i);

    uint64_t timestamp = output_frame->get_timestamp();
    uint64_t ts_expect = timestamp_first + i*2048;
    if (i>0 && timestamp != ts_expect) {
      int64_t ts_diff = timestamp - ts_expect;
      std::cout << " |______ TIMESTAMP WILL BE OVERWRITTEN! OLD: " << timestamp << ", NEW: " << ts_expect <<  std::endl;
      output_frame->set_timestamp(ts_expect);
    }

    // default pattern  
    if (cfg.count("patt_default") == 1 && cfg["patt_default"] == "true") {
      std::cout << "********** PATT_DEFAULT " << std::endl;
      for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          output_frame->set_adc(ch, itime, 0);
        }
        if (itime == 0 && i==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {	
          output_frame->set_adc(input_ch, itime, 666);      
        }
        uint16_t adc_val = output_frame->get_adc(input_ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
    }
    // golden pattern
    if (cfg.count("patt_golden") == 1 && cfg["patt_golden"] == "true") {
      std::cout << "********** PATT_GOLDEN " << std::endl;

      int patt_time[9]{0};
      int patt_adc[9]{500, 502, 504, 505, 506, 505, 504, 502, 500};
      patt_time[0] = patt_time_offset;
      int npatt = sizeof(patt_time) / sizeof(*patt_time);
      for (int ipatt=1; ipatt<npatt; ipatt++) {
	patt_time[ipatt] = patt_time_offset+ipatt < 64 ? patt_time_offset+ipatt : patt_time_offset-64+ipatt;
      } 

      for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          output_frame->set_adc(ch, itime, 0); // set pedestal value, make it configurable
        }
        //if (i==0 && itime < npatt && std::find(patt_time, patt_time + npatt, itime) == patt_time + npatt ) {
        if (i==0 && itime < patt_time[0]) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {
          for (int ipatt=0; ipatt<sizeof(patt_time)/sizeof(patt_time[0]); ipatt++) {
            if (itime == patt_time[ipatt]) output_frame->set_adc(input_ch, itime, patt_adc[ipatt]);
	  }
        }
        uint16_t adc_val = output_frame->get_adc(input_ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
    }
    // next pattern
    if (cfg.count("patt_edge_square") == 1 && cfg["patt_edge_square"] == "true") {
      std::cout << "********** PATT_EDGE_SQUARE " << std::endl;
      for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          output_frame->set_adc(ch, itime, 0);
        }
        if (itime >= 0 && itime<=62 && i==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {	
          if (itime == 0) output_frame->set_adc(input_ch, itime, 500);      
          if (itime == 63) output_frame->set_adc(input_ch, itime, 500);      
        }
        uint16_t adc_val = output_frame->get_adc(input_ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
    }
    // next
    if (cfg.count("patt_edge_left") == 1 && cfg["patt_edge_left"] == "true") {
      std::cout << "********** PATT_EDGE_LEFT " << std::endl;
      for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          output_frame->set_adc(ch, itime, 0);
        }
        if (itime >= 0 && itime<=62 && i==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {	
          if (itime == 0) output_frame->set_adc(input_ch, itime, 500);      
          if (itime == 63) output_frame->set_adc(input_ch, itime, 501);      
        }
        uint16_t adc_val = output_frame->get_adc(input_ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
    }
    // next pattern
    if (cfg.count("patt_edge_right") == 1 && cfg["patt_edge_right"] == "true") {
      std::cout << "********** PATT_EDGE_RIGHT " << std::endl;
      for (int itime=0; itime<64; ++itime) {
        for (int ch=0; ch<64; ++ch) {
          output_frame->set_adc(ch, itime, 0);
        }
        if (itime >= 0 && itime<=62 && i==0) {
          std::cout << "Nothing to do for first frame" << std::endl;
        } else {	
          if (itime == 0) output_frame->set_adc(input_ch, itime, 501);      
          if (itime == 63) output_frame->set_adc(input_ch, itime, 500);      
        }
        uint16_t adc_val = output_frame->get_adc(input_ch, itime);
        std::cout << "Output ADC value: " << adc_val << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime <<  std::endl;
      }
    }
    // next pattern

    // save pattern
    output_file.write(reinterpret_cast<char*>(output_frame), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );  
  }
  output_file.close();
 

  if (cfg.count("do_pedsub") == 1 && cfg["do_pedsub"] == "true") {
    std::cout << "Run pedestal subtraction " << std::endl;

    // pedsub as stream - and reused headers for writing binary file 
    dunedaq::fddetdataformats::WIBEthFrame* output_frame_pedsub; 
    FrameFile input_file_fake = FrameFile((patt_name+"_wibeth_output.bin").c_str()); 
    std::cout << "Size of the input file wibeth_output.bin: " << input_file_fake.length() << std::endl;
    //std::cout << "Number of frames " << input_file_fake.num_frames() << std::endl;
    std::cout << "Number of frames " << num_frames << std::endl;

    std::fstream output_file_pedsub;
    if (cfg.count("do_filter") == 1 && cfg["do_filter"] == "true") {
      name_suffix = "_fir";
    }
    output_file_pedsub.open(patt_name+"_wibeth_output_pedsub"+name_suffix+".bin", std::ios::app | std::ios::binary);
    int16_t median_ssr[64] = { 0 };
    int16_t accum_ssr[64] = { 0 };
    int16_t median = 0;
    int16_t accum = 0;
    // filtering
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
    for (size_t it = 0; it < tpg_taps.size(); ++it) {
      std::cout << "DBG fir filter taps " << it << " : " << tpg_taps[it] << std::endl;
      taps[it] = tpg_taps[it];
    }
    std::cout << "DBG fir filter consts: multiplier, adcMax, INT16_MAX: " << multiplier << ", " << adcMax << ", " << INT16_MAX << std::endl;
    
    for (size_t i=0; i<num_frames; i++) {
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
	  std::cout << "DBG pedsub write to file sample 1 " << sample << ", " << median << ", " << accum << std::endl;
	  if (ch == input_ch) {
	    //std::cout << "DBG frugal before " << median << ", " << accum << ", " << sample << std::endl;
            swtpg_wibeth::frugal_accum_update(median, sample, accum, 10);
	    //std::cout << "DBG frugal after " << median << ", " << accum << ", " << sample << std::endl;
            median_ssr[ch] = median;
            accum_ssr[ch] = accum;
	    sample -= median;
	    //std::cout << "DBG frugal after sample " << sample << std::endl;
	    // filtering
	    if (cfg.count("do_filter") == 1 && cfg["do_filter"] == "true") {
              int16_t filt = swtpg_wibeth::fir_filter(sample, adcMax, NTAPS, absTimeModNTAPS, taps, prev_samp); 
	      std::cout << "Filtering before / after: " << sample << " / " << filt << " where " << absTimeModNTAPS << std::endl;
	      sample = filt >> tap_exponent;
            }
	  }
	  std::cout << "DBG pedsub write to file sample 2: " << sample << ", " << median << ", " << accum << ", " << adcMax << std::endl;
	  // NB TDAQ ERROR ADC value out of range 
	  if (sample > 0 && sample < INT16_MAX) {
	    output_frame_pedsub->set_adc(ch, itime, sample);
	  } else {
            output_frame_pedsub->set_adc(ch, itime, 0);
          }
          std::cout << "Pedsub ADC value: " << sample << "\t\t\tFrame: " << i << " \t\tChannel: " << input_ch << " \t\tTimeSample: " << itime << " \t\tiCh: " << ch <<  std::endl;
        }
      }
      output_file_pedsub.write(reinterpret_cast<char*>(output_frame_pedsub), sizeof(dunedaq::fddetdataformats::WIBEthFrame) );  
    }
    output_file_pedsub.close();
  }

  // pedestal subtraction - ADCs stored in long vector - works for small number of frames -> text file 
  // todo


  if (cfg.count("do_hitfind") == 1 && cfg["do_hitfind"] == "true") {
    std::cout << "Run hit finding " << std::endl;

    // hit finding - ADCs stored in long vector - works for small number of frames -> text file
    dunedaq::fddetdataformats::WIBEthFrame* input_frame_pedsub; 
    std::string input_file_name = patt_name+"_wibeth_output.bin";

    std::string file_name_hits = patt_name+"_wibeth_output_hits.txt";
    if (cfg.count("do_pedsub") == 1 && cfg["do_pedsub"] == "true") {
      input_file_name = patt_name+"_wibeth_output_pedsub"+name_suffix+".bin";
      file_name_hits = patt_name+"_wibeth_output_pedsub"+name_suffix+"_hits.txt";
    }
    FrameFile input_file_pedsub = FrameFile(input_file_name.c_str()); 
    std::cout << "Size of the input file " << input_file_name << ": " << input_file_pedsub.length() << std::endl;
    //std::cout << "Number of frames " << input_file_pedsub.num_frames() << std::endl;
    std::cout << "Number of frames " << num_frames << std::endl;

    //std::fstream output_file_pedsub;
    //output_file_pedsub.open("wibeth_output_hits.txt", std::ios::app | std::ios::binary);
    std::ofstream output_file_pedsub;
    output_file_pedsub.open(file_name_hits.c_str(), std::ofstream::app);

    if (cfg.count("threshold") == 1) {	    
      threshold = stoi(cfg["threshold"]);
    }

    for (int ch=0; ch<64; ++ch) {
      if (ch != input_ch) continue;
      //std::vector<int16_t> adcs;
      std::vector<int> adcs;
      std::vector<std::vector<int>> tmp_out;
      uint64_t timestamp;
      for (size_t i=0; i<num_frames; i++) {
        std::cout << "========== FRAME_NUM " << i <<  std::endl;
        input_frame_pedsub = input_file_pedsub.frame(i);
	if (i==0) {
	  timestamp = input_frame_pedsub->get_timestamp();
	}
        for (int itime=0; itime<64; ++itime) {
	  //int16_t adc = input_frame_pedsub->get_adc(ch, itime);
	  int adc = (int)input_frame_pedsub->get_adc(ch, itime);
	  adcs.push_back(adc);
        }
      }
      hit_finder(adcs, tmp_out, ch, threshold, timestamp);
      for (auto& it : tmp_out) {
        // start[k], end[k], peak_times[k], channel, sums[k], peak_adcs[k], hitcontinue[k]
        uint64_t ts_tov = (it[1]-it[0])*32;
        uint64_t ts_start = it[0] * 32 + timestamp;
        uint64_t ts_peak = it[2] * 32 + timestamp;
        output_file_pedsub << it[3] << "," << ts_start << "," << ts_tov << "," << ts_peak << "," << it[4] << "," << it[5] << "\n";
      }
    }
    output_file_pedsub.close();
  }

  // end
}

/* patterns 
patt_default

patt_golden
-----------
500 502 504 505 506 505 504 502 500
threshold : 499
tov       : 8 * 32 = 256 (ticks, counted from 0)
peak_time : ts_0  + (time_offset + 4) * 32  e.g.  79554162068719943+(1+4)*32
peak_adc  : 506
sum_adc   : 4528

*/
