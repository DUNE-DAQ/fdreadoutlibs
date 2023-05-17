/**
 * @file WIB2AddFakeHits.cxx: quick and dirty solution to modify an input binary file and add 
 * fake hits at a predefined frame and channel numbers 
 * Usage: ./WIB2AddFakeHits /nfs/home/aabedabu/dunedaq-v3.1.1dunedaq-v3.1.1_pre_source_id/dev/test_files/output_0_3.out 
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


#include "hdf5libs/HDF5RawDataFile.hpp"
#include "hdf5libs/hdf5filelayout/Nljs.hpp"

#include "logging/Logging.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace dunedaq::hdf5libs;
using namespace dunedaq::daqdataformats;


void
print_usage()
{
  TLOG() << "Usage: WIB2AddFakeHits <input_file_name> <input_channel> <input_frame>";
}


class FrameFile
{
public:

    FrameFile(const char* filename)
        : m_file(filename, std::ifstream::binary),
          m_buffer(new char[sizeof(dunedaq::detdataformats::wib2::WIB2Frame)])
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
        //if(m_length%sizeof(dunedaq::detdataformats::wib2::WIB2Frame)!=0){
        //    throw std::runtime_error("File does not contain an integer number of frames");
        //}
        //m_n_frames=50;//m_length/sizeof(dunedaq::detdataformats::wib2::WIB2Frame);
        m_n_frames=m_length/sizeof(dunedaq::detdataformats::wib2::WIB2Frame);

	// Reinterpret the frame as WIB2Frame
	
	m_file.read(m_buffer, m_file.eof());
        m_wib2_frame = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(m_buffer);
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

    dunedaq::detdataformats::wib2::WIB2Frame* get_wib2_frame() const { return m_wib2_frame; }

    dunedaq::detdataformats::wib2::WIB2Frame* frame(size_t i)
    {
        if(i>=num_frames()) return nullptr;
        // Seek to the right place in the file
        m_file.seekg(i*sizeof(dunedaq::detdataformats::wib2::WIB2Frame));
        // Check we didn't go past the end
        if(m_file.bad() || m_file.eof()) return nullptr;
        // Actually read the fragment into the buffer
        m_file.read(m_buffer,sizeof(dunedaq::detdataformats::wib2::WIB2Frame));
        if(m_file.bad() || m_file.eof()) return nullptr;
        return reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(m_buffer);
    }
    

    

protected:
    std::ifstream m_file;
    char* m_buffer;
    dunedaq::detdataformats::wib2::WIB2Frame* m_wib2_frame = nullptr;
    size_t m_length;
    size_t m_n_frames;
};




int main(int argc, char** argv)
{
 

  if (argc != 4) {
    print_usage();
    return 1;
  }

  const std::string ifile_name = std::string(argv[1]);

  const size_t input_ch = atoi(argv[2]);
  const size_t input_fr = atoi(argv[3]);

  // Read file
  FrameFile input_file = FrameFile(ifile_name.c_str()); 

  std::cout << "Size of the input file " << input_file.length() << std::endl;
  std::cout << "Number of frames " << input_file.num_frames() << std::endl;
   
  size_t frame_number = input_fr; // frame chosen randomly
  size_t ch = input_ch;

  dunedaq::detdataformats::wib2::WIB2Frame* frame = input_file.frame(frame_number);
  
  auto val_to_check = frame->get_adc(ch);
  std::cout << "Output adc value: " << val_to_check << " for frame " << frame_number << " and channel number " << ch << ", and frame " << frame_number << std::endl;

 
  // Write output file 
  // Modified two channels in the input
  std::fstream output_file;
  output_file.open("modified_output.bin", std::ios::app | std::ios::binary);

  uint64_t timestamp = 0;
  bool pattern_square = false;
  bool pattern_triangle = true;

  dunedaq::detdataformats::wib2::WIB2Frame* output_frame; 
  for (size_t i=0; i<input_file.num_frames(); i++) {

     if (i%12 == 0) {
       timestamp += 12*32;
     }
     
     output_frame = input_file.frame(i);
     auto val_to_check2 = output_frame->get_adc(ch);
     std::cout << "Output adc value: " << val_to_check2 << " for frame " << i << " and channel number " << ch << std::endl;


     output_frame->set_timestamp(timestamp);

     if (pattern_square) {
     if (i==frame_number || i==frame_number+1 || i==frame_number+2 
         || i==frame_number+12 || i==frame_number+12+1 || i==frame_number+12+2) {
       auto adc_val = output_frame->get_adc(242);
       auto adc_val2 = output_frame->get_adc(12);
       output_frame->set_adc(242, adc_val+500);
       output_frame->set_adc(12, adc_val2+500);
     } 
     }

     if (pattern_triangle) {
       if (i==frame_number || i==frame_number+12) {
         auto adc_val = output_frame->get_adc(242);
         auto adc_val2 = output_frame->get_adc(12);
         output_frame->set_adc(242, adc_val+500);
         output_frame->set_adc(12, adc_val2+500);
       }
       if (i==frame_number+1 || i==frame_number+12+1) {
         auto adc_val = output_frame->get_adc(242);
         auto adc_val2 = output_frame->get_adc(12);
         output_frame->set_adc(242, adc_val+1000);
         output_frame->set_adc(12, adc_val2+1000);
       }
       if (i==frame_number+2 || i==frame_number+12+2) {
         auto adc_val = output_frame->get_adc(242);
         auto adc_val2 = output_frame->get_adc(12);
         output_frame->set_adc(242, adc_val+500);
         output_frame->set_adc(12, adc_val2+500);
       }
     }

     output_file.write(reinterpret_cast<char*>(output_frame), sizeof(dunedaq::detdataformats::wib2::WIB2Frame) );
  
  } 


}

