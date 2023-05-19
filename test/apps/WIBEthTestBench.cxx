/**
 * @file WIB2TestBench.cxx  Example of expanding a WIB2 frame using AVX2 
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/WIBEthFrame.hpp"
#include "fdreadoutlibs/wibeth/tpg/TPGConstants_wibeth.hpp"

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

// A little wrapper around an array of 256-bit registers, so that we
// can explicitly access it as an array of 256-bit registers or as an
// array of uint16_t
template<size_t N>
class RegisterArray
{
public:
    // Get the value at the ith position as a 256-bit register
    inline __m256i ymm(size_t i) const { return _mm256_lddqu_si256(reinterpret_cast<const __m256i*>(m_array)+i); }
    inline void set_ymm(size_t i, __m256i val) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_array)+i, val); }

    inline uint16_t uint16(size_t i) const { return m_array[i]; }
    inline void set_uint16(size_t i, uint16_t val) { m_array[i]=val; }

    // Access the jth entry in the ith register
    inline uint16_t uint16(size_t i, size_t j) const { return m_array[16*i+j]; }
    inline void set_uint16(size_t i, size_t j, uint16_t val) { m_array[16*i+j]=val; }

    inline uint16_t* data() { return m_array; }
    inline const uint16_t* data() const { return m_array; }

    inline size_t size() const { return N; }
private:
    alignas(32) uint16_t __restrict__ m_array[N*16];
};

//==============================================================================
// Print a 256-bit register interpreting it as packed 8-bit values
void print256(__m256i var)
{
    uint8_t *val = (uint8_t*) &var;
    for(int i=0; i<32; ++i){
        printf("%02x ", val[i]);
        if(i==15) printf("|| ");
        else{
            if(i%2) printf("| ");
            else    printf(" ");
        }
    }
}

//==============================================================================
// Print a 256-bit register interpreting it as packed 16-bit values
void print256_as16(__m256i var)
{
    uint16_t *val = (uint16_t*) &var;
    for(int i=0; i<32; ++i){
        printf("%02x ", val[i]);
        if(i==15) printf("|| ");
        else{
            if(i%2) printf("| ");
            else    printf(" ");
        }
    }
}

//==============================================================================
__m256i unpack_one_register( dunedaq::fddetdataformats::WIBEthFrame::word_t* ptr_to_first_word)
{
    __m256i reg=_mm256_lddqu_si256((__m256i*)ptr_to_first_word);
    //printf("Input:      ");
    //print256(reg);
    //printf("\n");

       // The register initially contains 18-and-a-bit 14-bit ADCs, but
    // we only have space for 16 after expansion, so the last 32-bit
    // word is unused. Copy word 3 so it appears twice, and move the
    // later words down one
    __m256i idx=_mm256_set_epi32(6, 5, 4, 3, 3, 2, 1, 0);
    __m256i shuf1=_mm256_permutevar8x32_epi32(reg, idx);
    //printf("shuf1:      ");
    //print256(shuf1);
    //printf("\n");

    // Each 32-bit word contains at least one full 14-bit ADC. Shift
    // the words by variable amounts so that the high 16 bits of each
    // word contains a 14-bit ADC at the right place (with the two
    // high bits still needing to be masked to zero). That result is
    // in `high_half`
    // __mmask8 mask=0xffu;
    // __m256i src=_mm256_set1_epi32(0);
    // The amounts by which we shift each 32-bit word
    __m256i count1=_mm256_set_epi32(12, 8, 4, 0, 14, 10, 6, 2);
    // __m256i high_half=_mm256_mask_sllv_epi32(src, mask, shuf1, count1);
    __m256i high_half=_mm256_sllv_epi32(shuf1, count1);
    // Mask out the low 16 bits, and the high two bits in the high half
    __m256i high_half_mask=_mm256_set1_epi32(0x3fff0000u);

    high_half=_mm256_and_si256(high_half, high_half_mask);
    // high_half2=_mm256_and_si256(high_half2, high_half_mask);

    // printf("high_half:  ");
    // print256(high_half);
    // printf("\n");

    //------------------------------------------------------------------
    // Now we start the process of setting the low 16 bits of each
    // word to the right value. This is trickier because now the bits
    // are spread across two words. First, left-shift each word so
    // that the higher bits of the ADC are in the right place
    __m256i count2=_mm256_set_epi32(10, 6, 2, 0, 12, 8, 4, 0);
    __m256i shift2=_mm256_sllv_epi32(shuf1, count2);
    // printf("shift2:     ");
    // print256(shift2);
    // printf("\n");

    // Next, permute the register so that the words containing the low
    // bits of the ADCs we want are in the same positions as the words
    // containing the corresponding high bits. This just amounts to
    // moving the words down by one
    __m256i idx2=_mm256_set_epi32(5, 4, 3, 2, 2, 1, 0, 0);
    __m256i shuf2=_mm256_permutevar8x32_epi32(reg, idx2);
    // printf("shuf2:      ");
    // print256(shuf2);
    // printf("\n");

    // Shift each word right by the amount that brings those low bits
    // into the right place, putting the result in `shift3`
    __m256i count3=_mm256_set_epi32(22, 26, 30, 0, 20, 24, 28, 0);
    __m256i shift3=_mm256_srlv_epi32(shuf2, count3);
    // printf("shift3:     ");
    // print256(shift3);
    // printf("\n");

    // OR together the registers containing the high and low bits of
    // the ADCs. At this point, the low 16 bits of each word should
    // contain the 14 bits of the ADCs in the right place (with the
    // two high bits still needing to be masked out)
    __m256i low_half=_mm256_or_si256(shift2, shift3);
    // Mask out the high 16 bits, and the high two bits in the high half
    __m256i low_half_mask=_mm256_set1_epi32(0x3fffu);
    low_half=_mm256_and_si256(low_half, low_half_mask);
    // printf("low_half:   ");
    // print256(low_half);
    // printf("\n");

    // Nearly there... Now we OR together the low and high halves
    __m256i both=_mm256_or_si256(low_half, high_half);
    // zero out the slot where we want to put the 16th value
    both=_mm256_andnot_si256(_mm256_set_epi32(0, 0, 0, 0xffffu, 0, 0, 0, 0), both);
    // printf("both:       ");
    // print256(both);
    // printf("\n");

    // We just missed the 16th value, and the lw 16 bits of the 8th
    // word are available, so shuffle it around to put it there
    __m256i shift4=_mm256_srli_epi32(reg, 18);
    // Mask so that's the only nonzero thing
    shift4=_mm256_and_si256(_mm256_set_epi32(0, 0x3fffu, 0, 0, 0, 0, 0, 0), shift4);
    // Move the word containing the value we want into the position we want
    __m256i idx3=_mm256_set_epi32(0, 0, 0, 6, 0, 0, 0, 0);
    __m256i shuf3=_mm256_permutevar8x32_epi32(shift4, idx3);

    both=_mm256_or_si256(both, shuf3);
    //printf("both':      ");
    //print256_as16(both);
    //printf("\n");
    //printf("both':      ");
    //print256(both);
    //printf("\n");

    return both;
}





RegisterArray<4*64> unpack_wibeth( dunedaq::fddetdataformats::WIBEthFrame& frame)
{

    // Number of time samples (TS) per frame
    int time_samples_per_frame = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;
   
    // Number of ADC words per TS
    int num_adc_words_per_ts = dunedaq::fddetdataformats::WIBEthFrame::s_num_adc_words_per_ts;

    
    RegisterArray<4*64> ret;

    // Define a pointer to walk the rows of the WIBEth frame which is 2D array.
    //dunedaq::fddetdataformats::WIBEthFrame::word_t (*ptr)[14] = frame.adc_words;
    auto frame_words_ptr = frame.adc_words;

    for (int i = 0; i < time_samples_per_frame; i++) {

      // The register index is used to decide on which of the 
      // 4 registers we want to unpack the ADC messages
      int reg_index = 0;

      for (int j = 0; j < num_adc_words_per_ts; j++) {

          // The words repeat every 7 iterations. 
          // In this way we can use the same unpacking 
          // function (unpack_one_register) as for the DUNE WIBs
          if (j%7 == 0 ) {
            dunedaq::fddetdataformats::WIBEthFrame::word_t * first_half = (*(frame_words_ptr + i) + j);
                               

            // Unpack one register and add it to the register array
            ret.set_ymm(i+reg_index*time_samples_per_frame, unpack_one_register(first_half));

            reg_index += 1;

            // Increment the cursor by 224 bits to get the second part of the first time sample
            // 224 corresponds to 16 (U blocks or ADCs) times 14 which are the bits per ADC in the frame. 
            // Check the WIBEth spreadsheet for further details
            char* cursor = (char*) first_half;
            cursor += 224 / 8; // divide by 8 to get the results in bytes
            dunedaq::fddetdataformats::WIBEthFrame::word_t * second_half = (dunedaq::fddetdataformats::WIBEthFrame::word_t*) cursor;
            // Unpack another register and add it to the register array
            ret.set_ymm(i+reg_index*time_samples_per_frame, unpack_one_register(second_half));

            reg_index += 1;

    
          }       
          
      } // loop over number of words 
    } // loop over time frames
          
        

    return ret;
}




int main()
{
    
    // Test WIBEth unpacking

    dunedaq::fddetdataformats::WIBEthFrame frame;

    // Zero it out first
    std::memset(&frame, 0, sizeof(dunedaq::fddetdataformats::WIBEthFrame));

    // Number of channels per frame
    int num_channels = dunedaq::fddetdataformats::WIBEthFrame::s_channels_per_half_femb;

    // Number of time samples (TS) per frame
    int time_samples_per_frame = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;
   


    // First, loop over time frames
    for (int itime=0; itime<time_samples_per_frame; ++itime) {
      // Loop over all the channels  
      for(int i=0; i<num_channels; ++i){

        frame.set_adc(i, itime, i);           
      }
    }  

    // Check that the ADC values have been properly set
    for (int itime=0; itime<time_samples_per_frame; ++itime) {
      for(int j=0; j<num_channels; ++j){
          uint16_t adc_val = frame.get_adc(j, itime);
          //std::cout << "Index " << j << " time sample " << itime << " : " << adc_val << std::endl;
      }
    }


    
    std::array<int, 16> indices{0, 1, 2, 3, 4, 5, 6, 7, 15, 8, 9, 10, 11, 12, 13, 14};

    // Unpack the WIB Eth frame
    RegisterArray<4*64> unpacked = unpack_wibeth(frame);

    // AAA: TODO: the following quantities are hardcode for the moment. 
    int NREGISTERS = swtpg_wibeth::NUM_REGISTERS_PER_FRAME;
    int SAMPLES_PER_REGISTER = swtpg_wibeth::SAMPLES_PER_REGISTER;
    int timeWindowNumFrames = time_samples_per_frame;
    
    bool success=true;

	    
    for (size_t j = 0; j < NREGISTERS * SAMPLES_PER_REGISTER; ++j) {
      // Index of the ADC frame 
      int in_index=16*(j/16)+indices[j%16];
      
      const size_t register_offset = j % SAMPLES_PER_REGISTER; 
      const size_t register_index = j / SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * timeWindowNumFrames;

      int16_t out_val;
      for (size_t itime = 0; itime < timeWindowNumFrames; ++itime) {
        const size_t msg_index = itime / timeWindowNumFrames;
        const size_t msg_time_offset = itime % timeWindowNumFrames;
        
        const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        // The index in uint16_t of the start of the message we want. 
        const size_t index = msg_start_index + offset_within_msg;

        //const uint16_t* input16 = unpacked.data(); // NOLINT
        //out_val = input16[index]
        out_val = unpacked.uint16(index);
        // Input value into the frame
        uint16_t in_val=frame.get_adc(in_index, itime);

        //if (itime  == 0) {
          std::cout << "Time sample: " << itime << " index OUTPUT: " << index << "    value:   " << out_val << std::endl;
          std::cout << "Time sample: " << itime << " index  INPUT: " << in_index << "    value:   " << in_val << std::endl;
          std::cout << "============" << std::endl;
        //}
        if(in_val!=out_val){
            success=false;  
        }
      }
    }

    if (success) {
        std::cout << "\n\nALL THE ADC VALUES ARE MATCHING!\n\n" << std::endl;
    }
    

    std::cout << "Finished" << std::endl;




}

