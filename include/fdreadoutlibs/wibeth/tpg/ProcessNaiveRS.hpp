/**
 * @file ProcessNaiveRS.hpp Non AVX implementation of AbsRS tpg algorithm
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIBETH_TPG_PROCESSNAIVERS_HPP_
#define READOUT_SRC_WIBETH_TPG_PROCESSNAIVERS_HPP_

#include "FrameExpand.hpp"
#include "ProcessingInfo.hpp"
#include "TPGConstants_wibeth.hpp"

#include <algorithm>
#include <inttypes.h>
#include <limits>

namespace swtpg_wibeth {


template<size_t NREGISTERS>
void
process_window_naive_RS(ProcessingInfo<NREGISTERS>& info)
{
  const float   R = 0.8; //"deweighting factor" for running sum
  //scaling factor to stop the ADCs from overflowing (may not needs this, depends on magnitude of FIR output) 
  const size_t  scale = 2; 

  uint16_t* output_loc = info.output;           // NOLINT
  const uint16_t* input16 = info.input->data(); // NOLINT
  int nhits = 0;

  for (size_t ichan = 0; ichan < NREGISTERS * SAMPLES_PER_REGISTER; ++ichan) {

    // AAA: debugging 
    //std::cout << "CHANNEL: " << ichan << std::endl;
    const size_t register_index = ichan / SAMPLES_PER_REGISTER;
    if (register_index < info.first_register || register_index >= info.last_register)
      continue;
    const size_t register_offset = ichan % SAMPLES_PER_REGISTER;

    const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * FRAMES_PER_MSG;

    // Get all the state variables by reference so they "automatically" get saved for the next go-round
    ChanState<NREGISTERS>& state = info.chanState;
    int16_t& median     = state.pedestals[ichan];
    int16_t& accum      = state.accum[ichan];
    
    int16_t& RS         = state.RS[ichan]; //value of the RS for the previous sample
    int16_t& medianRS   = state.pedestalsRS[ichan]; //median for the RS waveform needed for IQR & separate pedsub
    int16_t& accumRS    = state.accumRS[ichan];
    //IQR
    //int16_t& accum25    = state.accum25[ichan];
    //int16_t& accum75    = state.accum75[ichan];
    //int16_t& quantile25 = state.quantile25[ichan];
    //int16_t& quantile75 = state.quantile75[ichan];

    // Variables for hit finding
    uint16_t& prev_was_over = state.prev_was_over[ichan]; // was the previous sample over threshold?
    uint16_t& hit_charge = state.hit_charge[ichan];
    uint16_t& hit_tover = state.hit_tover[ichan]; // time over threshold
    uint16_t& hit_peak_adc = state.hit_peak_adc[ichan]; // time over threshold
    uint16_t& hit_peak_time = state.hit_peak_time[ichan]; // time over threshold    

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      const size_t msg_index = itime / swtpg_wibeth::FRAMES_PER_MSG;
      const size_t msg_time_offset = itime % swtpg_wibeth::FRAMES_PER_MSG;

      // The index in uint16_t of the start of the message we want // NOLINT
      const size_t msg_start_index = msg_index * swtpg_wibeth::ADCS_SIZE / sizeof(uint16_t); // NOLINT
      const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
      const size_t index = msg_start_index + offset_within_msg;

      // --------------------------------------------------------------
      // Pedestal finding/coherent noise removal
      // --------------------------------------------------------------
      
      int16_t sample = input16[index];

      std::stringstream ss;

      ss << "ADC value: " << sample;

      //frugal_accum_update(median, sample, accum, 10);
      frugal_accum_update((int16_t&)median, sample, (int16_t&)accum, 10);
      sample -= median;

      ss << "\tsample: " << sample;

      //--------------------------------------------------------------
      // Absolute Running Sum
      //--------------------------------------------------------------
      
      //RS = (R * RS) + std::abs(sample)/scale;       

      float first_part = R*RS;
      float second_part = (float)std::abs(sample)/scale; 
      //int16_t second_part = std::abs(sample); 

      // Round the RS result to the closest int16_t because
      // AVX code uses only int16_t
      RS = std::round(first_part+second_part); 

      ss << "  \tFirst part: " << first_part;
      ss << "  \tSecond part: " << second_part;
      ss << "  \tRS value: " << RS;

      // --------------------------------------------------------------
      // Second pedsub 
      // --------------------------------------------------------------      
      frugal_accum_update(medianRS, RS, accumRS, 10);
      RS -= medianRS;
      ss << "  \tMedianRS: " << medianRS ;

      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      bool is_over = RS > info.threshold;
      if (is_over) {
        // Simulate saturated add
        int32_t tmp_charge = hit_charge;
	tmp_charge += sample;
        tmp_charge = std::min(tmp_charge, (int32_t)std::numeric_limits<int16_t>::max());
        if (sample > hit_peak_adc) {
          hit_peak_adc = (uint16_t)sample;
          hit_peak_time = hit_tover;
        }
        hit_charge = (int16_t)tmp_charge;
        hit_tover++;
      }
      if (prev_was_over && !is_over) {

        ss << "\tis_over: " << is_over;

        // if(hit_tover==1){
        //     printf("% 5d % 5d % 5d % 5d\n", (uint16_t)ichan, (uint16_t)itime, hit_charge, hit_tover); // NOLINT
        // }

        // We reached the end of the hit: write it out
        (*output_loc++) = (uint16_t)ichan; // NOLINT
        (*output_loc++) = (uint16_t)itime;           // NOLINT
        (*output_loc++) = hit_charge;      // NOLINT
        (*output_loc++) = hit_tover;       // NOLINT
        (*output_loc++) = hit_peak_adc;    // NOLINT
        (*output_loc++) = hit_peak_time;   // NOLINT        

        hit_charge = 0;
        hit_tover = 0;
        hit_peak_adc = 0;
        hit_peak_time = 0;

        ++nhits;

      } // end if left hit
      prev_was_over = is_over;

      //std::cout << ss.str() << std::endl;      


    } // end loop over samples
  }   // end loop over channels

  // Write a magic "end-of-hits" value into the list of hits
  for (int i = 0; i < 6; ++i) {
    (*output_loc++) = MAGIC; // NOLINT
  }

  info.nhits += nhits;

  //if (nhits > 0) { 
  //  std::cout << "FOUND HITS: " << nhits << std::endl;
  //} 


  printf("Found %d hits\n", nhits);


}

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBETH_TPG_PROCESSNAIVERS_HPP_
