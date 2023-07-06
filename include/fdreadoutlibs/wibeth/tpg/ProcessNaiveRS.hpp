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
#include "TPGConstants_wib2.hpp"

#include <algorithm>
#include <inttypes.h>
#include <limits>

namespace swtpg_wib2 {


template<size_t NREGISTERS>
void
process_window_naive_RS(ProcessingInfo<NREGISTERS>& info, size_t channel_offset)
{
  // Start with taps as floats that add to 1. Multiply by some
  // power of two (2**N) and round to int. Before filtering, cap the
  // value of the input to INT16_MAX/(2**N)
  const size_t  NTAPS = 8;
  const int16_t adcMax = info.adcMax;
  const int16_t sigmaMax = (1 << 15) / (info.multiplier * 5);
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
    int16_t& accum25    = state.accum25[ichan];
    int16_t& accum75    = state.accum75[ichan];
    int16_t& quantile25 = state.quantile25[ichan];
    int16_t& quantile75 = state.quantile75[ichan];

    // Variables for filtering
    int16_t* prev_samp = state.prev_samp + NTAPS * ichan;

    // Variables for hit finding
    int16_t& prev_was_over = state.prev_was_over[ichan]; // was the previous sample over threshold?
    int16_t& hit_charge = state.hit_charge[ichan];
    int16_t& hit_tover = state.hit_tover[ichan]; // time over threshold

    uint16_t absTimeModNTAPS = info.absTimeModNTAPS; // NOLINT

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      const size_t msg_index = itime / 12;
      const size_t msg_time_offset = itime % 12;

      // The index in uint16_t of the start of the message we want // NOLINT
      const size_t msg_start_index = msg_index * swtpg_wib2::ADCS_SIZE / sizeof(uint16_t); // NOLINT
      const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
      const size_t index = msg_start_index + offset_within_msg;

      // --------------------------------------------------------------
      // Pedestal finding/coherent noise removal
      // --------------------------------------------------------------
      
      int16_t sample = input16[index];

      std::stringstream ss;

      ss << "ADC value: " << sample;

      frugal_accum_update(median, sample, accum, 10);
      sample -= median;

      ss << "\tsample: " << sample;

      // --------------------------------------------------------------
      // Filtering
      // --------------------------------------------------------------
      /*
      // Don't let the sample exceed adcMax, which is the value
      // at which its filtered version might overflow
      sample = std::min(sample, adcMax);
      int16_t filt_tmp = 0;
      for (size_t j = 0; j < NTAPS; ++j) {
        filt_tmp += info.taps[j] * prev_samp[(j + absTimeModNTAPS) % NTAPS];
      }
      prev_samp[absTimeModNTAPS % NTAPS] = sample;

      absTimeModNTAPS = (absTimeModNTAPS + 1) % NTAPS;
      int16_t filt = filt_tmp;
      
      ss << "\tFilt: " << filt;
      */
     
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
      // Pedsub & IQR
      // --------------------------------------------------------------
      ss << "  \t25: " << quantile25;
      ss << "  \t75: " << quantile75;

      if (RS < medianRS) {
        frugal_accum_update(quantile25, RS, accum25, 10);      
      }
      if (RS > medianRS) {
        frugal_accum_update(quantile75, RS, accum75, 10);
      }  
      frugal_accum_update(medianRS, RS, accumRS, 10);

      ss << "  \tMedianRS: " << medianRS ;


      int16_t sigma = quantile75 - quantile25;

     // The maximum value that sigma can have before the threshold overflows a 16-bit signed integer
      sigma = std::min(sigma, sigmaMax);


      ss << "  \tsigma: " << sigma;


      RS -= medianRS;

      ss << "  \tRS after pedsub: " << RS ;

      
      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      bool is_over = RS > 5 * sigma;
      if (is_over) {
        // Simulate saturated add
        int32_t tmp_charge = hit_charge;
        tmp_charge += sample >> info.tap_exponent;
        tmp_charge = std::min(tmp_charge, (int32_t)std::numeric_limits<int16_t>::max());
        hit_charge = (int16_t)tmp_charge;
        hit_tover++;
        prev_was_over = true;
      }
      if (prev_was_over && !is_over) {

        ss << "\tis_over: " << is_over;

        // if(hit_tover==1){
        //     printf("% 5d % 5d % 5d % 5d\n", (uint16_t)ichan, (uint16_t)itime, hit_charge, hit_tover); // NOLINT
        // }

        // We reached the end of the hit: write it out
        (*output_loc++) = (uint16_t)ichan; // NOLINT
        (*output_loc++) = itime;           // NOLINT
        (*output_loc++) = hit_charge;      // NOLINT
        (*output_loc++) = hit_tover;       // NOLINT

        hit_charge = 0;
        hit_tover = 0;

        ++nhits;
        prev_was_over = false;

      } // end if left hit

      //std::cout << ss.str() << std::endl;      


    } // end loop over samples
  }   // end loop over channels

  // Write a magic "end-of-hits" value into the list of hits
  for (int i = 0; i < 4; ++i) {
    (*output_loc++) = MAGIC; // NOLINT
  }

  info.nhits += nhits;

  //if (nhits > 0) { 
  //  std::cout << "FOUND HITS: " << nhits << std::endl;
  //} 

  info.absTimeModNTAPS = (info.absTimeModNTAPS + info.timeWindowNumFrames) % NTAPS;


  printf("Found %d hits\n", nhits);


}

} // namespace swtpg_wib2

#endif // READOUT_SRC_WIBETH_TPG_PROCESSNAIVERS_HPP_
