/**                                                                                                                                   
 * @file ProcessingInfoRS.hpp ProcessInfoRS struct                                                                                  
 *                                                                                                                                    
 * This is part of the DUNE DAQ , copyright 2020.                                                                                     
 * Licensing/copyright details are in the COPYING file that you should have                                                           
 * received with this code.                                                                                                           
 */

#ifndef READOUT_SRC_WIBEth_TPG_PROCESSINGINFO_HPP_
#define READOUT_SRC_WIBEth_TPG_PROCESSINGINFO_HPP_

#include "FrameExpand.hpp"

#define BOOST_THREAD_PROVIDES_FUTURE_CONTINUATION
#include <boost/thread/future.hpp>

namespace swtpg_wibeth {

// The state variables for each channel in the stream, saved from the last time
template<size_t NREGISTERS>
struct ChanState
{
  ChanState()
  {
    for (size_t i = 0; i < NREGISTERS * SAMPLES_PER_REGISTER; ++i) {
      pedestals[i] = 0;
      accum[i] = 0;
      RS[i] = 0; 
      pedestalsRS[i] = 0;      
      accumRS[i] = 0;
      RS_memory_factor[i] = 0;
      accum25[i] = 0;
      accum75[i] = 0;
      prev_was_over[i] = 0;
      hit_charge[i] = 0;
      hit_tover[i] = 0;
      hit_peak_time[i] = 0;
      hit_peak_adc[i] = 0;
    }
  }

  alignas(32) int16_t __restrict__ pedestals[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) int16_t __restrict__ accum[NREGISTERS * SAMPLES_PER_REGISTER];
  
  //Variables for running sum
  alignas(32) int16_t __restrict__ RS[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) int16_t __restrict__ pedestalsRS[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) int16_t __restrict__ accumRS[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) uint16_t __restrict__ RS_memory_factor[NREGISTERS * SAMPLES_PER_REGISTER];


  //Variables for IQR
  alignas(32) int16_t __restrict__ accum25[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) int16_t __restrict__ accum75[NREGISTERS * SAMPLES_PER_REGISTER];
  
  alignas(32) int16_t __restrict__ quantile25[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) int16_t __restrict__ quantile75[NREGISTERS * SAMPLES_PER_REGISTER];

  // Variables for hit finding
  alignas(32) uint16_t __restrict__ prev_was_over[NREGISTERS * SAMPLES_PER_REGISTER]; // was the previous sample over threshold?
  alignas(32) uint16_t __restrict__ hit_charge[NREGISTERS * SAMPLES_PER_REGISTER];
  alignas(32) uint16_t __restrict__ hit_tover[NREGISTERS * SAMPLES_PER_REGISTER]; // time over threshold

  alignas(32) uint16_t __restrict__ hit_peak_time[NREGISTERS * SAMPLES_PER_REGISTER]; // time peak time
  alignas(32) uint16_t __restrict__ hit_peak_adc[NREGISTERS * SAMPLES_PER_REGISTER]; // time peak adc
};

template<size_t NREGISTERS>
struct ProcessingInfo
{
  ProcessingInfo(const RegisterArray<NREGISTERS * FRAMES_PER_MSG>* __restrict__ input_,
                 size_t timeWindowNumFrames_,
                 uint8_t first_register_,           // NOLINT
                 uint8_t last_register_,            // NOLINT
                 uint16_t* __restrict__ output_,    // NOLINT
                 const uint8_t exponent_, // NOLINT
                 uint16_t threshold_,         // NOLINT
                 uint16_t rs_memory_factor_, // NOLINT
                 uint16_t rs_scale_factor_, // NOLINT
                 int16_t frugal_streaming_accumulator_limit_, // NOLINT 
                 size_t nhits_
                ) // NOLINT
    : input(input_)
    , timeWindowNumFrames(timeWindowNumFrames_)
    , first_register(first_register_)
    , last_register(last_register_)
    , output(output_)
    , exponent(exponent_)
    , threshold(threshold_)
    , rs_memory_factor(rs_memory_factor_)
    , rs_scale_factor(rs_scale_factor_)
    , frugal_streaming_accumulator_limit(frugal_streaming_accumulator_limit_)
    , multiplier(1 << exponent)
    , adcMax(INT16_MAX / multiplier)
    , nhits(nhits_)    
  {}


  // Set the initial state from the window starting at first_msg_p
  template<size_t N>
  void setState(const RegisterArray<N>& first_tick_registers, 
                std::array<uint16_t, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER>& register_memory_factor
               )
  {
    static_assert(N >= NREGISTERS, "Wrong array size");

    // AAA: Loop through all the registers, loop through all the channels, look at the 
    // first message of the superchunk and read the ADC value. This will be used as the 
    // pedestal for the channel state
    std::cout << "Printing values of the memory factor: ";
    for (size_t j = 0; j < NREGISTERS * SAMPLES_PER_REGISTER; ++j) {
      const size_t register_offset = j % SAMPLES_PER_REGISTER; 
      const size_t register_index = j / SAMPLES_PER_REGISTER;
      const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * FRAMES_PER_MSG;

      int16_t ped=0;
      for (size_t itime = 0; itime < timeWindowNumFrames; ++itime) {
        const size_t msg_index = itime / timeWindowNumFrames;
        const size_t msg_time_offset = itime % timeWindowNumFrames;
        // The index in uint16_t of the start of the message we want.         
        const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
        const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
        const size_t index = msg_start_index + offset_within_msg;

        const uint16_t* input16 = first_tick_registers.data(); // NOLINT

        ped = input16[index];

        break; // breaking in order to select only the first entry
      }

      // Set up the channel state for the memory factor
      chanState.RS_memory_factor[j] = register_memory_factor[j];
    
      // Set the pedestals and the 25/75-percentiles
      chanState.pedestals[j] = ped;
      chanState.pedestalsRS[j] = 0;
      chanState.RS[j] = 0;
      // AAA: Quantiles are set to the pedestal value +/- 20 so that the IQR 
      // becomes above the RMS value of the input ADCs. We use the frugal 
      // streaming on the 25th/75th quantiles so that the IQR becomes
      // a good estimate of the RMS of the input. 
      chanState.quantile25[j] = ped-20;
      chanState.quantile75[j] = ped+20;
    }
    std::cout << '\n' ;
    
  }  

  const RegisterArray<NREGISTERS * FRAMES_PER_MSG>* __restrict__ input;
  size_t timeWindowNumFrames;
  uint8_t first_register;        // NOLINT
  uint8_t last_register;         // NOLINT
  uint16_t* __restrict__ output; // NOLINT
  uint8_t exponent; // NOLINT
  uint16_t threshold;   // NOLINT
  uint16_t rs_memory_factor;   // NOLINT
  uint16_t rs_scale_factor;   // NOLINT
  int16_t frugal_streaming_accumulator_limit;   // NOLINT


  int16_t multiplier;
  int16_t adcMax;
  size_t nhits;
  ChanState<NREGISTERS> chanState;
};

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBEth_TPG_PROCESSINGINFO_HPP_
