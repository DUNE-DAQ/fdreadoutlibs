/**
 * @file ProcessNaive.hpp Non AVX implementation of sw tpg
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIBEth_TPG_PROCESSNAIVE_HPP_
#define READOUT_SRC_WIBEth_TPG_PROCESSNAIVE_HPP_

#include "FrameExpand.hpp"
#include "ProcessingInfo.hpp"
#include "TPGConstants_wibeth.hpp"

#include <algorithm>
#include <inttypes.h>
#include <limits>

namespace swtpg_wibeth {

void
frugal_accum_update(int16_t& m, const int16_t s, int16_t& acc, const int16_t acclimit)
{
  if (s > m)
    ++acc;
  if (s < m)
    --acc;

  if (acc > acclimit) {
    ++m;
    acc = 0;
  }

  if (acc < -1 * acclimit) {
    --m;
    acc = 0;
  }
}

template<size_t NREGISTERS>
void
process_window_naive(ProcessingInfo<NREGISTERS>& info)
{
  // Start with taps as floats that add to 1. Multiply by some
  // power of two (2**N) and round to int. Before filtering, cap the
  // value of the input to INT16_MAX/(2**N)
  const size_t NTAPS = 8;
  const int16_t adcMax = info.adcMax;

  uint16_t* output_loc = info.output;           // NOLINT
  uint32_t* output_loc32u = info.output32u;           // NOLINT
  int16_t* output_loc16i = info.output16i;           // NOLINT
  const uint16_t* input16 = info.input->data(); // NOLINT
  int nhits = 0;

  printf("DBG Info data size % 5d % 5d % 5d \n", info.input->size(), swtpg_wibeth::NUM_REGISTERS_PER_FRAME, swtpg_wibeth::FRAMES_PER_MSG);
  printf("DBG Info NREGISTERS, SAMPLES_PER_REGISTER % 5d % 5d\n", NREGISTERS, SAMPLES_PER_REGISTER);

  for (size_t ichan = 0; ichan < NREGISTERS * SAMPLES_PER_REGISTER; ++ichan) {
    const size_t register_index = ichan / SAMPLES_PER_REGISTER;
    if (register_index < info.first_register || register_index >= info.last_register)
      continue;
    const size_t register_offset = ichan % SAMPLES_PER_REGISTER;

    const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * FRAMES_PER_MSG;

    printf("DBG DBG DBG  %5d:register_index %5d:register_offset %5d:register_t0_start %5d:ichan \n", register_index, register_offset, register_t0_start, ichan);

    printf("DBG ChanState ---------------------------------------------------------------------------------------------\n");
    // Get all the state variables by reference so they "automatically" get saved for the next go-round
    ChanState<NREGISTERS>& state = info.chanState;
    int16_t& median = state.pedestals[ichan];
    int16_t& accum = state.accum[ichan];

    // Variables for filtering
    int16_t* prev_samp = state.prev_samp + NTAPS * ichan;

    // Variables for hit finding
    uint16_t& prev_was_over = state.prev_was_over[ichan]; // was the previous sample over threshold?
    uint16_t& hit_charge = state.hit_charge[ichan];
    uint16_t& hit_tover = state.hit_tover[ichan]; // time over threshold
    uint16_t& hit_peak_adc = state.hit_peak_adc[ichan]; // time over threshold
    uint16_t& hit_peak_time = state.hit_peak_time[ichan]; // time over threshold

    uint16_t absTimeModNTAPS = info.absTimeModNTAPS; // NOLINT

    printf("DBG ChanState % 5d:prev_was_over % 5d:nframes \n", prev_was_over, info.timeWindowNumFrames);

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      const size_t msg_index = itime / info.timeWindowNumFrames;
      const size_t msg_time_offset = itime % info.timeWindowNumFrames;
      // The index in uint16_t of the start of the message we want // NOLINT
      const size_t msg_start_index = msg_index * swtpg_wibeth::ADCS_SIZE / sizeof(uint16_t); // NOLINT
      //const size_t msg_start_index = msg_index * swtpg_wibeth::ADCS_SIZE / sizeof(uint32_t); // NOLINT
      const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;
      const size_t index = msg_start_index + offset_within_msg;

      // --------------------------------------------------------------
      // Pedestal finding/coherent noise removal
      // --------------------------------------------------------------
      int16_t sample = input16[index]; 
      //uint16_t sample = input16[index]; // NO
      //uint32_t sample = input16[index]; // NO 
      printf("DBG DBG DBG  %5d:msg_index %5d:msg_time_offset %5d:msg_start_index, %5d:offset_within_msg \n", msg_index, msg_time_offset, msg_start_index, offset_within_msg);
      printf("DBG DBG DBG  %5d:ichan %5d:itime %5d:index %5d:sample \n", ichan, itime, index, sample);

      //frugal_accum_update(median, sample, accum, 10);
      printf("DBG Sample before pedsub % 5d:itime % 5d:msg % 5d:sample % 5d:ichan \n", itime, msg_index, (int16_t)sample, (uint16_t)ichan);
      frugal_accum_update((int16_t&)median, sample, (int16_t&)accum, 10);
      printf("DBG Sample before pedsub % 5d:sample % 5d:median % 5d:accum \n", (int16_t)sample, median, accum);

      sample -= median;

      printf("DBG Sample after pedsub % 5d % 5d % 5d \n", itime, msg_index ,(int16_t)sample);
      printf("DBG Sample after shift % 5d % 5d % 5d \n", itime, msg_index ,(int16_t)sample);

      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      bool is_over = sample > info.threshold;
      if (is_over) {
        // Simulate saturated add
        uint32_t tmp_charge = hit_charge;
        tmp_charge += sample;
        tmp_charge = std::min(tmp_charge, (uint32_t)std::numeric_limits<uint32_t>::max());
        if (sample > hit_peak_adc) {
          //hit_peak_adc = (uint16_t)sample;
          hit_peak_adc = (uint32_t)sample;
          hit_peak_time = hit_tover;
        }
        hit_charge = tmp_charge;
        hit_tover++;
        prev_was_over = is_over;
      }
      if (prev_was_over && !is_over) {
        // if(hit_tover==1){
        //     printf("% 5d % 5d % 5d % 5d\n", (uint16_t)ichan, (uint16_t)itime, hit_charge, hit_tover); // NOLINT
        // }

        // We reached the end of the hit: write it out
	(*output_loc++) = (uint16_t)ichan; // NOLINT
        (*output_loc16i++) = int16_t(itime)-1;         // NOLINT 
        (*output_loc++) = hit_charge;      // NOLINT
        (*output_loc++) = hit_tover-1;     // NOLINT
        (*output_loc++) = hit_peak_adc;    // NOLINT
        (*output_loc++) = hit_peak_time;   // NOLINT

        hit_charge = 0;
        hit_tover = 0;
        hit_peak_adc = 0;
        hit_peak_time = 0;

        ++nhits;
        prev_was_over = is_over;

      } // end if left hit

    } // end loop over samples
  }   // end loop over channels

  // printf("Found %d hits\n", nhits);
  info.nhits = nhits;
  info.absTimeModNTAPS = (info.absTimeModNTAPS + info.timeWindowNumFrames) % NTAPS;

  // Write a magic "end-of-hits" value into the list of hits
  for (int i = 0; i < 5; ++i) {
    (*output_loc++) = MAGIC; // NOLINT
  }
  for (int i = 0; i < 1; ++i) {
    (*output_loc16i++) = MAGIC16i; // NOLINT
  }
}

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBEth_TPG_PROCESSNAIVE_HPP_
