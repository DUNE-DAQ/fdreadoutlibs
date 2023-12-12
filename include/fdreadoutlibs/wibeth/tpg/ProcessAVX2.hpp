/**
 * @file ProcessAVX2.hpp 
 * Simplified hit finding algorithm
 * Process frames with AVX2 registers and instructions
 * Does not compute the FIR, uses a configurable fixed threshold
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIBEth_TPG_PROCESSAVX2_HPP_
#define READOUT_SRC_WIBEth_TPG_PROCESSAVX2_HPP_

#include "FrameExpand.hpp"
#include "UtilsAVX2.hpp"
#include "ProcessingInfo.hpp"
#include "TPGConstants_wibeth.hpp"

#include <immintrin.h>

namespace swtpg_wibeth {

template<size_t NREGISTERS>
inline void
process_window_avx2(ProcessingInfo<NREGISTERS>& info)
{
  const __m256i overflowMax = _mm256_set1_epi16(INT16_MAX);

  // Pointer to keep track of where we'll write the next output hit
  __m256i* output_loc = (__m256i*)(info.output); // NOLINT(readability/casting)

  const __m256i iota = _mm256_set_epi16(14, 13, 12, 11, 10, 9, 8, 15, 7, 6, 5, 4, 3, 2, 1, 0);

  int nhits = 0;
 
  for (uint16_t ireg = info.first_register; ireg < info.last_register; ++ireg) { // NOLINT(build/unsigned)

    // ------------------------------------
    // Variables for pedestal subtraction

    // The current estimate of the pedestal in each channel: get
    // from the previous go-around.

    ChanState<NREGISTERS>& state = info.chanState;
    __m256i median = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.pedestals) + ireg);      // NOLINT
    // The accumulator that we increase/decrease when the current
    // sample is greater/less than the median
    __m256i accum = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.accum) + ireg);     // NOLINT

    // ------------------------------------
    // Variables for hit finding

    // Was the previous step over threshold?
    __m256i prev_was_over = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.prev_was_over) + ireg); // NOLINT
    // The integrated charge (so far) of the current hit
    __m256i hit_charge = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_charge) + ireg); // NOLINT
    // The time-over-threshold (so far) of the current hit
    __m256i hit_tover = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_tover) + ireg); // NOLINT

    // The peak adc value
    __m256i hit_peak_adc = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_peak_adc) + ireg); // NOLINT

    // The time of the peak of the current hit
    __m256i hit_peak_time = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_peak_time) + ireg); // NOLINT

    // The channel numbers in each of the slots in the register
    __m256i channel_base = _mm256_set1_epi16(ireg * SAMPLES_PER_REGISTER);
    __m256i channels = _mm256_add_epi16(channel_base, iota);

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      const size_t msg_index = itime / info.timeWindowNumFrames;
      const size_t msg_time_offset = itime % info.timeWindowNumFrames;
      const size_t index = msg_index * NREGISTERS * FRAMES_PER_MSG + FRAMES_PER_MSG * ireg + msg_time_offset;
      // const __m256i* rawp=reinterpret_cast<const __m256i*>(info.input)+index; // NOLINT

      // The current sample
      __m256i s = info.input->ymm(index);
      //printf("raw ADC:          "); print256_as16_dec(s);          printf("\n");

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverflow"
      swtpg_wibeth::frugal_accum_update_avx2(median, s, accum, 10, _mm256_set1_epi16(0xffff));
#pragma GCC diagnostic pop
      // Actually subtract the pedestal
      s = _mm256_sub_epi16(s, median);

      // Don't let the sample exceed adcMax, which is the value
      // at which its filtered version might overflow
      //s = _mm256_min_epi16(s, overflowMax);

      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      // Mask for channels that are over the threshold in this step
      
      // FIXED THRESHOLD
      __m256i threshold = _mm256_set1_epi16(info.threshold);
      __m256i is_over = _mm256_cmpgt_epi16(s, threshold);

      // Mask for channels that left "over threshold" state this step
      // Comparison with previous TP
      __m256i left = _mm256_andnot_si256(is_over, prev_was_over);

      //-----------------------------------------
      // Update hit start times for the channels where a hit started
      //const __m256i timenow = _mm256_set1_epi16(itime);
      const __m256i timenow = _mm256_set1_epi16(itime);
      //-----------------------------------------
      // Accumulate charge and time-over-threshold in the is_over channels

      // Really want an epi16 version of this, but the cmpgt and
      // cmplt functions set their epi16 parts to 0xff or 0x0,
      // so treating everything as epi8 works the same
      __m256i to_add_charge = _mm256_blendv_epi8(_mm256_set1_epi16(0), s, is_over);
      // Divide by the multiplier before adding (implemented as a shift-right)
      //The saturated value of 0x7FFF(max) or 0x8000(min), respectively, is written to the destination vector.
      //hit_charge = _mm256_adds_epi16(hit_charge, to_add_charge);
      //You must control the range of values operated upon to prevent undetected overflow conditions.
      hit_charge = _mm256_add_epi16(hit_charge, to_add_charge);

      // Avoid overflow of the hit charge, if needed in practice 
      hit_charge = _mm256_min_epi16(hit_charge, overflowMax);

      if(ireg==0){
      //     printf("itime=%ld\n", itime);
      //     printf("ADC value:             "); print256_as16_dec(s);             printf("\n\n");
      //     printf("median:        "); print256_as16_dec(median);        printf("\n");
      //     printf("sigma:         "); print256_as16_dec(sigma);         printf("\n");
      //     printf("threshold:         "); print256_as16_dec(threshold);         printf("\n");
      //     printf("to_add_charge: "); print256_as16_dec(to_add_charge); printf("\n");
      //     printf("hit_charge:    "); print256_as16_dec(hit_charge);    printf("\n");
      //     printf("channels:    "); print256_as16_dec(channels);    printf("\n");
      //     printf("is_over:          "); print256_as16_dec(is_over);          printf("\n");
      //     printf("left:          "); print256_as16_dec(left);          printf("\n");
      }

      // 1. Calculation of the hit peak time and ADC
      __m256i is_sample_over_adc_peak = _mm256_cmpgt_epi16(s, hit_peak_adc);
      hit_peak_adc = _mm256_blendv_epi8(hit_peak_adc, s, is_sample_over_adc_peak); 
      hit_peak_time = _mm256_blendv_epi8(hit_peak_time, hit_tover, is_sample_over_adc_peak);

      // 2. Update of the hit time over threshold  
      __m256i to_add_tover = _mm256_blendv_epi8(_mm256_set1_epi16(0), _mm256_set1_epi16(1), is_over);
      hit_tover = _mm256_adds_epi16(hit_tover, to_add_tover);

      // 3. Stop tover from overflowing if needed in practice 
      //hit_tover = _mm256_min_epi16(hit_tover, overflowMax);
    
      // Only store the values if there are >0 hits ending on
      // this sample. We have to save the entire 16-channel
      // register, which is inefficient, but whatever

      // Testing whether a whole register is zeroes turns out to be tricky. Here's a way:
      //
      // https://stackoverflow.com/questions/22674205/is-there-an-or-equivalent-to-ptest-in-x64-assembly
      //
      // In x64 assembly, PTEST %XMM0 -> %XMM1 sets the
      // zero-flag if none of the same bits are set in %XMM0 and
      // %XMM1, and sets the carry-flag if everything that is
      // set in %XMM0 is also set in %XMM1:
      const int no_hits_to_store = _mm256_testc_si256(_mm256_setzero_si256(), left);

      if (!no_hits_to_store) {

        ++nhits;
        // We have to save the whole register, including the
        // lanes that don't have anything interesting, but
        // we'll mask them to zero so they're easy to remove
        // in a later processing step.
        //
        // (TODO: Maybe we should do that processing step in this function?)
        //#define STORE_MASK(x) _mm256_storeu_si256(output_loc++, _mm256_blendv_epi8(_mm256_set1_epi16(0), x, left));

        _mm256_storeu_si256(output_loc++, channels); // NOLINT(runtime/increment_decrement)
        // Store the end time of the hit, not the start
        // time. Since we also have the time-over-threshold,
        // we can calculate the absolute 64-bit start time in
        // the caller. This saves faffing with hits that span
        // a message boundary, hopefully

        _mm256_storeu_si256(output_loc++, timenow); // NOLINT(runtime/increment_decrement)
        // STORE_MASK(hit_charge);
        //_mm256_storeu_si256(output_loc++, // NOLINT(runtime/increment_decrement)
        //                    _mm256_blendv_epi8(_mm256_set1_epi16(0), hit_charge, left));
        _mm256_storeu_si256(output_loc++, hit_charge);

        _mm256_storeu_si256(output_loc++, hit_tover); // NOLINT(runtime/increment_decrement)

        _mm256_storeu_si256(output_loc++, hit_peak_adc); // NOLINT(runtime/increment_decrement)

	_mm256_storeu_si256(output_loc++, hit_peak_time); // NOLINT(runtime/increment_decrement)

        // reset hit_start, hit_charge and hit_tover in the channels we saved
        const __m256i zero = _mm256_setzero_si256();
        hit_charge = _mm256_blendv_epi8(hit_charge, zero, left);
        hit_tover = _mm256_blendv_epi8(hit_tover, zero, left);
        hit_peak_adc = _mm256_blendv_epi8(hit_peak_adc, zero, left);
        hit_peak_time = _mm256_blendv_epi8(hit_peak_time, zero, left);

      } // end if(!no_hits_to_store)
      prev_was_over = is_over;

    } // end loop over itime (times for this register)

    // Store the state, ready for the next time round
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.pedestals) + ireg, median);      // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.accum) + ireg, accum);     // NOLINT

    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.prev_was_over) + ireg, prev_was_over); // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_charge) + ireg, hit_charge);       // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_tover) + ireg, hit_tover);         // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_peak_adc) + ireg, hit_peak_adc);         // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_peak_time) + ireg, hit_peak_time);         // NOLINT

  } // end loop over ireg (the swtpg_wibeth::SAMPLES_PER_REGISTER, e.g. 8, registers in this frame)
  // Arguably not needed, we can avoid using MAGIC
  for (int i = 0; i < 6; ++i) {
    _mm256_storeu_si256(output_loc++, _mm256_set1_epi16(swtpg_wibeth::MAGIC)); // NOLINT(runtime/increment_decrement)
  }

  info.nhits = nhits;

} // NOLINT(readability/fn_size)

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBEth_TPG_PROCESSAVX2_HPP_
