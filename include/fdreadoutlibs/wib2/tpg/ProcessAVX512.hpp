/**
 * @file ProcessAVX512.hpp Process frames with AVX512 registers and instructions
 * @author Philip Rodrigues (rodriges@fnal.gov)
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIB2_TPG_PROCESSAVX512_HPP_
#define READOUT_SRC_WIB2_TPG_PROCESSAVX512_HPP_

#include "FrameExpand.hpp"
#include "ProcessingInfo.hpp"
#include "TPGConstants_wib2.hpp"

#include <immintrin.h>

namespace swtpg_wib2 {

inline void
frugal_accum_update_avx512(__m512i& __restrict__ median,
                         const __m512i s,
                         __m512i& __restrict__ accum,
                         const int16_t acclimit,
                         const __mmask32 mask) __attribute__((always_inline));

inline void
frugal_accum_update_avx512(__m512i& __restrict__ median,
                         const __m512i s,
                         __m512i& __restrict__ accum,
                         const int16_t acclimit,
                         const __mmask32 mask)
{
  // if the sample is greater than the median, add one to the accumulator
  // if the sample is less than the median, subtract one from the accumulator.

  __mmask32 is_gt = _mm512_cmpgt_epi16_mask(s, median);
  __m512i to_add = _mm512_set1_epi16(1);
  is_gt = _kand_mask32(is_gt, mask);
  accum = _mm512_mask_add_epi16(accum, is_gt, accum, to_add);

  __mmask32 is_lt = _mm512_cmplt_epi16_mask(s, median);
  __m512i to_sub = _mm512_set1_epi16(-1);
  is_lt = _kand_mask32(is_lt, mask);
  accum = _mm512_mask_add_epi16(accum, is_lt, accum, to_sub);

  // if the accumulator is >10, add one to the median and
  // set the accumulator to zero. if the accumulator is
  // <-10, subtract one from the median and set the
  // accumulator to zero
  is_gt = _mm512_cmpgt_epi16_mask(accum, _mm512_set1_epi16(acclimit));
  is_gt = _kand_mask32(is_gt, mask);
  median = _mm512_mask_adds_epi16(median, is_gt, median, to_add);
  is_lt = _mm512_cmplt_epi16_mask(accum, _mm512_set1_epi16(-1 * acclimit));
  is_lt = _kand_mask32(is_lt, mask);
  median = _mm512_mask_adds_epi16(median, is_lt, median, to_sub);


  // Reset the unmasked channels that were >10 or <-10 to zero, leaving the others unchanged
  __mmask32 need_reset = _kor_mask32(is_gt, is_lt);
  need_reset = _kand_mask32(need_reset, mask);
  accum = _mm512_mask_set1_epi16(accum, need_reset, 0);
}

template<size_t NREGISTERS>
inline void
process_window_avx512(ProcessingInfo<NREGISTERS>& info)
{
  // Start with taps as floats that add to 1. Multiply by some
  // power of two (2**N) and round to int. Before filtering, cap the
  // value of the input to INT16_MAX/(2**N)
  const size_t NTAPS = 8;
  // int16_t taps[NTAPS] = {0};
  // for (size_t i = 0; i < std::min(NTAPS, info.tapsv.size()); ++i) {
  //  taps[i] = info.tapsv[i];
  //}

  const __m512i adcMax = _mm512_set1_epi16(info.adcMax);
  // The maximum value that sigma can have before the threshold overflows a 16-bit signed integer
  const __m512i sigmaMax = _mm512_set1_epi16((1 << 15) / (info.multiplier * 5));

  __m512i tap_512[NTAPS];
  for (size_t i = 0; i < NTAPS; ++i) {
    tap_512[i] = _mm512_set1_epi16(info.taps[i]);
  }
  // Pointer to keep track of where we'll write the next output hit
  __m512i* output_loc = (__m512i*)(info.output); // NOLINT(readability/casting)

  const __m512i iota = _mm512_set_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
		  15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

  int nhits = 0;

  for (uint16_t ireg = info.first_register; ireg < info.last_register; ++ireg) { // NOLINT(build/unsigned)
    std::cout << "ireg = " << ireg << std::endl;

    uint16_t absTimeModNTAPS = info.absTimeModNTAPS; // NOLINT(build/unsigned)

    // ------------------------------------
    // Variables for pedestal subtraction

    // The current estimate of the pedestal in each channel: get
    // from the previous go-around.

    ChanState<NREGISTERS>& state = info.chanState;
    __m512i median = _mm512_loadu_epi16(reinterpret_cast<__m512i*>(state.pedestals) + ireg);      // NOLINT

    // The accumulator that we increase/decrease when the current
    // sample is greater/less than the median
    __m512i accum = _mm512_loadu_epi16(reinterpret_cast<__m512i*>(state.accum) + ireg);     // NOLINT

    // ------------------------------------
    // Variables for filtering

    // The (unfiltered) samples `n` places before the current one
    __m512i prev_samp[NTAPS];
    for (size_t j = 0; j < NTAPS; ++j) {
      prev_samp[j] = _mm512_loadu_epi16(reinterpret_cast<__m512i*>(state.prev_samp) + NTAPS * ireg + j); // NOLINT
    }

    // ------------------------------------
    // Variables for hit finding

    // Was the previous step over threshold?
    __mmask32 prev_was_over = *(reinterpret_cast<__mmask32*>(state.prev_was_over) + ireg); // NOLINT
    ;
    // The integrated charge (so far) of the current hit
    __m512i hit_charge = _mm512_loadu_epi16(reinterpret_cast<__m512i*>(state.hit_charge) + ireg); // NOLINT
    ;
    // The time-over-threshold (so far) of the current hit
    __m512i hit_tover = _mm512_loadu_epi16(reinterpret_cast<__m512i*>(state.hit_tover) + ireg); // NOLINT
    ;

    // The channel numbers in each of the slots in the register
    __m512i channel_base = _mm512_set1_epi16(ireg * SAMPLES_PER_REGISTER);
    __m512i channels = _mm512_add_epi16(channel_base, iota);

      std::cout << "median :";
      print512_as16(median);
      std::cout << "\n";
      std::cout << "accum  :";
      print512_as16(accum);
      std::cout << "\n\n";

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      const size_t msg_index = itime / 12;
      const size_t msg_time_offset = itime % 12;
      const size_t index = msg_index * NREGISTERS * FRAMES_PER_MSG + FRAMES_PER_MSG * ireg + msg_time_offset;
      // const __m256i* rawp=reinterpret_cast<const __m256i*>(info.input)+index; // NOLINT

      // The current sample
      __m512i s = info.input->ymm(index);
      std::cout << "itime = " << itime << std::endl;
      std::cout << "sample :";
      print512_as16(s);
      std::cout << "\n";


      // Update the median itself in all channels
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverflow"
      frugal_accum_update_avx512(median, s, accum, 10, _cvtu32_mask32(0xffffffff));
#pragma GCC diagnostic pop
      // Actually subtract the pedestal
      s = _mm512_sub_epi16(s, median);


      // __m512i sigma = _mm512_set1_epi16(2000); // 20 ADC
      // --------------------------------------------------------------
      // Filtering
      // --------------------------------------------------------------

      // Don't let the sample exceed adcMax, which is the value
      // at which its filtered version might overflow
      s = _mm512_min_epi16(s, adcMax);



      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      // Mask for channels that are over the threshold in this step
      // const uint16_t threshold=2000; // NOLINT(build/unsigned)
      __m512i threshold = _mm512_set1_epi16(info.threshold);
      __mmask32 is_over = _mm512_cmpgt_epi16_mask(s, threshold);
      // Mask for channels that left "over threshold" state this step
      __mmask32 left = _kand_mask32(_knot_mask32(is_over), prev_was_over);

      //-----------------------------------------
      // Update hit start times for the channels where a hit started
      const __m512i timenow = _mm512_set1_epi16(itime);
      //-----------------------------------------
      // Accumulate charge and time-over-threshold in the is_over channels

      // Really want an epi16 version of this, but the cmpgt and
      // cmplt functions set their epi16 parts to 0xff or 0x0,
      // so treating everything as epi8 works the same
      //__m256i to_add_charge = _mm256_blendv_epi8(_mm256_set1_epi16(0), filt, is_over);
      // Divide by the multiplier before adding (implemented as a shift-right)
      __m512i to_add_charge = _mm512_maskz_srai_epi16(is_over, s, info.tap_exponent);
      hit_charge = _mm512_adds_epi16(hit_charge, to_add_charge);

      // if(ireg==2){
      //     printf("itime=%ld\n", itime);
      //     printf("s:             "); print256_as16_dec(s);             printf("\n");
      //     printf("median:        "); print256_as16_dec(median);        printf("\n");
      //     printf("sigma:         "); print256_as16_dec(sigma);         printf("\n");
      //     printf("filt:          "); print256_as16_dec(filt);          printf("\n");
      //     printf("to_add_charge: "); print256_as16_dec(to_add_charge); printf("\n");
      //     printf("hit_charge:    "); print256_as16_dec(hit_charge);    printf("\n");
      //     printf("left:          "); print256_as16_dec(left);          printf("\n");
      // }

      //__m256i to_add_tover = _mm256_blendv_epi8(_mm256_set1_epi16(0), _mm256_set1_epi16(1), is_over);
      //hit_tover = _mm256_adds_epi16(hit_tover, to_add_tover);
      hit_tover = _mm512_mask_add_epi16(hit_tover, is_over, hit_tover, _mm512_set1_epi16(1));

      // Only store the values if there are >0 hits ending on
      // this sample. We have to save the entire 16-channel
      // register, which is inefficient, but whatever
      if (_cvtmask32_u32(left) != 0) {

        ++nhits;
        // We have to save the whole register, including the
        // lanes that don't have anything interesting, but
        // we'll mask them to zero so they're easy to remove
        // in a later processing step.
        //
        // (TODO: Maybe we should do that processing step in this function?)
        //#define STORE_MASK(x) _mm256_storeu_si256(output_loc++, _mm256_blendv_epi8(_mm256_set1_epi16(0), x, left));

        _mm512_storeu_si512(output_loc++, channels); // NOLINT(runtime/increment_decrement)
        // Store the end time of the hit, not the start
        // time. Since we also have the time-over-threshold,
        // we can calculate the absolute 64-bit start time in
        // the caller. This saves faffing with hits that span
        // a message boundary, hopefully

        _mm512_storeu_si512(output_loc++, timenow); // NOLINT(runtime/increment_decrement)
        // STORE_MASK(hit_charge);
        const __m512i zero = _mm512_setzero_si512();
        _mm512_storeu_si512(output_loc++, // NOLINT(runtime/increment_decrement)
                            _mm512_maskz_add_epi16(_knot_mask32(left), hit_charge, zero));
        _mm512_storeu_si512(output_loc++, hit_tover); // NOLINT(runtime/increment_decrement)

        // reset hit_start, hit_charge and hit_tover in the channels we saved
        hit_charge = _mm512_maskz_add_epi16(left, hit_charge, zero);
        hit_tover = _mm512_maskz_add_epi16(left, hit_tover, zero);
      } // end if(!no_hits_to_store)

      prev_was_over = is_over;
      std::cout << "median :";
      print512_as16(median);
      std::cout << "\n";
      std::cout << "accum  :";
      print512_as16(accum);
      std::cout << "\n\n";
    } // end loop over itime (times for this register)
      std::cout << "\n";

    // Store the state, ready for the next time round
    _mm512_storeu_si512(reinterpret_cast<__m512i*>(state.pedestals) + ireg, median);      // NOLINT

    _mm512_storeu_si512(reinterpret_cast<__m512i*>(state.accum) + ireg, accum);     // NOLINT


    // Tricky saving mask rather than __m512i
    __mmask32* msk=reinterpret_cast<__mmask32*>(reinterpret_cast<__m512i*>(state.prev_was_over) + ireg);
   *msk = prev_was_over; // NOLINT
    _mm512_storeu_si512(reinterpret_cast<__m512i*>(state.hit_charge) + ireg, hit_charge);       // NOLINT
    _mm512_storeu_si512(reinterpret_cast<__m512i*>(state.hit_tover) + ireg, hit_tover);         // NOLINT

  } // end loop over ireg (the 8 registers in this frame)
      std::cout << "\n";

  info.absTimeModNTAPS = (info.absTimeModNTAPS + info.timeWindowNumFrames) % NTAPS;
  // Store the output
  for (int i = 0; i < 4; ++i) {
    _mm512_storeu_si512(output_loc++, _mm512_set1_epi16(swtpg_wib2::MAGIC)); // NOLINT(runtime/increment_decrement)
  }

  info.nhits = nhits;

} // NOLINT(readability/fn_size)

} // namespace swtpg_wib2

#endif // READOUT_SRC_WIB2_TPG_PROCESSAVX2_HPP_
