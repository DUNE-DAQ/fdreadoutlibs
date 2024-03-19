/**
 * @file ProcessRSAVX2.hpp Process frames with AVX2 registers and instructions
 * using the Running Sum algorithm
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIBEth_TPG_PROCESSRSAVX2_HPP_
#define READOUT_SRC_WIBEth_TPG_PROCESSRSAVX2_HPP_

#include "FrameExpand.hpp"
#include "UtilsAVX2.hpp"
#include "ProcessingInfo.hpp"
#include "TPGConstants_wibeth.hpp"

#include <immintrin.h>

namespace swtpg_wibeth {

template<size_t NREGISTERS>
inline void
process_window_rs_avx2(ProcessingInfo<NREGISTERS>& info)
{

  const __m256i overflowMax = _mm256_set1_epi16(INT16_MAX);
      

  // Running sum scaling factor
  //const __m256i R_factor = _mm256_set1_epi16(info.rs_memory_factor);

  // Scaling factor to stop the ADCs from overflowing 
  // (may not needs this, depends on magnitude of FIR output) 
  const __m256i scale_factor = _mm256_set1_epi16(info.rs_scale_factor);

  // The maximum value that sigma can have before the threshold overflows a 16-bit signed integer
  //const __m256i sigmaMax = _mm256_set1_epi16((1 << 15) / (info.multiplier * info.threshold));

  // Pointer to keep track of where we'll write the next output hit
  __m256i* output_loc = (__m256i*)(info.output); // NOLINT(readability/casting)

  const __m256i iota = _mm256_set_epi16(14, 13, 12, 11, 10, 9, 8, 15, 7, 6, 5, 4, 3, 2, 1, 0);

  int nhits = 0;

  for (uint16_t ireg = info.first_register; ireg < info.last_register; ++ireg) { // NOLINT(build/unsigned)

    //printf("ireg:          "); std::cout << (ireg) << std::endl;


    // ------------------------------------
    // Variables for pedestal subtraction

    // The current estimate of the pedestal in each channel: get
    // from the previous go-around.

    ChanState<NREGISTERS>& state = info.chanState;
    __m256i median = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.pedestals) + ireg);      // NOLINT
    //__m256i quantile25 = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.quantile25) + ireg); // NOLINT
    //__m256i quantile75 = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.quantile75) + ireg); // NOLINT

    // The accumulator that we increase/decrease when the current
    // sample is greater/less than the median
    __m256i accum = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.accum) + ireg);     // NOLINT
    //__m256i accum25 = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.accum25) + ireg); // NOLINT
    //__m256i accum75 = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.accum75) + ireg); // NOLINT

    // Runnin Sum variables

    __m256i RS = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.RS) + ireg);     // NOLINT
    __m256i medianRS = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.pedestalsRS) + ireg);     // NOLINT
    __m256i accumRS = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.accumRS) + ireg);     // NOLINT
    __m256i R_factor = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.RS_memory_factor) + ireg);     // NOLINT

    // ------------------------------------
    // Variables for hit finding

    // Was the previous step over threshold?
    __m256i prev_was_over = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.prev_was_over) + ireg); // NOLINT
    ;
    // The integrated charge (so far) of the current hit
    __m256i hit_charge = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_charge) + ireg); // NOLINT
    ;
    // The time-over-threshold (so far) of the current hit
    __m256i hit_tover = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_tover) + ireg); // NOLINT
    ;

    // The peak adc value
    __m256i hit_peak_adc = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_peak_adc) + ireg); // NOLINT

    // The time of the peak of the current hit
    __m256i hit_peak_time = _mm256_lddqu_si256(reinterpret_cast<__m256i*>(state.hit_peak_time) + ireg); // NOLINT    

    // The channel numbers in each of the slots in the register
    __m256i channel_base = _mm256_set1_epi16(ireg * SAMPLES_PER_REGISTER);
    __m256i channels = _mm256_add_epi16(channel_base, iota);

    for (size_t itime = 0; itime < info.timeWindowNumFrames; ++itime) {
      //printf("itime=%ld\n", itime);
      const size_t msg_index = itime / info.timeWindowNumFrames;
      const size_t msg_time_offset = itime % info.timeWindowNumFrames;
      const size_t index = msg_index * NREGISTERS * FRAMES_PER_MSG + FRAMES_PER_MSG * ireg + msg_time_offset;
      // const __m256i* rawp=reinterpret_cast<const __m256i*>(info.input)+index; // NOLINT


      // --------------------------------------------------------------
      // Pedestal finding/coherent noise removal
      // --------------------------------------------------------------


      // The current sample
      __m256i s = info.input->ymm(index);
      //printf("Input ADC value:\t\t\t\t"); print256_as16_dec(s);         printf("\n");

      // Update the median itself in all channels
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverflow"
      swtpg_wibeth::frugal_accum_update_avx2(median, s, accum, info.frugal_streaming_accumulator_limit, _mm256_set1_epi16(0xffff));
#pragma GCC diagnostic pop
      // Actually subtract the pedestal
      s = _mm256_sub_epi16(s, median);
      //printf("after pedestal:        "); print256_as16_dec(s);        printf("\n");


    


      //--------------------------------------------------------------
      // Absolute Running Sum
      //--------------------------------------------------------------
      
      // Naive: RS = (R_factor * RS) + std::abs(filt)/scale; 

      // Instead of using floats in the calcualation of the RS we multiply by 10 and 
      // do operations on the integers. In the end we divide by 10. 

     __m256i first_part = _mm256_mullo_epi16(RS, R_factor);
     //__m256i first_part_div = _mm256_div_epi16(RS, 10);

     //__m256i second_part = s;
     //__m256i second_part_div = _mm256_div_epi16(_mm256_abs_epi16(s), 10);

     //RS = _mm256_div_epi16(_mm256_add_epi16(first_part, second_part), 10);
     RS = swtpg_wibeth::_mm256_div_epi16(_mm256_add_epi16(first_part, s), 10);

     //printf("first_part:\t\t\t\t"); print256_as16_dec(first_part);         printf("\n"); 
     //printf("second_part:\t\t\t\t"); print256_as16_dec(second_part);         printf("\n"); 
     //printf("RS_value:\t\t\t\t"); print256_as16_dec(RS);         printf("\n"); 

      // Update the medianRS itself in all channels
      //printf("MedianRS:\t\t\t\t"); print256_as16_dec(medianRS);         printf("\n"); 

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverflow"
      swtpg_wibeth::frugal_accum_update_avx2(medianRS, RS, accumRS, info.frugal_streaming_accumulator_limit, _mm256_set1_epi16(0xffff));
#pragma GCC diagnostic pop

      // __m256i sigma = _mm256_set1_epi16(2000); // 20 ADC
      RS = _mm256_sub_epi16(RS, medianRS);

      //printf("RS_after_medianRS:\t\t\t\t"); print256_as16_dec(RS);         printf("\n"); 



      // --------------------------------------------------------------
      // Inter-quantile range
      // --------------------------------------------------------------

      // Find the interquartile range
      //__m256i sigma = _mm256_sub_epi16(quantile75, quantile25);
      
      // Clamp sigma to a range where it won't overflow when
      // multiplied by info.multiplier*5
      //sigma = _mm256_min_epi16(sigma, sigmaMax);



      // --------------------------------------------------------------
      // Hit finding
      // --------------------------------------------------------------
      // Mask for channels that are over the threshold in this step

      // IQR-based THRESHOLD
      //__m256i is_over = _mm256_cmpgt_epi16(RS, sigma * info.threshold);

      // FIXED THRESHOLD
      __m256i threshold = _mm256_set1_epi16(info.threshold);
      __m256i is_over = _mm256_cmpgt_epi16(RS, threshold);


      // Mask for channels that left "over threshold" state this step
      __m256i left = _mm256_andnot_si256(is_over, prev_was_over);

      //-----------------------------------------
      // Update hit start times for the channels where a hit started
      const __m256i timenow = _mm256_set1_epi16(itime);
      //-----------------------------------------
      // Accumulate charge and time-over-threshold in the is_over channels

      // Really want an epi16 version of this, but the cmpgt and
      // cmplt functions set their epi16 parts to 0xff or 0x0,
      // so treating everything as epi8 works the same
      __m256i to_add_charge = _mm256_blendv_epi8(_mm256_set1_epi16(0), s, is_over);
      hit_charge = _mm256_adds_epi16(hit_charge, to_add_charge);

      // Avoid overflow of the hit charge, if needed in practice 
      hit_charge = _mm256_min_epi16(hit_charge, overflowMax);


      //if(ireg==0){
      //     printf("itime=%ld\n", itime);
      //     printf("s:             "); print256_as16_dec(s);             printf("\n");
      //     printf("median:        "); print256_as16_dec(median);        printf("\n");
      //     printf("RS:         "); print256_as16_dec(RS);         printf("\n");
      //     printf("sigma:         "); print256_as16_dec(sigma);         printf("\n");
      //     printf("filt:          "); print256_as16_dec(filt);          printf("\n");
      //     printf("to_add_charge: "); print256_as16_dec(to_add_charge); printf("\n");
      //     printf("hit_charge:    "); print256_as16_dec(hit_charge);    printf("\n");
      //     printf("channels:    "); print256_as16_dec(channels);    printf("\n");     
      //     printf("is_over:          "); print256_as16_dec(is_over);          printf("\n");
      //     printf("left:          "); print256_as16_dec(left);          printf("\n");
      //}

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
      //printf("left:          "); print256_as16_dec(left);          printf("\n");


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
        //printf("channels:          "); print256_as16_dec(channels);          printf("\n");
        //printf("to_add_charge:    "); print256_as16_dec(to_add_charge);    printf("\n");
        //printf("hit_charge:    "); print256_as16_dec(hit_charge);    printf("\n");
        //printf("to_add_tover:    "); print256_as16_dec(to_add_tover);    printf("\n");      
        //printf("hit_tover:    "); print256_as16_dec(hit_tover);    printf("\n");      

        // Store the end time of the hit, not the start
        // time. Since we also have the time-over-threshold,
        // we can calculate the absolute 64-bit start time in
        // the caller. This saves faffing with hits that span
        // a message boundary, hopefully

        _mm256_storeu_si256(output_loc++, timenow); // NOLINT(runtime/increment_decrement)
        // STORE_MASK(hit_charge);
        _mm256_storeu_si256(output_loc++, // NOLINT(runtime/increment_decrement)
                            _mm256_blendv_epi8(_mm256_set1_epi16(0), hit_charge, left));
        _mm256_storeu_si256(output_loc++, hit_tover); // NOLINT(runtime/increment_decrement)

        _mm256_storeu_si256(output_loc++, hit_peak_adc); // NOLINT(runtime/increment_decrement)

	      _mm256_storeu_si256(output_loc++, hit_peak_time); // NOLINT(runtime/increment_decrement)	

        // Make sure to count only the channels that are above threshold
        // Store the flags per channel that indicated the waveform went below threshold
        // in order to recover the exact channel(s) that have complete hits. 	
	      _mm256_storeu_si256(output_loc++, left); // NOLINT(runtime/increment_decrement)

        // reset hit_start, hit_charge and hit_tover in the channels we saved
        const __m256i zero = _mm256_setzero_si256();
        hit_charge = _mm256_blendv_epi8(hit_charge, zero, left);
        hit_tover = _mm256_blendv_epi8(hit_tover, zero, left);
        hit_peak_adc = _mm256_blendv_epi8(hit_peak_adc, zero, left);
        hit_peak_time = _mm256_blendv_epi8(hit_peak_time, zero, left);	

      } // end if(!no_hits_to_store)
      //printf("nhits:          "); std::cout << (nhits) << std::endl;


      prev_was_over = is_over;

    } // end loop over itime (times for this register)


    // Store the state, ready for the next time round
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.pedestals) + ireg, median);      // NOLINT
    //_mm256_storeu_si256(reinterpret_cast<__m256i*>(state.quantile25) + ireg, quantile25); // NOLINT
    //_mm256_storeu_si256(reinterpret_cast<__m256i*>(state.quantile75) + ireg, quantile75); // NOLINT

    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.accum) + ireg, accum);     // NOLINT
    //_mm256_storeu_si256(reinterpret_cast<__m256i*>(state.accum25) + ireg, accum25); // NOLINT
    //_mm256_storeu_si256(reinterpret_cast<__m256i*>(state.accum75) + ireg, accum75); // NOLINT


    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.RS) + ireg, RS);     // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.pedestalsRS) + ireg, medianRS); // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.accumRS) + ireg, accumRS); // NOLINT
    
    // No need to update the memory factor as it is a fixed quantity
    //_mm256_storeu_si256(reinterpret_cast<__m256i*>(state.RS_memory_factor) + ireg, R_factor); // NOLINT

    

    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.prev_was_over) + ireg, prev_was_over); // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_charge) + ireg, hit_charge);       // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_tover) + ireg, hit_tover);         // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_peak_adc) + ireg, hit_peak_adc);         // NOLINT
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(state.hit_peak_time) + ireg, hit_peak_time);         // NOLINT    

  } // end loop over ireg (the 8 registers in this frame)


  // Store the output
  for (int i = 0; i < 7; ++i) {
    _mm256_storeu_si256(output_loc++, _mm256_set1_epi16(swtpg_wibeth::MAGIC)); // NOLINT(runtime/increment_decrement)
  }

  info.nhits = nhits;
  //printf("Found %d hits\n", nhits);


} // NOLINT(readability/fn_size)

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBEth_TPG_PROCESSRSAVX2_HPP_

