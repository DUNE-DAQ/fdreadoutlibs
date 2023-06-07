/**
 * @file UtilsAVX2.hpp 
 * Utility methods based on AVX2
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIB2_TPG_UTILSAVX2_HPP_
#define READOUT_SRC_WIB2_TPG_UTILSAVX2_HPP_


#include <immintrin.h>

namespace swtpg_wib2 {

inline void
frugal_accum_update_avx2(__m256i& __restrict__ median,
                         const __m256i s,
                         __m256i& __restrict__ accum,
                         const int16_t acclimit,
                         const __m256i mask) __attribute__((always_inline));

inline void
frugal_accum_update_avx2(__m256i& __restrict__ median,
                         const __m256i s,
                         __m256i& __restrict__ accum,
                         const int16_t acclimit,
                         const __m256i mask)
{
  // if the sample is greater than the median, add one to the accumulator
  // if the sample is less than the median, subtract one from the accumulator.

  // For reasons that I don't understand, there's no cmplt
  // for "compare less-than", so we have to compare greater,
  // compare equal, and take everything else to be compared
  // less-then
  __m256i is_gt = _mm256_cmpgt_epi16(s, median);
  __m256i is_eq = _mm256_cmpeq_epi16(s, median);

  __m256i to_add = _mm256_set1_epi16(-1);
  // Really want an epi16 version of this, but the cmpgt and
  // cmplt functions set their epi16 parts to 0xff or 0x0,
  // so treating everything as epi8 works the same
  to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(1), is_gt);
  to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(0), is_eq);

  // Don't add anything to the channels which are masked out
  to_add = _mm256_and_si256(to_add, mask);

  accum = _mm256_add_epi16(accum, to_add);

  // if the accumulator is >10, add one to the median and
  // set the accumulator to zero. if the accumulator is
  // <-10, subtract one from the median and set the
  // accumulator to zero
  is_gt = _mm256_cmpgt_epi16(accum, _mm256_set1_epi16(acclimit));
  __m256i is_lt =
    _mm256_cmpgt_epi16(_mm256_sign_epi16(accum, _mm256_set1_epi16(-1 * acclimit)), _mm256_set1_epi16(acclimit));

  to_add = _mm256_setzero_si256();
  to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(1), is_gt);
  to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(-1), is_lt);

  // Don't add anything to the channels which are masked out
  to_add = _mm256_and_si256(to_add, mask);

  median = _mm256_adds_epi16(median, to_add);

  // Reset the unmasked channels that were >10 or <-10 to zero, leaving the others unchanged
  __m256i need_reset = _mm256_or_si256(is_lt, is_gt);
  need_reset = _mm256_and_si256(need_reset, mask);
  accum = _mm256_blendv_epi8(accum, _mm256_setzero_si256(), need_reset);
}

// Perform the division of __m256i with a const int
inline __m256i _mm256_div_epi16 (const __m256i va, const int b)
{
    __m256i vb = _mm256_set1_epi16(32768 / b);
    return _mm256_mulhrs_epi16(va, vb);
}


}
#endif // READOUT_SRC_WIB2_TPG_UTILSAVX2_HPP_
