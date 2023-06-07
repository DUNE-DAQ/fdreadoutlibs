/**
 * @file frame_expand.h WIB specific frame expansion
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND_HPP_

#include "TPGConstants_wib2.hpp"
#include "fddetdataformats/WIB2Frame.hpp"
#include "fdreadoutlibs/ProtoWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"


#include <array>
#include <immintrin.h>

namespace swtpg_wib2 {

// A little wrapper around an array of 256-bit registers, so that we
// can explicitly access it as an array of 256-bit registers or as an
// array of uint16_t
template<size_t N>
class RegisterArray
{
public:
  // RegisterArray() = default;

  // RegisterArray(RegisterArray& other)
  // {
  //     memcpy(m_array, other.m_array, N*sizeof(uint16_t)); NOLINT(build/unsigned)
  // }

  // RegisterArray(RegisterArray&& other) = default;

  // Get the value at the ith position as a 256-bit register
  inline __m256i ymm(size_t i) const
  {
    return _mm256_lddqu_si256(reinterpret_cast<const __m256i*>(m_array) + i); // NOLINT
  }
  inline void set_ymm(size_t i, __m256i val)
  {
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_array) + i, val); // NOLINT
  }
  inline uint16_t uint16(size_t i) const { return m_array[i]; }        // NOLINT(build/unsigned)
  inline void set_uint16(size_t i, uint16_t val) { m_array[i] = val; } // NOLINT(build/unsigned)

  // Access the jth entry in the ith register
  inline uint16_t uint16(size_t i, size_t j) const { return m_array[16 * i + j]; }        // NOLINT(build/unsigned)
  inline void set_uint16(size_t i, size_t j, uint16_t val) { m_array[16 * i + j] = val; } // NOLINT(build/unsigned)

  inline uint16_t* data() { return m_array; }             // NOLINT(build/unsigned)
  inline const uint16_t* data() const { return m_array; } // NOLINT(build/unsigned)

  inline size_t size() const { return N; }

private:
  alignas(32) uint16_t __restrict__ m_array[N * 16]; // NOLINT(build/unsigned)
};

typedef RegisterArray<swtpg_wib2::NUM_REGISTERS_PER_FRAME> FrameRegisters;

typedef RegisterArray<swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::FRAMES_PER_MSG> MessageRegisters;



//==============================================================================
// Print a 256-bit register interpreting it as packed 8-bit values
void
print256(__m256i var);

//==============================================================================
// Print a 256-bit register interpreting it as packed 16-bit values
void
print256_as16(__m256i var);

//==============================================================================
// Print a 256-bit register interpreting it as packed 16-bit values
void
print256_as16_dec(__m256i var);

//==============================================================================
inline __m256i unpack_one_register(const dunedaq::fddetdataformats::WIB2Frame::word_t* first_word)
{
    __m256i reg=_mm256_lddqu_si256((__m256i*)first_word);
    // printf("Input:      ");
    // print256(reg);
    // printf("\n");

    // The register initially contains 18-and-a-bit 14-bit ADCs, but
    // we only have space for 16 after expansion, so the last 32-bit
    // word is unused. Copy word 3 so it appears twice, and move the
    // later words down one
    __m256i idx=_mm256_set_epi32(6, 5, 4, 3, 3, 2, 1, 0);
    __m256i shuf1=_mm256_permutevar8x32_epi32(reg, idx);
    // printf("shuf1:      ");
    // print256(shuf1);
    // printf("\n");

    // Each 32-bit word contains at least one full 14-bit ADC. Shift
    // the words by variable amounts s.t. the high 16 bits of each
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
    // printf("both':      ");
    // print256(both);
    // printf("\n");

    return both;
}




// Expand 14-bit ADCs to 16-bits using the WIB2 format
inline void
expand_wib2_adcs(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* __restrict__ ucs,
                            swtpg_wib2::MessageRegisters* __restrict__ register_array, 
                            int registers_selection
                            )
{
  for (size_t iframe = 0; iframe < swtpg_wib2::FRAMES_PER_MSG; ++iframe) {
    const dunedaq::fddetdataformats::WIB2Frame* frame =
      reinterpret_cast<const dunedaq::fddetdataformats::WIB2Frame*>(ucs) + iframe; // NOLINT

    for (size_t iblock = 0; iblock < swtpg_wib2::NUM_REGISTERS_PER_FRAME ; ++iblock) {
      register_array->set_ymm(iframe + iblock * swtpg_wib2::FRAMES_PER_MSG, swtpg_wib2::unpack_one_register(frame->adc_words+7*(iblock+registers_selection*swtpg_wib2::NUM_REGISTERS_PER_FRAME)));
    }
    

  }
}



} // namespace swtpg_wib2

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND_HPP_
