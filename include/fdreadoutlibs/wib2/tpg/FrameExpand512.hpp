/**
 * @file frame_expand.h WIB specific frame expansion
 * @author Philip Rodrigues (rodriges@fnal.gov)
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND512_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND512_HPP_

#include "TPGConstants_wib2.hpp"
#include "detdataformats/wib2/WIB2Frame.hpp"
#include "fdreadoutlibs/ProtoWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"


#include <array>
#include <immintrin.h>

namespace swtpg_wib2 {

// A little wrapper around an array of 512-bit registers, so that we
// can explicitly access it as an array of 512-bit registers or as an
// array of uint16_t
template<size_t N>
class RegisterArray512
{
public:
  // RegisterArray() = default;

  // RegisterArray(RegisterArray& other)
  // {
  //     memcpy(m_array, other.m_array, N*sizeof(uint16_t)); NOLINT(build/unsigned)
  // }

  // RegisterArray(RegisterArray&& other) = default;

  // Get the value at the ith position as a 512-bit register
  inline __m512i ymm(size_t i) const
  {
    return _mm512_loadu_epi16(reinterpret_cast<const __m512i*>(m_array) + i); // NOLINT
  }
  inline void set_ymm(size_t i, __m512i val)
  {
    _mm512_storeu_epi16(reinterpret_cast<__m512i*>(m_array) + i, val); // NOLINT
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
  alignas(64) uint16_t __restrict__ m_array[N * 32]; // NOLINT(build/unsigned)
};

typedef RegisterArray512<swtpg_wib2::NUM_REGISTERS_PER_FRAME> FrameRegisters512;

typedef RegisterArray512<swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::FRAMES_PER_MSG> MessageRegisters512;



//==============================================================================
// Print a 512-bit register interpreting it as packed 8-bit values
void
print512(__m512i var);

//==============================================================================
// Print a 512-bit register interpreting it as packed 16-bit values
void
print512_as16(__m512i var);
//==============================================================================
// Print a 512-bit register interpreting it as packed 32-bit values
void
print512_as32(__m512i var);

//==============================================================================
// Print a 512-bit register interpreting it as packed 16-bit values
void
print512_as16_dec(__m512i var);

//==============================================================================
inline __m512i unpack_one_register512(const dunedaq::detdataformats::wib2::WIB2Frame::word_t* first_word)
{
    __m512i reg=_mm512_loadu_si512((__m512i*)first_word);
    // printf("Input:      ");
    // print512_as32(reg);
    // printf("\n");

    // The register initially contains 36-and-a-bit 14-bit ADCs, but
    // we only have space for 32 after expansion, so the last 2 32-bit
    // words are unused. Copy words 7 and 3 so they appear twice, and move the
    // later words down one
//    15        14        13       12       11         10       9        8
//  8|14|10  4|14|14    14|14|4  10|14|8  6|14|12  2|14|14|2  12|14|6  8|14|10
//                      32 31  30   29  28   27  26  25 24  23   22  21  20
//     7         6        5        4          3          2       1        0
//  4|14|14   14|14|4  10|14|8   6|14|12   2|14|14|2  12|14|6  8|14|10  4|14|14
//19  18 17   16 15  14   13  12   11    10   9  8   7    6   5   4    3   2  1

    __m512i idx=_mm512_set_epi32(13,12,11,10,10,9,8,7,6,5,4,3,3,2,1,0);
    __m512i shuf1=_mm512_permutexvar_epi32(idx,reg);
    // printf("shuf1:      ");
    // print512_as32(shuf1);
    // printf("\n");

    // Each 32-bit word contains at least one full 14-bit ADC. Shift
    // the words by variable amounts s.t. the high 16 bits of each
    // word contains a 14-bit ADC at the right place (with the two
    // high bits still needing to be masked to zero). That result is
    // in `high_half`
     
    // The amounts by which we shift each 32-bit word
    __m512i count1=_mm512_set_epi32(
      12, 8, 4, 0, 14, 10, 6, 2,  12, 8, 4, 0, 14, 10, 6, 2);
    __m512i high_half=_mm512_sllv_epi32(shuf1, count1);
    // Mask out the low 16 bits, and the high two bits in the high half
    __m512i high_half_mask=_mm512_set1_epi32(0x3fff0000u);

    high_half=_mm512_and_si512(high_half, high_half_mask);
    // printf("high_half:  ");
    // print512_as32(high_half);
    // printf("\n");

    //------------------------------------------------------------------
    // Now we start the process of setting the low 16 bits of each
    // word to the right value. This is trickier because now the bits
    // are spread across two words. First, left-shift each word so
    // that the higher bits of the ADC are in the right place
    __m512i count2=_mm512_set_epi32(10, 6, 2, 0, 12, 8, 4, 0,
		    10, 6, 2, 0, 12, 8, 4, 0);
    __m512i shift2=_mm512_sllv_epi32(shuf1, count2);
    // printf("shift2:     ");
    // print512_as32(shift2);
    // printf("\n");

    // Next, permute the register so that the words containing the low
    // bits of the ADCs we want are in the same positions as the words
    // containing the corresponding high bits. This just amounts to
    // moving the words down by one
    __m512i idx2=_mm512_set_epi32(12, 11, 10, 10, 9, 8, 7, 7, 5,
			    4, 3, 2, 2, 1, 0, 0);
    __m512i shuf2=_mm512_permutexvar_epi32(idx2, reg);
    // printf("shuf2:      ");
    // print512_as32(shuf2);
    // printf("\n");

    // Shift each word right by the amount that brings those low bits
    // into the right place, putting the result in `shift3`
    __m512i count3=_mm512_set_epi32(22, 26, 30, 0, 20, 24, 28, 0,
	    22, 26, 30, 0, 20, 24, 28, 0);
    __m512i shift3=_mm512_srlv_epi32(shuf2, count3);
    // printf("shift3:     ");
    // print512_as32(shift3);
    // printf("\n");

    // OR together the registers containing the high and low bits of
    // the ADCs. At this point, the low 16 bits of each word should
    // contain the 14 bits of the ADCs in the right place (with the
    // two high bits still needing to be masked out)
    __m512i low_half=_mm512_or_si512(shift2, shift3);
    // Mask out the high 16 bits, and the high two bits in the high half
    __m512i low_half_mask=_mm512_set1_epi32(0x3fffu);
    low_half=_mm512_and_si512(low_half, low_half_mask);
    // printf("low_half:   ");
    // print512_as32(low_half);
    // printf("\n");

    // Nearly there... Now we OR together the low and high halves
    __m512i both=_mm512_or_si512(low_half, high_half);
    // zero out the slots where we want to put the 16th and 32nd values
    both=_mm512_andnot_si512(_mm512_set_epi32(0, 0, 0, 0xffffu, 0, 0, 0, 0,
			    0, 0, 0, 0xffffu, 0, 0, 0, 0), both);
    // printf("both:       ");
    // print512_as32(both);
    // printf("\n");

    // We just missed the 16th and 32nd values, and the lw 16 bits of the
    // ?? 8th ??
    // word are available, so shuffle it around to put it there
    __m512i shift4=_mm512_srli_epi32(reg, 18);
    // Mask so that's the only nonzero thing
    shift4=_mm512_and_si512(_mm512_set_epi32(0, 0, 0x3fffu, 0, 0, 0, 0, 0,
			    0, 0x3fffu, 0, 0, 0, 0, 0, 0), shift4);
    // Move the word containing the value we want into the position we want
    __m512i idx3=_mm512_set_epi32(0, 0, 0, 13, 0, 0, 0, 0,
		    0, 0, 0, 6, 0, 0, 0, 0);
    __m512i shuf3=_mm512_permutexvar_epi32(idx3, shift4);

    both=_mm512_or_si512(both, shuf3);
    // printf("both':      ");
    // print512_as32(both);
    // printf("\n");

    return both;
}




// Expand 14-bit ADCs to 16-bits using the WIB2 format
inline void
expand_wib2_adcs_512(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* __restrict__ ucs,
                            swtpg_wib2::MessageRegisters512* __restrict__ register_array, 
                            int registers_selection
                            )
{
  std::cout << "expand_wib2_adcs() registers_selection=" << registers_selection << std::endl;
  for (size_t iframe = 0; iframe < swtpg_wib2::FRAMES_PER_MSG; ++iframe) {
    const dunedaq::detdataformats::wib2::WIB2Frame* frame =
      reinterpret_cast<const dunedaq::detdataformats::wib2::WIB2Frame*>(ucs) + iframe; // NOLINT

    //std::cout  << " Unpacking frame " << iframe << ":\n";
    for (size_t iblock = 0; iblock < swtpg_wib2::NUM_REGISTERS_PER_FRAME ; ++iblock) {
      //std::cout << " block " << iblock << " from offset "
      //          << 14*(iblock+registers_selection*swtpg_wib2::NUM_REGISTERS_PER_FRAME)
      //          << " (0x" << std::hex << 14*(iblock+registers_selection*swtpg_wib2::NUM_REGISTERS_PER_FRAME) << ")"
      //          << " to offset " << std::dec << iframe + iblock * swtpg_wib2::FRAMES_PER_MSG
      //          << " (0x" << std::hex << iframe + iblock * swtpg_wib2::FRAMES_PER_MSG << ")"
      //          ;//<< std::dec << std::endl;

      register_array->set_ymm(
        iframe + iblock * swtpg_wib2::FRAMES_PER_MSG,
        swtpg_wib2::unpack_one_register512(
          frame->adc_words+14*(iblock+registers_selection*swtpg_wib2::NUM_REGISTERS_PER_FRAME)));
      //std::cout << "  first word " << register_array->uint16(iframe + iblock * swtpg_wib2::FRAMES_PER_MSG, 0)
	    //  << std::dec << std::endl;
      //print512_as16(register_array->ymm(iframe + iblock * swtpg_wib2::FRAMES_PER_MSG));
      //std::cout << std::endl;
    }
  }
}



} // namespace swtpg_wib2

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_FRAMEEXPAND512_HPP_


