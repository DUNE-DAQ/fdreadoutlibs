/**
 * @file TPGConstants_wibeth.hpp TPG specific constants
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef READOUT_SRC_WIBEth_TPG_TPGCONSTANTS_WIBEth_HPP_
#define READOUT_SRC_WIBEth_TPG_TPGCONSTANTS_WIBEth_HPP_

#include <cstdint>
#include <cstdlib>
#include <limits>

namespace swtpg_wibeth {

const constexpr std::size_t SUPERCHUNK_FRAME_SIZE = 5664;
struct SUPERCHUNK_CHAR_STRUCT
{
  char fragments[SUPERCHUNK_FRAME_SIZE];
};
static_assert(sizeof(struct SUPERCHUNK_CHAR_STRUCT) == 5664,"Check your assumptions on SUPERCHUNK_CHAR_STRUCT"); 

const constexpr std::uint16_t MAGIC = std::numeric_limits<std::uint16_t>::max(); // NOLINT

const constexpr std::int16_t THRESHOLD = 2000;

// How many frames are concatenated in one message. We have 64 time frames 
// AAA: TODO: remove this and use it from detector data formats
const constexpr std::size_t FRAMES_PER_MSG = 64;

// How many AVX2 registers are returned per frame.
// Maximum is 16
const constexpr std::size_t NUM_REGISTERS_PER_FRAME = 4;

// How many bytes are in an AVX2 register
const constexpr std::size_t BYTES_PER_REGISTER = 32;

// How many samples are in a register
const constexpr std::size_t SAMPLES_PER_REGISTER = 16;

// One netio message's worth of channel ADCs after
// expansion: 12 frames per message times 8 registers per frame times
// 32 bytes (256 bits) per register
const constexpr std::size_t ADCS_SIZE = BYTES_PER_REGISTER * NUM_REGISTERS_PER_FRAME * FRAMES_PER_MSG;

} // namespace swtpg_wibeth

#endif // READOUT_SRC_WIBEth_TPG_TPGCONSTANTS_WIBEth_HPP_
