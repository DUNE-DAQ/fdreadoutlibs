
#ifndef WIB2_TPG_CONSTANTS_HPP_
#define WIB2_TPG_CONSTANTS_HPP_


#include <limits>

const constexpr std::size_t SUPERCHUNK_FRAME_SIZE = 5664;
struct SUPERCHUNK_CHAR_STRUCT
{
  char fragments[SUPERCHUNK_FRAME_SIZE];
};
static_assert(sizeof(struct SUPERCHUNK_CHAR_STRUCT) == 5664,"Check your assumptions on SUPERCHUNK_CHAR_STRUCT"); 

const constexpr std::int16_t THRESHOLD = 2000;

// How many frames are concatenated in one netio message
const constexpr std::size_t FRAMES_PER_MSG = 1;

// How many AVX2 registers are returned per frame.
// Maximum is 16
const constexpr std::size_t NUM_REGISTERS_PER_FRAME = 16;

// How many bytes are in an AVX2 register
const constexpr std::size_t BYTES_PER_REGISTER = 32;

// How many samples are in a register
const constexpr std::size_t SAMPLES_PER_REGISTER = 16;

// One netio message's worth of channel ADCs after
// expansion: 12 frames per message times 8 registers per frame times
// 32 bytes (256 bits) per register
const constexpr std::size_t ADCS_SIZE = BYTES_PER_REGISTER * NUM_REGISTERS_PER_FRAME * FRAMES_PER_MSG;


const constexpr std::uint16_t MAGIC = std::numeric_limits<std::uint16_t>::max(); // NOLINT


#endif // WIB2_TPG_CONSTANTS_HPP_
