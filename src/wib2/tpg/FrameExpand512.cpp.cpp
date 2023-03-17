/**
 * @file WIBFrameProcessor.hpp WIB specific Task based raw processor
 * @author Adam Abed Abud & Gordon Crone (adam.abed.abud@cern.ch)
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/wib2/tpg/FrameExpand512.hpp"

namespace swtpg_wib2 {



//==============================================================================
// Print a 512-bit register interpreting it as packed 8-bit values
void
print512(__m512i var)
{
  uint8_t* val = (uint8_t*)&var; // NOLINT
  for (int i = 0; i < 64; ++i)
    printf("%02x ", val[i]); // NOLINT(runtime/output_format)
}


//==============================================================================
// Print a 512-bit register interpreting it as packed 16-bit values
void
print512_as16(__m512i var)
{
  uint16_t* val = (uint16_t*)&var; // NOLINT
  for (int i = 0; i < 32; ++i)
    printf("%04x ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 512-bit register interpreting it as packed 32-bit values
void
print512_as32(__m512i var)
{
  uint32_t* val = (uint32_t*)&var; // NOLINT
  for (int i = 0; i < 16; ++i)
    printf("%08x ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 512-bit register interpreting it as packed 16-bit values
void
print512_as16_dec(__m512i var)
{
  int16_t* val = (int16_t*)&var; // NOLINT
  for (int i = 0; i < 32; ++i)
    printf("%+6d ", val[i]); // NOLINT(runtime/output_format)
}

} // namespace swtpg_wib2
