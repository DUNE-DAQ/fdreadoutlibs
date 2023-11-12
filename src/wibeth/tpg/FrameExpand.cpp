/**
 * @file WIBFrameProcessor.hpp WIB specific Task based raw processor
 * @author Philip Rodrigues (rodriges@fnal.gov)
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/wibeth/tpg/FrameExpand.hpp"

namespace swtpg_wibeth {

//==============================================================================
// Print a 256-bit register interpreting it as packed 8-bit values
void
print256(__m256i var)
{
  uint8_t* val = (uint8_t*)&var; // NOLINT
  for (int i = 0; i < 32; ++i)
    printf("%02x ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 256-bit register interpreting it as packed 8-bit values, in reverse order
void
printr256(__m256i var)
{
  uint8_t* val = (uint8_t*)&var; // NOLINT
  for (int i = 31; i >= 0; --i)
    printf("%02x ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 256-bit register interpreting it as packed 16-bit values
void
print256_as16(__m256i var)
{
  uint16_t* val = (uint16_t*)&var; // NOLINT
  for (int i = 0; i < 16; ++i)
    printf("%04x ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 256-bit register interpreting it as packed 16-bit values
void
print256_as16_dec(__m256i var)
{
  int16_t* val = (int16_t*)&var; // NOLINT
  for (int i = 0; i < 16; ++i)
    printf("%+6d ", val[i]); // NOLINT(runtime/output_format)
}

//==============================================================================
// Print a 256-bit register interpreting it as packed 32-bit values
void
print256_as32_dec(__m256i var)
{
  int32_t* val = (int32_t*)&var; // NOLINT
  for (int i = 0; i < 8; ++i)
    printf("%+6d ", val[i]); // NOLINT(runtime/output_format)
}

} // namespace swtpg_wibeth
