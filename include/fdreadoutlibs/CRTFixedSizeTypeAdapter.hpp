#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTFIXEDSIZETYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTFIXEDSIZETYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/CRTFixedSizeFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

/**
 * @brief For CRT the numbers are different.
 * DAQEthHeader (16 bytes) + Header(16 bytes) + 64 ch * 4 bytes/ch  = 288[Bytes]
 * */
const constexpr std::size_t kCRTFixedSize = 288; 
struct CRTFixedSizeTypeAdapter
{
  using FrameType = dunedaq::fddetdataformats::CRTFixedSizeFrame;
  // data
  char data[kCRTFixedSize];
  // comparable based on first timestamp
  bool operator<(const CRTFixedSizeTypeAdapter& other) const
  {
    auto thisptr = reinterpret_cast<const FrameType*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const FrameType*>(&data)->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT
    frame->set_timestamp(ts);
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset = 2048*/ ) // NOLINT(build/unsigned)
  {
    auto wef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    wef->set_timestamp(first_timestamp);
  }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
        df->daq_header.crate_id = crate_id;
        df->daq_header.slot_id = slot_id;
        df->daq_header.stream_id = stream_id;
      }
  }

  void fake_adc_pattern(int channel) {
    auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT
    // Set the ADC to the int16 maximum value 
    frame->set_adc(channel, 8191); 
  }

  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT
  {
    // Set error bits in header
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + kCRTFixedSize); // NOLINT
  }

  size_t get_payload_size() { return kCRTFixedSize; }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return kCRTFixedSize; }

  static const constexpr size_t fixed_payload_size = 288;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kCRT;
  static const constexpr uint64_t expected_tick_difference = 16; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct dunedaq::detdataformats::DAQEthHeader) == 16,
              "Check your assumptions on DAQEthHeader");

} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTFIXEDSIZETYPEADAPTER_HPP_ */