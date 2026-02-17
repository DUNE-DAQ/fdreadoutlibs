#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTBERNTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTBERNTYPEADAPTER_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/CRTBernFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

const constexpr std::size_t kCRTBernFrameSize = sizeof(dunedaq::fddetdataformats::CRTBernFrame);

struct CRTBernTypeAdapter
{
  using FrameType = dunedaq::fddetdataformats::CRTBernFrame;
  // data
  char data[kCRTBernFrameSize];
  // comparable based on first timestamp
  bool operator<(const CRTBernTypeAdapter& other) const
  {
    auto thisptr = reinterpret_cast<const FrameType*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const FrameType*>(&data)->get_timestamp(); // NOLINT
  }

  void set_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT
    frame->set_timestamp(ts);
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset*/= expected_tick_difference) // NOLINT(build/unsigned)
  {
    set_timestamp(first_timestamp);
  }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      auto df = reinterpret_cast<FrameType*>(reinterpret_cast<uint8_t*>(&data));
      df->daq_header.slot_id = slot_id;
      df->daq_header.stream_id = stream_id;
      df->daq_header.crate_id = crate_id;
      }

  void fake_adc_pattern(int channel, int time_sample = 0) {
    auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT
    // CRT frames don't use time samples in the same way as TPC, only placeholder
    frame->set_adc(channel,0xbeef);
  }

  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT
  {
    // Set frame error bits in header
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + kCRTBernFrameSize); // NOLINT
  }

  size_t get_payload_size() { return get_num_frames() * get_frame_size(); }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return kCRTBernFrameSize; }

  static const constexpr size_t fixed_payload_size = kCRTBernFrameSize;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kCRTBern;
  static const constexpr uint64_t expected_tick_difference = 1; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct CRTBernTypeAdapter) == kCRTBernFrameSize,
              "Check your assumptions on CRTBernTypeAdapter");


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTBERNTYPEADAPTER_ */
