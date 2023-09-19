#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSPFRAMETYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSPFRAMETYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/SSPTypes.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

const constexpr std::size_t kSSPFrameSize = 1012;
struct SSPFrameTypeAdapter
{
  using FrameType = SSPFrameTypeAdapter;

  // header
  fddetdataformats::ssp::EventHeader header;

  // data
  char data[kSSPFrameSize];

  // comparable based on start timestamp
  bool operator<(const SSPFrameTypeAdapter& other) const
  {
    return this->get_first_timestamp() < other.get_first_timestamp() ? true : false;
  }

  uint64_t get_timestamp() const  // NOLINT(build/unsigned)
  {
    return get_first_timestamp();
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    auto ehptr = &header;
    unsigned long ts = 0; // NOLINT(runtime/int)
    for (unsigned int iword = 0; iword <= 3; ++iword) {
      ts += ((unsigned long)(ehptr->timestamp[iword])) << 16 * iword; //NOLINT(runtime/int)
    }
    return ts;
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    uint64_t bitmask = (1 << 16) - 1; // NOLINT(build/unsigned)
    for (unsigned int iword = 0; iword <= 3; ++iword) {
      header.timestamp[iword] = static_cast<uint16_t>((ts & bitmask)); // NOLINT(build/unsigned)
      ts = ts >> 16;
    }
  }

  void fake_timestamps(uint64_t /*first_timestamp*/, uint64_t /*offset = 25*/) // NOLINT(build/unsigned)
  {
    // tp.time_start = first_timestamp;
  }

  void fake_geoid(uint16_t /*crate_id*/, uint16_t /*slot_id*/, uint16_t /*link_id*/) {
  } 

  void fake_adc_pattern(int /*channel*/) {
  }

  FrameType* begin() { return this; }

  FrameType* end() { return (this + 1); } // NOLINT

  size_t get_payload_size() {
    return kSSPFrameSize;
  }

  size_t get_num_frames() {
    return 1;
  }

  size_t get_frame_size() {
    return kSSPFrameSize;
  }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kUnknown;
  static const constexpr uint64_t expected_tick_difference = 25; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct SSPFrameTypeAdapter) == sizeof(fddetdataformats::ssp::EventHeader) + kSSPFrameSize,
              "Check your assumptions on SSPFrameTypeAdapter");


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq
#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSPFRAMETYPEADAPTER_HPP_ */
