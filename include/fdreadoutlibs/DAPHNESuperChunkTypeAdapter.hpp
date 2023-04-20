#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/DAPHNEFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

/**
 * @brief For DAPHNE the numbers are different.
 * 12[DAPHNE frames] x 584[Bytes] = 7008[Bytes]
 * */
const constexpr std::size_t kDAPHNESuperChunkSize = 7008; // for 12: 7008
struct DAPHNESuperChunkTypeAdapter
{
  using FrameType = dunedaq::detdataformats::daphne::DAPHNEFrame;
  // data
  char data[kDAPHNESuperChunkSize];
  // comparable based on first timestamp
  bool operator<(const DAPHNESuperChunkTypeAdapter& other) const
  {
    auto thisptr = reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data)->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    auto frame = reinterpret_cast<dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data); // NOLINT
    frame->header.timestamp_wf_1 = ts;
    frame->header.timestamp_wf_2 = ts >> 32;
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 25) // NOLINT(build/unsigned)
  {
    uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)
    for (unsigned int i = 0; i < 12; ++i) {
      auto df = reinterpret_cast<dunedaq::detdataformats::daphne::DAPHNEFrame*>(((uint8_t*)(&data)) + i * 584); // NOLINT
      df->header.timestamp_wf_1 = ts_next;
      df->header.timestamp_wf_2 = ts_next >> 32;
      ts_next += offset;
    }
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
    return reinterpret_cast<FrameType*>(data + kDAPHNESuperChunkSize); // NOLINT
  }

  size_t get_payload_size() { return 7008; }

  size_t get_num_frames() { return 12; }

  size_t get_frame_size() { return 584; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kDAPHNE;
  static const constexpr uint64_t expected_tick_difference = 16; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct DAPHNESuperChunkTypeAdapter) == kDAPHNESuperChunkSize,
              "Check your assumptions on DAPHNESuperChunkTypeAdapter");


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_ */