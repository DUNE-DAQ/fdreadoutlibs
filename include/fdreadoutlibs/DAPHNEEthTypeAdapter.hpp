#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/DAPHNEEthFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {


/**
 * @brief For DAPHNEEth the numbers are different.
 * Header + (64 channels * 64 time slices) = 233[Bytes]
 * */
const constexpr std::size_t kDAPHNEEthSize = 1864;

struct DAPHNEEthTypeAdapter
{
  using FrameType = dunedaq::fddetdataformats::DAPHNEEthFrame;
  char data[kDAPHNEEthSize];

  bool operator<(const DAPHNEEthTypeAdapter& other) const
  {
    auto thisptr = reinterpret_cast<const dunedaq::fddetdataformats::DAPHNEEthFrame*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const dunedaq::fddetdataformats::DAPHNEEthFrame*>(&other.data); // NOLINT

    return std::forward_as_tuple(thisptr->get_timestamp(), thisptr->get_channel()) < std::forward_as_tuple(otherptr->get_timestamp(), otherptr->get_channel());
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

    void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset*/= 2048 ) // NOLINT(build/unsigned)
  {
    auto wef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    wef->set_timestamp(first_timestamp);
  }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->set_geoid(crate_id, slot_id, stream_id);
      }
  }

  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT
  {
    // Set error bits in header
  }

  void fake_adc_pattern(int /*channel*/) {
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + kDAPHNEEthSize); // NOLINT
  }

  size_t get_payload_size() { return kDAPHNEEthSize; }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return kDAPHNEEthSize; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kDAPHNEEth;
  static const constexpr uint64_t expected_tick_difference = 1; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct dunedaq::fddetdataformats::DAPHNEEthFrame) == kDAPHNEEthSize,
              "Check your assumptions on DAPHNEEthEthTypeAdapter");


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_ */
