#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DUNEWIBETHTYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DUNEWIBETHTYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "detdataformats/wibeth/WIBEthFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

/**
 * @brief For WIBEth the numbers are different.
 * Header + (64 channels * 64 time slices) = 7200[Bytes]
 * */
const constexpr std::size_t kDUNEWIBEthSize = 7200; 
struct DUNEWIBEthTypeAdapter
{
  using FrameType = dunedaq::detdataformats::wibeth::WIBEthFrame;
  // data
  char data[kDUNEWIBEthSize];
  // comparable based on first timestamp
  bool operator<(const DUNEWIBEthTypeAdapter& other) const
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
    return reinterpret_cast<FrameType*>(data + kDUNEWIBEthSize); // NOLINT
  }

  size_t get_payload_size() { return kDUNEWIBEthSize; }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return kDUNEWIBEthSize; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kWIBEth;
  static const constexpr uint64_t expected_tick_difference = 2048; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct dunedaq::detdataformats::wibeth::WIBEthFrame) == kDUNEWIBEthSize,
              "Check your assumptions on DUNEWIBEthTypeAdapter");


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DUNEWIBETHTYPEADAPTER_HPP_ */
