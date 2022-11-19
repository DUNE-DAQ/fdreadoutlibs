#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEAMCFRAMETYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEAMCFRAMETYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "detdataformats/tde/TDE16Frame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

const constexpr std::size_t kTDEAMCFrameSize = 1 * sizeof(dunedaq::detdataformats::tde::TDE16Frame);
struct TDEAMCFrameTypeAdapter

{
  using FrameType = dunedaq::detdataformats::tde::TDE16Frame;

  char data[kTDEAMCFrameSize];

  bool operator<(const TDEAMCFrameTypeAdapter
& other) const
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
    reinterpret_cast<FrameType*>(&data)->get_tde_header()->set_timestamp(ts); // NOLINT
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t) // NOLINT(build/unsigned)
  {
    auto tdef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 64; ++i) {
      auto wfh = const_cast<dunedaq::detdataformats::tde::TDE16Header*>(tdef->get_tde_header());
      wfh->set_timestamp(first_timestamp);
      tdef++;
    }
  }

  void fake_frame_errors(std::vector<uint16_t>* fake_errors) // NOLINT(build/unsigned)
  {
    auto tdef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 64; ++i) {
      tdef->set_tde_errors((*fake_errors)[i]);
      tdef++;
    }
  }

  size_t get_payload_size() { return 64 * sizeof(FrameType); }

  size_t get_num_frames() { return 64; }

  size_t get_frame_size() { return sizeof(FrameType); }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + kTDEAMCFrameSize); // NOLINT
  }

  // static const constexpr size_t fixed_payload_size = 5568;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kTDE_AMC;
  static const constexpr uint64_t expected_tick_difference = 1000; // NOLINT(build/unsigned)

};


} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEAMCFRAMETYPEADAPTER_HPP_ */