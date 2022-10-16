#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_PROTOWIBSUPERCHUNKTYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_PROTOWIBSUPERCHUNKTYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "detdataformats/wib/WIBFrame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

/**
 * @brief SuperChunk concept: The FELIX user payloads are called CHUNKs.
 * There is mechanism in firmware to aggregate WIB frames to a user payload
 * that is called a SuperChunk. Default mode is with 12 frames:
 * 12[WIB frames] x 464[Bytes] = 5568[Bytes]
 */
const constexpr std::size_t kProtoWIBSuperChunkSize = 5568; // for 12: 5568
struct ProtoWIBSuperChunkTypeAdapter
{
  using FrameType = dunedaq::detdataformats::wib::WIBFrame;

  // data
  char data[kProtoWIBSuperChunkSize];
  // comparable based on first timestamp
  bool operator<(const ProtoWIBSuperChunkTypeAdapter& other) const
  {
    // auto thisptr = reinterpret_cast<const dunedaq::detdataformats::WIBHeader*>(&data);        // NOLINT
    // auto otherptr = reinterpret_cast<const dunedaq::detdataformats::WIBHeader*>(&other.data); // NOLINT
    return this->get_first_timestamp() < other.get_first_timestamp();
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const dunedaq::detdataformats::wib::WIBFrame*>(&data)->get_wib_header()->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(&data)->get_wib_header()->set_timestamp(ts); // NOLINT
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 25) // NOLINT(build/unsigned)
  {
    uint64_t ts_next = first_timestamp;                                               // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(((uint8_t*)(&data))); // NOLINT
    for (unsigned int i = 0; i < 12; ++i) {
      auto wfh = const_cast<dunedaq::detdataformats::wib::WIBHeader*>(wf->get_wib_header());
      wfh->set_timestamp(ts_next);
      ts_next += offset;
      wf++;
    }
  }

  void fake_frame_errors(std::vector<uint16_t>* fake_errors) // NOLINT(build/unsigned)
  {
    auto wf = reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 12; ++i) {
      wf->set_wib_errors((*fake_errors)[i]);
      wf++;
    }
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + kProtoWIBSuperChunkSize); // NOLINT
  }

  size_t get_payload_size() { return 5568; }

  size_t get_num_frames() { return 12; }

  size_t get_frame_size() { return 464; }
  
  static const constexpr size_t fixed_payload_size = 5568;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kProtoWIB;
  static const constexpr uint64_t expected_tick_difference = 25; // 2 MHz@50MHz clock // NOLINT(build/unsigned)
};

static_assert(sizeof(struct dunedaq::detdataformats::wib::WIBFrame)*12 == kProtoWIBSuperChunkSize,
              "Check your assumptions on ProtoWIBSuperChunkTypeAdapter");

} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_PROTOWIBSUPERCHUNKTYPEADAPTER_HPP_ */