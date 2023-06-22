#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEFRAMETYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEFRAMETYPEADAPTER_HPP_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/TDE16Frame.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

const constexpr std::size_t kTDEFrameSize = sizeof(dunedaq::fddetdataformats::TDE16Frame);
struct TDEFrameTypeAdapter

{
  using FrameType = dunedaq::fddetdataformats::TDE16Frame;

  //char data[kTDEFrameSize];
  FrameType data;

  bool operator<(const TDEFrameTypeAdapter& other) const
  {
    uint64_t ts = data.get_timestamp();
    uint32_t ch = data.get_channel();
    uint64_t ots = other.data.get_timestamp();
    //uint32_t och = 0xffff;
    uint32_t och = other.data.get_channel();;

    return std::tie(ts,ch) < std::tie(ots,och);
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return data.get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    data.set_timestamp(ts); // NOLINT
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t) // NOLINT(build/unsigned)
  {
    data.set_timestamp(first_timestamp); 
  }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t /*link_id*/) {
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
        df->get_daq_header()->crate_id = crate_id;
        df->get_daq_header()->slot_id = slot_id;
        //df->get_daq_header()->stream_id = link_id;
	//df->get_tde_header()->channel = link_id;
      }
  }
  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT(build/unsigned)
  {
  }

  size_t get_payload_size() { return sizeof(FrameType); }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return sizeof(FrameType); }

  FrameType* begin()
  {
    return &data; // NOLINT
  }

  FrameType* end()
  {
    return (reinterpret_cast<FrameType*>(reinterpret_cast<char*>(&data)+kTDEFrameSize)); // NOLINT
  }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kTDE_AMC;
  static const constexpr uint64_t expected_tick_difference = dunedaq::fddetdataformats::ticks_between_adc_samples * dunedaq::fddetdataformats::tot_adc16_samples; // NOLINT(build/unsigned)

};

static_assert(sizeof(dunedaq::fddetdataformats::TDE16Frame) == kTDEFrameSize,
              "Check your assumptions on TDEFrameTypeAdapter");

} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif /* FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEFRAMETYPEADAPTER_HPP_ */
