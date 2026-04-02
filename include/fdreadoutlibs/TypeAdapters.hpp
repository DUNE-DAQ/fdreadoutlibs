#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS_

#include "fddetdataformats/FrameConcepts.hpp"

#include <vector>

namespace dunedaq::fdreadoutlibs::types {

  template <typename FrameType>
  requires fddetdataformats::AdaptableFrameConcept<FrameType>
  void fake_timestamps(FrameType* frame, uint64_t first_timestamp, uint64_t offset ) {};

  template <typename FrameType>
  requires fddetdataformats::AdaptableFrameConcept<FrameType>
  void fake_adc_pattern(FrameType* frame, int channel) {};

  template <typename FrameType>
  requires fddetdataformats::AdaptableFrameConcept<FrameType>
  void fake_frame_errors(FrameType* frame, std::vector<uint16_t>* fake_errors) {};


  
  template <typename FrameType,
	    int NumFrames>
  requires fddetdataformats::AdaptableFrameConcept<FrameType>
  struct TypeAdapter {

    static constexpr int s_num_frames { NumFrames };
    
    char data[sizeof(FrameType)*s_num_frames];

    bool operator<(const TypeAdapter& other) const {
      auto thisptr = reinterpret_cast<FrameType*>(&data); // NO LINT
      auto otherptr = reinterpret_cast<FrameType*>(&other.data); // NOLINT 

      return *thisptr < *otherptr;
    }

    uint64_t get_timestamp() const { // NOLINT(build/unsigned)
      return reinterpret_cast<const FrameType*>(&data)->get_timestamp();
    }

    void set_timestamp(uint64_t ts) { // NOLINT(build/unsigned)
      auto frame = reinterpret_cast<FrameType*>(&data);
      frame->set_timestamp(ts);
    }

    FrameType* begin() {
      return reinterpret_cast<FrameType*>(&data[0]); // NOLINT                                                             
    }

    FrameType* end() {
      return reinterpret_cast<FrameType*>(data + sizeof(FrameType)*s_num_frames); // NOLINT
    }

    constexpr size_t get_payload_size() { return sizeof(FrameType)*s_num_frames; }

    constexpr size_t get_num_frames() { return s_num_frames; }

    constexpr size_t get_frame_size() { return sizeof(FrameType); }
    
    void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      for (int i = 0; i < s_num_frames; ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * sizeof(FrameType));
        df->daq_header.crate_id = crate_id;
        df->daq_header.slot_id = slot_id;
        df->daq_header.stream_id = stream_id;
      }
    }

    void fake_timestamps(uint64_t first_timestamp, uint64_t offset ) { // NOLINT(build/unsigned)
      dunedaq::fdreadoutlibs::types::fake_timestamps<FrameType>(reinterpret_cast<FrameType*>(data), first_timestamp, offset);
    }

    void fake_adc_pattern(int channel) {
      dunedaq::fdreadoutlibs::types::fake_adc_pattern<FrameType>(reinterpret_cast<FrameType*>(data), channel);
    }

    // Why not make this a reference rather than a pointer?
    void fake_frame_errors(std::vector<uint16_t>* fake_errors) { // NOLINT
      dunedaq::fdreadoutlibs::types::fake_frame_errors<FrameType>(reinterpret_cast<FrameType*>(data), fake_errors);
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS
