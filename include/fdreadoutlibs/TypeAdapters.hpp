#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS_

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/FrameConcepts.hpp"

#include <vector>

namespace dunedaq::fdreadoutlibs::types {

  // The following two structs will make it a lot easier to understand
  // the meaning of the relevant values used to instantiate the
  // TypeAdapter template - e.g., "ExpectedTickDifference{2800}"
  // rather than just "2800"
  
  struct NumFrames {
    const int val;
  };

  struct ExpectedTickDifference {
    const uint64_t val;
  };
  
  template <typename FrameType,
	    NumFrames NFrames,
	    ExpectedTickDifference ExpectedTickDiff,
	    daqdataformats::SourceID::Subsystem SubSystem,
	    daqdataformats::FragmentType FragmentType
	    >
  requires fddetdataformats::AdaptableFrameConcept<FrameType>
  struct TypeAdapter {

    static constexpr int s_num_frames { NFrames.val };

    // Naming convention is wrong for consistency with existing code
    static constexpr uint64_t expected_tick_difference { ExpectedTickDiff.val };
    static constexpr daqdataformats::SourceID::Subsystem subsystem { SubSystem };
    static constexpr daqdataformats::FragmentType fragment_type { FragmentType };

    char data[sizeof(FrameType)*s_num_frames];

    bool operator<(const TypeAdapter& other) const {
      auto thisptr = reinterpret_cast<const FrameType*>(&data); // NO LINT
      auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT 

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
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->set_geoid(crate_id, slot_id, stream_id);
      }
    }

    void fake_timestamps(uint64_t first_timestamp, uint64_t offset ) {} // NOLINT(build/unsigned)

    void fake_adc_pattern(int channel) {}

    // Why not make this a reference rather than a pointer?
    void fake_frame_errors(std::vector<uint16_t>* fake_errors) {} // NOLINT
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TYPEADAPTERS
