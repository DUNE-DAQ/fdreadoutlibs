#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/DAPHNEFrame.hpp"


namespace dunedaq::fdreadoutlibs::types {

  constexpr std::size_t kDAPHNEFrameSize = 1864;
  constexpr std::size_t kDAPHNENumFrames = 3;

  class DAPHNESuperChunkTypeAdapter : public TypeAdapter<fddetdataformats::DAPHNEFrame,
					 NumFrames{3},
    ExpectedTickDifference{1024},
						     daqdataformats::SourceID::Subsystem::kDetectorReadout,
							    daqdataformats::FragmentType::kDAPHNE> {
  public:
    void fake_timestamps(uint64_t ts, uint64_t offset = DAPHNESuperChunkTypeAdapter::expected_tick_difference) {
      uint64_t ts_next = ts; // NOLINT(build/unsigned)

      for (auto i = 0; i < get_num_frames(); ++i) {
	auto df = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEFrame*>(data + i * get_frame_size()); // NOLINT
	df->set_timestamp(ts_next);
	ts_next += offset;
      }
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESUPERCHUNKTYPEADAPTER_
