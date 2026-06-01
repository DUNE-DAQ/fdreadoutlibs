#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/DAPHNEEthStreamFrame.hpp"


namespace dunedaq::fdreadoutlibs::types {

  constexpr std::size_t kDAPHNEEthStreamNumFrames = 1;
  constexpr std::size_t kDAPHNEEthStreamSize = kDAPHNEEthStreamNumFrames * sizeof(dunedaq::fddetdataformats::DAPHNEEthStreamFrame); 

  class DAPHNEEthStreamTypeAdapter : public TypeAdapter<fddetdataformats::DAPHNEEthStreamFrame,
								  NumFrames{1},
    ExpectedTickDifference{280},
						     daqdataformats::SourceID::Subsystem::kDetectorReadout,
							    daqdataformats::FragmentType::kDAPHNEEthStream> {
  public:
    void fake_timestamps(uint64_t first_timestamp,
			 uint64_t offset = DAPHNEEthStreamTypeAdapter::expected_tick_difference) {
      uint64_t ts_next = first_timestamp;
      
      for (auto i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<fddetdataformats::DAPHNEEthStreamFrame*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->set_timestamp(ts_next);
        ts_next += offset;
      }
    }

    void fake_adc_pattern(int channel) {
      auto frame = reinterpret_cast<fddetdataformats::DAPHNEEthStreamFrame*>(&data); // NOLINT

      // Set the ADC for the first sample to the 14-bit max value
      frame->set_adc(0, channel, 0x3FFF);
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_
