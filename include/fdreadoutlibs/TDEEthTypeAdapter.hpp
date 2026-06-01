#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/TDEEthFrame.hpp"

namespace dunedaq::fdreadoutlibs::types {

  class TDEEthTypeAdapter : public TypeAdapter<fddetdataformats::TDEEthFrame,
				      NumFrames{1},
    ExpectedTickDifference{2000},
    daqdataformats::SourceID::Subsystem::kDetectorReadout,
    daqdataformats::FragmentType::kTDEEth> {
  public:
    static constexpr size_t fixed_payload_size = 7200;
    static constexpr uint64_t samples_per_frame = 64; 
    static constexpr float samples_tick_difference = 31.25; 

    void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset*/ = 2048) {
      reinterpret_cast<fddetdataformats::TDEEthFrame*>(data)->set_timestamp(first_timestamp);
    }

    void fake_adc_pattern(int channel) {
      reinterpret_cast<fddetdataformats::TDEEthFrame*>(data)->set_adc(channel, 0, 16383);
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHTYPEADAPTER_
