#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETHTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/WIBEthFrame.hpp"

namespace dunedaq::fdreadoutlibs::types {

  /**
   * @brief For WIBEth the numbers are different.
   * Header + (64 channels * 64 time slices) = 7200[Bytes]
   * */
  constexpr std::size_t kDUNEWIBEthSize = 7200;

  class DUNEWIBEthTypeAdapter : public TypeAdapter<fddetdataformats::WIBEthFrame,
				      NumFrames{1},
    ExpectedTickDifference{2048},
    daqdataformats::SourceID::Subsystem::kDetectorReadout,
    daqdataformats::FragmentType::kWIBEth> {
  public:
    static constexpr size_t fixed_payload_size = 7200;
    static constexpr uint64_t samples_per_frame = 64; 
    static constexpr float samples_tick_difference = 32; 

    void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset*/ = 2048) {
      auto* wf = reinterpret_cast<fddetdataformats::WIBEthFrame*>(data);
      wf->set_timestamp(first_timestamp);
    }

    void fake_adc_pattern(int channel) {
      auto* wf = reinterpret_cast<fddetdataformats::WIBEthFrame*>(data);
      wf->set_adc(channel, 0, 16383);
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETHTYPEADAPTER_
