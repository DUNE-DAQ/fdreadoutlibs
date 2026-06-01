#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/DAPHNEEthFrame.hpp"


namespace dunedaq::fdreadoutlibs::types {

  constexpr std::size_t kDAPHNEEthSize = 1864;

  class DAPHNEEthTypeAdapter : public TypeAdapter<fddetdataformats::DAPHNEEthFrame,
					 NumFrames{1},
    ExpectedTickDifference{1},
						     daqdataformats::SourceID::Subsystem::kDetectorReadout,
							    daqdataformats::FragmentType::kDAPHNEEth> {
  public:
    void fake_timestamps(uint64_t ts, uint64_t /*offset*/= 2048) {
      reinterpret_cast<fddetdataformats::DAPHNEEthFrame*>(&data)->set_timestamp(ts);
    }
  };

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHTYPEADAPTER_
