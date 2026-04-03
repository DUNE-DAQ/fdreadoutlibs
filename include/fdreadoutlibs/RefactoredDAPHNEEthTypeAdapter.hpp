#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/DAPHNEEthFrame.hpp"


namespace dunedaq::fdreadoutlibs::types {

  using RefactoredDAPHNEEthTypeAdapter = TypeAdapter<fddetdataformats::DAPHNEEthFrame,
						     1,
						     1,
						     daqdataformats::SourceID::Subsystem::kDetectorReadout,
						     daqdataformats::FragmentType::kDAPHNEEth>;

  template <>
  void fake_timestamps<dunedaq::fddetdataformats::DAPHNEEthFrame>(
								  dunedaq::fddetdataformats::DAPHNEEthFrame* frame,
								  uint64_t ts,
								  uint64_t) {
    frame->set_timestamp(ts);
  }
  
} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_
