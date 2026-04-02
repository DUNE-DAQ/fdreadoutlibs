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

  void fake_timestamps(fddetdataformats::DAPHNEEthFrame* frame,
		       uint64_t first_timestamp, uint64_t /*ignored*/) {
    frame->set_timestamp(first_timestamp);
  }

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_
