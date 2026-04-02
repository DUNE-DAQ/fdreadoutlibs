#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_

#include "fddetdataformats/DAPHNEEthFrame.hpp"
#include "fdreadoutlibs/TypeAdapters.hpp"

namespace dunedaq::fdreadoutlibs::types {

  using RefactoredDAPHNEEthTypeAdapter = TypeAdapter<fddetdataformats::DAPHNEEthFrame, 1>;

  void fake_timestamps(fddetdataformats::DAPHNEEthFrame* frame,
		       uint64_t first_timestamp, uint64_t /*ignored*/) {
    frame->set_timestamp(first_timestamp);
  }

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDDAPHNEETHTYPEADAPTER_
