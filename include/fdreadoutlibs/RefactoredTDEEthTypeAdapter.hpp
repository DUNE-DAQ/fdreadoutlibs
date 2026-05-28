#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDTDEETHTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDTDEETHTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/SourceID.hpp"
#include "daqdataformats/FragmentHeader.hpp"  // for FragmentType
#include "fddetdataformats/TDEEthFrame.hpp"

namespace dunedaq::fdreadoutlibs::types {

  using RefactoredTDEEthTypeAdapter = TypeAdapter<fddetdataformats::TDEEthFrame,
						  1,
						  2000
						  daqdataformats::SourceID::Subsystem::kDetectorReadout,
						  daqdataformats::FragmentType::kTDEEth>;

  void fake_timestamps(fddetdataformats::TDEEthFrame* frame,
		       uint64_t first_timestamp, uint64_t /*ignored*/) {
    frame->set_timestamp(first_timestamp);
  }

  void fake_adc_pattern(fddetdataformats::TDEEthFrame* frame, int channel) {
    frame->set_adc(channel, 0, 16383);
  }

} // namespace dunedaq::fdreadoutlibs::types

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_REFACTOREDTDEETHTYPEADAPTER_
