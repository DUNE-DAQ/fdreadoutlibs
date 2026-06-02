#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_

#include "fdreadoutlibs/TypeAdapters.hpp"

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/DAPHNEStreamFrame.hpp"



namespace dunedaq::fdreadoutlibs::types {

  /**                                                                                                                       
   * @brief For DAPHNE Stream the numbers are similar to DUNE-WIB                                                           
   * 12[DAPHNE frames] x 472[Bytes] = 5664[Bytes]                                                                           
   * */
  constexpr std::size_t kDAPHNEStreamNumFrames = 12;
  constexpr std::size_t kDAPHNEStreamFrameSize = sizeof(dunedaq::fddetdataformats::DAPHNEStreamFrame);
  constexpr std::size_t kDAPHNEStreamSuperChunkSize = kDAPHNEStreamNumFrames * kDAPHNEStreamFrameSize; // for 12: 5664 

  class DAPHNEStreamSuperChunkTypeAdapter : public TypeAdapter<fddetdataformats::DAPHNEStreamFrame,
					    NumFrames{kDAPHNEStreamNumFrames},
    ExpectedTickDifference{64},
						     daqdataformats::SourceID::Subsystem::kDetectorReadout,
							    daqdataformats::FragmentType::kDAPHNEStream> {
  public:

    void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 64) // NOLINT(build/unsigned)                          
    {
      uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)                                                         
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->set_timestamp(ts_next);
        ts_next += offset;
      }
    }

    void fake_geoid(uint16_t /*crate_id*/, uint16_t /*slot_id*/, uint16_t /*link_id*/) {
    }
  };

  static_assert(sizeof(struct DAPHNEStreamSuperChunkTypeAdapter) == kDAPHNEStreamSuperChunkSize,
                "Check your assumptions on DAPHNESuperChunkTypeAdapter");

} // namespace dunedaq::fdreadoutlibs::types


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_ 
