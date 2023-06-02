#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/DAPHNEStreamFrame.hpp"



namespace dunedaq::fdreadoutlibs::types {

  /**                                                                                                                       
   * @brief For DAPHNE Stream the numbers are similar to DUNE-WIB                                                           
   * 12[DAPHNE frames] x 472[Bytes] = 5664[Bytes]                                                                           
   * */
  const constexpr std::size_t kDAPHNEStreamSuperChunkSize = 5664; // for 12: 5664 

  struct DAPHNEStreamSuperChunkTypeAdapter {

    using FrameType = dunedaq::fddetdataformats::DAPHNEStreamFrame;

    char data[kDAPHNEStreamSuperChunkSize];

    // comparable based on first timestamp
    bool operator<(const DAPHNEStreamSuperChunkTypeAdapter& other) const
    {
      auto thisptr = reinterpret_cast<const FrameType*>(&data);        // NOLINT
      auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT
      return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
    }

    uint64_t get_first_timestamp() const // NOLINT(build/unsigned)                                                          
    {
      return reinterpret_cast<const FrameType*>(&data)->daq_header.get_timestamp(); // NOLINT                                                                                                               
    }

    void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)                                                         
    {
      auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT                  
      frame->daq_header.timestamp_1 = ts;
      frame->daq_header.timestamp_2 = ts >> 32;
    }

    void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 64) // NOLINT(build/unsigned)                          
    {
      uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)                                                         
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
	auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->daq_header.timestamp_1 = ts_next;
	df->daq_header.timestamp_2 = ts_next >> 32;
	ts_next += offset;
      }
    }

    void fake_geoid(uint16_t /*crate_id*/, uint16_t /*slot_id*/, uint16_t /*link_id*/) {
    }

    void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT                                                
    {
      // Set frame error bits in header                                                                                     
    }

    FrameType* begin()
    {
      return reinterpret_cast<FrameType*>(&data[0]); // NOLINT                                                              
    }

    FrameType* end()
    {
      return reinterpret_cast<FrameType*>(data + kDAPHNEStreamSuperChunkSize); // NOLINT                                  
    }

    constexpr size_t get_payload_size() const { return get_num_frames() * get_frame_size(); } // 12*472 -> 5664

    constexpr size_t get_num_frames() const { return 12; }

    constexpr size_t get_frame_size() const { return 472; }

    static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
    static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kDAPHNE;
    static const constexpr uint64_t expected_tick_difference = 64; // NOLINT(build/unsigned)    
  };

  static_assert(sizeof(struct DAPHNEStreamSuperChunkTypeAdapter) == kDAPHNEStreamSuperChunkSize,
		"Check your assumptions on DAPHNESuperChunkTypeAdapter");


} // namespace dunedaq::fdreadoutlibs::types


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNESTREAMSUPERCHUNKTYPEADAPTER_ 
