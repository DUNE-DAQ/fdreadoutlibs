#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/DAPHNEEthStreamFrame.hpp"



namespace dunedaq::fdreadoutlibs::types {

  /**                                                                                                                       
   */
  const constexpr std::size_t kDAPHNEEthStreamNumFrames = 1;
  const constexpr std::size_t kDAPHNEEthStreamSize = kDAPHNEEthStreamNumFrames * sizeof(dunedaq::fddetdataformats::DAPHNEEthStreamFrame); 

  struct DAPHNEEthStreamTypeAdapter {

    using FrameType = dunedaq::fddetdataformats::DAPHNEEthStreamFrame;

    char data[kDAPHNEEthStreamSize];

    // comparable based on first timestamp
    bool operator<(const DAPHNEEthStreamTypeAdapter& other) const
    {
      auto thisptr = reinterpret_cast<const FrameType*>(&data);        // NOLINT
      auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT
      return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
    }

    uint64_t get_timestamp() const { // NOLINT(build/unsigned)
      return reinterpret_cast<const FrameType*>(&data)->get_timestamp(); // NOLINT
    }

    void set_timestamp(uint64_t ts) // NOLINT(build/unsigned)
    {
      auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT                  
      frame->set_timestamp(ts);
    }

    void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 280) // NOLINT(build/unsigned)                          
    {
      uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)                                                         
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
        df->set_timestamp(ts_next);
        ts_next += offset;
      }
    }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<FrameType*>((reinterpret_cast<uint8_t*>(&data)) + i * get_frame_size());
	df->set_geoid(crate_id, slot_id, stream_id);
      }
  }

  void fake_adc_pattern(int channel) {
    // Set the ADC for the first sample to the 14-bit max value 
    auto frame = reinterpret_cast<FrameType*>(&data); // NOLINT
    frame->set_adc(0, channel, 0x3FFF);
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
      return reinterpret_cast<FrameType*>(data + kDAPHNEEthStreamSize); // NOLINT                                  
    }

    constexpr size_t get_payload_size() const { return kDAPHNEEthStreamSize; }

    constexpr size_t get_num_frames() const { return kDAPHNEEthStreamNumFrames; }

    constexpr size_t get_frame_size() const { return sizeof(dunedaq::fddetdataformats::DAPHNEEthStreamFrame); }

    static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
    static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kDAPHNEEthStream;
    static const constexpr uint64_t expected_tick_difference = 280; // NOLINT(build/unsigned)    
  };

  static_assert(sizeof(struct DAPHNEEthStreamTypeAdapter) == kDAPHNEEthStreamSize,
                "Check your assumptions on DAPHNEEthStreamTypeAdapter");


} // namespace dunedaq::fdreadoutlibs::types


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETHSTREAMTYPEADAPTER_ 
