/**
 * @file CRTTypeAdapter.hpp Payload type structures for the DUNE Far Detector
 *
 * This is part of the DUNE DAQ , copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTTYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTTYPEADAPTER_HPP_

#include "iomanager/IOManager.hpp"
#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "fddetdataformats/CRTFrame.hpp"
#include "logging/Logging.hpp"
#include <cstdint> // uint_t types
#include <memory>  // unique_ptr

namespace dunedaq {

  ERS_DECLARE_ISSUE(fdreadoutlibs,
		    InvalidDataSize,
		    " Unable to load all data. Size of data (" << data_size 
		    << ") is greater than PACMAN frame size (" << crt_frame_size << ").",
		    ((uint64_t)data_size)((uint64_t)crt_frame_size)) // NOLINT

  namespace fdreadoutlibs {
    namespace types {

      /**
       * @brief CRT frame
       * Size = 288[Bytes] (16+16+4*64)
       * */
      const constexpr std::size_t CRT_FRAME_SIZE = 288;
      struct CRTTypeAdapter
      {
	using FrameType = CRTTypeAdapter;
	// data
	std::vector<char> data{std::vector<char>(CRT_FRAME_SIZE,0)};
	void load_message( const void * load_data, const unsigned int size ) {
	  if( size > data.size() ) {
	    ers::error(InvalidDataSize(ERS_HERE, size, data.size()));
	    return;
	  }
	  memcpy(&data[0], load_data, size);
	}

	// comparable based on first timestamp
	bool operator<(const CRTTypeAdapter& other) const
	{
	  auto thisptr = reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&data[0]);        // NOLINT
	  auto otherptr = reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&other.data[0]); // NOLINT
	  return (thisptr->get_timestamp()) <  (otherptr->get_timestamp()) ? true : false;
	}

	uint64_t get_first_timestamp() const 
        { 
          return reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&data[0])->get_timestamp(); // NOLINT
        }

        void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
        {
          auto frame = reinterpret_cast<dunedaq::fddetdataformats::CRTFrame*>(&data[0]); // NOLINT
          frame->set_timestamp(ts);
        }

  void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset = 2048*/ ) // NOLINT(build/unsigned)
  {
    auto wef = reinterpret_cast<dunedaq::fddetdataformats::CRTFrame*>(((uint8_t*)(&data[0]))); // NOLINT
    wef->set_timestamp(first_timestamp);
  }

  void fake_geoid(uint16_t crate_id, uint16_t slot_id, uint16_t stream_id) {
      for (unsigned int i = 0; i < get_num_frames(); ++i) {
        auto df = reinterpret_cast<dunedaq::fddetdataformats::CRTFrame*>((reinterpret_cast<uint8_t*>(&data[0])) + i * get_frame_size());
        df->daq_header->crate_id = crate_id;
        df->daq_header->slot_id = slot_id;
        df->daq_header->stream_id = stream_id;
      }   
  }

  void fake_adc_pattern(int channel) 
  {
    auto frame = reinterpret_cast<dunedaq::fddetdataformats::CRTFrame*>(&data[0]); // NOLINT
    // Set the ADC to the int16 maximum value 
    frame->set_adc(channel, 8191); 
  }
    
  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT
  {
    // Set error bits in header
  }
    
	void inspect_message() const
	{
	  TLOG_DEBUG(1) << "Frame timestamp: " << get_first_timestamp();

	  uint16_t moduleNum = //NOLINT
	     reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&data[0])
	     ->crt_header->module_number;

	  TLOG_DEBUG(1) << "Module number: " << moduleNum;

	  uint8_t numHits = // NOLINT
	    reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&data[0])
	    ->crt_header->hit_count; // NOLINT

	  TLOG_DEBUG(1) << "Num hits in frame: " << numHits;

	  for (unsigned int i = 0; i < numHits; i++) {
	    TLOG_DEBUG(1) << "Inspecting hit " << i;

	    dunedaq::fddetdataformats::CRTFrame::CRTHit* theHit =
	      reinterpret_cast<const dunedaq::fddetdataformats::CRTFrame*>(&data[0])->get_crt_hit((void*)&data[0], i); // NOLINT

	    TLOG_DEBUG(1) << "Hit magic number: " << (char)theHit->hit_magic_number; // NOLINT
	    TLOG_DEBUG(1) << "Hit channel: " << (uint8_t)theHit->hit_channel; // NOLINT
	    TLOG_DEBUG(1) << "Hit ADC: " << (int16_t)theHit->hit_adc; // NOLINT
	  }
	}

	FrameType* begin()
	{
	  return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
	}

	FrameType* end()
	{
	  return reinterpret_cast<FrameType*>(data[0] + data.size()); // NOLINT
	}

	size_t get_payload_size() { return data.size(); }

	size_t get_num_frames() { return 1; }

	size_t get_frame_size() { return data.size(); }

	static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
	static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kCRT;

	// Set the right value for this field
	static const constexpr uint64_t expected_tick_difference = 0; // NOLINT(build/unsigned)
      };
    } // namespace types
  } // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRTTYPEADAPTER_HPP_
