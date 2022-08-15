/**
 * @file FDReadoutTypes.hpp Common types in udaq-readout
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTTYPES_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTTYPES_HPP_

//#include "iomanager/Sender.hpp"
//#include "iomanager/Receiver.hpp"

#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "detdataformats/daphne/DAPHNEFrame.hpp"
#include "detdataformats/ssp/SSPTypes.hpp"
#include "detdataformats/wib/WIBFrame.hpp"
#include "detdataformats/wib2/WIB2Frame.hpp"
#include "detdataformats/tde/TDE16Frame.hpp"
#include "detdataformats/fwtp/RawTp.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

// Location of this struct is very bad very bad
struct TriggerPrimitive
{
  TriggerPrimitive(uint64_t messageTimestamp_,  // NOLINT(build/unsigned)
                   uint16_t channel_,           // NOLINT(build/unsigned)
                   uint16_t endTime_,           // NOLINT(build/unsigned)
                   uint16_t charge_,            // NOLINT(build/unsigned)
                   uint16_t timeOverThreshold_) // NOLINT(build/unsigned)
    : messageTimestamp(messageTimestamp_)
    , channel(channel_)
    , endTime(endTime_)
    , charge(charge_)
    , timeOverThreshold(timeOverThreshold_)
  {}

  // The timestamp of the netio message that this hit comes from
  uint64_t messageTimestamp; // NOLINT(build/unsigned)
  // The electronics channel number within the (crate, slot, fiber)
  uint16_t channel; // NOLINT(build/unsigned)
  // In TPC ticks relative to the start of the netio message
  uint16_t endTime; // NOLINT(build/unsigned)
  // In ADC
  uint16_t charge; // NOLINT(build/unsigned)
  // In *TPC* clock ticks
  uint16_t timeOverThreshold; // NOLINT(build/unsigned)
};

/**
 * @brief SuperChunk concept: The FELIX user payloads are called CHUNKs.
 * There is mechanism in firmware to aggregate WIB frames to a user payload
 * that is called a SuperChunk. Default mode is with 12 frames:
 * 12[WIB frames] x 464[Bytes] = 5568[Bytes]
 */
const constexpr std::size_t WIB_SUPERCHUNK_SIZE = 5568; // for 12: 5568
struct WIB_SUPERCHUNK_STRUCT
{
  using FrameType = dunedaq::detdataformats::wib::WIBFrame;

  // data
  char data[WIB_SUPERCHUNK_SIZE];
  // comparable based on first timestamp
  bool operator<(const WIB_SUPERCHUNK_STRUCT& other) const
  {
    // auto thisptr = reinterpret_cast<const dunedaq::detdataformats::WIBHeader*>(&data);        // NOLINT
    // auto otherptr = reinterpret_cast<const dunedaq::detdataformats::WIBHeader*>(&other.data); // NOLINT
    return this->get_first_timestamp() < other.get_first_timestamp();
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const dunedaq::detdataformats::wib::WIBFrame*>(&data)->get_wib_header()->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(&data)->get_wib_header()->set_timestamp(ts); // NOLINT
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 25) // NOLINT(build/unsigned)
  {
    uint64_t ts_next = first_timestamp;                                               // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(((uint8_t*)(&data))); // NOLINT
    for (unsigned int i = 0; i < 12; ++i) {
      auto wfh = const_cast<dunedaq::detdataformats::wib::WIBHeader*>(wf->get_wib_header());
      wfh->set_timestamp(ts_next);
      ts_next += offset;
      wf++;
    }
  }

  void fake_frame_errors(std::vector<uint16_t>* fake_errors) // NOLINT(build/unsigned)
  {
    auto wf = reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 12; ++i) {
      wf->set_wib_errors((*fake_errors)[i]);
      wf++;
    }
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + WIB_SUPERCHUNK_SIZE); // NOLINT
  }

  size_t get_payload_size() { return 5568; }

  size_t get_num_frames() { return 12; }

  size_t get_frame_size() { return 464; }
  
  static const constexpr size_t fixed_payload_size = 5568;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kProtoWIB;
  static const constexpr uint64_t expected_tick_difference = 25; // 2 MHz@50MHz clock // NOLINT(build/unsigned)
};

const constexpr std::size_t TDE_AMC_STRUCT_SIZE = 64 * sizeof(dunedaq::detdataformats::tde::TDE16Frame);
struct TDE_AMC_STRUCT
{
  using FrameType = dunedaq::detdataformats::tde::TDE16Frame;

  char data[TDE_AMC_STRUCT_SIZE];

  bool operator<(const TDE_AMC_STRUCT& other) const
  {
    auto thisptr = reinterpret_cast<const FrameType*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const FrameType*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const FrameType*>(&data)->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    reinterpret_cast<FrameType*>(&data)->get_tde_header()->set_timestamp(ts); // NOLINT
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t) // NOLINT(build/unsigned)
  {
    auto tdef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 64; ++i) {
      auto wfh = const_cast<dunedaq::detdataformats::tde::TDE16Header*>(tdef->get_tde_header());
      wfh->set_timestamp(first_timestamp);
      tdef++;
    }
  }

  void fake_frame_errors(std::vector<uint16_t>* fake_errors) // NOLINT(build/unsigned)
  {
    auto tdef = reinterpret_cast<FrameType*>(((uint8_t*)(&data))); // NOLINT
    for (int i = 0; i < 64; ++i) {
      tdef->set_tde_errors((*fake_errors)[i]);
      tdef++;
    }
  }

  size_t get_payload_size() { return 64 * sizeof(FrameType); }

  size_t get_num_frames() { return 64; }

  size_t get_frame_size() { return sizeof(FrameType); }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + TDE_AMC_STRUCT_SIZE); // NOLINT
  }

  // static const constexpr size_t fixed_payload_size = 5568;
  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kTDE_AMC;
  static const constexpr uint64_t expected_tick_difference = 1000; // NOLINT(build/unsigned)

};

static_assert(sizeof(struct dunedaq::detdataformats::wib::WIBFrame)*12 == WIB_SUPERCHUNK_SIZE,
              "Check your assumptions on WIB_SUPERCHUNK_STRUCT");

/**
 * @brief For WIB2 the numbers are different.
 * 12[WIB2 frames] x 468[Bytes] = 5616[Bytes]
 * */
const constexpr std::size_t WIB2_SUPERCHUNK_SIZE = 5664; // for 12: 5664
struct WIB2_SUPERCHUNK_STRUCT
{
  using FrameType = dunedaq::detdataformats::wib2::WIB2Frame;
  // data
  char data[WIB2_SUPERCHUNK_SIZE];
  // comparable based on first timestamp
  bool operator<(const WIB2_SUPERCHUNK_STRUCT& other) const
  {
    auto thisptr = reinterpret_cast<const dunedaq::detdataformats::wib2::WIB2Frame*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const dunedaq::detdataformats::wib2::WIB2Frame*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const dunedaq::detdataformats::wib2::WIB2Frame*>(&data)->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    auto frame = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(&data); // NOLINT
    frame->header.timestamp_1 = ts;
    frame->header.timestamp_2 = ts >> 32;
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 32) // NOLINT(build/unsigned)
  {
    uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)
    for (unsigned int i = 0; i < 12; ++i) {
      auto w2f = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(((uint8_t*)(&data)) + i * 472); // NOLINT
      w2f->header.timestamp_1 = ts_next;
      w2f->header.timestamp_2 = ts_next >> 32;
      ts_next += offset;
    }
  }

  void fake_frame_errors(std::vector<uint16_t>* /*fake_errors*/) // NOLINT
  {
    // Set error bits in header
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(&data[0]); // NOLINT
  }

  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(data + WIB2_SUPERCHUNK_SIZE); // NOLINT
  }

  size_t get_payload_size() { return 5664; }

  size_t get_num_frames() { return 12; }

  size_t get_frame_size() { return 472; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kWIB;
  static const constexpr uint64_t expected_tick_difference = 32; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct dunedaq::detdataformats::wib2::WIB2Frame)*12 == WIB2_SUPERCHUNK_SIZE,
              "Check your assumptions on WIB2_SUPERCHUNK_STRUCT");

/**
 * @brief For DAPHNE the numbers are different.
 * 12[DAPHNE frames] x 584[Bytes] = 7008[Bytes]
 * */
const constexpr std::size_t DAPHNE_SUPERCHUNK_SIZE = 7008; // for 12: 7008
struct DAPHNE_SUPERCHUNK_STRUCT
{
  using FrameType = dunedaq::detdataformats::daphne::DAPHNEFrame;
  // data
  char data[DAPHNE_SUPERCHUNK_SIZE];
  // comparable based on first timestamp
  bool operator<(const DAPHNE_SUPERCHUNK_STRUCT& other) const
  {
    auto thisptr = reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data);        // NOLINT
    auto otherptr = reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&other.data); // NOLINT
    return thisptr->get_timestamp() < otherptr->get_timestamp() ? true : false;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return reinterpret_cast<const dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data)->get_timestamp(); // NOLINT
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    auto frame = reinterpret_cast<dunedaq::detdataformats::daphne::DAPHNEFrame*>(&data); // NOLINT
    frame->header.timestamp_wf_1 = ts;
    frame->header.timestamp_wf_2 = ts >> 32;
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t offset = 25) // NOLINT(build/unsigned)
  {
    uint64_t ts_next = first_timestamp; // NOLINT(build/unsigned)
    for (unsigned int i = 0; i < 12; ++i) {
      auto df = reinterpret_cast<dunedaq::detdataformats::daphne::DAPHNEFrame*>(((uint8_t*)(&data)) + i * 584); // NOLINT
      df->header.timestamp_wf_1 = ts_next;
      df->header.timestamp_wf_2 = ts_next >> 32;
      ts_next += offset;
    }
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
    return reinterpret_cast<FrameType*>(data + DAPHNE_SUPERCHUNK_SIZE); // NOLINT
  }

  size_t get_payload_size() { return 7008; }

  size_t get_num_frames() { return 12; }

  size_t get_frame_size() { return 584; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kDAPHNE;
  static const constexpr uint64_t expected_tick_difference = 16; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct DAPHNE_SUPERCHUNK_STRUCT) == DAPHNE_SUPERCHUNK_SIZE,
              "Check your assumptions on DAPHNE_SUPERCHUNK_STRUCT");

const constexpr std::size_t TP_SIZE = sizeof(triggeralgs::TriggerPrimitive);
struct SW_WIB_TRIGGERPRIMITIVE_STRUCT
{
  using FrameType = SW_WIB_TRIGGERPRIMITIVE_STRUCT;
  // data
  triggeralgs::TriggerPrimitive tp;
  // comparable based on start timestamp
  bool operator<(const SW_WIB_TRIGGERPRIMITIVE_STRUCT& other) const
  {
    return this->tp.time_start < other.tp.time_start;
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return tp.time_start;
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    tp.time_start = ts;
  }

  uint64_t get_timestamp() const // NOLINT(build/unsigned)
  {
    return tp.time_start;
  }

  void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset = 25*/) // NOLINT(build/unsigned)
  {
    tp.time_start = first_timestamp;
  }

  FrameType* begin() { return this; }

  FrameType* end() { return (this + 1); } // NOLINT

  size_t get_payload_size() { return TP_SIZE; }

  size_t get_num_frames() { return 1; }

  size_t get_frame_size() { return TP_SIZE; }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kTrigger;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kSW_TriggerPrimitive;
  static const constexpr uint64_t expected_tick_difference = 25; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct SW_WIB_TRIGGERPRIMITIVE_STRUCT) == sizeof(triggeralgs::TriggerPrimitive),
              "Check your assumptions on SW_WIB_TRIGGERPRIMITIVE_STRUCT");

const constexpr std::size_t SSP_FRAME_SIZE = 1012;
struct SSP_FRAME_STRUCT
{
  using FrameType = SSP_FRAME_STRUCT;

  // header
  detdataformats::ssp::EventHeader header;

  // data
  char data[SSP_FRAME_SIZE];

  // comparable based on start timestamp
  bool operator<(const SSP_FRAME_STRUCT& other) const
  {
    return this->get_first_timestamp() < other.get_first_timestamp() ? true : false;
  }

  uint64_t get_timestamp() const  // NOLINT(build/unsigned)
  {
    return get_first_timestamp();
  }

  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    auto ehptr = &header;
    unsigned long ts = 0; // NOLINT(runtime/int)
    for (unsigned int iword = 0; iword <= 3; ++iword) {
      ts += ((unsigned long)(ehptr->timestamp[iword])) << 16 * iword; //NOLINT(runtime/int)
    }
    return ts;
  }

  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    uint64_t bitmask = (1 << 16) - 1; // NOLINT(build/unsigned)
    for (unsigned int iword = 0; iword <= 3; ++iword) {
      header.timestamp[iword] = static_cast<uint16_t>((ts & bitmask)); // NOLINT(build/unsigned)
      ts = ts >> 16;
    }
  }

  void fake_timestamps(uint64_t /*first_timestamp*/, uint64_t /*offset = 25*/) // NOLINT(build/unsigned)
  {
    // tp.time_start = first_timestamp;
  }

  FrameType* begin() { return this; }

  FrameType* end() { return (this + 1); } // NOLINT

  size_t get_payload_size() {
    return SSP_FRAME_SIZE;
  }

  size_t get_num_frames() {
    return 1;
  }

  size_t get_frame_size() {
    return SSP_FRAME_SIZE;
  }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kDetectorReadout;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kUnknown;
  static const constexpr uint64_t expected_tick_difference = 25; // NOLINT(build/unsigned)
};

static_assert(sizeof(struct SSP_FRAME_STRUCT) == sizeof(detdataformats::ssp::EventHeader) + SSP_FRAME_SIZE,
              "Check your assumptions on SSP_FRAME_STRUCT");

/**
 * @brief Convencience wrapper to take ownership over char pointers with
 * corresponding allocated memory size.
 * */
struct VariableSizePayloadWrapper
{
  VariableSizePayloadWrapper() {}
  VariableSizePayloadWrapper(size_t size, char* data)
    : size(size)
    , data(data)
  {}

  size_t size = 0;
  std::unique_ptr<char> data = nullptr;
};

// raw WIB TP
struct RAW_WIB_TRIGGERPRIMITIVE_STRUCT
{
  RAW_WIB_TRIGGERPRIMITIVE_STRUCT()
  {
    m_raw_tp_frame_chunksize = 0;
    m_first_timestamp = 0;
  }

  using FrameType = RAW_WIB_TRIGGERPRIMITIVE_STRUCT;

  // latency buffer 
  bool operator<(const RAW_WIB_TRIGGERPRIMITIVE_STRUCT& other) const 
  {
    return this->get_first_timestamp() < other.get_first_timestamp();
  }
  
  uint64_t get_timestamp() const  // NOLINT(build/unsigned)
  {
    return get_first_timestamp();
  }
  uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
  {
    return m_first_timestamp;
  }
  void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
  {
    m_first_timestamp = ts;
  }
  void fake_timestamps(uint64_t first_timestamp, uint64_t /*offset = 25*/) // NOLINT(build/unsigned)
  {
    m_first_timestamp = first_timestamp;
  }

  FrameType* begin()
  {
    return reinterpret_cast<FrameType*>(m_raw_tp_frame_chunk.data()); // NOLINT
  }
  FrameType* end()
  {
    return reinterpret_cast<FrameType*>(m_raw_tp_frame_chunk.data()+m_raw_tp_frame_chunksize); // NOLINT
  }

  static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kTrigger;
  static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kSW_TriggerPrimitive;
  static const constexpr uint64_t expected_tick_difference = 0; // 2 MHz@50MHz clock // NOLINT(build/unsigned)
  //static const constexpr size_t frame_size = TP_SIZE;
  //static const constexpr size_t element_size = TP_SIZE;
  static const constexpr uint64_t tick_dist = 25; // 2 MHz@50MHz clock // NOLINT(build/unsigned)
  static const constexpr uint8_t frames_per_element = 1; // NOLINT(build/unsigned)

  // raw WIB TP frames are variable size
  size_t get_payload_size() {
    return m_raw_tp_frame_chunksize;
  }

  size_t get_num_frames() {
    return 1;
  }

  size_t get_frame_size() {
    return m_raw_tp_frame_chunksize;
  }

  void set_raw_tp_frame_chunk(std::vector<char>& source)
  {
    int bsize = source.capacity();
    m_raw_tp_frame_chunk.reserve(bsize);
    ::memcpy(static_cast<void*>(m_raw_tp_frame_chunk.data()),
             static_cast<void*>(source.data()),
             bsize);
    m_raw_tp_frame_chunksize = bsize;

    unpack_timestamp();
  }

  std::vector<std::uint8_t>& get_data()  // NOLINT(build/unsigned)
  {
    return std::ref(m_raw_tp_frame_chunk);
  }

  void set_data_size(const int& bytes)
  {
    m_raw_tp_frame_chunksize = bytes;

    unpack_timestamp();
  }
  int get_raw_tp_frame_chunksize()
  {
    return m_raw_tp_frame_chunksize;
  }

private:
void unpack_timestamp() {
    std::unique_ptr<detdataformats::fwtp::RawTp> rwtp = 
                    std::make_unique<detdataformats::fwtp::RawTp>();
    ::memcpy(static_cast<void*>(&rwtp->m_head),
             static_cast<void*>(m_raw_tp_frame_chunk.data()),
             2*RAW_WIB_TP_SUBFRAME_SIZE); 
    m_first_timestamp = rwtp->m_head.get_timestamp();
} 

private:
  std::vector<std::uint8_t> m_raw_tp_frame_chunk;  // NOLINT(build/unsigned)
  int m_raw_tp_frame_chunksize;

  uint64_t m_first_timestamp;
  static const constexpr std::size_t RAW_WIB_TP_SUBFRAME_SIZE = 12;

};

struct TpSubframe
{
  uint32_t word1; // NOLINT(build/unsigned)
  uint32_t word2; // NOLINT(build/unsigned)
  uint32_t word3; // NOLINT(build/unsigned)
};

} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTTYPES_HPP_
