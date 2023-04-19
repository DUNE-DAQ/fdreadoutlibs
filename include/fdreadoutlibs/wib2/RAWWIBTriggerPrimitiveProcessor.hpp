/**
 * @file RAWWIBTriggerPrimitiveProcessor.hpp WIB TP specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_RAWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_RAWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_

#include "appfwk/DAQModuleHelper.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"


#include "detdataformats/fwtp/RawTp.hpp"
#include "logging/Logging.hpp"
#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutLogging.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"
#include "fdreadoutlibs/wib2/WIB2TPHandler.hpp"
#include "fdreadoutlibs/DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "rcif/cmd/Nljs.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"
#include "trigger/TPSet.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

using dunedaq::readoutlibs::logging::TLVL_WORK_STEPS;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

namespace dunedaq {
namespace fdreadoutlibs {

class RAWWIBTriggerPrimitiveProcessor
  : public readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>;
  using frame_ptr = types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter*;
  using rwtp_ptr = detdataformats::fwtp::RawTp*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Channel map function type
  typedef int (*chan_map_fn_t)(int);
  
  explicit RAWWIBTriggerPrimitiveProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>(error_registry)
    , m_fw_tpg_enabled(false)
  {}

  void conf(const nlohmann::json& args) override
  {
    auto config = args["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
    m_sourceid.id = config.source_id;
    m_sourceid.subsystem = types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter::subsystem;

    TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>::add_preprocess_task(
                std::bind(&RAWWIBTriggerPrimitiveProcessor::tp_unpack, this, std::placeholders::_1));
    TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>::conf(args);

    if (config.enable_firmware_tpg) {
      m_fw_tpg_enabled = true;
      daqdataformats::SourceID tpset_sourceid;
      tpset_sourceid.id = config.tpset_sourceid;
      tpset_sourceid.subsystem = daqdataformats::SourceID::Subsystem::kTrigger;
      m_tphandler.reset(
            new WIB2TPHandler(m_tp_sink, m_tpset_sink, config.tp_timeout, config.tpset_window_size, tpset_sourceid));
    }

    m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

    m_stitch_constant = config.fwtp_number_of_ticks * config.fwtp_tick_length;
    m_time_tick = config.fwtp_tick_length;
    m_enable_fake_timestamp = config.fwtp_fake_timestamp;
    TLOG_DEBUG(15) << "TP frame stitching parameters are ( " << m_stitch_constant << ", " << m_time_tick << ")";
  }

  void init(const nlohmann::json& args) override
  {
    m_fake_timestamp = 0;
    try {
      auto queue_index = appfwk::connection_index(args, {});
      if (queue_index.find("tp_out") != queue_index.end()) {
        m_tp_sink = get_iom_sender<types::TriggerPrimitiveTypeAdapter>(queue_index["tp_out"]);
      }
      if (queue_index.find("tpset_out") != queue_index.end()) {
        m_tpset_sink = get_iom_sender<trigger::TPSet>(queue_index["tpset_out"]);
      }
    } catch (const ers::Issue& excpt) {
      throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
    }

/*
    try {
      auto queue_index = appfwk::connection_index(args, {});
      if (queue_index.find("tp_out") != queue_index.end()) {
        m_tp_sink = get_iom_sender<types::TriggerPrimitiveTypeAdapter>(queue_index["tp_out"]);
      }
      if (queue_index.find("tpset_out") != queue_index.end()) {
        m_tpset_sink = get_iom_sender<trigger::TPSet>(queue_index["tpset_out"]);
      }
      m_err_frame_sink = get_iom_sender<detdataformats::wib::WIBFrame>(queue_index["errored_frames"]);
    } catch (const ers::Issue& excpt) {
      throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
    }
i*/

  }

  void start(const nlohmann::json& args) override
  {
    if (m_fw_tpg_enabled) {
      rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
      m_tphandler->set_run_number(start_params.run);
      m_tphandler->reset();
    }
    // Reset stats 
    m_tps_stitched = 0;
    m_tp_frames = 0;
    m_tp_hits = 0;
    m_tps_dropped = 0;
    m_sent_tps = 0;
    std::fill(m_nhits.begin(), m_nhits.end(), 0);
    m_total_hits_count.exchange(0);
  }

  void stop(const nlohmann::json& /*args*/) override
  {
    TLOG_DEBUG(20) << "Number of TP frames " << m_tp_frames;
    TLOG_DEBUG(20) << "Number of TPs stitched " << m_tps_stitched;
    TLOG_DEBUG(20) << "Number of TPs dropped " << m_tps_dropped; 

    for (size_t i = 0; i < m_nhits.size(); i++) {
      TLOG_DEBUG(20) << "Number of frames with hits " << i << ": " << m_nhits[i] << ", " 
        << static_cast<double>(static_cast<double>(m_nhits[i])/static_cast<double>(m_tp_frames)) << "\n";
    }
  }

  void scrap(const nlohmann::json& args) override
  {
    m_tphandler.reset();

    TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>::scrap(args);
  }

  void get_info(opmonlib::InfoCollector& ci, int level)
  {
    readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>::get_info(ci, level);

    readoutlibs::readoutinfo::RawDataProcessorInfo info;

    if (m_tphandler != nullptr) {
      info.num_tps_sent = m_tphandler->get_and_reset_num_sent_tps();
      info.num_tpsets_sent = m_tphandler->get_and_reset_num_sent_tpsets();
      info.num_tps_dropped = m_tps_dropped.exchange(0);
    }
    auto now = std::chrono::high_resolution_clock::now();
    if (m_fw_tpg_enabled) {
      int new_hits = m_total_hits_count.exchange(0);
      double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
      TLOG_DEBUG(TLVL_TAKE_NOTE) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
      TLOG_DEBUG(TLVL_TAKE_NOTE) << "Total new hits: " << new_hits;
      info.rate_tp_hits = new_hits / seconds / 1000.;
    }
    m_t0 = now;

    ci.add(info);
  }

void tp_stitch(rwtp_ptr rwtp)
{
  m_tp_frames++;
  uint64_t ts_0;
  if(m_enable_fake_timestamp == true)
  {
    m_fake_timestamp += 6400; 
    ts_0 = m_fake_timestamp; // NOLINT
  }
  else
  {
    ts_0 = rwtp->m_head.get_timestamp(); // NOLINT
  }
  int nhits = rwtp->m_head.get_nhits(); // NOLINT
  uint8_t m_channel_no = rwtp->m_head.m_wire_no; // NOLINT
  uint8_t m_fiber_no = rwtp->m_head.m_fiber_no; // NOLINT
  uint8_t m_crate_no = rwtp->m_head.m_crate_no; // NOLINT
  uint8_t m_slot_no = (rwtp->m_head.m_slot_no) & ((uint8_t) 0x7); // NOLINT
  uint offline_channel = m_channel_map->get_offline_channel_from_crate_slot_fiber_chan(m_crate_no, m_slot_no, m_fiber_no, m_channel_no);
 
  if (nhits < 8) { 
    m_nhits[nhits] += 1;
  }
  m_total_hits_count += nhits;

  m_tp_hits += nhits;

  for (int i = 0; i < nhits; i++) {

    triggeralgs::TriggerPrimitive trigprim;
    trigprim.time_start = ts_0 + rwtp->m_blocks[i].m_start_time * m_time_tick;
    trigprim.time_peak = ts_0 + rwtp->m_blocks[i].m_peak_time * m_time_tick;
    trigprim.time_over_threshold = (rwtp->m_blocks[i].m_end_time - rwtp->m_blocks[i].m_start_time) * m_time_tick;
    trigprim.channel = offline_channel; //offline_channel; // m_channel_no;
    trigprim.adc_integral = rwtp->m_blocks[i].m_sum_adc;
    trigprim.adc_peak = rwtp->m_blocks[i].m_peak_adc;
    trigprim.detid =
            m_fiber_no; // TODO: convert crate/slot/fiber to SourceID Roland Sipos rsipos@cern.ch July-22-2021
    trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
    trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
    trigprim.version = 1;

    // stitch current hit to previous hit
    if (m_A[m_channel_no][m_fiber_no].size() == 1) {
      if (static_cast<int>(rwtp->m_blocks[i].m_start_time) == 0
          && (
          static_cast<int>(trigprim.time_start) - static_cast<int>(m_A[m_channel_no][m_fiber_no][0].time_start)
             <= static_cast<int>(m_stitch_constant)
          || (static_cast<int>(trigprim.time_start) - static_cast<int>(m_T[m_channel_no][m_fiber_no][0])
             <= static_cast<int>(m_stitch_constant))
             )
          ) {
        // current hit is continuation of previous hit
        m_T[m_channel_no][m_fiber_no].clear();
        if (trigprim.adc_peak > m_A[m_channel_no][m_fiber_no][0].adc_peak) {
          m_A[m_channel_no][m_fiber_no][0].time_peak = trigprim.time_peak;
          m_A[m_channel_no][m_fiber_no][0].adc_peak = trigprim.adc_peak;
        }
        m_A[m_channel_no][m_fiber_no][0].time_over_threshold += trigprim.time_over_threshold;
        m_A[m_channel_no][m_fiber_no][0].adc_integral += trigprim.adc_integral;
        m_T[m_channel_no][m_fiber_no].push_back(trigprim.time_start);

      } else {
        // current hit is not continuation of previous hit
        // add previous hit to TriggerPrimitives
        if (m_tphandler != nullptr && !m_tphandler->add_tp(std::move(m_A[m_channel_no][m_fiber_no][0]), ts_0)) {
          m_tps_dropped++;
        }
        m_tps_stitched++;
        if (m_tphandler != nullptr) { m_tphandler->try_sending_tpsets(ts_0); } 
        m_A[m_channel_no][m_fiber_no].clear();
        m_T[m_channel_no][m_fiber_no].clear();
      }
    }

    // NB for TPSets: this assumes hits come ordered in time 
    // current hit (is, completes or starts) one TriggerPrimitive 
    uint8_t m_tp_continue = rwtp->m_blocks[i].m_hit_continue; // NOLINT
    uint8_t m_tp_end_time = rwtp->m_blocks[i].m_end_time; // NOLINT
 
    if (m_tp_continue == 0 && m_tp_end_time != 63) {
      if (m_A[m_channel_no][m_fiber_no].size() == 1) {
        // the current hit completes one stitched TriggerPrimitive
        if (m_tphandler != nullptr && !m_tphandler->add_tp(std::move(m_A[m_channel_no][m_fiber_no][0]), ts_0)) {
          m_tps_dropped++;
        }
        if (m_tphandler != nullptr) { m_tphandler->try_sending_tpsets(ts_0); } 
        m_tps_stitched++;
        m_A[m_channel_no][m_fiber_no].clear();
        m_T[m_channel_no][m_fiber_no].clear();
      } else {
        // the current hit is one TriggerTrimitive
        if (m_tphandler != nullptr && !m_tphandler->add_tp(std::move(trigprim), ts_0)) {
          m_tps_dropped++;
        }
        if (m_tphandler != nullptr) { m_tphandler->try_sending_tpsets(ts_0); } 
        m_tps_stitched++;      
      }
    } else {
      // the current hit starts one TriggerPrimitive
      if (m_A[m_channel_no][m_fiber_no].size() == 0) {
        m_A[m_channel_no][m_fiber_no].push_back(trigprim);
        m_T[m_channel_no][m_fiber_no].push_back(trigprim.time_start);
      } else { // decide to add long TriggerPrimitive even when it doesn't end properly
               // this is rare case and can be removed for efficiency    
        // the current hit is "bad"
        // add one TriggerPrimitive from previous stitched hits except the current hit  
        if ( m_tp_continue == 0 && m_tp_end_time == 63 &&
             static_cast<int>(trigprim.time_start) - static_cast<int>(m_T[m_channel_no][m_fiber_no][0])
             <= static_cast<int>(m_stitch_constant)) {
          if (m_tphandler != nullptr && !m_tphandler->add_tp(std::move(m_A[m_channel_no][m_fiber_no][0]), ts_0)) {
            m_tps_dropped++;
          }
          if (m_tphandler != nullptr) { m_tphandler->try_sending_tpsets(ts_0); }
          m_tps_stitched++;      
          m_A[m_channel_no][m_fiber_no].clear();
          m_T[m_channel_no][m_fiber_no].clear();
        }
      }
    }
  }
} // NOLINT (exceeding 80 lines)


void tp_unpack(frame_ptr fr)  
{
  auto& srcbuffer = fr->get_data();
  int num_elem = fr->get_raw_tp_frame_chunksize();

  if (num_elem == 0) {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "No raw WIB TP elements to read from buffer! ";
    return;
  }
  if (num_elem % RAW_WIB_TP_SUBFRAME_SIZE != 0) {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "Raw WIB TP elements not multiple of subframe size (3)! ";
    return;
  }

  double dbg_frames = 0;
  double dbg_microseconds = 0;
  int offset = 0;
  while (offset <= num_elem) {

    if (offset == num_elem) break;
 
    auto now = std::chrono::high_resolution_clock::now();
    // Count number of subframes in a TP frame
    int n;
    bool ped_found { false };
    for (n=2; offset + n*RAW_WIB_TP_SUBFRAME_SIZE<(size_t)num_elem; ++n) {
      if (reinterpret_cast<types::TpSubframe*>(((uint8_t*)srcbuffer.data()) // NOLINT
           + offset + n*RAW_WIB_TP_SUBFRAME_SIZE)->word3 == 0xDEADBEEF) {
        ped_found = true;
        break; 
      }  
    }
    // Found no pedestal block
    if (!ped_found) {
      TLOG() << "Debug message: Raw WIB TP chunk contains no TP frames! Chunk size / offset / subframes is " << num_elem << " / " << offset << " / " << n;
      return;
    }

    int bsize = n * RAW_WIB_TP_SUBFRAME_SIZE;
    std::vector<char> tmpbuffer;
    tmpbuffer.reserve(bsize);
    int nhits = n - 1;  // n is subframe counter (starting from 0, not 1)
    // add header block 
    ::memcpy(static_cast<void*>(tmpbuffer.data() + 0),
             static_cast<void*>(srcbuffer.data() + offset),
             RAW_WIB_TP_SUBFRAME_SIZE);

    // add pedinfo block 
    ::memcpy(static_cast<void*>(tmpbuffer.data() + RAW_WIB_TP_SUBFRAME_SIZE),
             static_cast<void*>(srcbuffer.data() + offset + (n-1)*RAW_WIB_TP_SUBFRAME_SIZE),
             RAW_WIB_TP_SUBFRAME_SIZE);

    // add TP hits
    ::memcpy(static_cast<void*>(tmpbuffer.data() + 2*RAW_WIB_TP_SUBFRAME_SIZE),
             static_cast<void*>(srcbuffer.data() + offset + RAW_WIB_TP_SUBFRAME_SIZE),
             nhits*RAW_WIB_TP_SUBFRAME_SIZE);

    auto heap_memory_block = malloc(
         sizeof(dunedaq::detdataformats::fwtp::TpHeader) + 
         nhits * sizeof(dunedaq::detdataformats::fwtp::TpData));
    rwtp_ptr rwtp =
         static_cast<dunedaq::detdataformats::fwtp::RawTp*>(heap_memory_block);

    ::memcpy(static_cast<void*>(&rwtp->m_head),
             static_cast<void*>(tmpbuffer.data() + 0),
             2*RAW_WIB_TP_SUBFRAME_SIZE);

    for (int i=0; i<nhits; i++) {
      ::memcpy(static_cast<void*>(&rwtp->m_blocks[i]),
               static_cast<void*>(tmpbuffer.data() + (2+i)*RAW_WIB_TP_SUBFRAME_SIZE),
               RAW_WIB_TP_SUBFRAME_SIZE);
    }
    dbg_frames++;
    dbg_microseconds += std::chrono::duration_cast<std::chrono::microseconds>(now - m_u_t0).count();
    m_u_t0 = now;   
 
    // old format lacks number of hits
    rwtp->set_nhits(nhits); // explicitly set number of hits in new format

    // stitch TP hits
    tp_stitch(rwtp);
    offset += (2+nhits)*RAW_WIB_TP_SUBFRAME_SIZE;
    free(heap_memory_block);
  }
}

protected:
  double m_time_tick { 1.0 }; // 

private:
  using source_t = iomanager::ReceiverConcept<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>;
  std::shared_ptr<source_t> m_tp_source;

  // unpacking
  static const constexpr std::size_t RAW_WIB_TP_SUBFRAME_SIZE = 12;

  // stitching algorithm
  std::vector<triggeralgs::TriggerPrimitive> m_A[256][10]; // keep track of TPs to stitch per channel
  std::vector<uint64_t> m_T[256][10]; // NOLINT // keep track of last stitched start time
  std::atomic<uint64_t> m_tps_stitched { 0 }; // NOLINT
  std::atomic<uint64_t> m_tp_frames  { 0 }; // NOLINT
  std::atomic<int> m_tp_hits { 0 };

  int m_stitch_constant { 2048 }; // number of ticks between WIB-to-TP packets
  std::vector<int> m_nhits { 0, 0, 0, 0, 0, 0, 0, 0};

  // interface to DS
  bool m_fw_tpg_enabled;
  bool m_enable_fake_timestamp;
  std::shared_ptr<iomanager::SenderConcept<types::TriggerPrimitiveTypeAdapter>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> m_tpset_sink;
  std::unique_ptr<WIB2TPHandler> m_tphandler;
  std::atomic<uint64_t> m_tps_dropped{ 0 }; // NOLINT
  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;
  uint64_t m_fake_timestamp { 0 }; // NOLINT

  // info
  std::atomic<uint64_t> m_sent_tps{ 0 }; // NOLINT(build/unsigned)
  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_u_t0;
  std::atomic<int> m_total_hits_count{ 0 };
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_RAWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_
