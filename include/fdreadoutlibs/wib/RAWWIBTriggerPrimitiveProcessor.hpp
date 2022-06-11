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

#include "fdreadoutlibs/FDReadoutTypes.hpp"
#include "detdataformats/wib/RawWIBTp.hpp"
#include "logging/Logging.hpp"
#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutLogging.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"
#include "fdreadoutlibs/wib/WIBTPHandler.hpp"
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

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using namespace dunedaq::readoutlibs::logging;

namespace dunedaq {
namespace fdreadoutlibs {

class RAWWIBTriggerPrimitiveProcessor
  : public readoutlibs::TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>;
  using frame_ptr = types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT*;
  using rwtp_ptr = detdataformats::wib::RawWIBTp*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Channel map function type
  typedef int (*chan_map_fn_t)(int);
  
  explicit RAWWIBTriggerPrimitiveProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>(error_registry)
    , m_fw_tpg_enabled(false)
  {}

  void conf(const nlohmann::json& args) override
  {
    TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>::add_preprocess_task(
                std::bind(&RAWWIBTriggerPrimitiveProcessor::tp_unpack, this, std::placeholders::_1));
    TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>::conf(args);

    auto config = args["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
    if (config.enable_firmware_tpg) {
      m_fw_tpg_enabled = true;

      TLOG(1) << "IRHRI fwTPG enabled -- before reset tp handler";
      m_tphandler.reset(
            new WIBTPHandler(*m_tp_sink, *m_tpset_sink, config.tp_timeout, config.tpset_window_size, m_geoid, config.tpset_topic));
      TLOG(1) << "IRHRI fwTPG enabled -- after reset tp handler";
    }

    m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);
  }

  void init(const nlohmann::json& args) override
  {
    m_fake_timestamp = 0;
    /*
    try {
      auto queue_index = appfwk::connection_index(args, {});
      if (queue_index.find("tp") != queue_index.end()) {
        m_tp_source = get_iom_receiver<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>(queue_index["tp"]);
      }
    } catch (const ers::Issue& excpt) {
      // error
    }
    */ 

    try {
      auto queue_index = appfwk::connection_index(args, {});
      //if (queue_index.find("tp_out") != queue_index.end()) {
      //  m_tp_sink = get_iom_sender<types::SW_WIB_TRIGGERPRIMITIVE_STRUCT>(queue_index["tp_out"]);
      //}
      if (queue_index.find("tpset_out") != queue_index.end()) {
        m_tpset_sink = get_iom_sender<trigger::TPSet>(queue_index["tpset_out"]);
      }
    } catch (const ers::Issue& excpt) {
      throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
    }
  }

  void start(const nlohmann::json& args) override
  {
    if (m_fw_tpg_enabled) {
      rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
      TLOG(1) << "IRHRI fwTPG enabled -- setting run number";
      m_tphandler->set_run_number(start_params.run);

      TLOG(1) << "IRHRI fwTPG enabled -- going to reset tp handler";
      m_tphandler->reset();
      m_tps_dropped = 0;
    }
  }

  void stop(const nlohmann::json& /*args*/) override
  {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "Number of TP frames " << m_tp_frames;
    TLOG_DEBUG(TLVL_WORK_STEPS) << "Number of TPs stitched " << m_tps_stitched;
    TLOG_DEBUG(TLVL_WORK_STEPS) << "Number of TPs dropped " << m_tps_dropped;
  }

  void scrap(const nlohmann::json& args) override
  {
    TLOG(1) << "IRHRI fwTPG enabled -- before reset tp handler";
    m_tphandler.reset();
    TLOG(1) << "IRHRI fwTPG enabled -- after reset tp handler";

    TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>::scrap(args);
  }

  void get_info(opmonlib::InfoCollector& ci, int level)
  {
    readoutlibs::readoutinfo::RawDataProcessorInfo info;

    if (m_tphandler != nullptr) {
      info.num_tps_sent = m_tphandler->get_and_reset_num_sent_tps();
      info.num_tpsets_sent = m_tphandler->get_and_reset_num_sent_tpsets();
      info.num_tps_dropped = m_tps_dropped.exchange(0);
    }

    readoutlibs::TaskRawDataProcessorModel<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>::get_info(ci, level);
    ci.add(info);
  }

void tp_stitch(rwtp_ptr rwtp)
{
  m_fake_timestamp += 6400;
  TLOG(1) << "IRHRI tp_stitch " << m_fake_timestamp; 
  m_tp_frames++;
  uint64_t ts_0 = rwtp->m_head.get_timestamp(); // NOLINT
  //uint64_t ts_0 = m_fake_timestamp; // NOLINT
  int nhits = rwtp->m_head.get_nhits(); // NOLINT
  uint8_t m_channel_no = rwtp->m_head.m_wire_no; // NOLINT
  uint8_t m_fiber_no = rwtp->m_head.m_fiber_no; // NOLINT
  uint8_t m_crate_no = rwtp->m_head.m_crate_no; // NOLINT
  uint8_t m_slot_no = rwtp->m_head.m_slot_no; // NOLINT
  //uint offline_channel = m_channel_map->get_offline_channel_from_crate_slot_fiber_chan(m_crate_no, m_slot_no, m_fiber_no, m_channel_no);

  TLOG(1) << "IRHRI fwTPG enabled -- will loop over " << nhits << " hits";
  for (int i = 0; i < nhits; i++) {

    triggeralgs::TriggerPrimitive trigprim;
    trigprim.time_start = ts_0 + rwtp->m_blocks[i].m_start_time * m_time_tick;
    trigprim.time_peak = ts_0 + rwtp->m_blocks[i].m_peak_time * m_time_tick;
    trigprim.time_over_threshold = (rwtp->m_blocks[i].m_end_time - rwtp->m_blocks[i].m_start_time) * m_time_tick;
    trigprim.channel = m_channel_no; //offline_channel; // m_channel_no;
    trigprim.adc_integral = rwtp->m_blocks[i].m_sum_adc;
    trigprim.adc_peak = rwtp->m_blocks[i].m_peak_adc;
    trigprim.detid =
            m_fiber_no; // TODO: convert crate/slot/fiber to GeoID Roland Sipos rsipos@cern.ch July-22-2021
    trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
    trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
    trigprim.version = 1;

    TLOG(1) << "IRHRI tp_stitch time_start " << trigprim.time_start; 
    TLOG(1) << "IRHRI tp_stitch time_peak " << trigprim.time_peak; 
    TLOG(1) << "IRHRI tp_stitch time_over_threshold " << trigprim.time_over_threshold; 
    TLOG(1) << "IRHRI tp_stitch channel " << trigprim.channel; 

    // stitch current hit to previous hit
    if (m_A[m_channel_no].size() == 1) {
      if (static_cast<int>(rwtp->m_blocks[i].m_start_time) == 0
          && (
          static_cast<int>(trigprim.time_start) - static_cast<int>(m_A[m_channel_no][0].time_start)
             <= static_cast<int>(m_stitch_constant)
          || (static_cast<int>(trigprim.time_start) - static_cast<int>(m_T[m_channel_no][0])
             <= static_cast<int>(m_stitch_constant))
             )
          ) {
        // current hit is continuation of previous hit
        m_T[m_channel_no].clear();
        if (trigprim.adc_peak > m_A[m_channel_no][0].adc_peak) {
          m_A[m_channel_no][0].time_peak = trigprim.time_peak;
          m_A[m_channel_no][0].adc_peak = trigprim.adc_peak;
        }
        m_A[m_channel_no][0].time_over_threshold += trigprim.time_over_threshold;
        m_A[m_channel_no][0].adc_integral += trigprim.adc_integral;
        m_T[m_channel_no].push_back(trigprim.time_start);

      } else {
        // current hit is not continuation of previous hit
        // add previous hit to TriggerPrimitives
        if (!m_tphandler->add_tp(std::move(m_A[m_channel_no][0]), ts_0)) {
          m_tps_dropped++;
        }
        m_tps_stitched++;
        m_tphandler->try_sending_tpsets(ts_0);
        TLOG(1) << "IRHRI tp_stitch try_sending_tpsets 0 " ; 
        m_A[m_channel_no].clear();
        m_T[m_channel_no].clear();
      }
    }

    // NB for TPSets: this assumes hits come ordered in time 
    // current hit (is, completes or starts) one TriggerPrimitive 
    uint8_t m_tp_continue = rwtp->m_blocks[i].m_hit_continue; // NOLINT
    uint8_t m_tp_end_time = rwtp->m_blocks[i].m_end_time; // NOLINT
 
    if (m_tp_continue == 0 && m_tp_end_time != 63) {
      if (m_A[m_channel_no].size() == 1) {
        // the current hit completes one stitched TriggerPrimitive
        if (!m_tphandler->add_tp(std::move(m_A[m_channel_no][0]), ts_0)) {
          m_tps_dropped++;
        }
        m_tphandler->try_sending_tpsets(ts_0);
        TLOG(1) << "IRHRI tp_stitch try_sending_tpsets 1 " ; 
        m_tps_stitched++;
        m_A[m_channel_no].clear();
        m_T[m_channel_no].clear();
      } else {
        // the current hit is one TriggerTrimitive
        if (!m_tphandler->add_tp(std::move(trigprim), ts_0)) {
          m_tps_dropped++;
        }
        m_tphandler->try_sending_tpsets(ts_0);
        TLOG(1) << "IRHRI tp_stitch try_sending_tpsets 2 " ; 
        TLOG(1) << "IRHRI tp_stitch current time " <<  ts_0; 
        m_tps_stitched++;      
      }
    } else {
      // the current hit starts one TriggerPrimitive
      if (m_A[m_channel_no].size() == 0) {
        m_A[m_channel_no].push_back(trigprim);
        m_T[m_channel_no].push_back(trigprim.time_start);
      } else { // decide to add long TriggerPrimitive even when it doesn't end properly
               // this is rare case and can be removed for efficiency    
        // the current hit is "bad"
        // add one TriggerPrimitive from previous stitched hits except the current hit  
        if ( m_tp_continue == 0 && m_tp_end_time == 63 &&
             static_cast<int>(trigprim.time_start) - static_cast<int>(m_T[m_channel_no][0])
             <= static_cast<int>(m_stitch_constant)) {
          if (!m_tphandler->add_tp(std::move(m_A[m_channel_no][0]), ts_0)) {
            m_tps_dropped++;
          }
          m_tphandler->try_sending_tpsets(ts_0);
          TLOG(1) << "IRHRI tp_stitch try_sending_tpsets 3 " ; 
          m_tps_stitched++;      
          m_A[m_channel_no].clear();
          m_T[m_channel_no].clear();
        }
      }
    }
  }
} // NOLINT (exceeding 80 lines)


//void unpack_tpframe_version_1(frame_ptr fr)
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

  
  TLOG(1) << "IRHRI fwTPG enabled -- tp unpack received " << num_elem << " vytes from FELIX";
  int offset = 0;
  while (offset <= num_elem) {

    if (offset == num_elem) break;

    // Count number of subframes in a TP frame
    int n = 1;
    while (reinterpret_cast<types::TpSubframe*>(((uint8_t*)srcbuffer.data()) // NOLINT
           + offset + (n-1)*RAW_WIB_TP_SUBFRAME_SIZE)->word3 != 0xDEADBEEF) {
      n++;
    }

    int bsize = n * RAW_WIB_TP_SUBFRAME_SIZE;
    std::vector<char> tmpbuffer;
    tmpbuffer.reserve(bsize);
    int nhits = n - 2;
    TLOG(1) << "IRHRI fwTPG enabled -- tp unpack TP frame size " << bsize;

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

    rwtp_ptr rwtp =
         static_cast<dunedaq::detdataformats::wib::RawWIBTp*>( malloc(
         sizeof(dunedaq::detdataformats::wib::TpHeader) + nhits * sizeof(dunedaq::detdataformats::wib::TpData)
         ));

    ::memcpy(static_cast<void*>(&rwtp->m_head),
             static_cast<void*>(tmpbuffer.data() + 0),
             2*RAW_WIB_TP_SUBFRAME_SIZE);

    for (int i=0; i<nhits; i++) {
      ::memcpy(static_cast<void*>(&rwtp->m_blocks[i]),
               static_cast<void*>(tmpbuffer.data() + (2+i)*RAW_WIB_TP_SUBFRAME_SIZE),
               RAW_WIB_TP_SUBFRAME_SIZE);
    }

    // old format lacks number of hits
    rwtp->set_nhits(nhits); // explicitly set number of hits in new format
    TLOG(1) << "IRHRI fwTPG enabled -- before stitch we found " << nhits << " hits";

    // stitch TP hits
    tp_stitch(rwtp);
    offset += (2+nhits)*RAW_WIB_TP_SUBFRAME_SIZE;
  }
}

protected:
  int m_time_tick = 25;

private:
  using source_t = iomanager::ReceiverConcept<types::RAW_WIB_TRIGGERPRIMITIVE_STRUCT>;
  std::shared_ptr<source_t> m_tp_source;

  // unpacking
  static const constexpr std::size_t RAW_WIB_TP_SUBFRAME_SIZE = 12;

  // stitching algorithm
  std::vector<triggeralgs::TriggerPrimitive> m_A[256]; // keep track of TPs to stitch per channel
  std::vector<uint64_t> m_T[256]; // NOLINT // keep track of last stitched start time
  std::atomic<uint64_t> m_tps_stitched { 0 }; // NOLINT
  std::atomic<uint64_t> m_tp_frames  { 0 }; // NOLINT
  uint64_t m_stitch_constant = 1600; // NOLINT  // one packet = 64 * 25 ns

  // interface to DS
  bool m_fw_tpg_enabled;
  std::shared_ptr<iomanager::SenderConcept<types::SW_WIB_TRIGGERPRIMITIVE_STRUCT>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> m_tpset_sink;
  std::unique_ptr<WIBTPHandler> m_tphandler;
  std::atomic<uint64_t> m_tps_dropped{ 0 }; // NOLINT
  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;
  uint64_t m_fake_timestamp{ 0 };

  // info
  std::atomic<uint64_t> m_sent_tps{ 0 }; // NOLINT(build/unsigned)
  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_RAWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_
