/**
 * @file DAPHNEEthFrameProcessor.hpp DAPHNE specific Task based raw processor
 * implementation
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/DAPHNEEthFrame.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "fdreadoutlibs/daphneeth/DAPHNEEthFrameProcessor.hpp"
#include "fdreadoutlibs/pds/DescriptorToTP.hpp"

#include "confmodel/GeoId.hpp"
#include "confmodel/Queue.hpp"
#include "appmodel/DataMoveCallbackConf.hpp"

#include <atomic>
#include <cstring>
#include <functional>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;
using dunedaq::datahandlinglibs::logging::TLVL_FRAME_RECEIVED;

DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")
DUNE_DAQ_TYPESTRING(std::vector<dunedaq::trigger::TriggerPrimitiveTypeAdapter>, "TriggerPrimitiveVector")

namespace dunedaq {
namespace fdreadoutlibs {

DAPHNEEthFrameProcessor::~DAPHNEEthFrameProcessor()
{
  pds::remove_descriptor_processor(m_descriptor_key, m_descriptor_processor);
}

void 
DAPHNEEthFrameProcessor::conf(const appmodel::DataHandlerModule* conf)
{
  pds::remove_descriptor_processor(m_descriptor_key, m_descriptor_processor);
  m_descriptor_processor.reset();
  bool separate_descriptors = false;
  m_tp_sink.reset();
  m_channel_map.reset();
  m_channel_mask_set.clear();
  TLOG() << "Looking for TP sink...";

  for (auto output : conf->get_outputs()) {
    TLOG() << "On outputs... (" << output->UID() << "," << output->get_data_type() << ")";
    try {
      if (output->get_data_type() == "TriggerPrimitiveVector") {
         TLOG() << "Found TP sink.";
         m_tp_sink = get_iom_sender<std::vector<trigger::TriggerPrimitiveTypeAdapter>>(output->UID());
         TLOG() << " SINK INITIALIZED for TriggerPrimitives with UID : " << output->UID();
      }
    } catch (const ers::Issue& excpt) {
      ers::error(datahandlinglibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
    }
  }

  auto dp = conf->get_module_configuration()->get_data_processor();
  if (dp == nullptr) {
    TLOG()<< " PDS Data processor does not exist.";
  } else {
    auto proc_conf = dp->cast<appmodel::PDSRawDataProcessor>();
    if (proc_conf == nullptr) {
      TLOG()<< "PDS RawDataProcessor does not exist.";
    } else { 
      m_def_adc_intg_thresh = proc_conf-> get_default_adc_intg_thresh();
      separate_descriptors = proc_conf->get_separate_descriptor_processing();
      
      auto geo_id = conf->get_geo_id();
      if (geo_id != nullptr) {
        m_det_id = geo_id->get_detector_id();
        m_crate_id = geo_id->get_crate_id();
        m_slot_id = geo_id->get_slot_id();
        m_stream_id = geo_id->get_stream_id();
      }    
    
      m_channel_map = dunedaq::detchannelmaps::make_pds_map(proc_conf->get_channel_map());
      const std::vector<unsigned int> channel_mask_vec = proc_conf->get_channel_mask();
    
      m_channel_mask_set.clear();
      m_channel_mask_set.insert(channel_mask_vec.begin(), channel_mask_vec.end());
      
    }
  }

  if (m_post_processing_enabled && (!m_tp_sink || !m_channel_map)) {
    throw std::runtime_error("DAPHNE descriptor TPG requires a PDS channel map and TriggerPrimitiveVector output");
  }
  if (m_post_processing_enabled && separate_descriptors) {
    for (auto output : conf->get_outputs()) {
      if (output->get_data_type() != "TriggerPrimitiveVector") continue;
      auto queue = output->cast<confmodel::Queue>();
      if (!queue || queue->get_queue_type() != confmodel::Queue::Queue_type::KFollyMPMCQueue) {
        throw std::runtime_error("Separate descriptor processing requires a kFollyMPMCQueue TP output");
      }
    }
    if (!conf->get_raw_data_callback()) {
      throw std::runtime_error("Separate descriptor processing requires a raw callback");
    }
  }
  TLOG() << "Calling parent conf.";
  inherited::conf(conf);
  // Register only after configuration succeeds, so a failed attempt leaves no callbacks.
  inherited::add_preprocess_task(std::bind(&DAPHNEEthFrameProcessor::timestamp_check, this, std::placeholders::_1));
  if (m_post_processing_enabled && separate_descriptors) {
    auto raw_callback = conf->get_raw_data_callback();
    m_descriptor_key = raw_callback->UID();
    auto map = m_channel_map;
    auto sink = m_tp_sink;
    m_descriptor_processor = std::make_shared<pds::DescriptorProcessor>(
      [map](const pds::Frame& input) {
        return map->get_offline_channel_from_det_crate_slot_stream_chan(
          input.daq_header.det_id, input.daq_header.crate_id, input.daq_header.slot_id,
          input.daq_header.stream_id, input.get_channel());
      }, m_channel_mask_set, m_def_adc_intg_thresh,
      [sink](std::vector<pds::TP>&& primitives) {
        std::vector<trigger::TriggerPrimitiveTypeAdapter> output;
        output.reserve(primitives.size());
        for (const auto& tp : primitives) {
          trigger::TriggerPrimitiveTypeAdapter adapter;
          adapter.tp = tp;
          output.push_back(adapter);
        }
        return sink->try_send(std::move(output), iomanager::Sender::s_no_block);
      });
    pds::register_descriptor_processor(m_descriptor_key, m_descriptor_processor);
  } else if (m_post_processing_enabled) {
    // Pre-process because SkipList latency buffers do not support this post-processing path.
    inherited::add_preprocess_task(std::bind(&DAPHNEEthFrameProcessor::extract_tps, this, std::placeholders::_1));
  }
}

void DAPHNEEthFrameProcessor::start(const appfwk::DAQModule::CommandData_t& args)
{
  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_ts_error_ctr = 0;

  // Reset stats
  m_t0 = std::chrono::high_resolution_clock::now();
  m_num_new_tps.exchange(0);
  m_tps_send_failed.exchange(0);
  m_descriptor_frames_rejected.exchange(0);

  if (m_descriptor_processor) {
    auto& c = m_descriptor_processor->counters;
    c.frames = 0; c.overflow = 0; c.malformed = 0; c.sent = 0; c.send_failed = 0;
  }
  inherited::start(args);
}
void
DAPHNEEthFrameProcessor::stop(const appfwk::DAQModule::CommandData_t& args)
{
  inherited::stop(args);
}
/**
 * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
 * */
void 
DAPHNEEthFrameProcessor::timestamp_check(frameptr fp)
{

/*
  for (size_t i=0; i<types::kDAPHNENumFrames; i++){
    auto df_ptr = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEEthFrame*>(fp);

    if(df_ptr[i].get_timestamp() > 0xFFFFFFFFFFFF0000 || df_ptr[i].get_timestamp() < 0xFFFF){
      ers::warning(PDSUnphysicalFrameTimestamp(ERS_HERE, df_ptr[i].get_timestamp(), df_ptr[i].get_channel(), i));
      // Force the TS to 0
      df_ptr[i].daq_header.timestamp_1 = df_ptr[i].daq_header.timestamp_2 = 0;
    }
  }
*/

  // Acquire timestamp
  m_current_ts = fp->get_timestamp();
  uint64_t k_clock_frequency = 62500000; // NOLINT(build/unsigned)
  TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Received DAPHNE frame timestamp value of " << m_current_ts << " ticks (..." << std::fixed << std::setprecision(8) << (static_cast<double>(m_current_ts % (k_clock_frequency*1000)) / static_cast<double>(k_clock_frequency)) << " sec)";// NOLINT


  if (m_ts_error_ctr > 1000) {
    if (!m_problem_reported) {
      std::cout << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!\n";
      m_problem_reported = true;
    }
  }

  m_previous_ts = m_current_ts;
  m_last_processed_daq_ts = m_current_ts;
}

/**
 * Pipeline Stage 2.: Check DAPHNE headers for error flags
 * */
void 
DAPHNEEthFrameProcessor::frame_error_check(frameptr /*fp*/)
{
  // check error fields
}


void DAPHNEEthFrameProcessor::extract_tps(constframeptr fp)
{
  if (!fp || !m_tp_sink || !m_channel_map) return;
  // The adapter's char buffer need not satisfy the frame's alignment.
  pds::Frame frame{};
  std::memcpy(&frame, fp->data, sizeof(frame));
  std::vector<trgdataformats::TriggerPrimitive> primitives;
  try {
    primitives = pds::descriptor_tps(frame, [this](const pds::Frame& input) {
      return m_channel_map->get_offline_channel_from_det_crate_slot_stream_chan(
        input.daq_header.det_id, input.daq_header.crate_id, input.daq_header.slot_id,
        input.daq_header.stream_id, input.get_channel());
    }, pds::DescriptorConfig{m_def_adc_intg_thresh, 1, true});
  } catch (const std::bad_alloc&) {
    throw; // Resource exhaustion is not a malformed detector frame.
  } catch (const std::exception& error) {
    // One warning per monitoring interval; persistent overflow must not flood ERS.
    if (++m_descriptor_frames_rejected == 1) {
      ers::warning(PDSDescriptorFrameRejected(ERS_HERE, frame.get_timestamp(), error.what()));
    }
    return;
  }
  std::vector<trigger::TriggerPrimitiveTypeAdapter> output;
  for (const auto& tp : primitives) {
    if (m_channel_mask_set.count(tp.channel)) continue;
    trigger::TriggerPrimitiveTypeAdapter adapter;
    adapter.tp = tp;
    output.push_back(adapter);
  }
  if (output.empty()) return;
  const auto count = output.size();
  const uint64_t start = output.front().tp.time_start;
  const uint64_t end = output.back().tp.time_start;
  const uint64_t first_channel = output.front().tp.channel;
  const uint64_t last_channel = output.back().tp.channel;
  if (!m_tp_sink->try_send(std::move(output), iomanager::Sender::s_no_block)) {
    m_tps_send_failed += count;
    ers::warning(FailedToSendTPVector(ERS_HERE, start, first_channel, end, last_channel));
  } else {
    m_num_new_tps += count;
  }
}

void
DAPHNEEthFrameProcessor::generate_opmon_data() {

  //right now, just fill some basic tp info...
  if (m_post_processing_enabled) {
    auto now = std::chrono::high_resolution_clock::now();
    uint64_t num_new_tps = m_descriptor_processor ? m_descriptor_processor->counters.sent.exchange(0) : m_num_new_tps.exchange(0);
    int num_new_tps_suppressed_too_long = 0; // not relevant for PDS TPs
    uint64_t num_new_tps_send_failed = m_descriptor_processor ? m_descriptor_processor->counters.send_failed.exchange(0) : m_tps_send_failed.exchange(0);
    const auto rejected_frames = m_descriptor_frames_rejected.exchange(0);
    double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "TP rate: " << std::to_string(num_new_tps / seconds / 1000.) << " [kHz]";
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new TPs: " << num_new_tps;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Rejected descriptor frames in monitoring interval: " << rejected_frames;
    
    datahandlinglibs::opmon::HitFindingInfo tp_info;
    tp_info.set_rate_tp_hits(num_new_tps / seconds / 1000.);
    
    tp_info.set_num_tps_sent(num_new_tps);
    tp_info.set_num_tps_suppressed_too_long(num_new_tps_suppressed_too_long);
    tp_info.set_num_tps_send_failed(num_new_tps_send_failed);
    
    publish(std::move(tp_info));

    m_t0 = now;

  }

 inherited::generate_opmon_data();
  
}

} // namespace fdreadoutlibs
} // namespace dunedaq
