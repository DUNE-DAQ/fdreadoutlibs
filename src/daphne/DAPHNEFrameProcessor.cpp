/**
 * @file DAPHNEFrameProcessor.hpp DAPHNE specific Task based raw processor
 * implementation
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fddetdataformats/DAPHNEFrame.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "fdreadoutlibs/daphne/DAPHNEFrameProcessor.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;
using dunedaq::datahandlinglibs::logging::TLVL_FRAME_RECEIVED;

DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")
DUNE_DAQ_TYPESTRING(std::vector<dunedaq::trigger::TriggerPrimitiveTypeAdapter>, "TriggerPrimitiveVector")

namespace dunedaq {
namespace fdreadoutlibs {

void 
DAPHNEFrameProcessor::conf(const appmodel::DataHandlerModule* conf)
{
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
  
  TLOG() << "Registering processing tasks...";
  inherited::add_preprocess_task(std::bind(&DAPHNEFrameProcessor::timestamp_check, this, std::placeholders::_1));
  
  if (m_post_processing_enabled) { 
    // Extract TPs back as a pre-processing task, due to LatencyBuffer post-proc issues using SkipList.
    inherited::add_preprocess_task(std::bind(&DAPHNEFrameProcessor::extract_tps, this, std::placeholders::_1));
  }

  TLOG() << "Calling parent conf.";
  inherited::conf(conf);
}

void DAPHNEFrameProcessor::start(const nlohmann::json& args)
{
  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_ts_error_ctr = 0;

  // Reset stats
  m_t0 = std::chrono::high_resolution_clock::now();
  m_new_hits = 0;
  m_new_tps = 0;
  //m_tpg_hits_count.exchange(0);

  inherited::start(args);
}
void DAPHNEFrameProcessor::stop(const nlohmann::json& args)
{
  inherited::stop(args);
}
/**
 * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
 * */
void 
DAPHNEFrameProcessor::timestamp_check(frameptr fp)
{
  // Let Source Emulator deal with this
  /*
  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) { // emulate perfectly incrementing timestamp
    // RS warning : not fixed rate!
    if (m_first_ts_fake) {
      fp->fake_timestamps(m_previous_ts, 16);
      m_first_ts_fake = false;
    } else {
      fp->fake_timestamps(m_previous_ts + 192, 16);
    }
  }*/

  // Acquire timestamp
  m_current_ts = fp->get_timestamp();
  uint64_t k_clock_frequency = 62500000; // NOLINT(build/unsigned)
  TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Received DAPHNE frame timestamp value of " << m_current_ts << " ticks (..." << std::fixed << std::setprecision(8) << (static_cast<double>(m_current_ts % (k_clock_frequency*1000)) / static_cast<double>(k_clock_frequency)) << " sec)";// NOLINT

  // Check timestamp
  // RS warning : not fixed rate!
  // if (m_current_ts - m_previous_ts != ???) {
  //  ++m_ts_error_ctr;
  //}

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
DAPHNEFrameProcessor::frame_error_check(frameptr /*fp*/)
{
  // check error fields
}


void DAPHNEFrameProcessor::extract_tps(constframeptr fp)
{

//  size_t nhits = 0;
  if (!fp || fp==nullptr)
    return;

  auto nonconstframeptr = const_cast<frameptr>(fp);
  auto df_ptr = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEFrame*>((uint8_t*)nonconstframeptr); // NOLINT

  std::vector<trigger::TriggerPrimitiveTypeAdapter> ttpp;
  for (size_t i=0; i<types::kDAPHNENumFrames; i++)
  {
    for(size_t j=0; j<fddetdataformats::DAPHNEFrame::PeakDescriptorData::max_peaks;j++)
    {
      if(df_ptr[i].peaks_data.is_found(j))
      {
        trigger::TriggerPrimitiveTypeAdapter tpa;
        tpa.tp = peak_to_tp(df_ptr[i],j);

	//check for timestamps that are due to frame timestamps ~ ts=0, and ignore these peaks
	if(tpa.tp.time_start > 0xFFFFFFFFFFFF0000 || tpa.tp.time_start < 0xFFFF){
	  ers::warning(PDSPeakIgnored(ERS_HERE, tpa.tp.time_start, tpa.tp.channel, i, j));
	  continue;
	}
	
        tpa.tp.detid = df_ptr->daq_header.det_id;
        ttpp.push_back(tpa);
      }
    }
  }

  int new_tps = ttpp.size();
  if (new_tps > 0) {

    const auto s_ts_begin = ttpp.front().tp.time_start;
    const auto channel_begin = ttpp.front().tp.channel;
    const auto s_ts_end = ttpp.back().tp.time_start;
    const auto channel_end = ttpp.back().tp.channel;      
    
    if (!m_tp_sink->try_send(std::move(ttpp), iomanager::Sender::s_no_block)) {
      ers::warning(FailedToSendTPVector(ERS_HERE, s_ts_begin, channel_begin, s_ts_end, channel_end));
      m_tps_send_failed += new_tps;
    } else {
      m_new_tps += new_tps;
      m_new_hits += new_tps;
    }
  }

  return;
}

dunedaq::trgdataformats::TriggerPrimitive 
DAPHNEFrameProcessor::peak_to_tp(dunedaq::fddetdataformats::DAPHNEFrame &frame, int i)
{
  dunedaq::trgdataformats::TriggerPrimitive tp;
  // TODO: add check on peak presence
  tp.version = frame.version;
  tp.time_start = frame.get_timestamp()+frame.peaks_data.get_sample_start(i);
  tp.samples_to_peak = frame.peaks_data.get_sample_max(i);
  tp.samples_over_threshold = frame.peaks_data.get_samples_over_baseline(i);
  // FIXME : hard-coded channel map
  // WARNING: slot ids in DAPHNEs are all 0!
  tp.channel = frame.daq_header.slot_id*100+frame.get_channel();
  tp.adc_integral = frame.peaks_data.get_adc_integral(i);
  tp.adc_peak = frame.peaks_data.get_adc_max(i);
  tp.detid = dunedaq::trgdataformats::INVALID_DETID;
  return tp;
}


void
DAPHNEFrameProcessor::generate_opmon_data()
{

  //right now, just fill some basic tp info...
  if (m_post_processing_enabled) {
    auto now = std::chrono::high_resolution_clock::now();
    int new_hits = m_new_hits.exchange(0);
    int new_tps = m_new_tps.exchange(0);
    int new_tps_suppressed_too_long = 0; // not relevant for PDS TPs
    int new_tps_send_failed = m_tps_send_failed.exchange(0);
    double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new hits: " << new_hits << " new TPs: " << new_tps;
    
    datahandlinglibs::opmon::HitFindingInfo tp_info;
    tp_info.set_rate_tp_hits(new_hits / seconds / 1000.);
    
    tp_info.set_num_tps_sent(new_tps);
    tp_info.set_num_tps_suppressed_too_long(new_tps_suppressed_too_long);
    tp_info.set_num_tps_send_failed(new_tps_send_failed);
    
    publish(std::move(tp_info));

    m_t0 = now;

  }
  
}
  
} // namespace fdreadoutlibs
} // namespace dunedaq
