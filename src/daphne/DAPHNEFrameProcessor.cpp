/**
 * @file DAPHNEFrameProcessor.hpp DAPHNE specific Task based raw processor
 * implementation
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fddetdataformats/DAPHNEFrame.hpp"
#include "fdreadoutlibs/daphne/DAPHNEFrameProcessor.hpp"
#include "fdreadoutlibs/TriggerPrimitivePDSTypeAdapter.hpp"
#include "fdreadoutlibs/FDReadoutIssues.hpp"

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_FRAME_RECEIVED;


//DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")

namespace dunedaq {
namespace fdreadoutlibs {

void 
DAPHNEFrameProcessor::conf(const nlohmann::json& args)
{
  readoutlibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>::add_preprocess_task(
    std::bind(&DAPHNEFrameProcessor::timestamp_check, this, std::placeholders::_1));
  // m_tasklist.push_back( std::bind(&DAPHNEFrameProcessor::frame_error_check, this, std::placeholders::_1) );
  TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>::conf(args);

  auto config = args["RawDataProcessorConf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();

  m_sourceid.id = config.source_id;
  //m_sourceid.subsystem = types::DUNEWIBEthTypeAdapter::subsystem;
  m_tpg_algorithm = config.tpg_algorithm;    
//  auto dp = conf->get_module_configuration()->get_data_processor();
//  m_channel_map = dunedaq::detchannelmaps::make_map(conf->get_channel_map());
   m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

  inherited::add_postprocess_task(std::bind(&DAPHNEFrameProcessor::extract_tps, this, std::placeholders::_1));
}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
 * */
void 
DAPHNEFrameProcessor::timestamp_check(frameptr fp)
{
  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) { // emulate perfectly incrementing timestamp
    // RS warning : not fixed rate!
    if (m_first_ts_fake) {
      fp->fake_timestamps(m_previous_ts, 16);
      m_first_ts_fake = false;
    } else {
      fp->fake_timestamps(m_previous_ts + 192, 16);
    }
  }

  // Acquire timestamp
  m_current_ts = fp->get_first_timestamp();
  uint64_t k_clock_frequency = 62500000; // NOLINT(build/unsigned)
  TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Received DAPHNE frame timestamp value of " << m_current_ts << " ticks (..." << std::fixed << std::setprecision(8) << (static_cast<double>(m_current_ts % (k_clock_frequency*1000)) / static_cast<double>(k_clock_frequency)) << " sec)"; // NOLINT

  // Check timestamp
  // RS warning : not fixed rate!
  // if (m_current_ts - m_previous_ts != ???) {
  //  ++m_ts_error_ctr;
  //}

  if (m_ts_error_ctr > 1000) {
    if (!m_problem_reported) {
      TLOG() << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!";
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

void
DAPHNEFrameProcessor::extract_tps(constframeptr fp)
{
  size_t nhits = 0;
  if (!fp)
    return;
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEFrame*>((uint8_t*)fp); // NOLINT

  // Check that the system is properly configured from the first hit.
  if (m_first_tp) {
//    if (wfptr->daq_header.crate_id != m_crate_id || wfptr->daq_header.slot_id != m_slot_id || wfptr->daq_header.stream_id != m_stream_id) {
//      ers::error(LinkMisconfiguration(ERS_HERE, wfptr->daq_header.crate_id, wfptr->daq_header.slot_id, wfptr->daq_header.stream_id, m_crate_id, m_slot_id, m_stream_id));
//    }
    m_first_tp = false;
  }

  std::vector<trgdataformats::TriggerPrimitivePDS> tps;

  for(size_t ii=0; ii<5; ii++) tps.push_back(wfptr->get_TP(ii));

  for (auto tp : tps) {
    // If this TP is on a masked channel, skip it -> We might add this in the future
//    if (std::binary_search(m_channel_mask_set.begin(), m_channel_mask_set.end(), tp.channel))
//      continue;
    // Need to move into a type adapter.
    fdreadoutlibs::types::TriggerPrimitivePDSTypeAdapter tpa;
    tpa.tp = tp;
    tpa.tp.detid = m_det_id;  // Last missing piece.
    tpa.tp.algorithm = m_tp_algo;
    if(!m_tp_sink->try_send(std::move(tpa), iomanager::Sender::s_no_block)) {
      ers::warning(FailedToSendTP(ERS_HERE, tp.time_start, tp.channel));
      m_tps_send_failed++;
    } else {
      m_new_tps++;
      ++nhits;
    }
  }
  m_tpg_hits_count += nhits;
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
