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
  
  // Extract TPs back as a pre-processing task, due to LatencyBuffer post-proc issues using SkipList.
  inherited::add_preprocess_task(std::bind(&DAPHNEFrameProcessor::extract_tps, this, std::placeholders::_1));

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

  //std::cout << "wfptr timestamp: " << fp->get_timestamp() << '\n';

  auto nonconstframeptr = const_cast<frameptr>(fp);
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEFrame*>((uint8_t*)nonconstframeptr); // NOLINT

  //std::cout << "wfptr timestamp: " << wfptr->get_timestamp() << '\n';

  //std::vector<trigger::TriggerPrimitiveTypeAdapter> ttpp;
  for (size_t i=0; i<dunedaq::fdreadoutlibs::types::kDAPHNENumFrames;i++)
  {
    for(size_t j=0; j<5;j++)
    {
      if(wfptr[i].get_da(j)==1)
      {
        trigger::TriggerPrimitiveTypeAdapter tpa;
        tpa.tp = get_TP(wfptr[i],j);

//        tpa.tp.detid = m_det_id;  // Missing piece.
//        tpa.tp.algorithm = m_tp_algo; // to be filled
        //ttpp.push_back(tpa);

        // 27-Mar-2025, KAB: this local vector is a temporary change to get things working!
        // I imagine that there can/should be better grouping of TPs into a vector.
        std::vector<trigger::TriggerPrimitiveTypeAdapter> tptav;
        tptav.push_back(tpa);
        if (!m_tp_sink->try_send(std::move(tptav), iomanager::Sender::s_no_block)) {
          //std::cout << "sind failed " << std::endl;
          //ers::warning(FailedToSendTP(ERS_HERE, s_ts_begin, channel_begin, s_ts_end, channel_end));
          m_tps_send_failed++;
        } else {
          //std::cout << "send success" << std::endl;
          m_new_tps++;
          m_new_hits++;
        }

      }
    }
  }

  /*
  int new_tps = ttpp.size();
  if (!m_tp_sink->try_send(std::move(ttpp), iomanager::Sender::s_no_block)) {
   //std::cout << "sind failed " << std::endl;
   //ers::warning(FailedToSendTP(ERS_HERE, s_ts_begin, channel_begin, s_ts_end, channel_end));
    m_tps_send_failed++;
  } else {
	  //std::cout << "send success" << std::endl;
    m_new_tps += new_tps;
    nhits += new_tps;
  }
  */
  return;
}

dunedaq::trgdataformats::TriggerPrimitive DAPHNEFrameProcessor::get_TP(dunedaq::fddetdataformats::DAPHNEFrame &frame, int i)
{
  dunedaq::trgdataformats::TriggerPrimitive tp;
  tp.version = frame.version;
  tp.time_start = frame.get_timestamp()+frame.get_time_start(i);
  tp.samples_to_peak = frame.get_time_peak(i);
  tp.samples_over_threshold = frame.get_time_over_baseline(i);
  tp.channel = frame.daq_header.slot_id*100+frame.get_channel();
  tp.adc_integral = frame.get_adc_integral(i);
  tp.adc_peak = frame.get_adc_peak(i);
  tp.detid = dunedaq::trgdataformats::INVALID_DETID;
  return tp;
}


void
DAPHNEFrameProcessor::generate_opmon_data()
{
}

} // namespace fdreadoutlibs
} // namespace dunedaq
