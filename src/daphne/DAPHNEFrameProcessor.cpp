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

namespace dunedaq {
namespace fdreadoutlibs {

void 
DAPHNEFrameProcessor::conf(const appmodel::DataHandlerModule* conf)
{
  TLOG() << "Looking for TP sink...";

  for (auto output : conf->get_outputs()) {
    TLOG() << "On outputs...";
    try {
      if (output->get_data_type() == "TriggerPrimitive") {
         TLOG() << "Found TP sink.";
         m_tp_sink = get_iom_sender<trigger::TriggerPrimitiveTypeAdapter>(output->UID());
         TLOG() << " SINK INITIALIZAED for TriggerPrimitives with UID : " << output->UID();
      }
    } catch (const ers::Issue& excpt) {
      ers::error(datahandlinglibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
    }
  }
  
  TLOG() << "Registering processing tasks...";
  inherited::add_preprocess_task(std::bind(&DAPHNEFrameProcessor::timestamp_check, this, std::placeholders::_1));
  inherited::add_postprocess_task(std::bind(&DAPHNEFrameProcessor::extract_tps, this, std::placeholders::_1));

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


void DAPHNEFrameProcessor::extract_tps(constframeptr fp)
{

  size_t nhits = 0;
  if (!fp || fp==nullptr)
    return;

  auto nonconstframeptr = const_cast<frameptr>(fp);
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::DAPHNEFrame*>((uint8_t*)nonconstframeptr); // NOLINT

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

        if (!m_tp_sink->try_send(std::move(tpa), iomanager::Sender::s_no_block)) {
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
  tp.time_start = frame.get_timestamp()+64;
  //std::cout << "TIME START: " << (unsigned)tp.time_start << '\n';
  tp.time_peak = frame.get_timestamp()+64+frame.get_time_peak(i);
  //std::cout << "TIME PEAK: " << (unsigned)tp.time_peak << '\n';
  tp.time_over_threshold = frame.get_timestamp()+frame.get_time_over_baseline(i);
  tp.channel = frame.daq_header.slot_id*100+frame.get_channel();
  tp.adc_integral = frame.get_adc_integral(i);
  tp.adc_peak = frame.get_adc_peak(i);
  tp.detid = dunedaq::trgdataformats::INVALID_DETID;
  tp.type = dunedaq::trgdataformats::TriggerPrimitive::Type::kPDS;
  tp.algorithm = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kUnknown;
  return tp;
}


void
DAPHNEFrameProcessor::generate_opmon_data()
{
}

} // namespace fdreadoutlibs
} // namespace dunedaq
