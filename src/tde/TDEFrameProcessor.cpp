/*
 * @file TDEFrameProcessor.cpp TDE specific Task based raw processor
 * implementation
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fddetdataformats/TDE16Frame.hpp"
#include "fdreadoutlibs/tde/TDEFrameProcessor.hpp"

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_FRAME_RECEIVED;

namespace dunedaq {
namespace fdreadoutlibs {

void 
TDEFrameProcessor::conf(const nlohmann::json& args)
{
  inherited::add_preprocess_task(
      std::bind(&TDEFrameProcessor::timestamp_check, this, std::placeholders::_1));
  // m_tasklist.push_back( std::bind(&TDEFrameProcessor::frame_error_check, this, std::placeholders::_1) );
  inherited::conf(args);

  auto config = args["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
  m_clock_frequency = config.clock_speed_hz;
}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in TDE frames
 * */
void 
TDEFrameProcessor::timestamp_check(frameptr fp)
{
  
  //auto tdef.= reinterpret_cast<dunedaq::fddetdataformats::TDE16Frame*>(fp); // NOLINT
  auto tdef = fp->data; // NOLINT
  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) {         // emulate perfectly incrementing timestamp
    if (m_previous_ts[tdef.get_channel()] == 0) 
	    m_previous_ts[tdef.get_channel()] = tdef.get_timestamp();
    auto ts_next = m_previous_ts[tdef.get_channel()] + (dunedaq::fddetdataformats::ticks_between_adc_samples * dunedaq::fddetdataformats::tot_adc16_samples); // NOLINT(build/unsigned)
    tdef.set_timestamp(ts_next);
  }

  // Acquire timestamp
  m_current_ts = tdef.get_timestamp();
  auto ch = tdef.get_channel();
  auto tdefh = tdef.get_daq_header();
  TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Checking TDE frame timestamp value of " << m_current_ts 
	  << " , crate " << tdefh->crate_id << ", slot " << tdefh->slot_id << ", stream " << tdefh->stream_id; // NOLINT

  // Check timestamp
  if (m_previous_ts[ch]!=0 && m_current_ts - m_previous_ts[ch] != dunedaq::fddetdataformats::ticks_between_adc_samples * dunedaq::fddetdataformats::tot_adc16_samples) {
    ++m_ts_error_ctr;
    m_error_registry->add_error("MISSING_FRAMES",
                                readoutlibs::FrameErrorRegistry::ErrorInterval(m_previous_ts[ch] + (dunedaq::fddetdataformats::ticks_between_adc_samples * dunedaq::fddetdataformats::tot_adc16_samples), m_current_ts));
    if (m_first_ts_missmatch) { // log once
      //TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISSMATCH! -> | previous: " << std::to_string(m_previous_ts[ch])
      TLOG() << "First timestamp MISSMATCH for channel " << ch<< "! -> | previous: " << std::to_string(m_previous_ts[ch])
                                   << " current: " + std::to_string(m_current_ts);
      m_first_ts_missmatch = false;
    }
  }

  if (m_ts_error_ctr > 1000) {
    if (!m_problem_reported) {
      TLOG() << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!";
      m_problem_reported = true;
    }
  }

  m_previous_ts[ch] = m_current_ts;
  m_last_processed_daq_ts = m_current_ts;
}

/**
 * Pipeline Stage 2.: Check TDE headers for error flags
 * */
void 
TDEFrameProcessor::frame_error_check(frameptr /*fp*/)
{
  // check error fields
}

} // namespace fdreadoutlibs
} // namespace dunedaq
