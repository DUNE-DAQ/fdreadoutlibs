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

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

void 
DAPHNEFrameProcessor::conf(const nlohmann::json& args)
{
  readoutlibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>::add_preprocess_task(
    std::bind(&DAPHNEFrameProcessor::timestamp_check, this, std::placeholders::_1));
  // m_tasklist.push_back( std::bind(&DAPHNEFrameProcessor::frame_error_check, this, std::placeholders::_1) );
  TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>::conf(args);
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

} // namespace fdreadoutlibs
} // namespace dunedaq
