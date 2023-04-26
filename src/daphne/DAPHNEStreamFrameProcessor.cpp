/**
 * @file DAPHNEStreamFrameProcessor.hpp DAPHNE specific Task based raw processor
 * implementation for streaming mode
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fddetdataformats/DAPHNEStreamFrame.hpp"
#include "fdreadoutlibs/daphne/DAPHNEStreamFrameProcessor.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

void 
DAPHNEStreamFrameProcessor::conf(const nlohmann::json& args)
{
  readoutlibs::TaskRawDataProcessorModel<types::DAPHNEStreamSuperChunkTypeAdapter>::add_preprocess_task(
    std::bind(&DAPHNEStreamFrameProcessor::timestamp_check, this, std::placeholders::_1));
  // m_tasklist.push_back( std::bind(&DAPHNEStreamFrameProcessor::frame_error_check, this, std::placeholders::_1) );
  TaskRawDataProcessorModel<types::DAPHNEStreamSuperChunkTypeAdapter>::conf(args);
}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
 * */
void 
DAPHNEStreamFrameProcessor::timestamp_check(frameptr fp)
{
  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) { // emulate perfectly incrementing timestamp
    uint64_t ts_next = m_previous_ts + 64; // NOLINT(build/unsigned)
    auto df = reinterpret_cast<daphneframeptr>(((uint8_t*)fp));  // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) { // NOLINT(build/unsigned)
      //auto wfh = const_cast<dunedaq::fddetdataformats::WIB2Header*>(wf->get_wib_header());
      df->set_timestamp(ts_next);
      ts_next += 64;
      df++;
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
DAPHNEStreamFrameProcessor::frame_error_check(frameptr /*fp*/)
{
  // check error fields
}

} // namespace fdreadoutlibs
} // namespace dunedaq
