/**
 * @file DAPHNEStreamFrameProcessor.hpp DAPHNE specific Task based raw processor
 * for DAPHNE Streaming mode
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNESTREAMFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNESTREAMFRAMEPROCESSOR_HPP_

#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fddetdataformats/DAPHNEStreamFrame.hpp"
#include "fdreadoutlibs/DAPHNEStreamSuperChunkTypeAdapter.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class DAPHNEStreamFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::DAPHNEStreamSuperChunkTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DAPHNEStreamSuperChunkTypeAdapter>;
  using frameptr = types::DAPHNEStreamSuperChunkTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::DAPHNEStreamFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit DAPHNEStreamFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::TaskRawDataProcessorModel<types::DAPHNEStreamSuperChunkTypeAdapter>(error_registry)
  {}

  // Override config for pipeline setup
  void conf(const appdal::ReadoutModule* conf) override;

protected:
  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
   * */
  void timestamp_check(frameptr /*fp*/);

  /**
   * Pipeline Stage 2.: Check DAPHNE headers for error flags
   * */
  void frame_error_check(frameptr /*fp*/);

  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_fake = true;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };

};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNESTREAMFRAMEPROCESSOR_HPP_
