/**
 * @file DAPHNEEthStreamFrameProcessor.hpp DAPHNE specific Task based raw processor
 * for DAPHNE Eth Streaming mode
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETH_DAPHNEETHSTREAMFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETH_DAPHNEETHSTREAMFRAMEPROCESSOR_HPP_

#include "logging/Logging.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"

#include "fddetdataformats/DAPHNEEthStreamFrame.hpp"
#include "fdreadoutlibs/DAPHNEEthStreamTypeAdapter.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class DAPHNEEthStreamFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNEEthStreamTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNEEthStreamTypeAdapter>;
  using frameptr = types::DAPHNEEthStreamTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::DAPHNEEthStreamFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit DAPHNEEthStreamFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool post_processing_enabled)
    : datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNEEthStreamTypeAdapter>(error_registry, post_processing_enabled)
  {}

  // Override config for pipeline setup
  void conf(const appmodel::DataHandlerModule* conf) override;

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

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNEETH_DAPHNEETHSTREAMFRAMEPROCESSOR_HPP_
