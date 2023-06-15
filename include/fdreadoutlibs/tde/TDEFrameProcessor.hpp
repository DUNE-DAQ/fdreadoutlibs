/*
 * @file TDEFrameProcessor.hpp TDE specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDE_TDEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDE_TDEFRAMEPROCESSOR_HPP_

#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fddetdataformats/TDE16Frame.hpp"

#include "fdreadoutlibs/TDEFrameTypeAdapter.hpp"

#include "logging/Logging.hpp"
#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutLogging.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class TDEFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::TDEFrameTypeAdapter
>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::TDEFrameTypeAdapter>;
  using frameptr = types::TDEFrameTypeAdapter*;
  using tdeframeptr = dunedaq::fddetdataformats::TDE16Frame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit TDEFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::TDEFrameTypeAdapter
>(error_registry)
  {}

  // Override configuration to register pipeline
  void conf(const nlohmann::json& args) override;

protected:
  /**
   * Pipeline Stage 1.: Check proper timestamp increments in TDE frames
   * */
  void timestamp_check(frameptr fp);

  /**
   * Pipeline Stage 2.: Check TDE headers for error flags
   * */
  void frame_error_check(frameptr /*fp*/);

  // Internals
  timestamp_t m_previous_ts[64] = {0};
  timestamp_t m_current_ts = 0;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };
  uint64_t m_clock_frequency; // NOLINT(build/unsigned)

private:
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDE_TDEFRAMEPROCESSOR_HPP_
