/**
 * @file CRTFrameProcessor.hpp CRT specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFRAMEPROCESSOR_HPP_

#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fddetdataformats/CRTFrame.hpp"
#include "logging/Logging.hpp"
#include "fdreadoutlibs/CRTTypeAdapter.hpp"
#include "readoutlibs/ReadoutLogging.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_FRAME_RECEIVED;

namespace dunedaq {
namespace fdreadoutlibs {

class CRTFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::CRTTypeAdapter>
{
public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::CRTTypeAdapter>;
  using frameptr = types::CRTTypeAdapter*;
  using crtframeptr = dunedaq::fddetdataformats::CRTFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit CRTFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::TaskRawDataProcessorModel<types::CRTTypeAdapter>(error_registry)
  {}

  // Custom pipeline registration
  void conf(const nlohmann::json& args) override;

protected:

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
   * */
  void timestamp_check(frameptr /*fp*/);

  /**
   * Pipeline Stage 2.: Check headers for error flags
   * */
  void frame_error_check(frameptr /*fp*/);

  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_fake = true;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };
private:
};

} // namespace ndreadoutlibs
} // namespace dunedaq


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFRAMEPROCESSOR_HPP_
