/**
 * @file CRTFixedSizeFrameProcessor.hpp CRT specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFIXEDSIZEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFIIXEDSIZEFRAMEPROCESSOR_HPP_

#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fddetdataformats/CRTFixedSizeFrame.hpp"

#include "fdreadoutlibs/CRTFixedSizeTypeAdapter.hpp"


#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class CRTFixedSizeFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::CRTFixedSizeTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::CRTFixedSizeTypeAdapter>;
  using frameptr = types::CRTFixedSizeTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::CRTFixedSizeFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Constructor
  explicit CRTFixedSizeFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::TaskRawDataProcessorModel<types::CRTFixedSizeTypeAdapter>(error_registry)
  {}

  // Override config for pipeline setup
  void conf(const nlohmann::json& args) override;

protected:
  /**
   * Pipeline Stage 1.: Check proper timestamp increments in CRT frame
   * */
  void timestamp_check(frameptr /*fp*/);

  /**
   * Pipeline Stage 2.: Check CRT headers for error flags
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

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTFIXEDSIZEFRAMEPROCESSOR_HPP_
