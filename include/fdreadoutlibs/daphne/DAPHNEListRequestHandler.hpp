/**
 * @file DAPHNEListRequestHandler.hpp Trigger matching mechanism for
 * DAPHNE frames with SkipList latency buffer
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNELISTREQUESTHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNELISTREQUESTHANDLER_HPP_

#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/DefaultRequestHandlerModel.hpp"
#include "readoutlibs/models/SkipListLatencyBufferModel.hpp"

#include "fddetdataformats/DAPHNEFrame.hpp"

#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"

#include <atomic>
#include <deque>
#include <functional>
#include <future>
#include <iomanip>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_WORK_STEPS;

namespace dunedaq {
namespace fdreadoutlibs {

class DAPHNEListRequestHandler
  : public readoutlibs::DefaultRequestHandlerModel<
      types::DAPHNESuperChunkTypeAdapter,
      readoutlibs::SkipListLatencyBufferModel<types::DAPHNESuperChunkTypeAdapter>>
{
public:
  // Using shorter typenames
  using inherited =
    readoutlibs::DefaultRequestHandlerModel<types::DAPHNESuperChunkTypeAdapter,
                                            readoutlibs::SkipListLatencyBufferModel<types::DAPHNESuperChunkTypeAdapter>>;
  using SkipListAcc = typename folly::ConcurrentSkipList<types::DAPHNESuperChunkTypeAdapter>::Accessor;
  using SkipListSkip = typename folly::ConcurrentSkipList<types::DAPHNESuperChunkTypeAdapter>::Skipper;

  // Constructor
  DAPHNEListRequestHandler(
    std::unique_ptr<readoutlibs::SkipListLatencyBufferModel<types::DAPHNESuperChunkTypeAdapter>>& latency_buffer,
    std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : DefaultRequestHandlerModel<types::DAPHNESuperChunkTypeAdapter,
                                 readoutlibs::SkipListLatencyBufferModel<types::DAPHNESuperChunkTypeAdapter>>(
        latency_buffer,
        error_registry)
  {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "DAPHNEListRequestHandler created...";
  }

protected:
  // Override cleanup with special implementation
  void cleanup() override { daphne_cleanup_request(); }

  // The DAPHNE with SkipList needs a special cleanup method.
  void daphne_cleanup_request();

private:
  // Constants
  static const constexpr uint64_t m_max_ts_diff = 10000000; // NOLINT(build/unsigned)

  // Stats
  std::atomic<int> m_found_requested_count{ 0 };
  std::atomic<int> m_bad_requested_count{ 0 };
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNELISTREQUESTHANDLER_HPP_
