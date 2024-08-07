/**
 * @file SSPFrameProcessor.hpp SSP specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSP_SSPFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSP_SSPFRAMEPROCESSOR_HPP_

//#include "appfwk/DAQModuleHelper.hpp"
#include "logging/Logging.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/IterableQueueModel.hpp"
#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"
#include "datahandlinglibs/utils/ReusableThread.hpp"

#include "fddetdataformats/SSPTypes.hpp"

#include "fdreadoutlibs/SSPFrameTypeAdapter.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class SSPFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::SSPFrameTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::SSPFrameTypeAdapter>;
  using frameptr = types::SSPFrameTypeAdapter*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Channel map funciton type
  typedef int (*chan_map_fn_t)(int);

  explicit SSPFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry)
    : datahandlinglibs::TaskRawDataProcessorModel<types::SSPFrameTypeAdapter>(error_registry)
  {
    // Setup pre-processing pipeline
    datahandlinglibs::TaskRawDataProcessorModel<types::SSPFrameTypeAdapter>::add_preprocess_task(
      std::bind(&SSPFrameProcessor::timestamp_check, this, std::placeholders::_1));
  }

  ~SSPFrameProcessor() {}

  void start(const nlohmann::json& args) override { inherited::start(args); }

  void stop(const nlohmann::json& args) override { inherited::stop(args); }

  void conf(const appmodel::DataHandlerModule* conf) override
  {
    // Setup pre-processing pipeline
    datahandlinglibs::TaskRawDataProcessorModel<types::SSPFrameTypeAdapter>::add_preprocess_task(
      std::bind(&SSPFrameProcessor::timestamp_check, this, std::placeholders::_1));

    inherited::conf(conf);
  }

  #warning MISSING OPMON
  // void get_info(opmonlib::InfoCollector& /*ci*/, int /*level*/) {}

  void timestamp_check(frameptr fp)
  {
    // TLOG() << "Got frame with timestamp: " << fp->get_timestamp();
    inherited::m_last_processed_daq_ts = fp->get_first_timestamp();
  }

protected:
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_SSP_SSPFRAMEPROCESSOR_HPP_
