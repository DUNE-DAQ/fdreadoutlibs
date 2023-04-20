/**
 * @file SWWIB2TriggerPrimitiveProcessor.hpp WIB TP specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_SWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_SWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_

#include "appfwk/DAQModuleHelper.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "trigger/TPSet.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include "fddetdataformats/WIB2Frame.hpp"


#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class SWWIB2TriggerPrimitiveProcessor
  : public readoutlibs::TaskRawDataProcessorModel<types::SW_WIB2_TRIGGERPRIMITIVE_STRUCT>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::SW_WIB2_TRIGGERPRIMITIVE_STRUCT>;
  using frameptr = types::SW_WIB2_TRIGGERPRIMITIVE_STRUCT*;
  using wibframeptr = dunedaq::fddetdataformats::WIB2Frame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit SWWIB2TriggerPrimitiveProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::SW_WIB2_TRIGGERPRIMITIVE_STRUCT>(error_registry)
  {}

private:
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_SWWIBTRIGGERPRIMITIVEPROCESSOR_HPP_
