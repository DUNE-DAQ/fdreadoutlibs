/**
 * @file CRTGrenobleFrameProcessor.hpp CRTGrenoble specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/crt/CRTGrenobleFrameProcessor.hpp" // NOLINT(build/include)

namespace dunedaq {
namespace fdreadoutlibs {

CRTGrenobleFrameProcessor::CRTGrenobleFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled)
  : TaskRawDataProcessorModel<types::CRTGrenobleTypeAdapter>(error_registry, processing_enabled)
{
}

} // namespace fdreadoutlibs
} // namespace dunedaq    
