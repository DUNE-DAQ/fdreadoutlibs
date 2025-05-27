/**
 * @file CRTGrenobleFrameProcessor.hpp CRTGrenoble specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTGRENOBLEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTGRENOBLEFRAMEPROCESSOR_HPP_

#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/CRTGrenobleTypeAdapter.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

class CRTGrenobleFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::CRTGrenobleTypeAdapter>
{
public:
    CRTGrenobleFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled);
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTGRENOBLEFRAMEPROCESSOR_HPP_

