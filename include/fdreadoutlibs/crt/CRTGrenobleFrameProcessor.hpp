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
    explicit CRTGrenobleFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool post_processing_enabled)
    : datahandlinglibs::TaskRawDataProcessorModel<types::CRTGrenobleTypeAdapter>(error_registry, post_processing_enabled)
    {}

    void conf(const appmodel::DataHandlerModule* conf) override;

protected:
    using dunedaq::datahandlinglibs::logging::TLVL_FRAME_RECEIVED;

    void timestamp_check(types::CRTGrenobleTypeAdapter* fp);
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTGRENOBLEFRAMEPROCESSOR_HPP_

