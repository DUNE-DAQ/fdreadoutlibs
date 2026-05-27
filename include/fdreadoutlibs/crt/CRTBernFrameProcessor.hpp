/**
 * @file CRTBernFrameProcessor.hpp CRTBern specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTBERNFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTBERNFRAMEPROCESSOR_HPP_

#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/CRTBernTypeAdapter.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

class CRTBernFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::CRTBernTypeAdapter>
{
public:
    explicit CRTBernFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool post_processing_enabled)
            : datahandlinglibs::TaskRawDataProcessorModel<types::CRTBernTypeAdapter>(error_registry, post_processing_enabled)
    {}

    void conf(const appmodel::DataHandlerModule* conf) override;

protected:
    using dunedaq::datahandlinglibs::logging::TLVL_FRAME_RECEIVED;

    void timestamp_check(types::CRTBernTypeAdapter* fp);
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_CRT_CRTBERNFRAMEPROCESSOR_HPP_

