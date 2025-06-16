/**
 * @file CRTBernFrameProcessor.hpp CRTBern specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/crt/CRTBernFrameProcessor.hpp" // NOLINT(build/include)

namespace dunedaq {
namespace fdreadoutlibs {

    void CRTBernFrameProcessor::conf(const appmodel::DataHandlerModule* /*conf*/)
    {
        TLOG() << "Registering processing tasks...";
        datahandlinglibs::TaskRawDataProcessorModel<types::CRTBernTypeAdapter>::add_preprocess_task(std::bind(&CRTBernFrameProcessor::timestamp_check, this, std::placeholders::_1));
    }

    void CRTBernFrameProcessor::timestamp_check(types::CRTBernTypeAdapter* fp)
    {
        static const uint64_t k_clock_frequency = 62500000; // NOLINT(build/unsigned)
        auto current_ts = fp->get_timestamp();
        TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Received CRTBern frame timestamp value of " << current_ts << " ticks (..." << std::fixed << std::setprecision(8) << (static_cast<double>(current_ts % (k_clock_frequency*1000)) / static_cast<double>(k_clock_frequency)) << " sec)";// NOLINT

        if(current_ts > m_last_processed_daq_ts) m_last_processed_daq_ts = current_ts;
    }


} // namespace fdreadoutlibs
} // namespace dunedaq    
