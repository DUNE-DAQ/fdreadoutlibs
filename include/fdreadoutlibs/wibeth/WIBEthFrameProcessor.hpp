/**
 * @file WIBEthProcessor.hpp WIBEth specific Task based raw processor
 * @author Giovanna Lehmann Miotto (giovanna.lehmann@cern.ch)
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETH_WIBETHFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETH_WIBETHFRAMEPROCESSOR_HPP_

#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/IterableQueueModel.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"
#include "readoutlibs/readoutconfig/Nljs.hpp"
#include "readoutlibs/readoutinfo/InfoNljs.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"


#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

//#include "fdreadoutlibs/wibeth/WIB2TPHandler.hpp"
#include "rcif/cmd/Nljs.hpp"
//#include "trigger/TPSet.hpp"
//#include "triggeralgs/TriggerPrimitive.hpp"


//#include "tpg/DesignFIR.hpp"
//#include "tpg/FrameExpand.hpp"
//#include "tpg/ProcessAVX2.hpp"
//#include "tpg/ProcessingInfo.hpp"
//#include "tpg/RegisterToChannelNumber.hpp"
//#include "tpg/TPGConstants_wibeth.hpp"

#include <atomic>
#include <bitset>
#include <functional>
#include <future>
#include <memory>
#include <pthread.h>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

namespace dunedaq {
namespace fdreadoutlibs {


class WIBEthFrameHandler {

public: 
  explicit WIBEthFrameHandler(int register_selector_params) {
    m_register_selector = register_selector_params;
  }
  WIBEthFrameHandler(const WIBEthFrameHandler&) = delete;
  WIBEthFrameHandler& operator=(const WIBEthFrameHandler&) = delete;
  ~WIBEthFrameHandler() { }
  void initialize() { }

private: 
  int m_register_selector;    
};




class WIBEthFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>;
  using frameptr = types::DUNEWIBEthTypeAdapter*;
  using constframeptr = const types::DUNEWIBEthTypeAdapter*;
  using wibframeptr = dunedaq::fddetdataformats::WIBEthFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  explicit WIBEthFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>(error_registry)
  {}

  ~WIBEthFrameProcessor(){}

  void start(const nlohmann::json& args) override
  {
    inherited::start(args);
  }

  void stop(const nlohmann::json& args) override
  {
    inherited::stop(args);
  }

  void init(const nlohmann::json& /*args*/) override
  {
  }

  void conf(const nlohmann::json& cfg) override
  {
    auto config = cfg["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
    m_sourceid.id = config.source_id;
    m_sourceid.subsystem = types::DUNEWIBEthTypeAdapter::subsystem;
    m_error_counter_threshold = config.error_counter_threshold;
    m_error_reset_freq = config.error_reset_freq;

    // Setup pre-processing pipeline
    TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>::add_preprocess_task(
      std::bind(&WIBEthFrameProcessor::timestamp_check, this, std::placeholders::_1));

    TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>::conf(cfg);
  }

  void scrap(const nlohmann::json& args) override
  {
    TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>::scrap(args);
  }

  void get_info(opmonlib::InfoCollector& ci, int level)
  {
    readoutlibs::readoutinfo::RawDataProcessorInfo info;

    info.num_frame_errors = m_frame_error_count.exchange(0);

    auto now = std::chrono::high_resolution_clock::now();
    m_t0 = now;

    readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>::get_info(ci, level);
    ci.add(info);
  }

protected:
  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };



  void postprocess_example(const types::DUNEWIBEthTypeAdapter* fp)
  {
    TLOG() << "Postprocessing: " << fp->get_first_timestamp();
  }

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
   * */
  void timestamp_check(frameptr fp)
  {

    uint16_t wibeth_tick_difference = types::DUNEWIBEthTypeAdapter::expected_tick_difference;


    auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(fp); // NOLINT

    // If EMU data, emulate perfectly incrementing timestamp
    if (inherited::m_emulator_mode) {                           // emulate perfectly incrementing timestamp
      uint64_t ts_next = m_previous_ts + wibeth_tick_difference;                   // NOLINT(build/unsigned)
      wfptr->set_timestamp(ts_next);
    }

    // Acquire timestamp
    m_current_ts = wfptr->get_timestamp();

    // Check timestamp
    if (m_current_ts - m_previous_ts != wibeth_tick_difference) {
      ++m_ts_error_ctr;
      m_error_registry->add_error("MISSING_FRAMES",
                                  readoutlibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + wibeth_tick_difference, m_current_ts));
      if (m_first_ts_missmatch) { // log once
        TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISSMATCH! -> | previous: " << std::to_string(m_previous_ts)
                                     << " current: " + std::to_string(m_current_ts);
        m_first_ts_missmatch = false;
      }
    }

    if (m_ts_error_ctr > 1000) {
      if (!m_problem_reported) {
        TLOG() << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
               << "Something is wrong with the FE source or with the configuration!";
        m_problem_reported = true;
      }
    }

    m_previous_ts = m_current_ts;
    m_last_processed_daq_ts = m_current_ts;
  }

  /**
   * Pipeline Stage 2.: Check WIB headers for error flags
   * */
  void frame_error_check(frameptr fp)
  {
    if (!fp)
      return;

    if (m_frames_processed % 10000 == 0) {
      for (int i = 0; i < m_num_frame_error_bits; ++i) {
        if (m_error_occurrence_counters[i])
          m_error_occurrence_counters[i]--;
      }
    }

    m_frames_processed++;
  }


private:

  // Frame error check
  int m_error_counter_threshold;
  const int m_num_frame_error_bits = 16;
  int m_error_occurrence_counters[16] = { 0 };
  int m_error_reset_freq;

  std::atomic<uint64_t> m_frame_error_count{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_frames_processed{ 0 };  // NOLINT(build/unsigned)
  daqdataformats::SourceID m_sourceid;

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBETH_WIBETHFRAMEPROCESSOR_HPP_
