/**
 * @file TDEEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETH_TDEETHFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETH_TDEETHFRAMEPROCESSOR_HPP_

// #include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/FDReadoutIssues.hpp"
#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"
#include "trigger/TriggerPrimitiveTypeAdapter.hpp"

#include "daqdataformats/Types.hpp"

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
#include <random>


namespace dunedaq {
namespace fdreadoutlibs {

class TDEEthFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::TDEEthTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::TDEEthTypeAdapter>;
  using frameptr = types::TDEEthTypeAdapter*;
  using constframeptr = const types::TDEEthTypeAdapter*;
  using tdeframeptr = dunedaq::fddetdataformats::TDEEthFrame*;

  explicit TDEEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry);

  ~TDEEthFrameProcessor();

  void start(const nlohmann::json& args) override;

  void stop(const nlohmann::json& args) override;

  // void init(const nlohmann::json& args) override;

  void conf(const appmodel::DataHandlerModule* conf) override;

  void get_info(opmonlib::InfoCollector& ci, int level) override;

protected:
  // Internals
  dunedaq::daqdataformats::timestamp_t m_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_current_ts = 0;

  uint16_t m_previous_seq_id = 0;
  uint16_t m_current_seq_id = 0;

  dunedaq::daqdataformats::timestamp_t m_pattern_generator_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_pattern_generator_current_ts = 0;

  bool m_first_ts_missmatch = true;
  bool m_ts_problem_reported = false;
  std::atomic<uint64_t> m_ts_error_ctr{ 0 };

  bool m_first_seq_id_mismatch = true;
  bool m_seq_id_problem_reported = false;
  std::atomic<uint64_t> m_seq_id_error_ctr{ 0 };
  std::atomic<int16_t> m_seq_id_min_jump{ 0 };
  std::atomic<int16_t> m_seq_id_max_jump{ 0 };

  /**
   * Pipeline Stage 0: Pattern generator for hit finding in emulated mode
   * */
  void use_pattern_generator(frameptr fp);

  /**
   * Pipeline Stage 1.: Check proper sequence id increments in DAQ Eth header
   * */

  void sequence_check(frameptr fp);

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAQ Eth header
   * */

  void timestamp_check(frameptr fp);

private:
  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_id; // NOLINT(build/unsigned)
  uint32_t m_slot_id;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)
  bool m_emulator_mode = false;

  std::shared_ptr<iomanager::SenderConcept<trigger::TriggerPrimitiveTypeAdapter>> m_tp_sink;

};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETH_TDEETHFRAMEPROCESSOR_HPP_
