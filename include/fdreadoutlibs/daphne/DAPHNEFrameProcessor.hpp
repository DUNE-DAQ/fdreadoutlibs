/**
 * @file DAPHNEFrameProcessor.hpp DAPHNE specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_

#include "logging/Logging.hpp"

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "fddetdataformats/DAPHNEFrame.hpp"
#include "detchannelmaps/TPCChannelMap.hpp"

#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"


#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class DAPHNEFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>;
  using frameptr = types::DAPHNESuperChunkTypeAdapter*;
  using constframeptr = const types::DAPHNESuperChunkTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::DAPHNEFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Constructor
  explicit DAPHNEFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>(error_registry)
  {}

  // Override config for pipeline setup
  void conf(const nlohmann::json& args) override;
  void extract_tps(constframeptr fp);

protected:
  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
   * */
  void timestamp_check(frameptr /*fp*/);

  /**
   * Pipeline Stage 2.: Check DAPHNE headers for error flags
   * */
  void frame_error_check(frameptr /*fp*/);

  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_fake = true;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };

  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_no; // NOLINT(build/unsigned)
  uint32_t m_slot_no;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)


  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;


private:
  bool m_first_tp = true;
  std::string m_tpg_algorithm;
  std::atomic<int> m_tpg_hits_count{ 0 };
  std::shared_ptr<iomanager::SenderConcept<fdreadoutlibs::types::TriggerPrimitiveTypeAdapter>> m_tp_sink;  // Algorithm used to form a trigger primitive
  dunedaq::trgdataformats::TriggerPrimitive::Algorithm m_tp_algo = trgdataformats::TriggerPrimitive::Algorithm::kUnknown; 

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_suppressed_too_long{ 0 };
  std::atomic<uint64_t> m_tps_send_failed{ 0 };

};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
