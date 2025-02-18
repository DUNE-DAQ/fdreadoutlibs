/**
 * @file DAPHNEFrameProcessor.hpp DAPHNE specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"
#include "appmodel/RawDataProcessor.hpp"

#include "fddetdataformats/DAPHNEFrame.hpp"
#include "detchannelmaps/TPCChannelMap.hpp"

#include "trigger/TriggerPrimitiveTypeAdapterPDS.hpp"
#include "fdreadoutlibs/FDReadoutIssues.hpp"


#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"


#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

class DAPHNEFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>;
  using frameptr = types::DAPHNESuperChunkTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::DAPHNEFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)
  using constframeptr = const types::DAPHNESuperChunkTypeAdapter*;

  // Constructor
  DAPHNEFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool post_processing_enabled)
    : datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>(error_registry, post_processing_enabled)
  {}

  // Override config for pipeline setup
  void conf(const appmodel::DataHandlerModule* conf) override;
  void extract_tps( constframeptr fp);
  void start(const nlohmann::json& args) override;
  void stop(const nlohmann::json& args) override; 
  // Algorithm used to form a trigger primitive
  dunedaq::trgdataformats::TriggerPrimitivePDS::Algorithm m_tp_algo = trgdataformats::TriggerPrimitivePDS::Algorithm::kUnknown; 



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
  std::shared_ptr<iomanager::SenderConcept<std::vector<trigger::TriggerPrimitiveTypeAdapterPDS>>> m_tp_sink;
  std::vector<trigger::TriggerPrimitiveTypeAdapterPDS> m_tpa_vector;

private:
  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;
  bool m_first_tp = true;
  std::atomic<int> m_tpg_hits_count{ 0 };
  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_send_failed{ 0 };
  std::atomic<uint64_t> m_frame_counter{ 0 };

  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_id; // NOLINT(build/unsigned)
  uint32_t m_slot_id;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)
  bool m_emulator_mode = false;


};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
