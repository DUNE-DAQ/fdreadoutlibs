/**
 * @file TPCEthFrameProcessor.hpp TPCEth generic task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCETHFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCETHFRAMEPROCESSOR_HPP_

#include "fdreadoutlibs/FDReadoutIssues.hpp"

#include "appmodel/DataHandlerModule.hpp"
#include "appmodel/ProcessingStep.hpp"
#include "appmodel/RawDataProcessor.hpp"
#include "appmodel/SamplesOverThresholdMinima.hpp"
#include "appmodel/TPCRawDataProcessor.hpp"

#include "confmodel/GeoId.hpp"

#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"
#include "datahandlinglibs/opmon/datahandling_info.pb.h"

#include "daqdataformats/Types.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"

#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "tpglibs/TPGenerator.hpp"
#include "trigger/TriggerPrimitiveTypeAdapter.hpp"
#include "trgdataformats/Types.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dunedaq {
namespace fdreadoutlibs {

template <class ReadoutTypeAdapter>
class TPCEthFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<ReadoutTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<ReadoutTypeAdapter>;
  using frameptr = ReadoutTypeAdapter*;
  using constframeptr = const ReadoutTypeAdapter*;
  using tpcframeptr = ReadoutTypeAdapter::FrameType*;

  explicit TPCEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled);

  void start(const appfwk::DAQModule::CommandData_t& args) override;

  void stop(const appfwk::DAQModule::CommandData_t& args) override;

  void conf(const appmodel::DataHandlerModule* conf) override;

  void scrap(const appfwk::DAQModule::CommandData_t& cfg) override;

protected:
  void generate_opmon_data() override;

  void configure_source_and_geo_ids(const appmodel::DataHandlerModule* conf);

  void configure_preprocessing(const appmodel::DataHandlerModule* conf);

  void configure_postprocessing(const appmodel::DataHandlerModule* conf);

  void configure_channel_plane_numbers(const appmodel::TPCRawDataProcessor* proc_conf);

  void configure_find_tps(const appmodel::DataHandlerModule* conf, const appmodel::TPCRawDataProcessor* proc_conf);

  void scrap_source_and_geo_ids();

  void scrap_preprocessing();

  void scrap_postprocessing();

  void scrap_channel_plane_numbers();

  void scrap_find_tps();

  /**
   * Publishes collected processor metrics to opmon, currently called in generate_opmon_data()
   * */
  void publish_processor_metric_to_opmon();

  /**
   * Publishes collected processor metrics to opmon, with aggregation of metrics to summary statistics across physical planes
   * */
  void publish_processor_metric_to_opmon_with_aggregation();

  /**
   * Optimized version that calculates all metric summaries across all planes in a single pass
   * Returns a map of plane_number -> map of metric_name -> summary statistics
   * */
  std::map<int16_t, std::map<std::string, std::tuple<float, int16_t, int16_t, float, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>>
  calculate_all_metric_summaries_across_planes(const std::unordered_map<dunedaq::trgdataformats::channel_t, std::vector<std::pair<std::string, int16_t>>>& metrics);
  /**
   * Pipeline Stage 1.: Check proper sequence id increments in DAQ Eth header
   * */

  void sequence_check(frameptr fp);

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAQ Eth header
   * */

  void timestamp_check(frameptr fp);

  /**
   * Pipeline Stage 2.: Do software TPG
   * */

  void find_tps(constframeptr fp);

  bool m_emulator_mode = false;
  bool m_first_frame = true;

  // Timestamp related variables.
  dunedaq::daqdataformats::timestamp_t m_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_current_ts = 0;

  dunedaq::daqdataformats::timestamp_t m_pattern_generator_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_pattern_generator_current_ts = 0;

  bool m_first_ts_missmatch = true;
  bool m_ts_problem_reported = false;
  bool m_ts_error_state = false;
  std::atomic<uint64_t> m_ts_error_ctr{ 0 };

  // Sequence ID related variables.
  uint16_t m_previous_seq_id = 0;
  uint16_t m_current_seq_id = 0;

  bool m_first_seq_id_mismatch = true;
  bool m_seq_id_problem_reported = false;
  bool m_seq_id_error_state = false;
  std::atomic<uint64_t> m_seq_id_error_ctr{ 0 };
  std::atomic<int16_t> m_seq_id_min_jump{ 0 };
  std::atomic<int16_t> m_seq_id_max_jump{ 0 };

  // TPG related variables.
  std::unique_ptr<tpglibs::TPGenerator> m_tp_generator;
  std::vector<std::pair<std::string, nlohmann::json>> m_tpg_configs;

  std::unordered_map<unsigned int, std::vector<trigger::TriggerPrimitiveTypeAdapter>> m_plane_to_tpa_vector_map;
  std::unordered_map<unsigned int, std::shared_ptr<iomanager::SenderConcept<std::vector<trigger::TriggerPrimitiveTypeAdapter>>>> m_plane_to_tp_sink_map;

  uint32_t m_tp_count_limit = 0;
  uint32_t m_frame_count_limit = 0;
  uint32_t m_current_tp_count = 0;
  uint32_t m_frame_count_at_last_send = 0;

  bool m_tp_limit_enabled = false;
  bool m_frame_limit_enabled = false;

  // TPG: channel variables.
  std::set<unsigned int> m_channel_mask_set;
  std::set<unsigned int> m_plane_numbers_set;
  std::vector<std::pair<trgdataformats::channel_t, int16_t>> m_channel_plane_numbers;
  std::unordered_map<trgdataformats::channel_t, unsigned int> m_channel_plane_map;

  // OpMon related variables.
  bool m_tpg_metric_collect_enabled{false};
  uint32_t m_metric_collect_opmon_period { 128 };

  std::map<uint, std::atomic<int>> m_tp_channel_rate_map;

  std::atomic<uint64_t> m_num_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_suppressed_too_long{ 0 };
  std::atomic<uint64_t> m_tps_send_failed{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;

  std::atomic<uint64_t> m_frame_counter{ 0 };

  // Source & Geo ID related variables.
  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_id; // NOLINT(build/unsigned)
  uint32_t m_slot_id;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)

  daqdataformats::SourceID m_sourceid;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#include "fdreadoutlibs/detail/TPCEthFrameProcessor.hxx"

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCETHFRAMEPROCESSOR_HPP_
