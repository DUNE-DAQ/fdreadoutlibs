/**
 * @file TDEEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_

#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

class TDEEthFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::TDEEthTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::TDEEthTypeAdapter>;
  using frameptr = types::TDEEthTypeAdapter*;
  using constframeptr = const types::TDEEthTypeAdapter*;
  using wibframeptr = dunedaq::fddetdataformats::TDEEthFrame*;
  // Channel map function type
  //typedef int (*chan_map_fn_t)(int);

  explicit TDEEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled);

  void start(const appfwk::DAQModule::CommandData_t& args) override;

  void stop(const appfwk::DAQModule::CommandData_t& args) override;

  void conf(const appmodel::DataHandlerModule* conf) override;

protected:
  virtual void generate_opmon_data() override;

  // Internals
  dunedaq::daqdataformats::timestamp_t m_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_current_ts = 0;

  uint16_t m_previous_seq_id = 0;
  uint16_t m_current_seq_id = 0;

  dunedaq::daqdataformats::timestamp_t m_pattern_generator_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_pattern_generator_current_ts = 0;

  bool m_first_ts_missmatch = true;
  bool m_ts_problem_reported = false;
  bool m_ts_error_state = false;
  std::atomic<uint64_t> m_ts_error_ctr{ 0 };

  bool m_first_seq_id_mismatch = true;
  bool m_seq_id_problem_reported = false;
  bool m_seq_id_error_state = false;
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

  /**
   * Pipeline Stage 2.: Do software TPG
   * */

  void find_hits(constframeptr fp);
  //void find_hits(constframeptr fp);


private:
  bool m_first_hit = true;
  std::unique_ptr<tpglibs::TPGenerator> m_tp_generator;
  std::vector<std::pair<std::string, nlohmann::json>> m_tpg_configs;
  uint32_t m_tp_max_width;
  std::set<unsigned int> m_channel_mask_set;
  uint16_t m_tpg_threshold_selected;

  std::map<uint, std::atomic<int>> m_tp_channel_rate_map;

  size_t m_num_msg = 0;
  size_t m_num_push_fail = 0;

  std::atomic<int> m_tpg_hits_count{ 0 };

  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_id; // NOLINT(build/unsigned)
  uint32_t m_slot_id;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)
  bool m_emulator_mode = false;

  uint32_t m_tp_count_limit = 0;
  uint32_t m_frame_count_limit = 0;
  size_t m_current_tp_count = 0;
  size_t m_current_frame_count = 0;
  bool m_tp_limit_enabled = false;
  bool m_frame_limit_enabled = false;

  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;

  // Mapping from expanded AVX register position to offline channel number
  //std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER> m_register_channels;
  std::vector<std::pair<trgdataformats::channel_t, int16_t>> m_channel_plane_numbers;
  std::vector<trigger::TriggerPrimitiveTypeAdapter> m_tpa_vectors[3];

  std::shared_ptr<iomanager::SenderConcept<std::vector<trigger::TriggerPrimitiveTypeAdapter>>> m_tp_sink[3];
  std::shared_ptr<iomanager::SenderConcept<fddetdataformats::TDEEthFrame>> m_err_frame_sink;

  //std::thread m_add_hits_tphandler_thread;

  daqdataformats::SourceID m_sourceid;

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_suppressed_too_long{ 0 };
  std::atomic<uint64_t> m_tps_send_failed{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_
