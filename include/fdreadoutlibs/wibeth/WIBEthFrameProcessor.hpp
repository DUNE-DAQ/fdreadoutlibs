/**
 * @file WIBEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_

// #include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "fdreadoutlibs/FDReadoutIssues.hpp"
//#include "fdreadoutlibs/wibeth/WIBEthTPHandler.hpp"
//#include "trigger/TPSet.hpp"

#include "daqdataformats/Types.hpp"

#include "tpg/ProcessingInfo.hpp"
#include "tpg/RegisterToChannelNumber.hpp"

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

class WIBEthFrameHandler {

public: 
  explicit WIBEthFrameHandler();
  ~WIBEthFrameHandler();
  std::unique_ptr<swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>> m_tpg_processing_info;

  // Map from expanded AVX register position to offline channel number
  swtpg_wibeth::RegisterChannelMap register_channel_map; 

  bool first_hit = true;                                                  
                                                  
  int get_registers_selector();

  void reset();

  void initialize(int threshold_value);
 
  uint16_t* get_hits_dest();
private: 
  int m_register_selector;    
  uint16_t* m_hits_dest ;
  uint16_t m_tpg_threshold;                    // units of sigma // NOLINT(build/unsigned)
  const uint8_t m_tpg_tap_exponent = 6;                  // NOLINT(build/unsigned)
  const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
  std::vector<int16_t> m_tpg_taps;                       // firwin_int(7, 0.1, multiplier);
  int16_t* m_tpg_taps_p = nullptr;
};

class WIBEthFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>;
  using frameptr = types::DUNEWIBEthTypeAdapter*;
  using constframeptr = const types::DUNEWIBEthTypeAdapter*;
  using wibframeptr = dunedaq::fddetdataformats::WIBEthFrame*;
  // Channel map function type
  typedef int (*chan_map_fn_t)(int);

  explicit WIBEthFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry);

  ~WIBEthFrameProcessor();

  void start(const nlohmann::json& args) override;

  void stop(const nlohmann::json& args) override;

  void init(const nlohmann::json& args) override;

  void conf(const nlohmann::json& cfg) override;

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

  /**
   * Pipeline Stage 2.: Do software TPG
   * */

  void find_hits(constframeptr fp, WIBEthFrameHandler* frame_handler);
  //void find_hits(constframeptr fp);


  void process_swtpg_hits(uint16_t* primfind_it, dunedaq::daqdataformats::timestamp_t timestamp);

private:
  bool m_tpg_enabled;
  std::string m_tpg_algorithm;
  uint32_t m_tp_max_width;
  std::vector<int> m_channel_mask_vec;
  std::set<uint> m_channel_mask_set;
  uint16_t m_tpg_threshold_selected;

  std::map<uint, std::atomic<int>> m_tp_channel_rate_map;

  size_t m_num_msg = 0;
  size_t m_num_push_fail = 0;

  std::atomic<int> m_tpg_hits_count{ 0 };

  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_no; // NOLINT(build/unsigned)
  uint32_t m_slot_no;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)

  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;

  // Mapping from expanded AVX register position to offline channel number
  std::array<uint, swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER> m_register_channels;

    std::function<void(swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>& info)> m_assigned_tpg_algorithm_function;

  std::shared_ptr<iomanager::SenderConcept<fdreadoutlibs::types::ExtendedTriggerPrimitiveTypeAdapter>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<fddetdataformats::WIBEthFrame>> m_err_frame_sink;
  std::unique_ptr<WIBEthFrameHandler> m_wibeth_frame_handler = std::make_unique<WIBEthFrameHandler>();

  //std::thread m_add_hits_tphandler_thread;

  daqdataformats::SourceID m_sourceid;

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_dropped{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_
