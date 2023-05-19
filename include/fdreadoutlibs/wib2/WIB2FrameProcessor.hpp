/**
 * @file WIB2FrameProcessor.hpp WIB2 specific Task based raw processor
 * @author Adam Abed Abud (adam.abed.abud@cern.ch)
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIBFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIBFRAMEPROCESSOR_HPP_

// #include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "fdreadoutlibs/wib2/WIB2TPHandler.hpp"
#include "trigger/TPSet.hpp"

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
ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPHandlerBacklog,
                  "Failed to push hits to TP handler " << sid,
                  ((int)sid))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGAlgorithmInexistent,
                  "The selected algorithm does not exist: " << algorithm_selection << " . Check your configuration file and seelect either SWTPG or AbsRS.",
                  ((std::string)algorithm_selection))



namespace fdreadoutlibs {


struct swtpg_output{
  uint16_t* output_location;
  uint64_t timestamp;
};

// Pattern generator class for creating different types of TP patterns
// AAA: at the moment the pattern generator is very simple: random source id and random channel
class WIB2PatternGenerator {

public:
  WIB2PatternGenerator() {
    m_size = 100000;
  };
  ~WIB2PatternGenerator() {};

  void generate(int source_id);

  std::vector<int> get_sourceids() {
    return m_sourceid;
  }

  std::vector<int> get_channels() {
    return m_channel;
  }

  int get_total_size() {
    return m_size;
  }
 
private: 
  int m_size;
  std::vector<int> m_sourceid;
  std::vector<int> m_channel;
  

};


class WIB2FrameHandler {

public: 
  explicit WIB2FrameHandler(int register_selector_params, iomanager::FollyMPMCQueue<uint16_t*> &dest_queue );
  ~WIB2FrameHandler();

  std::unique_ptr<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>> m_tpg_processing_info;

  // Map from expanded AVX register position to offline channel number
  swtpg_wib2::RegisterChannelMap register_channel_map; 

  bool first_hit = true;                                                  
                                                  
  int get_registers_selector();

  void reset();
   

  void initialize(int threshold_value);
  // Pop one destination ptr for the frame handler
  uint16_t* get_primfind_dest();


private: 
  int m_register_selector;    
  iomanager::FollyMPMCQueue<uint16_t*> &m_dest_queue_frame_handler;
  uint16_t m_tpg_threshold;                    // units of sigma // NOLINT(build/unsigned)
  const uint8_t m_tpg_tap_exponent = 6;                  // NOLINT(build/unsigned)
  const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
  std::vector<int16_t> m_tpg_taps;                       // firwin_int(7, 0.1, multiplier);
  int16_t* m_tpg_taps_p = nullptr;
};




class WIB2FrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>;
  using frameptr = types::DUNEWIBSuperChunkTypeAdapter*;
  using constframeptr = const types::DUNEWIBSuperChunkTypeAdapter*;
  using wibframeptr = dunedaq::fddetdataformats::WIB2Frame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Channel map function type
  typedef int (*chan_map_fn_t)(int);

  explicit WIB2FrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry);

  ~WIB2FrameProcessor();

  void start(const nlohmann::json& args) override;

  void stop(const nlohmann::json& args) override;

  void init(const nlohmann::json& args) override;

  void conf(const nlohmann::json& cfg) override;

  void scrap(const nlohmann::json& args) override;

  void get_info(opmonlib::InfoCollector& ci, int level);

protected:
  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };

  timestamp_t m_pattern_generator_previous_ts = 0;
  timestamp_t m_pattern_generator_current_ts = 0;

  /**
   * Pattern generator for hit finding in emulated mode
   * */
  void use_pattern_generator(frameptr fp);


  /**
   * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
   * */
  void timestamp_check(frameptr fp);

  /**
   * Pipeline Stage 2.: Do software TPG
   * */
  void find_hits(constframeptr fp, WIB2FrameHandler* frame_handler);


  unsigned int process_swtpg_hits(uint16_t* primfind_it, timestamp_t timestamp);

  // Push destination ptr to the common queue
  void free_primfind_dest(uint16_t* dest);


  // Function for the TPHandler threads. 
  // Reads the m_tpthread_queue and then process the hits
  void add_hits_to_tphandler();

private:
  bool m_sw_tpg_enabled;
  std::string m_tpg_algorithm;
  std::vector<int> m_channel_mask_vec;
  std::set<uint> m_channel_mask_set;
  uint16_t m_tpg_threshold_selected;
  bool m_emulator_mode_enabled = false;

  std::map<uint, std::atomic<int>> m_tp_channel_rate_map;

  size_t m_num_msg = 0;
  size_t m_num_push_fail = 0;

  std::atomic<int> m_swtpg_hits_count{ 0 };

  std::atomic<bool> m_add_hits_tphandler_thread_should_run;

  uint32_t m_link; // NOLINT(build/unsigned)
  uint32_t m_slot_no;  // NOLINT(build/unsigned)
  uint32_t m_crate_no; // NOLINT(build/unsigned)

  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;

  // Mapping from expanded AVX register position to offline channel number
  std::array<uint, 256> m_register_channels;



  std::shared_ptr<iomanager::SenderConcept<types::TriggerPrimitiveTypeAdapter>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> m_tpset_sink;
  std::shared_ptr<iomanager::SenderConcept<fddetdataformats::WIB2Frame>> m_err_frame_sink;

  std::unique_ptr<WIB2TPHandler> m_tphandler;

  // AAA: TODO: make selection of the initial capacity of the queue configurable
  size_t m_capacity_mpmc_queue = 900000; 
  iomanager::FollyMPMCQueue<swtpg_output> m_tphandler_queue{"tphandler_queue", m_capacity_mpmc_queue};


  // Destination queue represents the queue of primfind
  // destinations after the execution of the SWTPG. We 
  // are using queues because we want to avoid adding a 
  // memcpy. Previously the processing threads were also 
  // performing the TPHandler task so this was not needed
  // (look at ProtoWIB SWTPG implementation). The size of 
  // the queues was set big enough for the March coldbox 
  // tests. The size of each uint16_t buffer is 100 000, 
  // therefore the total impact on memory is: 
  // 100 000 * (16bit = 2 bytes) * 300 000 ~ 55 GB
  size_t m_capacity_dest_queue = m_capacity_mpmc_queue+4; 
  iomanager::FollyMPMCQueue<uint16_t*> m_dest_queue{"dest_queue", m_capacity_dest_queue};


  // Select the registers where to process the frame expansion. 
  // E.g.: 0 --> divide registers by 2 (= 16/2) and select the first half 
  // E.g.: 1 --> divide registers by 2 (= 16/2) and select the second half
  int selection_of_register = 0; 
  std::unique_ptr<WIB2FrameHandler> m_wib2_frame_handler = std::make_unique<WIB2FrameHandler>(selection_of_register, m_dest_queue);
  int selection_of_register_second_half = 1; 
  std::unique_ptr<WIB2FrameHandler> m_wib2_frame_handler_second_half = std::make_unique<WIB2FrameHandler>(selection_of_register_second_half, m_dest_queue);

  // Pattern generator configs
  WIB2PatternGenerator m_wib2_pattern_generator;
  std::vector<int> m_random_channels; 
  int m_pattern_index = 0;

  std::thread m_add_hits_tphandler_thread;

  daqdataformats::SourceID m_sourceid;

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_dropped{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIBFRAMEPROCESSOR_HPP_
