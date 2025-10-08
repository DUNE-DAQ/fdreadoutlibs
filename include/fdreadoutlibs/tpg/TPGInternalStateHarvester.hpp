/**
 * @file TPGInternalStateHarvester.hpp TPG internal state harvester
 *
 * This is part of the DUNE DAQ , copyright 2025.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <array>
#include <cstdint>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <chrono>


#include "tpglibs/AbstractProcessor.hpp"
#include "trgdataformats/Types.hpp" // or the header that defines trgdataformats::channel_t

namespace dunedaq {
namespace fdreadoutlibs {

class TPGInternalStateHarvester {
public:
  using ProcRef = std::pair<std::shared_ptr<tpglibs::AbstractProcessor<__m256i>>, int /*pipeline_id*/>;

  // Destructor ensures proper thread cleanup
  ~TPGInternalStateHarvester();

  void set_processor_references(std::vector<ProcRef> refs);
  const std::vector<ProcRef>& get_processor_references() const;

  /**
   * @brief Cuts a full list of (channel, plane) into per-pipeline lists of 16 lanes each
   * 
   * @param channel_plane_numbers 
   * @param num_channels_per_pipeline 
   * @param num_pipelines 
   */
  void update_channel_plane_numbers(const std::vector<std::pair<trgdataformats::channel_t,int16_t>>& channel_plane_numbers,
                                    uint8_t num_channels_per_pipeline,
                                    uint8_t num_pipelines);

  /**
   * @brief Harvest once, outputs channel -> [(metric_name, value)...]
   * 
   * @return std::unordered_map<trgdataformats::channel_t, std::vector<std::pair<std::string,int16_t>>>
   */
  std::unordered_map<trgdataformats::channel_t,
                     std::vector<std::pair<std::string,int16_t>>> harvest_once();

  // --- Multi-threaded interface ---
  
  /**
   * @brief Start the background collection thread
   * Must be called before using async collection features
   */
  void start_collection_thread();
  
  /**
   * @brief Stop the background collection thread
   * Blocks until thread is fully stopped
   */
  void stop_collection_thread();
  
  /**
   * @brief Signal the collection thread to perform one harvest cycle
   * Non-blocking - returns immediately
   */
  void trigger_harvest();
  
  /**
   * @brief Get the latest collected results (thread-safe, non-blocking)
   * Returns a copy of the most recent harvest results
   * 
   * @return std::unordered_map<trgdataformats::channel_t, std::vector<std::pair<std::string,int16_t>>>
   */
  std::unordered_map<trgdataformats::channel_t,
                     std::vector<std::pair<std::string,int16_t>>> get_latest_results() const;
  
  /**
   * @brief Check if collection thread is running
   * 
   * @return true if thread is active, false otherwise
   */
  bool is_collection_thread_running() const;

private:
  
  // --- Original data structures ---
  std::vector<ProcRef> m_processor_references;
  // index: pipeline_id -> vector of (channel, plane) for its 16 lanes
  std::vector<std::vector<std::pair<trgdataformats::channel_t,int16_t>>> m_channel_plane_numbers_per_pipeline;
  uint8_t m_num_channels_per_pipeline;
  uint8_t m_num_pipelines;

  // --- Preallocation caches (rebuilt when refs/channels change) ---
  std::vector<std::vector<std::string>> m_metric_items_per_proc; // size == m_processor_references.size()
  std::vector<size_t> m_expected_items_per_pipeline;             // size == m_num_pipelines
  size_t m_expected_total_channels = 0;                          // == m_num_channels_per_pipeline * m_num_pipelines

  // --- Multi-threaded data structures ---
  using ResultContainer = std::unordered_map<trgdataformats::channel_t,
                                            std::vector<std::pair<std::string,int16_t>>>;
  
  // Single result container with mutex protection
  mutable std::mutex m_results_mutex;
  ResultContainer m_latest_results;
  
  // Thread synchronization
  std::thread m_collection_thread;
  std::atomic<bool> m_thread_should_stop{false};
  std::atomic<bool> m_thread_running{false};
  std::atomic<bool> m_harvest_requested{false};
  
  // Synchronization primitives
  mutable std::mutex m_config_mutex;  // Protects configuration changes
  std::mutex m_collection_mutex;      // Protects collection process
  std::condition_variable m_collection_cv;  // Signals collection thread
  
  // Thread management
  void collection_thread_worker_();
  
  // Rebuild caches after refs or channel layout changes.
  void rebuild_prealloc_caches_();

};

} // namespace fdreadoutlibs
} // namespace dunedaq
