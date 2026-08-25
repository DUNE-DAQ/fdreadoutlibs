/**
 * @file TPGInternalStateHarvester.cpp TPG internal state harvester
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifdef TPGLIBS_ENABLE_STATE_MONITORING
#include "fdreadoutlibs/tpg/TPGInternalStateHarvester.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "logging/Logging.hpp"
#include <cassert>
#include <algorithm>
#include <iostream>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {
namespace fdreadoutlibs {

TPGInternalStateHarvester::~TPGInternalStateHarvester() {
  // Ensure thread is properly stopped
  stop_collection_thread();
}

void TPGInternalStateHarvester::set_processor_references(std::vector<ProcRef> refs) {
  m_processor_references = std::move(refs);
  rebuild_prealloc_caches_(); // if num_pipelines is not set, this will gracefully get empty per-pipeline statistics, which will be filled later by update
}

const std::vector<TPGInternalStateHarvester::ProcRef>&
TPGInternalStateHarvester::get_processor_references() const {
  return m_processor_references;
}

void TPGInternalStateHarvester::update_channel_plane_numbers(
    const std::vector<std::pair<trgdataformats::channel_t,int16_t>>& channel_plane_numbers,
    uint8_t num_channels_per_pipeline,
    uint8_t num_pipelines)
{
  m_channel_plane_numbers_per_pipeline.clear();
  m_channel_plane_numbers_per_pipeline.resize(num_pipelines);

  // Cut in order: each pipeline has num_channels_per_pipeline channels
  for (uint8_t p = 0; p < num_pipelines; ++p) {
    auto begin = channel_plane_numbers.begin() + p * num_channels_per_pipeline;
    auto end   = begin + num_channels_per_pipeline;
    m_channel_plane_numbers_per_pipeline[p] =
      std::vector<std::pair<trgdataformats::channel_t,int16_t>>(begin, end);
  }
  m_num_channels_per_pipeline = num_channels_per_pipeline;
  m_num_pipelines = num_pipelines;
  m_expected_total_channels = static_cast<size_t>(m_num_channels_per_pipeline) * m_num_pipelines;

  rebuild_prealloc_caches_(); // now pipelines number is known, can accurately calculate the expected number of metric items per pipeline
}

void TPGInternalStateHarvester::rebuild_prealloc_caches_()
{
  // pre-clean
  m_metric_items_per_proc.clear();
  m_metric_items_per_proc.reserve(m_processor_references.size());

  // first clear/reset the expected number of metric items per pipeline
  m_expected_items_per_pipeline.clear();
  m_expected_items_per_pipeline.resize(m_num_pipelines, 0);

  // cache the metric names for each processor and accumulate to the corresponding pipeline
  for (const auto& [proc, pipeline_id] : m_processor_references) {
    if (proc) {
      // Use the processor's interface method which delegates to the registry
      auto names = proc->get_requested_internal_state_names();
      if (static_cast<size_t>(pipeline_id) < m_expected_items_per_pipeline.size()) {
        m_expected_items_per_pipeline[pipeline_id] += names.size();
      }
      m_metric_items_per_proc.emplace_back(std::move(names));
    } else {
      m_metric_items_per_proc.emplace_back(); // empty
    }
  }
}


std::unordered_map<trgdataformats::channel_t,
                   std::vector<std::pair<std::string,int16_t>>>
TPGInternalStateHarvester::harvest_once()
{
  std::unordered_map<trgdataformats::channel_t,
                     std::vector<std::pair<std::string,int16_t>>> out;

  // 1. pre-estimate the capacity of the buckets (usually 64)
  if (m_expected_total_channels > 0) {
    out.reserve(m_expected_total_channels);
  } else {
    out.reserve(static_cast<size_t>(m_num_channels_per_pipeline) * m_num_pipelines);
  }

  // defensive: if the cache is not complete for the current pipelines, rebuild it
  if (m_expected_items_per_pipeline.size() != m_num_pipelines) {
    rebuild_prealloc_caches_();
  }
  
  for (size_t i = 0; i < m_processor_references.size(); ++i) {
    const auto& [proc, pipeline_id] = m_processor_references[i];
    if (!proc) {
      continue;
    }

    // 1) get the cached metric names; if empty, fall back to pulling directly from processor
    std::vector<std::string> metric_names_cached;
    if (m_metric_items_per_proc.size() > i) {
      metric_names_cached = m_metric_items_per_proc[i];
    } else {
      // Fallback: use processor's interface method
      metric_names_cached = proc->get_requested_internal_state_names();
    }

    // 2) get the current snapshot
    const auto arr = proc->read_internal_states_as_integer_array();

    // basic consistency: the number of items should be the same as the number of snapshot items
    if (metric_names_cached.size() != arr.m_size) {
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "Processor " << i << " size mismatch: metric_names=" 
                                   << metric_names_cached.size() << " vs array=" << arr.m_size;
      continue;
    }

    // the current pipeline's lane -> (channel, plane)
    assert(static_cast<size_t>(pipeline_id) < m_channel_plane_numbers_per_pipeline.size());
    const auto& chan_plane_vec = m_channel_plane_numbers_per_pipeline[pipeline_id];

    // 3) first reserve the capacity for all channels in the current pipeline in out
    //    so that subsequent emplace_back will not trigger allocation
    std::vector<decltype(out.begin())> iters_for_lanes;
    iters_for_lanes.resize(chan_plane_vec.size());

    const size_t expected_items_here =
      (pipeline_id < m_expected_items_per_pipeline.size())
      ? m_expected_items_per_pipeline[pipeline_id]
      : metric_names_cached.size(); // fallback, use the number of items of the current processor

    for (size_t lane = 0; lane < chan_plane_vec.size(); ++lane) {
      const auto ch = chan_plane_vec[lane].first; // offline channel id
      auto [it, inserted] = out.try_emplace(ch, std::vector<std::pair<std::string,int16_t>>{});
      if (inserted) {
        // only reserve the capacity for the first time the channel is encountered
        it->second.reserve(expected_items_here);
      }
      iters_for_lanes[lane] = it;
    }

    // 4) map the 16-lane array of each metric to the corresponding channel
    for (size_t item = 0; item < arr.m_size; ++item) {
      const auto& lanes = arr.m_data[item]; // std::array<int16_t, 16>
      const auto& name  = metric_names_cached[item];
      const size_t L    = lanes.size();     // e.g. 16

      const size_t up_to = std::min(L, chan_plane_vec.size());
      for (size_t lane = 0; lane < up_to; ++lane) {
        // directly use the iterator, avoid the repeated lookup of the map
        iters_for_lanes[lane]->second.emplace_back(name, lanes[lane]);
      }
    }
  }

  return out;
}

// --- Multi-threaded implementation ---

void TPGInternalStateHarvester::start_collection_thread() {
  std::lock_guard<std::mutex> lock(m_config_mutex);
  
  if (m_thread_running.load()) {
    return; // Already running
  }
  
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Starting internal state collection thread";
  
  // Initialize result container
  m_latest_results.clear();
  
  // Reset thread control flags
  m_thread_should_stop.store(false);
  m_harvest_requested.store(false);
  
  // Start the collection thread
  m_collection_thread = std::thread(&TPGInternalStateHarvester::collection_thread_worker_, this);
  m_thread_running.store(true);
}

void TPGInternalStateHarvester::stop_collection_thread() {
  {
    std::lock_guard<std::mutex> config_lock(m_config_mutex);
    
    if (!m_thread_running.load()) {
      return; // Already stopped
    }
    
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Stopping internal state collection thread";
    
    // Signal thread to stop
    m_thread_should_stop.store(true);
  }
  
  // Notify the collection thread using the correct mutex
  {
    std::lock_guard<std::mutex> collection_lock(m_collection_mutex);
    m_collection_cv.notify_all();
  }
  
  // Wait for thread to finish
  if (m_collection_thread.joinable()) {
    m_collection_thread.join();
  }
  
  m_thread_running.store(false);
  
  // Clear results
  {
    std::lock_guard<std::mutex> lock(m_results_mutex);
    m_latest_results.clear();
  }
}

void TPGInternalStateHarvester::trigger_harvest() {
  if (!m_thread_running.load()) {
    return; // Thread not running
  }
  
  m_harvest_requested.store(true);
  
  // Notify the collection thread using the correct mutex
  {
    std::lock_guard<std::mutex> lock(m_collection_mutex);
    m_collection_cv.notify_all();
  }
}

std::unordered_map<trgdataformats::channel_t,
                   std::vector<std::pair<std::string,int16_t>>>
TPGInternalStateHarvester::get_latest_results() const {
  std::lock_guard<std::mutex> lock(m_results_mutex);
  return m_latest_results; // Return copy under lock (blocking read is acceptable)
}

bool TPGInternalStateHarvester::is_collection_thread_running() const {
  return m_thread_running.load();
}

void TPGInternalStateHarvester::collection_thread_worker_() {
  while (!m_thread_should_stop.load()) {
    std::unique_lock<std::mutex> lock(m_collection_mutex);
    
    // Wait for harvest request or stop signal
    m_collection_cv.wait(lock, [this] {
      return m_harvest_requested.load() || m_thread_should_stop.load();
    });
    
    if (m_thread_should_stop.load()) {
      break;
    }
    
    // Reset the harvest request flag
    m_harvest_requested.store(false);
    
    // Release lock during collection to allow concurrent reads
    lock.unlock();
    
    // Perform the actual harvest (expensive operation in background)
    auto new_results = harvest_once();
    
    // Update results with simple mutex protection
    {
      std::lock_guard<std::mutex> results_lock(m_results_mutex);
      m_latest_results = std::move(new_results);
    }
  }
}

} // namespace fdreadoutlibs
} // namespace dunedaq
#endif // TPGLIBS_ENABLE_STATE_MONITORING
