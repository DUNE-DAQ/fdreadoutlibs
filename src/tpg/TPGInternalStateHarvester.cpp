/**
 * @file TPGInternalStateHarvester.cpp TPG internal state harvester
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/tpg/TPGInternalStateHarvester.hpp"
#include <cassert>

namespace dunedaq {
namespace fdreadoutlibs {

void TPGInternalStateHarvester::set_processor_references(std::vector<ProcRef> refs) {
  m_processor_references = std::move(refs);
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
}

std::unordered_map<trgdataformats::channel_t,
                   std::vector<std::pair<std::string,int16_t>>>
TPGInternalStateHarvester::harvest_once()
{
  std::unordered_map<trgdataformats::channel_t,
                     std::vector<std::pair<std::string,int16_t>>> out;

  // Pre-allocate (common: 64 channels/stream)
  out.reserve(m_num_channels_per_pipeline * m_num_pipelines);

  for (const auto& [proc, pipeline_id] : m_processor_references) {
    if (!proc) continue;
    // 1) Get metric names
    const auto metric_names = proc->get_metric_items();
    // 2) Get current snapshot
    const auto arr = proc->read_internal_states_as_integer_array();

    // Basic consistency: items number matches snapshot entries number
    if (metric_names.size() != arr.m_size) continue;
    assert(static_cast<size_t>(pipeline_id) < m_channel_plane_numbers_per_pipeline.size());
    const auto& chan_plane_vec = m_channel_plane_numbers_per_pipeline[pipeline_id];
    // Each metric has one 16-lane array
    for (size_t item = 0; item < arr.m_size; ++item) {
      const auto& lanes = arr.m_data[item]; // std::array<int16_t,16>
      const auto& name  = metric_names[item];
      const size_t L    = lanes.size();     // 16

      // Map lane values to the channel corresponding to this pipeline
      // Assumes: lane order in pipeline matches lane order in buffer
      for (size_t lane = 0; lane < L && lane < chan_plane_vec.size(); ++lane) {
        const auto ch = chan_plane_vec[lane].first; // offline channel id
        out[ch].emplace_back(name, lanes[lane]);
      }
    }
  }
  return out;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
