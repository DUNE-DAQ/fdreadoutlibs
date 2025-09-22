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

#include "tpglibs/AbstractProcessor.hpp"
#include "trgdataformats/Types.hpp" // or the header that defines trgdataformats::channel_t

namespace dunedaq {
namespace fdreadoutlibs {

class TPGInternalStateHarvester {
public:
  using ProcRef = std::pair<std::shared_ptr<tpglibs::AbstractProcessor<__m256i>>, int /*pipeline_id*/>;

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

private:
  std::vector<ProcRef> m_processor_references;
  // index: pipeline_id -> vector of (channel, plane) for its 16 lanes
  std::vector<std::vector<std::pair<trgdataformats::channel_t,int16_t>>> m_channel_plane_numbers_per_pipeline;
  uint8_t m_num_channels_per_pipeline;
  uint8_t m_num_pipelines;
};

} // namespace fdreadoutlibs
} // namespace dunedaq
