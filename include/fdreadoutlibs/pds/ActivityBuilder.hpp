/** Streaming, ordered-TP PDS coincidence and descriptor-light scaffolding. */
#ifndef FDREADOUTLIBS_PDS_ACTIVITYBUILDER_HPP_
#define FDREADOUTLIBS_PDS_ACTIVITYBUILDER_HPP_

#include "fdreadoutlibs/pds/DescriptorToTP.hpp"
#include "triggeralgs/TriggerActivity.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <optional>
#include <set>

namespace dunedaq::fdreadoutlibs::pds {

inline uint64_t microseconds_to_ticks(double us, double tick_hz = 62500000.0)
{
  const long double ticks = std::ceil(static_cast<long double>(us) * tick_hz / 1.e6L);
  if (!std::isfinite(us) || !std::isfinite(tick_hz) || us <= 0 || tick_hz <= 0 ||
      ticks >= static_cast<long double>(trgdataformats::INVALID_TIMESTAMP)) {
    throw std::invalid_argument("Window and clock rate must be finite, positive and representable");
  }
  return static_cast<uint64_t>(ticks);
}

// One builder per detector/region, fed by an at-most-once, globally ordered TP merger.
// No ordering inference from arrival time: late data is explicitly rejected.
class OrderedInput
{
public:
  explicit OrderedInput(uint32_t sample_ticks) : m_sample_ticks(sample_ticks)
  {
    if (!sample_ticks) throw std::invalid_argument("Zero sample period");
  }
  void validate(const TP& tp) const
  {
    if (tp.version != TP::s_trigger_primitive_version || tp.channel >= trgdataformats::INVALID_TP_CHANNEL ||
        tp.detid == trgdataformats::INVALID_DETID || !tp.samples_over_threshold ||
        tp.samples_over_threshold > Frame::s_num_adcs || tp.samples_to_peak >= tp.samples_over_threshold ||
        tp.time_start < m_watermark || (m_last && tp.time_start < *m_last) ||
        (m_detector && tp.detid != *m_detector)) {
      throw std::invalid_argument("Invalid, late, unordered or mixed-detector PDS TP");
    }
    checked_add(tp.time_start, uint64_t(tp.samples_over_threshold) * m_sample_ticks);
  }
  void accept(const TP& tp) { m_last = uint64_t(tp.time_start); m_detector = uint8_t(tp.detid); }
  void watermark(uint64_t timestamp)
  {
    if (timestamp < m_watermark || (m_last && timestamp < *m_last)) {
      throw std::invalid_argument("Watermark must advance monotonically past processed data");
    }
    m_watermark = timestamp; // Future inputs must have time_start >= watermark.
  }

protected:
  uint32_t m_sample_ticks;
  uint64_t m_watermark = 0;
  std::optional<uint64_t> m_last;
  std::optional<uint8_t> m_detector;
};

struct CoincidenceConfig
{
  uint64_t window_ticks = 625; // 10 us at 62.5 MHz; inclusive [t-window, t].
  size_t minimum_channels = 6;
  size_t maximum_inputs = 100000;
  uint32_t ticks_per_sample = 1;
};

class CoincidenceBuilder : private OrderedInput
{
public:
  explicit CoincidenceBuilder(CoincidenceConfig config = {})
    : OrderedInput(config.ticks_per_sample), m_config(config)
  {
    if (!config.window_ticks || !config.minimum_channels ||
        config.maximum_inputs < config.minimum_channels) {
      throw std::invalid_argument("Invalid coincidence configuration");
    }
  }

  std::optional<triggeralgs::TriggerActivity> push(const TP& tp)
  {
    validate(tp);
    expire(tp.time_start);
    if (m_inputs.size() >= m_config.maximum_inputs) {
      throw std::length_error("PDS coincidence input capacity exceeded");
    }
    accept(tp);
    m_inputs.push_back(tp);
    std::set<uint32_t> channels;
    for (const auto& input : m_inputs) channels.insert(input.channel);
    if (channels.size() < m_config.minimum_channels) return std::nullopt;

    triggeralgs::TriggerActivity ta;
    ta.type = trgdataformats::TriggerActivityData::Type::kPDS;
    // This release has no PDS coincidence algorithm ID. Do not claim a TPC algorithm.
    ta.algorithm = trgdataformats::TriggerActivityData::Algorithm::kUnknown;
    ta.detid = tp.detid;
    ta.time_start = m_inputs.front().time_start;
    ta.time_end = 0;
    ta.time_activity = tp.time_start;
    ta.channel_start = *channels.begin();
    ta.channel_end = *channels.rbegin();
    bool first = true;
    for (const auto& input : m_inputs) {
      ta.inputs.push_back(input);
      ta.adc_integral += input.adc_integral;
      ta.time_end = std::max(ta.time_end, checked_add(input.time_start,
                              uint64_t(input.samples_over_threshold) * m_sample_ticks));
      if (first || input.adc_peak > ta.adc_peak) {
        ta.adc_peak = input.adc_peak;
        ta.channel_peak = input.channel;
        ta.time_peak = checked_add(input.time_start, uint64_t(input.samples_to_peak) * m_sample_ticks);
        first = false;
      }
    }
    // Consume contributing TPs. Persistent light can produce another TA using new TPs.
    m_inputs.clear();
    return ta;
  }

  void advance(uint64_t timestamp) { watermark(timestamp); expire(timestamp); }

private:
  void expire(uint64_t timestamp)
  {
    while (!m_inputs.empty() && timestamp - m_inputs.front().time_start > m_config.window_ticks) {
      m_inputs.pop_front();
    }
  }
  CoincidenceConfig m_config;
  std::deque<TP> m_inputs;
};

} // namespace dunedaq::fdreadoutlibs::pds
#endif
