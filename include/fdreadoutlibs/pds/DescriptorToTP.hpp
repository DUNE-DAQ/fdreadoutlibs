/** DAPHNE format-4 fragment-local peak descriptors to native DAQ TPs. */
#ifndef FDREADOUTLIBS_PDS_DESCRIPTORTOTP_HPP_
#define FDREADOUTLIBS_PDS_DESCRIPTORTOTP_HPP_

#include "fddetdataformats/DAPHNEEthFrame.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"

#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

namespace dunedaq::fdreadoutlibs::pds {

using Frame = fddetdataformats::DAPHNEEthFrame;
using TP = trgdataformats::TriggerPrimitive;
using ChannelMap = std::function<uint32_t(const Frame&)>;

struct DescriptorConfig
{
  uint32_t minimum_integral = 0;
  uint32_t ticks_per_sample = 1; // Current DAPHNE: 62.5 MHz sample and timestamp clocks.
  bool reject_overflow = true;
};

inline uint64_t checked_add(uint64_t timestamp, uint64_t delta)
{
  if (timestamp >= trgdataformats::INVALID_TIMESTAMP - delta) {
    throw std::overflow_error("PDS timestamp reaches invalid sentinel or wraps; start a new timing epoch");
  }
  return timestamp + delta;
}

inline std::vector<TP> descriptor_tps(const Frame& frame, const ChannelMap& map,
                                      DescriptorConfig config = {})
{
  if (!map || config.ticks_per_sample == 0) {
    throw std::invalid_argument("PDS conversion requires a channel map and nonzero sample period");
  }
  if (frame.header.version != Frame::version || !frame.header.fragment_descriptor) {
    throw std::invalid_argument("Expected format-4 fragment-local DAPHNE descriptors");
  }
  if (config.reject_overflow && frame.header.descriptor_overflow) {
    throw std::invalid_argument("DAPHNE descriptor overflow: incomplete light information");
  }
  const auto channel = map(frame);
  if (channel >= trgdataformats::INVALID_TP_CHANNEL) {
    throw std::invalid_argument("Missing or out-of-range PDS offline channel mapping");
  }
  // Validate the entire frame before returning any TPs. No partial frame emission.
  checked_add(frame.get_timestamp(), uint64_t(Frame::s_num_adcs) * config.ticks_per_sample);
  std::vector<TP> result;
  unsigned previous_end = 0;
  bool absent_seen = false;
  for (int i = 0; i < Frame::s_max_peaks; ++i) {
    const auto& peak = frame.header.peaks_data.peaks[i];
    if (!peak.found) {
      if (frame.header.get_descriptor_word(i) != 0) {
        throw std::invalid_argument("Nonzero absent DAPHNE descriptor");
      }
      absent_seen = true;
      continue;
    }
    const unsigned duration = peak.get_duration();
    if (absent_seen || peak.reserved || !peak.adc_peak || peak.sample_start < previous_end ||
        peak.sample_start + duration > Frame::s_num_adcs || peak.time_peak >= duration ||
        peak.adc_integral < peak.adc_peak || peak.adc_integral > duration * peak.adc_peak) {
      throw std::invalid_argument("Malformed DAPHNE peak descriptor");
    }
    previous_end = peak.sample_start + duration;
    if (peak.adc_integral < config.minimum_integral) {
      continue;
    }
    TP tp;
    tp.detid = frame.daq_header.det_id;
    tp.channel = channel;
    tp.time_start = checked_add(frame.get_timestamp(), uint64_t(peak.sample_start) * config.ticks_per_sample);
    tp.samples_over_threshold = duration;
    tp.samples_to_peak = peak.time_peak; // Offset from excursion start, NOT from frame start.
    tp.adc_integral = peak.adc_integral; // Already baseline subtracted by firmware.
    tp.adc_peak = peak.adc_peak;
    result.push_back(tp);
  }
  return result;
}

} // namespace dunedaq::fdreadoutlibs::pds
#endif
