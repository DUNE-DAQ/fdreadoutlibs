#ifndef FDREADOUTLIBS_PDS_DESCRIPTORPROCESSOR_HPP_
#define FDREADOUTLIBS_PDS_DESCRIPTORPROCESSOR_HPP_

#include "fdreadoutlibs/pds/DescriptorToTP.hpp"
#include <atomic>
#include <algorithm>
#include <new>
#include <cstddef>
#include <memory>
#include <set>
#include <string>

namespace dunedaq::fdreadoutlibs::pds {

struct DescriptorFrame
{
  decltype(Frame::daq_header) daq_header;
  Frame::Header header;
};
static_assert(sizeof(DescriptorFrame) == 64);
static_assert(offsetof(Frame, adc_words) == sizeof(DescriptorFrame));

class DescriptorProcessor
{
public:
  using Emit = std::function<bool(std::vector<TP>&&)>;
  struct Counters {
    std::atomic<uint64_t> frames{0}, overflow{0}, malformed{0}, sent{0}, send_failed{0};
  } counters;

  DescriptorProcessor(ChannelMap map, std::set<unsigned int> mask, uint32_t threshold, Emit emit)
    : m_map(std::move(map)), m_mask(std::move(mask)), m_threshold(threshold), m_emit(std::move(emit))
  {
    if (!m_map || !m_emit) throw std::invalid_argument("Descriptor processor requires map and sink");
  }

  void process(const DescriptorFrame& descriptor)
  {
    ++counters.frames;
    Frame frame{};
    frame.daq_header = descriptor.daq_header;
    frame.header = descriptor.header;
    if (frame.header.version == Frame::version && frame.header.fragment_descriptor &&
        frame.header.descriptor_overflow) {
      ++counters.overflow;
      return;
    }
    std::vector<TP> output;
    try {
      output = descriptor_tps(frame, m_map, DescriptorConfig{m_threshold, 1, true});
    } catch (const std::bad_alloc&) { throw; }
    catch (const std::exception&) { ++counters.malformed; return; }
    output.erase(std::remove_if(output.begin(), output.end(), [this](const TP& tp) {
      return m_mask.count(tp.channel);
    }), output.end());
    if (output.empty()) return;
    auto count = output.size();
    if (m_emit(std::move(output))) counters.sent += count;
    else counters.send_failed += count;
  }

private:
  const ChannelMap m_map;
  const std::set<unsigned int> m_mask;
  const uint32_t m_threshold;
  const Emit m_emit;
};

// Defined in the shared library so the reader and raw-processor plugins share one registry.
void register_descriptor_processor(const std::string& key, const std::shared_ptr<DescriptorProcessor>& processor);
std::shared_ptr<DescriptorProcessor> get_descriptor_processor(const std::string& key);
void remove_descriptor_processor(const std::string& key, const std::shared_ptr<DescriptorProcessor>& processor);

} // namespace dunedaq::fdreadoutlibs::pds
#endif
