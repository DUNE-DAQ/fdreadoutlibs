#ifndef FDREADOUTLIBS_UNITTEST_PDSDESCRIPTORFIXTURES_HPP_
#define FDREADOUTLIBS_UNITTEST_PDSDESCRIPTORFIXTURES_HPP_
#include "fdreadoutlibs/pds/DescriptorToTP.hpp"

namespace pds_test {
using namespace dunedaq::fdreadoutlibs::pds;
inline Frame frame()
{
  Frame f{};
  f.header.version = 4;
  f.header.fragment_descriptor = 1;
  f.daq_header.det_id = 2;
  f.set_timestamp(1000);
  // Encode an independent wire descriptor, rather than filling the bitfield struct.
  f.header.set_descriptor_word(0, uint64_t(300) | (uint64_t(100) << 22) |
    (uint64_t(3) << 36) | (uint64_t(2) << 44) | (uint64_t(7) << 52) | (uint64_t(1) << 60));
  return f;
}
inline const ChannelMap mapping = [](const Frame& f) { return 1000 + f.get_channel() + 100 * f.daq_header.slot_id; };
inline TP hit(uint32_t channel, uint64_t timestamp, uint32_t integral = 100)
{
  TP tp;
  tp.channel = channel;
  tp.detid = 2;
  tp.time_start = timestamp;
  tp.samples_over_threshold = 4;
  tp.samples_to_peak = 1;
  tp.adc_integral = integral;
  tp.adc_peak = 50;
  return tp;
}
} // namespace pds_test

#endif
