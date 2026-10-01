#define BOOST_TEST_MODULE PDSDescriptorProcessor_test
#include "boost/test/unit_test.hpp"
#include "fdreadoutlibs/pds/DescriptorProcessor.hpp"

using namespace dunedaq::fdreadoutlibs::pds;

DescriptorFrame descriptor(unsigned channel = 2)
{
  Frame frame{};
  frame.daq_header.timestamp = 1000;
  frame.header.version = Frame::version;
  frame.header.fragment_descriptor = 1;
  frame.header.channel = channel;
  auto& peak = frame.header.peaks_data.peaks[0];
  peak.found = 1; peak.sample_start = 3; peak.duration_minus_one = 1;
  peak.adc_peak = 10; peak.adc_integral = 15; peak.time_peak = 1;
  return {frame.daq_header, frame.header};
}

BOOST_AUTO_TEST_CASE(converts_metadata_without_waveforms)
{
  std::vector<TP> got;
  DescriptorProcessor processor([](const Frame& f) { return f.get_channel() + 100; }, {}, 0,
    [&](std::vector<TP>&& tps) { got = std::move(tps); return true; });
  processor.process(descriptor());
  BOOST_REQUIRE_EQUAL(got.size(), 1);
  BOOST_CHECK_EQUAL(got[0].channel, 102);
  BOOST_CHECK_EQUAL(got[0].time_start, 1003);
  BOOST_CHECK_EQUAL(got[0].adc_integral, 15);
  BOOST_CHECK_EQUAL(processor.counters.sent.load(), 1);
}

BOOST_AUTO_TEST_CASE(counts_rejection_and_sink_failure)
{
  DescriptorProcessor processor([](const Frame& f) { return f.get_channel(); }, {}, 0,
    [](std::vector<TP>&&) { return false; });
  auto input = descriptor();
  processor.process(input);
  input.header.descriptor_overflow = 1;
  processor.process(input);
  input.header.descriptor_overflow = 0;
  input.header.peaks_data.peaks[0].reserved = 1;
  processor.process(input);
  BOOST_CHECK_EQUAL(processor.counters.frames.load(), 3);
  BOOST_CHECK_EQUAL(processor.counters.send_failed.load(), 1);
  BOOST_CHECK_EQUAL(processor.counters.overflow.load(), 1);
  BOOST_CHECK_EQUAL(processor.counters.malformed.load(), 1);
}

BOOST_AUTO_TEST_CASE(applies_mask_and_cut)
{
  unsigned calls = 0;
  DescriptorProcessor processor([](const Frame& f) { return f.get_channel(); }, {2}, 16,
    [&](std::vector<TP>&&) { ++calls; return true; });
  processor.process(descriptor(2));
  processor.process(descriptor(3));
  BOOST_CHECK_EQUAL(calls, 0);
}

BOOST_AUTO_TEST_CASE(registry_ownership_and_removal)
{
  auto processor = std::make_shared<DescriptorProcessor>([](const Frame&) { return 1; },
    std::set<unsigned>{}, 0, [](std::vector<TP>&&) { return true; });
  register_descriptor_processor("test", processor);
  BOOST_CHECK(get_descriptor_processor("test") == processor);
  BOOST_CHECK_THROW(register_descriptor_processor("test", processor), std::logic_error);
  remove_descriptor_processor("test", nullptr);
  BOOST_CHECK(get_descriptor_processor("test") == processor);
  remove_descriptor_processor("test", processor);
  BOOST_CHECK(!get_descriptor_processor("test"));
  register_descriptor_processor("test", processor);
  processor.reset();
  BOOST_CHECK(!get_descriptor_processor("test"));
}
