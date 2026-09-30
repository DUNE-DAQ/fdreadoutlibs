#define BOOST_TEST_MODULE PDSDescriptorToTP_test
#include "boost/test/unit_test.hpp"
#include "PDSDescriptorFixtures.hpp"

using namespace pds_test;

BOOST_AUTO_TEST_CASE(wire_descriptor_maps_to_native_tp)
{
  auto tps = descriptor_tps(frame(), mapping);
  BOOST_REQUIRE_EQUAL(tps.size(), 1);
  const auto& tp = tps.front();
  BOOST_CHECK_EQUAL(+tp.time_start, 1007);
  BOOST_CHECK_EQUAL(+tp.samples_to_peak, 2);
  BOOST_CHECK_EQUAL(+tp.samples_over_threshold, 4);
  BOOST_CHECK_EQUAL(+tp.adc_integral, 300);
  BOOST_CHECK_EQUAL(+tp.adc_peak, 100);
  BOOST_CHECK_EQUAL(+tp.channel, 1000);
  BOOST_CHECK_EQUAL(+tp.detid, 2);
  BOOST_CHECK_EQUAL(+tp.version, TP::s_trigger_primitive_version);
  BOOST_CHECK_EQUAL(+descriptor_tps(frame(), mapping, {0, 2, true})[0].time_start, 1014);
}

BOOST_AUTO_TEST_CASE(empty_threshold_and_channel_mapping)
{
  auto f = frame();
  BOOST_CHECK(descriptor_tps(f, mapping, {301, 1, true}).empty());
  BOOST_CHECK_EQUAL(descriptor_tps(f, mapping, {300, 1, true}).size(), 1);
  f.daq_header.slot_id = 2;
  BOOST_CHECK_EQUAL(+descriptor_tps(f, mapping)[0].channel, 1200);
  f.header.set_descriptor_word(0, 0);
  BOOST_CHECK(descriptor_tps(f, mapping).empty());
  BOOST_CHECK_THROW(descriptor_tps(f, {}), std::invalid_argument);
  BOOST_CHECK_THROW(descriptor_tps(f, [](const Frame&) { return 0xFFFFFF; }), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(reject_bad_format_overflow_and_malformed_fields)
{
  auto f = frame();
  f.header.version = 3;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.fragment_descriptor = 0;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.descriptor_overflow = 1;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  BOOST_CHECK_EQUAL(descriptor_tps(f, mapping, {0, 1, false}).size(), 1);
  f = frame(); f.header.peaks_data.peaks[0].time_peak = 4;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.peaks_data.peaks[0].sample_start = 254;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.peaks_data.peaks[0].reserved = 1;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.set_descriptor_word(1, 1);
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.peaks_data.peaks[0].adc_integral = 401;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f = frame(); f.header.peaks_data.peaks[0].adc_peak = 0; f.header.peaks_data.peaks[0].adc_integral = 0;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(full_duration_and_continuation)
{
  auto f = frame();
  f.header.continuation = 1;
  auto& p = f.header.peaks_data.peaks[0];
  p.sample_start = 0; p.duration_minus_one = 255; p.time_peak = 255;
  p.adc_peak = 16383; p.adc_integral = 4194048;
  const auto tps = descriptor_tps(f, mapping);
  BOOST_CHECK_EQUAL(+tps[0].samples_over_threshold, 256);
  BOOST_CHECK_EQUAL(+tps[0].samples_to_peak, 255);
  BOOST_CHECK_EQUAL(+tps[0].adc_integral, 4194048);
}

BOOST_AUTO_TEST_CASE(five_descriptors_and_atomic_frame_validation)
{
  auto f = frame();
  const auto word = f.header.get_descriptor_word(0);
  for (int i = 0; i < 5; ++i) {
    f.header.set_descriptor_word(i, word);
    f.header.peaks_data.peaks[i].sample_start = 7 + 10 * i;
  }
  BOOST_CHECK_EQUAL(descriptor_tps(f, mapping).size(), 5);
  f.header.peaks_data.peaks[4].sample_start = 38;
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
  f.header.set_descriptor_word(0, 0);
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(timestamp_wrap_is_rejected)
{
  auto f = frame(); f.set_timestamp(UINT64_MAX - 100);
  BOOST_CHECK_THROW(descriptor_tps(f, mapping), std::overflow_error);
  BOOST_CHECK_THROW(checked_add(UINT64_MAX - 1, 1), std::overflow_error);
}
