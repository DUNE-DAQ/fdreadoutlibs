#include "fdreadoutlibs/pds/ActivityBuilder.hpp"

#define BOOST_TEST_MODULE PDSDescriptorPipeline_test
#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::pds;

#include "PDSDescriptorFixtures.hpp"
using namespace pds_test;

BOOST_AUTO_TEST_CASE(six_distinct_channels_not_six_pulses)
{
  CoincidenceBuilder builder({10, 6, 100, 1});
  for (unsigned i = 0; i < 10; ++i) BOOST_CHECK(!builder.push(hit(0, 100 + i)));
  for (unsigned ch = 1; ch < 5; ++ch) BOOST_CHECK(!builder.push(hit(ch, 110)));
  auto ta = builder.push(hit(5, 110));
  BOOST_REQUIRE(ta);
  BOOST_CHECK(ta->type == dunedaq::trgdataformats::TriggerActivityData::Type::kPDS);
  BOOST_CHECK_EQUAL(ta->inputs.size(), 15);
  BOOST_CHECK_EQUAL(ta->adc_integral, 1500);
  BOOST_CHECK_EQUAL(ta->time_start, 100);
  BOOST_CHECK_EQUAL(ta->time_end, 114);
  BOOST_CHECK_EQUAL(ta->time_activity, 110);
  BOOST_CHECK_EQUAL(ta->channel_end, 5);
  BOOST_CHECK(!builder.push(hit(6, 111))); // Contributing inputs were consumed.
}

BOOST_AUTO_TEST_CASE(inclusive_coincidence_boundary_and_expiry)
{
  CoincidenceBuilder inclusive({10, 2, 100, 1});
  BOOST_CHECK(!inclusive.push(hit(0, 100)));
  BOOST_CHECK(inclusive.push(hit(1, 110)));
  CoincidenceBuilder outside({10, 2, 100, 1});
  BOOST_CHECK(!outside.push(hit(0, 100)));
  BOOST_CHECK(!outside.push(hit(1, 111)));
  outside.advance(200);
  BOOST_CHECK(!outside.push(hit(2, 200)));
}

BOOST_AUTO_TEST_CASE(late_data_mixed_detectors_and_capacity)
{
  CoincidenceBuilder builder({10, 2, 2, 1});
  builder.push(hit(0, 100));
  BOOST_CHECK_THROW(builder.push(hit(1, 99)), std::invalid_argument);
  auto other = hit(1, 100); other.detid = 3;
  BOOST_CHECK_THROW(builder.push(other), std::invalid_argument);
  builder.push(hit(0, 101));
  BOOST_CHECK_THROW(builder.push(hit(0, 102)), std::length_error);
  builder.advance(200);
  BOOST_CHECK_THROW(builder.push(hit(1, 199)), std::invalid_argument);
  BOOST_CHECK_THROW(builder.advance(199), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(prompt_half_open_gates_and_quiet_watermark)
{
  PromptLightBuilder builder({5, 20, 1});
  BOOST_CHECK(!builder.push(hit(0, 100, 100)));
  BOOST_CHECK(!builder.push(hit(1, 104, 200)));
  BOOST_CHECK(!builder.push(hit(2, 105, 300))); // Excluded from prompt, included in total.
  BOOST_CHECK(!builder.advance(119));
  auto light = builder.advance(120);
  BOOST_REQUIRE(light);
  BOOST_CHECK_EQUAL(light->prompt_integral, 300);
  BOOST_CHECK_EQUAL(light->total_integral, 600);
  BOOST_CHECK_CLOSE(*light->fraction(), 0.5, 1.e-9);
  BOOST_CHECK_EQUAL(light->pulse_count, 3);
  BOOST_CHECK(!builder.advance(130)); // No duplicate flush.
  BOOST_CHECK_THROW(builder.push(hit(0, 129)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(prompt_boundary_starts_next_gate_and_zero_charge)
{
  PromptLightBuilder builder({5, 20, 1});
  builder.push(hit(0, 100, 0));
  auto first = builder.push(hit(1, 120, 100));
  BOOST_REQUIRE(first);
  BOOST_CHECK(!first->fraction());
  auto second = builder.advance(140);
  BOOST_REQUIRE(second);
  BOOST_CHECK_EQUAL(second->total_integral, 100);
  BOOST_CHECK_EQUAL(second->time_start, 120);
}

BOOST_AUTO_TEST_CASE(window_units_and_configuration_validation)
{
  BOOST_CHECK_EQUAL(microseconds_to_ticks(10), 625);
  BOOST_CHECK_EQUAL(microseconds_to_ticks(0.1), 7);
  BOOST_CHECK_THROW(microseconds_to_ticks(0), std::invalid_argument);
  BOOST_CHECK_THROW(microseconds_to_ticks(-1), std::invalid_argument);
  BOOST_CHECK_THROW(microseconds_to_ticks(NAN), std::invalid_argument);
  BOOST_CHECK_THROW(microseconds_to_ticks(INFINITY), std::invalid_argument);
  BOOST_CHECK_THROW(CoincidenceBuilder(CoincidenceConfig{0, 6, 100, 1}), std::invalid_argument);
  BOOST_CHECK_THROW(PromptLightBuilder(PromptConfig{21, 20, 1}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(prompt_existing_gate_near_timestamp_limit)
{
  PromptLightBuilder builder({5, 60, 1});
  builder.push(hit(0, UINT64_MAX - 100));
  BOOST_CHECK_NO_THROW(builder.push(hit(1, UINT64_MAX - 50)));
  auto light = builder.finish();
  BOOST_REQUIRE(light);
  BOOST_CHECK_EQUAL(light->time_end, UINT64_MAX - 40);
  BOOST_CHECK_EQUAL(light->total_integral, 200);
  BOOST_CHECK(!builder.finish());
  BOOST_CHECK_THROW(builder.push(hit(2, UINT64_MAX - 41)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(prompt_rejected_new_gate_preserves_pending_light)
{
  PromptLightBuilder builder({5, 60, 1});
  builder.push(hit(0, UINT64_MAX - 100));
  BOOST_CHECK_THROW(builder.push(hit(1, UINT64_MAX - 30)), std::overflow_error);
  auto light = builder.finish();
  BOOST_REQUIRE(light);
  BOOST_CHECK_EQUAL(light->total_integral, 100);
}

BOOST_AUTO_TEST_CASE(end_to_end_descriptors_to_activity_and_light)
{
  CoincidenceBuilder coincidence;
  PromptLightBuilder prompt({10, 100, 1});
  unsigned count = 0;
  for (unsigned ch = 0; ch < 6; ++ch) {
    auto f = frame(); f.set_channel(ch);
    f.header.peaks_data.peaks[0].sample_start = 7 + ch;
    auto tp = descriptor_tps(f, mapping)[0];
    if (auto ta = coincidence.push(tp)) {
      ++count;
      BOOST_CHECK_EQUAL(ta->inputs.size(), 6);
      BOOST_CHECK_EQUAL(ta->adc_integral, 1800);
    }
    BOOST_CHECK(!prompt.push(tp));
  }
  BOOST_CHECK_EQUAL(count, 1);
  auto light = prompt.advance(1107);
  BOOST_REQUIRE(light);
  BOOST_CHECK_EQUAL(light->prompt_integral, 1800);
}
