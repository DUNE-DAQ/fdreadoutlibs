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
