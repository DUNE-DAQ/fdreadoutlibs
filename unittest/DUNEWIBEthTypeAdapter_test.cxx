/**
 * @file DUNEWIBEthTypeAdapter_test.cxx DUNEWIBEthTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"

#define BOOST_TEST_MODULE DUNEWIBEthTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(DUNEWIBEthTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  DUNEWIBEthTypeAdapter frame;

  BOOST_REQUIRE(static_cast<size_t>(frame.end() - frame.begin()) == frame.get_num_frames());

  BOOST_REQUIRE(static_cast<size_t>(reinterpret_cast<uint8_t*>(frame.end()) -
                  reinterpret_cast<uint8_t*>(frame.begin())) ==
                frame.get_payload_size());

  BOOST_REQUIRE(frame.get_payload_size() == kDUNEWIBEthSize);
  BOOST_REQUIRE(frame.get_frame_size() == kDUNEWIBEthSize);
  BOOST_REQUIRE(frame.get_num_frames() == 1);
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::fixed_payload_size == kDUNEWIBEthSize);
}

BOOST_AUTO_TEST_CASE(TimestampsAndOrdering)
{
  DUNEWIBEthTypeAdapter frame;
  DUNEWIBEthTypeAdapter later_frame;

  const uint64_t timestamp = 0x1234567800000000ULL;

  frame.set_timestamp(timestamp);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp);

  frame.fake_timestamps(timestamp + 10, 9999);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp + 10);

  later_frame.set_timestamp(timestamp + 11);

  BOOST_REQUIRE(frame < later_frame);
}

BOOST_AUTO_TEST_CASE(GeoID)
{
  DUNEWIBEthTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->daq_header.crate_id == 7);
  BOOST_REQUIRE(frame.begin()->daq_header.slot_id == 8);
  BOOST_REQUIRE(frame.begin()->daq_header.stream_id == 9);
}

BOOST_AUTO_TEST_CASE(AdcPattern)
{
  DUNEWIBEthTypeAdapter frame{};

  frame.fake_adc_pattern(2);

  BOOST_REQUIRE(frame.begin()->get_adc(2, 0) == 0x3FFF);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kWIBEth);
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::expected_tick_difference == 2048);
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::samples_per_frame == 64);
  BOOST_REQUIRE(DUNEWIBEthTypeAdapter::samples_tick_difference == 32.0f);
}

BOOST_AUTO_TEST_SUITE_END()