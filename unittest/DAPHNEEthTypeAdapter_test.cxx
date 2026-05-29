/**
 * @file DAPHNEEthTypeAdapter_test.cxx DAPHNEEthTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/DAPHNEEthTypeAdapter.hpp"

#define BOOST_TEST_MODULE DAPHNEEthTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(DAPHNEEthTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  DAPHNEEthTypeAdapter frame;

  BOOST_REQUIRE(static_cast<size_t>(frame.end() - frame.begin()) == frame.get_num_frames());

  BOOST_REQUIRE(static_cast<size_t>(reinterpret_cast<uint8_t*>(frame.end()) -
                  reinterpret_cast<uint8_t*>(frame.begin())) ==
                frame.get_payload_size());

  BOOST_REQUIRE(frame.get_payload_size() == kDAPHNEEthSize);
  BOOST_REQUIRE(frame.get_frame_size() == kDAPHNEEthSize);
  BOOST_REQUIRE(frame.get_num_frames() == 1);
}

BOOST_AUTO_TEST_CASE(TimestampsAndOrdering)
{
  DAPHNEEthTypeAdapter frame;
  DAPHNEEthTypeAdapter same_timestamp_higher_channel;
  DAPHNEEthTypeAdapter later_frame;

  const uint64_t timestamp = 0x1234567800000000ULL;

  frame.set_timestamp(timestamp);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp);

  frame.fake_timestamps(timestamp + 10, 9999);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp + 10);

  frame.begin()->set_channel(3);
  same_timestamp_higher_channel.set_timestamp(timestamp + 10);
  same_timestamp_higher_channel.begin()->set_channel(4);
  later_frame.set_timestamp(timestamp + 11);
  later_frame.begin()->set_channel(0);

  BOOST_REQUIRE(frame < same_timestamp_higher_channel);
  BOOST_REQUIRE(same_timestamp_higher_channel < later_frame);
}

BOOST_AUTO_TEST_CASE(GeoID)
{
  DAPHNEEthTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->get_daqheader().crate_id == 7);
  BOOST_REQUIRE(frame.begin()->get_daqheader().slot_id == 8);
  BOOST_REQUIRE(frame.begin()->get_daqheader().stream_id == 9);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(DAPHNEEthTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(DAPHNEEthTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kDAPHNEEth);
  BOOST_REQUIRE(DAPHNEEthTypeAdapter::expected_tick_difference == 1);
}

BOOST_AUTO_TEST_SUITE_END()
