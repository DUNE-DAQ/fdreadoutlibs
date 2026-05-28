/**
 * @file TDEEthTypeAdapter_test.cxx TDEEthTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"

#define BOOST_TEST_MODULE TDEEthTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(TDEEthTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  TDEEthTypeAdapter frame;

  BOOST_REQUIRE(static_cast<size_t>(frame.end() - frame.begin()) == frame.get_num_frames());

  BOOST_REQUIRE(static_cast<size_t>(reinterpret_cast<uint8_t*>(frame.end()) -
                  reinterpret_cast<uint8_t*>(frame.begin())) ==
                frame.get_payload_size());

  BOOST_REQUIRE(frame.get_payload_size() == TDEEthTypeAdapter::fixed_payload_size);
  BOOST_REQUIRE(frame.get_frame_size() == TDEEthTypeAdapter::fixed_payload_size);
  BOOST_REQUIRE(frame.get_num_frames() == 1);
}

BOOST_AUTO_TEST_CASE(TimestampsAndOrdering)
{
  TDEEthTypeAdapter frame;
  TDEEthTypeAdapter next_frame;

  const uint64_t timestamp = 0x1234567800000000ULL;

  frame.set_timestamp(timestamp);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp);

  frame.fake_timestamps(timestamp + 10, 9999);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp + 10);

  next_frame.set_timestamp(timestamp + 11);
  BOOST_REQUIRE(frame < next_frame);
}

BOOST_AUTO_TEST_CASE(GeoID)
{
  TDEEthTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->daq_header.crate_id == 7);
  BOOST_REQUIRE(frame.begin()->daq_header.slot_id == 8);
  BOOST_REQUIRE(frame.begin()->daq_header.stream_id == 9);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(TDEEthTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(TDEEthTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kTDEEth);
  BOOST_REQUIRE(TDEEthTypeAdapter::expected_tick_difference == 2000);
  BOOST_REQUIRE(TDEEthTypeAdapter::samples_per_frame == 64);
  BOOST_REQUIRE(TDEEthTypeAdapter::samples_tick_difference == 31.25f);
}

BOOST_AUTO_TEST_SUITE_END()
