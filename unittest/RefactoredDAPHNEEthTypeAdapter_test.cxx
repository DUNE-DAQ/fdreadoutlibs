/**
 * @file RefactoredDAPHNEEthTypeAdapter_test.cxx RefactoredDAPHNEEthTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/RefactoredDAPHNEEthTypeAdapter.hpp"
#include "fdreadoutlibs/DAPHNEEthTypeAdapter.hpp"

#define BOOST_TEST_MODULE RefactoredDAPHNEEthTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(RefactoredDAPHNEEthTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  RefactoredDAPHNEEthTypeAdapter frame;

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
  RefactoredDAPHNEEthTypeAdapter frame;
  RefactoredDAPHNEEthTypeAdapter same_timestamp_higher_channel;
  RefactoredDAPHNEEthTypeAdapter later_frame;

  const uint64_t timestamp = 0x1234567800000000ULL;
  std::cout << "timestamp is " << timestamp << "\n";
  
  frame.set_timestamp(timestamp);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp);

  frame.fake_timestamps(timestamp + 10, 9999);
  std::cout << "new timestamp is " << frame.get_timestamp() << "\n";
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
  RefactoredDAPHNEEthTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->get_daqheader().crate_id == 7);
  BOOST_REQUIRE(frame.begin()->get_daqheader().slot_id == 8);
  BOOST_REQUIRE(frame.begin()->get_daqheader().stream_id == 9);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(RefactoredDAPHNEEthTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(RefactoredDAPHNEEthTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kDAPHNEEth);
  BOOST_REQUIRE(RefactoredDAPHNEEthTypeAdapter::expected_tick_difference == 1);
}

BOOST_AUTO_TEST_SUITE_END()
