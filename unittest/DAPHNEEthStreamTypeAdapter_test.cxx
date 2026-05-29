/**
 * @file DAPHNEEthStreamTypeAdapter_test.cxx DAPHNEEthStreamTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/DAPHNEEthStreamTypeAdapter.hpp"

#define BOOST_TEST_MODULE DAPHNEEthStreamTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(DAPHNEEthStreamTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  DAPHNEEthStreamTypeAdapter frame;

  BOOST_REQUIRE(static_cast<size_t>(frame.end() - frame.begin()) == frame.get_num_frames());

  BOOST_REQUIRE(static_cast<size_t>(reinterpret_cast<uint8_t*>(frame.end()) -
                  reinterpret_cast<uint8_t*>(frame.begin())) ==
                frame.get_payload_size());

  BOOST_REQUIRE(frame.get_payload_size() == kDAPHNEEthStreamSize);
  BOOST_REQUIRE(frame.get_frame_size() == sizeof(DAPHNEEthStreamTypeAdapter::FrameType));
  BOOST_REQUIRE(frame.get_num_frames() == kDAPHNEEthStreamNumFrames);
}

BOOST_AUTO_TEST_CASE(TimestampsAndOrdering)
{
  DAPHNEEthStreamTypeAdapter frame;
  DAPHNEEthStreamTypeAdapter later_frame;

  const uint64_t timestamp = 0x1234567800000000ULL;

  frame.set_timestamp(timestamp);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp);

  frame.fake_timestamps(timestamp + 10, 280);
  BOOST_REQUIRE(frame.get_timestamp() == timestamp + 10);

  later_frame.set_timestamp(timestamp + 11);
  BOOST_REQUIRE(frame < later_frame);
}

BOOST_AUTO_TEST_CASE(GeoID)
{
  DAPHNEEthStreamTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->get_daqheader().crate_id == 7);
  BOOST_REQUIRE(frame.begin()->get_daqheader().slot_id == 8);
  BOOST_REQUIRE(frame.begin()->get_daqheader().stream_id == 9);
}

BOOST_AUTO_TEST_CASE(AdcPattern)
{
  DAPHNEEthStreamTypeAdapter frame{};

  frame.fake_adc_pattern(2);

  bool any_nonzero = false;
  for (char byte : frame.data) {
    if (byte != 0) {
      any_nonzero = true;
      break;
    }
  }

  BOOST_REQUIRE(any_nonzero);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(DAPHNEEthStreamTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(DAPHNEEthStreamTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kDAPHNEEthStream);
  BOOST_REQUIRE(DAPHNEEthStreamTypeAdapter::expected_tick_difference == 280);
}

BOOST_AUTO_TEST_SUITE_END()
