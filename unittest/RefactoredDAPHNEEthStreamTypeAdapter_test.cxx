/**
 * @file RefactoredDAPHNEEthStreamTypeAdapter_test.cxx RefactoredDAPHNEEthStreamTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework.
 */

#include "fdreadoutlibs/RefactoredDAPHNEEthStreamTypeAdapter.hpp"

#define BOOST_TEST_MODULE RefactoredDAPHNEEthStreamTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(RefactoredDAPHNEEthStreamTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  RefactoredDAPHNEEthStreamTypeAdapter adapter;
  
  BOOST_REQUIRE(static_cast<size_t>(adapter.end() - adapter.begin()) == adapter.get_num_frames());

  BOOST_REQUIRE(static_cast<size_t>(reinterpret_cast<uint8_t*>(adapter.end()) -
                  reinterpret_cast<uint8_t*>(adapter.begin())) ==
                adapter.get_payload_size());
  BOOST_REQUIRE(sizeof(adapter) == adapter.get_payload_size());
  BOOST_REQUIRE(adapter.get_frame_size() == sizeof(dunedaq::fddetdataformats::DAPHNEEthStreamFrame));
}

BOOST_AUTO_TEST_CASE(TimestampsAndOrdering)
{
  RefactoredDAPHNEEthStreamTypeAdapter frame;
  RefactoredDAPHNEEthStreamTypeAdapter later_frame;

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
  RefactoredDAPHNEEthStreamTypeAdapter frame;

  frame.fake_geoid(7, 8, 9);

  BOOST_REQUIRE(frame.begin()->get_daqheader().crate_id == 7);
  BOOST_REQUIRE(frame.begin()->get_daqheader().slot_id == 8);
  BOOST_REQUIRE(frame.begin()->get_daqheader().stream_id == 9);
}

BOOST_AUTO_TEST_CASE(AdcPattern)
{
  RefactoredDAPHNEEthStreamTypeAdapter adapter{};

  adapter.fake_adc_pattern(2);

  bool any_nonzero = false;
  for (char byte : adapter.data) {
    if (byte != 0) {
      any_nonzero = true;
      break;
    }
  }

  BOOST_REQUIRE(any_nonzero);
}

BOOST_AUTO_TEST_CASE(StaticMetadata)
{
  BOOST_REQUIRE(RefactoredDAPHNEEthStreamTypeAdapter::subsystem == dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout);
  BOOST_REQUIRE(RefactoredDAPHNEEthStreamTypeAdapter::fragment_type == dunedaq::daqdataformats::FragmentType::kDAPHNEEthStream);
  BOOST_REQUIRE(RefactoredDAPHNEEthStreamTypeAdapter::expected_tick_difference == 280);
}

BOOST_AUTO_TEST_SUITE_END()
