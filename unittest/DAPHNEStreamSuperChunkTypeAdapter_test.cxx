/**
 * @file DAPHNEStreamSuperChunkTypeAdapter_test.cxx DAPHNEStreamSuperChunkTypeAdapter class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/DAPHNEStreamSuperChunkTypeAdapter.hpp"

#define BOOST_TEST_MODULE DAPHNEStreamSuperChunkTypeAdapter_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(DAPHNEStreamSuperChunkTypeAdapter_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  DAPHNEStreamSuperChunkTypeAdapter blob_of_junk;
  
  BOOST_REQUIRE(blob_of_junk.end() - blob_of_junk.begin() \
		== blob_of_junk.get_num_frames());

  BOOST_REQUIRE(reinterpret_cast<uint8_t*>(blob_of_junk.end()) - \
		reinterpret_cast<uint8_t*>(blob_of_junk.begin()) \
		== blob_of_junk.get_payload_size());
}

BOOST_AUTO_TEST_CASE(Timestamps)
{
  DAPHNEStreamSuperChunkTypeAdapter blob_of_junk;
  DAPHNEStreamSuperChunkTypeAdapter blob_of_junk2;

  const uint64_t timestamp = 0xDEADBEEFA0B0C0D0;

  blob_of_junk.set_timestamp(timestamp);
  BOOST_REQUIRE(blob_of_junk.get_timestamp() == timestamp);

  blob_of_junk2.set_timestamp(timestamp + 1);

  BOOST_REQUIRE(blob_of_junk < blob_of_junk2);

  const uint64_t offset = 32;
  blob_of_junk.fake_timestamps(timestamp, offset);

  for (size_t i = 0; i < blob_of_junk.get_num_frames(); ++i) {
    const DAPHNEStreamSuperChunkTypeAdapter::FrameType* frame = blob_of_junk.begin() + i;
    BOOST_REQUIRE(frame->daq_header.get_timestamp() == timestamp + i*offset);
  }
  
}

BOOST_AUTO_TEST_SUITE_END()
