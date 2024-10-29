/**
 * @file FDTypeAdapters_test.cxx  Unittest for expanding the WIBEth frames
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#define BOOST_TEST_MODULE FDTypeAdaptersBuffers_test // NOLINT

#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/DAPHNEStreamSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"

#include "datahandlinglibs/testutils/TestUtilities.hpp"

#include "datahandlinglibs/models/FixedRateQueueModel.hpp"
#include "datahandlinglibs/models/BinarySearchQueueModel.hpp"
#include "datahandlinglibs/models/SkipListLatencyBufferModel.hpp"


#include "boost/test/unit_test.hpp"

#include <iostream>
#include <sstream>
#include <set>
#include <iterator>

BOOST_AUTO_TEST_SUITE(FDReadoutTypeAdaptersBuffers_test)

BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DUNEWIBEth)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::FixedRateQueueModel,dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DUNEWIBEth)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::BinarySearchQueueModel,dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DAPHNEStreamSuperChunk)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::FixedRateQueueModel,dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DAPHNEStreamSuperChunk)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::BinarySearchQueueModel,dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(SkipListLatencyBufferModel_DAPHNESuperChunk)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::SkipListLatencyBufferModel,dunedaq::fdreadoutlibs::types::DAPHNESuperChunkTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(FixedRateQueueModel_TDEEth)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::FixedRateQueueModel,dunedaq::fdreadoutlibs::types::TDEEthTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_TDEEth)
{
    dunedaq::datahandlinglibs::test::test_queue_model<dunedaq::datahandlinglibs::BinarySearchQueueModel,dunedaq::fdreadoutlibs::types::TDEEthTypeAdapter>();
}

BOOST_AUTO_TEST_SUITE_END()


