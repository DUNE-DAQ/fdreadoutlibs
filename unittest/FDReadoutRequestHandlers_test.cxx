/**
 * @file FDReadoutRequestHandlers_test.cxx  Unittest for expanding the WIBEth frames
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#define BOOST_TEST_MODULE FDReadoutRequestHandlers_test // NOLINT

#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/DAPHNEStreamSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"

#include "datahandlinglibs/models/FixedRateQueueModel.hpp"
#include "datahandlinglibs/models/BinarySearchQueueModel.hpp"
#include "datahandlinglibs/models/SkipListLatencyBufferModel.hpp"

#include "datahandlinglibs/models/DefaultRequestHandlerModel.hpp"
#include "datahandlinglibs/FrameErrorRegistry.hpp"

#include "boost/test/unit_test.hpp"

#include <iostream>
#include <sstream>
#include <set>
#include <iterator>

//general test buffer struct that can be reused
template <template<class> class BufferType, class TypeAdapter>
void fill_buffer( std::shared_ptr< BufferType<TypeAdapter> >& buffer,
		  uint64_t const init_timestamp=0,
		  size_t const n_obj=10,
		  std::set<size_t> const obj_to_skip = {})
{
  //allocate mem in buffer
  buffer->allocate_memory(n_obj);
  
  //some accounting
  size_t i=0,i_obj=0;
  uint64_t next_timestamp=init_timestamp;
  
  //while we want to add to the buffer
  while(i_obj<n_obj){
    
    //create a frame with appropriate timestamp
    TypeAdapter frame;
    frame.fake_timestamps(next_timestamp);
    
    //increment to the next timestamp
    next_timestamp += uint64_t(frame.get_num_frames()*TypeAdapter::expected_tick_difference);
    
    //only write in buffer if not in the skip list
    if(obj_to_skip.count(i)==0) {
      buffer->write(std::move(frame));
      ++i_obj;
    }
    ++i;
    
    
  } // end while(obj_in_buffer<n_obj)
}

template <template<class> class BufferType, class TypeAdapter>
void print_buffer(std::shared_ptr<BufferType<TypeAdapter>>& buffer, std::string desc)
{
  std::stringstream ss;
  ss << "Buffer (" << desc << "): ";
  typename BufferType<TypeAdapter>::Iterator iter=buffer->begin();
  while(iter!=buffer->end()){
    ss << iter->get_timestamp() << " ";
    ++iter;
  }
  BOOST_TEST_MESSAGE(ss.str());

}

template<class ReadoutType, class LatencyBufferType>
class TestDefaultRequestHandlerModel : public dunedaq::datahandlinglibs::DefaultRequestHandlerModel<ReadoutType, LatencyBufferType>
{
public:
    TestDefaultRequestHandlerModel(std::shared_ptr<LatencyBufferType>& latency_buffer,
                                   std::unique_ptr<dunedaq::datahandlinglibs::FrameErrorRegistry>& error_registry)
          : dunedaq::datahandlinglibs::DefaultRequestHandlerModel<ReadoutType, LatencyBufferType>(latency_buffer,error_registry) {}
    using dunedaq::datahandlinglibs::DefaultRequestHandlerModel<ReadoutType, LatencyBufferType>::get_fragment_pieces;
};


template <template<class> class BufferType, class TypeAdapter>
void test_request_model()
{
    using DefaultRequestHandler = TestDefaultRequestHandlerModel<TypeAdapter, BufferType<TypeAdapter> >;

  //some testing vars
  TypeAdapter test_element;
  uint64_t ticks_between = TypeAdapter::expected_tick_difference*test_element.get_num_frames();


  /*
   * Unskipped buffer should have elements with index [0, 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ]
   *                                   and timestamps [0,1*T,2*T,3*T,4*T,5*T,6*T,7*T,8*T,9*T]
   * where T = DTS ticks between successive elements (tick_diff_per_frame * n_frames_per_obj_in_buffer)
   */
  BOOST_TEST_MESSAGE("Testing buffer without skips...");
  auto buffer_noskip = std::make_shared< BufferType<TypeAdapter> >();
  fill_buffer<BufferType,TypeAdapter>(buffer_noskip,0,10);
  print_buffer<BufferType,TypeAdapter>(buffer_noskip,"noskip");

  //make the error registry we need
  auto errorRegistry_noskip = std::make_unique<dunedaq::datahandlinglibs::FrameErrorRegistry>();

  //create the request handler
  DefaultRequestHandler requestHandler(buffer_noskip,errorRegistry_noskip);

  //
  auto test_req_bounds = [&](uint64_t start_win, uint64_t end_win,
                             uint64_t expected_start, uint64_t expected_end){

      auto dfmessage = dunedaq::dfmessages::DataRequest();
      auto req_res = typename DefaultRequestHandler::RequestResult(DefaultRequestHandler::ResultCode::kUnknown,dfmessage);
      auto ret = requestHandler.get_fragment_pieces(start_win,end_win,req_res);

      //check that the return code is correct.
      BOOST_CHECK_EQUAL(req_res.result_code,DefaultRequestHandler::ResultCode::kFound);

      //grab the (first) timestamps of the first and last frames
      uint64_t first_ts = reinterpret_cast<const TypeAdapter::FrameType*>(ret.front().first)->get_timestamp();
      uint64_t last_ts = reinterpret_cast<const TypeAdapter::FrameType*>((char*)(ret.back().first)+(ret.back().second)-sizeof(typename TypeAdapter::FrameType))->get_timestamp();

      //general check:
      // first_ts <= start_win < first_ts + ticks_per_frame
      // last_ts < end_win <= last_ts + ticks_per_frame
      BOOST_CHECK_LE(first_ts,start_win);
      BOOST_CHECK_GT(first_ts+TypeAdapter::expected_tick_difference,start_win);

      BOOST_CHECK_GE(last_ts+TypeAdapter::expected_tick_difference,end_win);
      BOOST_CHECK_LT(last_ts,end_win);

      //specfic check: are values for this request what we expect
      BOOST_CHECK_EQUAL(first_ts,expected_start);
      BOOST_CHECK_EQUAL(last_ts+TypeAdapter::expected_tick_difference,expected_end);
  };
  test_req_bounds(ticks_between*2,ticks_between*5,ticks_between*2,ticks_between*5);
  test_req_bounds(ticks_between*3/2,ticks_between*9/2,ticks_between,ticks_between*5);

}

BOOST_AUTO_TEST_SUITE(FDReadoutRequestHandlers_test)

BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DUNEWIBEth)
{
  test_request_model<
          dunedaq::datahandlinglibs::FixedRateQueueModel,
          dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
/*
BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DUNEWIBEth)
{
  test_queue_model<dunedaq::datahandlinglibs::BinarySearchQueueModel,dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DAPHNEStreamSuperChunk)
{
  test_queue_model<dunedaq::datahandlinglibs::FixedRateQueueModel,dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DAPHNEStreamSuperChunk)
{
  test_queue_model<dunedaq::datahandlinglibs::BinarySearchQueueModel,dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(SkipListLatencyBufferModel_DAPHNESuperChunk)
{
  test_queue_model<dunedaq::datahandlinglibs::SkipListLatencyBufferModel,dunedaq::fdreadoutlibs::types::DAPHNESuperChunkTypeAdapter>();
}
*/

BOOST_AUTO_TEST_SUITE_END()


