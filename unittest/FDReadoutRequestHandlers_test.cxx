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
#include <cmath>
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
  uint64_t n_frames = test_element.get_num_frames();
  uint64_t ticks_per_frame = TypeAdapter::expected_tick_difference;
  uint64_t ticks_between = ticks_per_frame*n_frames;


  // function for testing get_fragment_pieces in cases where there are no skips in the buffer
  auto test_req_bounds = [&](std::shared_ptr<BufferType<TypeAdapter>> buffer,
                             uint64_t start_win, uint64_t end_win,
                             std::set<size_t> objects_skipped={}){

      //make the error registry we need
      auto errorRegistry = std::make_unique<dunedaq::datahandlinglibs::FrameErrorRegistry>();

      //create the request handler
      DefaultRequestHandler requestHandler(buffer,errorRegistry);

      //call the get_fragment_pieces
      auto dfmessage = dunedaq::dfmessages::DataRequest();
      auto req_res = typename DefaultRequestHandler::RequestResult(DefaultRequestHandler::ResultCode::kUnknown,dfmessage);
      auto ret = requestHandler.get_fragment_pieces(start_win,end_win,req_res);

      //check that the return code is correct.
      BOOST_CHECK_EQUAL(req_res.result_code,DefaultRequestHandler::ResultCode::kFound);

      //check that the first and last returned blocks are note empty
      BOOST_REQUIRE_GT(ret.front().second,0);
      BOOST_REQUIRE_GT(ret.back().second,0);

      //grab the (first) timestamps of the first and last frames
      uint64_t first_ts = reinterpret_cast<const TypeAdapter::FrameType*>(ret.front().first)->get_timestamp();
      uint64_t last_ts = reinterpret_cast<const TypeAdapter::FrameType*>((char*)(ret.back().first)+(ret.back().second)-sizeof(typename TypeAdapter::FrameType))->get_timestamp();

      //general check:
      // first_ts <= start_win < first_ts + ticks_per_frame
      // last_ts < end_win <= last_ts + ticks_per_frame
      if(objects_skipped.size()==0)
          BOOST_CHECK_MESSAGE(first_ts<=start_win,
                              "first_frame_ts{" << first_ts << "} <= start_win{" << start_win << "}");
      BOOST_CHECK_MESSAGE(start_win<first_ts+ticks_per_frame,
                          "start_win{" << start_win << "} < first_frame_ts+ticks_per_frame{" << first_ts+ticks_per_frame << "}");
      BOOST_CHECK_MESSAGE(last_ts<end_win,
                          "Check last_frame_ts{" << last_ts << "} < end_win{" << end_win << "}");
      if(objects_skipped.size()==0)
          BOOST_CHECK_MESSAGE(end_win<=last_ts+ticks_per_frame,
                              "end_win{" << end_win << "} <= last_frame_ts+ticks_per_frame{" << last_ts+ticks_per_frame << "}");

      //expected timestamps for begin of fragment and 'end' of fragment, assuming no skipping
      auto expected_start = TypeAdapter::expected_tick_difference *
              (uint64_t) std::floor((float) (start_win) / (float) (ticks_per_frame));
      auto expected_end = TypeAdapter::expected_tick_difference *
              (uint64_t) std::ceil((float) (end_win) / (float) (ticks_per_frame));

      //correct for the cases there that object has been skipped
      auto expected_start_obj = expected_start - expected_start%ticks_between;
      while(objects_skipped.count(expected_start_obj/ticks_between)>0) {
          expected_start += ticks_per_frame;
          expected_start_obj = expected_start - expected_start%ticks_between;
      }

      auto expected_end_obj = expected_end - ticks_per_frame - (expected_end-ticks_per_frame)%ticks_between;
      while(objects_skipped.count(expected_end_obj/ticks_between)>0) {
          expected_end -= ticks_per_frame;
          expected_start_obj = expected_end - ticks_per_frame - (expected_end-ticks_per_frame)%ticks_between;
      }

      //specfic check: are values for this request what we expect
      BOOST_CHECK_MESSAGE(first_ts == expected_start,
                          "Fragment start ts {" << first_ts << "} is expected value {" << expected_start << "}");
      BOOST_CHECK_MESSAGE((last_ts + TypeAdapter::expected_tick_difference) == expected_end,
                          "Fragment 'end' ts {" << last_ts + TypeAdapter::expected_tick_difference
                          << "} is expected value {" << expected_end << "}");
  };

  /*
   * Unskipped buffer should have elements with index [0, 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ]
   *                                   and timestamps [0,1*T,2*T,3*T,4*T,5*T,6*T,7*T,8*T,9*T]
   * where T = DTS ticks between successive elements (tick_diff_per_frame * n_frames_per_obj_in_buffer)
   */
  BOOST_TEST_MESSAGE("Testing buffer without skips...");
  auto buffer_noskip = std::make_shared< BufferType<TypeAdapter> >();
  fill_buffer<BufferType,TypeAdapter>(buffer_noskip,0,10);
  print_buffer<BufferType,TypeAdapter>(buffer_noskip,"noskip");

  test_req_bounds(buffer_noskip,ticks_between*2,ticks_between*5);
  test_req_bounds(buffer_noskip,ticks_between*3/2,ticks_between*9/2);
  test_req_bounds(buffer_noskip,ticks_between*11/5,ticks_between*21/5);
  test_req_bounds(buffer_noskip,ticks_between*2+1,ticks_between*5+1);

  /*
   * Skipped buffer should have elements with index [0, 1 , 2 , 3 , 4 , 5 , 6 , 7 ,  8 ,  9 ]
   *                                 and timestamps [0,1*T,2*T,5*T,6*T,7*T,8*T,9*T,10*T,11*T]
   * where T = DTS ticks between successive elements (tick_diff_per_frame * n_frames_per_obj_in_buffer)
  */

  BOOST_TEST_MESSAGE("Testing buffer with skips...");
  std::set<size_t> obj_to_skip = {2,3};
  auto buffer_skip = std::make_shared< BufferType<TypeAdapter> >();
  fill_buffer<BufferType,TypeAdapter>(buffer_skip,0,10,obj_to_skip);
  print_buffer<BufferType,TypeAdapter>(buffer_skip,"skip");

  test_req_bounds(buffer_skip,ticks_between*2,ticks_between*5,obj_to_skip);
  test_req_bounds(buffer_skip,ticks_between*3/2,ticks_between*9/2,obj_to_skip);
  test_req_bounds(buffer_skip,ticks_between*11/5,ticks_between*21/5,obj_to_skip);
  test_req_bounds(buffer_skip,ticks_between*2+1,ticks_between*5+1,obj_to_skip);

}

BOOST_AUTO_TEST_SUITE(FDReadoutRequestHandlers_test)

BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DUNEWIBEth)
{
    test_request_model<
            dunedaq::datahandlinglibs::FixedRateQueueModel,
            dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}

BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DUNEWIBEth)
{
    test_request_model<
            dunedaq::datahandlinglibs::BinarySearchQueueModel,
            dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DAPHNEStreamSuperChunk)
{
    test_request_model<
            dunedaq::datahandlinglibs::FixedRateQueueModel,
            dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}

BOOST_AUTO_TEST_CASE(BinarySearchQueueModel_DAPHNEStreamSuperChunk)
{
    test_request_model<
            dunedaq::datahandlinglibs::BinarySearchQueueModel,
            dunedaq::fdreadoutlibs::types::DAPHNEStreamSuperChunkTypeAdapter>();
}

BOOST_AUTO_TEST_CASE(SkipListLatencyBufferModel_DAPHNESuperChunk)
{
    test_request_model<
            dunedaq::datahandlinglibs::SkipListLatencyBufferModel,
            dunedaq::fdreadoutlibs::types::DAPHNESuperChunkTypeAdapter>();
}


BOOST_AUTO_TEST_SUITE_END()


