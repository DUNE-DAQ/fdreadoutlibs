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

#include "datahandlinglibs/models/FixedRateQueueModel.hpp"
#include "datahandlinglibs/models/BinarySearchQueueModel.hpp"
#include "datahandlinglibs/models/SkipListLatencyBufferModel.hpp"
#include "datahandlinglibs/concepts/LatencyBufferConcept.hpp"


#include "boost/test/unit_test.hpp"

#include <iostream>
#include <sstream>
#include <set>
#include <iterator>

//general test buffer struct that can be reused
template <template<class> class BufferType, class TypeAdapter>
void fill_buffer( BufferType<TypeAdapter>& buffer,
		  uint64_t const init_timestamp=0,
		  size_t const n_obj=10,
		  std::set<size_t> const obj_to_skip = {})
{
  //allocate mem in buffer
  buffer.allocate_memory(n_obj);
  
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
      buffer.write(std::move(frame));
      ++i_obj;
    }
    ++i;
    
    
  } // end while(obj_in_buffer<n_obj)
}

template <template<class> class BufferType, class TypeAdapter>
void test_lower_bound( BufferType<TypeAdapter>& buffer,
		       uint64_t test_ts,
		       uint32_t expected_idx,
		       bool with_errors=false)
{
  TypeAdapter test_element; test_element.set_timestamp(test_ts);
  typename BufferType<TypeAdapter>::Iterator expected_el =  buffer.begin();
  for(size_t i=0; i<expected_idx; ++i){
    ++expected_el;
  }

  //get our lower_bound call
  auto return_el = buffer.lower_bound(test_element,with_errors);

  //loop through to get the previous element to return_el
  typename BufferType<TypeAdapter>::Iterator scan_el =  buffer.begin();
  typename BufferType<TypeAdapter>::Iterator prev_el =  buffer.begin();
  while(scan_el!=return_el){
    ++scan_el;
    if(scan_el==return_el) break;
    ++prev_el;
  }
  
  //check that expected and return timestamps agree
  BOOST_CHECK_MESSAGE(expected_el->get_timestamp()==return_el->get_timestamp(),
		      "Expected ts{" << expected_el->get_timestamp() << "} == return ts{" << return_el->get_timestamp() << "} for test_ts=" << test_ts);
  
  //check that we satisfy the lower bound condition
  BOOST_CHECK_MESSAGE(return_el->get_timestamp()>=test_ts,
		      "Returned ts{" << return_el->get_timestamp() << "} is >= test_ts{" << test_ts << "}");
  BOOST_CHECK_MESSAGE((prev_el->get_timestamp()<test_ts || return_el==buffer.begin()),
		      "Prev ts{" << prev_el->get_timestamp() << "} is < test_ts{" << test_ts << "} (or lower bound is begin of buffer)");
}

template <template<class> class BufferType, class TypeAdapter>
void print_buffer(BufferType<TypeAdapter>& buffer, std::string desc)
{
  std::stringstream ss;
  ss << "Buffer (" << desc << "): ";
  typename BufferType<TypeAdapter>::Iterator iter=buffer.begin();
  while(iter!=buffer.end()){
    ss << iter->get_timestamp() << " ";
    ++iter;
  }
  BOOST_TEST_MESSAGE(ss.str());

}


template <template<class> class BufferType, class TypeAdapter>
void test_queue_model()
{

  //create our buffer
  BufferType<TypeAdapter> buffer_noskip,buffer_skip;

  //some testing vars
  TypeAdapter test_element;
  uint64_t ticks_between = TypeAdapter::expected_tick_difference*test_element.get_num_frames();

  /*
   * Unskipped buffer should have elements with index [0, 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ]
   *                                   and timestamps [0,1*T,2*T,3*T,4*T,5*T,6*T,7*T,8*T,9*T]
   * where T = DTS ticks between successive elements (tick_diff_per_frame * n_frames_per_obj_in_buffer)
   */
  BOOST_TEST_MESSAGE("Testing buffer without skips...");
  fill_buffer(buffer_noskip,0,10);
  print_buffer(buffer_noskip,"noskip");
  
  // get lower bound on aligned element
  // should return the exact value
  test_lower_bound<BufferType,TypeAdapter>(buffer_noskip,ticks_between*2,2);

  // get lower bound when maximally unaligned
  // should get the next element up
  test_lower_bound<BufferType,TypeAdapter>(buffer_noskip,ticks_between*5/2,3);
  
  // get lower bound when minimally unaligned, next instance
  test_lower_bound<BufferType,TypeAdapter>(buffer_noskip,ticks_between+1,2);

  /*
   * Skipped buffer should have elements with index [0, 1 , 2 , 3 , 4 , 5 , 6 , 7 ,  8 ,  9 ]
   *                                 and timestamps [0,1*T,2*T,5*T,6*T,7*T,8*T,9*T,10*T,11*T]
   * where T = DTS ticks between successive elements (tick_diff_per_frame * n_frames_per_obj_in_buffer)
  */
  BOOST_TEST_MESSAGE("Testing buffer with skips...");
  std::set<size_t> obj_to_skip = {2,3};
  fill_buffer(buffer_skip,0,10,obj_to_skip);
  print_buffer(buffer_skip,"skip");

  // get lower bound on aligned but skipped element
  // should return next available
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*2,2,true);
  // should be unaffected
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between,1,true);

  // get lower bound when maximally unaligned
  // should return next available
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*3/2,2,true);
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*5/2,2,true);
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*7/2,2,true);
  // should be unaffected
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*1/2,1,true);
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*9/2,3,true);
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*11/2,4,true);

  // get lower bound when minimally unaligned, next instance
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between+1,2,true);  
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,ticks_between*2+1,2,true);  
  // should be unaffected
  test_lower_bound<BufferType,TypeAdapter>(buffer_skip,1,1,true);
}

BOOST_AUTO_TEST_SUITE(FDReadoutTypeAdaptersBuffers_test)

BOOST_AUTO_TEST_CASE(FixedRateQueueModel_DUNEWIBEth)
{
  test_queue_model<dunedaq::datahandlinglibs::FixedRateQueueModel,dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter>();
}
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

BOOST_AUTO_TEST_SUITE_END()


