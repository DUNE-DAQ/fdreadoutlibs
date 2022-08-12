/**
 * @file RegisterToChannelNumber.hpp Convert from WIB1 data AVX register position to offline channel numbers
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_TPG_REGISTERTOCHANNELNUMBER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_TPG_REGISTERTOCHANNELNUMBER_HPP_

#include "detchannelmaps/TPCChannelMap.hpp"
#include "detdataformats/wib/WIBFrame.hpp"
#include "readoutlibs/ReadoutTypes.hpp"
#include "FrameExpand.hpp"
#include "TPGConstants.hpp"

#include <boost/chrono/duration.hpp>
#include <chrono>
#include <sys/types.h>
#include <vector>

namespace swtpg {

struct RegisterChannelMap
{
  // TODO: Make these the right size
  uint collection[256];
  uint induction[256];
};

/**
 * The code that expands ADCs from WIB frame format into AVX registers
 * for the software TP generation puts the channels in the output
 * registers in some order that is convenient for the expansion code,
 * but doesn't have any particular pattern to it. So we need to map
 * from position-in-register to offline channel number. This function
 * creates that map
 */
RegisterChannelMap
get_register_to_offline_channel_map(const dunedaq::detdataformats::wib::WIBFrame* frame,
                                    std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap>& ch_map)
{
  auto start_time = std::chrono::steady_clock::now();

  // Find the lowest offline channel number of all the channels in the input frame
  uint min_ch = UINT_MAX;
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib::WIBFrame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->get_wib_header()->crate_no, frame->get_wib_header()->slot_no, frame->get_wib_header()->fiber_no, ich);
    min_ch = std::min(min_ch, offline_ch);
  }
  TLOG_DEBUG(0) << "get_register_to_offline_channel_map for crate " << frame->get_wib_header()->crate_no << " slot "
                << frame->get_wib_header()->slot_no << " fiber " << frame->get_wib_header()->fiber_no << ". min_ch is "
                << min_ch;
  // Now set each of the channels in our test frame to their
  // corresponding offline channel number, minus the minimum channel
  // number we just calculated (so we don't overflow the 12 bits we
  // have available)
  dunedaq::fdreadoutlibs::types::WIB_SUPERCHUNK_STRUCT superchunk;
  memset(superchunk.data, 0, sizeof(dunedaq::fdreadoutlibs::types::WIB_SUPERCHUNK_STRUCT));

  dunedaq::detdataformats::wib::WIBFrame* test_frame =
    reinterpret_cast<dunedaq::detdataformats::wib::WIBFrame*>(&superchunk);
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib::WIBFrame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->get_wib_header()->crate_no, frame->get_wib_header()->slot_no, frame->get_wib_header()->fiber_no, ich);
    test_frame->set_channel(ich, offline_ch - min_ch);
  }

  // Expand the test frame, so the offline channel numbers are now in the relevant places in the output registers
  swtpg::MessageRegistersCollection collection_registers;
  swtpg::MessageRegistersInduction induction_registers;
  expand_message_adcs_inplace(&superchunk, &collection_registers, &induction_registers);

  RegisterChannelMap ret;
  for (size_t i = 0; i < 6 * SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram
    size_t index = (i/16)*16*12 + (i%16);
    ret.collection[i] = collection_registers.uint16(index) + min_ch;
  }
  for (size_t i = 0; i < 10 * SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram
    size_t index = (i/16)*16*12 + (i%16);
    ret.induction[i] = induction_registers.uint16(index) + min_ch;
  }

  auto end_time = std::chrono::steady_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
  TLOG_DEBUG(0) << "get_register_to_offline_channel_map built map in " << dur << "us";
  return ret;
}

RegisterChannelMap
get_register_to_offline_channel_map(const dunedaq::detdataformats::wib::WIBFrame* frame, std::string channel_map_name)
{
  auto ch_map = dunedaq::detchannelmaps::make_map(channel_map_name);
  return get_register_to_offline_channel_map(frame, ch_map);
}



/**
 * The code that expands ADCs from WIB2 frame format into AVX registers
 * for the software TP generation puts the channels in the output
 * registers in some order that is convenient for the expansion code,
 * but doesn't have any particular pattern to it. So we need to map
 * from position-in-register to offline channel number. This function
 * creates that map
 */
RegisterChannelMap
get_register_to_offline_channel_map_wib2(const dunedaq::detdataformats::wib2::WIB2Frame* frame,
                                    std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap>& ch_map)
{
  auto start_time = std::chrono::steady_clock::now();

  // Find the lowest offline channel number of all the channels in the input frame
  uint min_ch = UINT_MAX;
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib2::WIB2Frame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->get_wib2_header()->crate, frame->get_wib2_header()->slot, frame->get_wib2_header()->link, ich);
    min_ch = std::min(min_ch, offline_ch);
  }
  TLOG_DEBUG(0) << "get_register_to_offline_channel_map for crate " << frame->get_wib2_header()->crate << " slot "
                << frame->get_wib2_header()->slot << " link " << frame->get_wib2_header()->link << ". min_ch is "
                << min_ch;
  // Now set each of the channels in our test frame to their
  // corresponding offline channel number, minus the minimum channel
  // number we just calculated (so we don't overflow the 12 bits we
  // have available)
  dunedaq::fdreadoutlibs::types::WIB2_SUPERCHUNK_STRUCT superchunk;
  memset(superchunk.data, 0, sizeof(dunedaq::fdreadoutlibs::types::WIB2_SUPERCHUNK_STRUCT));

  dunedaq::detdataformats::wib2::WIB2Frame* test_frame =
    reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(&superchunk);
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib2::WIB2Frame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->get_wib2_header()->crate, frame->get_wib2_header()->slot, frame->get_wib2_header()->link, ich);
    test_frame->set_channel(ich, offline_ch - min_ch);
  }

  // Expand the test frame, so the offline channel numbers are now in the relevant places in the output registers
  swtpg::MessageRegistersCollection collection_registers;
  swtpg::MessageRegistersInduction induction_registers;
  expand_message_adcs_inplace_wib2(&superchunk, &collection_registers, &induction_registers);

  RegisterChannelMap ret;
  for (size_t i = 0; i < 6 * SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram
    size_t index = (i/16)*16*12 + (i%16);
    ret.collection[i] = collection_registers.uint16(index) + min_ch;
  }
  for (size_t i = 0; i < 10 * SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram
    size_t index = (i/16)*16*12 + (i%16);
    ret.induction[i] = induction_registers.uint16(index) + min_ch;
  }

  auto end_time = std::chrono::steady_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
  TLOG_DEBUG(0) << "get_register_to_offline_channel_map built map in " << dur << "us";
  return ret;
}

RegisterChannelMap
get_register_to_offline_channel_map(const dunedaq::detdataformats::wib2::WIB2Frame* frame, std::string channel_map_name)
{
  auto ch_map = dunedaq::detchannelmaps::make_map(channel_map_name);
  return get_register_to_offline_channel_map(frame, ch_map);
}



} // namespace swtpg

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_TPG_REGISTERTOCHANNELNUMBER_HPP_

