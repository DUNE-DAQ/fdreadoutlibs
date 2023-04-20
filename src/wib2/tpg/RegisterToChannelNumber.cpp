/**
 * @file RegisterToChannelNumber.hpp Convert from WIB1 data AVX register position to offline channel numbers
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/wib2/tpg/RegisterToChannelNumber.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/ReadoutTypes.hpp"
#include "fdreadoutlibs/wib2/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wib2/tpg/TPGConstants_wib2.hpp"

#include <boost/chrono/duration.hpp>
#include <chrono>
#include <sys/types.h>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

namespace swtpg_wib2 {

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
                                    std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap>& ch_map,
                                    int registers_selection
                                    )
{

  using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
  using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

  auto start_time = std::chrono::steady_clock::now();

  // Find the lowest offline channel number of all the channels in the input frame
  uint min_ch = UINT_MAX;
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib2::WIB2Frame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->header.crate, frame->header.slot, frame->header.link, ich);
    TLOG_DEBUG(TLVL_BOOKKEEPING) << " offline_ch " << offline_ch; 
    min_ch = std::min(min_ch, offline_ch);
  }
  TLOG() << "get_register_to_offline_channel_map_wib2 for crate " << frame->header.crate << " slot "
                << frame->header.slot << " link " << frame->header.link << ". min_ch is "
                << min_ch;
  // Now set each of the channels in our test frame to their
  // corresponding offline channel number, minus the minimum channel
  // number we just calculated (so we don't overflow the 12 bits we
  // have available)
  dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter superchunk;
  memset(superchunk.data, 0, sizeof(dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter));

  dunedaq::detdataformats::wib2::WIB2Frame* test_frame =
    reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(&superchunk);
  for (size_t ich = 0; ich < dunedaq::detdataformats::wib2::WIB2Frame::s_num_ch_per_frame; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_fiber_chan(
      frame->header.crate, frame->header.slot, frame->header.link, ich);
      test_frame->set_adc(ich, offline_ch - min_ch);
  }

  // Expand the test frame, so the offline channel numbers are now in the relevant places in the output registers
  swtpg_wib2::MessageRegisters register_array;
  expand_wib2_adcs(&superchunk, &register_array, registers_selection); 


  RegisterChannelMap ret;
  for (size_t i = 0;  i < swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram
    size_t index = (i/16)*16*12 + (i%16);
    ret.channel[i] = register_array.uint16(index) + min_ch;
  }

  auto end_time = std::chrono::steady_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "get_register_to_offline_channel_map_wib2 built map in " << dur << "us";
  return ret;
}

RegisterChannelMap
get_register_to_offline_channel_map_wib2(
  const dunedaq::detdataformats::wib2::WIB2Frame* frame,
  std::string channel_map_name,
  int registers_selection
  )
{
  auto ch_map = dunedaq::detchannelmaps::make_map(channel_map_name);
  return get_register_to_offline_channel_map_wib2(frame, ch_map, registers_selection);
}



} // namespace swtpg_wib2

