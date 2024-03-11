/**
 * @file RegisterToChannelNumber.hpp Convert from WIB1 data AVX register position to offline channel numbers
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fdreadoutlibs/wibeth/tpg/RegisterToChannelNumber.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/ReadoutTypes.hpp"
#include "fdreadoutlibs/wibeth/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wibeth/tpg/TPGConstants_wibeth.hpp"

#include <boost/chrono/duration.hpp>
#include <chrono>
#include <sys/types.h>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

namespace swtpg_wibeth {

/**
 * The code that expands ADCs from WIBEth frame format into AVX registers
 * for the software TP generation puts the channels in the output
 * registers in some order that is convenient for the expansion code,
 * but doesn't have any particular pattern to it. So we need to map
 * from position-in-register to offline channel number. This function
 * creates that map
 */
RegisterChannelMap
get_register_to_offline_channel_map_wibeth(const dunedaq::fddetdataformats::WIBEthFrame* frame,
                                    std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap>& ch_map
                                    )
{

  using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
  using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

  auto start_time = std::chrono::steady_clock::now();

  // Find the lowest offline channel number of all the channels in the input frame
  uint min_ch = UINT_MAX;
  for (size_t ich = 0; ich < dunedaq::fddetdataformats::WIBEthFrame::s_channels_per_half_femb; ++ich) {
    auto offline_ch = ch_map->get_offline_channel_from_crate_slot_stream_chan(
      frame->daq_header.crate_id, frame->daq_header.slot_id, frame->daq_header.stream_id, ich);
    //TLOG_DEBUG(TLVL_BOOKKEEPING) << " offline_ch " << offline_ch; 
    min_ch = std::min(min_ch, offline_ch);
  }
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "get_register_to_offline_channel_map_wibeth for crate " << frame->daq_header.crate_id << " slot "
                << frame->daq_header.slot_id << " stream " << frame->daq_header.stream_id << ". min_ch is "
                << min_ch;
  // Now set each of the channels in our test frame to their
  // corresponding offline channel number, minus the minimum channel
  // number we just calculated (so we don't overflow the 12 bits we
  // have available)
  dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter wibeth_frame;
  memset(wibeth_frame.data, 0, sizeof(dunedaq::fdreadoutlibs::types::DUNEWIBEthTypeAdapter));

  
  dunedaq::fddetdataformats::WIBEthFrame* test_frame =
    reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(&wibeth_frame);
  for (size_t time_sample = 0; time_sample < dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame; ++time_sample) {  
    for (size_t ich = 0; ich < dunedaq::fddetdataformats::WIBEthFrame::s_channels_per_half_femb; ++ich) { 
      //s_channels_per_half_femb is 64 because there is only one half FEMB
      auto offline_ch = ch_map->get_offline_channel_from_crate_slot_stream_chan(
        frame->daq_header.crate_id, frame->daq_header.slot_id, frame->daq_header.stream_id, ich);
        test_frame->set_adc(ich, time_sample, offline_ch - min_ch);
    }
  }

  // Expand the test frame, so the offline channel numbers are now in the relevant places in the output registers
  swtpg_wibeth::MessageRegisters register_array;
  expand_wibeth_adcs(&wibeth_frame, &register_array); 


  RegisterChannelMap ret;
 
  // Define the following variables for convenience
  //int NREGISTERS = swtpg_wibeth::NUM_REGISTERS_PER_FRAME;
  int SAMPLES_PER_REGISTER = swtpg_wibeth::SAMPLES_PER_REGISTER;
  int TIME_WINDOW_NUM_FRAMES = dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame;

  for (size_t i = 0;  i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
    // expand_message_adcs_inplace reorders the output so
    // adjacent-in-time registers are adjacent in memory, hence the
    // need for this indexing. See the comment in that function for a
    // diagram

    const size_t register_offset = i % SAMPLES_PER_REGISTER;
    const size_t register_index = i / SAMPLES_PER_REGISTER;
    const size_t register_t0_start = register_index * SAMPLES_PER_REGISTER * TIME_WINDOW_NUM_FRAMES;

    // AAA: in this function we only want the offline channel number
    // so we do not care about the other time samples in the frame
    // that's why itime is zero
    size_t itime = 0; 
    const size_t msg_index = itime / TIME_WINDOW_NUM_FRAMES;
    const size_t msg_time_offset = itime % TIME_WINDOW_NUM_FRAMES;
    const size_t msg_start_index = msg_index * (swtpg_wibeth::ADCS_SIZE) / sizeof(uint16_t); // NOLINT
    const size_t offset_within_msg = register_t0_start + SAMPLES_PER_REGISTER * msg_time_offset + register_offset;

    // The index in uint16_t of the start of the message we want. 
    const size_t index = msg_start_index + offset_within_msg;
    int16_t out_val = register_array.uint16(index);
    ret.channel[i] = out_val + min_ch;

    //std::cout << " index: " << index << "    value:   " << out_val << std::endl;



  }

  auto end_time = std::chrono::steady_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "get_register_to_offline_channel_map_wibeth built map in " << dur << "us";
  return ret;
}

RegisterChannelMap
get_register_to_offline_channel_map_wibeth(
  const dunedaq::fddetdataformats::WIBEthFrame* frame,
  std::string channel_map_name
  )
{
  auto ch_map = dunedaq::detchannelmaps::make_map(channel_map_name);
  return get_register_to_offline_channel_map_wibeth(frame, ch_map);
}



} // namespace swtpg_wibeth

