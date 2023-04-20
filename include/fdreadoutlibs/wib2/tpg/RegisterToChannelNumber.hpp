/**
 * @file RegisterToChannelNumber.hpp Convert from WIB1 data AVX register position to offline channel numbers
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_REGISTERTOCHANNELNUMBER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_REGISTERTOCHANNELNUMBER_HPP_

#include "detchannelmaps/TPCChannelMap.hpp"
#include "detdataformats/wib2/WIB2Frame.hpp"


namespace swtpg_wib2 {

struct RegisterChannelMap
{
  uint channel[swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER];
};

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
                                    );

RegisterChannelMap
get_register_to_offline_channel_map_wib2(
  const dunedaq::detdataformats::wib2::WIB2Frame* frame,
  std::string channel_map_name,
  int registers_selection
  );


} // namespace swtpg_wib2

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPG_REGISTERTOCHANNELNUMBER_HPP_

