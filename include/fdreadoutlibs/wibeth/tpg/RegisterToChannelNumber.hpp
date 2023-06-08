/**
 * @file RegisterToChannelNumber.hpp Convert from WIB1 data AVX register position to offline channel numbers
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_REGISTERTOCHANNELNUMBER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_REGISTERTOCHANNELNUMBER_HPP_

#include "detchannelmaps/TPCChannelMap.hpp"

#include "fddetdataformats/WIBEthFrame.hpp"
#include "fdreadoutlibs/wibeth/tpg/TPGConstants_wibeth.hpp"


namespace swtpg_wibeth {

struct RegisterChannelMap
{
  uint channel[swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER];
};

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
                                    );

RegisterChannelMap
get_register_to_offline_channel_map_wibeth(
  const dunedaq::fddetdataformats::WIBEthFrame* frame,
  std::string channel_map_name
  );


} // namespace swtpg_wibeth

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_REGISTERTOCHANNELNUMBER_HPP_

