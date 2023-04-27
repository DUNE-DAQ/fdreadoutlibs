/**
 * @file SwtpgBase.hpp  Base class for the SWTPG algorithms
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef TEST_SWTPG_BASE_HPPP_
#define TEST_SWTPG_BASE_HPPP_

#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib2/tpg/RegisterToChannelNumber.hpp"

class SwtpgBase {
    
public:
   virtual void reset(bool first_hit, int threshold_value, std::string channel_map_name) {};
   virtual void find_hits(const dunedaq::fdreadoutlibs::types::DUNEWIBSuperChunkTypeAdapter* fp, bool first_hit) {};

   // Total sum of the number of found hits. Just for debugging, to be removed later
   size_t m_total_swtpg_hits = 0;


protected:
  uint16_t m_tpg_threshold = 5; // units of sigma 
  const uint8_t m_tpg_tap_exponent = 6; 
  const int m_tpg_multiplier = 1 << m_tpg_tap_exponent;  // 64
  std::vector<int16_t> m_tpg_taps;  // firwin_int(7, 0.1, multiplier)
  int16_t* m_tpg_taps_p = nullptr;
  uint16_t* m_primfind_dest = nullptr;  

  // Map from expanded AVX register position to offline channel number
  swtpg_wib2::RegisterChannelMap m_register_channel_map; 
  std::shared_ptr<dunedaq::detchannelmaps::TPCChannelMap> m_channel_map;

  // AAA: Hardcoded value but for now it is fine
  std::string m_ch_map_name = "HDColdboxChannelMap";


  int m_link_no; 
  int m_crate_no;
  int m_slot_no;




};

#endif // TEST_SWTPG_BASE_HPPP_
