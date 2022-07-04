// dvargas 2022

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <atomic>
#include <chrono>
#include <random>
#include <vector>
#include <string>

#include "readoutlibs/utils/RateLimiter.hpp"
#include "detdataformats/tde/TDE16Frame.hpp"
#include "detdataformats/tde/TDE12Frame.hpp"
using namespace dunedaq::readoutlibs;
using namespace dunedaq::detdataformats::tde;
using namespace std;

static constexpr int rate_12bits = 328;                              // in Hz = 3.054 in ms
static constexpr int rate_16bits = 437;                              // in Hz = 2.290 in ms


int main()
{
  TDE16Header tdeheader;
  TDE16Frame tde16frame;
  ADC16Data adc16data;
  TDE12Frame tde12frame;
  ADC12Data adc12data;

  FILE *fp_16, *fp_12;

  bool killswitch16 { true };
  bool killswitch12 { true };
  uint16_t new_adc_val_12bits = 0x43;
  uint16_t new_adc_val_16bits = 0x72;
  uint64_t toptime = 0xfa;
  uint64_t time = 0;

  // Seting the header
  tdeheader.version = 0x1;
  tdeheader.det_id = 0x2e;
  tdeheader.crate = 0x339;
  tdeheader.slot = 0x7;
  tdeheader.link = 0x23;
  tdeheader.tde_header = 0x61;
  tdeheader.tde_errors = 0xfc5;

  for(int i=0; i<tot_adc16_samples; i++) { tde16frame.set_adc_samples(new_adc_val_16bits, i); }
  for(int i=0; i<tot_adc12_samples; i++) { tde12frame.set_adc(i, new_adc_val_12bits); }

  // using the RateLimiter in Hz for 16 and 12 bits
  auto limiter16 = RateLimiter(rate_16bits); 
  limiter16.init();
  auto limiter12 = RateLimiter(rate_12bits); 
  limiter12.init();

  while (killswitch16) 
  {
    for(uint64_t i = 0; i < toptime; i++)
    {
      time += i;
      tdeheader.set_timestamp(time);
      tdeheader.set_TAItime(time);
      limiter16.limit();

      std::string str = std::to_string(i);
      std::string filename16 = "/nfs/sw/TDE-testing/binary_16bits_" + str + ".bit";
      fp_16 = fopen(filename16.c_str(), "wb");
      if (fp_16 == NULL) { std::cout << "Could not open output file fp_16" << std::endl; }
      fwrite(&tdeheader, sizeof(tdeheader), 1, fp_16);
      fwrite(&adc16data, sizeof(adc16data), 1, fp_16);
      fclose(fp_16);
    }
    killswitch16 = false;
  }

  while (killswitch12) 
  {
    for(uint64_t i = 0; i < toptime; i++)
    {
      time += i;
      tdeheader.set_timestamp(time);
      tdeheader.set_TAItime(time);
      limiter12.limit();

      std::string str = std::to_string(i);
      std::string filename12 = "/nfs/sw/TDE-testing/binary_12bits_" + str + ".bit";
      fp_12 = fopen(filename12.c_str(), "wb");
      if (fp_12 == NULL) { std::cout << "Could not open output file fp_12" << std::endl; }
      fwrite(&tdeheader, sizeof(tdeheader), 1, fp_12);
      fwrite(&adc12data, sizeof(adc12data), 1, fp_12);
      fclose(fp_12);
    }
    killswitch12 = false;
  }

  std::cout << "Finish while stamenet" << std::endl;

  return 0;
}
