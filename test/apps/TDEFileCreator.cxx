/**
 * @file TDEFileCreator.cxx Create a binary file with TDE Frames for reading it from readout
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "fddetdataformats/TDE16Frame.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>

using namespace dunedaq;

int main()
{
  int num_batches = 5;

  std::ofstream out("frames.bin", std::ios::out | std::ios::binary);

  std::random_device rd;
  std::default_random_engine rng(rd());


  for (int batch = 0; batch < num_batches; batch++)
  {
    std::vector<fddetdataformats::TDE16Frame> v;
    uint64_t timestamp = batch;
    for (int amc = 0; amc < 12; amc++)
    {
      for (int i = 0; i < 64; i++)
      {
        fddetdataformats::TDE16Frame fr;
        // Timestamp
        fr.set_timestamp(timestamp);
        // Channel
        fr.get_tde_header()->slot = amc;
        fr.get_tde_header()->link = i;
        fr.set_adc_sample(batch,0);
        v.push_back(fr);
      }
    }
    std::shuffle(v.begin(), v.end(), rng);
    for (auto& fr: v) {
        out.write(reinterpret_cast<char*>(&fr), sizeof(fddetdataformats::TDE16Frame));
    }
  }
  out.close();

}
