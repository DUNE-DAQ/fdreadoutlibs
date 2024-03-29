/**
 * @file design_fir.h Functions to design a lowpass filter with the
 * given number of taps and cutoff, as a fraction of the
 * Nyquist frequency. Copied from the corresponding functions
 * in scipy and converted to C++
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_DESIGNFIR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_DESIGNFIR_HPP_

#include <cmath>
#include <cstdint>
#include <vector>

namespace swtpg_wibeth {

// constexpr double pi() { return 3.1416; }
constexpr double
pi()
{
  return std::atan(1) * 4;
}

std::vector<double>
hamming(int M);

double
sinc(double x);

std::vector<double>
firwin(int N, double cutoff);

std::vector<int16_t>
firwin_int(int N, double cutoff, const int multiplier);

} // namespace swtpg_wibeth

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_TPG_DESIGNFIR_HPP_
