#ifndef WIB2_TPG_UTILITIES_HPP_
#define WIB2_TPG_UTILITIES_HPP_



#include <immintrin.h>
//#include "TPGConstants_wib2.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>
#include <cstdio> // For printf
#include <array>
#include <chrono>
#include <stdio.h>
#include <stdint.h>
#include <memory>


constexpr double
pi()
{
  return std::atan(1) * 4;
}


std::vector<double>
hamming(int M)
{
  std::vector<double> ret;
  for (int n = 0; n < M; ++n) {
    ret.push_back(0.54 - 0.46 * cos(2.0 * pi() * n / (M - 1)));
  }
  return ret;
}

double
sinc(double x)
{
  if (x == 0) {
    return 1;
  }
  return sin(pi() * x) / (pi() * x);
}

std::vector<double>
firwin(int N, double cutoff)
{
  int alpha = N / 2;
  std::vector<double> window = hamming(N);
  std::vector<double> ret;
  double sum = 0;
  for (int m = 0; m < N; ++m) {
    double val = window[m] * sinc(cutoff * (m - alpha));
    ret.push_back(val);
    sum += val;
  }
  for (int m = 0; m < N; ++m) {
    ret[m] /= sum;
  }
  return ret;
}

std::vector<int16_t>
firwin_int(int N, double cutoff, const int multiplier)
{
  std::vector<double> coeffs_double(firwin(N, cutoff));
  // coeffs_double.push_back(0);

  std::vector<int16_t> coeffs(N, 0);
  for (int i = 0; i < N; ++i) {
    coeffs[i] = round(multiplier * coeffs_double[i]);
  }
  return coeffs;
}

// A little wrapper around an array of 256-bit registers, so that we
// can explicitly access it as an array of 256-bit registers or as an
// array of uint16_t
template<size_t N>
class RegisterArray
{
public:
    // Get the value at the ith position as a 256-bit register
    inline __m256i ymm(size_t i) const { return _mm256_lddqu_si256(reinterpret_cast<const __m256i*>(m_array)+i); }
    inline void set_ymm(size_t i, __m256i val) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_array)+i, val); }

    inline uint16_t uint16(size_t i) const { return m_array[i]; }
    inline void set_uint16(size_t i, uint16_t val) { m_array[i]=val; }

    // Access the jth entry in the ith register
    inline uint16_t uint16(size_t i, size_t j) const { return m_array[16*i+j]; }
    inline void set_uint16(size_t i, size_t j, uint16_t val) { m_array[16*i+j]=val; }

    inline uint16_t* data() { return m_array; }
    inline const uint16_t* data() const { return m_array; }

    inline size_t size() const { return N; }
private:
    alignas(32) uint16_t __restrict__ m_array[N*16];
};

#endif // WIB2_TPG_UTILITIES_HPP_
