/**
 * @file FDReadoutTypes.hpp Common types in udaq-readout
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_VARIABLESIZEPAYLOADTYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_VARIABLESIZEPAYLOADTYPEADAPTER_HPP_

#include <cstdint> // uint_t types
#include <cstring> // memcpy
#include <memory>  // unique_ptr
#include <tuple>   // tie
#include <vector>

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

/**
 * @brief Convencience wrapper to take ownership over char pointers with
 * corresponding allocated memory size.
 * */
struct VariableSizePayloadTypeAdapter
{
  VariableSizePayloadTypeAdapter() {}
  VariableSizePayloadTypeAdapter(size_t size, char* data)
    : size(size)
    , data(data)
  {
  }

  size_t size = 0;
  std::unique_ptr<char> data = nullptr;
};

} // namespace types
} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_VARIABLESIZEPAYLOADTYPEADAPTER_HPP_
