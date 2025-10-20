/**
 * @file WIBEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_

#include "fdreadoutlibs/TPCEthFrameProcessor.hpp"
#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

class WIBEthFrameProcessor : public TPCEthFrameProcessor<types::DUNEWIBEthTypeAdapter>
{
  public:
    using TPCEthFrameProcessor<types::DUNEWIBEthTypeAdapter>::TPCEthFrameProcessor;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIBEth_WIBFRAMEPROCESSOR_HPP_
