/**
 * @file TDEEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_

#include "fdreadoutlibs/TPCEthFrameProcessor.hpp"
#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

class TDEEthFrameProcessor : public TPCEthFrameProcessor<types::TDEEthTypeAdapter>
{
  public:
    using TPCEthFrameProcessor<types::TDEEthTypeAdapter>::TPCEthFrameProcessor;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TDEETHFRAMEPROCESSOR_HPP_
