
/**
 * @file TPCTPReadoutModel.hpp Variant of the default readout model which introduces
 * a SkipList after data reception and before insertion into the latency buffer
 * to reorder TPs
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "readoutlibs/models/ReadoutModel.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "readoutlibs/SkipListLatencyBufferModel.hpp"
#include "readoutlibs/TaskRawDataProcessorModel.hpp"
#include "fdreadoutlibs/TPCTPRequestHandler.hpp"

#include <map>

namespace dunedaq {
namespace fdreadoutlibs {

class TPCTPReadoutModel : public ReadoutModel<types::TriggerPrimitiveTypeAdapter,
      TPCTPRequestHandler,
      readoutlibs::SkipListLatencyBufferModel<types::TriggerPrimitiveTypeAdapter>,
      readoutlibs::TaskRawDataProcessorModel<types::TriggerPrimitiveTypeAdapter>>
{
private:
    void push_reordered_tp(trgfdataformats::TriggerPrimitive tp);
    void run_consume();
    std::map<std::pair<trgdataformats::timestamp_t, trgdataformats::channel_t>, types::TriggerPrimitive> m_tps_map;
    
}
}} //namespace dunedaq::fdreadoutlibs