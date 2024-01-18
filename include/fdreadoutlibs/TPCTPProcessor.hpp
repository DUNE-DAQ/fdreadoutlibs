/**
 * @file TPCTPProcessor.hpp TPC TP specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCTPPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCTPPROCESSOR_HPP_

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "fdreadoutlibs/FDReadoutIssues.hpp"

#include "daqdataformats/Types.hpp"


namespace dunedaq {
namespace fdreadoutlibs {

class TPCTPProcessor : public readoutlibs::TaskRawDataProcessorModel<types::TriggerPrimitiveTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::TriggerPrimitiveTypeAdapter>;
  using tpptr = types::TriggerPrimitiveTypeAdapter*;
  using consttpptr = const types::TriggerPrimitiveTypeAdapter*;


  explicit TPCTPProcessor();

  ~TPCTPProcessor();

  void start(const nlohmann::json& args) override;

  void stop(const nlohmann::json& args) override;

  void init(const nlohmann::json& args) override;

  void conf(const nlohmann::json& cfg) override;

  void get_info(opmonlib::InfoCollector& ci, int level) override;

protected:
  // Internals
  dunedaq::daqdataformats::timestamp_t m_previous_ts = 0;
  dunedaq::daqdataformats::timestamp_t m_current_ts = 0;

  /**
   * Pipeline Stage 2.: Do TA finding
   * */

  void find_ta(consttpptr fp);

private:

  std::vector<std::shared_ptr<triggeralgs::TriggerActivityMaker>> m_taas;
  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;


  std::shared_ptr<iomanager::SenderConcept<fdreadoutlibs::types::TriggerActivityTypeAdapter>> m_ta_sink;

  //std::thread m_add_hits_tphandler_thread;

  daqdataformats::SourceID m_sourceid;

  std::atomic<uint64_t> m_new_tas{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tas_dropped{ 0 };
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPCTPPROCESSOR_HPP_
