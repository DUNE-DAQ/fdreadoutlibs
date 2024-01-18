/**
 * @file TPCTPProcessor.hpp TPC TP specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/TPCTPProcessor.hpp" // NOLINT(build/include)

#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/IterableQueueModel.hpp"
#include "readoutlibs/readoutconfig/Nljs.hpp"
#include "readoutlibs/readoutinfo/InfoNljs.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "fdreadoutlibs/TriggerActivityTypeAdapter.hpp"

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

// THIS SHOULDN'T BE HERE!!!!! But it is necessary.....
DUNE_DAQ_TYPESTRING(dunedaq::fdreadoutlibs::types::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")
DUNE_DAQ_TYPESTRING(dunedaq::fdreadoutlibs::types::TriggerActivityTypeAdapter, "TriggerActivity")

namespace dunedaq {
namespace fdreadoutlibs {

TPCTPProcessor::TPCTPProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
  : TaskRawDataProcessorModel<types::TriggerPrimitiveTypeAdapter>(error_registry)
{
}

TPCTPProcessor::~TPCTPProcessor()
{}

void
TPCTPProcessor::start(const nlohmann::json& args)
{

  // Reset stats
  m_new_hits = 0;
  m_new_tps = 0;
  inherited::start(args);
}

void
TPCTPProcessor::stop(const nlohmann::json& args)
{
  inherited::stop(args);
}

void
TPCTPProcessor::init(const nlohmann::json& args)
{
//  inherited::init(args);

  try {
    auto queue_index = appfwk::connection_index(args, {});
    if (queue_index.find("ta_out") != queue_index.end()) {
      m_tp_sink = get_iom_sender<types::TriggerActivityTypeAdapter>(queue_index["ta_out"]);
    }
  } catch (const ers::Issue& excpt) {
    ers::error(readoutlibs::ResourceQueueError(ERS_HERE, "ta", "TPCTPProcessor", excpt));
  }

}

void
TPCTPProcessor::conf(const nlohmann::json& cfg)
{
  auto config = cfg["tpprocessorconf"].get<readoutlibs::readoutconfig::TPProcessorConf>(); //FIXME: put correct schema

  m_sourceid.id = config.source_id;
  m_sourceid.subsystem = types::TriggerPrimitiveTypeAdapter::subsystem;

  auto ta_algorithms = config.ta_algorithms;

// FIXME: instantiate here the different TA algorithms: cannot use a list of basic TAMakers
// FIXME: add each algo as a post processor task

  for (auto algo_name in ta_algorithms)  {
    TLOG() << "Selected TA algorithm: " << algo_name;
    if (algo_name == "Prescale") {
      auto taa = std::make_shared<triggeralgs::TriggerActivityMakerPrescale>();
      taa->config(cfg);
      inherited::add_postprocess_task(std::bind(&TPCProcessor::find_ta, this, std::placeholders::_1, taa);
      m_taas.push_back(taa);
    }
    else if (algo_name == "PlaneCoincidence") {
      auto taa = std::make_shared<triggeralgs::TriggerActivityMakerPlaneCoincidence>();
      taa->config(cfg);
      inherited::add_postprocess_task(std::bind(&TPCProcessor::find_ta, this, std::placeholders::_1, taa);
      m_taas.push_back(taa);
    }
    else if (algo_name == "HorizontalMuon") {
      auto taa = std::make_shared<triggeralgs::TriggerActivityMakerHorizontalMuon>();
      taa->config(cfg);
      inherited::add_postprocess_task(std::bind(&TPCProcessor::find_ta, this, std::placeholders::_1,taa);
      m_taas.push_back(taa);
    }
    else if (algo_name == "DBSCAN") {
      auto taa = std::make_shared<triggeralgs::TriggerActivityMakerDBSCAN>();
      taa->config(cfg);
      inherited::add_postprocess_task(std::bind(&TPCProcessor::find_ta, this, std::placeholders::_1,taa);
      m_taas.push_back(taa);
    }
    else {
      ers::error(UnknownAlgorithm(ERS_HERE, algo_name));
    }
  }
  inherited::conf(cfg);
}

void
TPCTPProcessor::get_info(opmonlib::InfoCollector& ci, int level)
{

  inherited::get_info(ci, level);
  ci.add(info);
}


/**
 * Pipeline Stage 2.: Do software TPG
 * */
void
TPCTPProcessor::find_ta(types::TriggerPrimitiveTypeAdapter& tp,  std::shared_ptr<triggeralgs::TriggerActivityMaker> taa)
{
  std::vector<types::TriggerActivityTypeAdapter> tas;
  taa->operator()(tp, tas);
  for (auto ta in tas) {
    if(!m_ta_sink->try_send(ta, iomanager::Sender::s_no_block)) {
                 ers::warning(TADropped(ERS_HERE, ta.time_start, m_sourceid);
                   m_tas_dropped++;
            }
    m_new_tas++;
  }
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq