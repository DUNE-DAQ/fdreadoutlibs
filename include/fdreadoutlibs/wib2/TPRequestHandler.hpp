/**
 * @file DefaultSkipListRequestHandler.hpp Trigger matching mechanism 
 * used for skip list based LBs in readout models
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPREQUESTHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPREQUESTHANDLER_HPP_

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"

#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"
#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/DefaultSkipListRequestHandler.hpp"

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"
#include "trigger/TPSet.hpp"
 
#include <atomic>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_WORK_STEPS;


namespace dunedaq {
ERS_DECLARE_ISSUE(fdreadoutlibs,
                  DroppedTPSet,
                  "Failed to send TPs from  " << s_ts << " to " << e_ts,
                  ((uint64_t)s_ts) 
		  ((uint64_t)e_ts))
ERS_DECLARE_ISSUE(fdreadoutlibs,
                   TPHandlerMsg,
                   infomsg,
                   ((std::string) infomsg))

namespace fdreadoutlibs {

//template<>	
class TPRequestHandler : public dunedaq::readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>
{
public:
  using inherited2 = readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>;
  using timestamp_t = std::uint64_t;

  // Constructor that binds LB and error registry

  TPRequestHandler(std::unique_ptr<readoutlibs::SkipListLatencyBufferModel<types::TriggerPrimitiveTypeAdapter>>& latency_buffer,
                                std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>(
        latency_buffer,
        error_registry), m_tp_set_sender_thread(0)
  {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "TPRequestHandler created...";
  }
 
  void init(const nlohmann::json& args) override;
  void conf(const nlohmann::json& args) override;
  void start(const nlohmann::json& args) override;
  void stop(const nlohmann::json& args) override;

protected:
  
private:
  std::shared_ptr<iomanager::SenderConcept<dunedaq::trigger::TPSet>> m_tpset_sink;
  int m_tp_set_sender_rate_hz;
  uint64_t m_ts_set_sender_offset_ticks;
  uint64_t m_run_number;
  dunedaq::readoutlibs::ReusableThread  m_tp_set_sender_thread;
  void send_tp_sets();
  uint64_t m_next_tpset_seqno;
};

} // namespace readoutlibs
} // namespace dunedaq


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPREQUESTHANDLER_HPP_
