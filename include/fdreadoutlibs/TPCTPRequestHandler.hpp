/**
 * @file TPCTPRequestHandler.hpp Trigger matching mechanism 
 * used for skip list based LBs in readout models
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_

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
class TPCTPRequestHandler : public dunedaq::readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>
{
public:
  using inherited2 = readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>;
  using timestamp_t = std::uint64_t;

  // Constructor that binds LB and error registry

  TPCTPRequestHandler(std::unique_ptr<readoutlibs::SkipListLatencyBufferModel<types::TriggerPrimitiveTypeAdapter>>& latency_buffer,
                                std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>(
        latency_buffer,
        error_registry), m_tp_set_sender_thread(0)
  {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "TPCTPRequestHandler created...";
  }
 
  void conf(const appdal::ReadoutModule* conf) override;
  void start(const nlohmann::json& args) override;
  void stop(const nlohmann::json& args) override;
  void get_info(opmonlib::InfoCollector& ci, int level) override;

protected:
  
private:
  std::shared_ptr<iomanager::SenderConcept<dunedaq::trigger::TPSet>> m_tpset_sink;
  int m_tp_set_sender_sleep_us;
  uint64_t m_ts_set_sender_offset_ticks;
  uint64_t m_run_number;
  dunedaq::readoutlibs::ReusableThread  m_tp_set_sender_thread;
  void send_tp_sets();
  uint64_t m_next_tpset_seqno;

  std::atomic<uint64_t> m_new_tps{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tpsets{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps_dropped{ 0 };
  std::atomic<uint64_t> m_new_heartbeats{ 0 };

};

} // namespace readoutlibs
} // namespace dunedaq


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_
