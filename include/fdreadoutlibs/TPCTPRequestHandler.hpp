/**
 * @file DefaultSkipListRequestHandler.hpp Trigger matching mechanism 
 * used for skip list based LBs in readout models
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_

#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "daqdataformats/Types.hpp"
#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"
#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/DefaultSkipListRequestHandler.hpp"
#include "trigger/TPSet.hpp"
 
#include <atomic>
#include <chrono>
#include <memory>
#include <string>

using dunedaq::readoutlibs::logging::TLVL_WORK_STEPS;


namespace dunedaq {
ERS_DECLARE_ISSUE(fdreadoutlibs,
                  FailedToSendTPSet,
                  "Failed to send TPs from start time " << s_ts << " to end time " << e_ts << " in run " << run,
                  ((dunedaq::daqdataformats::timestamp_t)s_ts) 
                  ((dunedaq::daqdataformats::timestamp_t)e_ts)
                  ((dunedaq::daqdataformats::run_number_t)run))
ERS_DECLARE_ISSUE(fdreadoutlibs,
                   TPHandlerMsg,
                   infomsg,
                   ((std::string) infomsg))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  DataPacketArrivedTooLate,
                  "SourceID[" << sid_subsystem << "," << sid_id << "] Received late data packet, "
                  << "TP channel " << tpchan << ", "
                  << msec_diff << " ms beyond the allowed latency ("
                  << msec_allowed << " ms, " << ticks_allowed << " ticks).",
                  ((std::string)sid_subsystem)((daqdataformats::SourceID::ID_t)sid_id)
                  ((trgdataformats::channel_t)tpchan)((double)msec_diff)
                  ((double)msec_allowed)((daqdataformats::timestamp_t)ticks_allowed))

namespace fdreadoutlibs {

//template<>	
class TPCTPRequestHandler : public dunedaq::readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>
{
public:
  using inherited2 = readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>;
  using timestamp_t = dunedaq::daqdataformats::timestamp_t;

  // Constructor that binds LB and error registry

  TPCTPRequestHandler(std::unique_ptr<readoutlibs::SkipListLatencyBufferModel<types::TriggerPrimitiveTypeAdapter>>& latency_buffer,
                                std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : readoutlibs::DefaultSkipListRequestHandler<types::TriggerPrimitiveTypeAdapter>(
        latency_buffer,
        error_registry), m_tp_set_sender_thread(0)
  {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "TPCTPRequestHandler created...";
  }
 
  void init(const nlohmann::json& args) override;
  void conf(const nlohmann::json& args) override;
  void start(const nlohmann::json& args) override;
  void stop(const nlohmann::json& args) override;
  void get_info(opmonlib::InfoCollector& ci, int level) override;

  timestamp_t get_cutoff_timestamp() override {return m_cutoff_timestamp.load();}
  bool supports_cutoff_timestamp() override {return true;}
  void report_tardy_packet(const types::TriggerPrimitiveTypeAdapter& packet, int64_t tardy_ticks) override;

protected:
  
private:
  std::shared_ptr<iomanager::SenderConcept<dunedaq::trigger::TPSet>> m_tpset_sink;
  int m_tp_set_sender_sleep_us;
  uint64_t m_ts_set_sender_offset_ticks;
  int m_tardy_tp_quiet_time_at_start_sec;
  dunedaq::daqdataformats::run_number_t m_run_number;
  dunedaq::readoutlibs::ReusableThread  m_tp_set_sender_thread;
  void send_tp_sets();
  uint64_t m_next_tpset_seqno;

  std::atomic<uint64_t> m_new_tps{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tpsets{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps_in_tpsets_send_failed{ 0 };
  std::atomic<uint64_t> m_new_tpsets_send_failed{ 0 };
  std::atomic<uint64_t> m_new_tps_suppressed_tardy{ 0 };
  std::atomic<uint64_t> m_new_heartbeats{ 0 };

  std::atomic<timestamp_t> m_cutoff_timestamp{ 0 }; // NOLINT(build/unsigned)
  std::chrono::time_point<std::chrono::high_resolution_clock> m_run_start_timepoint;
};

} // namespace readoutlibs
} // namespace dunedaq


#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPCTPREQUESTHANDLER_HPP_
