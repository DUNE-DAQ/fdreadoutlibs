/**
 * @file TPEmulatorModel.hpp Emulates a tp source
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPEMULATORMODEL_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPEMULATORMODEL_HPP_

#include "iomanager/Sender.hpp"
#include "iomanager/Receiver.hpp"

#include "logging/Logging.hpp"

#include "opmonlib/InfoCollector.hpp"

#include "readoutlibs/sourceemulatorconfig/Nljs.hpp"
#include "readoutlibs/sourceemulatorinfo/InfoNljs.hpp"

#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/concepts/SourceEmulatorConcept.hpp"
#include "readoutlibs/utils/FileSourceBuffer.hpp"
#include "readoutlibs/utils/RateLimiter.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"


#include "fdreadoutlibs/DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter.hpp"

#include "detdataformats/wib/RawWIBTp.hpp"

#include <functional>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace dunedaq {
namespace fdreadoutlibs {

class TPEmulatorModel : public readoutlibs::SourceEmulatorConcept
{
public:
  using sink_t = iomanager::SenderConcept<detdataformats::wib::RawWIBTp>;


  // Raw WIB TP
  static const constexpr std::size_t RAW_WIB2_TP_SUBFRAME_SIZE = 12;
  // same size for header, tp data, pedinfo: 3 words * 4 bytes/word

  explicit TPEmulatorModel(std::atomic<bool>& run_marker, double rate_khz);

  void init(const nlohmann::json& /*args*/);

  void set_sender(const std::string& sink_name);

  void conf(const nlohmann::json& args, const nlohmann::json& link_conf);

  bool is_configured() override;

  void scrap(const nlohmann::json& /*args*/);

  void start(const nlohmann::json& /*args*/);

  void stop(const nlohmann::json& /*args*/);

  void get_info(opmonlib::InfoCollector& ci, int /*level*/);

protected:
  void run_produce();

private:
  // Constuctor params
  std::atomic<bool>& m_run_marker;

  // CONFIGURATION
  uint32_t m_this_apa_number;  // NOLINT(build/unsigned)
  uint32_t m_this_link_number; // NOLINT(build/unsigned)

  bool m_found_tp_header{ false };

  // STATS
  std::atomic<int> m_packet_count{ 0 };
  std::atomic<int> m_packet_count_tot{ 0 };

  readoutlibs::sourceemulatorconfig::Conf m_cfg;

  // RAW SINK
  std::chrono::milliseconds m_sink_queue_timeout_ms;
  using raw_sink_qt = iomanager::SenderConcept<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>;
  std::shared_ptr<raw_sink_qt> m_raw_data_sink;

  bool m_sink_is_set = false;
  using module_conf_t = dunedaq::readoutlibs::sourceemulatorconfig::Conf;
  module_conf_t m_conf;
  using link_conf_t = dunedaq::readoutlibs::sourceemulatorconfig::LinkConfiguration;
  link_conf_t m_link_conf;

  std::unique_ptr<readoutlibs::RateLimiter> m_rate_limiter;
  std::unique_ptr<readoutlibs::FileSourceBuffer> m_file_source;

  readoutlibs::ReusableThread m_producer_thread;

  bool m_is_configured = false;
  double m_rate_khz;
  daqdataformats::SourceID m_sourceid;

  types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter m_payload_wrapper;

};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPEMULATORMODEL_HPP_
