/**
 * @file TPEmulatorModel.hpp Emulates a tp source
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/wib2/TPEmulatorModel.hpp" // NOLINT(build/include)

// #ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPEMULATORMODEL_HPP_
// #define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_TPEMULATORMODEL_HPP_

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"
// #include "iomanager/Receiver.hpp"

// #include "logging/Logging.hpp"

// #include "opmonlib/InfoCollector.hpp"

// #include "readoutlibs/sourceemulatorconfig/Nljs.hpp"
// #include "readoutlibs/sourceemulatorinfo/InfoNljs.hpp"

// #include "readoutlibs/ReadoutIssues.hpp"
// #include "readoutlibs/concepts/SourceEmulatorConcept.hpp"
// #include "readoutlibs/utils/FileSourceBuffer.hpp"
// #include "readoutlibs/utils/RateLimiter.hpp"
// #include "readoutlibs/utils/ReusableThread.hpp"


// #include "fdreadoutlibs/DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter.hpp"

// #include "detdataformats/wib/RawWIBTp.hpp"

// #include <functional>
// #include <memory>
// #include <random>
// #include <string>
// #include <utility>
// #include <vector>

using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;
using dunedaq::readoutlibs::logging::TLVL_WORK_STEPS;

namespace dunedaq {
namespace fdreadoutlibs {

TPEmulatorModel::TPEmulatorModel(std::atomic<bool>& run_marker, double rate_khz)
  : m_run_marker(run_marker)
  , m_packet_count{ 0 }
  , m_sink_queue_timeout_ms(0)
  , m_raw_data_sink(nullptr)
  , m_producer_thread(0)
  , m_rate_khz(rate_khz)
{
}

void
TPEmulatorModel::init(const nlohmann::json& /*args*/)
{
}

void
TPEmulatorModel::set_sender(const std::string& sink_name)
{
  if (!m_sink_is_set) {
    m_raw_data_sink = get_iom_sender<types::DUNEWIBFirmwareTriggerPrimitiveSuperChunkTypeAdapter>(sink_name);
    m_sink_is_set = true;
  } else {
    // ers::error();
  }
}

void
TPEmulatorModel::conf(const nlohmann::json& args, const nlohmann::json& link_conf)
{
  if (m_is_configured) {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "This emulator is already configured!";
  } else {
    m_conf = args.get<module_conf_t>();
    m_link_conf = link_conf.get<link_conf_t>();
    m_sink_queue_timeout_ms = std::chrono::milliseconds(m_conf.queue_timeout_ms);

    m_sourceid.id = m_link_conf.source_id;
    m_sourceid.subsystem = daqdataformats::SourceID::Subsystem::kTrigger;

    m_file_source = std::make_unique<readoutlibs::FileSourceBuffer>(m_link_conf.input_limit, RAW_WIB2_TP_SUBFRAME_SIZE);

    try {
      m_file_source->read(m_link_conf.tp_data_filename);
    } catch (const ers::Issue& ex) {
      ers::fatal(ex);
      throw readoutlibs::ConfigurationError(ERS_HERE, m_sourceid, "", ex);
    }

    m_is_configured = true;
  }
  // Configure thread:
  m_producer_thread.set_name("fakeprod-tp", m_link_conf.source_id);
}

bool
TPEmulatorModel::is_configured()
{
  return m_is_configured;
}

void
TPEmulatorModel::scrap(const nlohmann::json& /*args*/)
{
  m_is_configured = false;
}

void
TPEmulatorModel::start(const nlohmann::json& /*args*/)
{
  m_packet_count_tot = 0;
  TLOG_DEBUG(TLVL_WORK_STEPS) << "Starting threads...";
  m_rate_limiter = std::make_unique<readoutlibs::RateLimiter>(m_rate_khz / m_link_conf.slowdown);
  m_producer_thread.set_work(&TPEmulatorModel::run_produce, this);
}

void
TPEmulatorModel::stop(const nlohmann::json& /*args*/)
{
  while (!m_producer_thread.get_readiness()) {
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
}

void
TPEmulatorModel::get_info(opmonlib::InfoCollector& ci, int /*level*/)
{
  readoutlibs::sourceemulatorinfo::Info info;
  info.packets = m_packet_count_tot.load();
  info.new_packets = m_packet_count.exchange(0);

  ci.add(info);
}

void
TPEmulatorModel::run_produce()
{
  TLOG_DEBUG(TLVL_WORK_STEPS) << "Data generation thread " << m_this_link_number << " started";

  int offset = 0;
  auto& source = m_file_source->get();
  int num_elem = m_file_source->num_elements(); // bytes

  if (num_elem == 0) {
    TLOG_DEBUG(TLVL_WORK_STEPS) << "No raw WIB2 TP elements to read from buffer! Sleeping...";
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    num_elem = m_file_source->num_elements();
  }
  TLOG_DEBUG(TLVL_WORK_STEPS) << "Raw WIB2 TP bytes to read from buffer: " << num_elem * static_cast<int>(RAW_WIB2_TP_SUBFRAME_SIZE);
  while (m_run_marker.load()) {
    // Which element to push to the buffer
    if (offset == num_elem * static_cast<int>(RAW_WIB2_TP_SUBFRAME_SIZE)) { // NOLINT(build/unsigned)
      offset = 0;
    }
    int bsize = num_elem * static_cast<int>(RAW_WIB2_TP_SUBFRAME_SIZE);
    std::vector<char> tmpbuffer;
    tmpbuffer.reserve(bsize);
    ::memcpy(static_cast<void*>(tmpbuffer.data()), static_cast<void*>(source.data() + offset), bsize);
    m_payload_wrapper.set_raw_tp_frame_chunk(tmpbuffer);

    offset += bsize;

    // queue in to actual iomanager::Sender
    try {
      m_raw_data_sink->send(std::move(m_payload_wrapper), m_sink_queue_timeout_ms);
    } catch (const dunedaq::iomanager::TimeoutExpired& excpt) {
      // std::runtime_error("Queue timed out...");
    }

    // Count packet and limit rate if needed.
    ++m_packet_count;
    ++m_packet_count_tot;
    m_rate_limiter->limit();
  }
  TLOG_DEBUG(TLVL_WORK_STEPS) << "Data generation thread " << m_this_link_number << " finished";
}

} // namespace fdreadoutlibs
} // namespace dunedaq
