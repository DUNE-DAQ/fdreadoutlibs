/**
 * @file WIBTPHandler.hpp Buffer for TPSets
 *
 * This is part of the DUNE DAQ , copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBTPHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBTPHANDLER_HPP_

#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "trigger/TPSet.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include <queue>
#include <utility>
#include <vector>

namespace dunedaq {
namespace fdreadoutlibs {

class WIBTPHandler
{
public:
  explicit WIBTPHandler(iomanager::SenderConcept<types::SW_WIB_TRIGGERPRIMITIVE_STRUCT>& tp_sink,
                        iomanager::SenderConcept<trigger::TPSet>& tpset_sink,
                        uint64_t tp_timeout,        // NOLINT(build/unsigned)
                        uint64_t tpset_window_size, // NOLINT(build/unsigned)
                        daqdataformats::GeoID geoId,
                        std::string tpset_topic)
    : m_tp_sink(tp_sink)
    , m_tpset_sink(tpset_sink)
    , m_tp_timeout(tp_timeout)
    , m_tpset_window_size(tpset_window_size)
    , m_geoid(geoId)
    , m_tpset_topic(tpset_topic)
  {}

  void set_run_number(daqdataformats::run_number_t run_number)
  {
    m_run_number = run_number;
  }

  daqdataformats::run_number_t get_run_number()
  {
    return m_run_number;
  }
  
  bool add_tp(triggeralgs::TriggerPrimitive trigprim, uint64_t currentTime) // NOLINT(build/unsigned)
  {
    if (trigprim.time_start + m_tp_timeout > currentTime) {
      m_tp_buffer.push(trigprim);
      return true;
    } else {
      return false;
    }
  }

  void try_sending_tpsets(uint64_t currentTime) // NOLINT(build/unsigned)
  {
    if (!m_tp_buffer.empty() && m_tp_buffer.top().time_start + m_tpset_window_size + m_tp_timeout < currentTime) {
      trigger::TPSet tpset;
      tpset.run_number = m_run_number;
      tpset.start_time = (m_tp_buffer.top().time_start / m_tpset_window_size) * m_tpset_window_size;
      tpset.end_time = tpset.start_time + m_tpset_window_size;
      tpset.seqno = m_next_tpset_seqno++; // NOLINT(runtime/increment_decrement)
      tpset.type = trigger::TPSet::Type::kPayload;
      tpset.origin = m_geoid;
      
      while (!m_tp_buffer.empty() && m_tp_buffer.top().time_start < tpset.end_time) {
        triggeralgs::TriggerPrimitive tp = m_tp_buffer.top();
        types::SW_WIB_TRIGGERPRIMITIVE_STRUCT* tp_readout_type =
          reinterpret_cast<types::SW_WIB_TRIGGERPRIMITIVE_STRUCT*>(&tp); // NOLINT
        try {
            types::SW_WIB_TRIGGERPRIMITIVE_STRUCT tp_copy(*tp_readout_type);
          m_tp_sink.send(std::move(tp_copy), std::chrono::milliseconds(10));
          m_sent_tps++;
        } catch (const dunedaq::iomanager::TimeoutExpired& excpt) {
          ers::error(readoutlibs::CannotWriteToQueue(ERS_HERE, m_geoid, "m_tp_sink"));
        }
        tpset.objects.emplace_back(std::move(tp));
        m_tp_buffer.pop();
      }

      try {
        m_tpset_sink.send(std::move(tpset), std::chrono::milliseconds(10), m_tpset_topic);
        m_sent_tpsets++;
      } catch (const dunedaq::iomanager::TimeoutExpired& excpt) {
        ers::error(readoutlibs::CannotWriteToQueue(ERS_HERE, m_geoid, "m_tpset_sink"));
      }
    }
  }

  void reset()
  {
    while (!m_tp_buffer.empty()) {
      m_tp_buffer.pop();
    }
    m_next_tpset_seqno = 0;
    m_sent_tps = 0;
    m_sent_tpsets = 0;
  }

  size_t get_and_reset_num_sent_tps() { return m_sent_tps.exchange(0); }

  size_t get_and_reset_num_sent_tpsets() { return m_sent_tpsets.exchange(0); }

private:
  iomanager::SenderConcept<types::SW_WIB_TRIGGERPRIMITIVE_STRUCT>& m_tp_sink;
  iomanager::SenderConcept<trigger::TPSet>& m_tpset_sink;
  daqdataformats::run_number_t m_run_number{ daqdataformats::TypeDefaults::s_invalid_run_number };
  uint64_t m_tp_timeout;           // NOLINT(build/unsigned)
  uint64_t m_tpset_window_size;    // NOLINT(build/unsigned)
  uint64_t m_next_tpset_seqno = 0; // NOLINT(build/unsigned)
  daqdataformats::GeoID m_geoid;
  std::string m_tpset_topic;
  
  std::atomic<size_t> m_sent_tps{ 0 };    // NOLINT(build/unsigned)
  std::atomic<size_t> m_sent_tpsets{ 0 }; // NOLINT(build/unsigned)

  class TPComparator
  {
  public:
    bool operator()(triggeralgs::TriggerPrimitive& left, triggeralgs::TriggerPrimitive& right)
    {
      return left.time_start > right.time_start;
    }
  };
  std::priority_queue<triggeralgs::TriggerPrimitive, std::vector<triggeralgs::TriggerPrimitive>, TPComparator>
    m_tp_buffer;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBTPHANDLER_HPP_
