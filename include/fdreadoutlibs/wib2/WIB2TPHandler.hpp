/**
 * @file WIB2TPHandler.hpp Buffer for TPSets
 *
 * This is part of the DUNE DAQ , copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIB2TPHANDLER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIB2TPHANDLER_HPP_

#include "logging/Logging.hpp"
#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "trigger/TPSet.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include <queue>
#include <utility>
#include <vector>

namespace dunedaq {

ERS_DECLARE_ISSUE(fdreadoutlibs,
                   TPHandlerTimestampIssue,
                  "Continuity of timestamps broken. Start ts:  " << start_ts << " ts previous tpset: " << current_ts,
                   ((uint64_t)start_ts)((uint64_t)current_ts))


namespace fdreadoutlibs {

class WIB2TPHandler
{
public:
  explicit WIB2TPHandler(std::shared_ptr<iomanager::SenderConcept<types::TriggerPrimitiveTypeAdapter>> tp_sink,
                        std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> tpset_sink,
                        uint64_t tp_timeout,        // NOLINT(build/unsigned)
                        uint64_t tpset_window_size, // NOLINT(build/unsigned)
                        daqdataformats::SourceID sourceId);

  void set_run_number(daqdataformats::run_number_t run_number);

  daqdataformats::run_number_t get_run_number();
  
  bool add_tp(triggeralgs::TriggerPrimitive trigprim, uint64_t currentTime) ;

  void try_sending_tpsets(uint64_t currentTime);

  void reset();

  size_t get_and_reset_num_sent_tps();

  size_t get_and_reset_num_sent_tpsets();

private:
  std::shared_ptr<iomanager::SenderConcept<types::TriggerPrimitiveTypeAdapter>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> m_tpset_sink;
  daqdataformats::run_number_t m_run_number{ daqdataformats::TypeDefaults::s_invalid_run_number };
  uint64_t m_tp_timeout;           // NOLINT(build/unsigned)
  uint64_t m_tpset_window_size;    // NOLINT(build/unsigned)
  uint64_t m_next_tpset_seqno = 0; // NOLINT(build/unsigned)
  daqdataformats::SourceID m_sourceid; 
  uint64_t m_timestamp_counter = 0;

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

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB2_WIB2TPHANDLER_HPP_
