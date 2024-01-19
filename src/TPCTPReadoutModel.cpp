#include "fdreadoutlibs/TPCTPReadoutModel.hpp"

namespace dunedaq {
    namespace fdreadoutlibs {

void TPCTPReadoutModel::push_reordered_tp(types::TriggerPrimitiveTypeAdapter payload)
{

  m_raw_processor_impl->preprocess_item(&payload);
  if (!m_latency_buffer_impl->write(std::move(payload))) {
    TLOG_DEBUG(TLVL_TAKE_NOTE) << "***ERROR: Latency buffer is full and data was overwritten!";
    m_num_payloads_overwritten++;
  }
  m_raw_processor_impl->postprocess_item(m_latency_buffer_impl->back());
  ++m_num_payloads;
  ++m_sum_payloads;
  ++m_stats_packet_count;
}

void 
TPCTPReadoutModel::run_consume()
{
  m_rawq_timeout_count = 0;
  m_num_payloads = 0;
  m_sum_payloads = 0;
  m_stats_packet_count = 0;

  trgdataformats::timestamp_t oldest_ts = 0;
  trgdataformats::timestamp_t newest_ts = 0;
  trgdataformats::timestamp_diff_t offset = 312500; // FIXME: we probably want to parametrize this
  bool first_cycle = true;

  TLOG() << "TPCTP Consumer thread started...";
  while (m_run_marker.load()) {
    // Try to acquire data
    auto  opt_payload = m_raw_data_receiver->try_receive(m_raw_receiver_timeout_ms);
    // Here put TP in a containter to reorder them and let them age... can be a map since we are in a single threaded environment
    if (opt_payload) {
	types::TriggerPrimitiveTypeAdapter& payload = opt_payload.value();
      if (first_cycle) {
        first_cycle = false;
      }
      m_tps_map[std::pair<trgdataformats::timestamp_t, trgdataformats::channel_t>
      (payload.get_timestamp(), payload.tp.channel)] = payload;

    } else {
      ++m_rawq_timeout_count;
      // Protection against a zero sleep becoming a yield
      if ( m_raw_receiver_sleep_us != std::chrono::microseconds::zero())
        std::this_thread::sleep_for(m_raw_receiver_sleep_us);
    }
    if (!m_tps_map.empty()) {
        oldest_ts = m_tps_map.begin()->first.first;
        newest_ts = m_tps_map.rbegin()->first.first;
        if (newest_ts - offset > oldest_ts) {
          auto itlow = m_tps_map.lower_bound(std::pair<trgdataformats::timestamp_t, trgdataformats::channel_t>
                                            (newest_ts - offset, trgdataformats::INVALID_CHANNEL));
          for(auto it = m_tps_map.begin(); it!= itlow; ++it) {
            push_reordered_tp(it->second);
          }
          m_tps_map.erase(m_tps_map.begin(), itlow);
        }   
    }
  }
  TLOG_DEBUG(TLVL_WORK_STEPS) << "Consumer thread joins... ";
}

}}
