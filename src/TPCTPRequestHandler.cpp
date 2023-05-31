#include "fdreadoutlibs/TPCTPRequestHandler.hpp"
#include "appfwk/DAQModuleHelper.hpp"
#include "rcif/cmd/Nljs.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

void
TPCTPRequestHandler::init(const nlohmann::json& args) {
	inherited2::init(args);
   try {
       auto qi = dunedaq::appfwk::connection_index(args, {"tpset_out"});
       m_tpset_sink = iomanager::IOManager::get()->get_sender<dunedaq::trigger::TPSet>(qi["tpset_out"]);
   } catch (const ers::Issue& excpt) {
       throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
   }

}	
void
TPCTPRequestHandler::conf(const nlohmann::json& args) {
   auto conf = args["readoutmodelconf"].get<readoutlibs::readoutconfig::ReadoutModelConf>();
   m_tp_set_sender_thread.set_name("tpset", conf.source_id);
   inherited2::conf(args);
   m_tp_set_sender_sleep_us = 1000000/conf.tpset_transmission_rate_hz;
   m_ts_set_sender_offset_ticks = conf.tpset_min_latency_ticks;
}


void 
TPCTPRequestHandler::start(const nlohmann::json& args) {
   m_new_tps = 0;
   m_new_tpsets = 0;
   m_new_tps_dropped = 0;
   inherited2::start(args);

   rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
   m_run_number = start_params.run;

   m_tp_set_sender_thread.set_work(&TPCTPRequestHandler::send_tp_sets, this);
}

void
TPCTPRequestHandler::stop(const nlohmann::json& args) {
	m_run_marker.store(false);
	while(!m_tp_set_sender_thread.get_readiness()) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
	inherited2::stop(args);

}

void
TPCTPRequestHandler::get_info(opmonlib::InfoCollector& ci, int level)
{
  readoutlibs::readoutinfo::RawDataProcessorInfo info;

  auto now = std::chrono::high_resolution_clock::now();
  int new_tps = m_new_tps.exchange(0);
  int new_tpsets = m_new_tpsets.exchange(0);
  int new_tps_dropped = m_new_tps_dropped.exchange(0);
  int new_heartbeats = m_new_heartbeats.exchange(0);
  double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
  TLOG() << "TPSets rate: " << std::to_string(new_tpsets / seconds) << " [Hz], TP rate: " << std::to_string(new_tps / seconds) << ", heartbeats: " << std::to_string(new_heartbeats / seconds) << " [Hz]";
  //info.rate_tp_hits = new_hits / seconds / 1000.;
 
  info.num_tps_sent = new_tps;
  info.num_tpsets_sent = new_tpsets;
  info.num_tps_dropped = new_tps_dropped;
  info.num_heartbeats = new_heartbeats;
  m_t0 = now;
  inherited2::get_info(ci, level);
  ci.add(info);
}


void
TPCTPRequestHandler::send_tp_sets() {
   timestamp_t oldest_ts=0;
   timestamp_t newest_ts=0;
   timestamp_t start_win_ts=0;
   timestamp_t end_win_ts=0;
   bool first_cycle = true;
   dunedaq::dfmessages::DataRequest dr;

   while (m_run_marker.load()) {
   {
      std::unique_lock<std::mutex> lock(m_cv_mutex);
      m_cv.wait(lock, [&] { return !m_cleanup_requested; });
      m_requests_running++;
    }
    m_cv.notify_all();
    if(m_latency_buffer->occupancy() != 0) {
       // Prepare response
       RequestResult rres(ResultCode::kUnknown, dr);
       std::vector<std::pair<void*, size_t>> frag_pieces;

       // Get the newest TP
       SkipListAcc acc(inherited2::m_latency_buffer->get_skip_list());
       auto tail = acc.last();
       auto head = acc.first();
       newest_ts = (*tail).get_first_timestamp();
       oldest_ts = (*head).get_first_timestamp();
       
       if (first_cycle) {
    	  start_win_ts = oldest_ts;
	  first_cycle = false;
       }
       if (newest_ts - start_win_ts > m_ts_set_sender_offset_ticks) {
         end_win_ts = newest_ts - m_ts_set_sender_offset_ticks;
         frag_pieces = get_fragment_pieces(start_win_ts, end_win_ts, rres);
         auto num_tps = frag_pieces.size();
         /*
	 if (num_tps == 0) {
	       std::stringstream s;
               s << "No TPs in time interval " << start_win_ts << " " << end_win_ts;
	       ers::info(TPHandlerMsg(ERS_HERE, s.str()));
         }
	 
         else {
	 */
	 
         trigger::TPSet tpset;
         tpset.run_number = m_run_number;
         tpset.type = num_tps>0 ? trigger::TPSet::Type::kPayload : trigger::TPSet::Type::kHeartbeat;
         tpset.origin = m_sourceid;
         tpset.start_time = start_win_ts; // provisory timestamp, will be filled with first TP
         tpset.end_time = end_win_ts; // provisory timestamp, will be filled with last TP
         tpset.seqno = m_next_tpset_seqno++; // NOLINT(runtime/increment_decrement)
       // reserve the space for efficiency
         if (num_tps > 0) {    
            tpset.objects.reserve(frag_pieces.size());
            bool first_tp = true;
            for( auto f : frag_pieces) {
               trgdataformats::TriggerPrimitive tp = *(static_cast<trgdataformats::TriggerPrimitive*>(f.first));
	       if(first_tp) {
		       tpset.start_time = tp.time_start;
		       first_tp = false;
	       }
	       tpset.end_time = tp.time_start;
	       tpset.objects.emplace_back(std::move(tp)); 
            }
	 } 
         if(!m_tpset_sink->try_send(std::move(tpset), iomanager::Sender::s_no_block)) {
	    ers::warning(DroppedTPSet(ERS_HERE, start_win_ts, end_win_ts));
	    m_new_tps_dropped += num_tps;
         }
         m_new_tps += num_tps;
         m_new_tpsets++;

	 if (num_tps == 0) {
		 m_new_heartbeats++;
	 }
         //remember what we sent for the next loop
         start_win_ts = end_win_ts;
	 
       }
    }
    {
      std::lock_guard<std::mutex> lock(m_cv_mutex);
      m_requests_running--;
    }
    m_cv.notify_all();  
    std::this_thread::sleep_for(std::chrono::microseconds(m_tp_set_sender_sleep_us));
   }
   return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
