#include "fdreadoutlibs/TPCTARequestHandler.hpp"
#include "appdal/ReadoutModuleConf.hpp"
#include "appdal/RequestHandler.hpp"
#include "appdal/TriggerRequestHandler.hpp"

//#include "appfwk/DAQModuleHelper.hpp"
#include "rcif/cmd/Nljs.hpp"

namespace dunedaq {
namespace fdreadoutlibs {

void
TPCTARequestHandler::conf(const appdal::ReadoutModule* conf) {

   for (auto output : conf->get_outputs()) {
      if (output->get_data_type() == "TASet") {
         try {
            m_taset_sink = iomanager::IOManager::get()->get_sender<dunedaq::trigger::TASet>(output->UID());
         } catch (const ers::Issue& excpt) {
            throw readoutlibs::ResourceQueueError(ERS_HERE, "ta set out", "TPCTARequestHandlerModel", excpt);
         }
      }
   }
   m_ta_set_sender_thread.set_name("taset", conf->get_source_id());
  
   auto tah_conf = conf->get_module_configuration()->get_request_handler()->cast<appdal::TriggerRequestHandler>();
   if (tah_conf == nullptr) {
      throw readoutlibs::GenericConfigurationError(ERS_HERE, "The request handler for the TAHandlerModule is not of the right class.");
   }
   else {
      m_ta_set_sender_sleep_us = 1000000/tah_conf->get_max_transmission_rate_hz();
      //m_ts_set_sender_offset_ticks = tph_conf->get_min_latency_ticks();
   }
   inherited2::conf(conf);
}

void 
TPCTARequestHandler::start(const nlohmann::json& args) {
   m_new_tas = 0;
   m_new_tasets = 0;
   m_new_tas_dropped = 0;
   inherited2::start(args);

   rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
   m_run_number = start_params.run;

   m_ta_set_sender_thread.set_work(&TPCTARequestHandler::send_ta_sets, this);
}

void
TPCTARequestHandler::stop(const nlohmann::json& args) {
	m_run_marker.store(false);
	while(!m_ta_set_sender_thread.get_readiness()) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
	inherited2::stop(args);

}

void
TPCTARequestHandler::get_info(opmonlib::InfoCollector& ci, int level)
{
  readoutlibs::readoutinfo::RawDataProcessorInfo info;

  auto now = std::chrono::high_resolution_clock::now();
  int new_tas = m_new_tps.exchange(0);
  int new_tasets = m_new_tpsets.exchange(0);
  int new_tas_dropped = m_new_tps_dropped.exchange(0);
  int new_heartbeats = m_new_heartbeats.exchange(0);
  //double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
  //TLOG() << "TPSets rate: " << std::to_string(new_tpsets / seconds) << " [Hz], TP rate: " << std::to_string(new_tps / seconds) << ", heartbeats: " << std::to_string(new_heartbeats / seconds) << " [Hz]";
  //info.rate_tp_hits = new_hits / seconds / 1000.;
 
  info.num_tas_sent = new_tps;
  info.num_tasets_sent = new_tpsets;
  info.num_tas_dropped = new_tps_dropped;
  info.num_heartbeats = new_heartbeats;
  m_t0 = now;
  inherited2::get_info(ci, level);
  ci.add(info);
}


void
TPCTARequestHandler::send_ta_sets() {
   
   timestamp_t newest_ts=0;
   timestamp_t start_win_ts=0;
   
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
         //SkipListAcc acc(inherited2::m_latency_buffer->get_skip_list());
         
         auto tail = inherited2::m_latency_buffer->back();
         newest_ts = tail->get_first_timestamp();
         
         if (first_cycle) {
         auto head = inherited2::m_latency_buffer->front(); 
         start_win_ts = head->get_first_timestamp(); 
         first_cycle = false;
         }
         
         
         frag_pieces = get_fragment_pieces(start_win_ts, newest_ts, rres);
         auto num_tas = frag_pieces.size();
         trigger::TASet taset;
         taset.run_number = m_run_number;
         taset.type = num_tas>0 ? trigger::TASet::Type::kPayload : trigger::TASet::Type::kHeartbeat;
         taset.origin = m_sourceid;
         taset.start_time = start_win_ts; // provisory timestamp, will be filled with first TP
         taset.end_time = end_win_ts; // provisory timestamp, will be filled with last TP
         taset.seqno = m_next_tpset_seqno++; // NOLINT(runtime/increment_decrement)
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
