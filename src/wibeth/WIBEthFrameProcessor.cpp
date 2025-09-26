/**
 * @file WIBEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp" // NOLINT(build/include)
#include "confmodel/GeoId.hpp"
#include "appmodel/RawDataProcessor.hpp"
#include "appmodel/TPCRawDataProcessor.hpp"
#include "appmodel/ProcessingStep.hpp"
#include "appmodel/SamplesOverThresholdMinima.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/IterableQueueModel.hpp"

#include  "datahandlinglibs/opmon/datahandling_info.pb.h"

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;
using dunedaq::datahandlinglibs::logging::TLVL_TAKE_NOTE;

// THIS SHOULDN'T BE HERE!!!!! But it is necessary.....
DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")
DUNE_DAQ_TYPESTRING(std::vector<dunedaq::trigger::TriggerPrimitiveTypeAdapter>, "TriggerPrimitiveVector")

namespace dunedaq {
namespace fdreadoutlibs {

WIBEthFrameProcessor::WIBEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled)
  : TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>(error_registry, processing_enabled)
{
}

void
WIBEthFrameProcessor::start(const appfwk::DAQModule::CommandData_t& args)
{
  // Reset software TPG resources
  if (m_post_processing_enabled) {
    m_tps_suppressed_too_long = 0;
    m_tps_send_failed = 0;
  } 

  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_ts_problem_reported = false;
  m_ts_error_state = false;
  m_ts_error_ctr = 0;

  m_first_seq_id_mismatch = true;
  m_seq_id_problem_reported = false;
  m_seq_id_error_state = false;
  m_seq_id_error_ctr = 0;


  // Reset stats
  m_t0 = std::chrono::high_resolution_clock::now();
  m_new_hits = 0;
  m_new_tps = 0;
  m_tpg_hits_count.exchange(0);
  inherited::start(args);
}

void
WIBEthFrameProcessor::stop(const appfwk::DAQModule::CommandData_t& args)
{
  inherited::stop(args);
  if (m_post_processing_enabled) {
    if (m_tpg_metric_collect_enabled) {
      m_tp_generator->free_metric_collector();
    }
    // Clears the pipelines and resets with the given configs.
    m_tp_generator->set_metric_collector_enable_state(m_tpg_metric_collect_enabled);
    m_tp_generator->configure(m_tpg_configs, m_channel_plane_numbers, types::DUNEWIBEthTypeAdapter::samples_tick_difference);
  }
}

void
WIBEthFrameProcessor::conf(const appmodel::DataHandlerModule* conf)
{
  size_t idx = 0;
  for (auto output : conf->get_outputs()) {
    try {
      if (output->get_data_type() == "TriggerPrimitiveVector") {
         m_tp_sink[idx++] = get_iom_sender<std::vector<trigger::TriggerPrimitiveTypeAdapter>>(output->UID());
      }
    } catch (const ers::Issue& excpt) {
      ers::error(datahandlinglibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
    }
  }

  m_sourceid.id = conf->get_source_id();
  m_sourceid.subsystem = types::DUNEWIBEthTypeAdapter::subsystem;
  auto geo_id = conf->get_geo_id();
  if (geo_id != nullptr) {
    m_det_id = geo_id->get_detector_id();
    m_crate_id = geo_id->get_crate_id();
    m_slot_id = geo_id->get_slot_id();
    m_stream_id = geo_id->get_stream_id();
  }
  m_emulator_mode = conf->get_emulation_mode();

  // Setup pre-processing pipeline
  if (!m_emulator_mode)
    inherited::add_preprocess_task(std::bind(&WIBEthFrameProcessor::sequence_check, this, std::placeholders::_1));

  inherited::add_preprocess_task(std::bind(&WIBEthFrameProcessor::timestamp_check, this, std::placeholders::_1));

  // Check it post-processing is active
  auto dp = conf->get_module_configuration()->get_data_processor();
  if (dp != nullptr) {
    auto proc_conf = dp->cast<appmodel::TPCRawDataProcessor>();
    if (proc_conf != nullptr && m_post_processing_enabled) {
      m_tp_generator = std::make_unique<tpglibs::TPGenerator>();

      // Set the number of frames and TPs above which TPs are sent to sink.
      m_tp_count_limit = proc_conf->get_tp_count_limit();
      m_frame_count_limit = proc_conf->get_frame_count_limit();

      // Set the minimum TP samples over threshold.
      auto conf_sot_minima = proc_conf->get_sot_minima();
      std::vector<uint16_t> sot_minima{conf_sot_minima->get_sot_minimum_plane0(),
                                       conf_sot_minima->get_sot_minimum_plane1(),
                                       conf_sot_minima->get_sot_minimum_plane2()};
      m_tp_generator->set_sot_minima(sot_minima);

      const std::vector<unsigned int> channel_mask_vec = proc_conf->get_channel_mask();

      std::vector<const appmodel::ProcessingStep*> processing_steps = proc_conf->get_processing_steps();
      for (auto step : processing_steps) {
        m_tpg_configs.push_back(std::make_pair(step->class_name(), step->to_json(false).back()));
      }

      // Setup post-processing pipeline
      m_channel_map = dunedaq::detchannelmaps::make_tpc_map(proc_conf->get_channel_map());
      for (int chan = 0; chan < 64; chan++) {
        trgdataformats::channel_t off_channel = m_channel_map->get_offline_channel_from_det_crate_slot_stream_chan(m_det_id, m_crate_id, m_slot_id, m_stream_id, chan);
        int16_t plane = m_channel_map->get_plane_from_offline_channel(off_channel);
        m_channel_plane_numbers.push_back(std::make_pair(off_channel, plane));

        // This processor only needs to handle some (maybe 0) of the masked channels.
        // Only get those relevant channels for the later check.
        if (std::find(channel_mask_vec.begin(), channel_mask_vec.end(), off_channel) != channel_mask_vec.end())
          m_channel_mask_set.insert(off_channel);
      }

      m_metric_collect_opmon_period = proc_conf->get_metric_collect_opmon_rate();

      // Let the TPG generator configure

      m_tp_generator->configure(m_tpg_configs, m_channel_plane_numbers, types::DUNEWIBEthTypeAdapter::samples_tick_difference);
      
      // After it sees the configs, it will set the metric collector enable state
      
      m_tpg_metric_collect_enabled = m_tp_generator->get_metric_collector_enable_state();

      inherited::add_postprocess_task(std::bind(&WIBEthFrameProcessor::find_hits, this, std::placeholders::_1));
    }
  }
  inherited::conf(conf);
}

void
WIBEthFrameProcessor::generate_opmon_data()
{
   datahandlinglibs::opmon::FixedRateDataProcessorInfo info;

   info.set_num_seq_id_errors(m_seq_id_error_ctr.load());
   info.set_min_seq_id_jump(m_seq_id_min_jump.exchange(0));
   info.set_max_seq_id_jump(m_seq_id_max_jump.exchange(0));

   info.set_num_ts_errors(m_ts_error_ctr.load());
   
   publish(std::move(info));

   m_error_registry->log_registered_errors();

   if (m_post_processing_enabled) {
     auto now = std::chrono::high_resolution_clock::now();
     int new_hits = m_tpg_hits_count.exchange(0);
     int new_tps = m_new_tps.exchange(0);
     int new_tps_suppressed_too_long = m_tps_suppressed_too_long.exchange(0);
     int new_tps_send_failed = m_tps_send_failed.exchange(0);
     double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
     TLOG_DEBUG(TLVL_BOOKKEEPING) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
     TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new hits: " << new_hits << " new TPs: " << new_tps;

     datahandlinglibs::opmon::HitFindingInfo tp_info;
     tp_info.set_rate_tp_hits(new_hits / seconds / 1000.);

     tp_info.set_num_tps_sent(new_tps);
     tp_info.set_num_tps_suppressed_too_long(new_tps_suppressed_too_long);
     tp_info.set_num_tps_send_failed(new_tps_send_failed);

     publish(std::move(tp_info));
     // Find the channels with the top  TP rates
     // Create a vector of pairs to store the map elements
     std::vector<std::pair<uint, int>> channel_tp_rate_vec(m_tp_channel_rate_map.begin(), m_tp_channel_rate_map.end());
     // Sort the vector in descending order of the value of the pairs
     sort(channel_tp_rate_vec.begin(), channel_tp_rate_vec.end(), [](std::pair<uint, int>& a, std::pair<uint, int>& b) { return a.second > b.second; });
     // Add the metrics to opmon
     // For convenience we are selecting only the top 10 elements
     if (channel_tp_rate_vec.size() != 0) {
       int top_highest_values = 10;
       if (channel_tp_rate_vec.size() < 10) {
         top_highest_values = channel_tp_rate_vec.size();
       }
       //datahandlinglibs::opmon::TPChannelsInfo channels_info;
       for (int i = 0; i < top_highest_values; i++) {
         datahandlinglibs::opmon::TPChannelInfo tpc_info;
         tpc_info.set_number_of_tps(channel_tp_rate_vec[i].second);
          tpc_info.set_channel_id(channel_tp_rate_vec[i].first);
         publish(std::move(tpc_info), {{"channel", std::to_string(channel_tp_rate_vec[i].first)}});
       }
     }

     // Reset the counter in the channel rate map
     for (auto& el : m_tp_channel_rate_map) {
       el.second = 0;
     }
     m_t0 = now;

     if (m_tpg_metric_collect_enabled && m_tp_generator) {
       publish_processor_metric_to_opmon();
       publish_processor_metric_to_opmon_with_aggregation();
     }
   }
   
   inherited::generate_opmon_data();
 }

 void
WIBEthFrameProcessor::publish_processor_metric_to_opmon() {
  auto metrics = m_tp_generator->get_processor_metrics();
  for (const auto& [channel, vec] : metrics) {
    datahandlinglibs::opmon::TPGProcessorInfo tpg_proc_info;
    for (const auto& [name, val] : vec) {
      if (name == "m_pedestal") {
        tpg_proc_info.set_pedestal(val);
      } else if (name == "m_accum") {
        tpg_proc_info.set_accum(val);
      }
    }
    publish(std::move(tpg_proc_info), {{"channel", std::to_string(channel)}});
  }
}

std::map<int16_t, std::map<std::string, std::tuple<float, int16_t, int16_t, float, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>> 
WIBEthFrameProcessor::calculate_all_metric_summaries_across_planes(const std::unordered_map<dunedaq::trgdataformats::channel_t, std::vector<std::pair<std::string, int16_t>>>& metrics) {
    // Structure to hold all statistics: plane -> metric -> (mean, min, max, stddev, min_channel_id, max_channel_id)
    std::map<int16_t, std::map<std::string, std::tuple<float, int16_t, int16_t, float, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>> all_stats;
    
    // Structure to accumulate statistics: plane -> metric -> (count, mean, M2, min, max, min_channel_id, max_channel_id)
    std::map<int16_t, std::map<std::string, std::tuple<size_t, double, double, int16_t, int16_t, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>> accumulators;
    
    // Single pass through all metrics to collect data using Welford's online algorithm for variance
    for (const auto& [channel, vec] : metrics) {
        if (!m_channel_map) continue;
        
        int16_t plane = m_channel_map->get_plane_from_offline_channel(channel);
        
        for (const auto& [name, val] : vec) {
            auto& [count, mean, M2, min, max, min_channel_id, max_channel_id] = accumulators[plane][name];
            
            count++;
            
            if (count == 1 || val < min) {
                min = val;
                min_channel_id = channel;
            }
            if (count == 1 || val > max) {
                max = val;
                max_channel_id = channel;
            }
            
            // Welford's online algorithm for variance calculation
            if (count == 1) {
                // First value: initialize mean and M2
                mean = val;
                M2 = 0.0;
            } else {
                // Update mean and M2 using Welford's algorithm
                double delta = val - mean;
                mean += delta / count;
                double delta2 = val - mean;
                M2 += delta * delta2;
            }
        }
    }
    
    // Calculate final statistics from accumulated data
    for (const auto& [plane, metric_map] : accumulators) {
        for (const auto& [metric_name, acc_data] : metric_map) {
            const auto& [count, mean, M2, min, max, min_channel_id, max_channel_id] = acc_data;
            
            if (count == 0) continue;
            
            float stddev = 0.0f;
            
            // Calculate standard deviation using accumulated M2
            if (count > 1) {
                stddev = std::sqrt(M2 / (count - 1));
            }
            
            all_stats[plane][metric_name] = std::make_tuple(static_cast<float>(mean), min, max, stddev, min_channel_id, max_channel_id);
        }
    }
    
    return all_stats;
}

void
WIBEthFrameProcessor::publish_processor_metric_to_opmon_with_aggregation() {
  auto metrics = m_tp_generator->get_processor_metrics();
  
  // Use optimized single-pass calculation for all metrics across all planes
  auto all_stats = calculate_all_metric_summaries_across_planes(metrics);
  
  // Publish all calculated statistics
  for (const auto& [plane, metric_map] : all_stats) {
    for (const auto& [metric_name, stats] : metric_map) {
      const auto& [mean, min, max, stddev, min_channel_id, max_channel_id] = stats;
      
      datahandlinglibs::opmon::TPGProcessorReducedInfo info;
      info.set_average(mean);
      info.set_max(max);
      info.set_min(min);
      info.set_standard_dev(stddev);
      info.set_max_channel_id(max_channel_id);
      info.set_min_channel_id(min_channel_id);
      publish(std::move(info), {{"plane", std::to_string(plane)}, {"metric", metric_name}});
    }
  }
}


/**
 * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
 * */
void
WIBEthFrameProcessor::sequence_check(frameptr fp)
{
  // FIXME: Make source emulator deal with this! Hard to do since source emu is templated...
  /* If EMU data, emulate perfectly incrementing timestamp
  if (m_emulator_mode) {  
    // uint64_t ts_next = m_previous_seq_id + 1; // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));            // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
      //auto wfh = const_cast<dunedaq::fddetdataformats::WIBEthFrame*>(wf->header());
      wf->daq_header.crate_id = m_crate_id;
      wf->daq_header.slot_id = m_slot_id;
      wf->daq_header.stream_id = m_stream_id; 
      wf->daq_header.seq_id = (m_previous_seq_id+i) & 0xfff;
      wf++;
    }
  }
  */
          
  // Acquire timestamp
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(fp); // NOLINT
  m_current_seq_id = wfptr->daq_header.seq_id;

  // Check sequence id
  // Calculate the next sequence id (12 bits)
  uint16_t expected_seq_id = (m_previous_seq_id + fp->get_num_frames()) & 0xfff;
  int16_t delta_seq_id = m_current_seq_id-expected_seq_id;
  if ( delta_seq_id > 0x800) {
    delta_seq_id -= 0x1000;
  } else if ( delta_seq_id < -0x7ff) {
    delta_seq_id += 0x1000;
  }

  if (delta_seq_id == 0) {
    m_seq_id_error_state = false;
  } else {
    // uint16_t delta_seq_id = (m_current_seq_id-expected_seq_id);
    ++m_seq_id_error_ctr;
    m_seq_id_max_jump = std::max(delta_seq_id, m_seq_id_max_jump.load());
    m_seq_id_min_jump = std::min(delta_seq_id, m_seq_id_min_jump.load());

    if (m_first_seq_id_mismatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First sequence id MISMATCH! -> | previous: " << std::to_string(m_previous_seq_id) << " current: " + std::to_string(m_current_seq_id);
      m_first_seq_id_mismatch = false;
    } else {
      if (!m_seq_id_error_state) {
        m_error_registry->add_error("Sequence ID jump", datahandlinglibs::FrameErrorRegistry::ErrorInterval(expected_seq_id, m_current_seq_id));
        m_seq_id_error_state = true;
      }
    }    
  }

  if (m_seq_id_error_ctr > 1000) {
    if (!m_seq_id_problem_reported) {
      TLOG() << "*** Data Integrity ERROR *** Sequence ID continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!";
      m_seq_id_problem_reported = true;
    }
  }

  m_previous_seq_id = m_current_seq_id;

}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
 * */
void
WIBEthFrameProcessor::timestamp_check(frameptr fp)
{

  uint16_t wibeth_tick_difference = types::DUNEWIBEthTypeAdapter::expected_tick_difference;
  uint16_t wibeth_frame_tick_difference = wibeth_tick_difference * fp->get_num_frames();

  // FIXME: let source emulator deal with this!
  /* If EMU data, emulate perfectly incrementing timestamp
  if (m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
    uint64_t ts_next = m_previous_ts + wibeth_frame_tick_difference; // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));            // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
      //auto wfh = const_cast<dunedaq::fddetdataformats::WIBEthFrame*>(wf->header());
      wf->daq_header.crate_id = m_crate_id;
      wf->daq_header.slot_id = m_slot_id;
      wf->daq_header.stream_id = m_stream_id; 
      wf->set_timestamp(ts_next);
      ts_next += wibeth_tick_difference;
      wf++;
    }
  }*/

  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(fp); // NOLINT
  m_current_ts = wfptr->get_timestamp();

  // Check timestamp
  if (m_previous_ts > 0 &&
      m_current_ts - m_previous_ts != wibeth_frame_tick_difference) [[unlikely]] {
    ++m_ts_error_ctr;
    if (m_first_ts_missmatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISMATCH! -> | previous: " << std::to_string(m_previous_ts) << " current: " + std::to_string(m_current_ts);
      m_first_ts_missmatch = false;
    } else {
      if (!m_ts_error_state) {
        m_error_registry->add_error("Timestamp jump", datahandlinglibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + wibeth_frame_tick_difference, m_current_ts));
        m_ts_error_state = true;
      }
    }
  } else {
    m_ts_error_state = false;
  }

  if (m_ts_error_ctr > 1000) {
    if (!m_ts_problem_reported) {
      TLOG() << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!";
      m_ts_problem_reported = true;
    }
  }

  m_previous_ts = m_current_ts;
  m_last_processed_daq_ts = m_current_ts;
}

/**
 * Pipeline Stage 2.: Do software TPG
 * */
void
WIBEthFrameProcessor::find_hits(constframeptr fp)
{
  size_t nhits = 0;
  if (!fp)
    return;
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp); // NOLINT

  // Check that the system is properly configured from the first hit.
  if (m_first_hit) {
    if (wfptr->daq_header.crate_id != m_crate_id || wfptr->daq_header.slot_id != m_slot_id || wfptr->daq_header.stream_id != m_stream_id) {
      ers::error(LinkMisconfiguration(ERS_HERE, wfptr->daq_header.crate_id, wfptr->daq_header.slot_id, wfptr->daq_header.stream_id, m_crate_id, m_slot_id, m_stream_id));
    }

    m_first_hit = false;
  }

  std::vector<trgdataformats::TriggerPrimitive> tps = (*m_tp_generator)(wfptr);
  m_current_frame_count++;
  if (m_tpg_metric_collect_enabled && m_frame_counter.load(std::memory_order_relaxed) % m_metric_collect_opmon_period == 0) {
    m_tp_generator->signal_metric_collection();
  }

  for (const auto& tp : tps) {
    // If this TP is on a masked channel, skip it.
    if (std::binary_search(m_channel_mask_set.begin(), m_channel_mask_set.end(), tp.channel))
      continue;
    m_current_tp_count++;
    // Need to move into a type adapter.
    trigger::TriggerPrimitiveTypeAdapter tpa;
    tpa.tp = tp;
    tpa.tp.detid = m_det_id;  // Last missing piece.
    m_tpa_vectors[m_channel_map->get_plane_from_offline_channel(tp.channel)].push_back(tpa);
    m_tp_channel_rate_map[tp.channel]++;
  }

  // if (m_current_frame_count > m_frame_count_limit || m_current_tp_count > m_tp_count_limit) {
  if (m_current_frame_count > m_frame_count_limit == 0 || m_current_tp_count > m_tp_count_limit) {
    for (int i = 0; i < 3; i++) {
      int new_tps = m_tpa_vectors[i].size();
      if (new_tps == 0) {
        continue;
      } 
      const auto s_ts_begin = m_tpa_vectors[i].front().tp.time_start;
      const auto channel_begin = m_tpa_vectors[i].front().tp.channel;
      const auto s_ts_end = m_tpa_vectors[i].back().tp.time_start;
      const auto channel_end = m_tpa_vectors[i].back().tp.channel;      
      if (!m_tp_sink[i]->try_send(std::move(m_tpa_vectors[i]), iomanager::Sender::s_no_block)) {
        ers::warning(FailedToSendTPVector(ERS_HERE, s_ts_begin, channel_begin, s_ts_end, channel_end));
        m_tps_send_failed++;
      } else {
        m_new_tps += new_tps;
        nhits += new_tps;
      }
    }
    m_current_tp_count=0;
    m_current_frame_count = 0;
  }
  m_tpg_hits_count += nhits;
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
