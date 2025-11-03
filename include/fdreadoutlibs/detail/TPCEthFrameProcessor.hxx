DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")
DUNE_DAQ_TYPESTRING(std::vector<dunedaq::trigger::TriggerPrimitiveTypeAdapter>, "TriggerPrimitiveVector")

namespace dunedaq {
namespace fdreadoutlibs {

using datahandlinglibs::logging::TLVL_BOOKKEEPING;
using datahandlinglibs::logging::TLVL_TAKE_NOTE;

template <class ReadoutTypeAdapter>
TPCEthFrameProcessor<ReadoutTypeAdapter>::TPCEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool processing_enabled)
  : datahandlinglibs::TaskRawDataProcessorModel<ReadoutTypeAdapter>(error_registry, processing_enabled)
{
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::start(const appfwk::DAQModule::CommandData_t& args)
{
  // Reset software TPG resources
  if (this->m_post_processing_enabled) {
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
  m_num_new_tps.exchange(0);

      
  // Start the state harvester collection thread if enabled
  if (m_state_harvester && m_tpg_metric_collect_enabled) {
    m_state_harvester->start_collection_thread();
  }
  inherited::start(args);
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::stop(const appfwk::DAQModule::CommandData_t& args)
{ 
  inherited::stop(args);
  if (this->m_post_processing_enabled) {
    // Stop the state harvester collection thread if it exists
    if (m_state_harvester) {
      m_state_harvester->stop_collection_thread();
    }
    // Clears the pipelines and resets with the given configs.
    m_tp_generator->configure(m_tpg_configs, m_channel_plane_numbers, ReadoutTypeAdapter::samples_tick_difference);
  }
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::configure_source_and_geo_ids(const appmodel::DataHandlerModule* conf)
{
  m_sourceid.id = conf->get_source_id();
  m_sourceid.subsystem = ReadoutTypeAdapter::subsystem;
  auto geo_id = conf->get_geo_id();
  if (geo_id != nullptr) {
    m_det_id = geo_id->get_detector_id();
    m_crate_id = geo_id->get_crate_id();
    m_slot_id = geo_id->get_slot_id();
    m_stream_id = geo_id->get_stream_id();
  }
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::configure_preprocessing(const appmodel::DataHandlerModule* conf)
{
  m_emulator_mode = conf->get_emulation_mode();
  if (!m_emulator_mode) {
    inherited::add_preprocess_task(std::bind(&TPCEthFrameProcessor<ReadoutTypeAdapter>::sequence_check, this, std::placeholders::_1));
  }

  inherited::add_preprocess_task(std::bind(&TPCEthFrameProcessor<ReadoutTypeAdapter>::timestamp_check, this, std::placeholders::_1));
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::configure_channel_plane_numbers(const appmodel::TPCRawDataProcessor* proc_conf)
{
  const std::shared_ptr<detchannelmaps::TPCChannelMap> channel_map = dunedaq::detchannelmaps::make_tpc_map(proc_conf->get_channel_map());
  const std::vector<unsigned int> channel_mask_vec = proc_conf->get_channel_mask();

  for (int chan = 0; chan < 64; chan++) {
    trgdataformats::channel_t off_channel = channel_map->get_offline_channel_from_det_crate_slot_stream_chan(m_det_id, m_crate_id, m_slot_id, m_stream_id, chan);
    int16_t plane = channel_map->get_plane_from_offline_channel(off_channel);
    m_channel_plane_numbers.push_back(std::make_pair(off_channel, plane));

    // This processor only needs to handle some (maybe 0) of the masked channels.
    // Only get those relevant channels for the later check.
    // Only get the planes for the channels that are not masked.
    if (std::find(channel_mask_vec.begin(), channel_mask_vec.end(), off_channel) != channel_mask_vec.end()) {
      m_channel_mask_set.insert(off_channel);
    } else {
      m_channel_plane_map[off_channel] = plane;
      m_plane_numbers_set.insert(plane);
    }
  }
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::configure_find_tps(const appmodel::DataHandlerModule* conf, const appmodel::TPCRawDataProcessor* proc_conf)
{
  // Setting TP sinks.
  // Configurations currently have the sinks iterate in order, but there may be more sinks than planes.
  int plane_number = 0;
  for (auto output : conf->get_outputs()) {
    try {
      if (output->get_data_type() == "TriggerPrimitiveVector") {
         if (m_plane_numbers_set.contains(plane_number)) {
           m_plane_to_tp_sink_map[plane_number] = get_iom_sender<std::vector<trigger::TriggerPrimitiveTypeAdapter>>(output->UID());
         }
         plane_number++;
      }
    } catch (const ers::Issue& excpt) {
      ers::error(datahandlinglibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
    }
  }

  // We do need a coverage for all planes.
  if (m_plane_numbers_set.size() > m_plane_to_tp_sink_map.size()) {
      ers::error(DetectorPlaneToTPSinkMismatch(ERS_HERE, m_plane_numbers_set.size(), m_plane_to_tp_sink_map.size()));
  }

  m_tp_generator = std::make_unique<tpglibs::TPGenerator>();

  // Set the minimum TP samples over threshold.
  auto conf_sot_minima = proc_conf->get_sot_minima();
  std::vector<uint16_t> sot_minima{conf_sot_minima->get_sot_minimum_plane0(),
                                   conf_sot_minima->get_sot_minimum_plane1(),
                                   conf_sot_minima->get_sot_minimum_plane2()};
  m_tp_generator->set_sot_minima(sot_minima);

  std::vector<const appmodel::ProcessingStep*> processing_steps = proc_conf->get_processing_steps();
  for (auto step : processing_steps) {
    m_tpg_configs.push_back(std::make_pair(step->class_name(), step->to_json(false).back()));
  }

  // Let the TPG generator configure
  m_tp_generator->configure(m_tpg_configs, m_channel_plane_numbers, ReadoutTypeAdapter::samples_tick_difference);

  // Set the limits on when to send TPs and check that we can actually send on these limits.
  m_frame_count_limit = proc_conf->get_frame_count_limit();
  m_tp_count_limit = proc_conf->get_tp_count_limit();
  m_frame_limit_enabled = m_frame_count_limit > 0;
  m_tp_limit_enabled = m_tp_count_limit > 0;

  if (!m_frame_limit_enabled && !m_tp_limit_enabled) {
    throw FrameAndTPCountersDisabled(ERS_HERE);
  }

  m_metric_collect_opmon_period = proc_conf->get_metric_collect_opmon_period();

  // Check if metric collection is enabled in the configs
  m_tpg_metric_collect_enabled = false;
  for (const auto& name_config : m_tpg_configs) {
    if (name_config.second.contains("metric_collect_toggle_state") && 
        name_config.second["metric_collect_toggle_state"] == true) {
      m_tpg_metric_collect_enabled = true;
      break;
    }
  }

  // Only create and configure state harvester if metric collection is enabled
  if (m_tpg_metric_collect_enabled) {
    auto processsor_references = m_tp_generator->get_all_processor_references_with_pipeline_index();

    m_state_harvester = std::make_unique<fdreadoutlibs::TPGInternalStateHarvester>();

    const uint8_t channels_per_pipeline = 16;
    const uint8_t pipelines = static_cast<uint8_t>(m_channel_plane_numbers.size() / channels_per_pipeline);
    
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Configuring state harvester with " << static_cast<int>(channels_per_pipeline) 
                                  << " channels per pipeline, " << static_cast<int>(pipelines) << " pipelines, " 
                                  << processsor_references.size() << " processor references";
    
    m_state_harvester->update_channel_plane_numbers(m_channel_plane_numbers,
                                                    channels_per_pipeline, pipelines);
    m_state_harvester->set_processor_references(processsor_references);
    
    // Start the collection thread immediately after configuration
    m_state_harvester->start_collection_thread();
    
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "State harvester configured and started successfully";
  }

  inherited::add_postprocess_task(std::bind(&TPCEthFrameProcessor<ReadoutTypeAdapter>::find_tps, this, std::placeholders::_1));
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::configure_postprocessing(const appmodel::DataHandlerModule* conf)
{
  const appmodel::DataProcessor* dp = conf->get_module_configuration()->get_data_processor();
  if (dp == nullptr) {
    return;
  }

  const appmodel::TPCRawDataProcessor* proc_conf = dp->cast<appmodel::TPCRawDataProcessor>();
  if (proc_conf == nullptr) {
    return;
  }

  // Need TPCRawDataProcessor configurations to configure the following.
  configure_channel_plane_numbers(proc_conf);
  configure_find_tps(conf, proc_conf);
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::conf(const appmodel::DataHandlerModule* conf)
{
  configure_source_and_geo_ids(conf);

  configure_preprocessing(conf);

  if (this->m_post_processing_enabled) {
    configure_postprocessing(conf);
  }

  inherited::conf(conf);
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::scrap_source_and_geo_ids()
{
  m_sourceid = daqdataformats::SourceID();

  m_det_id = 0;
  m_crate_id = 0;
  m_slot_id = 0;
  m_stream_id = 0;
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::scrap_preprocessing()
{
  m_emulator_mode = false;
  m_first_frame = true;

  // Timestamps.
  m_previous_ts = 0;
  m_current_ts = 0;

  m_pattern_generator_previous_ts = 0;
  m_pattern_generator_current_ts = 0;

  m_first_ts_missmatch = true;
  m_ts_problem_reported = false;
  m_ts_error_state = false;
  m_ts_error_ctr = 0;

  // Sequence ID.
  m_previous_seq_id = 0;
  m_current_seq_id = 0;

  m_first_seq_id_mismatch = true;
  m_seq_id_problem_reported = false;
  m_seq_id_error_state = false;
  m_seq_id_error_ctr = 0;
  m_seq_id_min_jump = 0;
  m_seq_id_max_jump = 0;

  // The preprocessing tasks scrap is handled by inherited::scrap().
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::scrap_postprocessing()
{
  // Channel-plane variables
  m_channel_mask_set.clear();
  m_plane_numbers_set.clear();
  m_channel_plane_numbers.clear();
  m_channel_plane_map.clear();

  // TP variables
  m_tp_generator->reset();
  m_tpg_configs.clear();
  m_plane_to_tpa_vector_map.clear();
  m_plane_to_tp_sink_map.clear();

  m_frame_limit_enabled = false;
  m_tp_limit_enabled = false;
  m_current_tp_count = 0;
  m_tp_count_limit = 0;
  m_frame_count_at_last_send = 0;

  // OpMon variables
  m_tpg_metric_collect_enabled = false;
  m_metric_collect_opmon_period = 128;
  m_tp_channel_rate_map.clear();
  m_num_new_tps.exchange(0);
  m_tps_suppressed_too_long.exchange(0);
  m_tps_send_failed.exchange(0);
  m_frame_counter.exchange(0);
  m_t0 = std::chrono::high_resolution_clock::now();
}

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::scrap(const appfwk::DAQModule::CommandData_t& cfg)
{
  scrap_source_and_geo_ids();
  scrap_preprocessing();

  if (this->m_post_processing_enabled) {
    scrap_postprocessing();
  }

  inherited::scrap(cfg);
}


template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::generate_opmon_data()
{
   datahandlinglibs::opmon::FixedRateDataProcessorInfo info;

   info.set_num_seq_id_errors(m_seq_id_error_ctr.load());
   info.set_min_seq_id_jump(m_seq_id_min_jump.exchange(0));
   info.set_max_seq_id_jump(m_seq_id_max_jump.exchange(0));

   info.set_num_ts_errors(m_ts_error_ctr.load());

   this->publish(std::move(info));

   this->m_error_registry->log_registered_errors();

   if (this->m_post_processing_enabled) {
     auto now = std::chrono::high_resolution_clock::now();
     int num_new_tps = m_num_new_tps.exchange(0);
     int num_new_tps_suppressed_too_long = m_tps_suppressed_too_long.exchange(0);
     int num_new_tps_send_failed = m_tps_send_failed.exchange(0);
     double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
     TLOG_DEBUG(TLVL_BOOKKEEPING) << "TP rate: " << std::to_string(num_new_tps / seconds / 1000.) << " [kHz]";
     TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new TPs: " << num_new_tps;

     datahandlinglibs::opmon::HitFindingInfo tp_info;
     tp_info.set_rate_tp_hits(num_new_tps / seconds / 1000.);

     tp_info.set_num_tps_sent(num_new_tps);
     tp_info.set_num_tps_suppressed_too_long(num_new_tps_suppressed_too_long);
     tp_info.set_num_tps_send_failed(num_new_tps_send_failed);

     this->publish(std::move(tp_info));
     // Find the channels with the top TP rates
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
         this->publish(std::move(tpc_info), {{"channel", std::to_string(channel_tp_rate_vec[i].first)}});
       }
     }

     // Reset the counter in the channel rate map
     for (auto& el : m_tp_channel_rate_map) {
       el.second = 0;
     }
     m_t0 = now;

    if (m_tpg_metric_collect_enabled && m_state_harvester) {
      publish_processor_metric_to_opmon();
      publish_processor_metric_to_opmon_with_aggregation();
    }
   }

   inherited::generate_opmon_data();
 }

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::publish_processor_metric_to_opmon() {
  if (!m_state_harvester) {
    return;
  }
  
  // Get latest results from background collection thread
  auto metrics = m_state_harvester->get_latest_results();
  
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Publishing processor metrics for " << metrics.size() << " channels";
  
  int metrics_published = 0;
  
  // Publish per-channel metrics
  for (const auto& [channel, vec] : metrics) {
    datahandlinglibs::opmon::TPGProcessorInfo tpg_proc_info;
    bool has_valid_metrics = false;
    
    for (const auto& [name, val] : vec) {
      if (name == "pedestal") {
        tpg_proc_info.set_pedestal(val);
        has_valid_metrics = true;
      } else if (name == "accum") {
        tpg_proc_info.set_accum(val);
        has_valid_metrics = true;
      }
    }
    
    if (has_valid_metrics) {
      this->publish(std::move(tpg_proc_info), {{"channel", std::to_string(channel)}});
      metrics_published++;
    }
  }
  
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Published " << metrics_published << " channel metrics";
}

template <class ReadoutTypeAdapter>
std::map<int16_t, std::map<std::string, std::tuple<float, int16_t, int16_t, float, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>>
TPCEthFrameProcessor<ReadoutTypeAdapter>::calculate_all_metric_summaries_across_planes(const std::unordered_map<dunedaq::trgdataformats::channel_t, std::vector<std::pair<std::string, int16_t>>>& metrics) {
    // Structure to hold all statistics: plane -> metric -> (mean, min, max, stddev, min_channel_id, max_channel_id)
    std::map<int16_t, std::map<std::string, std::tuple<float, int16_t, int16_t, float, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>> all_stats;

    // Structure to accumulate statistics: plane -> metric -> (count, mean, M2, min, max, min_channel_id, max_channel_id)
    std::map<int16_t, std::map<std::string, std::tuple<size_t, double, double, int16_t, int16_t, dunedaq::trgdataformats::channel_t, dunedaq::trgdataformats::channel_t>>> accumulators;

    // Single pass through all metrics to collect data using Welford's online algorithm for variance
    for (const auto& [channel, vec] : metrics) {
        if (m_channel_plane_map.empty()) continue;

        int16_t plane = m_channel_plane_map[channel];

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

template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::publish_processor_metric_to_opmon_with_aggregation() {
  if (!m_state_harvester) {
    return;
  }
  
  // Get latest results from background collection thread
  auto metrics = m_state_harvester->get_latest_results();
  
  // Use optimized single-pass calculation for all metrics across all planes
  auto all_stats = calculate_all_metric_summaries_across_planes(metrics);
  
  TLOG_DEBUG(TLVL_BOOKKEEPING) << "Publishing aggregated metrics for " << all_stats.size() << " planes";
  
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
      this->publish(std::move(info), {{"plane", std::to_string(plane)}, {"metric", metric_name}});
    }
  }
}


/**
 * Pipeline Stage 1.: Check proper timestamp increments in TPC frame
 * */
template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::sequence_check(frameptr fp)
{
  // Acquire timestamp
  auto wfptr = reinterpret_cast<tpcframeptr>(fp); // NOLINT
  m_current_seq_id = wfptr->daq_header.seq_id;

  // Check that the system is properly configured from the first frame.
  if (m_first_frame) [[unlikely]] {
    if (wfptr->daq_header.crate_id != m_crate_id || wfptr->daq_header.slot_id != m_slot_id || wfptr->daq_header.stream_id != m_stream_id) {
      ers::error(LinkMisconfiguration(ERS_HERE, wfptr->daq_header.crate_id, wfptr->daq_header.slot_id, wfptr->daq_header.stream_id, m_crate_id, m_slot_id, m_stream_id));
    }

    m_first_frame = false;
  }

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
        this->m_error_registry->add_error("Sequence ID jump", datahandlinglibs::FrameErrorRegistry::ErrorInterval(expected_seq_id, m_current_seq_id));
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
 * Pipeline Stage 1.: Check proper timestamp increments in TPC frame
 * */
template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::timestamp_check(frameptr fp)
{

  uint16_t tpceth_tick_difference = ReadoutTypeAdapter::expected_tick_difference;
  uint16_t tpceth_frame_tick_difference = tpceth_tick_difference * fp->get_num_frames();

  auto wfptr = reinterpret_cast<tpcframeptr>(fp); // NOLINT
  m_current_ts = wfptr->get_timestamp();

  // Check timestamp
  if (m_previous_ts > 0 &&
      m_current_ts - m_previous_ts != tpceth_frame_tick_difference) [[unlikely]] {
    ++m_ts_error_ctr;
    if (m_first_ts_missmatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISMATCH! -> | previous: " << std::to_string(m_previous_ts) << " current: " + std::to_string(m_current_ts);
      m_first_ts_missmatch = false;
    } else {
      if (!m_ts_error_state) {
        this->m_error_registry->add_error("Timestamp jump", datahandlinglibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + tpceth_frame_tick_difference, m_current_ts));
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
  this->m_last_processed_daq_ts = m_current_ts;
}

/**
 * Pipeline Stage 2.: Do software TPG
 * */
template <class ReadoutTypeAdapter>
void
TPCEthFrameProcessor<ReadoutTypeAdapter>::find_tps(constframeptr fp)
{
  if (!fp)
    return;
  auto wfptr = reinterpret_cast<tpcframeptr>((uint8_t*)fp); // NOLINT

  std::vector<trgdataformats::TriggerPrimitive> tps = (*m_tp_generator)(wfptr);

  uint64_t current_frame_count = m_frame_counter.fetch_add(1, std::memory_order_relaxed) + 1;
  
  // Trigger asynchronous metric collection in background thread
  if (m_tpg_metric_collect_enabled && m_state_harvester && 
      current_frame_count % m_metric_collect_opmon_period == 0) {
    m_state_harvester->trigger_harvest();
  }

  for (const auto& tp : tps) {
    // If this TP is on a masked channel, skip it.
    if (std::binary_search(m_channel_mask_set.begin(), m_channel_mask_set.end(), uint32_t(tp.channel)))
      continue;
    // Need to move into a type adapter.
    trigger::TriggerPrimitiveTypeAdapter tpa;
    tpa.tp = tp;

    tpa.tp.detid = m_det_id;  // Last missing piece.
    m_plane_to_tpa_vector_map[m_channel_plane_map[uint32_t(tp.channel)]].push_back(tpa);
    m_tp_channel_rate_map[uint32_t(tp.channel)]++;
    m_current_tp_count++;
  }

  const bool frame_limit_reached = m_frame_limit_enabled && (current_frame_count - m_frame_count_at_last_send >= m_frame_count_limit);
  const bool tp_limit_reached = m_tp_limit_enabled && (m_current_tp_count >= m_tp_count_limit);

  if (frame_limit_reached || tp_limit_reached) [[unlikely]] {
    m_current_tp_count = 0;
    m_frame_count_at_last_send = current_frame_count;
    for (auto& [plane_num, tpa_vector] : m_plane_to_tpa_vector_map) {
      int num_new_tps = tpa_vector.size();
      if (num_new_tps == 0) {
        continue;
      }
      const auto ts_begin = tpa_vector.front().tp.time_start;
      const auto channel_begin = tpa_vector.front().tp.channel;
      const auto ts_end = tpa_vector.back().tp.time_start;
      const auto channel_end = tpa_vector.back().tp.channel;
      if (!m_plane_to_tp_sink_map[plane_num]->try_send(std::move(tpa_vector), iomanager::Sender::s_no_block)) {
        ers::warning(FailedToSendTPVector(ERS_HERE, ts_begin, channel_begin, ts_end, channel_end));
        m_tps_send_failed++;
      } else {
        m_num_new_tps += num_new_tps;
      }
    }
  }
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
