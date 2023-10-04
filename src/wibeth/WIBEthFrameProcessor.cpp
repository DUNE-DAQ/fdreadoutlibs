/**
 * @file WIBEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/wibeth/WIBEthFrameProcessor.hpp" // NOLINT(build/include)

#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/IterableQueueModel.hpp"
#include "readoutlibs/readoutconfig/Nljs.hpp"
#include "readoutlibs/readoutinfo/InfoNljs.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"
#include "fddetdataformats/WIBEthFrame.hpp"


#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "fdreadoutlibs/wibeth/tpg/DesignFIR.hpp"
#include "fdreadoutlibs/wibeth/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/ProcessRSAVX2.hpp"
#include "fdreadoutlibs/wibeth/tpg/TPGConstants_wibeth.hpp"

#include <atomic>
#include <bitset>
#include <functional>
#include <future>
#include <memory>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

// THIS SHOULDN'T BE HERE!!!!! But it is necessary.....
DUNE_DAQ_TYPESTRING(dunedaq::fdreadoutlibs::types::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")


namespace dunedaq {
namespace fdreadoutlibs {


WIBEthFrameHandler::WIBEthFrameHandler()
  : m_hits_dest(nullptr)
  , m_tpg_taps_p(nullptr)
{}

WIBEthFrameHandler::~WIBEthFrameHandler()
{
  if (m_tpg_taps_p) {
    delete[] m_tpg_taps_p;
  }
  if (m_hits_dest)
	  delete[] m_hits_dest;
}

void
WIBEthFrameHandler::reset()
{
  if (m_tpg_taps_p)
      	delete[] m_tpg_taps_p;
  m_tpg_taps_p = nullptr;
  if (m_hits_dest)
        delete[] m_hits_dest;
  m_hits_dest = nullptr;
  first_hit = true;
}

void
WIBEthFrameHandler::initialize(int threshold_value)
{
  m_tpg_taps = swtpg_wibeth::firwin_int(7, 0.1, m_tpg_multiplier);
  m_tpg_taps.push_back(0);

  m_tpg_threshold = threshold_value;

  if (m_tpg_taps_p == nullptr) {
    m_tpg_taps_p = new int16_t[m_tpg_taps.size()];
  }
  for (size_t i = 0; i < m_tpg_taps.size(); ++i) {
    m_tpg_taps_p[i] = m_tpg_taps[i];
  }

  if(m_hits_dest == nullptr) {
	m_hits_dest = new uint16_t[100000];
  }

  m_tpg_processing_info = std::make_unique<swtpg_wibeth::ProcessingInfo<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>>(nullptr,
                                                                                                            swtpg_wibeth::FRAMES_PER_MSG,
                                                                                                            0,
                                                                                                            swtpg_wibeth::NUM_REGISTERS_PER_FRAME,
                                                                                                            m_hits_dest,
                                                                                                            m_tpg_taps_p,
                                                                                                            (uint8_t)m_tpg_taps.size(), // NOLINT(build/unsigned)
                                                                                                            m_tpg_tap_exponent,
                                                                                                            m_tpg_threshold,
                                                                                                            0,
                                                                                                            0);
}

// Get destination ptr for the frame handler
uint16_t*
WIBEthFrameHandler::get_hits_dest()
{
  return m_hits_dest;
}


WIBEthFrameProcessor::WIBEthFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
  : TaskRawDataProcessorModel<types::DUNEWIBEthTypeAdapter>(error_registry)
  , m_tpg_enabled(false)
{
}

WIBEthFrameProcessor::~WIBEthFrameProcessor()
{
  m_wibeth_frame_handler->reset();
}

void
WIBEthFrameProcessor::start(const nlohmann::json& args)
{
  // Reset software TPG resources
  if (m_tpg_enabled) {
    m_tps_dropped = 0;

    m_wibeth_frame_handler->initialize(m_tpg_threshold_selected);
  } // end if(m_tpg_enabled)

  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_ts_problem_reported = false;
  m_ts_error_ctr = 0;

  m_first_seq_id_mismatch = true;
  m_seq_id_problem_reported = false;
  m_seq_id_error_ctr = 0;


  // Reset stats
  m_t0 = std::chrono::high_resolution_clock::now();
  m_new_hits = 0;
  m_new_tps = 0;
  m_tpg_hits_count.exchange(0);
  inherited::start(args);
}

void
WIBEthFrameProcessor::stop(const nlohmann::json& args)
{
  inherited::stop(args);
  if (m_tpg_enabled) {
    // Make temp. buffers reusable on next start.
    m_wibeth_frame_handler->reset();
  }
}

void
WIBEthFrameProcessor::init(const nlohmann::json& args)
{
//  inherited::init(args);

  try {
    auto queue_index = appfwk::connection_index(args, {});
    if (queue_index.find("tp_out") != queue_index.end()) {
      m_tp_sink = get_iom_sender<types::TriggerPrimitiveTypeAdapter>(queue_index["tp_out"]);
    }
  } catch (const ers::Issue& excpt) {
    ers::error(readoutlibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
  }

}

void
WIBEthFrameProcessor::conf(const nlohmann::json& cfg)
{
  auto config = cfg["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();

  m_sourceid.id = config.source_id;
  m_sourceid.subsystem = types::DUNEWIBEthTypeAdapter::subsystem;

  m_tpg_algorithm = config.tpg_algorithm;  
  TLOG() << "Selected software TPG algorithm: " << m_tpg_algorithm;
  if (m_tpg_algorithm == "SimpleThreshold") {
    m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
  } else if (m_tpg_algorithm == "AbsRS" ) {
    m_assigned_tpg_algorithm_function = &swtpg_wibeth::process_window_rs_avx2<swtpg_wibeth::NUM_REGISTERS_PER_FRAME>;
  } else {
    throw TPGAlgorithmInexistent(ERS_HERE, m_tpg_algorithm);
  }

  m_tp_max_width = config.tp_timeout;

  m_channel_mask_vec = config.tpg_channel_mask;
  // Converting the input vector of channels masks into an std::set
  // AAA: The set provides faster look up than a std::vector
  m_channel_mask_set.insert(m_channel_mask_vec.begin(), m_channel_mask_vec.end());

  m_tpg_threshold_selected = config.tpg_threshold;

  m_crate_no = config.crate_id;
  m_slot_no = config.slot_id;
  m_stream_id = config.link_id;
  // Setup pre-processing pipeline
  inherited::add_preprocess_task(std::bind(&WIBEthFrameProcessor::sequence_check, this, std::placeholders::_1));
  inherited::add_preprocess_task(std::bind(&WIBEthFrameProcessor::timestamp_check, this, std::placeholders::_1));
  if (config.enable_tpg) {
    m_tpg_enabled = true;
    m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

    inherited::add_postprocess_task(std::bind(&WIBEthFrameProcessor::find_hits, this, std::placeholders::_1, m_wibeth_frame_handler.get()));
  }

  inherited::conf(cfg);
}

void
WIBEthFrameProcessor::get_info(opmonlib::InfoCollector& ci, int level)
{
  readoutlibs::readoutinfo::RawDataProcessorInfo info;

  info.num_seq_id_errors = m_seq_id_error_ctr.load();
  info.min_seq_id_jump = m_seq_id_min_jump.exchange(0);
  info.max_seq_id_jump = m_seq_id_max_jump.exchange(0);

  info.num_ts_errors = m_ts_error_ctr.load();
  

  auto now = std::chrono::high_resolution_clock::now();
  if (m_tpg_enabled) {
    int new_hits = m_tpg_hits_count.exchange(0);
    int new_tps = m_new_tps.exchange(0);
    int new_dropped_tps = m_tps_dropped.exchange(0);
    double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
    //TLOG() << " Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz], dropped rate: " << std::to_string(new_dropped_tps / seconds / 1000.) << " [kHz]";;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new hits: " << new_hits << " new TPs: " << new_tps;
    info.rate_tp_hits = new_hits / seconds / 1000.;

    info.num_tps_sent = new_tps;
    info.num_tps_dropped = new_dropped_tps;
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
      for (int i = 0; i < top_highest_values; i++) {
        std::stringstream info_name;
        info_name << "channel_" << channel_tp_rate_vec[i].first;
        opmonlib::InfoCollector tmp_ic;
        readoutlibs::readoutinfo::TPChannelInfo tp_info;
        tp_info.num_tp = channel_tp_rate_vec[i].second;
        tmp_ic.add(tp_info);
        ci.add(info_name.str(), tmp_ic);
      }
    }

    // Reset the counter in the channel rate map
    for (auto& el : m_tp_channel_rate_map) {
      el.second = 0;
    }
  }
  m_t0 = now;
  inherited::get_info(ci, level);
  ci.add(info);
}


/**
 * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
 * */
void
WIBEthFrameProcessor::sequence_check(frameptr fp)
{

  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
    // uint64_t ts_next = m_previous_seq_id + 1; // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));            // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
      //auto wfh = const_cast<dunedaq::fddetdataformats::WIBEthFrame*>(wf->header());
      wf->daq_header.crate_id = m_crate_no;
      wf->daq_header.slot_id = m_slot_no;
      wf->daq_header.stream_id = m_stream_id; 
      wf->daq_header.seq_id = (m_previous_seq_id+i) & 0xfff;
      wf++;
    }
  }

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

  if (delta_seq_id != 0) {
    // uint16_t delta_seq_id = (m_current_seq_id-expected_seq_id);
    ++m_seq_id_error_ctr;
    m_seq_id_max_jump = std::max(delta_seq_id, m_seq_id_max_jump.load());
    m_seq_id_min_jump = std::min(delta_seq_id, m_seq_id_min_jump.load());

    m_error_registry->add_error("SEQUENCE_ID_JUMP", readoutlibs::FrameErrorRegistry::ErrorInterval(expected_seq_id, m_current_seq_id));
    if (m_first_seq_id_mismatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First sequence id MISSMATCH! -> | previous: " << std::to_string(m_previous_seq_id) << " current: " + std::to_string(m_current_seq_id);
      m_first_seq_id_mismatch = false;
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

  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
    uint64_t ts_next = m_previous_ts + wibeth_frame_tick_difference; // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));            // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
      //auto wfh = const_cast<dunedaq::fddetdataformats::WIBEthFrame*>(wf->header());
      wf->daq_header.crate_id = m_crate_no;
      wf->daq_header.slot_id = m_slot_no;
      wf->daq_header.stream_id = m_stream_id; 
      wf->set_timestamp(ts_next);
      ts_next += wibeth_tick_difference;
      wf++;
    }
  }

  // Acquire timestamp
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>(fp); // NOLINT
  m_current_ts = wfptr->get_timestamp();

  // Check timestamp
  if (m_current_ts - m_previous_ts != wibeth_frame_tick_difference) {
    ++m_ts_error_ctr;
    m_error_registry->add_error("MISSING_FRAMES", readoutlibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + wibeth_frame_tick_difference, m_current_ts));
    if (m_first_ts_missmatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISSMATCH! -> | previous: " << std::to_string(m_previous_ts) << " current: " + std::to_string(m_current_ts);
      m_first_ts_missmatch = false;
    }
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
WIBEthFrameProcessor::find_hits(constframeptr fp, WIBEthFrameHandler* frame_handler)
{
  if (!fp)
    return;
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBEthFrame*>((uint8_t*)fp); // NOLINT
  uint64_t timestamp = wfptr->get_timestamp();                                            // NOLINT(build/unsigned)

  // Frame expansion
  swtpg_wibeth::MessageRegisters registers_array;
  expand_wibeth_adcs(fp, &registers_array);

  // For debugging purposes you can check the single ADCs 
  //parse_wibeth_adcs(&registers_array);

  // Only for the first WIBEth frame, create an offline register map
  if (frame_handler->first_hit) {
    frame_handler->register_channel_map = swtpg_wibeth::get_register_to_offline_channel_map_wibeth(wfptr, m_channel_map);

    frame_handler->m_tpg_processing_info->setState(registers_array);

    m_det_id = wfptr->daq_header.det_id;
    if (wfptr->daq_header.crate_id != m_crate_no || wfptr->daq_header.slot_id != m_slot_no || wfptr->daq_header.stream_id != m_stream_id) {
      ers::error(LinkMisconfiguration(ERS_HERE, wfptr->daq_header.crate_id, wfptr->daq_header.slot_id, wfptr->daq_header.stream_id, m_crate_no, m_slot_no, m_stream_id));
    }
    // Add WIBEthFrameHandler channel map to the common m_register_channels.
    // Populate the array 
    for (size_t i = 0; i < swtpg_wibeth::NUM_REGISTERS_PER_FRAME * swtpg_wibeth::SAMPLES_PER_REGISTER; ++i) {
      m_register_channels[i] = frame_handler->register_channel_map.channel[i];

      //TLOG () << "Index number " << i << " offline channel " << frame_handler->register_channel_map.channel[i]; 

      // Set up a map of channels and number of TPs for monitoring/debug
      m_tp_channel_rate_map[frame_handler->register_channel_map.channel[i]] = 0;
    }

    //TLOG() << "Processed the first frame ";

    // Set first hit bool to false so that registration of channel map is not executed twice
    frame_handler->first_hit = false;

  } // end if (frame_handler->first_hit)


  // Execute the SWTPG algorithm
  frame_handler->m_tpg_processing_info->input = &registers_array;
  // Set the first word to "magic" indicating there is no hit, initially
  frame_handler->m_tpg_processing_info->output[0] = swtpg_wibeth::MAGIC;

  m_assigned_tpg_algorithm_function(*frame_handler->m_tpg_processing_info);

  //GLM: avoid the tp_handler queue/thread
  process_swtpg_hits(frame_handler->m_tpg_processing_info->output, timestamp);

}

void
WIBEthFrameProcessor::process_swtpg_hits(uint16_t* primfind_it, dunedaq::daqdataformats::timestamp_t timestamp)
{

  constexpr int clocksPerTPCTick = types::DUNEWIBEthTypeAdapter::samples_tick_difference;

  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16], hit_peak_time[16], hit_peak_adc[16]; // NOLINT(build/unsigned)
  unsigned int nhits = 0;

  while (*primfind_it != swtpg_wibeth::MAGIC) {
    // First, get all of the register values (including those with no hit) into local variables
    for (int i = 0; i < 16; ++i) {
      chan[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
    }
    for (int i = 0; i < 16; ++i) {
      hit_end[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
    }
    for (int i = 0; i < 16; ++i) {
      hit_charge[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      // TLOG() << "hit_charge:" << hit_charge[i];
    }
    for (int i = 0; i < 16; ++i) {
      // hit_tover[i] = static_cast<uint16_t>(*primfind_it++); // NOLINT(runtime/increment_decrement)
      hit_tover[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
    }
    for (int i = 0; i < 16; ++i) {
      hit_peak_time[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
    }
    for (int i = 0; i < 16; ++i) {
      hit_peak_adc[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
    }  

    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (int i = 0; i < 16; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wibeth::MAGIC) {

        uint64_t tp_t_begin =                                                           // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_tover[i])); // NOLINT(build/unsigned)
        uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);         // NOLINT(build/unsigned)
        //TLOG() << "Hit start " << tp_t_begin << ", end: " << tp_t_end << ", online channel: " << chan[i];
	
	uint64_t tp_t_peak =                                                           // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_peak_time[i])); // NOLINT(build/unsigned)

        // This channel had a hit ending here, so we can create and output the hit here
        const uint16_t offline_channel = m_register_channels[chan[i]];
        if (m_channel_mask_set.find(offline_channel) == m_channel_mask_set.end()) {
          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          // sed -n -e 's/.*Hit: \(.*\) \(.*\).*/\1 \2/p' log.txt  > hits.txt
          //

	  fdreadoutlibs::types::TriggerPrimitiveTypeAdapter tp;
          tp.tp.time_start = tp_t_begin;
          tp.tp.time_peak = tp_t_peak;
          tp.tp.time_over_threshold = int64_t(hit_tover[i]) * clocksPerTPCTick;
          tp.tp.channel = offline_channel;
          tp.tp.adc_integral = hit_charge[i];
          tp.tp.adc_peak = hit_peak_adc[i];
          tp.tp.detid =  m_det_id; // TODO: convert crate/slot/link to SourceID Roland Sipos rsipos@cern.ch July-22-2021
          tp.tp.type = trgdataformats::TriggerPrimitive::Type::kTPC;
          tp.tp.algorithm = trgdataformats::TriggerPrimitive::Algorithm::kTPCDefault;
          tp.tp.version = 1;
          if(tp.tp.time_over_threshold > m_tp_max_width) {
		  ers::warning(TPTooLong(ERS_HERE, tp.tp.time_over_threshold, tp.tp.channel));
		  m_tps_dropped++;
	  }
	  //Send the TP to the TP handler module
	  else if(!m_tp_sink->try_send(std::move(tp), iomanager::Sender::s_no_block)) {
		 ers::warning(TPDropped(ERS_HERE, tp.tp.time_start, tp.tp.channel));
		 m_tps_dropped++;
	  }	 
          m_new_tps++;
          ++nhits;

          // Update the channel/rate map. Increment the value associated with the TP channel
          m_tp_channel_rate_map[offline_channel]++;
        }
      }
    }
  }
  m_tpg_hits_count += nhits;
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
