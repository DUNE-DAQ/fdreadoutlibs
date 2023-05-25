/**
 * @file WIB2FrameProcessor.hpp WIB2 specific Task based raw processor
 * @author Adam Abed Abud (adam.abed.abud@cern.ch)
 *
 * This is part of the DUNE DAQ , copyright 2022.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/wib2/WIB2FrameProcessor.hpp" // NOLINT(build/include)

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
#include "fddetdataformats/WIB2Frame.hpp"


#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "fdreadoutlibs/wib2/tpg/DesignFIR.hpp"
#include "fdreadoutlibs/wib2/tpg/FrameExpand.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessAVX2.hpp"
#include "fdreadoutlibs/wib2/tpg/ProcessRSAVX2.hpp"
#include "fdreadoutlibs/wib2/tpg/TPGConstants_wib2.hpp"

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

void
WIB2PatternGenerator::generate(int source_id)
{
  //TLOG() << "Generate random ADC patterns" ;
  std::srand(source_id*12345678);
  m_channel.reserve(m_size);
  for (int i = 0; i < m_size; i++) {
      int random_ch = std::rand()%256;
      m_channel.push_back(random_ch);
  }
}

WIB2FrameHandler::WIB2FrameHandler(int register_selector_params)
  : m_register_selector(register_selector_params)
  , m_hits_dest(nullptr)
  , m_tpg_taps_p(nullptr)
{}

WIB2FrameHandler::~WIB2FrameHandler()
{
  if (m_tpg_taps_p) {
    delete[] m_tpg_taps_p;
  }
  if (m_hits_dest)
	  delete[] m_hits_dest;
}

int
WIB2FrameHandler::get_registers_selector()
{
  return m_register_selector;
}

void
WIB2FrameHandler::reset()
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
WIB2FrameHandler::initialize(int threshold_value)
{
  m_tpg_taps = swtpg_wib2::firwin_int(7, 0.1, m_tpg_multiplier);
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

  m_tpg_processing_info = std::make_unique<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>>(nullptr,
                                                                                                            swtpg_wib2::FRAMES_PER_MSG,
                                                                                                            0,
                                                                                                            swtpg_wib2::NUM_REGISTERS_PER_FRAME,
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
WIB2FrameHandler::get_hits_dest()
{
  return m_hits_dest;
}


WIB2FrameProcessor::WIB2FrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
  : TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>(error_registry)
  , m_sw_tpg_enabled(false)
{
}

WIB2FrameProcessor::~WIB2FrameProcessor()
{
  m_wib2_frame_handler->reset();
  m_wib2_frame_handler_second_half->reset();
}

void
WIB2FrameProcessor::start(const nlohmann::json& args)
{
  // Reset software TPG resources
  if (m_sw_tpg_enabled) {
    m_tps_dropped = 0;

    m_wib2_frame_handler->initialize(m_tpg_threshold_selected);
    m_wib2_frame_handler_second_half->initialize(m_tpg_threshold_selected);
  } // end if(m_sw_tpg_enabled)

  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_problem_reported = false;
  m_ts_error_ctr = 0;

  // Reset stats
  m_t0 = std::chrono::high_resolution_clock::now();
  m_new_hits = 0;
  m_new_tps = 0;
  m_swtpg_hits_count.exchange(0);
  inherited::start(args);
}

void
WIB2FrameProcessor::stop(const nlohmann::json& args)
{
  inherited::stop(args);
  if (m_sw_tpg_enabled) {
    // Make temp. buffers reusable on next start.
    m_wib2_frame_handler->reset();
    m_wib2_frame_handler_second_half->reset();

    auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_t0).count();
  }
}

void
WIB2FrameProcessor::init(const nlohmann::json& args)
{
//  inherited::init(args);

  try {
	  auto qi = appfwk::connection_index(args, {"tp_out"});
	  m_tp_sink = iomanager::IOManager::get()->get_sender<dunedaq::fdreadoutlibs::types::TriggerPrimitiveTypeAdapter>(qi["tp_out"]);
  } catch (const ers::Issue& excpt) {
	  ers::warning(readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt));
  }

}

void
WIB2FrameProcessor::conf(const nlohmann::json& cfg)
{
  auto config = cfg["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();

  m_sourceid.id = config.source_id;
  m_sourceid.subsystem = types::DUNEWIBSuperChunkTypeAdapter::subsystem;
  m_tpg_algorithm = config.software_tpg_algorithm;
  m_tp_max_width = config.tp_timeout;
  TLOG() << "Selected software TPG algorithm: " << m_tpg_algorithm;

  m_channel_mask_vec = config.software_tpg_channel_mask;
  // Converting the input vector of channels masks into an std::set
  // AAA: The set provides faster look up than a std::vector
  m_channel_mask_set.insert(m_channel_mask_vec.begin(), m_channel_mask_vec.end());

  m_tpg_threshold_selected = config.software_tpg_threshold;

  m_crate_no = config.crate_id;
  m_slot_no = config.slot_id;
  m_link = config.link_id;
  // Setup pre-processing pipeline
  inherited::add_preprocess_task(std::bind(&WIB2FrameProcessor::timestamp_check, this, std::placeholders::_1));
  if (config.enable_software_tpg) {
    m_sw_tpg_enabled = true;
    if (config.emulator_mode) {
      m_wib2_pattern_generator.generate(m_sourceid.id);
      m_random_channels = m_wib2_pattern_generator.get_channels();
      inherited::add_preprocess_task(std::bind(&WIB2FrameProcessor::use_pattern_generator, this, std::placeholders::_1));
    }

    m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

    inherited::add_postprocess_task(std::bind(&WIB2FrameProcessor::find_hits, this, std::placeholders::_1, m_wib2_frame_handler.get()));
    inherited::add_postprocess_task(std::bind(&WIB2FrameProcessor::find_hits, this, std::placeholders::_1, m_wib2_frame_handler_second_half.get()));
  }

  inherited::conf(cfg);
}

void
WIB2FrameProcessor::get_info(opmonlib::InfoCollector& ci, int level)
{
  readoutlibs::readoutinfo::RawDataProcessorInfo info;

  auto now = std::chrono::high_resolution_clock::now();
  if (m_sw_tpg_enabled) {
    int new_hits = m_swtpg_hits_count.exchange(0);
    int new_tps = m_new_tps.exchange(0);
    int new_dropped_tps = m_tps_dropped.exchange(0);
    double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
    TLOG() << " Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz], dropped rate: " << std::to_string(new_dropped_tps / seconds / 1000.) << " [kHz]";;
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
 * Add hits using the pattern generator only when in emulated mode
 * */
void
WIB2FrameProcessor::use_pattern_generator(frameptr fp)
{

  // If we are not in the first superchunk then we start applying the pattern generator
  // This is because we use the ADC values of the first wib frame as the pedestal baseline
  if (m_current_ts != 0) {
    auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>((uint8_t*)fp);

    m_pattern_generator_current_ts = wfptr->get_timestamp();

    // Adding a hit every 2442 gives a total Sent TP rate of approx 100 Hz/wire
    if (m_pattern_generator_current_ts - m_pattern_generator_previous_ts > 2442) {

      // Reset the pattern from the beginning if it reaches the maximum
      m_pattern_index++;
      if (m_pattern_index == m_wib2_pattern_generator.get_total_size()) {
        m_pattern_index = 0;
      }

      // Set the ADC to the uint16 maximum value
     
      wfptr->set_adc(m_random_channels[m_pattern_index], 16383);
      //TLOG() << "Lift channel " << m_random_channels[m_pattern_index] << " to " << wfptr->get_adc(m_random_channels[m_pattern_index]);
      // Update the previous timestamp of the pattern generator
      m_pattern_generator_previous_ts = m_pattern_generator_current_ts;
    } // timestamp difference
  } // if not first superchunk
}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
 * */
void
WIB2FrameProcessor::timestamp_check(frameptr fp)
{

  uint16_t wib2_tick_difference = types::DUNEWIBSuperChunkTypeAdapter::expected_tick_difference;
  uint16_t wib2_superchunk_tick_difference = wib2_tick_difference * fp->get_num_frames();

  // If EMU data, emulate perfectly incrementing timestamp
  if (inherited::m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
    uint64_t ts_next = m_previous_ts + wib2_superchunk_tick_difference; // NOLINT(build/unsigned)
    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));            // NOLINT
    for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
      // auto wfh = const_cast<dunedaq::fddetdataformats::WIB2Header*>(wf->get_wib_header());
      wf->header.crate = m_crate_no;
      wf->header.slot = m_slot_no;
      wf->header.link = m_link; 
      wf->set_timestamp(ts_next);
      ts_next += wib2_tick_difference;
      wf++;
    }
  }

  // Acquire timestamp
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>(fp); // NOLINT
  m_current_ts = wfptr->get_timestamp();

  // Check timestamp
  if (m_current_ts - m_previous_ts != wib2_superchunk_tick_difference) {
    ++m_ts_error_ctr;
    m_error_registry->add_error("MISSING_FRAMES", readoutlibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + wib2_superchunk_tick_difference, m_current_ts));
    if (m_first_ts_missmatch) { // log once
      TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISSMATCH! -> | previous: " << std::to_string(m_previous_ts) << " current: " + std::to_string(m_current_ts);
      m_first_ts_missmatch = false;
    }
  }

  if (m_ts_error_ctr > 1000) {
    if (!m_problem_reported) {
      TLOG() << "*** Data Integrity ERROR *** Timestamp continuity is completely broken! "
             << "Something is wrong with the FE source or with the configuration!";
      m_problem_reported = true;
    }
  }

  m_previous_ts = m_current_ts;
  m_last_processed_daq_ts = m_current_ts;
}

/**
 * Pipeline Stage 2.: Do software TPG
 * */
void
WIB2FrameProcessor::find_hits(constframeptr fp, WIB2FrameHandler* frame_handler)
{
  if (!fp)
    return;
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIB2Frame*>((uint8_t*)fp); // NOLINT
  uint64_t timestamp = wfptr->get_timestamp();                                            // NOLINT(build/unsigned)

  // Frame expansion
  swtpg_wib2::MessageRegisters registers_array;
  int register_selection = frame_handler->get_registers_selector();
  expand_wib2_adcs(fp, &registers_array, register_selection);

  // Only for the first superchunk, create an offline register map
  if (frame_handler->first_hit) {
    frame_handler->register_channel_map = swtpg_wib2::get_register_to_offline_channel_map_wib2(wfptr, m_channel_map, register_selection);

    frame_handler->m_tpg_processing_info->setState(registers_array);

    m_det_id = wfptr->header.detector_id;
    if (wfptr->header.crate != m_crate_no || wfptr->header.slot != m_slot_no || wfptr->header.link != m_link) {
      ers::error(LinkMisconfiguration(ERS_HERE, wfptr->header.crate, wfptr->header.slot, wfptr->header.link, m_crate_no, m_slot_no, m_link));
    }
    // Debugging statements
    m_crate_no = wfptr->header.crate;
    m_slot_no = wfptr->header.slot;
    m_link = wfptr->header.link;
    //TLOG() << "Got first item link/crate/slot=" << m_link << "/" << m_crate_no << "/" << m_slot_no;

    // Add WIB2FrameHandler channel map to the common m_register_channels.
    // Populate the array taking into account the position of the register selector
    for (size_t i = 0; i < swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER; ++i) {
      m_register_channels[i + register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER] = frame_handler->register_channel_map.channel[i];

      // Set up a map of channels and number of TPs for monitoring/debug
      m_tp_channel_rate_map[frame_handler->register_channel_map.channel[i]] = 0;
    }

    //TLOG() << "Processed the first superchunk ";

    // Set first hit bool to false so that registration of channel map is not executed twice
    frame_handler->first_hit = false;

  } // end if (frame_handler->first_hit)

  // Execute the SWTPG algorithm
  frame_handler->m_tpg_processing_info->input = &registers_array;
  // Set the first word to "magic" indicating there is no hit, initially
  frame_handler->m_tpg_processing_info->output[0] = swtpg_wib2::MAGIC;

  if (m_tpg_algorithm == "SWTPG") {
    swtpg_wib2::process_window_avx2(*frame_handler->m_tpg_processing_info, register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER);
  } else if (m_tpg_algorithm == "AbsRS") {
    swtpg_wib2::process_window_rs_avx2(*frame_handler->m_tpg_processing_info, register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER);
  } else {
    throw TPGAlgorithmInexistent(ERS_HERE, "m_tpg_algo");
  }

  //GLM: avoid the tp_handler queue/thread
  process_swtpg_hits(frame_handler->m_tpg_processing_info->output, timestamp);
}

void
WIB2FrameProcessor::process_swtpg_hits(uint16_t* primfind_it, dunedaq::daqdataformats::timestamp_t timestamp)
{

  constexpr int clocksPerTPCTick = types::DUNEWIBSuperChunkTypeAdapter::expected_tick_difference;;

  uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; // NOLINT(build/unsigned)
  unsigned int nhits = 0;

  while (*primfind_it != swtpg_wib2::MAGIC) {
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

    // Now that we have all the register values in local
    // variables, loop over the register index (ie, channel) and
    // find the channels which actually had a hit, as indicated by
    // nonzero value of hit_charge
    for (int i = 0; i < 16; ++i) {
      if (hit_charge[i] && chan[i] != swtpg_wib2::MAGIC) {

        uint64_t tp_t_begin =                                                           // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_tover[i])); // NOLINT(build/unsigned)
        uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);         // NOLINT(build/unsigned)
        //TLOG() << "Hit start " << tp_t_begin << ", end: " << tp_t_end << ", online channel: " << chan[i];
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
          tp.tp.time_peak = (tp_t_begin + tp_t_end) / 2;
          tp.tp.time_over_threshold = int64_t(hit_tover[i]) * clocksPerTPCTick;
          tp.tp.channel = offline_channel;
          tp.tp.adc_integral = hit_charge[i];
          tp.tp.adc_peak = hit_charge[i] / 20;
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
  m_swtpg_hits_count += nhits;
  return;
}

} // namespace fdreadoutlibs
} // namespace dunedaq
