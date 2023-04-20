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
#include "detdataformats/wib2/WIB2Frame.hpp"


#include "fdreadoutlibs/DUNEWIBSuperChunkTypeAdapter.hpp"

#include "fdreadoutlibs/wib2/WIB2TPHandler.hpp"
#include "rcif/cmd/Nljs.hpp"
#include "trigger/TPSet.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

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
#include <pthread.h>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using dunedaq::readoutlibs::logging::TLVL_BOOKKEEPING;
using dunedaq::readoutlibs::logging::TLVL_TAKE_NOTE;

namespace dunedaq {
namespace fdreadoutlibs {

WIB2FrameHandler::WIB2FrameHandler(int register_selector_params, iomanager::FollyMPMCQueue<uint16_t*>& dest_queue)
  : m_dest_queue_frame_handler(dest_queue)
{
  m_register_selector = register_selector_params;
}
WIB2FrameHandler::~WIB2FrameHandler()
{
  if (m_tpg_taps_p) {
    delete[] m_tpg_taps_p;
  }
}

int
WIB2FrameHandler::get_registers_selector()
{
  return m_register_selector;
}

void
WIB2FrameHandler::reset()
{
  delete[] m_tpg_taps_p;
  m_tpg_taps_p = nullptr;
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

  m_tpg_processing_info = std::make_unique<swtpg_wib2::ProcessingInfo<swtpg_wib2::NUM_REGISTERS_PER_FRAME>>(nullptr,
                                                                                                            swtpg_wib2::FRAMES_PER_MSG,
                                                                                                            0,
                                                                                                            swtpg_wib2::NUM_REGISTERS_PER_FRAME,
                                                                                                            nullptr,
                                                                                                            m_tpg_taps_p,
                                                                                                            (uint8_t)m_tpg_taps.size(), // NOLINT(build/unsigned)
                                                                                                            m_tpg_tap_exponent,
                                                                                                            m_tpg_threshold,
                                                                                                            0,
                                                                                                            0);
}

// Pop one destination ptr for the frame handler
uint16_t*
WIB2FrameHandler::get_primfind_dest()
{
  uint16_t* primfind_dest;

  while (!m_dest_queue_frame_handler.try_pop(primfind_dest, std::chrono::milliseconds(0))) {
    std::this_thread::sleep_for(std::chrono::microseconds(10));
  }
  return primfind_dest;
}

WIB2FrameProcessor::WIB2FrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
  : TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>(error_registry)
  , m_sw_tpg_enabled(false)
  , m_add_hits_tphandler_thread_should_run(false)
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

    rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
    m_tphandler->set_run_number(start_params.run);

    m_tphandler->reset();

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
    TLOG() << "Ran for " << runtime << "ms.";
  }
}

void
WIB2FrameProcessor::init(const nlohmann::json& args)
{
  try {
    auto queue_index = appfwk::connection_index(args, {});
    if (queue_index.find("tp_out") != queue_index.end()) {
      m_tp_sink = get_iom_sender<types::TriggerPrimitiveTypeAdapter>(queue_index["tp_out"]);
    }
    if (queue_index.find("tpset_out") != queue_index.end()) {
      m_tpset_sink = get_iom_sender<trigger::TPSet>(queue_index["tpset_out"]);
    }
  } catch (const ers::Issue& excpt) {
    throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
  }
}

void
WIB2FrameProcessor::conf(const nlohmann::json& cfg)
{
  // Populate the queue of primfind destinations
  for (size_t i = 0; i < m_capacity_dest_queue; ++i) {
    uint16_t* dest = new uint16_t[100000];
    m_dest_queue.push(std::move(dest), std::chrono::milliseconds(0));
  }

  auto config = cfg["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
  m_sourceid.id = config.source_id;
  m_sourceid.subsystem = types::DUNEWIBSuperChunkTypeAdapter::subsystem;
  m_tpg_algorithm = config.software_tpg_algorithm;
  TLOG() << "Selected software TPG algorithm: " << m_tpg_algorithm;

  m_channel_mask_vec = config.software_tpg_channel_mask;
  // Converting the input vector of channels masks into an std::set
  // AAA: The set provides faster look up than a std::vector
  m_channel_mask_set.insert(m_channel_mask_vec.begin(), m_channel_mask_vec.end());
  for (int el : m_channel_mask_set) {
    TLOG() << "Software TPG channel mask: " << el;
  }

  m_tpg_threshold_selected = config.software_tpg_threshold;
  TLOG() << "Selected threshold value: " << m_tpg_threshold_selected;

  if (config.enable_software_tpg) {
    m_sw_tpg_enabled = true;

    m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

    daqdataformats::SourceID tpset_sourceid;
    tpset_sourceid.id = config.tpset_sourceid;
    tpset_sourceid.subsystem = daqdataformats::SourceID::Subsystem::kTrigger;

    m_tphandler.reset(new WIB2TPHandler(m_tp_sink, m_tpset_sink, config.tp_timeout, config.tpset_window_size, tpset_sourceid));

    TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::add_postprocess_task(std::bind(&WIB2FrameProcessor::find_hits, this, std::placeholders::_1, m_wib2_frame_handler.get()));

    TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::add_postprocess_task(std::bind(&WIB2FrameProcessor::find_hits, this, std::placeholders::_1, m_wib2_frame_handler_second_half.get()));

    // Launch the thread for adding hits to tphandler
    m_add_hits_tphandler_thread_should_run.store(true);
    m_add_hits_tphandler_thread = std::thread(&WIB2FrameProcessor::add_hits_to_tphandler, this);
    TLOG() << "Launched thread for adding hits to tphandler";
  }

  // Setup pre-processing pipeline
  TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::add_preprocess_task(std::bind(&WIB2FrameProcessor::timestamp_check, this, std::placeholders::_1));

  TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::conf(cfg);
}

void
WIB2FrameProcessor::scrap(const nlohmann::json& args)
{
  if (m_sw_tpg_enabled) {
    TLOG() << "Waiting to join add_hits_tphandler_thread";
    m_add_hits_tphandler_thread_should_run.store(false);
    m_add_hits_tphandler_thread.join();
    TLOG() << "add_hits_tphandler_thread joined";
    m_tphandler.reset();

    // Pop and delete all the elements of the destination queues
    while (m_dest_queue.can_pop()) {
      uint16_t* dest;
      m_dest_queue.pop(dest, std::chrono::milliseconds(0));
      delete[] dest;
    }
  }
  TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::scrap(args);
}

void
WIB2FrameProcessor::get_info(opmonlib::InfoCollector& ci, int level)
{
  readoutlibs::readoutinfo::RawDataProcessorInfo info;

  if (m_tphandler != nullptr) {
    info.num_tps_sent = m_tphandler->get_and_reset_num_sent_tps();
    info.num_tpsets_sent = m_tphandler->get_and_reset_num_sent_tpsets();
    info.num_tps_dropped = m_tps_dropped.exchange(0);
  }

  auto now = std::chrono::high_resolution_clock::now();
  if (m_sw_tpg_enabled) {
    int new_hits = m_swtpg_hits_count.exchange(0);
    int new_tps = m_new_tps.exchange(0);
    double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
    TLOG_DEBUG(TLVL_BOOKKEEPING) << "Total new hits: " << new_hits << " new TPs: " << new_tps;
    info.rate_tp_hits = new_hits / seconds / 1000.;

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

  readoutlibs::TaskRawDataProcessorModel<types::DUNEWIBSuperChunkTypeAdapter>::get_info(ci, level);
  ci.add(info);
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
      // auto wfh = const_cast<dunedaq::detdataformats::wib2::WIB2Header*>(wf->get_wib_header());
      wf->set_timestamp(ts_next);
      ts_next += wib2_tick_difference;
      wf++;
    }
  }

  // Acquire timestamp
  auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>(fp); // NOLINT
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

  auto wfptr = reinterpret_cast<dunedaq::detdataformats::wib2::WIB2Frame*>((uint8_t*)fp); // NOLINT
  uint64_t timestamp = wfptr->get_timestamp();                                            // NOLINT(build/unsigned)

  // Frame expansion
  swtpg_wib2::MessageRegisters registers_array;
  int register_selection = frame_handler->get_registers_selector();
  expand_wib2_adcs(fp, &registers_array, register_selection);

  // Only for the first superchunk, create an offline register map
  if (frame_handler->first_hit) {

    // Print thread ID and PID of the executing thread
    std::thread::id thread_id = std::this_thread::get_id();
    pid_t tid;
    tid = syscall(SYS_gettid);
    TLOG_DEBUG(TLVL_BOOKKEEPING) << " Thread ID " << thread_id << " PID " << tid;

    frame_handler->register_channel_map = swtpg_wib2::get_register_to_offline_channel_map_wib2(wfptr, m_channel_map, register_selection);

    frame_handler->m_tpg_processing_info->setState(registers_array);

    // Debugging statements
    m_link = wfptr->header.link;
    m_crate_no = wfptr->header.crate;
    m_slot_no = wfptr->header.slot;
    TLOG() << "Got first item, link/crate/slot=" << m_link << "/" << m_crate_no << "/" << m_slot_no;

    // Add WIB2FrameHandler channel map to the common m_register_channels.
    // Populate the array taking into account the position of the register selector
    for (size_t i = 0; i < swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER; ++i) {
      m_register_channels[i + register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER] = frame_handler->register_channel_map.channel[i];

      // Set up a map of channels and number of TPs for monitoring/debug
      m_tp_channel_rate_map[frame_handler->register_channel_map.channel[i]] = 0;
    }

    TLOG() << "Processed the first superchunk ";

    // Set first hit bool to false so that registration of channel map is not executed twice
    frame_handler->first_hit = false;

  } // end if (frame_handler->first_hit)

  // Execute the SWTPG algorithm
  frame_handler->m_tpg_processing_info->input = &registers_array;
  uint16_t* destination_ptr = frame_handler->get_primfind_dest();
  *destination_ptr = swtpg_wib2::MAGIC;
  frame_handler->m_tpg_processing_info->output = destination_ptr;

  if (m_tpg_algorithm == "SWTPG") {
    swtpg_wib2::process_window_avx2(*frame_handler->m_tpg_processing_info, register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER);
  } else if (m_tpg_algorithm == "AbsRS") {
    swtpg_wib2::process_window_rs_avx2(*frame_handler->m_tpg_processing_info, register_selection * swtpg_wib2::NUM_REGISTERS_PER_FRAME * swtpg_wib2::SAMPLES_PER_REGISTER);
  } else {
    throw TPGAlgorithmInexistent(ERS_HERE, "m_tpg_algo");
  }

  // Insert output of the AVX processing into the swtpg_output
  swtpg_output swtpg_processing_result = { destination_ptr, timestamp };

  // Push to the MPMC tphandler queue only if it's possible, else drop the TPs.
  if (!m_tphandler_queue.try_push(std::move(swtpg_processing_result), std::chrono::milliseconds(0))) {
    // we're going to loose these hits
    ers::warning(TPHandlerBacklog(ERS_HERE, m_sourceid.id));
  }
}

unsigned int
WIB2FrameProcessor::process_swtpg_hits(uint16_t* primfind_it, timestamp_t timestamp)
{

  constexpr int clocksPerTPCTick = 32;

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

        // This channel had a hit ending here, so we can create and output the hit here
        const uint16_t offline_channel = m_register_channels[chan[i]];
        if (m_channel_mask_set.find(offline_channel) == m_channel_mask_set.end()) {

          uint64_t tp_t_begin =                                                           // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - int64_t(hit_tover[i])); // NOLINT(build/unsigned)
          uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);         // NOLINT(build/unsigned)

          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          // sed -n -e 's/.*Hit: \(.*\) \(.*\).*/\1 \2/p' log.txt  > hits.txt
          //
          // TLOG() << "Hit: " << tp_t_begin << " " << offline_channel;

          triggeralgs::TriggerPrimitive trigprim;
          trigprim.time_start = tp_t_begin;
          trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
          trigprim.time_over_threshold = int64_t(hit_tover[i]) * clocksPerTPCTick;
          trigprim.channel = offline_channel;
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_charge[i] / 20;
          trigprim.detid = m_link; // TODO: convert crate/slot/link to SourceID Roland Sipos rsipos@cern.ch July-22-2021
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kTPCDefault;
          trigprim.version = 1;

          if (!m_tphandler->add_tp(trigprim, timestamp)) {
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
  return nhits;
}

// Push destination ptr to the common queue
void
WIB2FrameProcessor::free_primfind_dest(uint16_t* dest)
{
  m_dest_queue.push(std::move(dest), std::chrono::milliseconds(0));
}

// Function for the TPHandler threads.
// Reads the m_tpthread_queue and then process the hits
void
WIB2FrameProcessor::add_hits_to_tphandler()
{

  std::stringstream thread_name;
  thread_name << "tphandler-" << m_sourceid.id;
  pthread_setname_np(pthread_self(), thread_name.str().c_str());

  while (m_add_hits_tphandler_thread_should_run.load()) {
    swtpg_output result_from_swtpg;
    while (m_tphandler_queue.can_pop()) {
      if (m_tphandler_queue.try_pop(result_from_swtpg, std::chrono::milliseconds(0))) {

        // Process the trigger primitve
        unsigned int nhits = process_swtpg_hits(result_from_swtpg.output_location, result_from_swtpg.timestamp);

        m_swtpg_hits_count += nhits;
        m_tphandler->try_sending_tpsets(result_from_swtpg.timestamp);

        // After sending the TPset, add back the dest pointer to the queue
        free_primfind_dest(result_from_swtpg.output_location);
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
}

} // namespace fdreadoutlibs
} // namespace dunedaq
