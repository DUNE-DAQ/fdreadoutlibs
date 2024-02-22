/**
 * @file WIBFrameProcessor.hpp WIB specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBFRAMEPROCESSOR_HPP_

#include "appfwk/DAQModuleHelper.hpp"
#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "readoutlibs/FrameErrorRegistry.hpp"
#include "readoutlibs/ReadoutIssues.hpp"
#include "readoutlibs/ReadoutLogging.hpp"
#include "readoutlibs/models/IterableQueueModel.hpp"
#include "readoutlibs/models/TaskRawDataProcessorModel.hpp"
#include "readoutlibs/readoutconfig/Nljs.hpp"
#include "readoutlibs/readoutinfo/InfoNljs.hpp"
#include "readoutlibs/utils/ReusableThread.hpp"

#include "detchannelmaps/TPCChannelMap.hpp"
#include "fddetdataformats/WIBFrame.hpp"


#include "fdreadoutlibs/TriggerPrimitiveTypeAdapter.hpp"

#include "fdreadoutlibs/ProtoWIBSuperChunkTypeAdapter.hpp"
#include "fdreadoutlibs/wib/WIBTPHandler.hpp"
#include "rcif/cmd/Nljs.hpp"
#include "trigger/TPSet.hpp"
#include "triggeralgs/TriggerPrimitive.hpp"

#include "tpg/DesignFIR.hpp"
#include "tpg/FrameExpand.hpp"
#include "tpg/ProcessAVX2.hpp"
#include "tpg/ProcessingInfo.hpp"
#include "tpg/RegisterToChannelNumber.hpp"
#include "tpg/TPGConstants.hpp"

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
using dunedaq::readoutlibs::logging::TLVL_FRAME_RECEIVED;

namespace {
enum CollectionOrInduction {
  kCollection,
  kInduction
};
}

namespace dunedaq {

namespace fdreadoutlibs {

class WIBFrameProcessor : public readoutlibs::TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>
{

public:
  using inherited = readoutlibs::TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>;
  using frameptr = types::ProtoWIBSuperChunkTypeAdapter*;
  using constframeptr = const types::ProtoWIBSuperChunkTypeAdapter*;
  using wibframeptr = dunedaq::fddetdataformats::WIBFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)

  // Channel map function type
  typedef int (*chan_map_fn_t)(int);

  explicit WIBFrameProcessor(std::unique_ptr<readoutlibs::FrameErrorRegistry>& error_registry)
    : TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>(error_registry)
    , m_sw_tpg_enabled(false)
    , m_ind_thread_should_run(false)
    , m_coll_primfind_dest(nullptr)
    , m_coll_taps_p(nullptr)
    , m_ind_primfind_dest(nullptr)
    , m_ind_taps_p(nullptr)
  {}

  ~WIBFrameProcessor()
  {
    if (m_coll_taps_p) {
      delete[] m_coll_taps_p;
    }
    if (m_coll_primfind_dest) {
      delete[] m_coll_primfind_dest;
    }
    if (m_ind_taps_p) {
      delete[] m_ind_taps_p;
    }
    if (m_ind_primfind_dest) {
      delete[] m_ind_primfind_dest;
    }
  }

  void start(const nlohmann::json& args) override
  {
    // Reset software TPG resources
    if (m_sw_tpg_enabled) {

      rcif::cmd::StartParams start_params = args.get<rcif::cmd::StartParams>();
      m_tphandler->set_run_number(start_params.run);
  
      m_tphandler->reset();
      m_tps_suppressed_too_long = 0;

      m_coll_taps = swtpg::firwin_int(7, 0.1, m_coll_multiplier);
      m_coll_taps.push_back(0);
      m_ind_taps = swtpg::firwin_int(7, 0.1, m_ind_multiplier);
      m_ind_taps.push_back(0);

      if (m_coll_taps_p == nullptr) {
        m_coll_taps_p = new int16_t[m_coll_taps.size()];
      }
      for (size_t i = 0; i < m_coll_taps.size(); ++i) {
        m_coll_taps_p[i] = m_coll_taps[i];
      }

      if (m_ind_taps_p == nullptr) {
        m_ind_taps_p = new int16_t[m_ind_taps.size()];
      }
      for (size_t i = 0; i < m_ind_taps.size(); ++i) {
        m_ind_taps_p[i] = m_ind_taps[i];
      }

      // Temporary place to stash the hits
      if (m_coll_primfind_dest == nullptr) {
        m_coll_primfind_dest = new uint16_t[100000]; // NOLINT(build/unsigned)
      }
      if (m_ind_primfind_dest == nullptr) {
        m_ind_primfind_dest = new uint16_t[100000]; // NOLINT(build/unsigned)
      }

      TLOG() << "COLL TAPS SIZE: " << m_coll_taps.size() << " threshold:" << m_coll_threshold
             << " exponent:" << m_coll_tap_exponent;

      m_coll_tpg_pi = std::make_unique<swtpg::ProcessingInfo<swtpg::COLLECTION_REGISTERS_PER_FRAME>>(
        nullptr,
        swtpg::FRAMES_PER_MSG,
        0,
        swtpg::COLLECTION_REGISTERS_PER_FRAME,
        m_coll_primfind_dest,
        m_coll_taps_p,
        (uint8_t)m_coll_taps.size(), // NOLINT(build/unsigned)
        m_coll_tap_exponent,
        m_coll_threshold,
        0,
        0);

      m_ind_tpg_pi = std::make_unique<swtpg::ProcessingInfo<swtpg::INDUCTION_REGISTERS_PER_FRAME>>(
        nullptr,
        swtpg::FRAMES_PER_MSG,
        0,
        swtpg::INDUCTION_REGISTERS_PER_FRAME,
        m_ind_primfind_dest,
        m_ind_taps_p,
        (uint8_t)m_ind_taps.size(), // NOLINT(build/unsigned)
        m_ind_tap_exponent,
        m_ind_threshold,
        0,
        0);

      //m_induction_thread = std::thread(&WIBFrameProcessor::find_induction_hits_thread, this);
    } // end if(m_sw_tpg_enabled)

    // Reset timestamp check
    m_previous_ts = 0;
    m_current_ts = 0;
    m_first_ts_missmatch = true;
    m_problem_reported = false;
    m_ts_error_ctr = 0;

    // Reset stats
    m_first_coll = true;
    m_t0 = std::chrono::high_resolution_clock::now();
    m_new_hits = 0;
    m_new_tps = 0;
    m_coll_hits_count.exchange(0);
    m_frame_error_count = 0;
    m_frames_processed = 0;

    inherited::start(args);
  }

  void stop(const nlohmann::json& args) override
  {
    inherited::stop(args);
    if (m_sw_tpg_enabled) {
      // Make temp. buffers reusable on next start.
      if (m_coll_taps_p) {
        delete[] m_coll_taps_p;
        m_coll_taps_p = nullptr;
      }
      if (m_coll_primfind_dest) {
        delete[] m_coll_primfind_dest;
        m_coll_primfind_dest = nullptr;
      }
      if (m_ind_taps_p) {
        delete[] m_ind_taps_p;
        m_ind_taps_p = nullptr;
      }
      if (m_ind_primfind_dest) {
        delete[] m_ind_primfind_dest;
        m_ind_primfind_dest = nullptr;
      }
      auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_t0).count();
      TLOG() << "Ran for " << runtime << "ms. Found " << m_num_hits_coll << " collection hits and " << m_num_hits_ind << " induction hits";
    }

  }

  void init(const nlohmann::json& args) override
  {
    try {
      auto queue_index = appfwk::connection_index(args, {});
      if (queue_index.find("tp_out") != queue_index.end()) {
        m_tp_sink = get_iom_sender<types::TriggerPrimitiveTypeAdapter>(queue_index["tp_out"]);
      }
      if (queue_index.find("tpset_out") != queue_index.end()) {
        m_tpset_sink = get_iom_sender<trigger::TPSet>(queue_index["tpset_out"]);
      }
      m_err_frame_sink = get_iom_sender<fddetdataformats::WIBFrame>(queue_index["errored_frames"]);
    } catch (const ers::Issue& excpt) {
      throw readoutlibs::ResourceQueueError(ERS_HERE, "tp queue", "DefaultRequestHandlerModel", excpt);
    }
  }

  void conf(const nlohmann::json& cfg) override
  {
    auto config = cfg["rawdataprocessorconf"].get<readoutlibs::readoutconfig::RawDataProcessorConf>();
    m_sourceid.id = config.source_id;
    m_sourceid.subsystem = types::ProtoWIBSuperChunkTypeAdapter::subsystem;
    m_error_counter_threshold = config.error_counter_threshold;
    m_error_reset_freq = config.error_reset_freq;
    m_clock_frequency = config.clock_speed_hz;

    if (config.enable_software_tpg) {
      m_sw_tpg_enabled = true;

      m_channel_map = dunedaq::detchannelmaps::make_map(config.channel_map_name);

      daqdataformats::SourceID tpset_sourceid;
      tpset_sourceid.id = config.tpset_sourceid;
      tpset_sourceid.subsystem = daqdataformats::SourceID::Subsystem::kTrigger;
      m_tphandler.reset(
        new WIBTPHandler(*m_tp_sink, *m_tpset_sink, config.tp_timeout, config.tpset_window_size, tpset_sourceid));

      // m_induction_items_to_process = std::make_unique<readoutlibs::IterableQueueModel<InductionItemToProcess>>(
      //   200000, false, 0, true, 64); // 64 byte aligned

      // Setup parallel post-processing
      TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::add_postprocess_task(
        std::bind(&WIBFrameProcessor::find_collection_hits, this, std::placeholders::_1));

      // start the thread for induction hit finding
      TLOG() << "Launch induction hit finding thread"; 
      m_ind_thread_should_run.store(true);
      m_induction_thread = std::thread(&WIBFrameProcessor::find_induction_hits_thread, this);
    }

    // Setup pre-processing pipeline
    TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::add_preprocess_task(
      std::bind(&WIBFrameProcessor::timestamp_check, this, std::placeholders::_1));
    //TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::add_preprocess_task(
      //std::bind(&WIBFrameProcessor::frame_error_check, this, std::placeholders::_1));

    TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::conf(cfg);
  }

  void scrap(const nlohmann::json& args) override
  {
   if(m_sw_tpg_enabled) {	  
      TLOG() << "Waiting to join induction thread";
      m_ind_thread_should_run.store(false);
      m_induction_thread.join();
      TLOG() << "Induction thread joined";
   }
      m_tphandler.reset();

    TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::scrap(args);
  }

  void get_info(opmonlib::InfoCollector& ci, int level)
  {
    readoutlibs::readoutinfo::RawDataProcessorInfo info;

    if (m_tphandler != nullptr) {
      info.num_tps_sent = m_tphandler->get_and_reset_num_sent_tps();
      info.num_tpsets_sent = m_tphandler->get_and_reset_num_sent_tpsets();
      info.num_tps_suppressed_too_long = m_tps_suppressed_too_long.exchange(0);
    }
    info.num_frame_errors = m_frame_error_count.exchange(0);

    auto now = std::chrono::high_resolution_clock::now();
    if (m_sw_tpg_enabled) {
      int new_hits = m_coll_hits_count.exchange(0);
      int new_tps = m_num_tps_pushed.exchange(0);
      double seconds = std::chrono::duration_cast<std::chrono::microseconds>(now - m_t0).count() / 1000000.;
      TLOG_DEBUG(TLVL_TAKE_NOTE) << "Hit rate: " << std::to_string(new_hits / seconds / 1000.) << " [kHz]";
      TLOG_DEBUG(TLVL_TAKE_NOTE) << "Total new hits: " << new_hits << " new pushes: " << new_tps;
      info.rate_tp_hits = new_hits / seconds / 1000.;
    }
    m_t0 = now;

    readoutlibs::TaskRawDataProcessorModel<types::ProtoWIBSuperChunkTypeAdapter>::get_info(ci, level);
    ci.add(info);
  }

protected:
  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };

  struct InductionItemToProcess
  {
    // Horribly, `registers` has to be the first item in the
    // struct, because the first item in the queue has to be
    // correctly aligned, and we're going to put this in an
    // AlignedProducerConsumerQueue, which aligns the *starts* of
    // the contained objects to 64-byte boundaries, not any later
    // items
    swtpg::MessageRegistersInduction registers;
    uint64_t timestamp; // NOLINT(build/unsigned)

    static constexpr uint64_t END_OF_MESSAGES = UINT64_MAX; // NOLINT(build/unsigned)
  };

  void postprocess_example(const types::ProtoWIBSuperChunkTypeAdapter* fp)
  {
    TLOG() << "Postprocessing: " << fp->get_first_timestamp();
  }

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
   * */
  void timestamp_check(frameptr fp)
  {
    // If EMU data, emulate perfectly incrementing timestamp
    if (inherited::m_emulator_mode) {                           // emulate perfectly incrementing timestamp
      uint64_t ts_next = m_previous_ts + 300;                   // NOLINT(build/unsigned)
      auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp));  // NOLINT
      for (unsigned int i = 0; i < fp->get_num_frames(); ++i) { // NOLINT(build/unsigned)
        auto wfh = const_cast<dunedaq::fddetdataformats::WIBHeader*>(wf->get_wib_header());
        wfh->set_timestamp(ts_next);
        ts_next += 25;
        wf++;
      }
    }

    // Acquire timestamp
    auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBFrame*>(fp); // NOLINT
    auto wfhdrptr = wfptr->get_wib_header();
    m_current_ts = wfhdrptr->get_timestamp();
    TLOG_DEBUG(TLVL_FRAME_RECEIVED) << "Received ProtoWIB frame timestamp value of " << m_current_ts << " ticks (..." << std::fixed << std::setprecision(8) << (static_cast<double>(m_current_ts % (m_clock_frequency*1000)) / static_cast<double>(m_clock_frequency)) << " sec), crate " << wfhdrptr->crate_no << ", slot " << wfhdrptr->slot_no; // NOLINT

    // Check timestamp
    if (m_current_ts - m_previous_ts != 300) {
      ++m_ts_error_ctr;
      m_error_registry->add_error("MISSING_FRAMES",
                                  readoutlibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + 300, m_current_ts));
      if (m_first_ts_missmatch) { // log once
        TLOG_DEBUG(TLVL_BOOKKEEPING) << "First timestamp MISSMATCH! -> | previous: " << std::to_string(m_previous_ts)
                                     << " current: " + std::to_string(m_current_ts);
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
   * Pipeline Stage 2.: Check WIB headers for error flags
   * */
  void frame_error_check(frameptr fp)
  {
    if (!fp)
      return;

    auto wf = reinterpret_cast<wibframeptr>(((uint8_t*)fp)); // NOLINT
    for (size_t i = 0; i < fp->get_num_frames(); ++i) {
      if (m_frames_processed % 10000 == 0) {
        for (int i = 0; i < m_num_frame_error_bits; ++i) {
          if (m_error_occurrence_counters[i])
            m_error_occurrence_counters[i]--;
        }
      }

      auto wfh = const_cast<dunedaq::fddetdataformats::WIBHeader*>(wf->get_wib_header());
      if (wfh->wib_errors) {
        m_frame_error_count += std::bitset<16>(wfh->wib_errors).count();
      }

      m_current_frame_pushed = false;
      for (int j = 0; j < m_num_frame_error_bits; ++j) {
        if (wfh->wib_errors & (1 << j)) {
          if (m_error_occurrence_counters[j] < m_error_counter_threshold) {
            m_error_occurrence_counters[j]++;
            if (!m_current_frame_pushed) {
              try {
                  dunedaq::fddetdataformats::WIBFrame wf_copy(*wf);
                m_err_frame_sink->send(std::move(wf_copy), std::chrono::milliseconds(10));
                m_current_frame_pushed = true;
              } catch (const ers::Issue& excpt) {
                ers::warning(readoutlibs::CannotWriteToQueue(ERS_HERE, m_sourceid, "Errored frame queue", excpt));
              }
            }
          }
        }
      }
      wf++;
      m_frames_processed++;
    }
  }

  /**
   * Pipeline Stage 3.: Do software TPG
   * */
  void find_collection_hits(constframeptr fp)
  {
    if (!fp)
      return;

    auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::WIBFrame*>((uint8_t*)fp); // NOLINT
    uint64_t timestamp = wfptr->get_wib_header()->get_timestamp();                        // NOLINT(build/unsigned)

    // First, expand the frame ADCs into 16-bit values in AVX2
    // registers. They're split into "collection" and "induction"
    // channels: this classification is correct for ProtoDUNE-I data,
    // but not for the VD coldbox, where the channel map is very
    // different. In the VD coldbox case, all that the division into
    // "collection" and "induction" registers does is split the
    // channels into two groups so we can process some of them on one
    // thread, and some on another, since a single thread can't keep
    // up with all channels
    swtpg::MessageRegistersCollection collection_registers;
    InductionItemToProcess ind_item;
    expand_message_adcs_inplace(fp, &collection_registers, &ind_item.registers);

    if (m_first_coll) {
      m_register_channel_map = swtpg::get_register_to_offline_channel_map(wfptr, m_channel_map);

      m_coll_tpg_pi->setState(collection_registers);

      m_fiber_no = wfptr->get_wib_header()->fiber_no;
      m_crate_no = wfptr->get_wib_header()->crate_no;
      m_slot_no = wfptr->get_wib_header()->slot_no;

      TLOG() << "Got first item, fiber/crate/slot=" << m_fiber_no << "/" << m_crate_no << "/" << m_slot_no;

      std::stringstream ss;
      ss << "Collection channels are:\n";
      for(size_t i=0; i<swtpg::COLLECTION_REGISTERS_PER_FRAME*swtpg::SAMPLES_PER_REGISTER; ++i){
        ss << i << "\t" << m_register_channel_map.collection[i] << "\n";
      }
      TLOG_DEBUG(2) << ss.str();

      std::stringstream ss2;
      ss2 << "Induction channels are:\n";
      for(size_t i=0; i<swtpg::INDUCTION_REGISTERS_PER_FRAME*swtpg::SAMPLES_PER_REGISTER; ++i){
        ss2 << i << "\t" << m_register_channel_map.induction[i] << "\n";
      }
      TLOG_DEBUG(2) << ss2.str();

    } // end if (m_first_coll)

    // Signal to the induction thread that there's an item ready
    m_induction_item_to_process = &ind_item;
    m_induction_item_ready.store(true);

    // Find the hits in the "collection" registers
    m_coll_tpg_pi->input = &collection_registers;
    *m_coll_primfind_dest = swtpg::MAGIC;
    swtpg::process_window_avx2(*m_coll_tpg_pi);

    unsigned int nhits = add_hits_to_tphandler(m_coll_primfind_dest, timestamp, kCollection);

    // if (nhits > 0) {
    // TLOG() << "NON null hits: " << nhits << " for ts: " << timestamp;
    // TLOG() << *wfptr;
    //}

    m_num_hits_coll += nhits;
    m_coll_hits_count += nhits;

    if (m_first_coll) {
      TLOG() << "Total hits in first superchunk: " << nhits;
      m_first_coll = false;
    }

    // Wait for the induction item to be done. We have to spin here,
    // and not sleep, because we only have 6 microseconds to process
    // each superchunk. It appears that anything that makes a system
    // call or puts the thread to sleep causes too much latency
    while (m_induction_item_ready.load()) {
      _mm_pause();
      // std::this_thread::sleep_for(std::chrono::microseconds(100));
    }

    m_num_hits_ind += add_hits_to_tphandler(m_ind_primfind_dest, timestamp, kInduction);

    m_tphandler->try_sending_tpsets(timestamp);
  }

  void find_induction_hits(InductionItemToProcess* induction_item_to_process)
  {
    if (m_first_ind) {
      m_ind_tpg_pi->setState(induction_item_to_process->registers);
      TLOG() << "Got first item, fiber/crate/slot=" << m_fiber_no << "/" << m_crate_no << "/" << m_slot_no;
    }

    m_ind_tpg_pi->input = &induction_item_to_process->registers;
    *m_ind_primfind_dest = swtpg::MAGIC;
    swtpg::process_window_avx2(*m_ind_tpg_pi);

    m_first_ind = false;

  }

  // Stage: induction hit finding port
  void find_induction_hits_thread()
  {
    std::stringstream thread_name;
    thread_name << "ind-hits-" << m_sourceid.id;
    pthread_setname_np(pthread_self(), thread_name.str().c_str());

    size_t n_items=0;
    while (m_ind_thread_should_run.load()) {
      // There must be a nicer way to write this
      while (!m_induction_item_ready.load()) {
        // PAR 2022-03-16 Empirically, we can't sleep here, or wait on
        // a condition variable: that causes everything to back up. I
        // think the reason is that we're processing superchunks
        // without any buffering, so we have to be done in 6
        // microseconds. When we put the thread to sleep, we don't
        // wake up quickly enough, and problems ensue.
        //
        // It would probably be nicer to have the collection and
        // induction threads talk to each other via a queue, so we get
        // buffering and can relax the latency requirement, but then
        // we would have to think carefully about how we pass hits to
        // m_tp_handler. So we do it this way and spin-wait

        // std::this_thread::sleep_for(std::chrono::microseconds(1));
        _mm_pause();
        if(!m_ind_thread_should_run.load()) break;
      }
      if(!m_ind_thread_should_run.load()) break;

      find_induction_hits(m_induction_item_to_process);

      // Signal back to the collection thread that we're done
      m_induction_item_ready.store(false);
      ++n_items;
    }
    // Make sure this gets set so the collection thread isn't waiting forever at stop
    m_induction_item_ready.store(false);

    TLOG() << "Induction hit-finding thread stopping after processing " << n_items << " frames";
  }

  unsigned int add_hits_to_tphandler(uint16_t* primfind_it, timestamp_t timestamp, CollectionOrInduction coll_or_ind)
  {
    constexpr int clocksPerTPCTick = 25;

    uint16_t chan[16], hit_end[16], hit_charge[16], hit_tover[16]; // NOLINT(build/unsigned)
    unsigned int nhits = 0;

    // process_window_avx2 stores its output in the buffer pointed to
    // by m_coll_primfind_dest in a (necessarily) complicated way: for
    // every set of 16 channels (one AVX2 register) that has at least
    // one hit which ends at this tick, the full 16-channel registers
    // of channel number, hit end time, hit charge and hit t-o-t are
    // stored. This is done for each of the (6 collection registers
    // per tick) x (12 ticks per superchunk), and the end of valid
    // hits is indicated by the presence of the value "MAGIC" (defined
    // in TPGConstants.h).
    //
    // Since not all channels in a register will have hits ending at
    // this tick, we look at the stored hit_charge to determine
    // whether this channel in the register actually had a hit ending
    // in it: for channels which *did* have a hit ending, the value of
    // hit_charge is nonzero.
    while (*primfind_it != swtpg::MAGIC) {
      // First, get all of the register values (including those with no hit) into local variables
      for (int i = 0; i < 16; ++i) {
        chan[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_end[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_charge[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }
      for (int i = 0; i < 16; ++i) {
        hit_tover[i] = *primfind_it++; // NOLINT(runtime/increment_decrement)
      }

      // Now that we have all the register values in local
      // variables, loop over the register index (ie, channel) and
      // find the channels which actually had a hit, as indicated by
      // nonzero value of hit_charge
      for (int i = 0; i < 16; ++i) {
        if (hit_charge[i] && chan[i] != swtpg::MAGIC) {
          // This channel had a hit ending here, so we can create and output the hit here
          const uint16_t offline_channel = (coll_or_ind == kCollection) ?
            m_register_channel_map.collection[chan[i]] : m_register_channel_map.induction[chan[i]];
          uint64_t tp_t_begin =                                                        // NOLINT(build/unsigned)
            timestamp + clocksPerTPCTick * (int64_t(hit_end[i]) - hit_tover[i]);       // NOLINT(build/unsigned)
          uint64_t tp_t_end = timestamp + clocksPerTPCTick * int64_t(hit_end[i]);      // NOLINT(build/unsigned)

          // May be needed for TPSet:
          // uint64_t tspan = clocksPerTPCTick * hit_tover[i]; // is/will be this needed?
          //

          // For quick n' dirty debugging: print out time/channel of hits.
          // Can then make a text file suitable for numpy plotting with, eg:
          //
          // sed -n -e 's/.*Hit: \(.*\) \(.*\).*/\1 \2/p' log.txt  > hits.txt
          //
          // TLOG() << "Hit: " << hit_start << " " << offline_channel;

          triggeralgs::TriggerPrimitive trigprim;
          trigprim.time_start = tp_t_begin;
          trigprim.time_peak = (tp_t_begin + tp_t_end) / 2;
          trigprim.time_over_threshold = hit_tover[i] * clocksPerTPCTick;
          trigprim.channel = offline_channel;
          trigprim.adc_integral = hit_charge[i];
          trigprim.adc_peak = hit_charge[i] / 20;
          trigprim.detid =
            m_fiber_no; // TODO: convert crate/slot/fiber to SourceID Roland Sipos rsipos@cern.ch July-22-2021
          trigprim.type = triggeralgs::TriggerPrimitive::Type::kTPC;
          trigprim.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kSimpleThreshold;
          trigprim.version = 1;

          if (m_first_coll) {
            TLOG() << "TP makes sense? -> hit_t_begin:" << tp_t_begin << " hit_t_end:" << tp_t_end
                   << " time_peak:" << (tp_t_begin + tp_t_end) / 2;
          }

          if (!m_tphandler->add_tp(trigprim, timestamp)) {
            m_tps_suppressed_too_long++;
          }

          m_new_tps++;
          ++nhits;
        }
      }
    }

    return nhits;
  }

private:
  bool m_sw_tpg_enabled;
  std::atomic<bool> m_ind_thread_should_run;

  InductionItemToProcess* m_induction_item_to_process;
  std::atomic<bool> m_induction_item_ready{false};

  size_t m_num_msg = 0;
  size_t m_num_push_fail = 0;

  size_t m_num_hits_coll = 0;
  size_t m_num_hits_ind = 0;

  std::atomic<int> m_coll_hits_count{ 0 };
  std::atomic<int> m_indu_hits_count{ 0 };
  std::atomic<int> m_num_tps_pushed{ 0 };

  bool m_first_coll = true;
  bool m_first_ind = true;

  InductionItemToProcess m_dummy_induction_item;

  uint8_t m_fiber_no; // NOLINT(build/unsigned)
  uint8_t m_slot_no;  // NOLINT(build/unsigned)
  uint8_t m_crate_no; // NOLINT(build/unsigned)

  std::shared_ptr<detchannelmaps::TPCChannelMap> m_channel_map;
  swtpg::RegisterChannelMap m_register_channel_map; // Map from
                                                    // expanded AVX
                                                    // register
                                                    // position to
                                                    // offline channel
                                                    // number

  uint32_t m_offline_channel_base;           // NOLINT(build/unsigned)
  uint32_t m_offline_channel_base_induction; // NOLINT(build/unsigned)

  // Frame error check
  bool m_current_frame_pushed = false;
  int m_error_counter_threshold;
  const int m_num_frame_error_bits = 16;
  int m_error_occurrence_counters[16] = { 0 };
  int m_error_reset_freq;
  uint64_t m_clock_frequency; // NOLINT(build/unsigned)

  // Collection
  const uint16_t m_coll_threshold = 5;                    // units of sigma // NOLINT(build/unsigned)
  const uint8_t m_coll_tap_exponent = 6;                  // NOLINT(build/unsigned)
  const int m_coll_multiplier = 1 << m_coll_tap_exponent; // 64
  std::vector<int16_t> m_coll_taps;                       // firwin_int(7, 0.1, multiplier);
  uint16_t* m_coll_primfind_dest;                         // NOLINT(build/unsigned)
  int16_t* m_coll_taps_p;
  std::unique_ptr<swtpg::ProcessingInfo<swtpg::COLLECTION_REGISTERS_PER_FRAME>> m_coll_tpg_pi;

  // Induction
  const uint16_t m_ind_threshold = 5;                   // units of sigma // NOLINT(build/unsigned)
  const uint8_t m_ind_tap_exponent = 6;                 // NOLINT(build/unsigned)
  const int m_ind_multiplier = 1 << m_ind_tap_exponent; // 64
  std::vector<int16_t> m_ind_taps;                      // firwin_int(7, 0.1, multiplier);
  uint16_t* m_ind_primfind_dest;                        // NOLINT(build/unsigned)
  int16_t* m_ind_taps_p;
  std::unique_ptr<swtpg::ProcessingInfo<swtpg::INDUCTION_REGISTERS_PER_FRAME>> m_ind_tpg_pi;
  std::thread m_induction_thread;

  std::shared_ptr<iomanager::SenderConcept<types::TriggerPrimitiveTypeAdapter>> m_tp_sink;
  std::shared_ptr<iomanager::SenderConcept<trigger::TPSet>> m_tpset_sink;
  std::shared_ptr<iomanager::SenderConcept<fddetdataformats::WIBFrame>> m_err_frame_sink;

  std::unique_ptr<WIBTPHandler> m_tphandler;

  std::atomic<uint64_t> m_frame_error_count{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_frames_processed{ 0 };  // NOLINT(build/unsigned)
  daqdataformats::SourceID m_sourceid;

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_suppressed_too_long{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_WIB_WIBFRAMEPROCESSOR_HPP_
