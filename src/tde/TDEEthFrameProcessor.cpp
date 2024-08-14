/**
 * @file TDEEthFrameProcessor.hpp WIBEth specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "fdreadoutlibs/tde/TDEEthFrameProcessor.hpp" // NOLINT(build/include)
#include "confmodel/GeoId.hpp"
#include "appmodel/RawDataProcessor.hpp"

#include "iomanager/Sender.hpp"
#include "logging/Logging.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"
#include "datahandlinglibs/models/IterableQueueModel.hpp"
// #include "datahandlinglibs/readoutconfig/Nljs.hpp"
//#include "datahandlinglibs/readoutinfo/InfoNljs.hpp"
#include "datahandlinglibs/utils/ReusableThread.hpp"

#include "fddetdataformats/TDEEthFrame.hpp"


#include "fdreadoutlibs/TDEEthTypeAdapter.hpp"
#include "trigger/TriggerPrimitiveTypeAdapter.hpp"

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

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;
using dunedaq::datahandlinglibs::logging::TLVL_TAKE_NOTE;

// THIS SHOULDN'T BE HERE!!!!! But it is necessary.....
DUNE_DAQ_TYPESTRING(dunedaq::trigger::TriggerPrimitiveTypeAdapter, "TriggerPrimitive")


namespace dunedaq {
namespace fdreadoutlibs {

TDEEthFrameProcessor::TDEEthFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry)
  : TaskRawDataProcessorModel<types::TDEEthTypeAdapter>(error_registry)
{
}

TDEEthFrameProcessor::~TDEEthFrameProcessor()
{
}

void
TDEEthFrameProcessor::start(const nlohmann::json& args)
{
  // Reset timestamp check
  m_previous_ts = 0;
  m_current_ts = 0;
  m_first_ts_missmatch = true;
  m_ts_problem_reported = false;
  m_ts_error_ctr = 0;

  m_first_seq_id_mismatch = true;
  m_seq_id_problem_reported = false;
  m_seq_id_error_ctr = 0;

  inherited::start(args);
}

void
TDEEthFrameProcessor::stop(const nlohmann::json& args)
{
  inherited::stop(args);
}

void
TDEEthFrameProcessor::conf(const appmodel::DataHandlerModule* conf)
{
  // auto config = cfg["rawdataprocessorconf"].get<datahandlinglibs::readoutconfig::RawDataProcessorConf>();

  for (auto output : conf->get_outputs()) {
    try {
      if (output->get_data_type() == "TriggerPrimitive") {
         m_tp_sink = get_iom_sender<trigger::TriggerPrimitiveTypeAdapter>(output->UID());
      }
    } catch (const ers::Issue& excpt) {
      ers::error(datahandlinglibs::ResourceQueueError(ERS_HERE, "tp", "DefaultRequestHandlerModel", excpt));
    }
  }

  m_sourceid.id = conf->get_source_id();
  m_sourceid.subsystem = types::TDEEthTypeAdapter::subsystem;
  auto geo_id = conf->get_geo_id();
  if (geo_id != nullptr) {
    m_det_id = geo_id->get_detector_id();
    m_crate_id = geo_id->get_crate_id();
    m_slot_id = geo_id->get_slot_id();
    m_stream_id = geo_id->get_stream_id();
  }
  m_emulator_mode = conf->get_emulation_mode();

  inherited::conf(conf);
}

void
TDEEthFrameProcessor::generate_opmon_data()
{
   datahandlinglibs::opmon::FixedRateDataProcessorInfo info;

   info.set_num_seq_id_errors(m_seq_id_error_ctr.load());
   info.set_min_seq_id_jump(m_seq_id_min_jump.exchange(0));
   info.set_max_seq_id_jump(m_seq_id_max_jump.exchange(0));

   info.set_num_ts_errors(m_ts_error_ctr.load());

   publish(std::move(info));

   inherited::generate_opmon_data();
}

/**
 * Pipeline Stage 1.: Check proper timestamp increments in WIB frame
 * */
void
TDEEthFrameProcessor::sequence_check(frameptr fp)
{

  // // If EMU data, emulate perfectly incrementing timestamp
  // if (inherited::m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
  //   // uint64_t ts_next = m_previous_seq_id + 1; // NOLINT(build/unsigned)
  //   auto tf = reinterpret_cast<tdeframeptr>(((uint8_t*)fp));            // NOLINT
  //   for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
  //     //auto wfh = const_cast<tdeframeptr>(tf->header());
  //     tf->daq_header.crate_id = m_crate_no;
  //     tf->daq_header.slot_id = m_slot_no;
  //     tf->daq_header.stream_id = m_stream_id; 
  //     tf->daq_header.seq_id = (m_previous_seq_id+i) & 0xfff;
  //     tf++;
  //   }
  // }

  // Acquire timestamp
  auto wfptr = reinterpret_cast<dunedaq::fddetdataformats::TDEEthFrame*>(fp); // NOLINT
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

    m_error_registry->add_error("SEQUENCE_ID_JUMP", datahandlinglibs::FrameErrorRegistry::ErrorInterval(expected_seq_id, m_current_seq_id));
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
TDEEthFrameProcessor::timestamp_check(frameptr fp)
{

  uint16_t tdeeth_tick_difference = types::TDEEthTypeAdapter::expected_tick_difference;
  uint16_t tdeeth_frame_tick_difference = tdeeth_tick_difference * fp->get_num_frames();

  // If EMU data, emulate perfectly incrementing timestamp
  // if (inherited::m_emulator_mode) {                                     // emulate perfectly incrementing timestamp
  //   uint64_t ts_next = m_previous_ts + tdeeth_frame_tick_difference; // NOLINT(build/unsigned)
  //   auto tf = reinterpret_cast<tdeframeptr>(((uint8_t*)fp));            // NOLINT
  //   for (unsigned int i = 0; i < fp->get_num_frames(); ++i) {           // NOLINT(build/unsigned)
  //     //auto wfh = const_cast<tdeframeptr>(tf->header());
  //     tf->daq_header.crate_id = m_crate_no;
  //     tf->daq_header.slot_id = m_slot_no;
  //     tf->daq_header.stream_id = m_stream_id; 
  //     tf->set_timestamp(ts_next);
  //     ts_next += tdeeth_tick_difference;
  //     tf++;
  //   }
  // }

  // Acquire timestamp
  auto wfptr = reinterpret_cast<tdeframeptr>(fp); // NOLINT
  m_current_ts = wfptr->get_timestamp();

  // Check timestamp
  if (m_current_ts - m_previous_ts != tdeeth_frame_tick_difference) {
    ++m_ts_error_ctr;
    m_error_registry->add_error("MISSING_FRAMES", datahandlinglibs::FrameErrorRegistry::ErrorInterval(m_previous_ts + tdeeth_frame_tick_difference, m_current_ts));
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

} // namespace fdreadoutlibs
} // namespace dunedaq
