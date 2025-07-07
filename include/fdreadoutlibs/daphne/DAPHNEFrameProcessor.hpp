/**
 * @file DAPHNEFrameProcessor.hpp DAPHNE specific Task based raw processor
 *
 * This is part of the DUNE DAQ , copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_

#include "logging/Logging.hpp"

#include "datahandlinglibs/FrameErrorRegistry.hpp"
#include "datahandlinglibs/DataHandlingIssues.hpp"
#include "datahandlinglibs/ReadoutLogging.hpp"

#include "iomanager/IOManager.hpp"
#include "iomanager/Sender.hpp"

#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"
#include "trigger/TriggerPrimitiveTypeAdapter.hpp"
#include "fdreadoutlibs/FDReadoutIssues.hpp"

#include "fddetdataformats/DAPHNEFrame.hpp"
#include "trgdataformats/TriggerPrimitive.hpp"
#include "fdreadoutlibs/DAPHNESuperChunkTypeAdapter.hpp"

#include "appmodel/TPCRawDataProcessor.hpp"
#include "appmodel/PDSRawDataProcessor.hpp"


#include "detchannelmaps/PDSChannelMap.hpp"


#include "appmodel/DataHandlerModule.hpp"
#include "confmodel/Connection.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

using dunedaq::datahandlinglibs::logging::TLVL_BOOKKEEPING;

namespace dunedaq {

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  PDSPeakIgnored,
                  "Ignoring PDS Peak with ts=" << timestamp << ", ch=" << channel << ", sc_iframe=" << superchunk_iframe << ", ipeak=" << ipeak,
                  ((uint64_t)timestamp) ((uint64_t)channel) ((size_t)superchunk_iframe) ((size_t)ipeak))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  PDSUnphysicalFrameTimestamp,
                  "PDS Frame with unphysical timestamp detected with ts=" << timestamp << ", ch=" << channel << ", sc_iframe=" << superchunk_iframe,
                  ((uint64_t)timestamp) ((uint64_t)channel) ((size_t)superchunk_iframe))
                  
  
namespace fdreadoutlibs {

class DAPHNEFrameProcessor : public datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>
{

public:
  using inherited = datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>;
  using frameptr = types::DAPHNESuperChunkTypeAdapter*;
  using daphneframeptr = dunedaq::fddetdataformats::DAPHNEFrame*;
  using timestamp_t = std::uint64_t; // NOLINT(build/unsigned)
  using constframeptr = const types::DAPHNESuperChunkTypeAdapter*;

  // Constructor
  explicit DAPHNEFrameProcessor(std::unique_ptr<datahandlinglibs::FrameErrorRegistry>& error_registry, bool post_processing_enabled)
    : datahandlinglibs::TaskRawDataProcessorModel<types::DAPHNESuperChunkTypeAdapter>(error_registry, post_processing_enabled)
  {}

 

  // Override config for pipeline setup
  void conf(const appmodel::DataHandlerModule* conf) override;

  void start(const nlohmann::json& args) override;
  void stop(const nlohmann::json& args) override;

protected:
  virtual void generate_opmon_data() override;

  /**
   * Pipeline Stage 1.: Check proper timestamp increments in DAPHNE frame
   * */
  void timestamp_check(frameptr /*fp*/);

  /**
   * Pipeline Stage 2.: Check DAPHNE headers for error flags
   * */
  void frame_error_check(frameptr /*fp*/);

  // Internals
  timestamp_t m_previous_ts = 0;
  timestamp_t m_current_ts = 0;
  bool m_first_ts_fake = true;
  bool m_first_ts_missmatch = true;
  bool m_problem_reported = false;
  std::atomic<int> m_ts_error_ctr{ 0 };

  void extract_tps( constframeptr fp);
  dunedaq::trgdataformats::TriggerPrimitive peak_to_tp( dunedaq::fddetdataformats::DAPHNEFrame &frame, int i);
  
private:

  //PDSChannelMap
  std::shared_ptr<detchannelmaps::PDSChannelMap> m_channel_map;
  std::vector<std::pair<trgdataformats::channel_t, int16_t>> m_channel_plane_numbers;

  uint32_t m_det_id; // NOLINT(build/unsigned)
  uint32_t m_crate_id; // NOLINT(build/unsigned)
  uint32_t m_slot_id;  // NOLINT(build/unsigned)
  uint32_t m_stream_id; // NOLINT(build/unsigned)

  std::set<unsigned int> m_channel_mask_set;
  uint32_t m_def_adc_intg_thresh = 0;



  std::shared_ptr<iomanager::SenderConcept<std::vector<trigger::TriggerPrimitiveTypeAdapter>>> m_tp_sink;

  std::atomic<uint64_t> m_new_hits{ 0 }; // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_new_tps{ 0 };  // NOLINT(build/unsigned)
  std::atomic<uint64_t> m_tps_suppressed_too_long{ 0 };
  std::atomic<uint64_t> m_tps_send_failed{ 0 };
  std::atomic<uint64_t> m_frame_counter{ 0 };

  std::chrono::time_point<std::chrono::high_resolution_clock> m_t0;
  

};

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_DAPHNE_DAPHNEFRAMEPROCESSOR_HPP_
