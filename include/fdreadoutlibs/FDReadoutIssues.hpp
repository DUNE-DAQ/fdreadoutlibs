/**
 * @file FDReadoutIssues.hpp Readout system related 
 * ERS issues for fdreadoutlibs
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTISSUES_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTISSUES_HPP_

#include "daqdataformats/Types.hpp"

#include <ers/Issue.hpp>
#include "logging/Logging.hpp" // NOTE: if ISSUES ARE DECLARED BEFORE include logging/Logging.hpp, TLOG_DEBUG<<issue wont work.
#include <string>




namespace dunedaq {
/*
ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPHandlerBacklog,
                  "Failed to push hits to TP handler " << sid,
                  ((int)sid))
*/

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGAlgorithmInexistent,
                  "The selected algorithm does not exist: " << algorithm_selection << " . Check your configuration file and seelect either SimpleThreshold or AbsRS.",
                  ((std::string)algorithm_selection))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  FrameAndTPCountersDisabled,
                  "Both frame_count_limit and tp_count_limit were set to 0 (disabled) in the TPCRawDataProcessor config. TPs will not send.",
                  ) 

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPTooLong,
                  "TP with SOT " << width << " for channel " << channel,
                  ((uint64_t)width) ((uint64_t)channel))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  FailedToSendTPVector,
                  "Failed to send TP vector beginning with start time " << s_ts_begin << " and channel number " << channel_begin << ", ending with start time " << s_ts_end << " and channel number " << channel_end,
                  ((daqdataformats::timestamp_t)s_ts_begin) ((uint64_t)channel_begin) ((daqdataformats::timestamp_t)s_ts_end) ((uint64_t)channel_end))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  DetectorPlaneToTPSinkMismatch,
                  "There are more detector planes " << num_planes << " than available TP sinks " << num_tp_sinks << ".",
                  ((size_t) num_planes) ((size_t) num_tp_sinks))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  LinkMisconfiguration,
                  "WIB data have crate/slot/link " << wcrate << "/" << wslot << "/" << wlink << " while this readout link is configured for " << crate << "/" << slot << "/" << link,
                  ((uint32_t)wcrate) ((uint32_t)wslot) ((uint32_t)wlink) ((uint32_t)crate) ((uint32_t)slot) ((uint32_t)link))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  PDSPeakIgnored,
                  "Ignoring PDS Peak with ts=" << timestamp << ", ch=" << channel << ", sc_iframe=" << superchunk_iframe << ", ipeak=" << ipeak,
                  ((uint64_t)timestamp) ((uint64_t)channel) ((size_t)superchunk_iframe) ((size_t)ipeak))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  PDSUnphysicalFrameTimestamp,
                  "PDS Frame with unphysical timestamp detected with ts=" << timestamp << ", ch=" << channel << ", sc_iframe=" << superchunk_iframe,
                  ((uint64_t)timestamp) ((uint64_t)channel) ((size_t)superchunk_iframe))
ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGStateMonitoringDisabledAtBuildTime,
                  "TPG state monitoring config flags are set but tpglibs was built with "
                  "TPGLIBS_ENABLE_STATE_MONITORING=OFF. No processor metrics will be "
                  "collected. To enable: export TPGLIBS_ENABLE_STATE_MONITORING=ON, then "
                  "rebuild from a clean build directory (dbt-build -c).",
                  )

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGToggleStateDeprecated,
                  "ProcessingStep attribute 'metric_collect_toggle_state' is set but this "
                  "flag is deprecated. State monitoring is now controlled at build time: "
                  "export TPGLIBS_ENABLE_STATE_MONITORING=ON, then rebuild from a clean "
                  "build directory (dbt-build -c). This attribute will be removed in a "
                  "future release.",
                  )

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGStateMonitoringConfigIgnored,
                  "ProcessingStep attribute '" << attribute_name << "' is configured but "
                  "has no effect because TPGLIBS_ENABLE_STATE_MONITORING=OFF. To enable: "
                  "export TPGLIBS_ENABLE_STATE_MONITORING=ON, then rebuild from a clean "
                  "build directory (dbt-build -c).",
                  ((std::string)attribute_name))

} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTISSUES_HPP_
