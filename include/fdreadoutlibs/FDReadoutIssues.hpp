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

#include <ers/Issue.hpp>
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
                  "The selected algorithm does not exist: " << algorithm_selection << " . Check your configuration file and seelect either SWTPG or AbsRS.",
                  ((std::string)algorithm_selection))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPTooLong,
		  "TP with ToT " << width << " for channel " << channel,
		  ((uint64_t)width) ((uint64_t)channel))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPDropped,
                  "TP with ToT " << width << " for channel " << channel,
                  ((uint64_t)width) ((uint64_t)channel))


ERS_DECLARE_ISSUE(fdreadoutlibs,
                  LinkMisconfiguration,
                  "WIB data have crate/slot/link " << wcrate << "/" << wslot << "/" << wlink << " while this readout link is configured for " << crate << "/" << slot << "/" << link,
                  ((uint32_t)wcrate) ((uint32_t)wslot) ((uint32_t)wlink) ((uint32_t)crate) ((uint32_t)slot) ((uint32_t)link))


} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_FDREADOUTISSUES_HPP_
