/**
 * @file tpg_apps_issues.hpp
 * 
 * ERS issues for TPG applications
 *
 * This is part of the DUNE DAQ , copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef FDREADOUTLIBS_APPS_TPG_ISSUES_HPP_
#define FDREADOUTLIBS_APPS_TPG_ISSUES_HPP_

#include <ers/Issue.hpp>
#include <string>




namespace tpgtools {

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  TPGAlgorithmInexistent,
                  "The selected algorithm does not exist: " << algorithm_selection << " . Check your command line options and select either SimpleThreshold or AbsRS.",
                  ((std::string)algorithm_selection))

ERS_DECLARE_ISSUE(fdreadoutlibs,
                  FileInexistent,
                  "The selected input file does not exist. Input file: " << input_file_path << "  Check the path of the input file.",
                  ((std::string)input_file_path))



} // namespace dunedaq

#endif // FDREADOUTLIBS_APPS_TPG_ISSUES_HPP_