/**
 * This class is a hack to match configured TPG processing steps to the legacy, defined
 * algorithms.
 */

#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPGALGORITHMCLASSIFIER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPGALGORITHMCLASSIFIER_HPP_

#include "trgdataformats/TriggerPrimitive.hpp"

#include <unordered_map>
#include <vector>

namespace dunedaq {
namespace fdreadoutlibs {

class TPGAlgorithmClassifier {
  std::vector<int> m_step_order{};
  const std::unordered_map<std::string, int> m_tp_alg_id_map = {
    {"AVXFrugalPedestalSubtractProcessor", 1},
    {"AVXAbsRunSumProcessor", 2},
    {"AVXRunSumProcessor", 3},
    {"AVXThresholdProcessor", 4}
  };

  public:
    /**
     * @brief Take the processing step name and append the int value to m_step_order.
     */
    void
    append_processing_step(const std::string& step_name) {
      auto search = m_tp_alg_id_map.find(step_name);
      if (search != m_tp_alg_id_map.end())
        m_step_order.push_back(search->second);
    }

    /**
     * @brief Check against the defined algorithms if these steps match.
     */
    trgdataformats::TriggerPrimitive::Algorithm
    get_tpg_algorithm() {
      const std::vector<int> abs_rs_order = {1, 2, 1, 4}; // Order that defines kAbsRunningSum.
      const std::vector<int> rs_order = {1, 3, 1, 4}; // Order that defines kRunningSum.
      const std::vector<int> st_order = {1, 4}; // Order that defines kSimpleThreshold.

      if (m_step_order == abs_rs_order)
        return trgdataformats::TriggerPrimitive::Algorithm::kAbsRunningSum;
      if (m_step_order == rs_order)
        return trgdataformats::TriggerPrimitive::Algorithm::kRunningSum;
      if (m_step_order == st_order)
        return trgdataformats::TriggerPrimitive::Algorithm::kSimpleThreshold;

      // If it is none of the above, then it is unknown.
      return trgdataformats::TriggerPrimitive::Algorithm::kUnknown;
    }
  };

} // namespace fdreadoutlibs
} // namespace dunedaq

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TPGALGORITHMCLASSIFIER_HPP_
