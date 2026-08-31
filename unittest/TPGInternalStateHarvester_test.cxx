/**
 * @file TPGInternalStateHarvester_test.cxx TPGInternalStateHarvester class Unit Tests
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2025.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifdef TPGLIBS_ENABLE_STATE_MONITORING

#include "fdreadoutlibs/tpg/TPGInternalStateHarvester.hpp"
#include "tpglibs/AbstractProcessor.hpp"
#include "tpglibs/ProcessorMetricArray.hpp"

#define BOOST_TEST_MODULE TPGInternalStateHarvester_test // NOLINT

#include "boost/test/unit_test.hpp"
#include <thread>
#include <chrono>
#include <random>

using namespace dunedaq::fdreadoutlibs;
using namespace dunedaq::trgdataformats;

// Mock processor for testing
class MockProcessor : public tpglibs::AbstractProcessor<__m256i> {
public:
  MockProcessor(const std::vector<std::string>& metric_names, 
                const std::vector<std::array<int16_t, 16>>& metric_values)
    : m_metric_names(metric_names), m_metric_values(metric_values) {}

  void configure(const nlohmann::json& config, const int16_t* plane_numbers) override {
    // Mock implementation - do nothing
  }

  std::vector<std::string> get_requested_internal_state_names() const override {
    return m_metric_names;
  }

  tpglibs::ProcessorMetricArray<std::array<int16_t, 16>> read_internal_states_as_integer_array() override {
    tpglibs::ProcessorMetricArray<std::array<int16_t, 16>> result;
    result.m_size = m_metric_values.size();
    result.m_data = m_metric_values.data();
    return result;
  }

private:
  std::vector<std::string> m_metric_names;
  std::vector<std::array<int16_t, 16>> m_metric_values;
};

BOOST_AUTO_TEST_SUITE(TPGInternalStateHarvester_test)

BOOST_AUTO_TEST_CASE(ConstructorDestructor)
{
  TPGInternalStateHarvester harvester;
  // Test that constructor and destructor work without issues
  BOOST_REQUIRE(true);
}

BOOST_AUTO_TEST_CASE(SetGetProcessorReferences)
{
  TPGInternalStateHarvester harvester;
  
  // Create mock processors
  std::vector<std::string> metric_names = {"baseline", "accumulator"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}},
    {{200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215}}
  };
  
  auto proc1 = std::make_shared<MockProcessor>(metric_names, metric_values);
  auto proc2 = std::make_shared<MockProcessor>(metric_names, metric_values);
  
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {
    {proc1, 0}, {proc2, 1}
  };
  
  harvester.set_processor_references(refs);
  
  const auto& retrieved_refs = harvester.get_processor_references();
  BOOST_REQUIRE_EQUAL(retrieved_refs.size(), 2);
  BOOST_REQUIRE_EQUAL(retrieved_refs[0].second, 0);
  BOOST_REQUIRE_EQUAL(retrieved_refs[1].second, 1);
}

BOOST_AUTO_TEST_CASE(UpdateChannelPlaneNumbers)
{
  TPGInternalStateHarvester harvester;
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers = {
    {100, 0}, {101, 0}, {102, 1}, {103, 1}, {104, 2}, {105, 2}
  };
  
  harvester.update_channel_plane_numbers(channel_plane_numbers, 2, 3);
  
  // Test that the method completes without throwing
  BOOST_REQUIRE(true);
}

BOOST_AUTO_TEST_CASE(HarvestOnceBasic)
{
  TPGInternalStateHarvester harvester;
  
  // Setup mock processors
  std::vector<std::string> metric_names = {"baseline", "accumulator"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}},
    {{200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215}}
  };
  
  auto proc1 = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc1, 0}};
  harvester.set_processor_references(refs);
  
  // Setup channel-plane mapping
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers = {
    {100, 0}, {101, 0}, {102, 0}, {103, 0}, {104, 0}, {105, 0}, {106, 0}, {107, 0},
    {108, 0}, {109, 0}, {110, 0}, {111, 0}, {112, 0}, {113, 0}, {114, 0}, {115, 0}
  };
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Perform harvest
  auto results = harvester.harvest_once();
  
  // Verify results
  BOOST_REQUIRE_EQUAL(results.size(), 16); // 16 channels
  
  // Check first channel
  auto it = results.find(100);
  BOOST_REQUIRE(it != results.end());
  BOOST_REQUIRE_EQUAL(it->second.size(), 2); // 2 metrics
  
  // Check metric values
  BOOST_REQUIRE_EQUAL(it->second[0].first, "baseline");
  BOOST_REQUIRE_EQUAL(it->second[0].second, 100);
  BOOST_REQUIRE_EQUAL(it->second[1].first, "accumulator");
  BOOST_REQUIRE_EQUAL(it->second[1].second, 200);
}

BOOST_AUTO_TEST_CASE(HarvestOnceMultiplePipelines)
{
  TPGInternalStateHarvester harvester;
  
  // Setup mock processors for 2 pipelines
  std::vector<std::string> metric_names = {"baseline"};
  std::vector<std::array<int16_t, 16>> metric_values_pipeline0 = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}}
  };
  std::vector<std::array<int16_t, 16>> metric_values_pipeline1 = {
    {{200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215}}
  };
  
  auto proc0 = std::make_shared<MockProcessor>(metric_names, metric_values_pipeline0);
  auto proc1 = std::make_shared<MockProcessor>(metric_names, metric_values_pipeline1);
  
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {
    {proc0, 0}, {proc1, 1}
  };
  harvester.set_processor_references(refs);
  
  // Setup channel-plane mapping for 2 pipelines of 16 channels each
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 32; ++i) {
    channel_plane_numbers.push_back({100 + i, i / 16}); // First 16 in plane 0, next 16 in plane 1
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 2);
  
  // Perform harvest
  auto results = harvester.harvest_once();
  
  // Verify results
  BOOST_REQUIRE_EQUAL(results.size(), 32); // 32 channels total
  
  // Check first channel from pipeline 0
  auto it0 = results.find(100);
  BOOST_REQUIRE(it0 != results.end());
  BOOST_REQUIRE_EQUAL(it0->second[0].second, 100);
  
  // Check first channel from pipeline 1
  auto it1 = results.find(116);
  BOOST_REQUIRE(it1 != results.end());
  BOOST_REQUIRE_EQUAL(it1->second[0].second, 200);
}

BOOST_AUTO_TEST_CASE(ThreadManagement)
{
  TPGInternalStateHarvester harvester;
  
  // Test thread state before starting
  BOOST_REQUIRE(!harvester.is_collection_thread_running());
  
  // Start collection thread
  harvester.start_collection_thread();
  BOOST_REQUIRE(harvester.is_collection_thread_running());
  
  // Give the thread a moment to fully start
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  // Stop collection thread
  harvester.stop_collection_thread();
  BOOST_REQUIRE(!harvester.is_collection_thread_running());
}

BOOST_AUTO_TEST_CASE(TriggerHarvest)
{
  TPGInternalStateHarvester harvester;
  
  // Setup mock processor
  std::vector<std::string> metric_names = {"baseline"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({100 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Start collection thread
  harvester.start_collection_thread();
  
  // Give thread time to start
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  // Trigger harvest
  harvester.trigger_harvest();
  
  // Wait a bit for processing
  std::this_thread::sleep_for(std::chrono::milliseconds(50));
  
  // Get results
  auto results = harvester.get_latest_results();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Stop collection thread
  harvester.stop_collection_thread();
}

BOOST_AUTO_TEST_CASE(GetLatestResultsEmpty)
{
  TPGInternalStateHarvester harvester;
  
  // Get results without starting thread or triggering harvest
  auto results = harvester.get_latest_results();
  BOOST_REQUIRE(results.empty());
}

BOOST_AUTO_TEST_CASE(ConcurrentAccess)
{
  TPGInternalStateHarvester harvester;
  
  // Setup mock processor
  std::vector<std::string> metric_names = {"baseline"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({100 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Start collection thread
  harvester.start_collection_thread();
  
  // Give thread time to start
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  // Spawn multiple threads to trigger harvests and read results
  std::vector<std::thread> threads;
  std::atomic<int> success_count{0};
  
  for (int i = 0; i < 3; ++i) {
    threads.emplace_back([&harvester, &success_count]() {
      for (int j = 0; j < 5; ++j) {
        harvester.trigger_harvest();
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        auto results = harvester.get_latest_results();
        if (!results.empty()) {
          success_count++;
        }
      }
    });
  }
  
  // Wait for all threads to complete
  for (auto& t : threads) {
    t.join();
  }
  
  // Stop collection thread
  harvester.stop_collection_thread();
  
  // Verify that at least some operations succeeded
  BOOST_REQUIRE(success_count.load() > 0);
}

BOOST_AUTO_TEST_CASE(EmptyProcessorReferences)
{
  TPGInternalStateHarvester harvester;
  
  // Set empty processor references
  std::vector<TPGInternalStateHarvester::ProcRef> empty_refs;
  harvester.set_processor_references(empty_refs);
  
  // Setup channel-plane mapping
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({100 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Perform harvest - should handle empty processors gracefully
  auto results = harvester.harvest_once();
  BOOST_REQUIRE(results.empty());
}

BOOST_AUTO_TEST_CASE(MismatchedMetricSizes)
{
  TPGInternalStateHarvester harvester;
  
  // Create processor with mismatched metric names and values (3 names, 2 values)
  std::vector<std::string> metric_names = {"baseline", "accumulator", "threshold"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}},
    {{200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215}}
    // Missing third metric value - this causes size mismatch
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({100 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Perform harvest - should skip processor when sizes don't match
  auto results = harvester.harvest_once();
  // When metric names and values don't match, the entire processor is skipped
  BOOST_REQUIRE_EQUAL(results.size(), 0);
}

BOOST_AUTO_TEST_CASE(StressTest)
{
  TPGInternalStateHarvester harvester;
  
  // Create multiple processors with random data
  std::vector<std::string> metric_names = {"baseline", "accumulator", "threshold"};
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 1000);
  
  std::vector<TPGInternalStateHarvester::ProcRef> refs;
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  
  for (int pipeline = 0; pipeline < 4; ++pipeline) {
    std::vector<std::array<int16_t, 16>> metric_values;
    for (const auto& name : metric_names) {
      std::array<int16_t, 16> values;
      for (int i = 0; i < 16; ++i) {
        values[i] = dis(gen);
      }
      metric_values.push_back(values);
    }
    
    auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
    refs.push_back({proc, pipeline});
    
    for (int i = 0; i < 16; ++i) {
      channel_plane_numbers.push_back({100 + pipeline * 16 + i, pipeline % 3});
    }
  }
  
  harvester.set_processor_references(refs);
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 4);
  
  // Perform multiple harvests
  for (int i = 0; i < 100; ++i) {
    auto results = harvester.harvest_once();
    BOOST_REQUIRE_EQUAL(results.size(), 64); // 4 pipelines * 16 channels
    
    // Verify some random channels
    for (int j = 0; j < 10; ++j) {
      channel_t test_channel = 100 + (j % 64);
      auto it = results.find(test_channel);
      BOOST_REQUIRE(it != results.end());
      BOOST_REQUIRE_EQUAL(it->second.size(), 3); // 3 metrics
    }
  }
}

BOOST_AUTO_TEST_CASE(ExactValueCollectionVerification)
{
  TPGInternalStateHarvester harvester;
  
  // Create processor with known, specific values
  std::vector<std::string> metric_names = {"baseline", "accumulator", "threshold"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{-100, -99, -98, -97, -96, -95, -94, -93, -92, -91, -90, -89, -88, -87, -86, -85}}, // baseline
    {{1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015}}, // accumulator  
    {{50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65}} // threshold
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  // Setup channel mapping: channels 100-115 map to lanes 0-15
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({100 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Harvest and verify exact values
  auto results = harvester.harvest_once();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Verify each channel has exactly 3 metrics
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 100 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 3);
    
    // Verify metric names and exact values
    const auto& metrics = it->second;
    
    // Find each metric by name and verify its value
    for (const auto& [name, value] : metrics) {
      if (name == "baseline") {
        BOOST_REQUIRE_EQUAL(value, -100 + i); // lane i corresponds to channel 100+i
      } else if (name == "accumulator") {
        BOOST_REQUIRE_EQUAL(value, 1000 + i);
      } else if (name == "threshold") {
        BOOST_REQUIRE_EQUAL(value, 50 + i);
      } else {
        BOOST_FAIL("Unexpected metric name: " << name);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(BoundaryValueCollection)
{
  TPGInternalStateHarvester harvester;
  
  // Test with boundary values (min/max int16_t)
  std::vector<std::string> metric_names = {"extreme_values"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{INT16_MIN, INT16_MIN + 1, -1, 0, 1, INT16_MAX - 1, INT16_MAX, 
      INT16_MIN, INT16_MIN + 1, -1, 0, 1, INT16_MAX - 1, INT16_MAX, INT16_MIN, INT16_MAX}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({200 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  auto results = harvester.harvest_once();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Verify boundary values are preserved exactly
  std::array<int16_t, 16> expected_values = {
    INT16_MIN, INT16_MIN + 1, -1, 0, 1, INT16_MAX - 1, INT16_MAX,
    INT16_MIN, INT16_MIN + 1, -1, 0, 1, INT16_MAX - 1, INT16_MAX, INT16_MIN, INT16_MAX
  };
  
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 200 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 1);
    BOOST_REQUIRE_EQUAL(it->second[0].first, "extreme_values");
    BOOST_REQUIRE_EQUAL(it->second[0].second, expected_values[i]);
  }
}

BOOST_AUTO_TEST_CASE(MultiPipelineValueMapping)
{
  TPGInternalStateHarvester harvester;
  
  // Create two processors with different values
  std::vector<std::string> metric_names = {"pipeline_id"};
  
  // Pipeline 0: values 100-115
  std::vector<std::array<int16_t, 16>> metric_values_pipeline0 = {
    {{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115}}
  };
  
  // Pipeline 1: values 200-215  
  std::vector<std::array<int16_t, 16>> metric_values_pipeline1 = {
    {{200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215}}
  };
  
  auto proc0 = std::make_shared<MockProcessor>(metric_names, metric_values_pipeline0);
  auto proc1 = std::make_shared<MockProcessor>(metric_names, metric_values_pipeline1);
  
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {
    {proc0, 0}, {proc1, 1}
  };
  harvester.set_processor_references(refs);
  
  // Setup channels: 300-315 for pipeline 0, 400-415 for pipeline 1
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({300 + i, 0}); // pipeline 0
  }
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({400 + i, 1}); // pipeline 1
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 2);
  
  auto results = harvester.harvest_once();
  BOOST_REQUIRE_EQUAL(results.size(), 32);
  
  // Verify pipeline 0 values (channels 300-315)
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 300 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 1);
    BOOST_REQUIRE_EQUAL(it->second[0].first, "pipeline_id");
    BOOST_REQUIRE_EQUAL(it->second[0].second, 100 + i);
  }
  
  // Verify pipeline 1 values (channels 400-415)
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 400 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 1);
    BOOST_REQUIRE_EQUAL(it->second[0].first, "pipeline_id");
    BOOST_REQUIRE_EQUAL(it->second[0].second, 200 + i);
  }
}

BOOST_AUTO_TEST_CASE(ValueConsistencyAcrossHarvests)
{
  TPGInternalStateHarvester harvester;
  
  // Create processor with fixed values
  std::vector<std::string> metric_names = {"stable_metric"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({500 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Perform multiple harvests and verify values remain consistent
  for (int harvest = 0; harvest < 10; ++harvest) {
    auto results = harvester.harvest_once();
    BOOST_REQUIRE_EQUAL(results.size(), 16);
    
    // Verify all values are exactly as expected
    for (int i = 0; i < 16; ++i) {
      channel_t channel = 500 + i;
      auto it = results.find(channel);
      BOOST_REQUIRE(it != results.end());
      BOOST_REQUIRE_EQUAL(it->second.size(), 1);
      BOOST_REQUIRE_EQUAL(it->second[0].first, "stable_metric");
      BOOST_REQUIRE_EQUAL(it->second[0].second, 42 + i);
    }
  }
}

BOOST_AUTO_TEST_CASE(MetricNameValuePairOrdering)
{
  TPGInternalStateHarvester harvester;
  
  // Create processor with multiple metrics in specific order
  std::vector<std::string> metric_names = {"first", "second", "third", "fourth"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},   // first
    {{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},   // second
    {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}},   // third
    {{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}}    // fourth
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({600 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  auto results = harvester.harvest_once();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Verify metric ordering is preserved
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 600 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 4);
    
    // Verify order and values
    BOOST_REQUIRE_EQUAL(it->second[0].first, "first");
    BOOST_REQUIRE_EQUAL(it->second[0].second, 1);
    BOOST_REQUIRE_EQUAL(it->second[1].first, "second");
    BOOST_REQUIRE_EQUAL(it->second[1].second, 2);
    BOOST_REQUIRE_EQUAL(it->second[2].first, "third");
    BOOST_REQUIRE_EQUAL(it->second[2].second, 3);
    BOOST_REQUIRE_EQUAL(it->second[3].first, "fourth");
    BOOST_REQUIRE_EQUAL(it->second[3].second, 4);
  }
}

BOOST_AUTO_TEST_CASE(AsyncValueCollectionVerification)
{
  TPGInternalStateHarvester harvester;
  
  // Setup processor with known values
  std::vector<std::string> metric_names = {"async_metric"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({700 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  // Start collection thread
  harvester.start_collection_thread();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  
  // Trigger harvest and wait for results
  harvester.trigger_harvest();
  std::this_thread::sleep_for(std::chrono::milliseconds(50));
  
  auto results = harvester.get_latest_results();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Verify exact values were collected asynchronously
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 700 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 1);
    BOOST_REQUIRE_EQUAL(it->second[0].first, "async_metric");
    BOOST_REQUIRE_EQUAL(it->second[0].second, 777 + i);
  }
  
  harvester.stop_collection_thread();
}

BOOST_AUTO_TEST_CASE(ZeroValueCollection)
{
  TPGInternalStateHarvester harvester;
  
  // Test with all zero values
  std::vector<std::string> metric_names = {"zero_metric"};
  std::vector<std::array<int16_t, 16>> metric_values = {
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}
  };
  
  auto proc = std::make_shared<MockProcessor>(metric_names, metric_values);
  std::vector<TPGInternalStateHarvester::ProcRef> refs = {{proc, 0}};
  harvester.set_processor_references(refs);
  
  std::vector<std::pair<channel_t, int16_t>> channel_plane_numbers;
  for (int i = 0; i < 16; ++i) {
    channel_plane_numbers.push_back({800 + i, 0});
  }
  harvester.update_channel_plane_numbers(channel_plane_numbers, 16, 1);
  
  auto results = harvester.harvest_once();
  BOOST_REQUIRE_EQUAL(results.size(), 16);
  
  // Verify all values are exactly zero
  for (int i = 0; i < 16; ++i) {
    channel_t channel = 800 + i;
    auto it = results.find(channel);
    BOOST_REQUIRE(it != results.end());
    BOOST_REQUIRE_EQUAL(it->second.size(), 1);
    BOOST_REQUIRE_EQUAL(it->second[0].first, "zero_metric");
    BOOST_REQUIRE_EQUAL(it->second[0].second, 0);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif // TPGLIBS_ENABLE_STATE_MONITORING
