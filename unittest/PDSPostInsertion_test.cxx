#define BOOST_TEST_MODULE PDSPostInsertion_test
#include "boost/test/unit_test.hpp"
#include "fdreadoutlibs/DAPHNEEthTypeAdapter.hpp"
#include "fdreadoutlibs/pds/DescriptorToTP.hpp"
#include "datahandlinglibs/models/DataHandlingModel.hpp"
#include "datahandlinglibs/models/DefaultSkipListRequestHandler.hpp"
#include "datahandlinglibs/models/SkipListLatencyBufferModel.hpp"
#include "datahandlinglibs/models/TaskRawDataProcessorModel.hpp"
#include "PDSDescriptorFixtures.hpp"
#include <future>
#include <cstring>
#include <thread>
#include <vector>

using Adapter = dunedaq::fdreadoutlibs::types::DAPHNEEthTypeAdapter;
namespace dhl = dunedaq::datahandlinglibs;
using Buffer = dhl::SkipListLatencyBufferModel<Adapter>;

class Processor : public dhl::TaskRawDataProcessorModel<Adapter>
{
public:
  using TaskRawDataProcessorModel::TaskRawDataProcessorModel;
  std::promise<void> stop_entered;
  void stop(const dunedaq::appfwk::DAQModule::CommandData_t& args) override
  {
    stop_entered.set_value();
    TaskRawDataProcessorModel::stop(args);
  }
  void prepare(std::function<void(const Adapter*)> callback)
  {
    add_postprocess_task(std::move(callback));
    m_items_to_postprocess_queues.push_back(
      std::make_unique<folly::ProducerConsumerQueue<const Adapter*>>(16));
  }
};

class Model : public dhl::DataHandlingModel<Adapter, dhl::DefaultSkipListRequestHandler<Adapter>, Buffer, Processor>
{
public:
  using DataHandlingModel::DataHandlingModel;
  void prepare(std::shared_ptr<Buffer> buffer, std::shared_ptr<Processor> processor)
  {
    m_latency_buffer_impl = std::move(buffer);
    m_raw_processor_impl = std::move(processor);
    m_request_handler_supports_cutoff_timestamp = false;
    m_processing_delay_ticks = 0;
  }
  void receive(Adapter frame) { process_item(std::move(frame)); }
  void finish() { stop_processing_and_flush(dunedaq::appfwk::DAQModule::CommandData_t{}); }
};

BOOST_AUTO_TEST_CASE(out_of_order_frames_run_once_on_worker_and_drain_before_flush)
{
  std::atomic<bool> running{false};
  Model model(running);
  auto buffer = std::make_shared<Buffer>();
  auto errors = std::make_unique<dhl::FrameErrorRegistry>();
  auto processor = std::make_shared<Processor>(errors, true);
  std::promise<void> release;
  auto gate = release.get_future().share();
  std::vector<uint64_t> starts;
  bool buffer_alive = true;
  bool on_worker = true;
  const auto producer = std::this_thread::get_id();
  processor->prepare([&](const Adapter* input) {
    gate.wait();
    buffer_alive = buffer_alive && buffer->occupancy() == 3;
    if (!buffer_alive) return;
    on_worker = on_worker && std::this_thread::get_id() != producer;
    dunedaq::fdreadoutlibs::pds::Frame frame{};
    std::memcpy(&frame, input->data, sizeof(frame));
    for (const auto& tp : dunedaq::fdreadoutlibs::pds::descriptor_tps(frame, pds_test::mapping)) {
      starts.push_back(tp.time_start);
    }
  });
  model.prepare(buffer, processor);
  processor->start(dunedaq::appfwk::DAQModule::CommandData_t{});
  for (uint64_t ts : {3000, 1000, 2000}) {
    auto frame = pds_test::frame();
    frame.set_timestamp(ts);
    Adapter input{};
    std::memcpy(input.data, &frame, sizeof(frame));
    model.receive(std::move(input));
  }
  BOOST_CHECK(starts.empty()); // Worker is held; reception did not run conversion inline.
  BOOST_CHECK_EQUAL(buffer->occupancy(), 3);
  auto stopping = processor->stop_entered.get_future();
  std::thread stopper([&] { model.finish(); });
  const bool reached_stop = stopping.wait_for(std::chrono::seconds(5)) == std::future_status::ready;
  const auto retained_at_stop = buffer->occupancy();
  release.set_value();
  stopper.join();
  BOOST_CHECK(reached_stop);
  BOOST_CHECK_EQUAL(retained_at_stop, 3);
  BOOST_REQUIRE_EQUAL(starts.size(), 3);
  const std::vector<uint64_t> expected{3007, 1007, 2007};
  BOOST_CHECK_EQUAL_COLLECTIONS(starts.begin(), starts.end(), expected.begin(), expected.end());
  BOOST_CHECK(buffer_alive);
  BOOST_CHECK(on_worker);
  BOOST_CHECK_EQUAL(buffer->occupancy(), 0);
}
