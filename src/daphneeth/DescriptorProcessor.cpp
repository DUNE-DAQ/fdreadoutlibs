#include "fdreadoutlibs/pds/DescriptorProcessor.hpp"
#include <map>
#include <mutex>

namespace dunedaq::fdreadoutlibs::pds {
namespace {
std::mutex mutex;
std::map<std::string, std::weak_ptr<DescriptorProcessor>> processors;
}
void register_descriptor_processor(const std::string& key, const std::shared_ptr<DescriptorProcessor>& processor)
{
  std::lock_guard lock(mutex);
  if (processors[key].lock()) throw std::logic_error("Duplicate descriptor processor: " + key);
  processors[key] = processor;
}
std::shared_ptr<DescriptorProcessor> get_descriptor_processor(const std::string& key)
{
  std::lock_guard lock(mutex);
  auto found = processors.find(key);
  return found == processors.end() ? nullptr : found->second.lock();
}
void remove_descriptor_processor(const std::string& key, const std::shared_ptr<DescriptorProcessor>& processor)
{
  std::lock_guard lock(mutex);
  auto found = processors.find(key);
  if (found != processors.end() && found->second.lock() == processor) processors.erase(found);
}
} // namespace dunedaq::fdreadoutlibs::pds
