#ifndef FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TRIGGERACTIVITYTYPEADAPTER_HPP_
#define FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TRIGGERACTIVITYTYPEADAPTER_HPP_


#include "daqdataformats/FragmentHeader.hpp"
#include "daqdataformats/SourceID.hpp"
#include "trgdataformats/TriggerActivity.hpp"
#include "trgdataformats/Types.hpp"

#include <cstdint> // uint_t types
#include <memory>  // unique_ptr
#include <vector>
#include <cstring> // memcpy
#include <tuple> // tie

namespace dunedaq {
namespace fdreadoutlibs {
namespace types {

//const constexpr std::size_t kTriggerActivitySize = sizeof(trgdataformats::TriggerActivity);  

struct TriggerActivityTypeAdapter
{
    triggeralgs::TriggerActivity activity;
    std::vector<uint8_t> activity_overlay_buffer;
    
    // Don't really want this default ctor, but IterableQueueModel requires it
    TriggerActivityTypeAdapter() {}
    
    TriggerActivityTypeAdapter(triggeralgs::TriggerActivity a)
      : activity(a)
    {
      populate_buffer();
    }

    void populate_buffer()
    {
      activity_overlay_buffer.resize(triggeralgs::get_overlay_nbytes(activity));
      triggeralgs::write_overlay(activity, activity_overlay_buffer.data());
    }
    
    // comparable based on first timestamp
    bool operator<(const TriggerActivityTypeAdapter& other) const
    {
      return this->activity.time_start < other.activity.time_start;
    }

    uint64_t get_first_timestamp() const // NOLINT(build/unsigned)
    {
      return activity.time_start;
    }

    void set_first_timestamp(uint64_t ts) // NOLINT(build/unsigned)
    {
      activity.time_start = ts;
    }

    uint64_t get_timestamp() const // NOLINT(build/unsigned)
    {
      return activity.time_start;
    }

    size_t get_payload_size() { return activity_overlay_buffer.size(); }

    size_t get_num_frames() { return 1; }

    size_t get_frame_size() { return get_payload_size(); }

    uint8_t* begin()
    {
      return activity_overlay_buffer.data();
    }
    
    uint8_t* end()
    {
      return activity_overlay_buffer.data()+activity_overlay_buffer.size();
    }

    //static const constexpr size_t fixed_payload_size = 5568;
    static const constexpr daqdataformats::SourceID::Subsystem subsystem = daqdataformats::SourceID::Subsystem::kTrigger;
    static const constexpr daqdataformats::FragmentType fragment_type = daqdataformats::FragmentType::kTriggerActivity;
    // No idea what this should really be set to
    static const constexpr uint64_t expected_tick_difference = 16; // NOLINT(build/unsigned)

#endif // FDREADOUTLIBS_INCLUDE_FDREADOUTLIBS_TRIGGERACTIVITYTYPEADAPTER_HPP_
