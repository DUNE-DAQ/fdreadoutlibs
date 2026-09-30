/** Bounded offline replay of concatenated DAQ DAPHNE 512-byte frames. */
#include "fdreadoutlibs/pds/ActivityBuilder.hpp"
#include "CLI/CLI.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <tuple>

using namespace dunedaq::fdreadoutlibs::pds;

namespace {
using Address = std::array<uint32_t, 5>;

void print_light(const std::optional<PromptLight>& light, bool incomplete)
{
  if (!light) return;
  std::cout << "{\"kind\":\"descriptor_light_proxy\",\"start\":" << light->time_start
            << ",\"end_exclusive\":" << light->time_end
            << ",\"prompt_integral\":" << light->prompt_integral
            << ",\"total_integral\":" << light->total_integral << ",\"fraction\":";
  if (light->fraction()) std::cout << *light->fraction();
  else std::cout << "null";
  std::cout << ",\"pulse_count\":" << light->pulse_count
            << ",\"incomplete_input\":" << (incomplete ? "true" : "false") << "}\n";
}
} // namespace

int main(int argc, char** argv)
{
  CLI::App app{"PDS descriptor TP/TA replay; finite-file sorting, not a live link merger"};
  std::string input_path, map_path;
  bool demo = false;
  bool allow_overflow = false;
  size_t overflow_frames = 0;
  double window_us = 10, prompt_us = 0.1, total_us = 10;
  size_t minimum_channels = 6, max_tps = 1000000;
  uint32_t minimum_integral = 0;
  app.add_option("--input", input_path, "Concatenated DAQ-format 512-byte frames (not PCAP/HDF5/504-byte payloads)");
  app.add_option("--channel-map", map_path, "CSV: det,crate,slot,stream,hardware_channel,offline_channel");
  app.add_flag("--demo", demo, "Six synthetic channels; no hardware access");
  app.add_flag("--allow-descriptor-overflow", allow_overflow, "Inspect retained peaks from incomplete frames; all output is marked incomplete");
  app.add_option("--window-us", window_us, "Inclusive coincidence window");
  app.add_option("--min-channels", minimum_channels, "Distinct offline channels (default 6)");
  app.add_option("--prompt-us", prompt_us, "Prompt gate width");
  app.add_option("--total-us", total_us, "Total light gate width");
  app.add_option("--minimum-integral", minimum_integral, "Per-descriptor ADC integral cut");
  app.add_option("--max-tps", max_tps, "Bound replay memory; fail instead of truncating");
  CLI11_PARSE(app, argc, argv);
  try {
    if (!max_tps || (demo && (!input_path.empty() || !map_path.empty())) ||
        (!demo && (input_path.empty() || map_path.empty()))) {
      throw std::invalid_argument("Use either --demo or both --input and --channel-map; max-tps must be positive");
    }
    CoincidenceBuilder coincidence({microseconds_to_ticks(window_us), minimum_channels, max_tps, 1});
    const auto total_ticks = microseconds_to_ticks(total_us);
    PromptLightBuilder prompt({microseconds_to_ticks(prompt_us), total_ticks, 1});
    std::map<Address, uint32_t> channels;
    if (!demo) {
      std::ifstream mapping(map_path);
      if (!mapping) throw std::runtime_error("Cannot open channel map");
      std::string line;
      while (std::getline(mapping, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        Address address{};
        uint32_t offline;
        std::string extra;
        if (!(fields >> address[0] >> address[1] >> address[2] >> address[3] >> address[4] >> offline) ||
            fields >> extra || address[0] >= 64 || address[1] >= 1024 || address[2] >= 16 ||
            address[3] >= 256 || address[4] >= 256 || offline >= dunedaq::trgdataformats::INVALID_TP_CHANNEL ||
            !channels.emplace(address, offline).second) {
          throw std::invalid_argument("Malformed or duplicate channel-map row");
        }
      }
      // Aliases would silently merge physically distinct channels in coincidence counting.
      std::set<uint32_t> offline_channels;
      for (const auto& entry : channels) {
        if (!offline_channels.insert(entry.second).second) throw std::invalid_argument("Duplicate offline channel mapping");
      }
    }
    const ChannelMap map = [&channels, demo](const Frame& frame) -> uint32_t {
      if (demo) return frame.get_channel(); // Explicitly synthetic single-board numbering.
      Address address{uint32_t(frame.daq_header.det_id), uint32_t(frame.daq_header.crate_id),
                      uint32_t(frame.daq_header.slot_id), uint32_t(frame.daq_header.stream_id), frame.get_channel()};
      return channels.at(address);
    };
    std::vector<TP> tps;
    auto append = [&](const Frame& frame) {
      auto decoded = descriptor_tps(frame, map, {minimum_integral, 1, !allow_overflow});
      if (frame.header.descriptor_overflow) ++overflow_frames;
      if (decoded.size() > max_tps - tps.size()) throw std::length_error("Replay max-tps exceeded");
      tps.insert(tps.end(), decoded.begin(), decoded.end());
    };
    if (demo) {
      for (unsigned channel = 0; channel < 6; ++channel) {
        Frame frame{};
        frame.header.version = Frame::version;
        frame.header.fragment_descriptor = 1;
        frame.daq_header.det_id = 2;
        frame.set_timestamp(100000);
        frame.set_channel(channel);
        auto& peak = frame.header.peaks_data.peaks[0];
        peak.found = 1;
        peak.sample_start = channel * 2;
        peak.time_peak = 1;
        peak.duration_minus_one = 3;
        peak.adc_peak = 100;
        peak.adc_integral = 250;
        append(frame);
      }
    } else {
      const uint16_t endian = 1;
      if (*reinterpret_cast<const uint8_t*>(&endian) != 1) throw std::runtime_error("Replay requires a little-endian host");
      std::ifstream input(input_path, std::ios::binary);
      if (!input) throw std::runtime_error("Cannot open frame file");
      Frame frame{};
      while (input.read(reinterpret_cast<char*>(&frame), sizeof(frame))) append(frame);
      if (input.gcount() || !input.eof()) throw std::runtime_error("Truncated frame or input read error");
    }
    std::sort(tps.begin(), tps.end(), [](const TP& a, const TP& b) {
      return std::make_tuple(a.time_start, a.channel) < std::make_tuple(b.time_start, b.channel);
    });
    // Reject exact repeated records, including replayed network packets. At-most-once contract.
    for (size_t i = 1; i < tps.size(); ++i) {
      if (tps[i].time_start == tps[i-1].time_start && tps[i].channel == tps[i-1].channel) {
        throw std::invalid_argument("Duplicate/overlapping TP start for one channel; deduplicate source frames first");
      }
    }
    size_t activities = 0;
    for (const auto& tp : tps) {
      std::cout << "{\"kind\":\"tp\",\"channel\":" << tp.channel << ",\"start\":" << tp.time_start
                << ",\"samples\":" << tp.samples_over_threshold << ",\"samples_to_peak\":" << tp.samples_to_peak
                << ",\"integral\":" << tp.adc_integral << ",\"peak\":" << tp.adc_peak
                << ",\"incomplete_input\":" << (overflow_frames ? "true" : "false") << "}\n";
      if (auto ta = coincidence.push(tp)) {
        ++activities;
        std::set<uint32_t> distinct;
        for (const auto& hit : ta->inputs) distinct.insert(hit.channel);
        std::cout << "{\"kind\":\"pds_ta\",\"algorithm\":\"pds_distinct_channels_prototype\",\"start\":"
                  << ta->time_start << ",\"end_exclusive\":" << ta->time_end << ",\"activity_time\":" << ta->time_activity
                  << ",\"distinct_channels\":" << distinct.size() << ",\"integral\":" << ta->adc_integral
                  << ",\"incomplete_input\":" << (overflow_frames ? "true" : "false") << "}\n";
      }
      print_light(prompt.push(tp), overflow_frames != 0);
    }
    // Explicit end-of-file watermark closes the final quiet gate.
    print_light(prompt.finish(), overflow_frames != 0);
    std::cerr << "Decoded " << tps.size() << " TPs; emitted " << activities << " PDS TAs; "
              << overflow_frames << " incomplete descriptor-overflow frames\n";
  } catch (const std::exception& error) {
    std::cerr << "PDS replay failed: " << error.what() << '\n';
    return 1;
  }
}
