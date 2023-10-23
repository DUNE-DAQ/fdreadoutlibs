/**
 * @file test_bufferedfilewriter_app.cxx Test application for
 * BufferedFileWriter implementation
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */
#include "readoutlibs/utils/BufferedFileWriter.hpp"
#include "readoutlibs/utils/RateLimiter.hpp"

#include "logging/Logging.hpp"
#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"

#include <atomic>
#include <chrono>
#include <memory>
#include <string>

using namespace dunedaq::fdreadoutlibs;

int
main(int argc, char* argv[])
{
  if (argc < 3 || argc > 5 || argc==4 || (argc == 5 && strcmp(argv[3], "-L") != 0)) {
    TLOG() << "usage: readoutlibs_test_bufferedfilewriter filename buffer_size <-L rate_limiter_frequency>" << std::endl;
    TLOG() << "-L frequency parameter is optional. Limiter will be disable" << std::endl;    
    exit(1);
  }
  remove(argv[1]); // NOLINT
  std::string filename(argv[1]);
  int buffer_size = std::stoi(argv[2]);
  TLOG() << "buffer_size: " << buffer_size << " Btyes";
  double limiter_freq = 0;
  bool use_limiter = false;
  if (argc == 5 && std::strcmp(argv[3], " -L")) {
    use_limiter = true;
    limiter_freq = std::stod(argv[4]);
    TLOG() << "ratelimiter freq: " << limiter_freq << "kHz";
  }

  types::DUNEWIBEthTypeAdapter chunk;

  dunedaq::readoutlibs::BufferedFileWriter writer(filename, buffer_size);
  for (uint i = 0; i < sizeof(chunk); ++i) {
    (reinterpret_cast<char*>(&chunk))[i] = static_cast<char>(i); // NOLINT
  }

  std::atomic<int64_t> bytes_written_total = 0;
  std::atomic<int64_t> bytes_written_since_last_statistics = 0;
  std::chrono::steady_clock::time_point time_point_last_statistics = std::chrono::steady_clock::now();

  auto statistics_thread = std::thread([&]() {
    while (true) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      double time_diff = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() -
                                                                                   time_point_last_statistics)
                           .count();
      TLOG() << "Bytes written: " << bytes_written_total << ", Throughput: "
             << static_cast<double>(bytes_written_since_last_statistics) / ((int64_t)1 << 20) / time_diff << " MiB/s"
             << std::endl;
      time_point_last_statistics = std::chrono::steady_clock::now();
      bytes_written_since_last_statistics = 0;
    }
  });

  // Initializing limiter
  if (use_limiter){
    TLOG() << "Starting with ratelimiter at " << limiter_freq << "kHz";
    auto limiter = dunedaq::readoutlibs::RateLimiter(limiter_freq);
    limiter.init();

    while (true) {
      if (!writer.write(reinterpret_cast<char*>(&chunk), sizeof(chunk))) {
        TLOG() << "Could not write to file" << std::endl;
        exit(1);
      }
      bytes_written_total += sizeof(chunk);
      bytes_written_since_last_statistics += sizeof(chunk);
      limiter.limit();
    }
  }
  else {
    while (true) {
      if (!writer.write(reinterpret_cast<char*>(&chunk), sizeof(chunk))) {
        TLOG() << "Could not write to file" << std::endl;
        exit(1);
      }
      bytes_written_total += sizeof(chunk);
      bytes_written_since_last_statistics += sizeof(chunk);
    }
  }

}
