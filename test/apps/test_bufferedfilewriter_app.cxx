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
//#include "readoutlibs/ReadoutTypes.hpp"
#include "fdreadoutlibs/DUNEWIBEthTypeAdapter.hpp"

#include <atomic>
#include <chrono>
#include <memory>
#include <string>

//using namespace dunedaq::readoutlibs;
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
  double buffer_size = std::stod(argv[2]);
  //types::DUMMY_FRAME_STRUCT chunk;
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
  double limiter_freq = 0;
  if (argc == 5) limiter_freq = std::stod(argv[4]);
  auto limiter = dunedaq::readoutlibs::RateLimiter(limiter_freq);

  if (argc == 5) {
    TLOG() << "Starting with ratelimiter at " << limiter_freq << "kHz";
    limiter.init();
  }
  

  while (true) {
    if (!writer.write(reinterpret_cast<char*>(&chunk), sizeof(chunk))) {
      TLOG() << "Could not write to file" << std::endl;
      exit(1);
    }
    bytes_written_total += sizeof(chunk);
    bytes_written_since_last_statistics += sizeof(chunk);
    if (argc == 4) limiter.limit();
  }
}
