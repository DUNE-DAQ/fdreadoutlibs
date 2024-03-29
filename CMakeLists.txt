cmake_minimum_required(VERSION 3.12)
project(fdreadoutlibs VERSION 1.9.0)

find_package(daq-cmake REQUIRED)

daq_setup_environment()

find_package(ers REQUIRED)
find_package(iomanager REQUIRED)
find_package(appfwk REQUIRED)
find_package(logging REQUIRED)
find_package(opmonlib REQUIRED)
find_package(readoutlibs REQUIRED)
find_package(daqdataformats REQUIRED)
find_package(detchannelmaps REQUIRED)
find_package(detdataformats REQUIRED)
find_package(fddetdataformats REQUIRED)
find_package(trigger REQUIRED)
find_package(triggeralgs REQUIRED)
find_package(folly REQUIRED)
find_package(Boost COMPONENTS iostreams unit_test_framework REQUIRED)
set(BOOST_LIBS Boost::iostreams ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_LIBRARIES})

#daq_codegen( readoutconfig.jsonnet datalinkhandler.jsonnet  datarecorder.jsonnet  sourceemulatorconfig.jsonnet TEMPLATES Structs.hpp.j2 Nljs.hpp.j2 )
#daq_codegen( *info.jsonnet DEP_PKGS opmonlib TEMPLATES opmonlib/InfoStructs.hpp.j2 opmonlib/InfoNljs.hpp.j2 )

##############################################################################
# Dependency sets
set(FDREADOUTLIBS_DEPENDENCIES
#tools
  Folly::folly
  ers::ers
  logging::logging
#dunedaq
  appfwk::appfwk
  readoutlibs::readoutlibs
  opmonlib::opmonlib
  daqdataformats::daqdataformats
  detdataformats::detdataformats
  fddetdataformats::fddetdataformats
  trigger::trigger
  triggeralgs::triggeralgs
  detchannelmaps::detchannelmaps
)

##############################################################################
# Extra options and tweaks
set(FDREADOUTLIBS_USE_INTRINSICS ON)

if(${FDREADOUTLIBS_USE_INTRINSICS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()

##############################################################################
# Main library
daq_add_library(
  *.cpp
  #  wib/*.cpp
  #wib/tpg/*.cpp
  wib2/*.cpp
  wib2/tpg/*.cpp
  wibeth/*.cpp
  wibeth/tpg/*.cpp
  daphne/*.cpp
  tde/*.cpp
  LINK_LIBRARIES ${FDREADOUTLIBS_DEPENDENCIES}
)

##############################################################################
# Applications



##############################################################################
# Test Apps

daq_add_application(tde_file_creator tde_file_creator.cxx TEST LINK_LIBRARIES fdreadoutlibs)
daq_add_application(wib2_test_bench wib2_test_bench.cxx TEST LINK_LIBRARIES fdreadoutlibs)
daq_add_application(wibeth_test_bench wibeth_test_bench.cxx TEST LINK_LIBRARIES fdreadoutlibs)

##############################################################################

daq_add_unit_test(DAPHNEStreamSuperChunkTypeAdapter_test LINK_LIBRARIES fdreadoutlibs)
daq_add_unit_test(WIBEthFrameExpansion_test LINK_LIBRARIES fdreadoutlibs)

##############################################################################
# Installation
daq_install()
