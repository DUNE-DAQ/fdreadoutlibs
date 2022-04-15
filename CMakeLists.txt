cmake_minimum_required(VERSION 3.12)
project(fdreadoutlibs VERSION 1.1.0)

find_package(daq-cmake REQUIRED)

daq_setup_environment()

find_package(ers REQUIRED)
find_package(appfwk REQUIRED)
find_package(logging REQUIRED)
find_package(opmonlib REQUIRED)
find_package(readoutlibs REQUIRED)
find_package(daqdataformats REQUIRED)
find_package(detchannelmaps REQUIRED)
find_package(detdataformats REQUIRED)
find_package(trigger REQUIRED)
find_package(triggeralgs REQUIRED)
find_package(folly REQUIRED)
find_package(Boost COMPONENTS iostreams REQUIRED)
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
  wib/tpg/*.cpp
  LINK_LIBRARIES ${FDREADOUTLIBS_DEPENDENCIES}
)

##############################################################################
# Plugins
#daq_add_plugin(FakeCardReader duneDAQModule LINK_LIBRARIES readoutmodules)
#daq_add_plugin(DataLinkHandler duneDAQModule LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_plugin(DataRecorder duneDAQModule LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_plugin(FragmentConsumer duneDAQModule LINK_LIBRARIES readoutmodules)
#daq_add_plugin(TimeSyncConsumer duneDAQModule LINK_LIBRARIES readoutmodules)

##############################################################################
# Integration tests
#daq_add_application(readoutmodules_test_ratelimiter test_ratelimiter_app.cxx TEST LINK_LIBRARIES readoutmodules)
#daq_add_application(readoutmodules_test_lb_allocation test_lb_allocation_app.cxx TEST LINK_LIBRARIES readoutmodules)
#daq_add_application(readoutmodules_test_bufferedfilewriter test_bufferedfilewriter_app.cxx TEST LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_application(readoutmodules_test_bufferedfilereader test_bufferedfilereader_app.cxx TEST LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_application(readoutmodules_test_skiplist test_skiplist_app.cxx TEST LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_application(readoutmodules_test_fast_expand_wib2frame test_fast_expand_wib2frame_app.cxx TEST LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#
#
###############################################################################
## Unit Tests
#daq_add_unit_test(RawWIBTp_test                LINK_LIBRARIES readoutmodules)
#daq_add_unit_test(BufferedReadWrite_test       LINK_LIBRARIES readoutmodules ${BOOST_LIBS})
#daq_add_unit_test(VariableSizeElementQueue_test LINK_LIBRARIES readoutmodules)

##############################################################################
# Installation
daq_install()
