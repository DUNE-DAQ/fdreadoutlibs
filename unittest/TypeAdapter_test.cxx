
#include "fdreadoutlibs/RefactoredDAPHNEEthTypeAdapter.hpp"

#define BOOST_TEST_MODULE TypeAdapters_test // NOLINT

#include "boost/test/unit_test.hpp"

using namespace dunedaq::fdreadoutlibs::types;

BOOST_AUTO_TEST_SUITE(TypeAdapters_test)

BOOST_AUTO_TEST_CASE(Sizes)
{
  RefactoredDAPHNEEthTypeAdapter daphne_eth_type_adapter {};
  uint64_t a = 0;
  uint64_t b = 1;
  daphne_eth_type_adapter.fake_timestamps(a, b);
}

BOOST_AUTO_TEST_SUITE_END()
