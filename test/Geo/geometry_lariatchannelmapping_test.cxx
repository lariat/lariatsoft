/**
 * @file   geometry_lariatchannelmapping_test.cxx
 * @brief  Unit test of channel mapping on LArIAT detector
 * @date   October 12th, 2016
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * Usage: just run the executable.
 * Or plug a FHiCL file in the command line.
 * 
 * This is pretty much the copy of the test for the "standard" channel mapping.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryLArIATChannelMappingTest

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
#include "test/Geometry/ChannelMapStandardTestAlg.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "Geo/ChannelMapLArIATAlg.h"

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; BasicGeometryEnvironmentConfiguration can read the
// configuration file name from command line, and
// BoostCommandLineConfiguration<> makes it initialize in time for Boost
// to catch it when instanciating the fixture.
struct LArIATGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    testing::BasicGeometryEnvironmentConfiguration<geo::ChannelMapLArIATAlg>
    >
{
  /// Constructor: overrides the application name
  LArIATGeometryConfiguration()
    { SetApplicationName("GeometryLArIATChannelMappingTest"); }
}; // class LArIATGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterEnvironment, configured with the object
 * above.
 * It provides to the testing environment:
 * - `Tester()`, returning a configured instance of the test algorithm;
 * - `GlobalTester()`, (static) returning a global configured instance of the
 *   test algorithm.
 * 
 * The testing::TestSharedGlobalResource<> facility provides a singleton
 * instance of the tester algorithm, shared through the whole program and with
 * this object too.
 * 
 * This sharing allows the fixture to be used either as global or as per-suite.
 * In the former case, the BOOST_AUTO_TEST_CASE's will access the global test
 * algorithm instance through the static call to
 * `GeometryLArIATChannelMappingTestFixture::GlobalTester()`; in the latter
 * case, it will access the local tester via the member function `Tester()`.
 * In this case, whether `GlobalTester()` and `Tester()` point to the same
 * tester depends on Boost unit test implementation.
 */
class GeometryLArIATChannelMappingTestFixture:
  private testing::GeometryTesterEnvironment<LArIATGeometryConfiguration>
{
  using Tester_t = geo::ChannelMapStandardTestAlg;
  
  using TesterRegistry_t = testing::TestSharedGlobalResource<Tester_t>;

    public:
  
  /// Constructor: initialize the tester with the Geometry from base class
  GeometryLArIATChannelMappingTestFixture()
    {
      // create a new tester
      tester_ptr = std::make_shared<Tester_t>(TesterParameters());
      tester_ptr->Setup(*(Provider<geo::GeometryCore>()));
      // if no tester is default yet, share ours:
      TesterRegistry_t::ProvideDefaultSharedResource(tester_ptr);
    }
  
  /// Retrieves the local tester
  Tester_t& Tester() { return *(tester_ptr.get()); }
  
  /// Retrieves the global tester
  static Tester_t& GlobalTester() { return TesterRegistry_t::Resource(); }
  
    private:
  std::shared_ptr<Tester_t> tester_ptr; ///< our tester (may be shared)
}; // class GeometryLArIATChannelMappingTestFixture



//------------------------------------------------------------------------------
//---  The tests
//---
//
// Note on Boost fixture options:
// - BOOST_FIXTURE_TEST_SUITE will provide the tester as a always-present data
//   member in the environment, as "Tester()"; but a new fixture, with a new
//   geometry and a new tester, is initialized on each test case
// - BOOST_GLOBAL_FIXTURE does not provide tester access, so one has to get it
//   as GeometryLArIATChannelMappingTestFixture::GlobalTester(); on the other
//   hand, the fixture is initialised only when a new global one is explicitly
//   created.
//

BOOST_GLOBAL_FIXTURE(GeometryLArIATChannelMappingTestFixture);

BOOST_AUTO_TEST_CASE( TPCsetMappingTestCase )
{
  GeometryLArIATChannelMappingTestFixture::GlobalTester().TPCsetMappingTest();
} // BOOST_AUTO_TEST_CASE( TPCsetMappingTestCase )

BOOST_AUTO_TEST_CASE( ROPMappingTestCase )
{
  GeometryLArIATChannelMappingTestFixture::GlobalTester().ROPMappingTest();
} // BOOST_AUTO_TEST_CASE( ROPMappingTestCase )

BOOST_AUTO_TEST_CASE( ChannelMappingTestCase )
{
  GeometryLArIATChannelMappingTestFixture::GlobalTester()
    .ChannelMappingTest();
} // BOOST_AUTO_TEST_CASE( ChannelMappingTestCase )


// BOOST_AUTO_TEST_SUITE_END()
