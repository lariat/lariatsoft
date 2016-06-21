/**
 * @file   geometry_iterator_lariat_test.cxx
 * @brief  Unit test for geometry iterators on LArIAT detector
 * @date   May 20th, 2015
 * @author petrillo@fnal.gov
 * 
 * Usage: just run the executable.
 * Boost unit testing environment keeps the arguments secret anyway.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by geometry_boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryIteratorTestLArIAT

// LArSoft libraries
#include "test/Geo/geometry_unit_test_lariat.h"
#include "larcore/TestUtils/boost_unit_test_base.h"
#include "test/Geometry/GeometryIteratorTestAlg.h"
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
// used for output only; we use LArIATGeometryFixtureConfigurer as base
// class, that is already configured to use LArIAT geometry.
struct LArIATGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    lariat::testing::LArIATGeometryEnvironmentConfiguration
      <geo::ChannelMapLArIATAlg>
    >
{
  /// Constructor: overrides the application name
  LArIATGeometryConfiguration()
    { SetApplicationName("GeometryIteratorUnitTest"); }
}; // class LArIATGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterFixture, configured with the object
 * above.
 * It provides:
 * - `Tester`, a configured instance of the test algorithm.
 */
class LArIATGeometryIteratorTestFixture:
  private testing::GeometryTesterEnvironment<LArIATGeometryConfiguration>
{
    public:
  geo::GeometryIteratorTestAlg Tester;
  
  /// Constructor: initialize the tester with the Geometry from base class
  LArIATGeometryIteratorTestFixture(): Tester(TesterParameters())
    { Tester.Setup(*Geometry()); }
  
}; // class LArIATGeometryIteratorTestFixture



//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_FIXTURE_TEST_SUITE
  (GeometryIteratorsLArIAT, LArIATGeometryIteratorTestFixture)
// BOOST_GLOBAL_FIXTURE(LArIATGeometryIteratorTestFixture)


BOOST_AUTO_TEST_CASE( AllTests )
{
  Tester.Run();
} // BOOST_AUTO_TEST_CASE( AllTests )

/*
BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )
{
  Tester.CryostatIteratorsTest();
} // BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )



BOOST_AUTO_TEST_CASE( TPCIteratorsTest )
{
  Tester.TPCIteratorsTest();
} // BOOST_AUTO_TEST_CASE( TPCIteratorsTest )



BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )
{
  Tester.PlaneIteratorsTest();
} // BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )



BOOST_AUTO_TEST_CASE( WireIteratorsTest )
{
  Tester.WireIteratorsTest();
} // BOOST_AUTO_TEST_CASE( WireIteratorsTest )
*/

BOOST_AUTO_TEST_SUITE_END()

