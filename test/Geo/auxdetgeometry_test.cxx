/**
 * @file   geometry_test.cxx
 * @brief  Unit test for geometry functionalities on a standard detector
 * @date   May 5th, 2015
 * @author petrillo@fnal.gov
 * 
 * Usage:
 *   geometry_test  [ConfigurationFile [GeometryTestParameterSet]]
 * 
 * By default, GeometryTestParameterSet is set to "physics.analysers.geotest".
 * 
 * This unit test uses geometry_unit_test_base.h to build an environment with a
 * geometry set up.
 * For an example of use with Boost unit test module, see
 * geometry_iterator_test.cxx .
 */


// LArSoft libraries
#include "test/Geometry/GeometryTestAlg.h"
#include "test/Geometry/geometry_unit_test_base.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// we use an existing class provided for this purpose, since our test
// environment allows us to tailor it at run time.
using StandardGeometryConfiguration
  = testing::BasicGeometryEnvironmentConfiguration<geo::ChannelMapStandardAlg>;

/*
 * GeometryTesterFixture, configured with the object above, is used in a
 * non-Boost-unit-test context.
 * It provides:
 * - `geo::GeometryCore const* Geometry()`
 * - `geo::GeometryCore const* GlobalGeometry()` (static member)
 */
using StandardGeometryTestEnvironment
  = testing::GeometryTesterEnvironment<StandardGeometryConfiguration>;


//------------------------------------------------------------------------------
//---  The tests
//---

/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 * 
 * The arguments in argv are:
 * 0. name of the executable ("Geometry_test")
 * 1. path to the FHiCL configuration file
 * 2. FHiCL path to the configuration of the geometry test
 *    (default: physics.analysers.geotest)
 * 3. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 * 
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv) {
  
  StandardGeometryConfiguration config("geometry_test");
  
  //
  // parameter parsing
  //
  int iParam = 0;
  
  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);
  
  // second argument: path of the parameter set for geometry test configuration
  // (optional; default: "physics.analysers.geotest")
  config.SetMainTesterParameterSetPath
    ((++iParam < argc)? argv[iParam]: "physics.analyzers.geotest");
  
  // third argument: path of the parameter set for geometry configuration
  // (optional; default: "services.Geometry" from the inherited object)
  if (++iParam < argc) config.SetGeometryParameterSetPath(argv[iParam]);
  
  //
  // testing environment setup
  //
  StandardGeometryTestEnvironment TestEnvironment(config);
  
  //
  // run the test algorithm
  //
  
  // 1. we initialize it from the configuration in the environment,
  geo::GeometryTestAlg Tester(TestEnvironment.TesterParameters());
  
  // 2. we set it up with the geometry from the environment
  Tester.Setup(*TestEnvironment.Geometry());
  
  // 3. then we run it!
  unsigned int nErrors = Tester.Run();
  
  // 4. And finally we cross fingers.
  if (nErrors > 0) {
    mf::LogError("geometry_test") << nErrors << " errors detected!";
  }
  
  return nErrors;
} // main()
