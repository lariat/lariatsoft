/**
 * @file   geometry_unit_test_lariat.h
 * @brief  Class for objects initializing LArIAT geometry
 * @date   May 20th, 2015
 * @author petrillo@fnal.gov
 * 
 * Provides an environment for easy set up of LArIAT-aware tests.
 * Keep in mind that the channel mapping algorithm must be hard-coded and, if
 * using Boost unit test, the configuration file location must be hard coded too
 * (or you can use the provided configuration).
 * 
 * For an example of usage, see larcore/test/Geometry/geometry_iterator_test.cxx
 */

#ifndef TEST_GEO_UNIT_TEST_LARIAT_H
#define TEST_GEO_UNIT_TEST_LARIAT_H

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"

// C/C++ standard libraries
#include <string>


namespace geo {
  class ChannelMapLArIATAlg ;
} // namespace geo

namespace lariat {
  
  /// Namespace including LArIAT-specific testing
  namespace testing {
    
    /** ************************************************************************
     * @brief Class holding the configuration for a LArIAT fixture
     * @tparam CHANNELMAP the class used for channel mapping
     * @see BasicGeometryEnvironmentConfiguration, GeometryTesterFixture
     *
     * This class needs to be fully constructed by the default constructor
     * in order to be useful as Boost unit test fixture.
     * It is supposed to be passed as a template parameter to another class
     * that can store an instance of it and extract configuration information
     * from it.
     * 
     * This class should be used with ChannelMapLArIATAlg .
     * 
     * We reuse BasicGeometryEnvironmentConfiguration as base class
     * and then we fix its setup.
     */
    template <typename CHANNELMAP = geo::ChannelMapLArIATAlg >
    struct LArIATGeometryEnvironmentConfiguration:
      public ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>
    {
      // remember that BasicGeometryEnvironmentConfiguration is not polymorphic
      using base_t
        = ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>;
      
      /// Default constructor; this is what is used in Boost unit test
      LArIATGeometryEnvironmentConfiguration() { LArIATdefaultInit(); }
      
      /// Constructor; accepts the name as parameter
      LArIATGeometryEnvironmentConfiguration(std::string name):
        LArIATGeometryEnvironmentConfiguration()
        { base_t::SetApplicationName(name); }
      
        private:
      
      void LArIATdefaultInit()
        {
          // overwrite the configuration that happened in the base class:
          base_t::SetApplicationName("LArIATGeometryTest");
          base_t::SetDefaultGeometryConfiguration(R"(
              services: {
                Geometry: {
                  SurfaceY: 130.0e2 #in cm, vertical distance to the surface
                  Name:     "lariat"
                  GDML:     "lariat.gdml"
                  ROOT:     "lariat.gdml"
                  SortingParameters: {}
                } # Geometry
              } # services
            )");
        }
    }; // class LArIATGeometryEnvironmentConfiguration<>
    
    
  } // namespace testing
} // namespace lariat

#endif // TEST_GEO_UNIT_TEST_LARIAT_H
