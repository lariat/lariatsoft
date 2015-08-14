/**
 * @file   AuxDetGeometryTestAlg.h
 * @brief  Unit test for auxiliary detector geometry functionalities
 * @date   2015/06/24
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryTestAlg.cxx
 * 
 */

#ifndef GEO_AUXDETGEOMETRYTESTALG_H
#define GEO_AUXDETGEOMETRYTESTALG_H

// LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h"

// C/C++ standard libraries
#include <string>
#include <set>
#include <vector>
#include <array>
#include <memory> // std::unique_ptr<>

// forward declarations
namespace fhicl {
  class ParameterSet;
}


namespace geo {
  
  // forward declarations
  class AuxDetGeometryCore;
  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  
  namespace details {
    class TestTrackerClassBase;
  }
  
  
  /** **************************************************************************
   * @brief Performs tests on the auxiliary detector geometry as seen by AuxDetGeometry service
   * 
   * Configuration parameters
   * =========================
   * 
   */
  class AuxDetGeometryTestAlg {
      public:
    explicit AuxDetGeometryTestAlg(fhicl::ParameterSet const& pset);
    
    /// Virtual destructor
    virtual ~AuxDetGeometryTestAlg() = default;
    
    /// Runs the test
    virtual void Setup(geo::AuxDetGeometryCore const& new_geo) { geom = &new_geo; }

    /// Runs the test, returns a number of errors (very unlikely!)
    virtual unsigned int Run();


  private:
    geo::AuxDetGeometryCore const* geom; ///< pointer to geometry service provider

    std::set<std::string> fNonFatalExceptions;
    
    std::unique_ptr<details::TestTrackerClassBase> fRunTests; ///< test filter
    
    void printAuxDetSummary();
    void testAuxDet();
 
    bool shouldRunTests(std::string test_name) const;
       
  };
  
  
  namespace details {
    /// Class telling whether a test needs to be run
    class TestTrackerClassBase {
        public:
      using TestList_t = std::set<std::string>;
      
      virtual ~TestTrackerClassBase() = default;
      
      /// Returns whether the specified test should run
      virtual bool ShouldRun(std::string test_name) const = 0;
      
      /// Checks the test and records the request
      bool operator() (std::string test_name);
      
      /// Allow the specified test to run
      virtual void PleaseRunAlso(std::string test_name) = 0;
      
      /// Returns the tests that have been run
      TestList_t const& RunTests() const { return run; }
      
      /// Returns the tests that have been skipped
      TestList_t const& SkippedTests() const { return skipped; }
      
      /// Returns the tests that have been queried
      TestList_t QueriedTests() const;
      
      /// Checks that the validity of the configuration (after the fact)
      virtual bool CheckQueriesRegistry() const;
      
      /// Prints information about the configuration of the filter
      virtual void PrintConfiguration(std::ostream&) const;
      
        protected:
      TestList_t run; ///< requested tests that should be run
      TestList_t skipped; ///< requested tests that should be skipped
      
      virtual void RecordRequest(std::string test_name, bool bRun);
      
      /// Checks the test and records the request
      virtual bool Query(std::string test_name);
      
      /// Adds a vector of tests into a test set
      static void CopyList
        (TestList_t& dest, std::vector<std::string> const& from);
    }; // class TestTrackerClassBase
    
  } // namespace details
  
  
} // namespace geo

#endif // GEO_AUXDETGEOMETRYTESTALG_H
