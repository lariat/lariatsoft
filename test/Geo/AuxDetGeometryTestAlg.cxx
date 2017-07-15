/**
 * @file   AuxDetGeometryTestAlg.cxx
 * @brief  Unit test for geometry functionalities: implementation file
 * @date   2011/02/17
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryTestAlg.h
 */

// our header
#include "test/Geometry/AuxDetGeometryTestAlg.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi<>
#include "Geo/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoManager.h"
#include "TStopwatch.h"

// C/C++ standard libraries
#include <cmath>
#include <vector>
#include <iterator> // std::inserter()
#include <algorithm> // std::copy()
#include <set>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <limits> // std::numeric_limits<>


namespace {
  template <typename T>
  inline T sqr(T v) { return v*v; }
  
  template <typename T>
  std::string to_string(const T& v) {
    std::ostringstream sstr;
    sstr << v;
    return sstr.str();
  } // ::to_string()
  
  
  /// Returns whether the CET exception e contains the specified category cat
  bool hasCategory(cet::exception const& e, std::string const& cat) {
    for (auto const& e_category: e.history())
      if (e_category == cat) return true;
    return false;
  } // hasCategory()
  
  
} // local namespace


namespace simple_geo {
  
  struct Point2D {
    double y = 0.;
    double z = 0.;
    
    Point2D() = default;
    Point2D(double new_y, double new_z): y(new_y), z(new_z) {}
  }; // struct Point2D
  
  Point2D operator+ (Point2D const& a, Point2D const& b)
    { return { a.y + b.y, a.z + b.z }; }
  Point2D operator* (Point2D const& p, double f)
    { return { p.y * f, p.z * f }; }
  Point2D operator/ (Point2D const& p, double f)
    { return { p.y / f, p.z / f }; }
  template <typename Stream>
  Stream& operator<< (Stream& out, Point2D const& p)
    { out << "( " << p.y << " ; " << p.z << " )"; return out; }
  
  class Area {
      public:
    Area() = default;
    
    Area(Point2D const& a, Point2D const& b)
      {
        set_sorted(min.y, max.y, a.y, b.y);
        set_sorted(min.z, max.z, a.z, b.z);
      } // Area(Point2D x2)
    
    Point2D const& Min() const { return min; }
    Point2D const& Max() const { return max; }
    Point2D Center() const { return (min + max) / 2; }
    double DeltaY() const { return Max().y - Min().y; }
    double DeltaZ() const { return Max().z - Min().z; }
    bool isEmpty() const { return (DeltaY() == 0) || (DeltaZ() == 0); }
    
    void IncludePoint(Point2D const& point)
      {
        set_min_max(min.y, max.y, point.y);
        set_min_max(min.z, max.z, point.z);
      } // Include()
    
    void Include(Area const& area)
      { IncludePoint(area.min); IncludePoint(area.max); }
    
    void Intersect(Area const& area)
      {
        set_max(min.y, area.min.y);
        set_min(max.y, area.max.y);
        set_max(min.z, area.min.z);
        set_min(max.z, area.max.z);
      }
    
      protected:
    Point2D min, max;
    
    void set_min(double& var, double val) { if (val < var) var = val; }
    void set_max(double& var, double val) { if (val > var) var = val; }
    void set_min_max(double& min_var, double& max_var, double val)
      { set_min(min_var, val); set_max(max_var, val); }
    void set_sorted(double& min_var, double& max_var, double a, double b)
      {
        if (a > b) { min_var = b; max_var = a; }
        else       { min_var = a; max_var = b; }
      }
  }; // class Area
  
  
} // namespace simple_geo


namespace geo{
  
  
  namespace details {
    
    /// Checks the test and records the request
    bool TestTrackerClassBase::operator() (std::string test_name)
      { return Query(test_name); }
    
    TestTrackerClassBase::TestList_t TestTrackerClassBase::QueriedTests() const
    {
      TestList_t all;
      std::set_union(SkippedTests().begin(), SkippedTests().end(),
        RunTests().begin(), RunTests().end(),
        std::inserter(all, all.end())
        );
      return all;
    } // QueriedTests()
    
    bool TestTrackerClassBase::CheckQueriesRegistry() const
      { return true; /* all fine */ }
    
    void TestTrackerClassBase::PrintConfiguration(std::ostream&) const {}
    
    void TestTrackerClassBase::RecordRequest(std::string test_name, bool bRun)
      { (bRun? run: skipped).insert(test_name); }
    
    bool TestTrackerClassBase::Query(std::string test_name) {
      bool bRun = ShouldRun(test_name);
      RecordRequest(test_name, bRun);
      return bRun;
    }
    
    /// Adds a vector of tests into a test set
    void TestTrackerClassBase::CopyList
      (TestList_t& dest, std::vector<std::string> const& from)
      { std::copy(from.begin(), from.end(), std::inserter(dest, dest.end())); }
    
    
    /// Asks to run all the tests
    class PassAllTestTrackerClass: public TestTrackerClassBase {
        public:
      
      /// Returns whether the specified test should run
      virtual bool ShouldRun(std::string test_name) const override
        { return true; }
      
      // everything always runs already
      virtual void PleaseRunAlso(std::string /* test_name */) override {}
      
    }; // class PassAllTestTrackerClass
    
    /// Asks to run only tests in a list
    class WhiteListTestTrackerClass: public TestTrackerClassBase {
        public:
      using TestList_t = TestTrackerClassBase::TestList_t;
      
      //@{
      /// Constructor: takes the list of tests to be skipped
      WhiteListTestTrackerClass(TestList_t skip_these):
        to_be_skipped(skip_these) {}
      WhiteListTestTrackerClass(std::vector<std::string> const& skip_these):
        to_be_skipped()
        { CopyList(to_be_skipped, skip_these); }
      //@}
      
      /// Returns whether the specified test should run
      virtual bool ShouldRun(std::string test_name) const override
        { return to_be_skipped.count(test_name) > 0; }
      
      // everything always runs already
      virtual void PleaseRunAlso(std::string test_name) override
        { to_be_skipped.erase(test_name); }
      
      virtual bool CheckQueriesRegistry() const override
        {
          TestList_t not_registered, queried = QueriedTests();
          std::set_difference(
            to_be_skipped.cbegin(), to_be_skipped.cend(),
            queried.cbegin(), queried.cend(),
            std::inserter(not_registered, not_registered.end())
            );
          if (!not_registered.empty()) {
            auto iTest = not_registered.cbegin(), tend = not_registered.cend();
            mf::LogError error("AuxDetGeometryTestAlg");
            error
              << "The configuration presents " << not_registered.size()
              << " tests that are not supported: " << *iTest;
            while (++iTest != tend) error << ", " << *iTest;
            return false;
          }
          return true;
        } // CheckQueriesRegistry()
      
      /// Prints information about the configuration of the filter
      virtual void PrintConfiguration(std::ostream& out) const override
        {
          auto iTest = to_be_skipped.cbegin(), tend = to_be_skipped.cend();
          if (iTest == tend) {
            out << "Will skip no tests.";
            return;
          }
          out << "Will skip " << to_be_skipped.size() << " tests: " << *iTest;
          while (++iTest != tend) out << ", " << *iTest;
        } // PrintConfiguration()
      
        protected:
      TestList_t to_be_skipped; ///< tests that should be skipped
      
    }; // class WhiteListTestTrackerClass
    
    /// Asks to run only tests in a list
    class BlackListTestTrackerClass: public TestTrackerClassBase {
        public:
      using TestList_t = TestTrackerClassBase::TestList_t;
      
      //@{
      /// Constructor: takes the list of tests to be skipped
      BlackListTestTrackerClass(TestList_t run_these): to_be_run(run_these) {}
      BlackListTestTrackerClass(std::vector<std::string> const& run_these):
        to_be_run()
        { CopyList(to_be_run, run_these); }
      //@}
      
      /// Returns whether the specified test should run
      virtual bool ShouldRun(std::string test_name) const override
        { return to_be_run.count(test_name) == 0; }
      
      // everything always runs already
      virtual void PleaseRunAlso(std::string test_name) override
        { to_be_run.insert(test_name); }
      
      virtual bool CheckQueriesRegistry() const override
        {
          TestList_t not_registered, queried = QueriedTests();
          std::set_difference(
            to_be_run.cbegin(), to_be_run.cend(),
            queried.cbegin(), queried.cend(),
            std::inserter(not_registered, not_registered.end())
            );
          if (!not_registered.empty()) {
            auto iTest = not_registered.cbegin(), tend = not_registered.cend();
            mf::LogError error("AuxDetGeometryTestAlg");
            error
              << "The configuration presents " << not_registered.size()
              << " tests that are not supported: " << *iTest;
            while (++iTest != tend) error << ", " << *iTest;
            return false;
          }
          return true;
        } // CheckQueriesRegistry()
      
      /// Prints information about the configuration of the filter
      virtual void PrintConfiguration(std::ostream& out) const override
        {
          auto iTest = to_be_run.cbegin(), tend = to_be_run.cend();
          if (iTest == tend) {
            out << "Will run no tests.";
            return;
          }
          out << "Will run only " << to_be_run.size() << " tests: " << *iTest;
          while (++iTest != tend) out << ", " << *iTest;
        } // PrintConfiguration()
      
        protected:
      TestList_t to_be_run; ///< tests that should be run
      
    }; // class BlackListTestTrackerClass
    
  } // namespace details
  
  
  
  //......................................................................
  AuxDetGeometryTestAlg::AuxDetGeometryTestAlg(fhicl::ParameterSet const& pset) 
    : geom(nullptr)
  {
    // initialize the list of non-fatal exceptions
    std::vector<std::string> NonFatalErrors(pset.get<std::vector<std::string>>
      ("ForgiveExceptions", std::vector<std::string>()));
    std::copy(NonFatalErrors.begin(), NonFatalErrors.end(),
      std::inserter(fNonFatalExceptions, fNonFatalExceptions.end()));
    
    // initialize the list of tests to be run
    std::vector<std::string> RunTests(pset.get<std::vector<std::string>>
      ("RunTests", std::vector<std::string>()));
    std::vector<std::string> SkipTests(pset.get<std::vector<std::string>>
      ("SkipTests", std::vector<std::string>()));
    if (!RunTests.empty() && !SkipTests.empty()) {
      throw cet::exception("AuxDetGeometryTestAlg") << "Configuration error: "
        "'RunTests' and 'SkipTests' can't be specified together.\n";
    }
    
    if (!RunTests.empty())
      fRunTests.reset(new details::WhiteListTestTrackerClass(RunTests));
    else if (!SkipTests.empty())
      fRunTests.reset(new details::BlackListTestTrackerClass(SkipTests));
    else
      fRunTests.reset(new details::PassAllTestTrackerClass());
    
    if (pset.get<bool>("CheckForOverlaps", false))
      fRunTests->PleaseRunAlso("CheckOverlaps");
    
    if (pset.get<bool>("PrintSummary", false))
      fRunTests->PleaseRunAlso("PrintSummary");
    
    std::ostringstream sstr;
    fRunTests->PrintConfiguration(sstr);
    mf::LogInfo("AuxDetGeometryTestAlg") << sstr.str();
    
  } // AuxDetGeometryTestAlg::AuxDetGeometryTestAlg()

  //......................................................................
  unsigned int AuxDetGeometryTestAlg::Run()
  {
    
    if (!geom) {
      throw cet::exception("AuxDetGeometryTestAlg") << "AuxDetGeometryTestAlg not configured: no valid geometry provided.\n";
    }
    
    unsigned int nErrors = 0; // currently unused
    
    // change the printed version number when changing the "GeometryTest" output
    LOG_VERBATIM("AuxDetGeometryTest") << "AuxDetGeometryTest version 1.0";
    
    LOG_VERBATIM("AuxDetGeometryTestInfo") << "Running on detector: '" << geom->DetectorName() << "'";
    
    try{

      if (shouldRunTests("CheckOverlaps")) {
        LOG_INFO("AuxDetGeometryTest") << "test for overlaps ...";
        gGeoManager->CheckOverlaps(1e-5);
        gGeoManager->PrintOverlaps();
        LOG_INFO("AuxDetGeometryTest") << "complete.";
      }

      if(shouldRunTests("PrintSummary") ){
	LOG_VERBATIM("AuxDetGeometryTest") << "AuxDet Summary ...";
	this->printAuxDetSummary();
	LOG_VERBATIM("AuxDetGeometryTest") << "finished";
      }

      if(shouldRunTests("FindAtPosition") ){
	LOG_VERBATIM("AuxDetGeometryTest") << "test find at position methods ...";
	this->testFindAtPosition();
	LOG_VERBATIM("AuxDetGeometryTest") << "finished";
      }

      if(shouldRunTests("ChannelMethods") ){
	LOG_VERBATIM("AuxDetGeometryTest") << "test channel mapping methods ...";
	this->testChannelMethods();
	LOG_VERBATIM("AuxDetGeometryTest") << "finished";
      }

    }
    catch (cet::exception &e) {
      mf::LogWarning("AuxDetGeometryTest") << "exception caught: \n" << e;
      if (fNonFatalExceptions.count(e.category()) == 0) throw;
    }
    
    if (!fRunTests->CheckQueriesRegistry()) {
      throw cet::exception("AuxDetGeometryTest") << "(postumous) configuration error detected!\n";
    }
    
    return nErrors;
  } // AuxDetAuxDetGeometryTestAlg::Run()


  //......................................................................  
  void AuxDetGeometryTestAlg::printAuxDetSummary()
  {
    // get the AuxDetGeo objects
    auto adGeoVec = geom->AuxDetGeoVec();

    for(auto ad : adGeoVec){
      LOG_VERBATIM("AuxDetGeometryTestAlg") << "AuxDetGeo: "                  << ad->TotalVolume()->GetName()  
					    << "\n\t length: "       	      << ad->Length()		      
					    << "\n\t half width 1: " 	      << ad->HalfWidth1()	      
					    << "\n\t half width 2: " 	      << ad->HalfWidth2()	      
					    << "\n\t half height: "  	      << ad->HalfHeight()              
					    << "\n\t num sensitive volumes: " << ad->NSensitiveVolume();
      for(size_t sv = 0; sv < ad->NSensitiveVolume(); ++sv){
	LOG_VERBATIM("AuxDetGeometryTestAlg") << "\t\t Sensitive volume: " << sv
					      << "\t\t\t length: "         << ad->SensitiveVolume(sv).Length()		      
					      << "\t\t\t half width 1: "   << ad->SensitiveVolume(sv).HalfWidth1()	      
					      << "\t\t\t half width 2: "   << ad->SensitiveVolume(sv).HalfWidth2()	      
					      << "\t\t\t half height: "    << ad->SensitiveVolume(sv).HalfHeight();            
      }  // end loop over sensitive volumes
    } // end loop over auxdetgeos

    return;
  }

  //......................................................................  
  void AuxDetGeometryTestAlg::testFindAtPosition()
  {
    // loop over all the auxiliary detectors and grab the center of each
    // detector.  For each center, call FindAuxDetAtPosition and see if 
    // it returns the correct one

    // get the AuxDetGeo objects
    auto adGeoVec    = geom->AuxDetGeoVec();
    double origin[3] = {0.};
    double world[3]  = {0.};
    size_t adTmp     = 0;
    size_t svTmp     = 0;

    for(size_t ad = 0; ad < adGeoVec.size(); ++ad){
      
      adGeoVec[ad]->LocalToWorld(origin, world);

      LOG_VERBATIM("AuxDetGeometryTestAlg") << "AuxDetGeometryCore::FindAuxDetAtPosition for "
					    << "AuxDetGeo " << ad << " with center ("
					    << world[0] << ", " << world[1] << ", " << world[2] << ")";

      if(ad != geom->FindAuxDetAtPosition(world) )
	throw cet::exception("FailAuxDetGeometry") << "AuxDetGeometryCore::FindAuxDetAtPosition failed";

      // now test the sensitive volumes for this aux det
      for(size_t s = 0; s < adGeoVec[ad]->NSensitiveVolume(); ++s){

	adGeoVec[ad]->SensitiveVolume(s).LocalToWorld(origin, world);
	
	LOG_VERBATIM("AuxDetGeometryTestAlg") << "AuxDetGeometryCore::FindAuxDetSensitiveAtPosition for "
					      << "AuxDetSensitiveGeo " << s << " with center ("
					      << world[0] << ", " << world[1] << ", " << world[2] << ")";
	
	geom->FindAuxDetSensitiveAtPosition(world, adTmp, svTmp);

	if(ad != adTmp || s != svTmp)
	  throw cet::exception("FailAuxDetGeometry") << "AuxDetGeometryCore::FindAuxDetSensitiveAtPosition failed";
		
      } // end loop over sensitive volumes

    } // end loop over auxdetgeos

    return;
  }

  //......................................................................  
  void AuxDetGeometryTestAlg::testChannelMethods()
  {
    // loop over all the auxiliary detectors and grab the center of each
    // detector.  For each center, call FindAuxDetAtPosition and see if 
    // it returns the correct one

    // get the AuxDetGeo objects
    auto     adGeoVec  = geom->AuxDetGeoVec();
    double   origin[3] = {0.};
    double   world[3]  = {0.};
    size_t   adTmp     = 0;
    size_t   svTmp     = 0;
    uint32_t chan      = 0;    

    for(size_t ad = 0; ad < adGeoVec.size(); ++ad){
      
      for(size_t s = 0; s < adGeoVec[ad]->NSensitiveVolume(); ++s){

	adGeoVec[ad]->SensitiveVolume(s).LocalToWorld(origin, world);
	
	LOG_VERBATIM("AuxDetGeometryTestAlg") << "AuxDetGeo " << ad  
					      << " AuxDetSensitiveGeo " << s << " with center ("
					      << world[0] << ", " << world[1] << ", " << world[2] << ")";
	
	chan = geom->PositionToAuxDetChannel(world, adTmp, svTmp);

	if(ad != adTmp || s != svTmp)
	  throw cet::exception("FailAuxDetGeometry") << "AuxDetGeometryCore::PositionToAuxDetChannel failed";
		
	auto adRet = geom->ChannelToAuxDet(adGeoVec[ad]->TotalVolume()->GetName(), chan);
	
	if( adRet.TotalVolume()->GetName() != adGeoVec[ad]->TotalVolume()->GetName() )
	  throw cet::exception("FailAuxDetGeometry") << "AuxDetGeometryCore::ChannelToAuxDet failed";

	auto svRet = geom->ChannelToAuxDetSensitive(adGeoVec[ad]->TotalVolume()->GetName(), chan);

	if( svRet.TotalVolume()->GetName() != adGeoVec[ad]->SensitiveVolume(sv).TotalVolume()->GetName() )
	  throw cet::exception("FailAuxDetGeometry") << "AuxDetGeometryCore::ChannelToAuxDetSensitive failed";

      } // end loop over sensitive volumes

    } // end loop over auxdetgeos

    return;
  }

  //......................................................................  
  inline bool AuxDetGeometryTestAlg::shouldRunTests(std::string test_name) const {
    return (*fRunTests.get())(test_name);
  } // AuxDetAuxDetGeometryTestAlg::shouldRunTests()

}//end namespace
