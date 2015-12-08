//////////////////////////////////////////////////////////////
// Name:      DatabaseExample_module.cc
// Date:      13 July 2015
// Author:    Everybody is an author! Except for Andrzej.
//////////////////////////////////////////////////////////////
// Example module that uses the DatabaseUtilityT1034 service.
//////////////////////////////////////////////////////////////

#ifndef DatabaseExample_Module
#define DatabaseExample_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// LArIATSoft includes
#include "Utilities/DatabaseUtilityT1034.h"

// C++ includes
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace DatabaseExample {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class DatabaseExample : public art::EDAnalyzer 
  {
  public:

    // Standard constructor and destructor for an ART module.
    explicit DatabaseExample(fhicl::ParameterSet const& pset);
    virtual ~DatabaseExample();

    // This method is called once, at the start of the job.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method is called once, at the start of each sub-run. It's a
    // good place to read databases or files that may have
    // sub-run-dependent information.
    void beginSubRun(const art::SubRun& subrun);

    // This method reads in any parameters from the .fcl files.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

  private:

    // meh
    int fEvent;
    int fRun;
    int fSubRun;

    // run timestamp
    std::uint32_t fRunTimestamp;  // Unix timestamp
    std::string   fRunDateTime;   // datetime string

    // convert Unix timestamp to string with format 'YYYY-MM-DD HH24:MI:SS'
    std::string TimestampToString(std::time_t const& Timestamp);

    // DatabaseUtilityT1034 service handle
    art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

    // single DAQ XML configuration parameter from the
    // lariat_xml_database table for a specified run
    int fV1740RecordLength;

    // single IFBeam parameter from the lariat_ifbeam_database
    // table for a specified sub-run
    double fCathodeVoltage;

    // NOTE: Check the LArIAT Run Summary page for a list
    //       of parameters that are stored in the database
    //       tables:
    //
    //           http://lariat-wbm.fnal.gov/wbm/servlet/LariatRunSummary
    //
    //       You will have to select a run to list the
    //       available parameters.

    // IMPORTANT: The lariat_xml_database table has one row
    //            for each run, so it only makes sense to
    //            query this table in the beginRun() method.
    //            The lariat_ifbeam_database table has one row
    //            for each SUB-RUN, so this table should be
    //            queried in the beginSubRun() method.

    // Vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run.
    // This vector is filled in the beginJob() method.
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

    //////////////////////////////////////////////////////////
    // (key, value) pairs for database query results
    //////////////////////////////////////////////////////////
    // container for entire row from the lariat_xml_database
    // table for a specified run
    std::map< std::string, std::string > fAllConfigValues;
    // containter for entire row from the lariat_ifbeam_database
    // table for a specified sub-run
    std::map< std::string, std::string > fAllIFBeamValues;
    // containter for entire row from the lariat_hardware_connections
    // table for a specified datetime
    std::map< std::string, std::string > fHardwareConnectionsValues;
    //////////////////////////////////////////////////////////

  }; // class DatabaseExample


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  DatabaseExample::DatabaseExample(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    // read in the parameters from the .fcl file
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  // destructor
  DatabaseExample::~DatabaseExample() 
  {}

  //-----------------------------------------------------------------------
  void DatabaseExample::beginJob()
  {
    // add parameters to be queried from the lariat_xml_database table
    fConfigParams.push_back("larasic_config_larasic_induction_gain");
    fConfigParams.push_back("larasic_config_larasic_collection_gain");

    for (size_t i = 0; i < 16; ++i) {
      std::string index = std::to_string(i);
      fConfigParams.push_back("v1495_config_v1495_trigpat" + index + "_on");
      fConfigParams.push_back("v1495_config_v1495_trigpat" + index + "_off");
      fConfigParams.push_back("v1495_config_v1495_in" + index + "_name");
    }
  }

  //-----------------------------------------------------------------------
  void DatabaseExample::beginRun(const art::Run& run)
  {
    // query lariat_xml_database table here

    //========================================================
    // Start database shenanigans.
    //========================================================

    fRun = run.run();

    ///////////////////////////////////////////////////////////////////
    // get single parameter value from the lariat_xml_database table
    ///////////////////////////////////////////////////////////////////

    // results from the database are returned as strings
    // cast/convert result from string to int
    fV1740RecordLength = std::stoi(
      fDatabaseUtility->GetConfigValue("v1740_config_caen_recordlength", fRun));

    // print out the result
    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nv1740_config_caen_recordlength: " << fV1740RecordLength << "\n"
      << "///////////////////////////////////////////////////////////////////";

    ///////////////////////////////////////////////////////////////////
    // dump specified columns of lariat_xml_database for specified run
    ///////////////////////////////////////////////////////////////////

    // get results as a map where the key is the parameter name and the
    // mapped value is the parameter value
    fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams, fRun);

    // let's iterate over the results and print them out
    std::map< std::string, std::string >::const_iterator config_iter;

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nQuery result from lariat_xml_database table for run " << fRun << ".\n"
      << "///////////////////////////////////////////////////////////////////";

    for (config_iter = fConfigValues.begin();
         config_iter != fConfigValues.end();
         ++config_iter) {

      mf::LogVerbatim("DatabaseExample")
        << "Column: " << config_iter->first << "; Value: " << config_iter->second;

    }

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////";

    ///////////////////////////////////////////////////////////////////
    // dump an entire row of lariat_xml_database for specified run
    ///////////////////////////////////////////////////////////////////

    // get results as a map where the key is the parameter name and the
    // mapped value is the parameter value
    fAllConfigValues = fDatabaseUtility->GetAllConfigValues(fRun);

    // iterate over the results and print them out
    std::map< std::string, std::string >::const_iterator all_config_iter;

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nBegin dump of lariat_xml_database table for run " << fRun << "\n"
      << "///////////////////////////////////////////////////////////////////";

    for (all_config_iter = fAllConfigValues.begin();
         all_config_iter != fAllConfigValues.end();
         ++all_config_iter) {

      mf::LogVerbatim("DatabaseExample")
        << "Column: " << all_config_iter->first << "; Value: " << all_config_iter->second;

    }

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nEnd dump.\n"
      << "///////////////////////////////////////////////////////////////////";

    ///////////////////////////////////////////////////////////////////
    // dump an entire row of lariat_hardware_connections for closest
    // datetime before specified datetime
    ///////////////////////////////////////////////////////////////////

    // Unix timestamp
    fRunTimestamp = run.beginTime().timeLow();  // Unix time
    // convert Unix timestamp to string with format 'YYYY-MM-DD HH24:MI:SS'
    fRunDateTime = this->TimestampToString(fRunTimestamp);

    std::cout << "fRunTimestamp: " << fRunTimestamp << std::endl;
    std::cout << "fRunDateTime: "  << fRunDateTime  << std::endl;

    // get results as a map where the key is the parameter name and the
    // mapped value is the parameter value
    fHardwareConnectionsValues = fDatabaseUtility->GetHardwareConnections(fRunDateTime);
    // fRunDateTime format should be 'YYYY-MM-DD HH24:MI:SS'
    // e.g., '2015-06-17 14:16:00'

    // iterate over the results and print them out
    std::map< std::string, std::string >::const_iterator all_connections_iter;

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nBegin dump of lariat_hardware_connections table for " << fRunDateTime << "\n"
      << "///////////////////////////////////////////////////////////////////";

    for (all_connections_iter = fHardwareConnectionsValues.begin();
         all_connections_iter != fHardwareConnectionsValues.end();
         ++all_connections_iter) {

      mf::LogVerbatim("DatabaseExample")
        << "Column: " << all_connections_iter->first << "; Value: " << all_connections_iter->second;

    }

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nEnd dump.\n"
      << "///////////////////////////////////////////////////////////////////";

    //========================================================
    // End database shenanigans.
    //========================================================
  }

  //-----------------------------------------------------------------------
  void DatabaseExample::beginSubRun(const art::SubRun& subrun)
  {
    // query lariat_ifbeam_database table here

    //========================================================
    // Start database shenanigans.
    //========================================================

    fRun = subrun.run();
    fSubRun = subrun.subRun();

    ///////////////////////////////////////////////////////////////////
    // get single parameter value from the lariat_ifbeam_database table
    ///////////////////////////////////////////////////////////////////

    // results from the database are returned as strings
    // cast/convert result from string to double
    fCathodeVoltage = std::stod(fDatabaseUtility->GetIFBeamValue("mid_e_gmv", fRun, fSubRun));

    // print out the result
    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nmid_e_gmv: " << fCathodeVoltage << "\n"
      << "///////////////////////////////////////////////////////////////////";

    ///////////////////////////////////////////////////////////////////
    // dump an entire row of lariat_ifbeam_database for specified sub-run
    ///////////////////////////////////////////////////////////////////

    // get results as a map where the key is the parameter name and the
    // mapped value is the parameter value
    fAllIFBeamValues = fDatabaseUtility->GetAllIFBeamValues(fRun, fSubRun);

    // iterate over the results and print them out
    std::map< std::string, std::string >::const_iterator all_ifbeam_iter;

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nBegin dump of lariat_ifbeam_database table for run " << fRun << ", sub-run " << fSubRun << ".\n"
      << "///////////////////////////////////////////////////////////////////";

    for (all_ifbeam_iter = fAllIFBeamValues.begin();
         all_ifbeam_iter != fAllIFBeamValues.end();
         ++all_ifbeam_iter) {

      mf::LogVerbatim("DatabaseExample")
        << "Column: " << all_ifbeam_iter->first << "; Value: " << all_ifbeam_iter->second;

    }

    mf::LogVerbatim("DatabaseExample")
      << "///////////////////////////////////////////////////////////////////"
      << "\nEnd dump.\n"
      << "///////////////////////////////////////////////////////////////////";

    //========================================================
    // End database shenanigans.
    //========================================================
  }

  //-----------------------------------------------------------------------
  void DatabaseExample::reconfigure(fhicl::ParameterSet const& p)
  {
    return;
  }

  //-----------------------------------------------------------------------
  void DatabaseExample::analyze(const art::Event& event) 
  {
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    // do stuff with the art::Event here

    return;
  }

  //-----------------------------------------------------------------------
  std::string DatabaseExample::TimestampToString(std::time_t const& Timestamp) {
    struct tm * TimeInfo;
    char Buffer[30];
    TimeInfo = std::localtime(&Timestamp);
    std::strftime(Buffer, 30, "%Y-%m-%d %H:%M", TimeInfo);
    return std::string(Buffer);
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(DatabaseExample)

} // namespace DatabaseExample

#endif // DatabaseExample_Module
