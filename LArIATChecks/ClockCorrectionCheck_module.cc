//////////////////////////////////////////////////////////////
// Name:      ClockCorrectionCzech_module.cc
// Date:      27 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef ClockCorrectionCheck_Module
#define ClockCorrectionCheck_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// LArIATSoft includes
#include "RawDataUtilities/ClockCorrectionAlg.h"
#include "RawDataUtilities/FragmentUtility.h"

// C++ includes
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ClockCorrectionCheck {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class ClockCorrectionCheck : public art::EDAnalyzer 
  {

   public:

    // Standard constructor and destructor for an ART module.
    explicit ClockCorrectionCheck(fhicl::ParameterSet const& pset);
    virtual ~ClockCorrectionCheck();

    // This method is called once, at the start of the job.
    void beginJob();

    // This method is called once, at the start of each run.
    void beginRun(const art::Run& run);

    // This method is called once, at the start of each sub-run
    void beginSubRun(const art::SubRun& subrun);

    // This method reads in any parameters from the .fcl files.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

   private:

    std::string fRawFragmentLabel;    ///< label for module producing artdaq fragments
    std::string fRawFragmentInstance; ///< instance label for artdaq fragments

    rdu::ClockCorrectionAlg fClockCorrectionAlg;

  }; // class ClockCorrectionCheck


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  ClockCorrectionCheck::ClockCorrectionCheck(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fClockCorrectionAlg(pset.get<fhicl::ParameterSet>("ClockCorrectionAlg"))
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  ClockCorrectionCheck::~ClockCorrectionCheck() 
  {}

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::beginJob()
  {
    fClockCorrectionAlg.hello_world();
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::beginRun(const art::Run& run)
  {}

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::beginSubRun(const art::SubRun& subrun)
  {}

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClockCorrectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ClockCorrectionAlg"));
    fRawFragmentLabel    = pset.get< std::string >("RawFragmentLabel",    "daq"  );
    fRawFragmentInstance = pset.get< std::string >("RawFragmentInstance", "SPILL");
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::analyze(const art::Event& event) 
  {

    // make the utility to access the fragments from the event record
    rdu::FragmentUtility fragUtil(event, fRawFragmentLabel, fRawFragmentInstance);

    std::vector<rdu::DataBlockCollection> collections;

    collections = fClockCorrectionAlg.GroupCollections(&fragUtil.DAQFragment());

    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(ClockCorrectionCheck)

} // namespace ClockCorrectionCheck

#endif // ClockCorrectionCheck_Module
