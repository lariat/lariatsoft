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
#include "RawDataUtilities/SpillWrapper.h"

// ROOT includes
#include "TTree.h"

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
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
    void analyze(const art::Event& evt); 

   private:

    std::string fRawFragmentLabel;    ///< label for module producing artdaq fragments
    std::string fRawFragmentInstance; ///< instance label for artdaq fragments

    // run number
    //int fRun;

    // number of fragments in raw data file
    size_t fNumberFragmentsPerSpill;

    // unique pointer to SpillWrapper
    std::unique_ptr<rdu::SpillWrapper> fSpillWrapper;

    // complete LariatFragment for event record
    LariatFragment * fLariatFragment;

    // clock correction algorithm
    rdu::ClockCorrectionAlg fClockCorrectionAlg;

    // TTree for TFile... T, T, T. T is the fruit of the ROOT.
    // TKabobs, TCreole. TGumbo... PineappleT, LemonT, CoconutT,
    // PepperT. TSoup, TStew, TSalad, T and Potatoes, TBurger,
    // TSandwich. Th-that's about it.
    TTree * fClockCorrectionTree;

    // variables that will go into the n-tuples
    int    fRun;
    int    fSubRun;
    int    fDeviceID;
    int    fReferenceClockDeviceID;
    double fSlope;
    double fIntercept;
    double fClockDrift;

  }; // class ClockCorrectionCheck


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  ClockCorrectionCheck::ClockCorrectionCheck(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fNumberFragmentsPerSpill(pset.get<std::size_t>("NumberFragmentsPerSpill", 4))
    , fSpillWrapper(new rdu::SpillWrapper(fNumberFragmentsPerSpill))
    , fLariatFragment(nullptr)
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

    art::ServiceHandle<art::TFileService> tfs;

    fClockCorrectionTree = tfs->make<TTree>("ClockCorrectionTree", "ClockCorrectionTree");

    fClockCorrectionTree->Branch("Run",        &fRun,        "Run/I");
    fClockCorrectionTree->Branch("SubRun",     &fSubRun,     "SubRun/I");
    fClockCorrectionTree->Branch("Slope",      &fSlope,      "Slope/D");
    fClockCorrectionTree->Branch("Intercept",  &fIntercept,  "Intercept/D");
    fClockCorrectionTree->Branch("ClockDrift", &fClockDrift, "ClockDrift/D");
    fClockCorrectionTree->Branch("DeviceID",   &fDeviceID,   "DeviceID/I");
    fClockCorrectionTree->Branch("ReferenceClockDeviceID", &fReferenceClockDeviceID, "ReferenceClockDeviceID/I");

    return;
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

    //fNumberFragmentsPerSpill = pset.get< size_t >("NumberFragmentsPerSpill", 4);

    return;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::analyze(const art::Event& event) 
  {
    fRun = event.run();
    fSubRun = event.subRun();

    // access fragments from the event record
    if (!fSpillWrapper->ready()) {
        std::cout << "Adding event to spill wrapper." << std::endl;
        fSpillWrapper->add(event);
        std::cout << "Done!" << std::endl;
    }
    if (!fSpillWrapper->ready()) {
        return;
    }

    const uint8_t * SpillDataPtr(fSpillWrapper->get());

    mf::LogInfo("EventBuilderCheck") << "Spill is " << fSpillWrapper->size() <<
      " bytes, starting at " << static_cast<const void*>(SpillDataPtr);

    //const char * bytePtr = reinterpret_cast<const char *> (fSpillWrapper->get());
    //fLariatFragment = new LariatFragment((char *) bytePtr, fSpillWrapper->size());
    fLariatFragment = new LariatFragment((char *) SpillDataPtr, fSpillWrapper->size());

    fSpillWrapper.reset(nullptr);

    std::vector< rdu::DataBlock > DataBlocks;
    std::map< unsigned int, std::vector< double > > TimeStampMap;
    fClockCorrectionAlg.GetDataBlocksTimeStampMap(fLariatFragment, DataBlocks, TimeStampMap);

    for (std::map< unsigned int, std::vector< double> >::const_iterator
         iter = TimeStampMap.begin(); iter != TimeStampMap.end(); ++iter) {
      mf::LogVerbatim("ClockCorrectionAlg") << "Device ID: "
                                            << iter->first
                                            << "; number of data blocks: "
                                            << iter->second.size();
    }

    // get clock correction parameters
    std::map< unsigned int, std::pair< double, double > > ClockCorrectionParameters;
    unsigned int ReferenceClockDeviceID;
    fClockCorrectionAlg.GetClockCorrectionParameters(TimeStampMap, ClockCorrectionParameters, ReferenceClockDeviceID);

    fReferenceClockDeviceID = static_cast <int> (ReferenceClockDeviceID);

    std::cout << "Reference clock device ID: " << ReferenceClockDeviceID << std::endl;

    for (std::map< unsigned int, std::pair< double, double > >::const_iterator
         iter = ClockCorrectionParameters.begin(); iter != ClockCorrectionParameters.end(); ++iter) {

      fDeviceID = iter->first;
      fSlope = iter->second.first;
      fIntercept = iter->second.second;
      fClockDrift = 1 - fSlope;

      std::cout << "Device ID: " << fDeviceID << std::endl;
      std::cout << "  Slope: " << fSlope << std::endl;
      std::cout << "  Intercept: " << fIntercept << std::endl;
      std::cout << "  Clock drift [s/s]: " << fClockDrift << std::endl;

      fClockCorrectionTree->Fill();
    }

    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(ClockCorrectionCheck)

} // namespace ClockCorrectionCheck

#endif // ClockCorrectionCheck_Module
