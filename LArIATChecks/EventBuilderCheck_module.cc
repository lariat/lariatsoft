//////////////////////////////////////////////////////////////
// Name:      EventBuilderCzech_module.cc
// Date:      22 August 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef EventBuilderCheck_Module
#define EventBuilderCheck_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// LArIATFragments
#include "LArIATFragments/LariatFragment.h"

// LArIATSoft includes
#include "RawDataUtilities/EventBuilderAlg.h"
#include "RawDataUtilities/SpillWrapper.h"
#include "Utilities/DatabaseUtilityT1034.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace EventBuilderCheck {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class EventBuilderCheck : public art::EDAnalyzer 
  {

   public:

    // Standard constructor and destructor for an ART module.
    explicit EventBuilderCheck(fhicl::ParameterSet const& pset);
    virtual ~EventBuilderCheck();

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
    int fRun;

    // DatabaseUtilityT1034 service handle
    art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

    // number of fragments in raw data file
    size_t fNumberFragmentsPerSpill;

    // unique pointer to SpillWrapper
    std::unique_ptr<rdu::SpillWrapper> fSpillWrapper;

    // complete LariatFragment for event record
    LariatFragment * fLariatFragment;

    // event builder algorithm
    rdu::EventBuilderAlg fEventBuilderAlg;

    // vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

    size_t CastToSizeT(std::string const& String);

    // V1495 delay for the delayed trigger
    size_t fV1495DelayTicks;
    double fV1495Delay;

    // sample reduction of the CAEN V1740 digitizers
    size_t fV1740SampleReduction;
    size_t fV1740BSampleReduction;

    // post percent of CAEN digitizers
    double fV1740PostPercent;
    double fV1740BPostPercent;
    double fV1751PostPercent;

    // record lengths (number of time ticks) of the CAEN
    // V1740 and V1751 digitizers
    size_t fV1740RecordLength;
    size_t fV1740BRecordLength;
    size_t fV1751RecordLength;

    // sampling rate in MHz
    double fV1740SamplingRate;
    double fV1740BSamplingRate;
    double fV1751SamplingRate;

    // readout window lengths in microseconds
    double fV1740ReadoutWindow;
    double fV1740BReadoutWindow;
    double fV1751ReadoutWindow;
    double fTDCReadoutWindow;

    // TDC parameters
    double fTDCPipelineDelay;
    double fTDCGateWidth;

    // pre-/post-acquisition and acquisition windows
    double fV1740PreAcquisitionWindow;
    double fV1740PostAcquisitionWindow;
    double fV1740AcquisitionWindow;
    double fV1740BPreAcquisitionWindow;
    double fV1740BPostAcquisitionWindow;
    double fV1740BAcquisitionWindow;
    double fV1751PreAcquisitionWindow;
    double fV1751PostAcquisitionWindow;
    double fV1751AcquisitionWindow;
    double fTDCPreAcquisitionWindow;
    double fTDCPostAcquisitionWindow;
    double fTDCAcquisitionWindow;

    // pointer to histograms
    TH1D * fIntervalsDeltaTHistogram;
    TH1D * fIntervalsDeltaTZHistogram;
    TH1D * fNumberTPCReadoutsHistogram;
    TH1D * fTPCIntervalsDeltaTHistogram;
    TH1D * fTPCIntervalsDeltaTZHistogram;

    // pointer to n-tuples
    TTree * fEventBuilderTree;
    TTree * fTPCTree;

    // variables that will go into the n-tuples
    int fSubRun;
    double fIntervalsDeltaT;
    double fIntervalLength;
    double fIntervalStart;
    double fIntervalStop;
    int fNumberTPCReadouts;
    int fNumberCAENBlocks;
    int fNumberCAENBoard0Blocks;
    int fNumberCAENBoard1Blocks;
    int fNumberCAENBoard2Blocks;
    int fNumberCAENBoard3Blocks;
    int fNumberCAENBoard4Blocks;
    int fNumberCAENBoard5Blocks;
    int fNumberCAENBoard6Blocks;
    int fNumberCAENBoard7Blocks;
    int fNumberCAENBoard8Blocks;
    int fNumberCAENBoard9Blocks;
    int fNumberCAENBoard24Blocks;
    int fNumberTDCBlocks;
    double fTPCIntervalsDeltaT;
    std::vector<double> fCAENBoard0TimeStamps;
    std::vector<double> fCAENBoard1TimeStamps;
    std::vector<double> fCAENBoard2TimeStamps;
    std::vector<double> fCAENBoard3TimeStamps;
    std::vector<double> fCAENBoard4TimeStamps;
    std::vector<double> fCAENBoard5TimeStamps;
    std::vector<double> fCAENBoard6TimeStamps;
    std::vector<double> fCAENBoard7TimeStamps;
    std::vector<double> fCAENBoard8TimeStamps;
    std::vector<double> fCAENBoard9TimeStamps;
    std::vector<double> fCAENBoard24TimeStamps;
    std::vector<double> fTDCTimeStamps;

  }; // class EventBuilderCheck


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  EventBuilderCheck::EventBuilderCheck(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fNumberFragmentsPerSpill(pset.get<std::size_t>("NumberFragmentsPerSpill", 4))
    , fSpillWrapper(new rdu::SpillWrapper(fNumberFragmentsPerSpill))
    , fLariatFragment(nullptr)
    , fEventBuilderAlg(pset.get<fhicl::ParameterSet>("EventBuilderAlg"))
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  EventBuilderCheck::~EventBuilderCheck() 
  {}

  //-----------------------------------------------------------------------
  void EventBuilderCheck::beginJob()
  {
    fEventBuilderAlg.hello_world();

    fConfigParams.push_back("v1495_config_v1495_delay_ticks");
    fConfigParams.push_back("v1740_config_caen_postpercent");
    fConfigParams.push_back("v1740_config_caen_recordlength");
    fConfigParams.push_back("v1740b_config_caen_postpercent");
    fConfigParams.push_back("v1740b_config_caen_recordlength");
    fConfigParams.push_back("v1751_config_caen_postpercent");
    fConfigParams.push_back("v1751_config_caen_recordlength");
    fConfigParams.push_back("v1740_config_caen_v1740_samplereduction");
    fConfigParams.push_back("v1740b_config_caen_v1740_samplereduction");
    fConfigParams.push_back("tdc_config_tdc_pipelinedelay");
    fConfigParams.push_back("tdc_config_tdc_gatewidth");

    art::ServiceHandle<art::TFileService> tfs;

    fEventBuilderTree = tfs->make<TTree>("EventBuilderTree", "EventBuilderTree");

    fEventBuilderTree->Branch("Run",                         &fRun,                         "Run/I");
    fEventBuilderTree->Branch("SubRun",                      &fSubRun,                      "SubRun/I");
    fEventBuilderTree->Branch("IntervalsDeltaT",             &fIntervalsDeltaT,             "IntervalsDeltaT/D");
    fEventBuilderTree->Branch("IntervalLength",              &fIntervalLength,              "IntervalLength/D");
    fEventBuilderTree->Branch("IntervalStart",               &fIntervalStart,               "IntervalStart/D");
    fEventBuilderTree->Branch("IntervalStop",                &fIntervalStop,                "IntervalStop/D");
    fEventBuilderTree->Branch("NumberTPCReadouts",           &fNumberTPCReadouts,           "NumberTPCReadouts/I");
    fEventBuilderTree->Branch("NumberCAENBlocks",            &fNumberCAENBlocks,            "NumberCAENBlocks/I");
    fEventBuilderTree->Branch("NumberTDCBlocks",             &fNumberTDCBlocks,             "NumberTDCBlocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard0Blocks",      &fNumberCAENBoard0Blocks,      "NumberCAENBoard0Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard1Blocks",      &fNumberCAENBoard1Blocks,      "NumberCAENBoard1Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard2Blocks",      &fNumberCAENBoard2Blocks,      "NumberCAENBoard2Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard3Blocks",      &fNumberCAENBoard3Blocks,      "NumberCAENBoard3Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard4Blocks",      &fNumberCAENBoard4Blocks,      "NumberCAENBoard4Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard5Blocks",      &fNumberCAENBoard5Blocks,      "NumberCAENBoard5Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard6Blocks",      &fNumberCAENBoard6Blocks,      "NumberCAENBoard6Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard7Blocks",      &fNumberCAENBoard7Blocks,      "NumberCAENBoard7Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard8Blocks",      &fNumberCAENBoard8Blocks,      "NumberCAENBoard8Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard9Blocks",      &fNumberCAENBoard9Blocks,      "NumberCAENBoard9Blocks/I");
    fEventBuilderTree->Branch("NumberCAENBoard24Blocks",     &fNumberCAENBoard24Blocks,     "NumberCAENBoard24Blocks/I");
    fEventBuilderTree->Branch("V1740PreAcquisitionWindow",   &fV1740PreAcquisitionWindow,   "V1740PreAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1740PostAcquisitionWindow",  &fV1740PostAcquisitionWindow,  "V1740PostAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1740AcquisitionWindow",      &fV1740AcquisitionWindow,      "V1740AcquisitionWindow/D");
    fEventBuilderTree->Branch("V1740BPreAcquisitionWindow",  &fV1740BPreAcquisitionWindow,  "V1740BPreAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1740BPostAcquisitionWindow", &fV1740BPostAcquisitionWindow, "V1740BPostAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1740BAcquisitionWindow",     &fV1740BAcquisitionWindow,     "V1740BAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1751PreAcquisitionWindow",   &fV1751PreAcquisitionWindow,   "V1751PreAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1751PostAcquisitionWindow",  &fV1751PostAcquisitionWindow,  "V1751PostAcquisitionWindow/D");
    fEventBuilderTree->Branch("V1751AcquisitionWindow",      &fV1751AcquisitionWindow,      "V1751AcquisitionWindow/D");
    fEventBuilderTree->Branch("TDCPreAcquisitionWindow",     &fTDCPreAcquisitionWindow,     "TDCPreAcquisitionWindow/D");
    fEventBuilderTree->Branch("TDCPostAcquisitionWindow",    &fTDCPostAcquisitionWindow,    "TDCPostAcquisitionWindow/D");
    fEventBuilderTree->Branch("TDCAcquisitionWindow",        &fTDCAcquisitionWindow,        "TDCAcquisitionWindow/D");
    fEventBuilderTree->Branch("CAENBoard0TimeStamps",        &fCAENBoard0TimeStamps);
    fEventBuilderTree->Branch("CAENBoard1TimeStamps",        &fCAENBoard1TimeStamps);
    fEventBuilderTree->Branch("CAENBoard2TimeStamps",        &fCAENBoard2TimeStamps);
    fEventBuilderTree->Branch("CAENBoard3TimeStamps",        &fCAENBoard3TimeStamps);
    fEventBuilderTree->Branch("CAENBoard4TimeStamps",        &fCAENBoard4TimeStamps);
    fEventBuilderTree->Branch("CAENBoard5TimeStamps",        &fCAENBoard5TimeStamps);
    fEventBuilderTree->Branch("CAENBoard5TimeStamps",        &fCAENBoard5TimeStamps);
    fEventBuilderTree->Branch("CAENBoard6TimeStamps",        &fCAENBoard6TimeStamps);
    fEventBuilderTree->Branch("CAENBoard7TimeStamps",        &fCAENBoard7TimeStamps);
    fEventBuilderTree->Branch("CAENBoard8TimeStamps",        &fCAENBoard8TimeStamps);
    fEventBuilderTree->Branch("CAENBoard9TimeStamps",        &fCAENBoard9TimeStamps);
    fEventBuilderTree->Branch("CAENBoard24TimeStamps",       &fCAENBoard24TimeStamps);
    fEventBuilderTree->Branch("TDCTimeStamps",               &fTDCTimeStamps);

    fTPCTree = tfs->make<TTree>("TPCTree", "TPCTree");
    fTPCTree->Branch("Run",                &fRun,                "Run/I");
    fTPCTree->Branch("SubRun",             &fSubRun,             "SubRun/I");
    fTPCTree->Branch("TPCIntervalsDeltaT", &fTPCIntervalsDeltaT, "TPCIntervalsDeltaT/D");

    fIntervalsDeltaTHistogram      = tfs->make<TH1D>("IntervalsDeltaT",     ";#Delta t [ms];Entries per ms", 10000, -0.5, 9999.5);
    fIntervalsDeltaTZHistogram     = tfs->make<TH1D>("IntervalsDeltaTZ",    ";#Delta t [ms];Entries per 0.001 ms", 1000, 0, 1);
    fNumberTPCReadoutsHistogram    = tfs->make<TH1D>("NumberTPCReadouts",   ";Number of TPC readouts;Entries", 10, -0.5, 9.5);
    fTPCIntervalsDeltaTHistogram   = tfs->make<TH1D>("TPCIntervalsDeltaT",  ";#Delta t [ms];Entries per ms", 10000, -0.5, 9999.5);
    fTPCIntervalsDeltaTZHistogram  = tfs->make<TH1D>("TPCIntervalsDeltaTZ", ";#Delta t [ms];Entries per 0.001 ms", 1000, 0, 1);

    return;
  }

  //-----------------------------------------------------------------------
  void EventBuilderCheck::beginRun(const art::Run& run)
  {
    fRun = run.run();
    fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams, fRun);

    std::cout << "//////////////////////////////////////////////"                                       << std::endl;
    std::cout << "V1495DelayTicks:       " << fConfigValues["v1495_config_v1495_delay_ticks"]           << std::endl;
    std::cout << "V1740PostPercent:      " << fConfigValues["v1740_config_caen_postpercent"]            << std::endl;
    std::cout << "V1740BPostPercent:     " << fConfigValues["v1740b_config_caen_postpercent"]           << std::endl;
    std::cout << "V1751PostPercent:      " << fConfigValues["v1751_config_caen_postpercent"]            << std::endl;
    std::cout << "V1740RecordLength:     " << fConfigValues["v1740_config_caen_recordlength"]           << std::endl;
    std::cout << "V1740BRecordLength:    " << fConfigValues["v1740b_config_caen_recordlength"]          << std::endl;
    std::cout << "V1751RecordLength:     " << fConfigValues["v1751_config_caen_recordlength"]           << std::endl;
    std::cout << "V1740SampleReduction:  " << fConfigValues["v1740_config_caen_v1740_samplereduction"]  << std::endl;
    std::cout << "V1740BSampleReduction: " << fConfigValues["v1740b_config_caen_v1740_samplereduction"] << std::endl;
    std::cout << "fTDCPipelineDelay:     " << fConfigValues["tdc_config_tdc_pipelinedelay"]             << std::endl;
    std::cout << "fTDCGateWidth:         " << fConfigValues["tdc_config_tdc_gatewidth"]                 << std::endl;
    std::cout << "//////////////////////////////////////////////"                                       << std::endl;

    //fV1495DelayTicks       = static_cast <size_t> (std::stoi(fConfigValues["v1495_config_v1495_delay_ticks"]));
    //fV1740PostPercent      = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_postpercent"]));
    //fV1740BPostPercent     = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_postpercent"]));
    //fV1751PostPercent      = static_cast <size_t> (std::stoi(fConfigValues["v1751_config_caen_postpercent"]));
    //fV1740RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_recordlength"]));
    //fV1740BRecordLength    = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_recordlength"]));
    //fV1751RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1751_config_caen_recordlength"]));
    //fV1740SampleReduction  = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_v1740_samplereduction"]));
    //fV1740BSampleReduction = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_v1740_samplereduction"]));
    //fTDCPipelineDelay      = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_pipelinedelay"]));
    //fTDCGateWidth          = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_gatewidth"]));

    fV1495DelayTicks       = this->CastToSizeT(fConfigValues["v1495_config_v1495_delay_ticks"]);
    fV1740PostPercent      = this->CastToSizeT(fConfigValues["v1740_config_caen_postpercent"]);
    fV1740BPostPercent     = this->CastToSizeT(fConfigValues["v1740b_config_caen_postpercent"]);
    fV1751PostPercent      = this->CastToSizeT(fConfigValues["v1751_config_caen_postpercent"]);
    fV1740RecordLength     = this->CastToSizeT(fConfigValues["v1740_config_caen_recordlength"]);
    fV1740BRecordLength    = this->CastToSizeT(fConfigValues["v1740b_config_caen_recordlength"]);
    fV1751RecordLength     = this->CastToSizeT(fConfigValues["v1751_config_caen_recordlength"]);
    fV1740SampleReduction  = this->CastToSizeT(fConfigValues["v1740_config_caen_v1740_samplereduction"]);
    fV1740BSampleReduction = this->CastToSizeT(fConfigValues["v1740b_config_caen_v1740_samplereduction"]);
    fTDCPipelineDelay      = this->CastToSizeT(fConfigValues["tdc_config_tdc_pipelinedelay"]);
    fTDCGateWidth          = this->CastToSizeT(fConfigValues["tdc_config_tdc_gatewidth"]);

    // V1495 delay
    fV1495Delay = static_cast <double> (fV1495DelayTicks) * 0.01;

    // sampling rate in MHz
    fV1740SamplingRate  = 62.5 / fV1740SampleReduction;
    fV1740BSamplingRate = 62.5 / fV1740BSampleReduction;
    fV1751SamplingRate  = 1e3;

    // readout window length in microseconds
    fV1740ReadoutWindow  = fV1740RecordLength  / fV1740SamplingRate;
    fV1740BReadoutWindow = fV1740BRecordLength / fV1740BSamplingRate;
    fV1751ReadoutWindow  = fV1751RecordLength  / fV1751SamplingRate;
    fTDCReadoutWindow    = fTDCGateWidth * 8 * 0.001177;

    // pre-/post-acquisition and acquisition windows
    if (fV1740PreAcquisitionWindow   < 0) fV1740PreAcquisitionWindow   = fV1740ReadoutWindow;
    if (fV1740PostAcquisitionWindow  < 0) fV1740PostAcquisitionWindow  = 0.128;
    if (fV1740AcquisitionWindow      < 0) fV1740AcquisitionWindow      = fV1740ReadoutWindow;
    if (fV1740BPreAcquisitionWindow  < 0) fV1740BPreAcquisitionWindow  = fV1740BReadoutWindow;
    if (fV1740BPostAcquisitionWindow < 0) fV1740BPostAcquisitionWindow = 0.128;
    if (fV1740BAcquisitionWindow     < 0) fV1740BAcquisitionWindow     = fV1740BReadoutWindow;
    if (fV1751PreAcquisitionWindow   < 0) fV1751PreAcquisitionWindow   = 0.128;
    if (fV1751PostAcquisitionWindow  < 0) fV1751PostAcquisitionWindow  = 0.128;
    if (fV1751AcquisitionWindow      < 0) fV1751AcquisitionWindow      = fV1751ReadoutWindow;
    if (fTDCPreAcquisitionWindow     < 0) fTDCPreAcquisitionWindow     = fTDCReadoutWindow;
    if (fTDCPostAcquisitionWindow    < 0) fTDCPostAcquisitionWindow    = 0;
    if (fTDCAcquisitionWindow        < 0) fTDCAcquisitionWindow        = fTDCReadoutWindow;

    std::cout << "//////////////////////////////////////////////"                << std::endl;
    std::cout << "V1495DelayTicks:             " << fV1495DelayTicks             << std::endl;
    std::cout << "V1495Delay:                  " << fV1495Delay                  << std::endl;
    std::cout << "V1740PreAcquisitionWindow:   " << fV1740PreAcquisitionWindow   << std::endl;
    std::cout << "V1740PostAcquisitionWindow:  " << fV1740PostAcquisitionWindow  << std::endl;
    std::cout << "V1740AcquisitionWindow:      " << fV1740AcquisitionWindow      << std::endl;
    std::cout << "V1740BPreAcquisitionWindow:  " << fV1740BPreAcquisitionWindow  << std::endl;
    std::cout << "V1740BPostAcquisitionWindow: " << fV1740BPostAcquisitionWindow << std::endl;
    std::cout << "V1740BAcquisitionWindow:     " << fV1740BAcquisitionWindow     << std::endl;
    std::cout << "V1751PreAcquisitionWindow:   " << fV1751PreAcquisitionWindow   << std::endl;
    std::cout << "V1751PostAcquisitionWindow:  " << fV1751PostAcquisitionWindow  << std::endl;
    std::cout << "V1751AcquisitionWindow:      " << fV1751AcquisitionWindow      << std::endl;
    std::cout << "TDCPreAcquisitionWindow:     " << fTDCPreAcquisitionWindow     << std::endl;
    std::cout << "TDCPostAcquisitionWindow:    " << fTDCPostAcquisitionWindow    << std::endl;
    std::cout << "TDCAcquisitionWindow:        " << fTDCAcquisitionWindow        << std::endl;
    std::cout << "//////////////////////////////////////////////"                << std::endl;

    return;
  }

  //-----------------------------------------------------------------------
  void EventBuilderCheck::beginSubRun(const art::SubRun& subrun)
  {}

  //-----------------------------------------------------------------------
  void EventBuilderCheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    fEventBuilderAlg.reconfigure(pset.get<fhicl::ParameterSet>("EventBuilderAlg"));
    fRawFragmentLabel    = pset.get< std::string >("RawFragmentLabel",    "daq"  );
    fRawFragmentInstance = pset.get< std::string >("RawFragmentInstance", "SPILL");

    //fNumberFragmentsPerSpill = pset.get< size_t >("NumberFragmentsPerSpill", 4);

    fV1740PreAcquisitionWindow   = pset.get< double >("V1740PreAcquisitionWindow",   -1);
    fV1740PostAcquisitionWindow  = pset.get< double >("V1740PostAcquisitionWindow",  -1);
    fV1740AcquisitionWindow      = pset.get< double >("V1740AcquisitionWindow",      -1);
    fV1740BPreAcquisitionWindow  = pset.get< double >("V1740BPreAcquisitionWindow",  -1);
    fV1740BPostAcquisitionWindow = pset.get< double >("V1740BPostAcquisitionWindow", -1);
    fV1740BAcquisitionWindow     = pset.get< double >("V1740BAcquisitionWindow",     -1);
    fV1751PreAcquisitionWindow   = pset.get< double >("V1751PreAcquisitionWindow",   -1);
    fV1751PostAcquisitionWindow  = pset.get< double >("V1751PostAcquisitionWindow",  -1);
    fV1751AcquisitionWindow      = pset.get< double >("V1751AcquisitionWindow",      -1);
    fTDCPreAcquisitionWindow     = pset.get< double >("TDCPreAcquisitionWindow",     -1);
    fTDCPostAcquisitionWindow    = pset.get< double >("TDCPostAcquisitionWindow",    -1);
    fTDCAcquisitionWindow        = pset.get< double >("TDCAcquisitionWindow",        -1);

    return;
  }

  //-----------------------------------------------------------------------
  void EventBuilderCheck::analyze(const art::Event& event) 
  {
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

    // configure the event builder algorithm
    fEventBuilderAlg.Configure(fV1740PreAcquisitionWindow,
                               fV1740PostAcquisitionWindow,
                               fV1740AcquisitionWindow,
                               fV1740BPreAcquisitionWindow,
                               fV1740BPostAcquisitionWindow,
                               fV1740BAcquisitionWindow,
                               fV1751PreAcquisitionWindow,
                               fV1751PostAcquisitionWindow,
                               fV1751AcquisitionWindow,
                               fTDCPreAcquisitionWindow,
                               fTDCPostAcquisitionWindow,
                               fTDCAcquisitionWindow);

    // group data blocks into collections
    std::vector< rdu::DataBlockCollection > Collections;
    Collections = fEventBuilderAlg.Build(fLariatFragment);

    for (size_t i = 0; i < Collections.size(); ++i) {
      rdu::DataBlockCollection const& Collection = Collections[i];

      size_t const& NumberCaenBlocks = Collection.caenBlocks.size();
      size_t const& NumberTdcBlocks  = Collection.tdcBlocks.size();

      fNumberCAENBlocks = NumberCaenBlocks;
      fNumberTDCBlocks  = NumberTdcBlocks;

      fIntervalsDeltaT = -1;
      fIntervalLength  = -1;
      fIntervalStart   = -1;
      fIntervalStop    = -1;
      fNumberTPCReadouts = 0;

      fNumberCAENBoard0Blocks  = 0;
      fNumberCAENBoard1Blocks  = 0;
      fNumberCAENBoard2Blocks  = 0;
      fNumberCAENBoard3Blocks  = 0;
      fNumberCAENBoard4Blocks  = 0;
      fNumberCAENBoard5Blocks  = 0;
      fNumberCAENBoard6Blocks  = 0;
      fNumberCAENBoard7Blocks  = 0;
      fNumberCAENBoard8Blocks  = 0;
      fNumberCAENBoard9Blocks  = 0;
      fNumberCAENBoard24Blocks = 0;

      int NumberCAENBlocksArray[32] = {};

      fCAENBoard0TimeStamps.clear();
      fCAENBoard1TimeStamps.clear();
      fCAENBoard2TimeStamps.clear();
      fCAENBoard3TimeStamps.clear();
      fCAENBoard4TimeStamps.clear();
      fCAENBoard5TimeStamps.clear();
      fCAENBoard6TimeStamps.clear();
      fCAENBoard7TimeStamps.clear();
      fCAENBoard8TimeStamps.clear();
      fCAENBoard9TimeStamps.clear();
      fCAENBoard24TimeStamps.clear();
      fTDCTimeStamps.clear();

      std::map< unsigned int, std::vector< double > > CAENTimeStamps;

      fNumberTPCReadouts = Collection.numberTPCReadouts;
      fIntervalLength = Collection.interval.second - Collection.interval.first;
      fIntervalStart  = Collection.interval.first;
      fIntervalStop   = Collection.interval.second;

      std::cout << "Collection: " << i << std::endl;
      std::cout << "  Number of CAEN data blocks: " << NumberCaenBlocks << std::endl;
      std::cout << "  Number of TDC data blocks:  " << NumberTdcBlocks << std::endl;
      std::cout << "  Number of TPC readouts:     " << fNumberTPCReadouts << std::endl;

      for (size_t j = 0; j < NumberCaenBlocks; ++j) {

        CAENFragment const& caenFrag = Collection.caenBlocks[j];
        double const& timestamp = Collection.caenBlockTimeStamps[j];
        unsigned int boardId = caenFrag.header.boardId;

        std::cout << "    CAEN block: " << j << std::endl;
        std::cout << "      Board ID: " << boardId << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;

        NumberCAENBlocksArray[boardId] += 1;
        CAENTimeStamps[boardId].push_back(timestamp);

      }

      //

      fNumberCAENBoard0Blocks  = NumberCAENBlocksArray[0];
      fNumberCAENBoard1Blocks  = NumberCAENBlocksArray[1];
      fNumberCAENBoard2Blocks  = NumberCAENBlocksArray[2];
      fNumberCAENBoard3Blocks  = NumberCAENBlocksArray[3];
      fNumberCAENBoard4Blocks  = NumberCAENBlocksArray[4];
      fNumberCAENBoard5Blocks  = NumberCAENBlocksArray[5];
      fNumberCAENBoard6Blocks  = NumberCAENBlocksArray[6];
      fNumberCAENBoard7Blocks  = NumberCAENBlocksArray[7];
      fNumberCAENBoard8Blocks  = NumberCAENBlocksArray[8];
      fNumberCAENBoard9Blocks  = NumberCAENBlocksArray[9];
      fNumberCAENBoard24Blocks = NumberCAENBlocksArray[24];

      fCAENBoard0TimeStamps  = CAENTimeStamps[0];
      fCAENBoard1TimeStamps  = CAENTimeStamps[1];
      fCAENBoard2TimeStamps  = CAENTimeStamps[2];
      fCAENBoard3TimeStamps  = CAENTimeStamps[3];
      fCAENBoard4TimeStamps  = CAENTimeStamps[4];
      fCAENBoard5TimeStamps  = CAENTimeStamps[5];
      fCAENBoard6TimeStamps  = CAENTimeStamps[6];
      fCAENBoard7TimeStamps  = CAENTimeStamps[7];
      fCAENBoard8TimeStamps  = CAENTimeStamps[8];
      fCAENBoard9TimeStamps  = CAENTimeStamps[9];
      fCAENBoard24TimeStamps = CAENTimeStamps[24];

      for (size_t j = 0; j < NumberTdcBlocks; ++j) {

        //std::vector<TDCFragment::TdcEventData> const& tdcEvents = Collection.tdcBlocks[j];
        double const& timestamp = Collection.tdcBlockTimeStamps[j];

        std::cout << "    TDC block: " << j << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
        //std::cout << "      TDC events: " << tdcEvents.size() << std::endl;

        fTDCTimeStamps.push_back(timestamp);
      }

      if (fNumberTPCReadouts > 0) {
        if (i < (Collections.size() - 1)) {
          //fIntervalsDeltaT = (Collections[i+1].interval.first - Collections[i].interval.second) / 1000.0;
          fIntervalsDeltaT = (Collections[i+1].interval.first - Collections[i].interval.second);
          std::cout << "fIntervalsDeltaT [usec]: " << fIntervalsDeltaT << std::endl;
          fIntervalsDeltaTHistogram->Fill(fIntervalsDeltaT / 1000.0);
          fIntervalsDeltaTZHistogram->Fill(fIntervalsDeltaT / 1000.0);
        }
      }

      fNumberTPCReadoutsHistogram->Fill(fNumberTPCReadouts);

      fEventBuilderTree->Fill();
    }

    std::vector< rdu::DataBlock > DataBlocks;
    DataBlocks = fEventBuilderAlg.GetDataBlocks(fLariatFragment);
    std::vector< std::pair< double, double > > CAENBoard0Intervals;
    //CAENBoard0Intervals = fEventBuilderAlg.CreateIntervals(DataBlocks, 0, fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow, fV1740AcquisitionWindow);
    CAENBoard0Intervals = fEventBuilderAlg.CreateIntervals(DataBlocks, 0, 0.128, 0.128, fV1740AcquisitionWindow);
    CAENBoard0Intervals = fEventBuilderAlg.IntervalsSelfMerge(CAENBoard0Intervals);

    //fTPCIntervalsDeltaT = -999999999999999999;
    fTPCIntervalsDeltaT = -1;

    if (CAENBoard0Intervals.size() > 1) {
      for (size_t i = 0; i < (CAENBoard0Intervals.size() - 1); ++i) {
        fTPCIntervalsDeltaT = -1;
        //fTPCIntervalsDeltaT = (CAENBoard0Intervals.at(i+1).first - CAENBoard0Intervals.at(i).second) / 1000.0;
        fTPCIntervalsDeltaT = (CAENBoard0Intervals.at(i+1).first - CAENBoard0Intervals.at(i).second);
        std::cout << "fTPCIntervalsDeltaT [usec]: " << fTPCIntervalsDeltaT << std::endl;
        fTPCIntervalsDeltaTHistogram->Fill(fTPCIntervalsDeltaT / 1000.0);
        fTPCIntervalsDeltaTZHistogram->Fill(fTPCIntervalsDeltaT / 1000.0);
        fTPCTree->Fill();
      }
    }

    return;
  }

  //-----------------------------------------------------------------------
  size_t EventBuilderCheck::CastToSizeT(std::string const& String) 
  {
    size_t SizeT;

    if (!String.empty()) {
        SizeT = static_cast <size_t> (std::stoi(String));
    }
    else {
      SizeT = 0;
    }

    return SizeT;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(EventBuilderCheck)

} // namespace EventBuilderCheck

#endif // EventBuilderCheck_Module
