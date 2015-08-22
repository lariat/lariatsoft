//////////////////////////////////////////////////////////////
// Name:      SlicerCzech_module.cc
// Date:      22 August 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef SlicerCheck_Module
#define SlicerCheck_Module

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
//#include "RawDataUtilities/ClockCorrectionAlg.h"
#include "RawDataUtilities/FragmentUtility.h"
#include "RawDataUtilities/SlicerAlg.h"
#include "Utilities/DatabaseUtilityT1034.h"

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace SlicerCheck {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class SlicerCheck : public art::EDAnalyzer 
  {

   public:

    // Standard constructor and destructor for an ART module.
    explicit SlicerCheck(fhicl::ParameterSet const& pset);
    virtual ~SlicerCheck();

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

    // generate a vector of intervals from a vector of data blocks
    std::vector< std::pair< double, double> > CreateIntervals(std::vector< rdu::DataBlock > const& DataBlocks,
                                                              unsigned int                  const& DeviceID,
                                                              double                        const& PreTriggerWindow,
                                                              double                        const& PostTriggerWindow);

    // merge overlapping intervals in a vector of intervals
    std::vector< std::pair< double, double > > IntervalsSelfMerge(std::vector< std::pair< double, double > > const& Intervals);

    // merge overlapping intervals between two vectors of intervals
    std::vector< std::pair< double, double > > MergeIntervals(std::vector< std::pair< double, double > > const& IntervalsA,
                                                              std::vector< std::pair< double, double > > const& IntervalsB);

    std::string fRawFragmentLabel;    ///< label for module producing artdaq fragments
    std::string fRawFragmentInstance; ///< instance label for artdaq fragments

    // run number
    int fRun;

    // DatabaseUtilityT1034 service handle
    art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

    // clock correction algorithm
    //rdu::ClockCorrectionAlg fClockCorrectionAlg;

    // slicer algorithm
    rdu::SlicerAlg fSlicerAlg;

    // vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

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

    // pre-/post-trigger windows
    double fV1740PreTriggerWindow;
    double fV1740PostTriggerWindow;
    double fV1740BPreTriggerWindow;
    double fV1740BPostTriggerWindow;
    double fV1751PreTriggerWindow;
    double fV1751PostTriggerWindow;
    double fTDCPreTriggerWindow;
    double fTDCPostTriggerWindow;

  }; // class SlicerCheck


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  SlicerCheck::SlicerCheck(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    //, fClockCorrectionAlg(pset.get<fhicl::ParameterSet>("ClockCorrectionAlg"))
    , fSlicerAlg(pset.get<fhicl::ParameterSet>("SlicerAlg"))
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  SlicerCheck::~SlicerCheck() 
  {}

  //-----------------------------------------------------------------------
  void SlicerCheck::beginJob()
  {
    fSlicerAlg.hello_world();

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

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerCheck::beginRun(const art::Run& run)
  {
    fRun = run.run();
    fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams, fRun);

    fV1495DelayTicks       = static_cast <size_t> (std::stoi(fConfigValues["v1495_config_v1495_delay_ticks"]));
    fV1740PostPercent      = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_postpercent"]));
    fV1740BPostPercent     = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_postpercent"]));
    fV1751PostPercent      = static_cast <size_t> (std::stoi(fConfigValues["v1751_config_caen_postpercent"]));
    fV1740RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_recordlength"]));
    fV1740BRecordLength    = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_recordlength"]));
    fV1751RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1751_config_caen_recordlength"]));
    fV1740SampleReduction  = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_v1740_samplereduction"]));
    fV1740BSampleReduction = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_v1740_samplereduction"]));
    fTDCPipelineDelay      = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_pipelinedelay"]));
    fTDCGateWidth          = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_gatewidth"]));

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

    // Pre-/post-trigger windows
    if (fV1740PreTriggerWindow   < 0) fV1740PreTriggerWindow   = fV1740ReadoutWindow;
    if (fV1740PostTriggerWindow  < 0) fV1740PostTriggerWindow  = fV1740ReadoutWindow;
    if (fV1740BPreTriggerWindow  < 0) fV1740BPreTriggerWindow  = fV1740BReadoutWindow;
    if (fV1740BPostTriggerWindow < 0) fV1740BPostTriggerWindow = fV1740BReadoutWindow;
    if (fV1751PreTriggerWindow   < 0) fV1751PreTriggerWindow   = 0.64;
    if (fV1751PostTriggerWindow  < 0) fV1751PostTriggerWindow  = fV1751ReadoutWindow + 0.64;
    if (fTDCPreTriggerWindow     < 0) fTDCPreTriggerWindow     = fTDCReadoutWindow;
    if (fTDCPostTriggerWindow    < 0) fTDCPostTriggerWindow    = fTDCReadoutWindow;

    std::cout << "//////////////////////////////////////////////"        << std::endl;
    std::cout << "V1495DelayTicks:         " << fV1495DelayTicks         << std::endl;
    std::cout << "V1495Delay:              " << fV1495Delay              << std::endl;
    std::cout << "V1740PreTriggerWindow:   " << fV1740PreTriggerWindow   << std::endl;
    std::cout << "V1740PostTriggerWindow:  " << fV1740PostTriggerWindow  << std::endl;
    std::cout << "V1740BPreTriggerWindow:  " << fV1740BPreTriggerWindow  << std::endl;
    std::cout << "V1740BPostTriggerWindow: " << fV1740BPostTriggerWindow << std::endl;
    std::cout << "V1751PreTriggerWindow:   " << fV1751PreTriggerWindow   << std::endl;
    std::cout << "V1751PostTriggerWindow:  " << fV1751PostTriggerWindow  << std::endl;
    std::cout << "TDCPreTriggerWindow:     " << fTDCPreTriggerWindow     << std::endl;
    std::cout << "TDCPostTriggerWindow:    " << fTDCPostTriggerWindow    << std::endl;
    std::cout << "//////////////////////////////////////////////"        << std::endl;

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerCheck::beginSubRun(const art::SubRun& subrun)
  {}

  //-----------------------------------------------------------------------
  void SlicerCheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSlicerAlg.reconfigure(pset.get<fhicl::ParameterSet>("SlicerAlg"));
    fRawFragmentLabel    = pset.get< std::string >("RawFragmentLabel",    "daq"  );
    fRawFragmentInstance = pset.get< std::string >("RawFragmentInstance", "SPILL");

    fV1740PreTriggerWindow   = pset.get< double >("V1740PreTriggerWindow",   -1);
    fV1740PostTriggerWindow  = pset.get< double >("V1740PostTriggerWindow",  -1);
    fV1740BPreTriggerWindow  = pset.get< double >("V1740BPreTriggerWindow",  -1);
    fV1740BPostTriggerWindow = pset.get< double >("V1740BPostTriggerWindow", -1);
    fV1751PreTriggerWindow   = pset.get< double >("V1751PreTriggerWindow",   -1);
    fV1751PostTriggerWindow  = pset.get< double >("V1751PostTriggerWindow",  -1);
    fTDCPreTriggerWindow     = pset.get< double >("TDCPreTriggerWindow",     -1);
    fTDCPostTriggerWindow    = pset.get< double >("TDCPostTriggerWindow",    -1);

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerCheck::analyze(const art::Event& event) 
  {
    // make the utility to access the fragments from the event record
    rdu::FragmentUtility fragUtil(event, fRawFragmentLabel, fRawFragmentInstance);

    // configure the slicer algorithm
    fSlicerAlg.Configure(fV1740PreTriggerWindow,
                         fV1740PostTriggerWindow,
                         fV1740BPreTriggerWindow,
                         fV1740BPostTriggerWindow,
                         fV1751PreTriggerWindow,
                         fV1751PostTriggerWindow,
                         fTDCPreTriggerWindow,
                         fTDCPostTriggerWindow);

    // group data blocks into collections
    std::vector< rdu::DataBlockCollection > Collections;
    Collections = fSlicerAlg.Slice(&fragUtil.DAQFragment());

    for (size_t i = 0; i < Collections.size(); ++i) {
      rdu::DataBlockCollection const& Collection = Collections[i];

      size_t const& NumberCaenBlocks = Collection.caenBlocks.size();
      size_t const& NumberTdcBlocks = Collection.tdcBlocks.size();

      std::cout << "Collection: " << i << std::endl;
      std::cout << "  Number of CAEN data blocks: " << NumberCaenBlocks << std::endl;
      std::cout << "  Number of TDC data blocks:  " << NumberTdcBlocks << std::endl;

      for (size_t j = 0; j < NumberCaenBlocks; ++j) {

        std::pair< double, const CAENFragment * > caenBlock = Collection.caenBlocks[j];
        double const& timestamp = caenBlock.first;
        const CAENFragment * caenFrag = caenBlock.second;
        unsigned int boardId = caenFrag->header.boardId;

        std::cout << "    CAEN block: " << j << std::endl;
        std::cout << "      Board ID: " << boardId << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
      }

      for (size_t j = 0; j < NumberTdcBlocks; ++j) {

        std::pair< double, const std::vector<TDCFragment::TdcEventData> * > tdcBlock = Collection.tdcBlocks[j];
        double const& timestamp = tdcBlock.first;
        //const std::vector<TDCFragment::TdcEventData> * tdcEvents = tdcBlock.second;

        std::cout << "    TDC block: " << j << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
        //std::cout << "      TDC events: " << tdcEvents->size() << std::endl;
      }
    }

    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(SlicerCheck)

} // namespace SlicerCheck

#endif // SlicerCheck_Module
