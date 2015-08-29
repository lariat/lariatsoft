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
#include "RawDataUtilities/FragmentUtility.h"
#include "RawDataUtilities/SlicerAlg.h"
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

    std::string fRawFragmentLabel;    ///< label for module producing artdaq fragments
    std::string fRawFragmentInstance; ///< instance label for artdaq fragments

    // run number
    int fRun;

    // DatabaseUtilityT1034 service handle
    art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

    // slicer algorithm
    rdu::SlicerAlg fSlicerAlg;

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

    // pre-/post-trigger windows
    double fV1740PreAcquisitionWindow;
    double fV1740AcquisitionWindow;
    double fV1740BPreAcquisitionWindow;
    double fV1740BAcquisitionWindow;
    double fV1751PreAcquisitionWindow;
    double fV1751AcquisitionWindow;
    double fTDCPreAcquisitionWindow;
    double fTDCAcquisitionWindow;

    // pointer to histograms
    TH1D * fIntervalsDeltaTHistogram;
    TH1D * fIntervalsDeltaTZHistogram;
    TH1D * fNumberTPCReadoutsHistogram;
    TH1D * fTPCIntervalsDeltaTHistogram;
    TH1D * fTPCIntervalsDeltaTZHistogram;

    // pointer to n-tuples
    TTree * fSlicerTree;
    TTree * fTPCTree;

    // variables that will go into the n-tuples
    int fSubRun;
    double fIntervalsDeltaT;
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

  }; // class SlicerCheck


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  SlicerCheck::SlicerCheck(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
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

    art::ServiceHandle<art::TFileService> tfs;

    fSlicerTree = tfs->make<TTree>("SlicerTree", "SlicerTree");

    fSlicerTree->Branch("Run",               &fRun,               "Run/I");
    fSlicerTree->Branch("SubRun",            &fSubRun,            "SubRun/I");
    fSlicerTree->Branch("IntervalsDeltaT",   &fIntervalsDeltaT,   "IntervalsDeltaT/D");
    fSlicerTree->Branch("NumberTPCReadouts", &fNumberTPCReadouts, "NumberTPCReadouts/I");
    fSlicerTree->Branch("NumberCAENBlocks",  &fNumberCAENBlocks,  "NumberCAENBlocks/I");
    fSlicerTree->Branch("NumberTDCBlocks",   &fNumberTDCBlocks,   "NumberTDCBlocks/I");
    fSlicerTree->Branch("NumberCAENBoard0Blocks",  &fNumberCAENBoard0Blocks,  "NumberCAENBoard0Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard1Blocks",  &fNumberCAENBoard1Blocks,  "NumberCAENBoard1Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard2Blocks",  &fNumberCAENBoard2Blocks,  "NumberCAENBoard2Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard3Blocks",  &fNumberCAENBoard3Blocks,  "NumberCAENBoard3Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard4Blocks",  &fNumberCAENBoard4Blocks,  "NumberCAENBoard4Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard5Blocks",  &fNumberCAENBoard5Blocks,  "NumberCAENBoard5Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard6Blocks",  &fNumberCAENBoard6Blocks,  "NumberCAENBoard6Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard7Blocks",  &fNumberCAENBoard7Blocks,  "NumberCAENBoard7Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard8Blocks",  &fNumberCAENBoard8Blocks,  "NumberCAENBoard8Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard9Blocks",  &fNumberCAENBoard9Blocks,  "NumberCAENBoard9Blocks/I");
    fSlicerTree->Branch("NumberCAENBoard24Blocks", &fNumberCAENBoard24Blocks, "NumberCAENBoard24Blocks/I");

    fTPCTree = tfs->make<TTree>("TPCTree", "TPCTree");
    fTPCTree->Branch("Run",                &fRun,                "Run/I");
    fTPCTree->Branch("SubRun",             &fSubRun,             "SubRun/I");
    fTPCTree->Branch("TPCIntervalsDeltaT", &fTPCIntervalsDeltaT, "TPCIntervalsDeltaT/D");

    fIntervalsDeltaTHistogram      = tfs->make<TH1D>("IntervalsDeltaT",     ";#Delta t [ms];Entries per ms", 5000, -0.5, 4999.5);
    fIntervalsDeltaTZHistogram     = tfs->make<TH1D>("IntervalsDeltaTZ",    ";#Delta t [ms];Entries per 0.001 ms", 1000, 0, 1);
    fNumberTPCReadoutsHistogram    = tfs->make<TH1D>("NumberTPCReadouts",   ";Number of TPC readouts;Entries", 10, -0.5, 9.5);
    fTPCIntervalsDeltaTHistogram   = tfs->make<TH1D>("TPCIntervalsDeltaT",  ";#Delta t [ms];Entries per ms", 5000, -0.5, 4999.5);
    fTPCIntervalsDeltaTZHistogram  = tfs->make<TH1D>("TPCIntervalsDeltaTZ", ";#Delta t [ms];Entries per 0.001 ms", 1000, 0, 1);

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerCheck::beginRun(const art::Run& run)
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

    // Pre-acquisition and acquisition windows
    if (fV1740PreAcquisitionWindow  < 0) fV1740PreAcquisitionWindow  = fV1740ReadoutWindow;
    if (fV1740AcquisitionWindow     < 0) fV1740AcquisitionWindow     = fV1740ReadoutWindow;
    if (fV1740BPreAcquisitionWindow < 0) fV1740BPreAcquisitionWindow = fV1740BReadoutWindow;
    if (fV1740BAcquisitionWindow    < 0) fV1740BAcquisitionWindow    = fV1740BReadoutWindow;
    if (fV1751PreAcquisitionWindow  < 0) fV1751PreAcquisitionWindow  = 0.64;
    if (fV1751AcquisitionWindow     < 0) fV1751AcquisitionWindow     = fV1751ReadoutWindow + 0.64;
    if (fTDCPreAcquisitionWindow    < 0) fTDCPreAcquisitionWindow    = fTDCReadoutWindow;
    if (fTDCAcquisitionWindow       < 0) fTDCAcquisitionWindow       = fTDCReadoutWindow;

    std::cout << "//////////////////////////////////////////////"              << std::endl;
    std::cout << "V1495DelayTicks:            " << fV1495DelayTicks            << std::endl;
    std::cout << "V1495Delay:                 " << fV1495Delay                 << std::endl;
    std::cout << "V1740PreAcquisitionWindow:  " << fV1740PreAcquisitionWindow  << std::endl;
    std::cout << "V1740AcquisitionWindow:     " << fV1740AcquisitionWindow     << std::endl;
    std::cout << "V1740BPreAcquisitionWindow: " << fV1740BPreAcquisitionWindow << std::endl;
    std::cout << "V1740BAcquisitionWindow:    " << fV1740BAcquisitionWindow    << std::endl;
    std::cout << "V1751PreAcquisitionWindow:  " << fV1751PreAcquisitionWindow  << std::endl;
    std::cout << "V1751AcquisitionWindow:     " << fV1751AcquisitionWindow     << std::endl;
    std::cout << "TDCPreAcquisitionWindow:    " << fTDCPreAcquisitionWindow    << std::endl;
    std::cout << "TDCAcquisitionWindow:       " << fTDCAcquisitionWindow       << std::endl;
    std::cout << "//////////////////////////////////////////////"              << std::endl;

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

    fV1740PreAcquisitionWindow  = pset.get< double >("V1740PreAcquisitionWindow",  -1);
    fV1740AcquisitionWindow     = pset.get< double >("V1740AcquisitionWindow",     -1);
    fV1740BPreAcquisitionWindow = pset.get< double >("V1740BPreAcquisitionWindow", -1);
    fV1740BAcquisitionWindow    = pset.get< double >("V1740BAcquisitionWindow",    -1);
    fV1751PreAcquisitionWindow  = pset.get< double >("V1751PreAcquisitionWindow",  -1);
    fV1751AcquisitionWindow     = pset.get< double >("V1751AcquisitionWindow",     -1);
    fTDCPreAcquisitionWindow    = pset.get< double >("TDCPreAcquisitionWindow",    -1);
    fTDCAcquisitionWindow       = pset.get< double >("TDCAcquisitionWindow",       -1);

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerCheck::analyze(const art::Event& event) 
  {
    fSubRun = event.subRun();

    // make the utility to access the fragments from the event record
    rdu::FragmentUtility fragUtil(event, fRawFragmentLabel, fRawFragmentInstance);

    // configure the slicer algorithm
    fSlicerAlg.Configure(fV1740PreAcquisitionWindow,
                         fV1740AcquisitionWindow,
                         fV1740BPreAcquisitionWindow,
                         fV1740BAcquisitionWindow,
                         fV1751PreAcquisitionWindow,
                         fV1751AcquisitionWindow,
                         fTDCPreAcquisitionWindow,
                         fTDCAcquisitionWindow);

    // group data blocks into collections
    std::vector< rdu::DataBlockCollection > Collections;
    Collections = fSlicerAlg.Slice(&fragUtil.DAQFragment());

    for (size_t i = 0; i < Collections.size(); ++i) {
      rdu::DataBlockCollection const& Collection = Collections[i];

      size_t const& NumberCaenBlocks = Collection.caenBlocks.size();
      size_t const& NumberTdcBlocks  = Collection.tdcBlocks.size();

      fNumberCAENBlocks = NumberCaenBlocks;
      fNumberTDCBlocks  = NumberTdcBlocks;

      fIntervalsDeltaT = -1;
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

        if (boardId == 0) {
          fNumberTPCReadouts += 1;
        }

        NumberCAENBlocksArray[boardId] += 1;

      }

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

      for (size_t j = 0; j < NumberTdcBlocks; ++j) {

        std::pair< double, const std::vector<TDCFragment::TdcEventData> * > tdcBlock = Collection.tdcBlocks[j];
        double const& timestamp = tdcBlock.first;
        //const std::vector<TDCFragment::TdcEventData> * tdcEvents = tdcBlock.second;

        std::cout << "    TDC block: " << j << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
        //std::cout << "      TDC events: " << tdcEvents->size() << std::endl;
      }

      if (fNumberTPCReadouts > 0) {
        if (i < (Collections.size() - 1)) {
          fIntervalsDeltaT = (Collections[i+1].interval.second - Collections[i].interval.first) / 1000.0;
          std::cout << "fIntervalsDeltaT [ms]: " << fIntervalsDeltaT << std::endl;
          fIntervalsDeltaTHistogram->Fill(fIntervalsDeltaT);
          fIntervalsDeltaTZHistogram->Fill(fIntervalsDeltaT);
        }
      }

      fNumberTPCReadoutsHistogram->Fill(fNumberTPCReadouts);

      fSlicerTree->Fill();
    }

    std::vector< rdu::DataBlock > DataBlocks;
    DataBlocks = fSlicerAlg.GetDataBlocks(&fragUtil.DAQFragment());
    std::vector< std::pair< double, double > > CAENBoard0Intervals;
    CAENBoard0Intervals = fSlicerAlg.CreateIntervals(DataBlocks, 0, fV1740PreAcquisitionWindow, fV1740AcquisitionWindow);
    CAENBoard0Intervals = fSlicerAlg.IntervalsSelfMerge(CAENBoard0Intervals);

    fTPCIntervalsDeltaT = -1;

    if (CAENBoard0Intervals.size() > 1) {
      for (size_t i = 0; i < (CAENBoard0Intervals.size() - 1); ++i) {
        fTPCIntervalsDeltaT = -1;
        fTPCIntervalsDeltaT = (CAENBoard0Intervals.at(i+1).second - CAENBoard0Intervals.at(i).first) / 1000.0;
        std::cout << "fTPCIntervalsDeltaT [ms]: " << fTPCIntervalsDeltaT << std::endl;
        fTPCIntervalsDeltaTHistogram->Fill(fTPCIntervalsDeltaT);
        fTPCIntervalsDeltaTZHistogram->Fill(fTPCIntervalsDeltaT);
        fTPCTree->Fill();
      }
    }

    return;
  }

  //-----------------------------------------------------------------------
  size_t SlicerCheck::CastToSizeT(std::string const& String) 
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
  DEFINE_ART_MODULE(SlicerCheck)

} // namespace SlicerCheck

#endif // SlicerCheck_Module
