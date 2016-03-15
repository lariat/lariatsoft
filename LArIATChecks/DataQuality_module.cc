//////////////////////////////////////////////////////////////
// Name:      DataQuality_module.cc
// Date:      4 February 2016
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef DataQuality_Module
#define DataQuality_Module

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

// LArSoft includes
#include "lardata/RawData/AuxDetDigit.h"
#include "lardata/RawData/OpDetPulse.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/TriggerData.h"

// LArIATSoft includes
#include "LArIATRecoAlg/TOFBuilderAlg.h"
#include "RawDataUtilities/FragmentToDigitAlg.h"
#include "RawDataUtilities/EventBuilderAlg.h"
#include "RawDataUtilities/SpillWrapper.h"
#include "Utilities/DatabaseUtilityT1034.h"

// ROOT includes
#include "TH1.h"
#include "TTree.h"

// C++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

// enumerated types for hardware
enum {
  V1740_N_BOARDS = 8,
  V1740_N_CHANNELS = 64,
  V1740_N_SAMPLES = 192 * 16,
  V1740B_N_BOARDS = 1,
  V1740B_N_CHANNELS = 64,
  V1740B_N_SAMPLES = 192 * 16,
  V1751_N_BOARDS = 2,
  V1751_N_CHANNELS = 8,
  V1751_N_SAMPLES = 1792 * 16,
  WUT_MAX_HITS = 128,
};

namespace DataQuality {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class DataQuality : public art::EDAnalyzer 
  {

   public:

    // Standard constructor and destructor for an ART module.
    explicit DataQuality(fhicl::ParameterSet const& pset);
    virtual ~DataQuality();

    // This method is called once, at the start of the job.
    void beginJob();

    // This method is called once, at the start of each run.
    void beginRun(const art::Run& run);

    // This method is called once, at the start of each sub-run
    void beginSubRun(const art::SubRun& subrun);

    // This method reads in any parameters from the .fcl files.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event.
    void analyze(const art::Event& event);

   private:

    std::string fRawFragmentLabel;    ///< label for module producing artdaq fragments
    std::string fRawFragmentInstance; ///< instance label for artdaq fragments

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

    // fragment-to-digit algorithm
    FragmentToDigitAlg fFragmentToDigitAlg;

    // TOF builder algorithm
    TOFBuilderAlg fTOFBuilderAlg;

    // timestamp from SpillTrailer
    std::uint64_t fTimestamp;

    // vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

    // flag for PEDESTALON gate
    bool fPedestalOn;

    size_t castToSizeT_(std::string const& String);

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

    TH1I * fUSTOFHitsHistogram;
    TH1I * fDSTOFHitsHistogram;
    TH1D * fTOFHistogram;

    std::vector< std::vector< TH1I * > > fCAENPedestalHistograms;
    std::vector< std::vector< TH1I * > > fCAENADCHistograms;
    std::vector< std::vector< TH1I * > > fCAENMinADCHistograms;
    std::vector< std::vector< TH1I * > > fCAENMaxADCHistograms;

    std::vector< TH1D * > fCAENPedestalTimeStampHistograms;
    std::vector< TH1D * > fCAENTimeStampHistograms;
    TH1D * fMWPCTDCTimeStampHistograms;
    TH1D * fWUTTimeStampHistograms;

    TH1I * fTDCTimeBitMismatchHistogram;

    // pointer to n-tuples
    TTree * fEventBuilderTree;    ///< Tree holding variables on the performance of the event builder
    TTree * fTPCTree;             ///< Tree holding variables on the performance of the event builder on TPC events
    TTree * fEventRecord;         ///< Tree holding some data from art::Event
    TTree * fSpillTrailerTree;    ///< Tree holding the data from the SpillTrailer fragments
    TTree * fCaenV1740DataTree;   ///< Tree holding the data from the CAEN V1740 fragments
    TTree * fCaenV1740BDataTree;  ///< Tree holding the data from the CAEN V1740B fragments
    TTree * fCaenV1751DataTree;   ///< Tree holding the data from the CAEN V1751 fragments
    TTree * fMwpcTdcDataTree;     ///< Tree holding the data from the MWPC TDC fragments
    TTree * fWutDataTree;         ///< Tree holding the data from the Wave Union TDC fragments

    // variables that will go into the n-tuples
    // variables that will go into fEventBuilderTree and/or fTPCTree
    int                 fRun;
    int                 fSubRun;
    int                 fEvent;
    //int                 fEventCounter;
    double              fIntervalsDeltaT;
    double              fIntervalsDeltaTBeginToBegin;
    double              fIntervalLength;
    double              fIntervalStart;
    double              fIntervalStop;
    int                 fNumberTPCReadouts;
    int                 fNumberCAENBlocks;
    int                 fNumberCAENBoard0Blocks;
    int                 fNumberCAENBoard1Blocks;
    int                 fNumberCAENBoard2Blocks;
    int                 fNumberCAENBoard3Blocks;
    int                 fNumberCAENBoard4Blocks;
    int                 fNumberCAENBoard5Blocks;
    int                 fNumberCAENBoard6Blocks;
    int                 fNumberCAENBoard7Blocks;
    int                 fNumberCAENBoard8Blocks;
    int                 fNumberCAENBoard9Blocks;
    int                 fNumberCAENBoard24Blocks;
    int                 fNumberTDCBlocks;
    double              fTPCIntervalsDeltaT;
    double              fTPCIntervalsDeltaTBeginToBegin;
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

    uint32_t            fTimeStampLow;
    uint32_t            fTimeStampHigh;
    uint32_t            fSpillTrailerRunNumber;
    uint32_t            fSpillTrailerSpillNumber;
    uint32_t            fSpillTrailerTimeStamp;

    // variables that will go into the CAEN trees
    uint32_t                               fCaenFragment;
    uint32_t                               fCaenBoardId;
    uint32_t                               fCaenEventCounter;    // The data block collection index will be used in place of the CAEN Event Counter
    uint32_t                               fCaenTriggerTimeTag;  // Each count in the CAEN Trigger Time Tag is 8 ns
    uint32_t                               fCaenNumberSamples;   // Number of samples
    std::vector< std::vector< uint16_t > > fCaenV1740Waveform;   // 128 ns per sample, 3072 samples per trigger
    std::vector< std::vector< uint16_t > > fCaenV1740BWaveform;  // 128 ns per sample, 3072 samples per trigger
    std::vector< std::vector< uint16_t > > fCaenV1751Waveform;   // 1 ns per sample

    // variables that will go into fMwpcTdcDataTree
    uint32_t              fMwpcEventCounter;         // The data block collection index will be used here
    uint32_t              fMwpcTriggerCounter;
    uint16_t              fMwpcControllerTimeStamp;
    uint32_t              fMwpcTdcTimeStamp;         // Each TDC timestamp count is 1/106.208e6 seconds
    uint32_t              fMwpcNumberHits;
    std::vector<uint16_t> fMwpcTdcNumber;
    std::vector<uint16_t> fMwpcHitChannel;
    std::vector<uint16_t> fMwpcHitTimeBin;           // Each time tick is 1.177 ns

    // variables that will go into fWutDataTree
    uint32_t              fWutTimeHeader;  // Each count in the time header is 16 us
    uint32_t              fWutNumberHits;
    std::vector<uint16_t> fWutHitChannel;
    std::vector<uint32_t> fWutHitTimeBin;  // Each time tick is 15.625 ps

  }; // class DataQuality


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  DataQuality::DataQuality(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fNumberFragmentsPerSpill(pset.get<std::size_t>("NumberFragmentsPerSpill", 4))
    , fSpillWrapper(new rdu::SpillWrapper(fNumberFragmentsPerSpill))
    , fLariatFragment(nullptr)
    , fEventBuilderAlg(pset.get<fhicl::ParameterSet>("EventBuilderAlg"))
    , fFragmentToDigitAlg(pset.get<fhicl::ParameterSet>("FragmentToDigitAlg"))
    , fTOFBuilderAlg(pset)
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  DataQuality::~DataQuality() 
  {}

  //-----------------------------------------------------------------------
  void DataQuality::beginJob()
  {
    fEventBuilderAlg.hello_world();

    // parameters to fetch from the lariat_xml_database
    // table of the lariat_prd database
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

    // resize vectors for TTrees
    fCaenV1740Waveform.resize(V1740_N_CHANNELS);
    fCaenV1740BWaveform.resize(V1740B_N_CHANNELS);
    fCaenV1751Waveform.resize(V1751_N_CHANNELS);

    // resize vectors for pedestal TH1 objects
    //  0 ...  7  CAEN V1740
    //  8 ... 15  CAEN V1751
    // 15 ... 23
    // 24 ... 31  CAEN V1740
    fCAENPedestalHistograms.resize(32);
    fCAENADCHistograms.resize(32);
    fCAENMinADCHistograms.resize(32);
    fCAENMaxADCHistograms.resize(32);

    fCAENPedestalTimeStampHistograms.resize(32);
    fCAENTimeStampHistograms.resize(32);

    // TFile service
    art::ServiceHandle<art::TFileService> tfs;

    // create sub-directory for TOF
    art::TFileDirectory timestampDir = tfs->mkdir("timestamps");

    // create sub-directory for TOF
    art::TFileDirectory tofDir = tfs->mkdir("tof");

    // create sub-directories for pedestal and ADC histograms
    art::TFileDirectory pedestalDir = tfs->mkdir("pedestal");
    art::TFileDirectory adcDir      = tfs->mkdir("adc");
    art::TFileDirectory minADCDir   = tfs->mkdir("min_adc");
    art::TFileDirectory maxADCDir   = tfs->mkdir("max_adc");

    // create TH1 objects
    fIntervalsDeltaTHistogram     = tfs->make<TH1D>("IntervalsDeltaT",     ";#Delta t [ms];Entries per ms",       10000, -0.5, 9999.5);
    fIntervalsDeltaTZHistogram    = tfs->make<TH1D>("IntervalsDeltaTZ",    ";#Delta t [ms];Entries per 0.001 ms",  1000,    0,    1);
    fNumberTPCReadoutsHistogram   = tfs->make<TH1D>("NumberTPCReadouts",   ";Number of TPC readouts;Entries",        10, -0.5,    9.5);
    fTPCIntervalsDeltaTHistogram  = tfs->make<TH1D>("TPCIntervalsDeltaT",  ";#Delta t [ms];Entries per ms",       10000, -0.5, 9999.5);
    fTPCIntervalsDeltaTZHistogram = tfs->make<TH1D>("TPCIntervalsDeltaTZ", ";#Delta t [ms];Entries per 0.001 ms",  1000,    0,    1);

    // TH1 objects for timestamps
    fMWPCTDCTimeStampHistograms = timestampDir.make<TH1D>("mwpc_tdc_timestamps", ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
    fWUTTimeStampHistograms     = timestampDir.make<TH1D>("wut_timestamps", ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);

    for (size_t i = 0; i < V1740_N_BOARDS; ++i) {
      std::string th1Title = "caen_board_" + std::to_string(i);
      fCAENPedestalTimeStampHistograms[i] = timestampDir.make<TH1D>((th1Title + "_pedestal_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
      fCAENTimeStampHistograms[i]         = timestampDir.make<TH1D>((th1Title + "_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
    }

    for (size_t i = 0; i < V1751_N_BOARDS; ++i) {
      size_t offset = 8;
      std::string th1Title = "caen_board_" + std::to_string(i + offset);
      fCAENPedestalTimeStampHistograms[i + offset] = timestampDir.make<TH1D>((th1Title + "_pedestal_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
      fCAENTimeStampHistograms[i + offset]         = timestampDir.make<TH1D>((th1Title + "_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
    }

    for (size_t i = 0; i < V1740B_N_BOARDS; ++i) {
      size_t offset = 24;
      std::string th1Title = "caen_board_" + std::to_string(i + offset);
      fCAENPedestalTimeStampHistograms[i + offset] = timestampDir.make<TH1D>((th1Title + "_pedestal_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
      fCAENTimeStampHistograms[i + offset]         = timestampDir.make<TH1D>((th1Title + "_timestamps").c_str(), ";Timestamp [s];Entries per 0.2 s", 300, 0, 60);
    }

    // create TH1 objects in pedestal and ADC sub-directories
    for (size_t i = 0; i < V1740_N_BOARDS; ++i) {
      fCAENPedestalHistograms[i].resize(V1740_N_CHANNELS);
      fCAENADCHistograms[i].resize(V1740_N_CHANNELS);
      fCAENMinADCHistograms[i].resize(V1740_N_CHANNELS);
      fCAENMaxADCHistograms[i].resize(V1740_N_CHANNELS);
      for (size_t j = 0; j < V1740_N_CHANNELS; ++j) {
        std::string th1Title = "caen_board_" + std::to_string(i) + "_channel_" + std::to_string(j);
        fCAENPedestalHistograms[i][j] = pedestalDir.make<TH1I>((th1Title + "_pedestal").c_str(),
                                                               ";ADC;Entries per ADC",
                                                               4096, 0, 4096);
        fCAENADCHistograms[i][j] = adcDir.make<TH1I>((th1Title + "_adc").c_str(),
                                                     ";ADC;Entries per ADC",
                                                     4096, 0, 4096);
        fCAENMinADCHistograms[i][j] = minADCDir.make<TH1I>((th1Title + "_min_adc").c_str(),
                                                           ";ADC;Entries per ADC",
                                                           4096, 0, 4096);
        fCAENMaxADCHistograms[i][j] = maxADCDir.make<TH1I>((th1Title + "_max_adc").c_str(),
                                                           ";ADC;Entries per ADC",
                                                           4096, 0, 4096);
      }
    }

    for (size_t i = 0; i < V1751_N_BOARDS; ++i) {
      size_t offset = 8;
      fCAENPedestalHistograms[i + offset].resize(V1751_N_CHANNELS);
      fCAENADCHistograms[i + offset].resize(V1751_N_CHANNELS);
      fCAENMinADCHistograms[i + offset].resize(V1751_N_CHANNELS);
      fCAENMaxADCHistograms[i + offset].resize(V1751_N_CHANNELS);
      for (size_t j = 0; j < V1751_N_CHANNELS; ++j) {
        std::string th1Title = "caen_board_" + std::to_string(i + offset) + "_channel_" + std::to_string(j);
        fCAENPedestalHistograms[i + offset][j] = pedestalDir.make<TH1I>((th1Title + "_pedestal").c_str(),
                                                                        ";ADC;Entries per ADC",
                                                                        1024, 0, 1024);
        fCAENADCHistograms[i + offset][j] = adcDir.make<TH1I>((th1Title + "_adc").c_str(),
                                                              ";ADC;Entries per ADC",
                                                              1024, 0, 1024);
        fCAENMinADCHistograms[i + offset][j] = minADCDir.make<TH1I>((th1Title + "_min_adc").c_str(),
                                                                    ";ADC;Entries per ADC",
                                                                    1024, 0, 1024);
        fCAENMaxADCHistograms[i + offset][j] = maxADCDir.make<TH1I>((th1Title + "_max_adc").c_str(),
                                                                    ";ADC;Entries per ADC",
                                                                    1024, 0, 1024);
      }
    }

    for (size_t i = 0; i < V1740B_N_BOARDS; ++i) {
      size_t offset = 24;
      fCAENPedestalHistograms[i + offset].resize(V1740B_N_CHANNELS);
      fCAENADCHistograms[i + offset].resize(V1740B_N_CHANNELS);
      fCAENMinADCHistograms[i + offset].resize(V1740B_N_CHANNELS);
      fCAENMaxADCHistograms[i + offset].resize(V1740B_N_CHANNELS);
      for (size_t j = 0; j < V1740B_N_CHANNELS; ++j) {
        std::string th1Title = "caen_board_" + std::to_string(i + offset) + "_channel_" + std::to_string(j);
        fCAENPedestalHistograms[i + offset][j] = pedestalDir.make<TH1I>((th1Title + "_pedestal").c_str(),
                                                                        ";ADC;Entries per ADC",
                                                                        4096, 0, 4096);
        fCAENADCHistograms[i + offset][j] = adcDir.make<TH1I>((th1Title + "_adc").c_str(),
                                                              ";ADC;Entries per ADC",
                                                              4096, 0, 4096);
        fCAENMinADCHistograms[i + offset][j] = minADCDir.make<TH1I>((th1Title + "_min_adc").c_str(),
                                                                    ";ADC;Entries per ADC",
                                                                    4096, 0, 4096);
        fCAENMaxADCHistograms[i + offset][j] = maxADCDir.make<TH1I>((th1Title + "_max_adc").c_str(),
                                                                    ";ADC;Entries per ADC",
                                                                    4096, 0, 4096);
      }
    }

    // TH1 objects for TDC time bit mismatches
    fTDCTimeBitMismatchHistogram = tfs->make<TH1I>("TDCTimeBitMismatch", ";TDC with time bit mismatch;Entries per bin", 16, 0+1, 16+1);

    // TH1 objects for TOF
    fUSTOFHitsHistogram = tofDir.make<TH1I>("USTOFHits", ";Clock tick;Entries per clock tick", V1751_N_SAMPLES, 0, V1751_N_SAMPLES);
    fDSTOFHitsHistogram = tofDir.make<TH1I>("DSTOFHits", ";Clock tick;Entries per clock tick", V1751_N_SAMPLES, 0, V1751_N_SAMPLES);
    fTOFHistogram = tofDir.make<TH1D>("TOF", ";TOF [ns];Entries per ns", 500, 0, 500);

    // create TTree objects
    fEventRecord        = tfs->make<TTree>("artEventRecord",  "artEventRecord");
    fSpillTrailerTree   = tfs->make<TTree>("spillTrailer",    "spillTrailer");
    fEventBuilderTree   = tfs->make<TTree>("EventBuilderTree", "EventBuilderTree");
    fTPCTree            = tfs->make<TTree>("TPCTree",          "TPCTree");
    fCaenV1740DataTree  = tfs->make<TTree>("v1740",            "v1740");
    fCaenV1740BDataTree = tfs->make<TTree>("v1740b",           "v1740b");
    fCaenV1751DataTree  = tfs->make<TTree>("v1751",            "v1751");
    fMwpcTdcDataTree    = tfs->make<TTree>("mwpc",             "mwpc");
    fWutDataTree        = tfs->make<TTree>("wut",              "wut");

    // fEventBuilderTree branches
    fEventBuilderTree->Branch("Run",                         &fRun,                         "Run/I");
    fEventBuilderTree->Branch("SubRun",                      &fSubRun,                      "SubRun/I");
    fEventBuilderTree->Branch("IntervalsDeltaT",             &fIntervalsDeltaT,             "IntervalsDeltaT/D");
    fEventBuilderTree->Branch("IntervalsDeltaTBeginToBegin", &fIntervalsDeltaTBeginToBegin, "IntervalsDeltaTBeginToBegin/D");
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

    // fTPCTree branches
    fTPCTree->Branch("Run",                            &fRun,                            "Run/I");
    fTPCTree->Branch("SubRun",                         &fSubRun,                         "SubRun/I");
    fTPCTree->Branch("TPCIntervalsDeltaT",             &fTPCIntervalsDeltaT,             "TPCIntervalsDeltaT/D");
    fTPCTree->Branch("TPCIntervalsDeltaTBeginToBegin", &fTPCIntervalsDeltaTBeginToBegin, "TPCIntervalsDeltaTBeginToBegin/D");

    // fEventRecord branches
    fEventRecord->Branch("run_number",      &fRun,           "run_number/I");
    fEventRecord->Branch("sub_run_number",  &fSubRun,        "sub_run_number/I");
    fEventRecord->Branch("event_number",    &fEvent,         "event_number/I");
    fEventRecord->Branch("time_stamp_low",  &fTimeStampLow,  "time_stamp_low/i");
    fEventRecord->Branch("time_stamp_high", &fTimeStampHigh, "time_stamp_high/i");

    // fSpillTrailerTree branches
    fSpillTrailerTree->Branch("runNumber",   &fSpillTrailerRunNumber,   "runNumber/i");
    fSpillTrailerTree->Branch("spillNumber", &fSpillTrailerSpillNumber, "spillNumber/i");
    fSpillTrailerTree->Branch("timeStamp",   &fSpillTrailerTimeStamp,   "timeStamp/i");

    // fCaenV1740DataTree branches
    fCaenV1740DataTree->Branch("run",              &fRun,                "run/I");
    fCaenV1740DataTree->Branch("spill",            &fSubRun,             "spill/I");
    fCaenV1740DataTree->Branch("fragment",         &fCaenFragment,       "fragment/i");
    fCaenV1740DataTree->Branch("event_counter",    &fCaenEventCounter,   "event_counter/i");
    fCaenV1740DataTree->Branch("board_id",         &fCaenBoardId,        "board_id/i");
    fCaenV1740DataTree->Branch("trigger_time_tag", &fCaenTriggerTimeTag, "trigger_time_tag/i");
    fCaenV1740DataTree->Branch("number_samples",   &fCaenNumberSamples,  "number_samples/i");

    for (size_t i = 0; i < V1740_N_CHANNELS; ++i) {
      std::string branchName = "channel_" + std::to_string(i);
      fCaenV1740DataTree->Branch(branchName.c_str(), &fCaenV1740Waveform[i]);
    }

    // fCaenV1740BDataTree branches
    fCaenV1740BDataTree->Branch("run",              &fRun,                "run/I");
    fCaenV1740BDataTree->Branch("spill",            &fSubRun,             "spill/I");
    fCaenV1740BDataTree->Branch("fragment",         &fCaenFragment,       "fragment/i");
    fCaenV1740BDataTree->Branch("event_counter",    &fCaenEventCounter,   "event_counter/i");
    fCaenV1740BDataTree->Branch("board_id",         &fCaenBoardId,        "board_id/i");
    fCaenV1740BDataTree->Branch("trigger_time_tag", &fCaenTriggerTimeTag, "trigger_time_tag/i");
    fCaenV1740BDataTree->Branch("number_samples",   &fCaenNumberSamples,  "number_samples/i");

    for (size_t i = 0; i < V1740B_N_CHANNELS; ++i) {
      std::string branchName = "channel_" + std::to_string(i);
      fCaenV1740BDataTree->Branch(branchName.c_str(), &fCaenV1740BWaveform[i]);
    }

    // fCaenV1751DataTree branches
    fCaenV1751DataTree->Branch("run",              &fRun,                "run/I");
    fCaenV1751DataTree->Branch("spill",            &fSubRun,             "spill/I");
    fCaenV1751DataTree->Branch("fragment",         &fCaenFragment,       "fragment/i");
    fCaenV1751DataTree->Branch("event_counter",    &fCaenEventCounter,   "event_counter/i");
    fCaenV1751DataTree->Branch("board_id",         &fCaenBoardId,        "board_id/i");
    fCaenV1751DataTree->Branch("trigger_time_tag", &fCaenTriggerTimeTag, "trigger_time_tag/i");
    fCaenV1751DataTree->Branch("number_samples",   &fCaenNumberSamples,  "number_samples/i");

    for (size_t i = 0; i < V1751_N_CHANNELS; ++i) {
      std::string branchName = "channel_" + std::to_string(i);
      fCaenV1751DataTree->Branch(branchName.c_str(), &fCaenV1751Waveform[i]);
    }

    // fMwpcTdcDataTree branches
    fMwpcTdcDataTree->Branch("run",                   &fRun,                     "run/I");
    fMwpcTdcDataTree->Branch("spill",                 &fSubRun,                  "spill/I");
    fMwpcTdcDataTree->Branch("event_counter",         &fMwpcEventCounter,        "event_counter/i");
    fMwpcTdcDataTree->Branch("trigger_counter",       &fMwpcTriggerCounter,      "trigger_counter/i");
    fMwpcTdcDataTree->Branch("controller_time_stamp", &fMwpcControllerTimeStamp, "controller_time_stamp/s");
    fMwpcTdcDataTree->Branch("tdc_time_stamp",        &fMwpcTdcTimeStamp,        "tdc_time_stamp/i");
    fMwpcTdcDataTree->Branch("number_hits",           &fMwpcNumberHits,          "number_hits/i");
    fMwpcTdcDataTree->Branch("tdc_number",            &fMwpcTdcNumber);
    fMwpcTdcDataTree->Branch("hit_channel",           &fMwpcHitChannel);
    fMwpcTdcDataTree->Branch("hit_time_bin",          &fMwpcHitTimeBin);

    // fWutDataTree branches
    fWutDataTree->Branch("run",          &fRun,           "run/I");
    fWutDataTree->Branch("spill",        &fSubRun,        "spill/I");
    fWutDataTree->Branch("time_header",  &fWutTimeHeader, "time_header/i");
    fWutDataTree->Branch("number_hits",  &fWutNumberHits, "number_hits/i");
    fWutDataTree->Branch("hit_channel",  &fWutHitChannel);
    fWutDataTree->Branch("hit_time_bin", &fWutHitTimeBin);

    return;
  }

  //-----------------------------------------------------------------------
  void DataQuality::beginRun(const art::Run& run)
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

    fV1495DelayTicks       = this->castToSizeT_(fConfigValues["v1495_config_v1495_delay_ticks"]);
    fV1740PostPercent      = this->castToSizeT_(fConfigValues["v1740_config_caen_postpercent"]);
    fV1740BPostPercent     = this->castToSizeT_(fConfigValues["v1740b_config_caen_postpercent"]);
    fV1751PostPercent      = this->castToSizeT_(fConfigValues["v1751_config_caen_postpercent"]);
    fV1740RecordLength     = this->castToSizeT_(fConfigValues["v1740_config_caen_recordlength"]);
    fV1740BRecordLength    = this->castToSizeT_(fConfigValues["v1740b_config_caen_recordlength"]);
    fV1751RecordLength     = this->castToSizeT_(fConfigValues["v1751_config_caen_recordlength"]);
    fV1740SampleReduction  = this->castToSizeT_(fConfigValues["v1740_config_caen_v1740_samplereduction"]);
    fV1740BSampleReduction = this->castToSizeT_(fConfigValues["v1740b_config_caen_v1740_samplereduction"]);
    fTDCPipelineDelay      = this->castToSizeT_(fConfigValues["tdc_config_tdc_pipelinedelay"]);
    fTDCGateWidth          = this->castToSizeT_(fConfigValues["tdc_config_tdc_gatewidth"]);

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
  void DataQuality::beginSubRun(const art::SubRun& subrun)
  {}

  //-----------------------------------------------------------------------
  void DataQuality::reconfigure(fhicl::ParameterSet const& pset)
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
  void DataQuality::analyze(const art::Event& event) 
  {
    //fRun           = event.run();
    fSubRun        = event.subRun();
    fEvent         = event.event();
    fTimeStampLow  = event.time().timeLow();
    fTimeStampHigh = event.time().timeHigh();

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

    // get SpillTrailer
    LariatFragment::SpillTrailer const& spillTrailer = fLariatFragment->spillTrailer;
    fSpillTrailerRunNumber = spillTrailer.runNumber;
    fSpillTrailerSpillNumber = spillTrailer.spillNumber;
    fSpillTrailerTimeStamp = spillTrailer.timeStamp;

    // get timestamp from SpillTrailer, cast as uint64_t
    fTimestamp = (static_cast <std::uint64_t> (spillTrailer.timeStamp));

    // initialize run for the fragment-to-digit algorithm
    fFragmentToDigitAlg.InitializeRun(fRun, fTimestamp);

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

    // flag for PEDESTALON gate
    // PEDESTALON: $00 to $00+999 milliseconds 
    // BEAMON:     $30+1 second to $36 
    fPedestalOn = false;

    // loop through collections
    for (size_t i = 0; i < Collections.size(); ++i) {

      rdu::DataBlockCollection const& Collection = Collections[i];

      size_t const& numberCaenBlocks = Collection.caenBlocks.size();
      size_t const& numberTdcBlocks  = Collection.tdcBlocks.size();

      fNumberCAENBlocks = numberCaenBlocks;
      fNumberTDCBlocks  = numberTdcBlocks;

      // reset interval variables
      fIntervalsDeltaT             = -1;
      fIntervalsDeltaTBeginToBegin = -1;
      fIntervalLength              = -1;
      fIntervalStart               = -1;
      fIntervalStop                = -1;
      fNumberTPCReadouts           = 0;

      // reset data block counters
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

      // clear timestamp vectors
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

      // reset PEDESTALON flag
      fPedestalOn = false;

      std::map< unsigned int, std::vector< double > > CAENTimeStamps;

      fNumberTPCReadouts = Collection.numberTPCReadouts;
      fIntervalLength    = Collection.interval.second - Collection.interval.first;
      fIntervalStart     = Collection.interval.first;
      fIntervalStop      = Collection.interval.second;

      // for pedestal readout, TriggerTimeTag should be
      // 0 < TTT < 124875000 (in between 0 and 0.999 seconds)
      // intervals are in microseconds

      if ((fIntervalStop > 8) and (fIntervalStop < 999000)) {
        fPedestalOn = true;
      }

      //fPedestalOn = true;

      std::cout << "Collection: " << i << std::endl;
      std::cout << "  Number of CAEN data blocks: " << numberCaenBlocks << std::endl;
      std::cout << "  Number of TDC data blocks:  " << numberTdcBlocks << std::endl;
      std::cout << "  Number of TPC readouts:     " << fNumberTPCReadouts << std::endl;

      for (size_t j = 0; j < numberCaenBlocks; ++j) {

        CAENFragment const& caenFrag = Collection.caenBlocks[j];
        double const& timestamp = Collection.caenBlockTimeStamps[j];
        unsigned int boardId = caenFrag.header.boardId;

        std::cout << "    CAEN block: " << j << std::endl;
        std::cout << "      Board ID: " << boardId << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;

        NumberCAENBlocksArray[boardId] += 1;
        CAENTimeStamps[boardId].push_back(timestamp);

        fCaenFragment       = j;
        fCaenEventCounter   = caenFrag.header.eventCounter;
        fCaenEventCounter   = static_cast <uint32_t> (i);
        fCaenBoardId        = caenFrag.header.boardId;
        fCaenTriggerTimeTag = caenFrag.header.triggerTimeTag;
        fCaenNumberSamples  = caenFrag.header.nSamples;

        // fill timestamp histograms
        if (fPedestalOn) {
          fCAENPedestalTimeStampHistograms[boardId]->Fill(timestamp * 1e-6);
        }
        fCAENTimeStampHistograms[boardId]->Fill(timestamp * 1e-6);

        if (fCaenNumberSamples < 1) continue;

        // fill pedestal and ADC histograms
        for (size_t k = 0; k < caenFrag.waveForms.size(); ++k) {
          for (size_t sample = 0; sample < fCaenNumberSamples; ++sample) {
            if (fPedestalOn) {
              fCAENPedestalHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            }
            else {
              fCAENADCHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            }
          }
          short unsigned int minADC = * std::min_element(std::begin(caenFrag.waveForms[k].data),
                                                         std::end(caenFrag.waveForms[k].data));
          short unsigned int maxADC = * std::max_element(std::begin(caenFrag.waveForms[k].data),
                                                         std::end(caenFrag.waveForms[k].data));
          fCAENMinADCHistograms[boardId][k]->Fill(minADC);
          fCAENMaxADCHistograms[boardId][k]->Fill(maxADC);
        }

        if (boardId == 0 or boardId == 1 or boardId == 2 or
            boardId == 3 or boardId == 4 or boardId == 5 or
            boardId == 6 or boardId == 7) {

          for (size_t k = 0; k < V1740_N_CHANNELS; ++k) {
            fCaenV1740Waveform[k].clear();
            fCaenV1740Waveform[k] = caenFrag.waveForms[k].data;

            //for (size_t sample = 0; sample < fCaenNumberSamples; ++sample) {
            //  if (fPedestalOn) {
            //    fCAENPedestalHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //  else {
            //    fCAENADCHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //}
          }

          fCaenV1740DataTree->Fill();

          //if (fCaenNumberSamples != fV1740RecordLength) {
          //  std::cout << "CAEN record length mismatch!"                     << std::endl
          //            << "  CAEN board:             " << boardId            << std::endl
          //            << "  nSamples:               " << fCaenNumberSamples << std::endl
          //            << "  Expected record length: " << fV1740RecordLength << std::endl;
          //}
        }

        else if (boardId == 24) {

          for (size_t k = 0; k < V1740B_N_CHANNELS; ++k) {
            fCaenV1740BWaveform[k].clear();
            fCaenV1740BWaveform[k] = caenFrag.waveForms[k].data;

            //for (size_t sample = 0; sample < fCaenNumberSamples; ++sample) {
            //  if (fPedestalOn) {
            //    fCAENPedestalHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //  else {
            //    fCAENADCHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //}
          }

          fCaenV1740BDataTree->Fill();

          //if (fCaenNumberSamples != fV1740BRecordLength) {
          //  std::cout << "CAEN record length mismatch!"                      << std::endl
          //            << "  CAEN board:             " << boardId             << std::endl
          //            << "  nSamples:               " << fCaenNumberSamples  << std::endl
          //            << "  Expected record length: " << fV1740BRecordLength << std::endl;
          //}
        }

        else if (boardId == 8 or boardId == 9) {

          for (size_t k = 0; k < V1751_N_CHANNELS; ++k) {
            fCaenV1751Waveform[k].clear();
            fCaenV1751Waveform[k] = caenFrag.waveForms[k].data;

            //for (size_t sample = 0; sample < fCaenNumberSamples; ++sample) {
            //  if (fPedestalOn) {
            //    fCAENPedestalHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //  else {
            //    fCAENADCHistograms[boardId][k]->Fill(caenFrag.waveForms[k].data[sample]);
            //  }
            //}
          }

          fCaenV1751DataTree->Fill();

          //if (fCaenNumberSamples != fV1751RecordLength) {
          //  std::cout << "CAEN record length mismatch!"                     << std::endl
          //            << "  CAEN board:             " << boardId            << std::endl
          //            << "  nSamples:               " << fCaenNumberSamples << std::endl
          //            << "  Expected record length: " << fV1751RecordLength << std::endl;
          //}
        }
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

      for (size_t j = 0; j < numberTdcBlocks; ++j) {

        fMwpcTdcNumber.clear();
        fMwpcHitChannel.clear();
        fMwpcHitTimeBin.clear();

        fMwpcNumberHits = 0;
        fMwpcEventCounter = static_cast <uint32_t> (i);

        //std::vector<TDCFragment::TdcEventData> const& tdcEvents = Collection.tdcBlocks[j];
        double const& timestamp = Collection.tdcBlockTimeStamps[j];

        std::cout << "    TDC block: " << j << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
        //std::cout << "      TDC events: " << tdcEvents.size() << std::endl;

        fTDCTimeStamps.push_back(timestamp);

        fMWPCTDCTimeStampHistograms->Fill(timestamp * 1e-6);

        ///////////////////////////////////////////////////////////////////////
        // Loop over TDCs to get a histogram of the TDC timestamps. The
        // timestamp with the most counts is taken as the TDC timestamp of
        // the TDC data block. We are doing this to work around the mismatches
        // of the TDC timestamps.
        ///////////////////////////////////////////////////////////////////////
        std::map<uint32_t, unsigned int> tdcTimeStampCounts;

        std::vector<TDCFragment::TdcEventData> const& tdcDataBlock = Collection.tdcBlocks[j];

        for (size_t tdcIndex = 0; tdcIndex < TDCFragment::MAX_TDCS; ++tdcIndex) {
          TDCFragment::TdcEventData tdcEventData = tdcDataBlock.at(tdcIndex);

          fMwpcTriggerCounter      = tdcEventData.tdcEventHeader.triggerCounter;
          fMwpcTdcTimeStamp        = tdcEventData.tdcEventHeader.tdcTimeStamp;
          fMwpcNumberHits         += tdcEventData.tdcEventHeader.nHits;
          fMwpcControllerTimeStamp = tdcEventData.tdcEventHeader.controllerTimeStamp;

          tdcTimeStampCounts[tdcEventData.tdcEventHeader.tdcTimeStamp] += 1;

          //std::cout << "tdcEventData.tdcEventHeader.tdcTimeStamp:        " << tdcEventData.tdcEventHeader.tdcTimeStamp << std::endl;
          //std::cout << "tdcEventData.tdcEventHeader.controllerTimeStamp: " << tdcEventData.tdcEventHeader.controllerTimeStamp << std::endl;

          for (size_t hitIndex = 0; hitIndex < tdcEventData.tdcHits.size(); ++hitIndex) {
            TDCFragment::TdcHit const& hit = tdcEventData.tdcHits[hitIndex];
            fMwpcHitChannel.push_back(static_cast <uint16_t> (hit.channel));
            fMwpcHitTimeBin.push_back(hit.timeBin);
            fMwpcTdcNumber.push_back(tdcIndex + 1);
          }
        }

        if (fMwpcNumberHits == 0) continue;

        fMwpcTdcDataTree->Fill();

        ///////////////////////////////////////////////////////////////////////
        // Here, we try to figure out which TDCs have mismatch in time bits
        ///////////////////////////////////////////////////////////////////////
        unsigned int counts = 0;
        uint32_t tdcTimeStamp = 0;

        for (auto const& k : tdcTimeStampCounts) {
          if (k.second > counts) {
            tdcTimeStamp = k.first;
            counts = k.second;
          }
        }

        for (size_t tdcIndex = 0; tdcIndex < TDCFragment::MAX_TDCS; ++tdcIndex) {
          TDCFragment::TdcEventData tdcEventData = tdcDataBlock.at(tdcIndex);

          if (tdcEventData.tdcEventHeader.tdcTimeStamp != tdcTimeStamp) {
            std::cout << "Mismatch in TDC time bits!"
                      << "\n  TDC:                                      " << tdcIndex + 1
                      << "\n  tdcEventData.tdcEventHeader.tdcTimeStamp: " << tdcEventData.tdcEventHeader.tdcTimeStamp
                      << "\n  Reference tdcTimeStamp:                   " << tdcTimeStamp
                      << "\n  Difference in tdcTimeStamp:               " << tdcEventData.tdcEventHeader.tdcTimeStamp - tdcTimeStamp
                      << std::endl;
            fTDCTimeBitMismatchHistogram->Fill(tdcIndex + 1);
          }
        }
      }

      if (fNumberTPCReadouts > 0) {
        if (i < (Collections.size() - 1)) {
          //fIntervalsDeltaT = (Collections[i+1].interval.first - Collections[i].interval.second) / 1000.0;
          fIntervalsDeltaT             = (Collections[i+1].interval.first - Collections[i].interval.second);
          fIntervalsDeltaTBeginToBegin = (Collections[i+1].interval.first - Collections[i].interval.first);
          std::cout << "fIntervalsDeltaT [usec]: " << fIntervalsDeltaT << std::endl;
          fIntervalsDeltaTHistogram->Fill(fIntervalsDeltaT / 1000.0);
          fIntervalsDeltaTZHistogram->Fill(fIntervalsDeltaT / 1000.0);
        }
      }

      fNumberTPCReadoutsHistogram->Fill(fNumberTPCReadouts);

      fEventBuilderTree->Fill();

      ////////////////////////////////////////////////////////
      // add digits
      ////////////////////////////////////////////////////////

      std::vector<raw::AuxDetDigit> auxDetDigits;
      std::vector<raw::RawDigit>    rawDigits;
      std::vector<raw::OpDetPulse>  opDetPulses;

      //std::cout << "Collection.caenBlocks.size(): " << Collection.caenBlocks.size() << std::endl;
      //std::cout << "Collection.tdcBlocks.size():  " << Collection.tdcBlocks.size()  << std::endl;
      //std::cout << "auxDetDigits.size(): " << auxDetDigits.size() << std::endl;
      //std::cout << "rawDigits.size():    " << rawDigits.size()    << std::endl;
      //std::cout << "opDetPulses.size():  " << opDetPulses.size()  << std::endl;

      fFragmentToDigitAlg.makeTheDigits(Collection.caenBlocks,
                                        Collection.tdcBlocks,
                                        auxDetDigits,
                                        rawDigits,
                                        opDetPulses);

      ////////////////////////////////////////////////////////
      // begin time of flight shenanigans
      ////////////////////////////////////////////////////////

      std::vector< const raw::AuxDetDigit * > ustofDigits;
      std::vector< const raw::AuxDetDigit * > dstofDigits;

      for (size_t j = 0; j < auxDetDigits.size(); ++j) {
        //std::cout << "AuxDetName(): " << auxDetDigits[j].AuxDetName() << std::endl;
        if (auxDetDigits[j].AuxDetName() == "TOFUS")
          ustofDigits.push_back(&(auxDetDigits[j]));
        if (auxDetDigits[j].AuxDetName() == "TOFDS")
          dstofDigits.push_back(&(auxDetDigits[j]));
      }

      if (ustofDigits.size() == 2 and dstofDigits.size() == 2) {

        // vector for TOF waveforms
        std::vector<short> ustofWaveformA;
        std::vector<short> ustofWaveformB;
        std::vector<short> dstofWaveformA;
        std::vector<short> dstofWaveformB;

        // fill vectors with TOF waveforms
        for (size_t j = 0; j < ustofDigits[0]->NADC(); ++j) {
          ustofWaveformA.push_back(ustofDigits[0]->ADC(j));
        }
        for (size_t j = 0; j < ustofDigits[1]->NADC(); ++j) {
          ustofWaveformB.push_back(ustofDigits[1]->ADC(j));
        }
        for (size_t j = 0; j < dstofDigits[0]->NADC(); ++j) {
          dstofWaveformA.push_back(dstofDigits[0]->ADC(j));
        }
        for (size_t j = 0; j < dstofDigits[1]->NADC(); ++j) {
          dstofWaveformB.push_back(dstofDigits[1]->ADC(j));
        }

        // find hits from TOF waveforms
        std::vector<short> ustofAHits = fTOFBuilderAlg.find_hits(ustofWaveformA);
        std::vector<short> ustofBHits = fTOFBuilderAlg.find_hits(ustofWaveformB);
        std::vector<short> dstofAHits = fTOFBuilderAlg.find_hits(dstofWaveformA);
        std::vector<short> dstofBHits = fTOFBuilderAlg.find_hits(dstofWaveformB);

        // match hits between TOF counters;
        // USTOF1 matched with USTOF2, DSTOF1 matched with DSTOF2
        std::vector<short> ustofHits = fTOFBuilderAlg.match_hits(ustofAHits, ustofBHits);
        std::vector<short> dstofHits = fTOFBuilderAlg.match_hits(dstofAHits, dstofBHits);

        //std::cout << "ustofHits.size(): " << ustofHits.size() << std::endl;
        //std::cout << "dstofHits.size(): " << dstofHits.size() << std::endl;

        for (size_t j = 0; j < ustofHits.size(); ++j) {
          fUSTOFHitsHistogram->Fill(ustofHits[j]);
        }

        for (size_t j = 0; j < dstofHits.size(); ++j) {
          fDSTOFHitsHistogram->Fill(dstofHits[j]);
        }

        std::pair< std::vector<short>, std::vector<long> > tofAndTimeStamp;
        tofAndTimeStamp = fTOFBuilderAlg.get_TOF_and_TimeStamp(ustofDigits, dstofDigits);

        //std::cout << "tofAndTimeStamp.first.size(): " << tofAndTimeStamp.first.size() << std::endl;

        for (size_t j = 0; j < tofAndTimeStamp.first.size(); ++j) {
          //std::cout << "tofAndTimeStamp.first.at(j): " << tofAndTimeStamp.first.at(j) << std::endl;
          fTOFHistogram->Fill(tofAndTimeStamp.first.at(j));
        }

      }

    }

    std::vector< rdu::DataBlock > DataBlocks;
    DataBlocks = fEventBuilderAlg.GetDataBlocks(fLariatFragment);
    std::vector< std::pair< double, double > > CAENBoard0Intervals;
    //CAENBoard0Intervals = fEventBuilderAlg.CreateIntervals(DataBlocks, 0, fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow, fV1740AcquisitionWindow);
    CAENBoard0Intervals = fEventBuilderAlg.CreateIntervals(DataBlocks, 0, 0.128, 0.128, fV1740AcquisitionWindow);
    CAENBoard0Intervals = fEventBuilderAlg.IntervalsSelfMerge(CAENBoard0Intervals);

    //fTPCIntervalsDeltaT = -999999999999999999;
    fTPCIntervalsDeltaT             = -1;
    fTPCIntervalsDeltaTBeginToBegin = -1;

    if (CAENBoard0Intervals.size() > 1) {
      for (size_t i = 0; i < (CAENBoard0Intervals.size() - 1); ++i) {
        fTPCIntervalsDeltaT = -1;
        //fTPCIntervalsDeltaT = (CAENBoard0Intervals.at(i+1).first - CAENBoard0Intervals.at(i).second) / 1000.0;
        fTPCIntervalsDeltaT             = (CAENBoard0Intervals.at(i+1).first - CAENBoard0Intervals.at(i).second);
        fTPCIntervalsDeltaTBeginToBegin = (CAENBoard0Intervals.at(i+1).first - CAENBoard0Intervals.at(i).first);
        //std::cout << "fTPCIntervalsDeltaT [usec]: " << fTPCIntervalsDeltaT << std::endl;
        fTPCIntervalsDeltaTHistogram->Fill(fTPCIntervalsDeltaT / 1000.0);
        fTPCIntervalsDeltaTZHistogram->Fill(fTPCIntervalsDeltaT / 1000.0);
        fTPCTree->Fill();
      }
    }

    size_t const& numberWutFrags = fLariatFragment->wutFrags.size();
    for (size_t i = 0; i < numberWutFrags; ++i) {
      WUTFragment const& wutFrag = fLariatFragment->wutFrags[i];
      std::vector<WUTFragment::WutHit> const& wutHits = wutFrag.hits;
      fWutNumberHits = wutFrag.header.nHits;
      fWutHitChannel.clear();
      fWutHitTimeBin.clear();
      fWutTimeHeader = wutFrag.header.timeHeader;
      fWUTTimeStampHistograms->Fill(fWutTimeHeader * 16e-6);
      for (size_t j = 0; j < wutFrag.hits.size(); ++j) {
        WUTFragment::WutHit const& wutHit = wutHits[j];
        fWutHitChannel.push_back(static_cast <uint16_t> (wutHit.channel));
        fWutHitTimeBin.push_back(wutHit.timeBin);
      }
      fWutDataTree->Fill();
    }

    fEventRecord->Fill();
    fSpillTrailerTree->Fill();

    return;
  }

  //-----------------------------------------------------------------------
  size_t DataQuality::castToSizeT_(std::string const& string)
  {
    size_t sizeT;

    if (!string.empty()) {
        sizeT = static_cast <size_t> (std::stoi(string));
    }
    else {
      sizeT = 0;
    }

    return sizeT;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(DataQuality)

} // namespace DataQuality

#endif // DataQuality_Module
