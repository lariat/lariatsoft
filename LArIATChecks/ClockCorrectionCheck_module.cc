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
    rdu::ClockCorrectionAlg fClockCorrectionAlg;

    // vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

    // sample reduction of the CAEN V1740 digitizers
    size_t fV1740SampleReduction;
    size_t fV1740BSampleReduction;

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

    // TDC parameters
    double fTDCPipelineDelay;
    double fTDCGateWidth;

    double fV1740PreTriggerWindow;
    double fV1740PostTriggerWindow;
    double fV1740BPreTriggerWindow;
    double fV1740BPostTriggerWindow;
    double fV1751PreTriggerWindow;
    double fV1751PostTriggerWindow;
    double fTDCPreTriggerWindow;
    double fTDCPostTriggerWindow;

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

    fConfigParams.push_back("v1740_config_caen_recordlength");
    fConfigParams.push_back("v1740b_config_caen_recordlength");
    fConfigParams.push_back("v1751_config_caen_recordlength");
    fConfigParams.push_back("v1740_config_caen_v1740_samplereduction");
    fConfigParams.push_back("v1740b_config_caen_v1740_samplereduction");
    fConfigParams.push_back("tdc_config_tdc_pipelinedelay");
    fConfigParams.push_back("tdc_config_tdc_gatewidth");
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::beginRun(const art::Run& run)
  {
    fRun = run.run();
    fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams, fRun);

    fV1740RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_recordlength"]));
    fV1740BRecordLength    = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_recordlength"]));
    fV1751RecordLength     = static_cast <size_t> (std::stoi(fConfigValues["v1751_config_caen_recordlength"]));
    fV1740SampleReduction  = static_cast <size_t> (std::stoi(fConfigValues["v1740_config_caen_v1740_samplereduction"]));
    fV1740BSampleReduction = static_cast <size_t> (std::stoi(fConfigValues["v1740b_config_caen_v1740_samplereduction"]));
    fTDCPipelineDelay      = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_pipelinedelay"]));
    fTDCGateWidth          = static_cast <size_t> (std::stoi(fConfigValues["tdc_config_tdc_gatewidth"]));

    // sampling rate in MHz
    fV1740SamplingRate  = 62.5 / fV1740SampleReduction;
    fV1740BSamplingRate = 62.5 / fV1740BSampleReduction;
    fV1751SamplingRate  = 1e3;

    // readout window length in microseconds
    fV1740ReadoutWindow  = fV1740RecordLength  / fV1740SamplingRate;
    fV1740BReadoutWindow = fV1740BRecordLength / fV1740BSamplingRate;
    fV1751ReadoutWindow  = fV1751RecordLength  / fV1751SamplingRate;

    // Pre-/post-trigger windows
    if (fV1740PreTriggerWindow   < 0) fV1740PreTriggerWindow   = fV1740ReadoutWindow;
    if (fV1740PostTriggerWindow  < 0) fV1740PostTriggerWindow  = fV1740ReadoutWindow;
    if (fV1740BPreTriggerWindow  < 0) fV1740BPreTriggerWindow  = fV1740BReadoutWindow;
    if (fV1740BPostTriggerWindow < 0) fV1740BPostTriggerWindow = fV1740BReadoutWindow;
    if (fV1751PreTriggerWindow   < 0) fV1751PreTriggerWindow   = 0.64;
    if (fV1751PostTriggerWindow  < 0) fV1751PostTriggerWindow  = fV1751ReadoutWindow + 0.64;
    if (fTDCPreTriggerWindow     < 0) fTDCPreTriggerWindow     = 1.196;
    if (fTDCPostTriggerWindow    < 0) fTDCPostTriggerWindow    = 1.196;

    std::cout << "//////////////////////////////////////////////"        << std::endl;
    std::cout << "V1740PreTriggerWindow:   " << fV1740PreTriggerWindow   << std::endl;
    std::cout << "V1740PostTriggerWindow:  " << fV1740PostTriggerWindow  << std::endl;
    std::cout << "V1740BPreTriggerWindow:  " << fV1740BPreTriggerWindow  << std::endl;
    std::cout << "V1740BPostTriggerWindow: " << fV1740BPostTriggerWindow << std::endl;
    std::cout << "V1751PreTriggerWindow:   " << fV1751PreTriggerWindow   << std::endl;
    std::cout << "V1751PostTriggerWindow:  " << fV1751PostTriggerWindow  << std::endl;
    std::cout << "TDCPreTriggerWindow:     " << fTDCPreTriggerWindow     << std::endl;
    std::cout << "TDCPostTriggerWindow:    " << fTDCPostTriggerWindow    << std::endl;
    std::cout << "//////////////////////////////////////////////"        << std::endl;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::beginSubRun(const art::SubRun& subrun)
  {}

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClockCorrectionAlg.reconfigure(pset.get<fhicl::ParameterSet>("ClockCorrectionAlg"));
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
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionCheck::analyze(const art::Event& event) 
  {
    // make the utility to access the fragments from the event record
    rdu::FragmentUtility fragUtil(event, fRawFragmentLabel, fRawFragmentInstance);

    std::vector< rdu::DataBlock > DataBlocks;
    std::map< unsigned int, std::vector< double > > TimeStampMap;
    fClockCorrectionAlg.GetDataBlocksTimeStampMap(&fragUtil.DAQFragment(), DataBlocks, TimeStampMap);

    for (std::map< unsigned int, std::vector< double> >::const_iterator
         iter = TimeStampMap.begin(); iter != TimeStampMap.end(); ++iter) {
      mf::LogVerbatim("ClockCorrectionAlg") << "Device ID: "
                                            << iter->first
                                            << "; number of data blocks: "
                                            << iter->second.size();
    }

    // get clock correction parameters
    std::map< unsigned int, std::pair< double, double > > ClockCorrectionParameters;
    fClockCorrectionAlg.GetClockCorrectionParameters(TimeStampMap, ClockCorrectionParameters);

    // apply clock correction to DataBlock timestamps
    for (size_t block_idx = 0; block_idx < DataBlocks.size(); ++block_idx) {
      rdu::DataBlock & block = DataBlocks.at(block_idx);

      unsigned int DeviceID = block.deviceId;

      // correction parameters
      double Slope = ClockCorrectionParameters[DeviceID].first;
      double Intercept = ClockCorrectionParameters[DeviceID].second;

      // apply correction
      block.correctedTimestamp = (block.timestamp - Intercept) / Slope;
    } // end loop over data blocks

    // sort data blocks by their corrected timestamps
    std::sort(DataBlocks.begin(), DataBlocks.end(),
              [] (rdu::DataBlock const& a, rdu::DataBlock const& b) {
                return a.correctedTimestamp < b.correctedTimestamp;
              });

    //for (size_t block_idx = 0; block_idx < DataBlocks.size(); ++block_idx) {
    //  rdu::DataBlock const& block = DataBlocks.at(block_idx);
    //  std::cout << block.timestamp << ", " << block.correctedTimestamp << ", " << block.deviceId << std::endl;
    //} // end loop over data blocks

    //for (std::vector< rdu::DataBlock >::const_reverse_iterator
    //     block = DataBlocks.rbegin(); block != DataBlocks.rend(); ++block) {
    //  std::cout << std::setfill(' ') << std::setw(3)
    //            << block->deviceId << "; " << block->correctedTimestamp
    //            << std::endl;
    //} // end loop over data blocks

    //for (size_t block_idx = DataBlocks.size(); block_idx-- > 0;) { // What is this sorcery?
    //  rdu::DataBlock const& block = DataBlocks.at(block_idx);
    //  //std::cout << block_idx << std::endl;
    //  std::cout << std::setfill(' ') << std::setw(3)
    //            << block.deviceId << "; " << block.correctedTimestamp
    //            << std::endl;
    //} // end loop over data blocks

    // initialize vector of intervals
    std::vector< std::pair< double, double > > CAENBoard0Intervals;
    std::vector< std::pair< double, double > > CAENBoard1Intervals;
    std::vector< std::pair< double, double > > CAENBoard2Intervals;
    std::vector< std::pair< double, double > > CAENBoard3Intervals;
    std::vector< std::pair< double, double > > CAENBoard4Intervals;
    std::vector< std::pair< double, double > > CAENBoard5Intervals;
    std::vector< std::pair< double, double > > CAENBoard6Intervals;
    std::vector< std::pair< double, double > > CAENBoard7Intervals;
    std::vector< std::pair< double, double > > CAENBoard8Intervals;
    std::vector< std::pair< double, double > > CAENBoard9Intervals;
    std::vector< std::pair< double, double > > CAENBoard24Intervals;
    std::vector< std::pair< double, double > > TDCIntervals;

    // CAEN V1740 digitizers
    CAENBoard0Intervals  = this->CreateIntervals(DataBlocks,  0,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard1Intervals  = this->CreateIntervals(DataBlocks,  1,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard2Intervals  = this->CreateIntervals(DataBlocks,  2,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard3Intervals  = this->CreateIntervals(DataBlocks,  3,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard4Intervals  = this->CreateIntervals(DataBlocks,  4,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard5Intervals  = this->CreateIntervals(DataBlocks,  5,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard6Intervals  = this->CreateIntervals(DataBlocks,  6,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    CAENBoard7Intervals  = this->CreateIntervals(DataBlocks,  7,  fV1740PreTriggerWindow, fV1740PostTriggerWindow);
    // CAEN V1751 digitizers
    CAENBoard8Intervals  = this->CreateIntervals(DataBlocks,  8,  fV1751PreTriggerWindow, fV1751PostTriggerWindow);
    CAENBoard9Intervals  = this->CreateIntervals(DataBlocks,  9,  fV1751PreTriggerWindow, fV1751PostTriggerWindow);
    // CAEN V1740B digitizer ("spare" CAEN V1740 digitizer)
    CAENBoard24Intervals = this->CreateIntervals(DataBlocks, 24, fV1740BPreTriggerWindow, fV1740BPostTriggerWindow);

    // WC TDC
    // see TDC readout documentation here:
    // https://cdcvs.fnal.gov/redmine/projects/lariat-online/wiki/TDC_Readout_Documentation
    // TODO: This is temporary. We will need to figure out
    //       how to convert tdc_config_tdc_gatewidth and
    //       tdc_config_tdc_pipelinedelay into a TDC
    //       readout window length.
    TDCIntervals = this->CreateIntervals(DataBlocks, 32, fTDCPreTriggerWindow, fTDCPostTriggerWindow);

    // self-merge intervals
    CAENBoard0Intervals  = this->IntervalsSelfMerge(CAENBoard0Intervals);
    CAENBoard1Intervals  = this->IntervalsSelfMerge(CAENBoard1Intervals);
    CAENBoard2Intervals  = this->IntervalsSelfMerge(CAENBoard2Intervals);
    CAENBoard3Intervals  = this->IntervalsSelfMerge(CAENBoard3Intervals);
    CAENBoard4Intervals  = this->IntervalsSelfMerge(CAENBoard4Intervals);
    CAENBoard5Intervals  = this->IntervalsSelfMerge(CAENBoard5Intervals);
    CAENBoard6Intervals  = this->IntervalsSelfMerge(CAENBoard6Intervals);
    CAENBoard7Intervals  = this->IntervalsSelfMerge(CAENBoard7Intervals);
    CAENBoard8Intervals  = this->IntervalsSelfMerge(CAENBoard8Intervals);
    CAENBoard9Intervals  = this->IntervalsSelfMerge(CAENBoard9Intervals);
    CAENBoard24Intervals = this->IntervalsSelfMerge(CAENBoard24Intervals);
    TDCIntervals         = this->IntervalsSelfMerge(TDCIntervals);

    std::cout << "//////////////////////////////////////////////" << std::endl;
    std::cout << "Pre-merge stage" << std::endl;
    std::cout << "//////////////////////////////////////////////" << std::endl;

    std::cout << "CAEN digitizer board 0" << std::endl;
    std::cout << "//////////////////////////////////////////////" << std::endl;

    for (std::vector< std::pair< double, double > >::const_iterator
         iter = CAENBoard0Intervals.begin(); iter != CAENBoard0Intervals.end(); ++iter) {
      std::cout << iter->first << ", " << iter->second << std::endl;
    }

    std::cout << "//////////////////////////////////////////////" << std::endl;
    std::cout << "CAEN digitizer board 8" << std::endl;
    std::cout << "//////////////////////////////////////////////" << std::endl;

    for (std::vector< std::pair< double, double > >::const_iterator
         iter = CAENBoard8Intervals.begin(); iter != CAENBoard8Intervals.end(); ++iter) {
      std::cout << iter->first << ", " << iter->second << std::endl;
    }

    std::cout << "//////////////////////////////////////////////" << std::endl;
    std::cout << "TDC" << std::endl;
    std::cout << "//////////////////////////////////////////////" << std::endl;

    for (std::vector< std::pair< double, double > >::const_iterator
         iter = TDCIntervals.begin(); iter != TDCIntervals.end(); ++iter) {
      std::cout << iter->first << ", " << iter->second << std::endl;
    }

    // vector of merged intervals
    std::vector< std::pair< double, double > > MergedIntervals;

    // merge the intervals!
    MergedIntervals = this->MergeIntervals(CAENBoard0Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard1Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard2Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard3Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard4Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard5Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard6Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard7Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard8Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard9Intervals,  MergedIntervals);
    MergedIntervals = this->MergeIntervals(CAENBoard24Intervals, MergedIntervals);
    MergedIntervals = this->MergeIntervals(TDCIntervals,         MergedIntervals);

    std::cout << "//////////////////////////////////////////////" << std::endl;
    std::cout << "Post-merge stage" << std::endl;
    std::cout << "//////////////////////////////////////////////" << std::endl;

    for (std::vector< std::pair< double, double > >::const_iterator
         iter = MergedIntervals.begin(); iter != MergedIntervals.end(); ++iter) {
      std::cout << iter->first << ", " << iter->second << std::endl;
    }

    std::cout << "//////////////////////////////////////////////" << std::endl;

    //////////////////////////////////////////////////////////
    // Cue Queen & David Bowie's Under Pressure:
    //   Slice, slice, baby
    //////////////////////////////////////////////////////////

    // group data blocks into collections
    std::vector< rdu::DataBlockCollection > Collections;

    for (std::vector< std::pair< double, double > >::const_iterator
         iter = MergedIntervals.begin(); iter != MergedIntervals.end(); ++iter) {

      // get interval
      double const& t_a = iter->first;
      double const& t_b = iter->second;

      rdu::DataBlockCollection Collection;
      Collection.interval = std::make_pair(t_a, t_b);

      // number of data blocks in this interval
      size_t NumberDataBlocks = 0;

      for (size_t block_idx = 0; block_idx < DataBlocks.size(); ++block_idx) {

        rdu::DataBlock const& block              = DataBlocks.at(block_idx);
        //unsigned int   const& DeviceID           = block.deviceId;
        double         const& correctedTimestamp = block.correctedTimestamp;

        if ((correctedTimestamp > t_a) and (correctedTimestamp < t_b)) {
          // there should be at most one data block for each
          // DataBlock struct

          for (size_t i = 0; i < block.caenBlocks.size(); ++i) {
            const CAENFragment * caenFrag = block.caenBlocks[i];
            //Collection.caenBlocks.push_back(caenFrag);
            Collection.caenBlocks.push_back(std::make_pair(correctedTimestamp, caenFrag));

            // increment number of data blocks
            ++NumberDataBlocks;
          } // end loop over CAEN data blocks

          for (size_t i = 0; i < block.tdcBlocks.size(); ++i) {
            const std::vector<TDCFragment::TdcEventData> * tdcEvents = block.tdcBlocks[i];
            //Collection.tdcBlocks.push_back(tdcEvents);
            Collection.tdcBlocks.push_back(std::make_pair(correctedTimestamp, tdcEvents));

            // increment number of data blocks
            ++NumberDataBlocks;
          } // end loop over TDC data blocks
        }

      } // end loop through data blocks

      // add collection to vector of collections only if
      // there are data blocks present
      if (NumberDataBlocks > 0) Collections.push_back(Collection);

    } // end loop through merged intervals

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

  //-----------------------------------------------------------------------
  std::vector< std::pair< double, double> > ClockCorrectionCheck::CreateIntervals(std::vector< rdu::DataBlock > const& DataBlocks,
                                                                                  unsigned int                  const& DeviceID,
                                                                                  double                        const& PreTriggerWindow,
                                                                                  double                        const& PostTriggerWindow)
  {
    std::vector< std::pair< double, double > > Intervals;

    for (std::vector< rdu::DataBlock >::const_iterator
         block = DataBlocks.begin(); block != DataBlocks.end(); ++block) {
      if (block->deviceId != DeviceID) continue;
      //std::cout << block->timestamp << ", " << block->correctedTimestamp << ", " << block->deviceId << std::endl;
      double IntervalLow = block->correctedTimestamp - PreTriggerWindow;
      double IntervalHigh = block->correctedTimestamp + PostTriggerWindow;
      Intervals.push_back(std::make_pair(IntervalLow, IntervalHigh));
    }

    return Intervals;
  }

  //-----------------------------------------------------------------------
  std::vector< std::pair< double, double > > ClockCorrectionCheck::IntervalsSelfMerge(std::vector< std::pair< double, double > > const& Intervals)
  {
    // vector merged intervals
    std::vector< std::pair< double, double > > MergedIntervals;

    // bookkeeping for keeping track of which intervals have
    // already been merged
    std::vector<size_t> MergedIndices;

    for (size_t i = 0; i < Intervals.size(); ++i) {

      // skip if interval has already been merged
      if (std::find(MergedIndices.begin(), MergedIndices.end(), i) != MergedIndices.end())
        continue;

      // get interval
      double t_a = Intervals[i].first;
      double t_b = Intervals[i].second;

      for (size_t j = 0; j < Intervals.size(); ++j) {

        // skip if interval has already been merged or if
        // this interval is the same as the interval in
        // the parent loop
        if ((i == j) or
            (std::find(MergedIndices.begin(), MergedIndices.end(), j) != MergedIndices.end()))
          continue;

        // get interval
        double t_c = Intervals[j].first;
        double t_d = Intervals[j].second;

        //std::cout << t_a << ", " << t_b << "; " << t_c << ", " << t_d << std::endl;

        if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_b)) {
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndices.push_back(j);
          //std::cout << "t_a = t_a; t_b = t_d" << std::endl;
        }
        else if ((t_c <= t_a) and (t_d >= t_a) and (t_d <= t_b)) {
          t_a = t_c;
          // add interval index for bookkeeping
          MergedIndices.push_back(j);
          //std::cout << "t_a = t_c; t_b = t_b" << std::endl;
        }
        else if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d <= t_b)) {
          // add interval index for bookkeeping
          MergedIndices.push_back(j);
          //std::cout << "t_a = t_a; t_b = t_b" << std::endl;
        }
        else if ((t_c <= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d >= t_b)) {
          t_a = t_c;
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndices.push_back(j);
          //std::cout << "t_a = t_c; t_b = t_d" << std::endl;
        }

      } // end loop over intervals

      // add interval to vector of intervals
      MergedIntervals.push_back(std::make_pair(t_a, t_b));

    } // end loop over intervals

    // sort intervals in ascending order
    std::sort(MergedIntervals.begin(), MergedIntervals.end());

    return MergedIntervals;
  }

  //-----------------------------------------------------------------------
  std::vector< std::pair< double, double > > ClockCorrectionCheck::MergeIntervals(std::vector< std::pair< double, double > > const& IntervalsA,
                                                                                  std::vector< std::pair< double, double > > const& IntervalsB)
  {
    // vector merged intervals
    std::vector< std::pair< double, double > > MergedIntervals;

    // bookkeeping for keeping track of which intervals have
    // already been merged
    std::vector<size_t> MergedIndicesA;
    std::vector<size_t> MergedIndicesB;

    for (size_t i = 0; i < IntervalsA.size(); ++i) {

      // skip if interval has already been merged
      if (std::find(MergedIndicesA.begin(), MergedIndicesA.end(), i) != MergedIndicesA.end())
        continue;

      // get interval
      double t_a = IntervalsA[i].first;
      double t_b = IntervalsA[i].second;

      for (size_t j = 0; j < IntervalsB.size(); ++j) {

        // skip if interval has already been merged
        if (std::find(MergedIndicesB.begin(), MergedIndicesB.end(), j) != MergedIndicesB.end())
          continue;

        // get interval
        double t_c = IntervalsB[j].first;
        double t_d = IntervalsB[j].second;

        //std::cout << t_a << ", " << t_b << "; " << t_c << ", " << t_d << std::endl;

        if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_b)) {
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndicesB.push_back(j);
          //std::cout << "t_a = t_a; t_b = t_d" << std::endl;
        }
        else if ((t_c <= t_a) and (t_d >= t_a) and (t_d <= t_b)) {
          t_a = t_c;
          // add interval index for bookkeeping
          MergedIndicesB.push_back(j);
          //std::cout << "t_a = t_c; t_b = t_b" << std::endl;
        }
        else if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d <= t_b)) {
          // add interval index for bookkeeping
          MergedIndicesB.push_back(j);
          //std::cout << "t_a = t_a; t_b = t_b" << std::endl;
        }
        else if ((t_c <= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d >= t_b)) {
          t_a = t_c;
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndicesB.push_back(j);
          //std::cout << "t_a = t_c; t_b = t_d" << std::endl;
        }

      } // end loop over IntervalsB

      // add interval index for bookkeeping
      MergedIndicesA.push_back(i);

      // add interval to vector of intervals
      MergedIntervals.push_back(std::make_pair(t_a, t_b));

    } // end loop over IntervalsA

    // TODO: TURN THIS INTO A FUNCTION. MAYBE. IDK.

    for (size_t i = 0; i < IntervalsB.size(); ++i) {

      // skip if interval has already been merged
      if (std::find(MergedIndicesB.begin(), MergedIndicesB.end(), i) != MergedIndicesB.end())
        continue;

      double t_a = IntervalsB[i].first;
      double t_b = IntervalsB[i].second;

      for (size_t j = 0; j < IntervalsA.size(); ++j) {

        // skip if interval has already been merged
        if (std::find(MergedIndicesA.begin(), MergedIndicesA.end(), j) != MergedIndicesA.end())
          continue;

        double t_c = IntervalsA[j].first;
        double t_d = IntervalsA[j].second;

        if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_b)) {
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndicesA.push_back(j);
        }
        else if ((t_c <= t_a) and (t_d >= t_a) and (t_d <= t_b)) {
          t_a = t_c;
          // add interval index for bookkeeping
          MergedIndicesA.push_back(j);
        }
        else if ((t_c >= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d <= t_b)) {
          // add interval index for bookkeeping
          MergedIndicesA.push_back(j);
        }
        else if ((t_c <= t_a) and (t_c <= t_b) and (t_d >= t_a) and (t_d >= t_b)) {
          t_a = t_c;
          t_b = t_d;
          // add interval index for bookkeeping
          MergedIndicesA.push_back(j);
        }

      } // end loop over IntervalsA

      // add interval index for bookkeeping
      MergedIndicesB.push_back(i);

      // add interval to vector of intervals
      MergedIntervals.push_back(std::make_pair(t_a, t_b));

    } // end loop over IntervalsB

    // sort intervals in ascending order
    std::sort(MergedIntervals.begin(), MergedIntervals.end());

    return MergedIntervals;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(ClockCorrectionCheck)

} // namespace ClockCorrectionCheck

#endif // ClockCorrectionCheck_Module
