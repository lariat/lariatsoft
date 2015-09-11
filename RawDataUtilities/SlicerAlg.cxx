//////////////////////////////////////////////////////////////
// Name:      SlicerAlg.cxx
// Date:      20 August 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Event slicing algorithm for LArIAT.
//////////////////////////////////////////////////////////////

// Class include
#include "RawDataUtilities/SlicerAlg.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>

namespace rdu {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  SlicerAlg::SlicerAlg(fhicl::ParameterSet const& pset)
    : fClockCorrectionAlg(pset.get<fhicl::ParameterSet>("ClockCorrectionAlg"))
  {
    // read in parameters from .fcl files
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  SlicerAlg::~SlicerAlg() {}

  //-----------------------------------------------------------------------
  void SlicerAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTPCReadoutBufferLow  = pset.get< double >("TPCReadoutBufferLow",  0.128);
    fTPCReadoutBufferHigh = pset.get< double >("TPCReadoutBufferHigh", 0.128);

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerAlg::Configure(double V1740PreAcquisitionWindow,
                            double V1740PostAcquisitionWindow,
                            double V1740AcquisitionWindow,
                            double V1740BPreAcquisitionWindow,
                            double V1740BPostAcquisitionWindow,
                            double V1740BAcquisitionWindow,
                            double V1751PreAcquisitionWindow,
                            double V1751PostAcquisitionWindow,
                            double V1751AcquisitionWindow,
                            double TDCPreAcquisitionWindow,
                            double TDCPostAcquisitionWindow,
                            double TDCAcquisitionWindow)
  {
    fV1740PreAcquisitionWindow   = V1740PreAcquisitionWindow;
    fV1740PostAcquisitionWindow  = V1740PostAcquisitionWindow;
    fV1740AcquisitionWindow      = V1740AcquisitionWindow;
    fV1740BPreAcquisitionWindow  = V1740BPreAcquisitionWindow;
    fV1740BPostAcquisitionWindow = V1740BPostAcquisitionWindow;
    fV1740BAcquisitionWindow     = V1740BAcquisitionWindow;
    fV1751PreAcquisitionWindow   = V1751PreAcquisitionWindow;
    fV1751PostAcquisitionWindow  = V1751PostAcquisitionWindow;
    fV1751AcquisitionWindow      = V1751AcquisitionWindow;
    fTDCPreAcquisitionWindow     = TDCPreAcquisitionWindow;
    fTDCPostAcquisitionWindow    = TDCPostAcquisitionWindow;
    fTDCAcquisitionWindow        = TDCAcquisitionWindow;

    return;
  }

  //-----------------------------------------------------------------------
  std::vector< std::pair< double, double> > SlicerAlg::CreateIntervals(std::vector< rdu::DataBlock > const& DataBlocks,
                                                                       unsigned int                  const& DeviceID,
                                                                       double                        const& PreAcquisitionWindow,
                                                                       double                        const& PostAcquisitionWindow,
                                                                       double                        const& AcquisitionWindow)
  {
    std::vector< std::pair< double, double > > Intervals;

    for (std::vector< rdu::DataBlock >::const_iterator
         block = DataBlocks.begin(); block != DataBlocks.end(); ++block) {
      if (block->deviceId != DeviceID) continue;
      if (block->correctedTimestamp < 0) continue;
      //std::cout << block->timestamp << ", " << block->correctedTimestamp << ", " << block->deviceId << std::endl;
      double IntervalLow = block->correctedTimestamp - PreAcquisitionWindow;
      double IntervalHigh = block->correctedTimestamp + AcquisitionWindow + PostAcquisitionWindow;
      Intervals.push_back(std::make_pair(IntervalLow, IntervalHigh));
    }

    return Intervals;
  }

  //-----------------------------------------------------------------------
  std::vector< std::pair< double, double > > SlicerAlg::IntervalsSelfMerge(std::vector< std::pair< double, double > > const& Intervals)
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
  std::vector< std::pair< double, double > > SlicerAlg::MergeIntervals(std::vector< std::pair< double, double > > const& IntervalsA,
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

    // self-merging of intervals
    MergedIntervals = this->IntervalsSelfMerge(MergedIntervals);

    // sort intervals in ascending order
    std::sort(MergedIntervals.begin(), MergedIntervals.end());

    return MergedIntervals;
  }

  //-----------------------------------------------------------------------
  std::vector< rdu::DataBlock > SlicerAlg::GetDataBlocks(const LariatFragment * data)
  //void SlicerAlg::GetDataBlocks(const LariatFragment * data, std::vector< rdu::DataBlock > & DataBlocks)
  {

    std::vector< rdu::DataBlock > DataBlocks;
    std::map< unsigned int, std::vector< double > > TimeStampMap;
    fClockCorrectionAlg.GetDataBlocksTimeStampMap(data, DataBlocks, TimeStampMap);

    //std::cout << "DataBlocks.size():   " << DataBlocks.size()   << std::endl;
    //std::cout << "TimeStampMap.size(): " << TimeStampMap.size() << std::endl;

    for (std::map< unsigned int, std::vector< double> >::const_iterator
         iter = TimeStampMap.begin(); iter != TimeStampMap.end(); ++iter) {
      mf::LogVerbatim("SlicerAlg") << "Device ID: "
                                   << iter->first
                                   << "; number of data blocks: "
                                   << iter->second.size();
    }

    // get clock correction parameters
    std::map< unsigned int, std::pair< double, double > > ClockCorrectionParameters;
    unsigned ReferenceClockDeviceID;

    fClockCorrectionAlg.GetClockCorrectionParameters(TimeStampMap,
                                                     ClockCorrectionParameters,
                                                     ReferenceClockDeviceID);

    // apply clock correction to DataBlock timestamps
    for (size_t block_idx = 0; block_idx < DataBlocks.size(); ++block_idx) {
      rdu::DataBlock & block = DataBlocks.at(block_idx);

      unsigned int DeviceID = block.deviceId;

      if (ClockCorrectionParameters.count(DeviceID) < 1) continue;

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

    return DataBlocks;
    //return;
  }

  //-----------------------------------------------------------------------
  std::vector< rdu::DataBlockCollection > SlicerAlg::Slice(const LariatFragment * data)
  //void SlicerAlg::Slice(const LariatFragment * data, std::vector< rdu::DataBlockCollection > & Collections)
  {

    std::vector< rdu::DataBlock > DataBlocks;
    DataBlocks = this->GetDataBlocks(data);
    //this->GetDataBlocks(data, DataBlocks);

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
    CAENBoard0Intervals  = this->CreateIntervals(DataBlocks,  0,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard1Intervals  = this->CreateIntervals(DataBlocks,  1,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard2Intervals  = this->CreateIntervals(DataBlocks,  2,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard3Intervals  = this->CreateIntervals(DataBlocks,  3,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard4Intervals  = this->CreateIntervals(DataBlocks,  4,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard5Intervals  = this->CreateIntervals(DataBlocks,  5,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard6Intervals  = this->CreateIntervals(DataBlocks,  6,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    CAENBoard7Intervals  = this->CreateIntervals(DataBlocks,  7,  fV1740PreAcquisitionWindow, fV1740PostAcquisitionWindow,  fV1740AcquisitionWindow);
    // CAEN V1751 digitizers
    CAENBoard8Intervals  = this->CreateIntervals(DataBlocks,  8,  fV1751PreAcquisitionWindow, fV1751PostAcquisitionWindow,  fV1751AcquisitionWindow);
    CAENBoard9Intervals  = this->CreateIntervals(DataBlocks,  9,  fV1751PreAcquisitionWindow, fV1751PostAcquisitionWindow,  fV1751AcquisitionWindow);
    // CAEN V1740B digitizer ("spare" CAEN V1740 digitizer)
    CAENBoard24Intervals = this->CreateIntervals(DataBlocks, 24, fV1740BPreAcquisitionWindow, fV1740BPostAcquisitionWindow, fV1740BAcquisitionWindow);

    // WC TDC
    // see TDC readout documentation here:
    // https://cdcvs.fnal.gov/redmine/projects/lariat-online/wiki/TDC_Readout_Documentation
    // TODO: This is temporary. We will need to figure out
    //       how to convert tdc_config_tdc_gatewidth and
    //       tdc_config_tdc_pipelinedelay into a TDC
    //       readout window length.
    TDCIntervals = this->CreateIntervals(DataBlocks, rdu::TDC_DEVICE_ID, fTDCPreAcquisitionWindow, fTDCPostAcquisitionWindow, fTDCAcquisitionWindow);

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
      rdu::DataBlockCollection & Collection = Collections[i];

      size_t const& NumberCaenBlocks = Collection.caenBlocks.size();
      size_t const& NumberTdcBlocks = Collection.tdcBlocks.size();

      std::cout << "Collection: " << i << std::endl;
      std::cout << "  Number of CAEN data blocks: " << NumberCaenBlocks << std::endl;
      std::cout << "  Number of TDC data blocks:  " << NumberTdcBlocks << std::endl;

      // used for counting the number of TPC readouts
      std::vector< std::pair< double, double > > TPCReadout;
      size_t NumberTPCReadouts = 0;

      for (size_t j = 0; j < NumberCaenBlocks; ++j) {

        std::pair< double, const CAENFragment * > caenBlock = Collection.caenBlocks[j];
        double const& timestamp = caenBlock.first;
        const CAENFragment * caenFrag = caenBlock.second;
        unsigned int boardId = caenFrag->header.boardId;

        std::cout << "    CAEN block: " << j << std::endl;
        std::cout << "      Board ID: " << boardId << std::endl;
        //std::cout << "           TTT: " << caenFrag->header.triggerTimeTag << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;

        // a +- 128-ns buffer time should be sufficient
        // in a TPCReadout "interval"
        double TPCReadoutLow  = timestamp - fTPCReadoutBufferLow;
        double TPCReadoutHigh = timestamp + fTPCReadoutBufferHigh;

        TPCReadout.push_back(std::make_pair(TPCReadoutLow, TPCReadoutHigh));
      }

      // these TPCReadout "intervals" are used for counting
      // the number of TPC readouts
      TPCReadout = this->IntervalsSelfMerge(TPCReadout);
      std::map< size_t, std::vector< unsigned int > > TPCReadoutCAENBoardIDCountMap;

      for (size_t j = 0; j < NumberCaenBlocks; ++j) {

        std::pair< double, const CAENFragment * > caenBlock = Collection.caenBlocks[j];
        double const& timestamp = caenBlock.first;
        const CAENFragment * caenFrag = caenBlock.second;
        unsigned int boardId = caenFrag->header.boardId;

        if (boardId < 8) {
          for (size_t k = 0; k < TPCReadout.size(); ++k) {
            if ((timestamp >= TPCReadout[k].first) and (timestamp <= TPCReadout[k].second)) {
              TPCReadoutCAENBoardIDCountMap[k].push_back(boardId);
            }
          }
        }

      }

      // count the number of TPC readouts
      for (std::map< size_t, std::vector< unsigned int > >::const_iterator
           iter = TPCReadoutCAENBoardIDCountMap.begin();
           iter != TPCReadoutCAENBoardIDCountMap.end();
           ++iter) {
        // require that boards 0, 1, 2, 3, 4, 5, 6, and 7
        // are present in a TPCReadout "interval" before
        // incrementing fNumberTPCReadouts
        if ((std::find(iter->second.begin(), iter->second.end(), 0) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 1) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 2) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 3) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 4) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 5) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 6) != iter->second.end()) and
            (std::find(iter->second.begin(), iter->second.end(), 7) != iter->second.end())) {
          NumberTPCReadouts += 1;
        }
      }

      Collection.numberTPCReadouts = NumberTPCReadouts;

      for (size_t j = 0; j < NumberTdcBlocks; ++j) {

        std::pair< double, const std::vector<TDCFragment::TdcEventData> * > tdcBlock = Collection.tdcBlocks[j];
        double const& timestamp = tdcBlock.first;
        //const std::vector<TDCFragment::TdcEventData> * tdcEvents = tdcBlock.second;

        std::cout << "    TDC block: " << j << std::endl;
        std::cout << "      Timestamp: " << timestamp << std::endl;
        //std::cout << "      TDC events: " << tdcEvents->size() << std::endl;
      }
    }

    return Collections;
  }

  //-----------------------------------------------------------------------
  void SlicerAlg::hello_world()
  {
    mf::LogVerbatim("SlicerAlg")
        << "\n///////////////////////////////////////////////////////////\n"
        << "Hello, World!\n"
        << "///////////////////////////////////////////////////////////\n";

    this->hello_kitty();

    return;
  }

  //-----------------------------------------------------------------------
  void SlicerAlg::hello_kitty()
  {
    return;
  }

}
