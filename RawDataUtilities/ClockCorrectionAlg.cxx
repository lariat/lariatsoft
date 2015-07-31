//////////////////////////////////////////////////////////////
// Name:      ClockCorrectionAlg.cxx
// Date:      15 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Clock correction algorithm for LArIAT. This corrects the
// timestamps of data blocks from the CAEN digitizers and
// wire chambers.
//////////////////////////////////////////////////////////////

// Class include
#include "RawDataUtilities/ClockCorrectionAlg.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <limits>

namespace rdu {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  ClockCorrectionAlg::ClockCorrectionAlg(fhicl::ParameterSet const& pset) {

    // since we can't use infinity here, this will have to do
    fMaxDouble = std::numeric_limits<double>::max();
    fMaxSize_T = std::numeric_limits<size_t>::max();

    // read in parameters from .fcl files
    this->reconfigure(pset);

  }

  //-----------------------------------------------------------------------
  // destructor
  ClockCorrectionAlg::~ClockCorrectionAlg() {}

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMinSamples                   = pset.get< size_t >("MinSamples",                    2);
    fResidualThreshold            = pset.get< double >("ResidualThreshold",             1);
    fMaxTrials                    = pset.get< size_t >("MaxTrials",                     10000);
    fStopSampleNumber             = pset.get< size_t >("StopSampleNumber",              fMaxSize_T);
    fStopResidualsSum             = pset.get< double >("StopResidualsSum",              0);
    fStopProbability              = pset.get< double >("StopProbability",               1);
    fTimeStampDifferenceThreshold = pset.get< double >("TimeStampDifferenceThresdhold", 1e6);
    fSampleSlopeCutLower          = pset.get< double >("SampleSlopeCutLower",           0.5);
    fSampleSlopeCutUpper          = pset.get< double >("SampleSlopeCutUpper",           1.5);
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::UnweightedLinearFit(std::vector< std::pair< double, double> > const& Data,
                                               double                                         & Slope,
                                               double                                         & Intercept)
  {
    const size_t NumberDataPoints = Data.size();

    if (NumberDataPoints < 2) {
      throw cet::exception("ClockCorrectionAlg") << "Not enough data points to perform fit: "
                                                 << NumberDataPoints;
    }

    double sum_x = 0;
    double sum_y = 0;
    double sum_xy = 0;
    double sum_xx = 0;
    double sum_yy = 0;

    for (size_t i = 0; i < NumberDataPoints; ++i) {
      double x = Data.at(i).first;
      double y = Data.at(i).second;
      sum_x += x;
      sum_y += y;
      sum_xy += x*y;
      sum_xx += x*x;
      sum_yy += y*y;
    }

    mf::LogDebug("ClockCorrectionAlg") << "sum_x: "    << sum_x
                                       << "\nsum_y: "  << sum_y
                                       << "\nsum_xy: " << sum_xy
                                       << "\nsum_xx: " << sum_xx
                                       << "\nsum_yy: " << sum_yy;

    const double d = NumberDataPoints * sum_xx - sum_x * sum_x;
    Slope = (NumberDataPoints * sum_xy - sum_x * sum_y) / d;
    Intercept = (sum_y * sum_xx - sum_x * sum_xy) / d;

    return;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::GetLinearResiduals(std::vector< std::pair<double, double> > const& Data,
                                              double                                   const& Slope,
                                              double                                   const& Intercept,
                                              std::vector<double>                           & Residuals,
                                              double                                        & ResidualsSum,
                                              bool                                     const& Orthogonal)
  {
    ResidualsSum = 0;

    for (size_t datum_idx = 0; datum_idx < Data.size(); ++datum_idx) {
      double x = Data.at(datum_idx).first;
      double y = Data.at(datum_idx).second;
      double Residual;
      if (Orthogonal) {
        // uses the shortest distance from data point to line
        Residual = (y - Slope * x - Intercept) / std::sqrt(1 + Slope * Slope);
      }
      else {
        Residual = y - Intercept - Slope * x;
      }
      Residuals.push_back(Residual);
      ResidualsSum += Residual * Residual;
    }

    if (Data.size() != Residuals.size()) {
      throw cet::exception("ClockCorrectionAlg")
        << "Number of residuals is not equal to number of data points!"
        << " Number of residuals: " << Residuals.size()
        << " Number of data points: " << Data.size();
    }

    return;
  }

  //-----------------------------------------------------------------------
  size_t ClockCorrectionAlg::DynamicMaxTrials(size_t const& NumberInliers,
                                              size_t const& NumberSamples,
                                              size_t const& MinSamples,
                                              double const& Probability)
  {
    // Determine number trials such that at least one outlier-free subset is
    // sampled for the given inlier/outlier ratio.

    if (NumberInliers == 0) return fMaxSize_T;

    double Numerator = 1.0 - Probability;

    if (Numerator == 0) return fMaxSize_T;

    double InlierRatio = NumberInliers / (double) NumberSamples;

    double Denominator = 1.0 - std::pow((double) InlierRatio, MinSamples);

    if (Denominator == 0) { return 1; }
    else if (Denominator == 1) { return fMaxSize_T; };

    Numerator = log(Numerator);
    Denominator = log(Denominator);

    if (Denominator == 0) return 0;

    return (size_t) std::ceil(Numerator / Denominator);
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::LinearRANSAC(std::vector< std::pair<double, double> > const& Data,
                                        size_t                                   const& MinSamples,
                                        double                                   const& ResidualThreshold,
                                        size_t                                   const& MaxTrials,
                                        size_t                                   const& StopSampleNumber,
                                        double                                   const& StopResidualsSum,
                                        double                                   const& StopProbability,
                                        double                                        & Slope,
                                        double                                        & Intercept)
  {
    // random seed for shuffle
    std::srand(std::time(0));

    // get number of data points
    size_t NumberDataPoints = Data.size();

    // initialize the parameters of the best linear model
    size_t BestInlierNumber = 0;
    double BestInlierResidualsSum = fMaxDouble;
    double BestSlope = 0;
    double BestIntercept = 0;
    std::vector< std::pair<double, double> > BestInliers;

    // get vector of datum indices for shuffling
    std::vector<size_t> DatumIndices;
    for (size_t idx = 0; idx < NumberDataPoints; ++idx) {
      DatumIndices.push_back(idx);
    }

    // number of steps taken to convergence
    size_t ConvergenceSteps = 0;
    // trial number of convergence
    size_t ConvergenceTrialNumber = 0;

    // initialize iteration bookkeeper
    size_t TrialNumber = 0;

    // commence iterations
    while (TrialNumber < MaxTrials) {

      // copy vector of datum indices to shuffle
      std::vector<size_t> RandomIndices(DatumIndices);

      // shuffle the datum indices
      std::random_shuffle(RandomIndices.begin(), RandomIndices.end());

      // initialize vector of sample data points
      std::vector< std::pair< double, double > > Sample;

      // get sample of size MinSamples from the shuffled datum indices
      for (size_t rand_idx = 0; rand_idx < MinSamples; ++rand_idx) {
        //std::cout << RandomIndices.at(rand_idx) << std::endl;
        size_t DatumIndex = RandomIndices.at(rand_idx);
        Sample.push_back(Data.at(DatumIndex));
      }

      // get sample size
      //size_t NumberSamples = Sample.size();

      // initialize variables for linear fit of sample data
      double SampleSlope = 0;
      double SampleIntercept = 0;
      std::vector<double> SampleResiduals;
      double SampleResidualsSum = fMaxDouble;
      std::vector< std::pair<double, double> > SampleInliers;

      try {
        // attempt linear fit
        this->UnweightedLinearFit(Sample, SampleSlope, SampleIntercept);
        // get residuals
        this->GetLinearResiduals(Data, SampleSlope, SampleIntercept, SampleResiduals, SampleResidualsSum, true);
      }
      catch (cet::exception &e) {
        // if fit fails, continue
        continue;
      }

      // only allow real, finite values
      if (std::isnan(SampleSlope) or std::isinf(SampleSlope)
          or std::isnan(SampleIntercept) or std::isinf(SampleIntercept)
          or std::isnan(SampleResidualsSum) or std::isinf(SampleResidualsSum)) {
        continue;
      }

      // only allow slope in range (fSampleSlopeCutLower, fSampleSlopeCutUpper)
      if (SampleSlope < fSampleSlopeCutLower
          or SampleSlope > fSampleSlopeCutUpper) continue;

      // get sample inliers that is within ResidualThreshold of model
      for (size_t datum_idx = 0; datum_idx < NumberDataPoints; ++datum_idx) {
        if (std::abs(SampleResiduals.at(datum_idx)) < ResidualThreshold) {
          SampleInliers.push_back(Data.at(datum_idx));
        }
      }

      // get number of sample inliers
      size_t SampleInlierNumber = SampleInliers.size();

      mf::LogVerbatim("ClockCorrectionAlg")
        << "/////////////////////////////////////////////"
        << "\nTrial:                       " << TrialNumber
        << "\nSample slope:                " << SampleSlope
        << "\nSample intercept:            " << SampleIntercept
        << "\nSample inlier residuals sum: " << SampleResidualsSum
        << "\nNumber of data points:       " << NumberDataPoints
        << "\nNumber of sample residuals:  " << SampleResiduals.size()
        << "\nNumber of sample inliers:    " << SampleInlierNumber;

      // only allow real, finite values
      if (std::isnan(SampleInlierNumber)
          or std::isinf(SampleInlierNumber)) continue;

      // choose as new best model if number of inliers is maximal
      if (SampleInlierNumber > BestInlierNumber
          or (SampleInlierNumber == BestInlierNumber
              and SampleResidualsSum < BestInlierResidualsSum)) {

        BestInlierNumber = SampleInlierNumber;
        BestInlierResidualsSum = SampleResidualsSum;
        BestSlope = SampleSlope;
        BestIntercept = SampleIntercept;
        BestInliers = SampleInliers;

        ConvergenceTrialNumber = TrialNumber;
        ++ConvergenceSteps;

        if ((BestInlierNumber >= StopSampleNumber)
            or (BestInlierResidualsSum <= StopResidualsSum)
            or (TrialNumber >= this->DynamicMaxTrials(BestInlierNumber, NumberDataPoints, MinSamples, StopProbability))) {
          mf::LogInfo("ClockCorrectionAlg") << "DynamicMaxTrials reached! Stopping after " << TrialNumber << " trials.";
          break;
        }
      } // if sample model is current best model

      // increment trial number
      ++TrialNumber;
    } // end while loop

    // get best model after trials are done
    Slope = BestSlope;
    Intercept = BestIntercept;

    mf::LogVerbatim("ClockCorrectionAlg")
      << "\n/////////////////////////////////////////////"
      << "\nNumber of trials:            " << TrialNumber
      << "\nNumber of convergence steps: " << ConvergenceSteps
      << "\nTrial number of convergence: " << ConvergenceTrialNumber
      << "\nBest slope:                  " << BestSlope
      << "\nBest intercept:              " << BestIntercept
      << "\nBestInlierNumber:            " << BestInlierNumber
      << "\nBest inlier residuals sum:   " << BestInlierResidualsSum
      << "\n/////////////////////////////////////////////"
      << "\n1 - (best slope):            " << 1 - BestSlope
      << "\n/////////////////////////////////////////////";

    return;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::ClockCorrection(std::map< unsigned int, std::vector< double > > const& TimeStampMap)
  {

    // Container for clock correction parameters. The key is the device
    // ID, the mapped value is a pair. The first value of the pair is the
    // slope and the second value is the intercept.
    std::map< unsigned int, std::pair< double, double > > ClockCorrectionParameters;

    // container of device IDs that have at least one data block
    std::vector<unsigned int> DeviceIDs;

    // device ID of reference clock
    unsigned int ReferenceClockDeviceID = 999;

    // get device IDs that have at least one data block and get the
    // device ID of the reference clock
    for (std::map< unsigned int, std::vector< double > >::const_iterator
        iter = TimeStampMap.begin(); iter != TimeStampMap.end(); ++iter) {

      unsigned int DeviceID = iter->first;
      size_t NumberDataBlocks = iter->second.size();

      DeviceIDs.push_back(DeviceID);

      if (NumberDataBlocks > 0) {
        if (DeviceID == 8) {
          ReferenceClockDeviceID = DeviceID;
        }
        else if (DeviceID == 0 and ReferenceClockDeviceID != 8) {
          ReferenceClockDeviceID = DeviceID;
        }
      }

    } // end loop over TimeStampMap

    // get reference timestamps
    std::vector<double> const& TimeStampsA = TimeStampMap.at(ReferenceClockDeviceID);

    // collect the data points
    for (size_t i = 0; i < DeviceIDs.size(); ++i) {

      unsigned int const& DeviceID = DeviceIDs[i];

      if (DeviceID == ReferenceClockDeviceID) continue;

      std::vector<double> const& TimeStampsB = TimeStampMap.at(DeviceID);

      std::vector< std::pair< double, double > > Data;

      for (size_t m = 0; m < TimeStampsA.size(); ++m) {
        for (size_t n = 0; n < TimeStampsB.size(); ++n) {

          if (std::abs(TimeStampsA.at(m) - TimeStampsB.at(n)) > fTimeStampDifferenceThreshold) continue;

          Data.push_back(std::make_pair(TimeStampsA.at(m), TimeStampsB.at(n)));
          mf::LogDebug("ClockCorrectionAlg")
              << "(" << TimeStampsA.at(m) << ", " << TimeStampsB.at(n) << ")";

        } // end loop over TimeStampsB
      } // end loop over TimeStampsA

      // clock correction parameters
      double Slope = 0;
      double Intercept = 0;

      // RANSAC
      this->LinearRANSAC(Data,
                         fMinSamples,
                         fResidualThreshold,
                         fMaxTrials,
                         fStopSampleNumber,
                         fStopResidualsSum,
                         fStopProbability,
                         Slope,
                         Intercept);

      // copy clock correction parameters to container
      ClockCorrectionParameters[DeviceID] = std::make_pair(Slope, Intercept);

      mf::LogInfo("ClockCorrectionAlg")
          << "Device IDs: (" << ReferenceClockDeviceID << ", " << DeviceID << ")";

    } // end loop over device IDs

    // clock correction parameters for reference clock
    ClockCorrectionParameters[ReferenceClockDeviceID] = std::make_pair(1, 0);

    return;
  }

  //-----------------------------------------------------------------------
  std::map< unsigned int, std::vector< double > > ClockCorrectionAlg::GetTimeStampMap(const LariatFragment * data)
  {

    //std::map< unsigned int, std::map< unsigned int, double > > TimeStampMap;
    std::map< unsigned int, std::vector< double > > TimeStampMap;

    std::map< unsigned int, unsigned int > numberCaenDataBlocks;

    // get number of CAEN fragments / data blocks
    const size_t numberCaenFrags = data->caenFrags.size();

    mf::LogInfo("ClockCorrectionAlg") << "Found "
                                      << numberCaenFrags
                                      << " CAEN fragments";

    if (numberCaenFrags > 0)
      mf::LogInfo("ClockCorrectionAlg") << "Looking at CAEN fragments...";

    for (size_t i = 0; i < numberCaenFrags; ++i) {

      // get CAEN fragments
      CAENFragment const& caenFrag = data->caenFrags[i];

      // get board ID of CAEN fragment
      unsigned int boardId = static_cast <unsigned int> (caenFrag.header.boardId);

      if (std::find(fDeviceId.begin(), fDeviceId.end(), boardId) == fDeviceId.end()) {
        fDeviceId.push_back(boardId);
      }

      // CAEN fragment index for this boardId
      unsigned int index = numberCaenDataBlocks[boardId];

      mf::LogDebug("ClockCorrectionAlg") << "CAENFragment: "
                                         << i
                                         << "; boardId "
                                         << boardId
                                         << " index: "
                                         << index;

      // each CAEN Trigger Time Tag count is 8 ns
      double timestamp = caenFrag.header.triggerTimeTag * 0.008;  // convert to microseconds

      //TimeStampMap[boardId][i] = timestamp;
      TimeStampMap[boardId].push_back(timestamp);

      // add count to number of CAEN data blocks map
      numberCaenDataBlocks[boardId] += 1;

    } // end CAENFragment loop

    for (std::map< unsigned int, unsigned int >::const_iterator
         iter = numberCaenDataBlocks.begin();
         iter != numberCaenDataBlocks.end();
         ++iter) {

      mf::LogDebug("ClockCorrectionAlg") << "board ID: "
                                         << iter->first
                                         << "; number of CAEN data blocks: "
                                         << iter->second;
    }

    unsigned int numberTdcDataBlocks;

    // get number of TDC fragments
    const size_t numberTdcFrags = data->tdcFrags.size();

    mf::LogInfo("ClockCorrectionAlg") << "Found "
                                      << numberTdcFrags
                                      << " TDC fragments";

    if (numberTdcFrags > 0)
      mf::LogInfo("ClockCorrectionAlg") << "Looking at TDC fragments...";

    for (size_t i = 0; i < numberTdcFrags; ++i) {

      // get TDC fragment
      TDCFragment const& tdcFrag = data->tdcFrags[i];

      // get TDC data blocks
      std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcEvents = tdcFrag.tdcEvents;

      // get number of TDC data blocks
      numberTdcDataBlocks = tdcEvents.size();

      for (size_t j = 0; j < tdcEvents.size(); ++j) {

        if (tdcFrag.controllerHeader.nTDCs != tdcEvents[j].size()) {
          mf::LogError("ClockCorrectionAlg") << "*** Fatal nTDCs mismatch: " << tdcEvents[j].size()
              << " != " << tdcFrag.controllerHeader.nTDCs<< " "<< j;
          continue;
        }

        if (std::find(fDeviceId.begin(), fDeviceId.end(), 32) == fDeviceId.end()) {
          fDeviceId.push_back(32);
        }

        mf::LogDebug("ClockCorrectionAlg") << "TDC event: " << j;

        ///////////////////////////////////////////////////////////////////////
        // Loop over TDCs to get a histogram of the TDC timestamps. The
        // timestamp with the most counts is taken as the TDC timestamp of
        // the TDC data block. We are doing this to work around the mismatches
        // of the TDC timestamps.
        ///////////////////////////////////////////////////////////////////////

        std::map<int, unsigned int> tdcTimeStampCounts;

        for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS; ++tdc_index) {
          TDCFragment::TdcEventData tdcEventData = tdcEvents[j].at(tdc_index);

          // each TDC Time Stamp count is 1/106.208 microseconds
          int tdcTimeStamp = static_cast <int> (tdcEventData.tdcEventHeader.tdcTimeStamp);

          tdcTimeStampCounts[tdcTimeStamp] += 1;

        } // end loop over TDCs

        unsigned int counts = 0;
        int tdcTimeStamp = 0;

        for (auto const& k : tdcTimeStampCounts) {
          if (k.second > counts) {
            tdcTimeStamp = k.first;
            counts = k.second;
          }
        }

        mf::LogDebug("ClockCorrectionAlg") << "tdcTimeStamp: " << tdcTimeStamp
                                           << "\ncounts: " << counts;

        ///////////////////////////////////////////////////////////////////////

        // each TDC Time Stamp count is 1/106.208 microseconds
        double timestamp = tdcTimeStamp / 106.208;  // convert to microseconds

        //TimeStampMap[32][j] = timestamp;
        TimeStampMap[32].push_back(timestamp);

      } // TDC data block loop
    } // TDCFragment loop

    mf::LogInfo("ClockCorrectionAlg") << "Number of TDC data blocks: "
                                      << numberTdcDataBlocks;

    mf::LogInfo("ClockCorrectionAlg")
        << "//////////////////////////////////////////"
        << "\n Device IDs with at least one data block"
        << "\n//////////////////////////////////////////";

    for (std::vector<unsigned int>::const_iterator
         iter = fDeviceId.begin(); iter != fDeviceId.end(); ++iter) {
      mf::LogInfo("ClockCorrectionAlg") << *iter;
    }

    mf::LogInfo("ClockCorrectionAlg")
        << "\n//////////////////////////////////////////";

    return TimeStampMap;
  }

  //-----------------------------------------------------------------------
  std::vector< DataBlockCollection > ClockCorrectionAlg::GroupCollections(const LariatFragment * data)
  {

    std::vector< DataBlockCollection > collections;

    //std::map< unsigned int, std::map< unsigned int, double > > TimeStampMap;
    std::map< unsigned int, std::vector< double > > TimeStampMap;

    TimeStampMap = this->GetTimeStampMap(data);

    //for (std::map< unsigned int, std::map< unsigned int, double> >::const_iterator
    for (std::map< unsigned int, std::vector< double> >::const_iterator
         iter = TimeStampMap.begin(); iter != TimeStampMap.end(); ++iter) {

      mf::LogVerbatim("ClockCorrectionAlg") << "Device ID: "
                                            << iter->first
                                            << "; number of data blocks: "
                                            << iter->second.size();
    }

    mf::LogVerbatim("ClockCorrectionAlg") << "blah";
    this->ClockCorrection(TimeStampMap);
    mf::LogVerbatim("ClockCorrectionAlg") << "blarg";

    return collections;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::hello_world()
  {

    mf::LogVerbatim("ClockCorrectionAlg")
        << "\n///////////////////////////////////////////////////////////\n"
        << "Hello, World!\n"
        << "///////////////////////////////////////////////////////////\n";

    this->hello_kitty();

    return;
  }

  //-----------------------------------------------------------------------
  void ClockCorrectionAlg::hello_kitty()
  {
    return;
  }

}
