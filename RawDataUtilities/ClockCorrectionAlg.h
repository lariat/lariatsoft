//////////////////////////////////////////////////////////////
// Name:      ClockCorrectionAlg.h
// Date:      15 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
//
// Clock correction algorithm for LArIAT. This corrects the
// timestamps of data blocks from the CAEN digitizers and
// wire chambers.
//
// The algorithm "corrects" the timestamps of all the devices
// accordingly to a reference device. Given two devices, A
// and B, each timestamp of device A is plotted against each
// timestamp of device B.
//
//////////////////////////////////////////////////////////////
//
// RANSAC algorithm:
//
//   1. Randomly choose a minimal subset of data
//   2. Use this subset to estimate the parameters of the
//      model
//   3. Compute the number of inliers for this model
//   4. Repeat steps 1-3 a fixed number of times
//   5. Re-estimate the model using inliers from the best
//      fit
//
//////////////////////////////////////////////////////////////
//
// NOTE:
// -----
// This algorithm does not care about the absolute time of
// when a trigger signal is received by a device. This means
// that the time difference between the fast and delayed
// triggers, along with the pre-/post-trigger windows of the
// CAEN V1740 and V1751 acquisition windows, does not affect
// the clock-drift correction.
//
//////////////////////////////////////////////////////////////

#ifndef CLOCKCORRECTIONALG_H
#define CLOCKCORRECTIONALG_H

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// LArIATFragments includes
#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/TDCFragment.h"
//#include "LArIATFragments/V1495Fragment.h"
//#include "LArIATFragments/WUTFragment.h"
//#include "LArIATFragments/LARASICFragment.h"
//#include "LArIATFragments/ReadoutError.h"

// LArIATSoft includes
#include "RawDataUtilities/FragmentUtility.h"

// C++ includes
#include <utility>
#include <vector>

namespace rdu {

  enum {
    TDC_DEVICE_ID = 32,
  };

  struct DataBlock {

    unsigned int deviceId = 999;
    double timestamp = -1;
    double correctedTimestamp = -1;

    // there should only be one data block in only one of these vectors
    std::vector< const CAENFragment * > caenBlocks;
    std::vector< const std::vector< TDCFragment::TdcEventData > * > tdcBlocks;
    //std::vector< const V1494Fragment * > v1495Blocks;
    //std::vector< const WUTFragment * > wutBlocks;
    //std::vector< const LARASICFragment * > larasicBlocks;
    //std::vector< const ReadoutError * > errorBlocks;

    void clear()
    {
      timestamp = -1;
      correctedTimestamp = -1;
      deviceId = 999;
      caenBlocks.clear();
      tdcBlocks.clear();
      //v1495Blocks.clear();
      //wutBlocks.clear();
      //larasicBlocks.clear();
      //errorBlocks.clear();
    }

  };

  struct DataBlockCollection {

    std::pair< double, double > interval;
    size_t numberTPCReadouts;
    int caenEventCounter = -1;

    // vector of pairs where the first is the corrected timestamp
    // and the second is a pointer to the data block
    //std::vector< std::pair< double, const CAENFragment * > > caenBlocks;
    //std::vector< std::pair< double, const std::vector<TDCFragment::TdcEventData> * > > tdcBlocks;
    //std::vector< const CAENFragment * > caenBlocks;
    //std::vector< const std::vector<TDCFragment::TdcEventData> * > tdcBlocks;
    //std::vector< const V1494Fragment * > v1495Blocks;
    //std::vector< const WUTFragment * > wutBlocks;
    //std::vector< const LARASICFragment * > larasicBlocks;
    //std::vector< const ReadoutError * > errorBlocks;

    // vectors of data blocks
    std::vector< CAENFragment > caenBlocks;
    std::vector< std::vector<TDCFragment::TdcEventData> > tdcBlocks;

    // vectors of timestamps of data blocks
    std::vector< double > caenBlockTimeStamps;
    std::vector< double > tdcBlockTimeStamps;

    void clear()
    {
      interval = std::make_pair(0.0, 0.0);
      numberTPCReadouts = 0;
      caenEventCounter = -1;
      caenBlocks.clear();
      tdcBlocks.clear();
      caenBlockTimeStamps.clear();
      tdcBlockTimeStamps.clear();
      //v1495Blocks.clear();
      //wutBlocks.clear();
      //larasicBlocks.clear();
      //errorBlocks.clear();
    }

  };

  class ClockCorrectionAlg {

   public:

    // constructor
    ClockCorrectionAlg(fhicl::ParameterSet const& pset);

    // destructor
    ~ClockCorrectionAlg();

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset);

    // RANSAC with a linear model
    void LinearRANSAC(std::vector< std::pair<double, double> > const& Data,
                      size_t                                   const& MinSamples,
                      double                                   const& ResidualThreshold,
                      size_t                                   const& MaxTrials,
                      size_t                                   const& StopSampleNumber,
                      double                                   const& StopResidualsSum,
                      double                                   const& StopProbability,
                      double                                        & Slope,
                      double                                        & Intercept);

    // get vector of DataBlocks and map of timestamps
    void GetDataBlocksTimeStampMap(const LariatFragment                            * data,
                                   std::vector< DataBlock >                        & DataBlocks,
                                   std::map< unsigned int, std::vector< double > > & TimeStampMap);

    // get clock correction parameters
    void GetClockCorrectionParameters(std::map< unsigned int, std::vector< double > >       const& TimeStampMap,
                                      std::map< unsigned int, std::pair< double, double > >      & ClockCorrectionParameters,
                                      unsigned int                                               & ReferenceClockDeviceID);

    // this method is used for testing porpoises
    void hello_world();

   private:

    // this method is used for testing dolphins
    void hello_kitty();

    // maximum values of different types
    double fMaxDouble;
    size_t fMaxSize_T;

    ///////////////////////////////////////////////////////////////////////////
    // Begin RANSAC helper functions
    ///////////////////////////////////////////////////////////////////////////

    // unweighted linear fit
    void UnweightedLinearFit(std::vector< std::pair<double, double> > const& Data,
                             double                                        & Slope,
                             double                                        & Intercept);

    // get residuals and residuals sum
    void GetLinearResiduals(std::vector< std::pair<double, double> > const& Data,
                            double                                   const& Slope,
                            double                                   const& Intercept,
                            std::vector<double>                           & Residuals,
                            double                                        & ResidualsSum,
                            bool                                     const& Orthogonal);

    // determine number trials such that at least one outlier-free subset is
    // sampled for the given inlier/outlier ratio
    size_t DynamicMaxTrials(size_t const& NumberInliers,
                            size_t const& NumberSamples,
                            size_t const& MinSamples,
                            double const& Probability);

    ///////////////////////////////////////////////////////////////////////////
    // End RANSAC helper functions
    ///////////////////////////////////////////////////////////////////////////

    // RANSAC parameters
    size_t fMinSamples;
    double fResidualThreshold;
    size_t fMaxTrials;
    size_t fStopSampleNumber;
    double fStopResidualsSum;
    double fStopProbability;
    size_t fMaxContinuations;
    size_t fMinNumberTimeStamps;
    size_t fMinNumberDataPoints;

    // restrict the linear fit for RANSAC
    double fTimeStampDifferenceThreshold;  // microseconds
    double fSampleSlopeCutLower;
    double fSampleSlopeCutUpper;

    // store device IDs with at least one data block
    std::vector<unsigned int> fDeviceId;

  };

}

#endif
