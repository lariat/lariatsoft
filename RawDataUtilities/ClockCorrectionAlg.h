//////////////////////////////////////////////////////////////
// Name:      ClockCorrectionAlg.h
// Date:      15 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Clock correction algorithm for LArIAT. This corrects the
// timestamps of data blocks from the CAEN digitizers and
// wire chambers.
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

enum {
  TDC_DEVICE_ID = 32,
};

namespace rdu {

  struct DataBlock {

    unsigned int deviceId;
    double timestamp;
    double correctedTimestamp;

    // there should only be one data block in only one of these vectors
    std::vector< const CAENFragment * > caenBlocks;
    std::vector< const std::vector< TDCFragment::TdcEventData > * > tdcBlocks;
    //std::vector< const V1494Fragment* > v1495Blocks;
    //std::vector< const WUTFragment* > wutBlocks;
    //std::vector< const LARASICFragment* > larasicBlocks;
    //std::vector< const ReadoutError* > errorBlocks;

    void clear()
    {
      timestamp = 0;
      correctedTimestamp = 0;
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
  //struct Slice {

    double timestamp;
    std::vector< const CAENFragment * > caenBlocks;
    std::vector< const std::vector<TDCFragment::TdcEventData> * > tdcBlocks;
    //std::vector< const V1494Fragment * > v1495Blocks;
    //std::vector< const WUTFragment * > wutBlocks;
    //std::vector< const LARASICFragment * > larasicBlocks;
    //std::vector< const ReadoutError * > errorBlocks;

    void clear()
    {
      timestamp = 0;
      caenBlocks.clear();
      tdcBlocks.clear();
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
    void GetClockCorrectionParameters(std::map< unsigned int, std::vector< double > >  const& TimeStampMap,
                                      std::map< unsigned int, std::pair< double, double > > & ClockCorrectionParameters);

    // return a vector of data block collections
    std::vector< DataBlockCollection > GroupCollections(const LariatFragment * data);

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
    void UnweightedLinearFit(std::vector< std::pair< double, double> > const& Data,
                             double                                         & Slope,
                             double                                         & Intercept);

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

    // restrict the linear fit for RANSAC
    double fTimeStampDifferenceThreshold;
    double fSampleSlopeCutLower;
    double fSampleSlopeCutUpper;

    // store device IDs with at least one data block
    std::vector<unsigned int> fDeviceId;

  };

}

#endif
