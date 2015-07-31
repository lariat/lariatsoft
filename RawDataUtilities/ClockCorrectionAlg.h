//////////////////////////////////////////////////////////////
// Name:      ClockCorrectionAlg.h
// Date:      15 July 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Relativistic clock correction algorithm for LArIAT.
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
#include <vector>

namespace rdu {

  //struct DataBlock {

  //  int deviceId;
  //  double timestamp;
  //  double correctedTimestamp;

  //  // there should only be one data block in only one of these vectors
  //  std::vector< const CAENFragment * > caenBlocks;
  //  std::vector< const std::vector< TDCFragment::TdcEventData > * > tdcBlocks;
  //  //std::vector< const V1494Fragment* > v1495Blocks;
  //  //std::vector< const WUTFragment* > wutBlocks;
  //  //std::vector< const LARASICFragment* > larasicBlocks;
  //  //std::vector< const ReadoutError* > errorBlocks;

  //  void clear()
  //  {
  //    timestamp = 0;
  //    correctedTimestamp = 0;
  //    deviceId = -1;
  //    caenBlocks.clear();
  //    tdcBlocks.clear();
  //    //v1495Blocks.clear();
  //    //wutBlocks.clear();
  //    //larasicBlocks.clear();
  //    //errorBlocks.clear();
  //  }

  //};

  struct DataBlockCollection {
  //struct Slice {

    double timestamp;
    std::vector< const CAENFragment* > caenBlocks;
    std::vector< const std::vector<TDCFragment::TdcEventData>* > tdcBlocks;
    //std::vector< const V1494Fragment* > v1495Blocks;
    //std::vector< const WUTFragment* > wutBlocks;
    //std::vector< const LARASICFragment* > larasicBlocks;
    //std::vector< const ReadoutError* > errorBlocks;

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

    //typedef std::map<size_t, size_t> inlier_number_map;
    //typedef std::map<size_t, double> residuals_sum_map;
    //typedef std::map<size_t, double> slope_map;
    //typedef std::map<size_t, double> intercept_map;
    //typedef std::map<size_t, std::vector< std::pair< double, double > > inliers_map;

    // BestInlierNumber = SampleInlierNumber;
    // BestInlierResidualsSum = SampleResidualsSum;
    // BestSlope = SampleSlope;
    // BestIntercept = SampleIntercept;
    // BestInliers = SampleInliers;

    //  << "\nNumber of trials:            " << TrialNumber
    //  << "\nNumber of convergence steps: " << ConvergenceSteps
    //  << "\nTrial number of convergence: " << ConvergenceTrialNumber
    //  << "\nBest slope:                  " << BestSlope
    //  << "\nBest intercept:              " << BestIntercept
    //  << "\nBestInlierNumber:            " << BestInlierNumber
    //  << "\nBest inlier residuals sum:   " << BestInlierResidualsSum

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

    // return a map of timestamps
    //std::map< unsigned int, std::map< unsigned int, double > > GetTimeStampMap(const LariatFragment * data);
    std::map< unsigned int, std::vector< double > > GetTimeStampMap(const LariatFragment * data);

    // apply clock correction
    void ClockCorrection(std::map< unsigned int, std::vector< double > > const& TimeStampMap);

    // return a vector of data block collections
    std::vector< DataBlockCollection > GroupCollections(const LariatFragment * data);

    // return a vector of Slice structs
    //std::vector< Slice > slice(LariatFragment * data);

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
