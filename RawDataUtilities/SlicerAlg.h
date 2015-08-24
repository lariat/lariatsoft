//////////////////////////////////////////////////////////////
// Name:      SlicerAlg.h
// Date:      20 August 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Event slicing algorithm for LArIAT.
//////////////////////////////////////////////////////////////

#ifndef SLICERALG_H
#define SLICERALG_H

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
#include "RawDataUtilities/ClockCorrectionAlg.h"

// C++ includes
#include <utility>
#include <vector>

namespace rdu {

  class SlicerAlg {

   public:

    // constructor
    SlicerAlg(fhicl::ParameterSet const& pset);

    // destructor
    ~SlicerAlg();

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset);

    // configure interval pre-/post-trigger window lengths
    void Configure(double V1740PreTriggerWindow,
                   double V1740PostTriggerWindow,
                   double V1740BPreTriggerWindow,
                   double V1740BPostTriggerWindow,
                   double V1751PreTriggerWindow,
                   double V1751PostTriggerWindow,
                   double TDCPreTriggerWindow,
                   double TDCPostTriggerWindow);

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

    // slice
    //void Slice(const LariatFragment * data);
    std::vector< rdu::DataBlockCollection > Slice(const LariatFragment * data);

    // this method is used for testing porpoises
    void hello_world();

   private:

    // this method is used for testing dolphins
    void hello_kitty();

    // parameters for generating intervals
    double fV1740PreTriggerWindow;
    double fV1740PostTriggerWindow;
    double fV1740BPreTriggerWindow;
    double fV1740BPostTriggerWindow;
    double fV1751PreTriggerWindow;
    double fV1751PostTriggerWindow;
    double fTDCPreTriggerWindow;
    double fTDCPostTriggerWindow;

    // clock correction algorithm
    rdu::ClockCorrectionAlg fClockCorrectionAlg;

  };

}

#endif
