//////////////////////////////////////////////////////////////
// Name:      EventBuilderAlg.h
// Date:      20 August 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////
// Event slicing algorithm for LArIAT.
//////////////////////////////////////////////////////////////

#ifndef EVENTBUILDERALG_H
#define EVENTBUILDERALG_H

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// LArIATFragments includes
#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/BERNFragment.h"
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/TDCFragment.h"
#include "LArIATFragments/TriggerFragment.h"
#include "LArIATFragments/V1495Fragment.h"
//#include "LArIATFragments/WUTFragment.h"
//#include "LArIATFragments/LARASICFragment.h"
//#include "LArIATFragments/ReadoutError.h"

// LArIATSoft includes
#include "RawDataUtilities/ClockCorrectionAlg.h"

// C++ includes
#include <utility>
#include <vector>

namespace rdu {

  class EventBuilderAlg {

   public:

    // constructor
    EventBuilderAlg(fhicl::ParameterSet const& pset);

    // destructor
    ~EventBuilderAlg();

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset);

    // configure interval pre-acquisition and acquisition window lengths
    void Configure(double V1740PreAcquisitionWindow,
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
                   double TDCAcquisitionWindow);

    // generate a vector of intervals from a vector of data blocks
    std::vector< std::pair< double, double> > CreateIntervals(std::vector< rdu::DataBlock > const& DataBlocks,
                                                              unsigned int                  const& DeviceID,
                                                              double                        const& PreAcquisitionWindow,
                                                              double                        const& PostAcquisitionWindow,
                                                              double                        const& AcquisitionWindow);

    // merge overlapping intervals in a vector of intervals
    std::vector< std::pair< double, double > > IntervalsSelfMerge(std::vector< std::pair< double, double > > const& Intervals);

    // merge overlapping intervals between two vectors of intervals
    std::vector< std::pair< double, double > > MergeIntervals(std::vector< std::pair< double, double > > const& IntervalsA,
                                                              std::vector< std::pair< double, double > > const& IntervalsB);

    void ReadV1495Fragments(const LariatFragment * data);
    void ReadTriggerFragments(const LariatFragment * data);
    void ReadBERNFragments(const LariatFragment * data);

    // get data blocks with corrected timestamps
    std::vector< rdu::DataBlock > GetDataBlocks(const LariatFragment * data);

    // build
    std::vector< rdu::DataBlockCollection > Build(const LariatFragment * data);

    // this method is used for testing porpoises
    void hello_world();

   private:

    // this method is used for testing dolphins
    void hello_kitty();

    // parameters for generating intervals
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

    // parameters for counting the number of TPC readouts in an interval
    double fTPCReadoutBufferLow;
    double fTPCReadoutBufferHigh;

    // flag for sorting data block collections by their CAEN event counter
    bool fSortCAENEventCounter;

    // clock correction algorithm
    rdu::ClockCorrectionAlg fClockCorrectionAlg;

  };

}

#endif
