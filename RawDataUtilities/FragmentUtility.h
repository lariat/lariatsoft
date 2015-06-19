////////////////////////////////////////////////////////////////////////////////
// Utility to grab fragments from the DAQ output
//
// Do not attempt to write this utility into the event record
//
// 
// FragmentUtility.h
//
// \author Brian Rebel brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef RAWDATAUTILITIES_FRAGMENTUTILITY_H
#define RAWDATAUTILITIES_FRAGMENTUTILITY_H

#include <iostream>

// ART Framework Includes
#include "art/Framework/Principal/Event.h"

// LArIATFragments
#include "LArIATFragments/LariatFragment.h"

namespace rdu{

  // An object to provide access to triggers and digits
  class FragmentUtility{

  public: 
    
    // The triggerModuleLabel is the label of the module creating both the triggers
    // and the associations to the triggers
    explicit FragmentUtility(art::Event  const& evt,
			     std::string const& daqModuleLabel,
			     std::string const& daqInstanceLabel);
    virtual  ~FragmentUtility();

    LariatFragment const& DAQFragment();

    std::vector<WUTFragment>     const& WUTFragments();
    std::vector<V1495Fragment>   const& V1495Fragments();
    std::vector<CAENFragment>    const& CAENFragments();
    std::vector<TDCFragment>     const& TDCFragments();
    std::vector<LARASICFragment> const& LARASICFragments();
    std::vector<ReadoutError>    const& ReadoutErrors();

  private:

    LariatFragment *fLariatFragment; ///< total fragment for the event record

  };   

}// end namespace

inline LariatFragment               const& rdu::FragmentUtility::DAQFragment()      { return *fLariatFragment;              }
inline std::vector<WUTFragment>     const& rdu::FragmentUtility::WUTFragments()     { return fLariatFragment->wutFrags;     }
inline std::vector<V1495Fragment>   const& rdu::FragmentUtility::V1495Fragments()   { return fLariatFragment->v1495Frags;   }
inline std::vector<CAENFragment>    const& rdu::FragmentUtility::CAENFragments()    { return fLariatFragment->caenFrags;    }
inline std::vector<TDCFragment>     const& rdu::FragmentUtility::TDCFragments()     { return fLariatFragment->tdcFrags;     }
inline std::vector<LARASICFragment> const& rdu::FragmentUtility::LARASICFragments() { return fLariatFragment->larasicFrags; }
inline std::vector<ReadoutError>    const& rdu::FragmentUtility::ReadoutErrors()    { return fLariatFragment->errorFrags;   }

#endif //RAWDATAUTILITIES_FRAGMENTUTILITY_H
 
