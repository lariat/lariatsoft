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

// ART Framework Includes
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "artdaq-core/Data/Fragment.hh"

#include "RawDataUtilities/FragmentUtility.h"

namespace rdu{

  //----------------------------------------------------------------------------
  // The triggerModuleLabel is the label of the module creating both the triggers
  // and the associations to the triggers
  FragmentUtility::FragmentUtility(art::Event  const& evt,
				   std::string const& daqModuleLabel,
				   std::string const& daqInstanceLabel)
    : fLariatFragment(nullptr)
  {

  art::Handle< std::vector<artdaq::Fragment> > fragments;
  evt.getByLabel(daqModuleLabel, daqInstanceLabel, fragments);

  if ( !fragments.isValid() )
      throw cet::exception("FragmentUtility") << "artdaq::Fragment handle is not valid, bail";
  if ( fragments->size() != 1 )
      throw cet::exception("FragmentUtility") << "artdaq::Fragment handle contains more than one fragment, bail";

  // get the fragments we are interested in
  const auto& frag((*fragments)[0]);

  const char * bytePtr = reinterpret_cast<const char *> (&*frag.dataBegin());
  fLariatFragment = new LariatFragment((char *) bytePtr, frag.dataSize() * sizeof(unsigned long long));
  LOG_VERBATIM("FragmentToDigit") << "Have data fragment "
				  << frag.dataSize() * sizeof(unsigned long long);
  fLariatFragment->print();
  fLariatFragment->printSpillTrailer();

  LariatFragment::SpillTrailer const& spillTrailer = fLariatFragment->spillTrailer;

  LOG_VERBATIM("FragmentUtility") << "evt.run(): "               << evt.run()   
				  << "; evt.subRun(): " 	 << evt.subRun()
				  << "; evt.event(): "  	 << evt.event() 
				  << "; evt.time().timeLow(): "  << evt.time().timeLow()
				  << "; evt.time().timeHigh(): " << evt.time().timeHigh()
				  << "\nrunNumber: "             << spillTrailer.runNumber  
				  << "; spillNumber: " 		 << spillTrailer.spillNumber
				  << "; timeStamp: "   		 << spillTrailer.timeStamp; 
    return;
  }

  //----------------------------------------------------------------------------
  FragmentUtility::~FragmentUtility()
  {
    delete fLariatFragment;
  }



} // end namespace
