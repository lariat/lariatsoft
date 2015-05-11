////////////////////////////////////////////////////////////////////////////////
// Utility to encapsulate code that grabs digits (RawDigit/AuxDetDigit/OpDetPulse)
// for a single trigger and returns them to the user.  Provides methods to grab
// digits from each auxiliary detector individually.
//
// Do not attempt to write this utility into the event record
//
// Use the versions that return art::PtrVector if you need to write a new Association 
// using the given data type into the event record, otherwise use the std::vector versions
// 
// TriggerDigitUtility.h
//
// \author Brian Rebel brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef RAWDATAUTILITIES_TRIGGERDIGITUTILITY_H
#define RAWDATAUTILITIES_TRIGGERDIGITUTILITY_H

#include <iostream>

// ART Framework Includes
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Persistency/Common/Assns.h"

// LArSoft Includes
#include "RawData/RawDigit.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/TriggerData.h"

namespace rdu{

  // An object to provide access to triggers and digits
  class TriggerDigitUtility{

  public: 
    
    // The triggerModuleLabel is the label of the module creating both the triggers
    // and the associations to the triggers
    explicit TriggerDigitUtility(art::Event  const& evt,
				 std::string const& triggerModuleLabel);
    virtual  ~TriggerDigitUtility();

    // Methods for accessing Triggers
    size_t                                    NTriggers()                                        const;
    std::vector<const raw::Trigger*>          EventTriggers() 				      	 const;			  
    art::PtrVector<raw::Trigger>       const& EventTriggersPtr() 			      	 const;

    // Methods for accessing RawDigits
    std::vector<const raw::RawDigit*>         EventRawDigits()   		   	         const;			   	 
    art::PtrVector<raw::RawDigit>             EventRawDigitsPtr()   		   	         const;                  	 
    std::vector<const raw::RawDigit*>         TriggerRawDigits(size_t const& t) 	      	 const;			 
    art::PtrVector<raw::RawDigit>      const& TriggerRawDigitsPtr(size_t const& t) 	      	 const;
												       
    // Methods for accessing OpDetPulses							       
    std::vector<const raw::OpDetPulse*>       EventOpDetPulses()   	       		      	 const;			       
    art::PtrVector<raw::OpDetPulse>           EventOpDetPulsesPtr()   	       		         const;                               
    std::vector<const raw::OpDetPulse*>       TriggerOpDetPulses(size_t const& t) 	      	 const;			       
    art::PtrVector<raw::OpDetPulse>    const& TriggerOpDetPulsesPtr(size_t const& t) 	      	 const;
												       
    // Methods for accessing AuxDetDigits from the Muon Range Stack				       
    std::vector<const raw::AuxDetDigit*>      EventMuonRangeStackDigits()   		      	 const;			
    art::PtrVector<raw::AuxDetDigit>          EventMuonRangeStackDigitsPtr()   		         const;                   
    std::vector<const raw::AuxDetDigit*>      TriggerMuonRangeStackDigits   (size_t const& t) 	 const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerMuonRangeStackDigitsPtr(size_t const& t) 	 const;

    // Methods for accessing AuxDetDigits from the TOF
    std::vector<const raw::AuxDetDigit*>      EventUpStreamTOFDigits()   	      	         const;   				      
    std::vector<const raw::AuxDetDigit*>      EventDownStreamTOFDigits()        	         const;			      
    art::PtrVector<raw::AuxDetDigit>          EventUpStreamTOFDigitsPtr()          	         const;                    	     
    art::PtrVector<raw::AuxDetDigit>          EventDownStreamTOFDigitsPtr()                      const;                              
    std::vector<const raw::AuxDetDigit*>      TriggerUpStreamTOFDigits     (size_t const& t)  	 const;		      
    std::vector<const raw::AuxDetDigit*>      TriggerDownStreamTOFDigits   (size_t const& t)  	 const;		      
    art::PtrVector<raw::AuxDetDigit>   const& TriggerUpStreamTOFDigitsPtr  (size_t const& t)  	 const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerDownStreamTOFDigitsPtr(size_t const& t)  	 const;

    // Methods for accessing AuxDetDigits from the AeroGel Counters
    std::vector<const raw::AuxDetDigit*>      EventUpStreamAeroGelDigits()   		         const; 				 
    std::vector<const raw::AuxDetDigit*>      EventDownStreamAeroGelDigits()   		         const;			
    art::PtrVector<raw::AuxDetDigit>          EventUpStreamAeroGelDigitsPtr()   		 const;                  
    art::PtrVector<raw::AuxDetDigit>          EventDownStreamAeroGelDigitsPtr()   	         const;                  
    std::vector<const raw::AuxDetDigit*>      TriggerUpStreamAeroGelDigits     (size_t const& t) const;
    std::vector<const raw::AuxDetDigit*>      TriggerDownStreamAeroGelDigits   (size_t const& t) const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerUpStreamAeroGelDigitsPtr  (size_t const& t) const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerDownStreamAeroGelDigitsPtr(size_t const& t) const;

    // Methods for accessing AuxDetDigits from the MWPCs
    std::vector<const raw::AuxDetDigit*>      EventMWPC1Digits()   			         const;			       	     
    std::vector<const raw::AuxDetDigit*>      EventMWPC2Digits()   			    	 const;		     
    std::vector<const raw::AuxDetDigit*>      EventMWPC3Digits()   			    	 const;		     
    std::vector<const raw::AuxDetDigit*>      EventMWPC4Digits()   			    	 const;		     
    art::PtrVector<raw::AuxDetDigit>          EventMWPC1DigitsPtr()   		                 const;                    
    art::PtrVector<raw::AuxDetDigit>          EventMWPC2DigitsPtr()   			         const;                 
    art::PtrVector<raw::AuxDetDigit>          EventMWPC3DigitsPtr()   			         const;                 
    art::PtrVector<raw::AuxDetDigit>          EventMWPC4DigitsPtr()   			         const;                 
    std::vector<const raw::AuxDetDigit*>      TriggerMWPC1Digits   (size_t const& t) 	    	 const;		     
    std::vector<const raw::AuxDetDigit*>      TriggerMWPC2Digits   (size_t const& t) 	    	 const;		     
    std::vector<const raw::AuxDetDigit*>      TriggerMWPC3Digits   (size_t const& t) 	    	 const;		     
    std::vector<const raw::AuxDetDigit*>      TriggerMWPC4Digits   (size_t const& t) 	    	 const;		     
    art::PtrVector<raw::AuxDetDigit>   const& TriggerMWPC1DigitsPtr(size_t const& t) 	    	 const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerMWPC2DigitsPtr(size_t const& t) 	    	 const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerMWPC3DigitsPtr(size_t const& t) 	    	 const;
    art::PtrVector<raw::AuxDetDigit>   const& TriggerMWPC4DigitsPtr(size_t const& t) 	    	 const;

  private:

    // The following methods work for any of the AuxDetDigit/Trigger vectors below and allow us
    // to not repeat the same code for the different detectors
    std::vector<const raw::AuxDetDigit*>       EventAuxDetDigits     (std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const;
    art::PtrVector<raw::AuxDetDigit>           EventAuxDetDigitsPtr  (std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const;
    std::vector<const raw::AuxDetDigit*>       TriggerAuxDetDigits   (size_t const& t,
								      std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const;
    art::PtrVector<raw::AuxDetDigit>    const& TriggerAuxDetDigitsPtr(size_t const& t,
								      std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const;


    art::PtrVector<raw::Trigger>                    fTriggers;
    std::vector< art::PtrVector<raw::RawDigit>    > fTriggerRawDigits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerUpStreamTOFDigits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerDownStreamTOFDigits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerUpStreamAGDigits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerDownStreamAGDigits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerMWPC1Digits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerMWPC2Digits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerMWPC3Digits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerMWPC4Digits;
    std::vector< art::PtrVector<raw::AuxDetDigit> > fTriggerMuonRangeDigits;
    std::vector< art::PtrVector<raw::OpDetPulse>  > fTriggerOpDetPulses;
 
  };   

}// end namespace

#endif //RAWDATAUTILITIES_TRIGGERDIGITUTILITY_H
 
