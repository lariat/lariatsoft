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
// TriggerDigitUtility.cxx
//
// \author Brian Rebel brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// ART Framework Includes
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Framework/Core/FindManyP.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "RawDataUtilities/TriggerDigitUtility.h"

namespace rdu{

  //----------------------------------------------------------------------------
  // The triggerModuleLabel is the label of the module creating both the triggers
  // and the associations to the triggers
  TriggerDigitUtility::TriggerDigitUtility(art::Event  const& evt,
					   std::string const& triggerModuleLabel)
  {

    // get the handle to the triggers out of the event record
    art::Handle< std::vector<raw::Trigger> > trigHandle;
    evt.getByLabel(triggerModuleLabel, trigHandle);

    if( !trigHandle.isValid() )
      throw cet::exception("TriggerDigitUtility") << "trigger handle from module "
						  << triggerModuleLabel
						  << " is not valid.";

    // loop over all the triggers and get the digits associated with each
    size_t numTrigs = trigHandle->size();
    
    // declare the FindManyP objects to get the associated digits for the triggers
    art::FindManyP<raw::RawDigit>    fmprd (trigHandle, evt, triggerModuleLabel);
    art::FindManyP<raw::AuxDetDigit> fmpadd(trigHandle, evt, triggerModuleLabel);
    art::FindManyP<raw::OpDetPulse>  fmpodp(trigHandle, evt, triggerModuleLabel);

    // resize the vectors first
    fTriggerRawDigits          .resize(numTrigs);
    fTriggerUpStreamTOFDigits  .resize(numTrigs);
    fTriggerDownStreamTOFDigits.resize(numTrigs);
    fTriggerUpStreamAGDigits   .resize(numTrigs);
    fTriggerDownStreamAGDigits .resize(numTrigs);
    fTriggerMuonRangeDigits    .resize(numTrigs);
    fTriggerMWPC1Digits        .resize(numTrigs);
    fTriggerMWPC2Digits        .resize(numTrigs);
    fTriggerMWPC3Digits        .resize(numTrigs);
    fTriggerMWPC4Digits        .resize(numTrigs);
    fTriggerTriggerWaveForms   .resize(numTrigs);
    fTriggerOpDetPulses        .resize(numTrigs);

    for(size_t t = 0; t < numTrigs; ++t){

      fTriggers.push_back(art::Ptr<raw::Trigger>(trigHandle, t));

      std::vector<art::Ptr<raw::RawDigit>    > rd  = fmprd .at(t);
      std::vector<art::Ptr<raw::AuxDetDigit> > add = fmpadd.at(t);
      std::vector<art::Ptr<raw::OpDetPulse>  > odp = fmpodp.at(t);

      for(auto rdp  : rd ) fTriggerRawDigits  [t].push_back(rdp);
      for(auto odpp : odp) fTriggerOpDetPulses[t].push_back(odpp); 
      for(auto addp : add){
	
	std::string const& detName = addp->AuxDetName();
	
	if     (detName.find("TOFUS")          != std::string::npos) fTriggerUpStreamTOFDigits[t]  .push_back(addp);
	else if(detName.find("TOFDS")          != std::string::npos) fTriggerDownStreamTOFDigits[t].push_back(addp);
	else if(detName.find("AeroGelUS")      != std::string::npos) fTriggerUpStreamAGDigits[t]   .push_back(addp);
	else if(detName.find("AeroGelDS")      != std::string::npos) fTriggerDownStreamAGDigits[t] .push_back(addp);
	else if(detName.find("MuonRangeStack") != std::string::npos) fTriggerMuonRangeDigits[t]    .push_back(addp);
	else if(detName.find("MWPC1")          != std::string::npos) fTriggerMWPC1Digits[t]        .push_back(addp);
	else if(detName.find("MWPC2")          != std::string::npos) fTriggerMWPC1Digits[t]        .push_back(addp);
	else if(detName.find("MWPC3")          != std::string::npos) fTriggerMWPC1Digits[t]        .push_back(addp);
	else if(detName.find("MWPC4")          != std::string::npos) fTriggerMWPC1Digits[t]        .push_back(addp);
	else if(detName.find("WC1")            != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("WC2")            != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);   
	else if(detName.find("WC3")            != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("WC4")            != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("BEAMON")         != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp); 
	else if(detName.find("USTOF")          != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp); 
	else if(detName.find("DSTOF")          != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("PUNCH")          != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("HALO")           != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);   
	else if(detName.find("PULSER")         != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("COSMICON")       != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("COSMIC")         != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("PILEUP")         != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("MICHEL")         != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("LARSCINT")       != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);
	else if(detName.find("MuRS")           != std::string::npos) fTriggerTriggerWaveForms[t]   .push_back(addp);

      } // end loop over AuxDetDigits for this trigger

    } // end loop over triggers

    return;
  }

  //----------------------------------------------------------------------------
  TriggerDigitUtility::~TriggerDigitUtility()
  {
  }

  //----------------------------------------------------------------------------
  size_t TriggerDigitUtility::NTriggers() const
  {
    return fTriggers.size();
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::Trigger*> TriggerDigitUtility::EventTriggers() const
  {
    // loop over the PtrVector and make the vector of raw::Triggers
    std::vector<const raw::Trigger*> trig;
    trig.reserve(fTriggers.size());

    for(size_t t = 0; t < fTriggers.size(); ++t) trig.push_back(fTriggers[t].get());
    
    return trig;
  }
				  
  //----------------------------------------------------------------------------
  art::PtrVector<raw::Trigger> const& TriggerDigitUtility::EventTriggersPtr() const
  {
    return fTriggers;
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::RawDigit*> TriggerDigitUtility::EventRawDigits() const
  {
    std::vector<const raw::RawDigit*> rdvec;
    for(auto rdpv : fTriggerRawDigits){
    std::cout<<std::endl;
    std::cout<<"rdpv.size() = "<<rdpv.size()<<std::endl;
    std::cout<<std::endl;
      for(size_t r = 0; r < rdpv.size(); ++r) rdvec.push_back(rdpv[r].get());
    }

    return rdvec;
  }
			   	 
  //----------------------------------------------------------------------------
  art::PtrVector<raw::RawDigit> TriggerDigitUtility::EventRawDigitsPtr() const
  {
    art::PtrVector<raw::RawDigit> rdPtrVec;

    for(auto rdpv : fTriggerRawDigits){
      for(size_t r = 0; r < rdpv.size(); ++r) rdPtrVec.push_back(rdpv[r]);
    }
    
    return rdPtrVec;
  }
 
  //----------------------------------------------------------------------------
  std::vector<const raw::RawDigit*> TriggerDigitUtility::TriggerRawDigits(size_t const& t) const
  {
    std::vector<const raw::RawDigit*> rdvec;

    if(t > fTriggerRawDigits.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested raw digits from trigger "
						  << t << " but only " << fTriggerRawDigits.size()
						  << " triggers in this event";

    for(size_t r = 0; r < fTriggerRawDigits[t].size(); ++r) rdvec.push_back(fTriggerRawDigits[t][r].get());
    std::cout<<"rdvec.size() = "<<rdvec.size()<<std::endl;
    return rdvec;
  }
			   	 
  //----------------------------------------------------------------------------
  art::PtrVector<raw::RawDigit> const& TriggerDigitUtility::TriggerRawDigitsPtr(size_t const& t) const
  {
  
    art::PtrVector<raw::RawDigit> rdPtrVec; 
    
    if(t > fTriggerRawDigits.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested raw digits from trigger "
						  << t << " but only " << fTriggerRawDigits.size()
						  << " triggers in this event";
    return fTriggerRawDigits[t];
    
    //Edit JAsaadi
    /*for(size_t r = 0; r < fTriggerRawDigits[t].size(); ++r) rdPtrVec.push_back(fTriggerRawDigits[r].get());
    std::cout<<"rdPtrVec.size() = "<<rdPtrVec.size()<<std::endl;
    return rdPtrVec;*/
    
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::OpDetPulse*> TriggerDigitUtility::EventOpDetPulses() const
  {
    std::vector<const raw::OpDetPulse*> odpvec;

    for(auto odpv : fTriggerOpDetPulses){
      for(size_t o = 0; o < odpv.size(); ++o) odpvec.push_back(odpv[o].get());
    }

    return odpvec;
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::OpDetPulse> TriggerDigitUtility::EventOpDetPulsesPtr() const
  {
    art::PtrVector<raw::OpDetPulse> odpPtrVec;

    for(auto odpv : fTriggerOpDetPulses){
      for(size_t o = 0; o < odpv.size(); ++o) odpPtrVec.push_back(odpv[o]);
    }
    
    return odpPtrVec;
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::OpDetPulse*> TriggerDigitUtility::TriggerOpDetPulses(size_t const& t) const
  {
    std::vector<const raw::OpDetPulse*> odpvec;

    if(t > fTriggerOpDetPulses.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested optical detector pulses from trigger "
						  << t << " but only " << fTriggerOpDetPulses.size()
						  << " triggers in this event";

    for(size_t o = 0; o < fTriggerOpDetPulses[t].size(); ++o) odpvec.push_back(fTriggerOpDetPulses[t][o].get());

    return odpvec;
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::OpDetPulse> const& TriggerDigitUtility::TriggerOpDetPulsesPtr(size_t const& t) const
  {
    if(t > fTriggerOpDetPulses.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested optical detector pulses from trigger "
						  << t << " but only " << fTriggerOpDetPulses.size()
						  << " triggers in this event";
    return fTriggerOpDetPulses[t];
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventAuxDetDigits(std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const
  {
    std::vector<const raw::AuxDetDigit*> addvec;

    for(auto addv : trigDigVec){
      for(size_t a = 0; a < addv.size(); ++a) addvec.push_back(addv[a].get());
    }

    return addvec;
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventAuxDetDigitsPtr(std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const
  {
    art::PtrVector<raw::AuxDetDigit> addPtrVec;

    for(auto addv : trigDigVec){
      for(size_t a = 0; a < addv.size(); ++a) addPtrVec.push_back(addv[a]);
    }
    
    return addPtrVec;
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerAuxDetDigits(size_t const& t,
									       std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const
  {
    std::vector<const raw::AuxDetDigit*> addvec;

    if(t > trigDigVec.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested aux detector digits from trigger "
						  << t << " but only " << trigDigVec.size()
						  << " triggers in this event";

    for(size_t a = 0; a < trigDigVec[t].size(); ++a) addvec.push_back(trigDigVec[t][a].get());

    return addvec;
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerAuxDetDigitsPtr(size_t const& t,
										      std::vector< art::PtrVector<raw::AuxDetDigit> > const& trigDigVec) const
  {
    if(t > trigDigVec.size()) 
      throw cet::exception("TriggerDigitUtility") << "requested aux detector digits from trigger "
						  << t << " but only " << trigDigVec.size()
						  << " triggers in this event";
    return trigDigVec[t];
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventMuonRangeStackDigits() const
  {
    return this->EventAuxDetDigits(fTriggerMuonRangeDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventMuonRangeStackDigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerMuonRangeDigits);
  }
                     
  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerMuonRangeStackDigits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerMuonRangeDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerMuonRangeStackDigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerMuonRangeDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventUpStreamTOFDigits() const
  {
    return this->EventAuxDetDigits(fTriggerUpStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventDownStreamTOFDigits() const
  {
    return this->EventAuxDetDigits(fTriggerDownStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventUpStreamTOFDigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerUpStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventDownStreamTOFDigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerDownStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerUpStreamTOFDigits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerUpStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerDownStreamTOFDigits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerDownStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerUpStreamTOFDigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerUpStreamTOFDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerDownStreamTOFDigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerDownStreamTOFDigits);
  }

  
  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventUpStreamAeroGelDigits() const
  {
    return this->EventAuxDetDigits(fTriggerUpStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventDownStreamAeroGelDigits() const
  {
    return this->EventAuxDetDigits(fTriggerDownStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventUpStreamAeroGelDigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerUpStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventDownStreamAeroGelDigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerDownStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerUpStreamAeroGelDigits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerUpStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerDownStreamAeroGelDigits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerDownStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerUpStreamAeroGelDigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerUpStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerDownStreamAeroGelDigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerDownStreamAGDigits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventMWPC1Digits() const
  {
    return this->EventAuxDetDigits(fTriggerMWPC1Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventMWPC2Digits() const
  {
    return this->EventAuxDetDigits(fTriggerMWPC2Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventMWPC3Digits() const
  {
    return this->EventAuxDetDigits(fTriggerMWPC3Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventMWPC4Digits() const
  {
    return this->EventAuxDetDigits(fTriggerMWPC4Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventMWPC1DigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerMWPC1Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventMWPC2DigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerMWPC2Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventMWPC3DigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerMWPC3Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventMWPC4DigitsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerMWPC4Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerMWPC1Digits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerMWPC1Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerMWPC2Digits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerMWPC2Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerMWPC3Digits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerMWPC3Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerMWPC4Digits(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerMWPC4Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerMWPC1DigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerMWPC1Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerMWPC2DigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerMWPC2Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerMWPC3DigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerMWPC3Digits);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerMWPC4DigitsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerMWPC4Digits);
  }

  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::EventTriggerWaveForms() const
  {
    return this->EventAuxDetDigits(fTriggerTriggerWaveForms);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> TriggerDigitUtility::EventTriggerWaveFormsPtr() const
  {
    return this->EventAuxDetDigitsPtr(fTriggerTriggerWaveForms);
  }
                     
  //----------------------------------------------------------------------------
  std::vector<const raw::AuxDetDigit*> TriggerDigitUtility::TriggerTriggerWaveForms(size_t const& t) const
  {
    return this->TriggerAuxDetDigits(t, fTriggerTriggerWaveForms);
  }

  //----------------------------------------------------------------------------
  art::PtrVector<raw::AuxDetDigit> const& TriggerDigitUtility::TriggerTriggerWaveFormsPtr(size_t const& t) const
  {
    return this->TriggerAuxDetDigitsPtr(t, fTriggerTriggerWaveForms);
  }


} // end namespace
