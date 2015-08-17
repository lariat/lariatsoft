////////////////////////////////////////////////////////////////////////
// Class:       AerogelCherenkovCounterSlicing
// Module Type: producer
// File:        AerogelCherenkovCounterSlicing_module.cc
//
// Generated at Tue Aug 11 18:44:09 2015 by Dung Phan using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

//ROOT Includes
#include <TH1F.h>

// LArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "RawData/TriggerData.h"

//LAriatSoft Includes
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/AGCounterAlg.h"

namespace lrm {
class AerogelCherenkovCounterSlicing;
}

class lrm::AerogelCherenkovCounterSlicing : public art::EDProducer {
public:
  explicit AerogelCherenkovCounterSlicing(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AerogelCherenkovCounterSlicing(AerogelCherenkovCounterSlicing const &) = delete;
  AerogelCherenkovCounterSlicing(AerogelCherenkovCounterSlicing &&) = delete;
  AerogelCherenkovCounterSlicing & operator = (AerogelCherenkovCounterSlicing const &) = delete;
  AerogelCherenkovCounterSlicing & operator = (AerogelCherenkovCounterSlicing &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  AGCounterAlg	fAGCounterAlg;

};


lrm::AerogelCherenkovCounterSlicing::AerogelCherenkovCounterSlicing(fhicl::ParameterSet const & p) : fAGCounterAlg(p)
{
  // Call appropriate produces<>() functions here.
  // Configures param set
  this->reconfigure(p);  

  // Produces LArIATSoft data object 
  produces<std::vector<ldp::AGCounter> >();
}

void lrm::AerogelCherenkovCounterSlicing::produce(art::Event & e)
{
  // Implementation of required member function here.
  std::unique_ptr<std::vector<ldp::AGCounter> > AGCHitsCol(new std::vector<ldp::AGCounter>);

  // Event Handling
  art::Handle< std::vector<raw::AuxDetDigit> > AuxDetDigitHandle;
  e.getByLabel("SlicerInput", AuxDetDigitHandle);
  
  // Grab Aerogel Counter Digits
  std::vector<raw::AuxDetDigit> AGCUSEDigits; // Channel 4 - Board 8 - CALLID_0
  std::vector<raw::AuxDetDigit> AGCUSWDigits; // Channel 5 - Board 8 - CALLID_1
  std::vector<raw::AuxDetDigit> AGCDS1Digits; // Channel 6 - Board 8 - CALLID_0 - Photonis Square 3in
  std::vector<raw::AuxDetDigit> AGCDS2Digits; // Channel 7 - Board 8 - CALLID_1 - Hamamatsu Round 2in
  for(size_t iDig = 0; iDig < AuxDetDigitHandle->size(); ++iDig){
    if((AuxDetDigitHandle->at(iDig).AuxDetName() == "AeroGelUS")&&(AuxDetDigitHandle->at(iDig).Channel() == 0)) AGCUSEDigits.push_back(AuxDetDigitHandle->at(iDig));
    if((AuxDetDigitHandle->at(iDig).AuxDetName() == "AeroGelUS")&&(AuxDetDigitHandle->at(iDig).Channel() == 1)) AGCUSWDigits.push_back(AuxDetDigitHandle->at(iDig));
    if((AuxDetDigitHandle->at(iDig).AuxDetName() == "AeroGelDS")&&(AuxDetDigitHandle->at(iDig).Channel() == 0)) AGCDS1Digits.push_back(AuxDetDigitHandle->at(iDig));
    if((AuxDetDigitHandle->at(iDig).AuxDetName() == "AeroGelDS")&&(AuxDetDigitHandle->at(iDig).Channel() == 1)) AGCDS2Digits.push_back(AuxDetDigitHandle->at(iDig));
  }
  
  fAGCounterAlg.ImportWaveform("USE", AGCUSEDigits);
  fAGCounterAlg.ImportWaveform("USW", AGCUSWDigits);
  fAGCounterAlg.ImportWaveform("DS1", AGCDS1Digits);
  fAGCounterAlg.ImportWaveform("DS2", AGCDS2Digits);
  
  std::vector<std::vector<AGCHits> > AllAGCHits = fAGCounterAlg.AGCHitsWrapper();
  // Linearization
  std::vector<AGCHits> AllAGCHitsLinearized;
  AllAGCHitsLinearized.clear();
  for (size_t i = 0; i < AllAGCHits.size(); i++) {
	AllAGCHitsLinearized.insert(AllAGCHitsLinearized.end(), AllAGCHits.at(i).begin(), AllAGCHits.at(i).end());
  }
  
  // Wrap-Up
  ldp::AGCounter AGCounterObj(AllAGCHitsLinearized);
  (*AGCHitsCol).push_back(AGCounterObj);
  
  // Push AGCHitsCol to the event
  e.put(std::move(AGCHitsCol));
}

void lrm::AerogelCherenkovCounterSlicing::beginJob()
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::endJob()
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::AerogelCherenkovCounterSlicing::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(lrm::AerogelCherenkovCounterSlicing)
