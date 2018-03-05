////////////////////////////////////////////////////////////////////////
// Class:       KaonFilter
// Module Type: filter
// File:        KaonFilter_module.cc
//
// Generated at Fri Sep 25 16:49:10 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AuxDetParticleID.h"
#include <TH2F.h>
#include "art/Framework/Services/Optional/TFileService.h"

#include <iostream>
#include <memory>

class KaonFilter;

class KaonFilter : public art::EDFilter {
public:
  explicit KaonFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  KaonFilter(KaonFilter const &) = delete;
  KaonFilter(KaonFilter &&) = delete;
  KaonFilter & operator = (KaonFilter const &) = delete;
  KaonFilter & operator = (KaonFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const  &fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  TH2F* fKaonPzVsTOF;
  
  std::string fParticleIDModuleLabel;
  std::string fWCTrackModuleLabel;
  std::string fTOFModuleLabel;
  
  double fKaonProbCutOff;

};


KaonFilter::KaonFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool KaonFilter::filter(art::Event & e)
{
  // Implementation of required member function here.
  //Retrieving the Particle IDs from the event record
  art::Handle< std::vector<ldp::AuxDetParticleID> > particleIDCol;
  e.getByLabel(fParticleIDModuleLabel,particleIDCol);

  //Get the collection of WCTracks produced by the WCTrackBuilder module
  art::Handle< std::vector<ldp::WCTrack> > WCTrackColHandle;
  e.getByLabel(fWCTrackModuleLabel,WCTrackColHandle);
  
  //Get the collection of TOF objects produced by the TOF module
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  e.getByLabel(fTOFModuleLabel,TOFColHandle);
  
  //Finding best-guess kaons
  if( particleIDCol->at(0).PDGCode() == 321 || particleIDCol->at(0).PDGCode() == -321 ){
    if( particleIDCol->at(0).KaonProbability() > fKaonProbCutOff ){
      std::cout << "Selected Possible Kaon Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
      fKaonPzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
      return true;
    }
  }
  return false;
   
}

void KaonFilter::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fKaonPzVsTOF = tfs->make<TH2F>("KaonPzVsTOF","Kaon Pz Vs. TOF",160,0,1600,70,10,80);  
}

void KaonFilter::endJob()
{
  // Implementation of optional member function here.
}

void KaonFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fKaonProbCutOff = p.get<float>("KaonProbabilityThreshold",0.5);
  fParticleIDModuleLabel = p.get<std::string>("ParticleIDModuleLabel");
  fWCTrackModuleLabel = p.get<std::string>("WCTrackModuleLabel");
  fTOFModuleLabel = p.get<std::string>("TOFModuleLabel");
  

}

void KaonFilter::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void KaonFilter::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void KaonFilter::respondToOpenInputFile(art::FileBlock const  &fb)
{
  // Implementation of optional member function here.
}

void KaonFilter::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(KaonFilter)
