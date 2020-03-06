////////////////////////////////////////////////////////////////////////
// Class:       HitNumberFilter
// Module Type: filter
// File:        HitNumberFilter_module.cc
//
// Rejects events with too many wire hits. Helpful for placing between
// hit-finder and clustering/tracking modules to save time in processing.
//
// Configured for 2 planes by default (LArIAT). 
//
// Generated at Wed Oct 12 15:23:51 2016 by William Foreman using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include <TH1F.h>
#include <TH2F.h>

// C++ includes
#include <iostream>
#include <memory>

class HitNumberFilter;

class HitNumberFilter : public art::EDFilter {
public:
  explicit HitNumberFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitNumberFilter(HitNumberFilter const &) = delete;
  HitNumberFilter(HitNumberFilter &&) = delete;
  HitNumberFilter & operator = (HitNumberFilter const &) = delete;
  HitNumberFilter & operator = (HitNumberFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  std::string fHitsModuleLabel;
  std::string fHitsInstance;
  std::vector<int> fMaxNumHits;
  
  TH1F* hNumHits[2];
  TH1F* hNumHits_pass[2];
  TH1F* hEventPass;

};


HitNumberFilter::HitNumberFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  art::ServiceHandle<art::TFileService> tfs;
  hNumHits[0]       = tfs->make<TH1F>("NumHits_Ind",";Number of reconstructed hits (induction plane)",200,0,1000);
  hNumHits[1]       = tfs->make<TH1F>("NumHits_Col",";Number of reconstructed hits (collection plane)",200,0,1000);
  hNumHits_pass[0]       = tfs->make<TH1F>("NumHits_Ind_pass",";Number of reconstructed hits for events passing filter (induction plane)",200,0,1000);
  hNumHits_pass[1]       = tfs->make<TH1F>("NumHits_Col_pass",";Number of reconstructed hits for evtents passing filter (collection plane)",200,0,1000);
  hEventPass        = tfs->make<TH1F>("EventPass","0 = rejected events, 1 = passed events",2,0,2);
}

bool HitNumberFilter::filter(art::Event & e)
{
  bool pass = true;

  // Get the hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (e.getByLabel(fHitsModuleLabel,fHitsInstance,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}
  int nHits[2]={0};
  for(size_t i=0; i<hitlist.size(); i++){
    int hitPlane = hitlist[i]->WireID().Plane;
    if( hitPlane < 2 ) nHits[hitPlane]++;
  }
  hNumHits[0]->Fill(nHits[0]);
  hNumHits[1]->Fill(nHits[1]);
  
  LOG_VERBATIM("HitNumberFilter")
  <<"---- HitNumberFilter -----";
  for(size_t i=0; i<2; i++){
  LOG_VERBATIM("HitNumberFilter")
  <<"#hits on plane "<<i<<": "<<nHits[i];
    if( fMaxNumHits[i] > 0 && nHits[i] > fMaxNumHits[i] )
      pass = false;
  }
  
  LOG_VERBATIM("HitNumberFilter")
  <<"Does event pass? "<<pass<<"\n"
  <<"--------------------------";
 
  if( pass ) {
    hNumHits_pass[0]->Fill(nHits[0]);
    hNumHits_pass[1]->Fill(nHits[1]);
  }
  
  hEventPass->Fill((int)pass);
  return pass;
  
}

void HitNumberFilter::beginJob()
{
}


void HitNumberFilter::endJob()
{
}


void HitNumberFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fMaxNumHits	      = p.get< std::vector<int> > ("MaxNumHits",{-999,-999});
  fHitsModuleLabel    = p.get< std::string >      ("HitsModuleLabel","gaushit");
  fHitsInstance       = p.get< std::string >      ("HitsInstance","");
}


DEFINE_ART_MODULE(HitNumberFilter)
