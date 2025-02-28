////////////////////////////////////////////////////////////////////////
// Class:       NNShowerFilter
// Module Type: filter
// File:        NNShowerFilter_module.cc
//
// Daniel Smith
// Boston University
// 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AuxDetParticleID.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h" 

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"

#include <iostream>
#include <memory>
#include <TH2F.h>

#define MVA_LENGTH 4

class NNShowerFilter;

class NNShowerFilter : public art::EDFilter {
public:

  explicit NNShowerFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NNShowerFilter(NNShowerFilter const &) = delete;
  NNShowerFilter(NNShowerFilter &&) = delete;
  NNShowerFilter & operator = (NNShowerFilter const &) = delete;
  NNShowerFilter & operator = (NNShowerFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  // Declare member data here.
  
  std::string fNNetModuleLabel;
  std::string fHitModuleLabel;
  std::string fTrackModuleLabel;
  std::string fWCTrackLabel;

  TH1F* fTrackProb; 
  TH1F* fShowerProb;
  TH1F* fMichelProb;
  TH1F* fEmptyProb;

  TH1F* fProbBoxShower;
  TH1F* fProbBoxTrack;
  TH1F* fProbBoxShowerMult;
  TH2F* fProbBoxTrackVsMomo;
  
  TH1F* fRejectionModes;
  TH1F* fPassingEvents;
  TH2F* fPassingEventsVsPz;
  TH1F* fTrackReNormProb; 
  TH1F* fShowerReNormProb;

  //TH2F* fReNormTrackVsEnergy;
  //TH2F* fReNormTrackVsWireID;

  TH1F* fMomo;

  TH1F* fXface;
  TH1F* fYface;
  TH1F* fDistface;
  TH1F* fAlpha;

  float fShowerThresh;
  bool bData;
  bool bVerbose;
  bool bRemoveOrSelect;
  bool DoWCTPCMatch;

};


NNShowerFilter::NNShowerFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
  std::cout << "All fhicle file stuff done" << std::endl; 
}

bool NNShowerFilter::filter(art::Event & e)
{

  // NN data product
  anab::MVAReader<recob::Hit, 4> hitResults(e, fNNetModuleLabel);
  std::vector<anab::FeatureVector<4>> featVec = hitResults.outputs();

  // Hit info
  art::Handle< std::vector<recob::Hit> > HitHandle;
  std::vector< art::Ptr<recob::Hit> > Hitlist;
  if(e.getByLabel(fHitModuleLabel, HitHandle))
    {//std::cout<<" Filling Hit List."<<std::endl;
     art::fill_ptr_vector(Hitlist, HitHandle);
    }

  // Track info
  art::Handle< std::vector<recob::Track> > TrackHandle;
  std::vector< art::Ptr<recob::Track> > Tracklist;

  if(e.getByLabel(fTrackModuleLabel, TrackHandle))
    {//std::cout<<" Filling Track List."<<std::endl;
    art::fill_ptr_vector(Tracklist, TrackHandle);
    }
  // Association between Tracks and 2d Hits
  art::FindManyP<recob::Track> ass_trk_hits(HitHandle,   e, fTrackModuleLabel); 
  art::FindManyP<recob::Hit> HitsInTrackAssn(TrackHandle,   e, fTrackModuleLabel); 



  // Get some basic histograms of what is going onss

  if(bVerbose) {
    std::cout << "Hitlist.size()=" << Hitlist.size()
              << " Tracklist.size()=" << Tracklist.size()  
	      << " featVec.size()=" << featVec.size() << std::endl;
  }

  for(size_t iHit = 0; iHit < Hitlist.size(); iHit++) {
    // Only take from collection view, only exists for col
    if(Hitlist[iHit]->View() != 1) { continue; } 

    fTrackProb->Fill(featVec[iHit][0]);
    fShowerProb->Fill(featVec[iHit][1]);
    fMichelProb->Fill(featVec[iHit][2]);
    fEmptyProb->Fill(featVec[iHit][3]);

    double renormTrack = featVec[iHit][0] / (featVec[iHit][0] + featVec[iHit][1]);
    double renormShower = featVec[iHit][1] / (featVec[iHit][0] + featVec[iHit][1]);

    fTrackReNormProb->Fill(renormTrack);
    fShowerReNormProb->Fill(renormShower);

   // fReNormTrackVsEnergy->Fill(renormTrack, Hitlist[iHit]->SummedADC());
   // fReNormTrackVsWireID->Fill(renormTrack, (float)Hitlist[iHit]->WireID().Wire);
  }



  // Need the WC projection only the TPC face
  double xface = 0.0; double yface = 0.0;
  double momo = 0.0;

  std::vector<double> WC_vec = {-1.0, -1.0, -1.0};

  if(bData) {
    if(bVerbose) { std::cout << " Data event. " << std::endl; }

    // WC info
    art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
    std::vector<art::Ptr<ldp::WCTrack> > wctrack;
   
    if(e.getByLabel(fWCTrackLabel, wctrackHandle))
      art::fill_ptr_vector(wctrack, wctrackHandle);

    // requires single WC track
    if(bVerbose) { std::cout << "wcsize = " << wctrack.size() << std::endl; }
    if(wctrack.size() != 1) { fRejectionModes->Fill(1.0); return false; }

    xface = wctrack[0]->XYFace(0);
    yface = wctrack[0]->XYFace(1);

    momo = wctrack[0]->Momentum();
    fMomo->Fill(wctrack[0]->Momentum());

    WC_vec[0] = wctrack[0]->HitPosition(3, 0) - wctrack[0]->HitPosition(2, 0);
    WC_vec[1] = wctrack[0]->HitPosition(3, 1) - wctrack[0]->HitPosition(2, 1);
    WC_vec[2] = wctrack[0]->HitPosition(3, 2) - wctrack[0]->HitPosition(2, 2);

    if(bVerbose) { std::cout << "Data WC projection (" << xface << ", " << yface << ")" << std::endl; }

  } else {
    if(bVerbose) { std::cout << " MC event. " << std::endl; }

    // ParticleInventoryService service ===
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();

    std::vector<const simb::MCParticle* > geant_part;     
    for(size_t p = 0; p < plist.size(); ++p) { geant_part.push_back(plist.Particle(p)); }


    fPassingEvents->Fill(0); // start of all 


    bool found = false;
    for(unsigned int i = 0; i < geant_part.size(); ++i ) {
      if(geant_part[i]->Process() != "primary") { continue; }

      // <x, y, z> = <sx, sy, sz> * t + <startx, starty, startz>
      // (z - startz) / pz = t
      // x = px * t + startx
      // y = px * t + srarty
      
      // intercept time

      double t = (0.0 - geant_part[i]->Vz()) / geant_part[i]->Pz(); 

      xface = geant_part[i]->Px() * t + geant_part[i]->Vx(); 
      yface = geant_part[i]->Py() * t + geant_part[i]->Vy();
    
      TVector3 wcVec;
      wcVec.SetX(geant_part[i]->Px());
      wcVec.SetY(geant_part[i]->Py());
      wcVec.SetZ(geant_part[i]->Pz());

      momo = wcVec.Mag()*1000.;
      fMomo->Fill(momo);

      WC_vec[0] = wcVec[0];
      WC_vec[1] = wcVec[1];
      WC_vec[2] = wcVec[2];

      found = true;
    }

    if(!found) { 
      if(bVerbose) { std::cout << "Could not find an MC WC-TPC match!" << std::endl; } 
      fRejectionModes->Fill(1.0); 
      return false; 
    }

    if(bVerbose) { std::cout << "MC TPC face intercept - " << xface << " " << yface << std::endl; }

  }

  fPassingEvents->Fill(1); // Pass WC existing
  //fPassingEventsVsPz->Fill(1, momo);

  // ----- 
  // Need to create box on Wire / TDC view starting from where WC projects into the event
  // Complicated because we are going to 3D space to Wire / TDC, reverse to normal progression
  //
  // First, start with tracks (3D object)
  //   find track closest to the projected WC track
  // If track is close enough, find associated hit that is lowest wire number
  // This hit is start point of box, then run for every hit in box
  // ---

  // First, find closest track to the projected WC track
  double closest_track_index = 0;
  double closest_dist=99999;
  int nMatched = 0.0;

  for(size_t i = 0; i < Tracklist.size(); i++) {

    if(bVerbose) { 
      std::cout << "Start of track: " << Tracklist[i]->Start().X() << " " << Tracklist[i]->Start().Y() << " " << Tracklist[i]->Start().Z() << std::endl;
      std::cout << "End of track: "   << Tracklist[i]->End().X()   << " " << Tracklist[i]->End().Y()   << " " << Tracklist[i]->End().Z()   << std::endl;
    }
    
    // Require track starts (or ends if backwards) within the first 10 cm of TPC
    if(Tracklist[i]->Start().Z() > 5.0) { continue; }

    auto larStartEnd = Tracklist[i]->Direction();

    TVector3 wcz(WC_vec[0], WC_vec[1], WC_vec[2]);
    TVector3 tpcz(larStartEnd.first.X(), larStartEnd.first.Y(), larStartEnd.first.Z());

    //larStartEnd.first.X();

    // track length requirement ... we wills ee if this does anything good for us

    /*
    if(TMath::Sqrt(TMath::Power(Tracklist[i]->Start().X() - Tracklist[i]->End().X(), 2) + 
		   TMath::Power(Tracklist[i]->Start().Y() - Tracklist[i]->End().Y(), 2) +
		   TMath::Power(Tracklist[i]->Start().Z() - Tracklist[i]->End().Z(), 2)) < 5.0) { continue; }
    */

    fXface->Fill(Tracklist[i]->Start().X() - xface);
    fYface->Fill(Tracklist[i]->Start().Y() - yface);
    fDistface->Fill(TMath::Sqrt(TMath::Power(Tracklist[i]->Start().X() - xface, 2) + 
				TMath::Power(Tracklist[i]->Start().Y() - yface, 2)));
    fAlpha->Fill(acos(tpcz.Dot(wcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159));


    // Track spatial match using wc2tpc filer values dx=[-2,6]cm dy=[-3,6]cm

   
    // Alpha
    std::cout << "alpha =" << acos(tpcz.Dot(wcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159) << std::endl;
    if(DoWCTPCMatch)
    {
    if(acos(wcz.Dot(tpcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159) > 20.0) { continue; }
    if(
      (Tracklist[i]->Start().X() - xface) < -2 || 
      (Tracklist[i]->Start().X() - xface) > 6  ||
      (Tracklist[i]->Start().Y() - xface) > -3 ||
      (Tracklist[i]->Start().Y() - xface) > 6 ||
       acos(wcz.Dot(tpcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159) > 20.0) {continue;}
    std::cout << "space match " << std::endl;
    closest_track_index = i;
    
    nMatched += 1;
    }
   if(!DoWCTPCMatch)
   {
    double temp_dist=TMath::Sqrt(TMath::Power(Tracklist[i]->Start().X() - xface, 2) + TMath::Power(Tracklist[i]->Start().Y() - yface, 2));
    if(temp_dist<closest_dist)
    {
      closest_dist=temp_dist;
      closest_track_index=i;
    }
   }
  }


  // If it didn't change, couldn't find a track match, so we Reject this event
  if(nMatched != 1) {
    if(bVerbose) { std::cout << "Found " <<nMatched<< " tracks. There can only be one!"<<std::endl; }
    fRejectionModes->Fill(2.0);
    return false; 
  }

  std::cout << "Matched track!" << std::endl;


  // Now, we find the lowest Hit on that very same track (hit with lowest wireid)

  // This is wrong! Ha. What if wires aren't sequential like this?

  if(!HitsInTrackAssn.isValid()) { fRejectionModes->Fill(3.0); return false; }


  int lowest_hit = -1;
  std::vector<art::Ptr<recob::Hit> > trackhits = HitsInTrackAssn.at(closest_track_index);
  for (size_t iHit =0 ; iHit < trackhits.size(); ++iHit) {
    if(trackhits[iHit]->View() != 1) { continue; }    
    if(lowest_hit == -1) { lowest_hit = iHit; }       
    if(trackhits[iHit]->WireID().Wire < trackhits[lowest_hit]->WireID().Wire) { lowest_hit = iHit; }
  }


  //if(bVerbose) { std::cout << "lowest hit number = " << lowest_hit << " at wire ID " << Hitlist[lowest_hit]->WireID() << std::endl; }
  
  // If we for some reason didn't find any hits on this track, we reject this event
  if(lowest_hit == -1) {   
    if(bVerbose) { std::cout << "lowest_hit never initialized" << std::endl; } 
    fRejectionModes->Fill(4.0);
    return false; 
  }

  float inter = trackhits[lowest_hit]->PeakTime(); //Hitlist[lowest_hit]->PeakTime();
  float offset = trackhits[lowest_hit]->WireID().Wire;

  if(bVerbose) { std::cout << " Network choose to start the box at: (wire, tdc) = " << offset << ", " << inter << std::endl; } 

  if(offset > 100.0) { return false; } // make sure things don't get out of hand

  // If inter is greater than the tdcs, reject it
  if(inter > 3000.) { 
    if(bVerbose) { std::cout << " Failed match ... " << std::endl; } 
    fRejectionModes->Fill(5.0);
    return false; 
  }


  double total_shower_prob = 0.0;

  double mult_shower = -1.0;
  double trac_shower = -1.0;

  int num_hits = 0;

  for(size_t jHit = 0; jHit < Hitlist.size(); ++jHit) {
    //std::cout << jHit << " of " << Hitlist.size() << " ";

    if(Hitlist[jHit]->View() != 1) { continue; }    

    int wireID = Hitlist[jHit]->WireID().Wire;
    int hitTime = Hitlist[jHit]->PeakTime();

    // Box --- --- 
    if(wireID > (offset+100.0) or wireID < offset) { continue; }
    if(hitTime > inter + 200. or hitTime < inter - 200.) { continue; }
    // --- --- ---

    if(hitTime > 3000.) {  if(bVerbose) { std::cout << "too hit high in time ... " << std::endl; } continue; }
    if(wireID > 240.)   {  if(bVerbose) { std::cout << "too hit high in wire ... " << std::endl; } continue; }

    // if(bVerbose) { std::cout << " wireID, hitTime about to do vector " << wireID << " " << hitTime << std::endl; }

    std::vector<float> results;
    results.push_back(featVec[jHit][0]);
    results.push_back(featVec[jHit][1]);
    results.push_back(featVec[jHit][2]);
    results.push_back(featVec[jHit][3]);

    if(bVerbose) { std::cout << "CNN results " << results[0] << " " << results[1] << " " << results[2] << " " << std::endl; }

    //results[0] += results[3];

    std::vector<float> temp_results = results;

    results[0] = temp_results[0] / (temp_results[0] + temp_results[1]);
    results[1] = temp_results[1] / (temp_results[0] + temp_results[1]);

    if(bVerbose) { std::cout << "Normalized CNN results " << results[0] << " " << results[1] << " " << std::endl; }

    total_shower_prob += results[1];

    if(trac_shower == -1.0) { trac_shower = results[0]; }
    if(mult_shower == -1.0) { mult_shower = results[1]; }

    trac_shower *= results[0];
    mult_shower *= results[1];

    num_hits += 1;
  }

  total_shower_prob = total_shower_prob / double(num_hits);

  fProbBoxShower->Fill(total_shower_prob);
  fProbBoxTrack->Fill(1-total_shower_prob);
  fProbBoxShowerMult->Fill(mult_shower / (trac_shower + mult_shower));
  fProbBoxTrackVsMomo->Fill(1-total_shower_prob, momo);

  if(bVerbose) { std::cout << "shower prob = " << total_shower_prob << std::endl; }

  if(total_shower_prob < fShowerThresh) {
    fPassingEvents->Fill(6); // My new EM cut 
    fPassingEventsVsPz->Fill(6, momo);
  }  

  // ifbRemoveOrSelect == true, Select things above the threshold
  // else, Remove things above the threshold
  if(bRemoveOrSelect) { return (total_shower_prob > fShowerThresh); }
  else { return (total_shower_prob < fShowerThresh); }


}

void NNShowerFilter::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTrackProb = tfs->make<TH1F>("TrackProb","TrackProb",100,0.0,1.0);
  fShowerProb = tfs->make<TH1F>("ShowerProb","ShowerProb",100,0.0,1.0);
  fMichelProb = tfs->make<TH1F>("MichelProb","MichelProb",100,0.0,1.0);
  fEmptyProb = tfs->make<TH1F>("EmptyProb","EmptyProb",100,0.0,1.0);

  fTrackReNormProb = tfs->make<TH1F>("TrackReNormProb","TrackReNormProb",100,0.0,1.0);
  fShowerReNormProb = tfs->make<TH1F>("ShowerReNormProb","ShowerReNormProb",100,0.0,1.0);

  fProbBoxShower = tfs->make<TH1F>("fProbBoxShower","fProbBoxShower",100, 0.0, 1.0); //, 500, 0, 500);
  fProbBoxTrack = tfs->make<TH1F>("fProbBoxTrack","fProbBoxTrack",100, 0.0, 1.0); //, 500, 0, 500);
  fProbBoxShowerMult = tfs->make<TH1F>("fProbBoxShowerMult","fProbBoxShowerMult",200, 0.0, 1.0); //, 500, 0, 500);
  fProbBoxTrackVsMomo  = tfs->make<TH2F>("fProbBoxTrackVsMomo","fProbBoxTrackVsMomo",100, 0.0, 1.0, 200, 0, 2000); //, 500, 0, 500);

  //fReNormTrackVsEnergy = tfs->make<TH2F>("ReNormTrackVsEnergy","ReNormTrackVsEnergy",100,0.0,1.0, 100, 0.0, 10000.0);
  //fReNormTrackVsWireID = tfs->make<TH2F>("ReNormTrackVsWireID","ReNormTrackVsWireID",100,0.0,1.0, 240, 0.0, 240.0);

  fRejectionModes = tfs->make<TH1F>("RejectionModes","RejectionModes",10, 0.0, 10.0);

  fPassingEvents = tfs->make<TH1F>("PassingEvents","PassingEvents",15, 0.0, 15.0);
  fPassingEventsVsPz = tfs->make<TH2F>("PassingEventsVsPz","PassingEventsVsPz",15, 0.0, 15.0, 200, 0, 2000);

  fMomo = tfs->make<TH1F>("Momo","Momo",200, 0.0, 2000.0);

  fXface = tfs->make<TH1F>("XFace","XFace",100, -20.0, 20.0);
  fYface = tfs->make<TH1F>("YFace","YFace",100, -20.0, 20.0);
  fDistface = tfs->make<TH1F>("DistFace","DistFace",100, 0.0, 50.0);
  fAlpha = tfs->make<TH1F>("Alpha","Alpha",100, 0.0, 50.0);


}

void NNShowerFilter::endJob()
{
  // Implementation of optional member function here.
  //  TH2F* fProbBoxShowerVsMomo;
  //TH1F* fMomo;

/*   for(int x = 1; x < fPassingEventsVsPz->GetNbinsX()+1; x++) {
    for(int y = 1; y < fPassingEventsVsPz->GetNbinsY()+1; y++) {
      //std::cout << fProbBoxShowerVsMomo->GetBinContent(x, y) << std::endl;
      if(fMomo->GetBinContent(y) != 0) {
	/fPassingEventsVsPz->SetBinContent(x, y, fPassingEventsVsPz->GetBinContent(x, y) / fMomo->GetBinContent(y));
      }
    }
  } */


}


void NNShowerFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fNNetModuleLabel = p.get<std::string>("NNetModuleLabel"); 
  fHitModuleLabel = p.get<std::string>("HitModuleLabel"); 
  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
  fWCTrackLabel = p.get<std::string>("WCTrackLabel");
  DoWCTPCMatch=p.get<bool>("DoWCTPCMatch",true);
  fShowerThresh = p.get<float>("ShowerThresh");
  bData = p.get<bool>("Data");
  bVerbose = p.get<bool>("Verbose");
  bRemoveOrSelect = p.get<bool>("RemoveOrSelect");

}

DEFINE_ART_MODULE(NNShowerFilter)
