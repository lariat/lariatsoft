////////////////////////////////////////////////////////////////////////
// Class:       StoppingTracks
// Module Type: filter
// File:        StoppingTracks_module.cc
//
// Generated at Wed Jan  6 11:59:49 2016 by Elena Gramellini using artmod
// from cetpkgsupport v1_08_06.
// Calorimetry and tagging for stopping tracks added by Irene Nutini on January 2016

// This filter works as follows:
// The filter rejects an event if all the "incoming" 
// tracks are stopping inside the TPC.
// 
// ### Failure mode in view of pion analysis ###
// Say that you have 2 incoming tracks
// Track1 is a stopping particle, Track2 is not.
// This filter would keep this event.
// You later find out that Track1 is the best match to wcTrack.
// You don't want to keep this track for the PionXSAnalysis!
////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////

// ########################## 
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// ########################### 
// ### LArIATsoft Includes ###
// ########################### 
#include "LArIATDataProducts/WCTrack.h"

// ######################## 
// ### LArsoft Includes ### 
// ########################
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"


// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

// #####################
// ### ROOT includes ###
// #####################
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>

//const int kMaxTrack      = 1000; // unused  //maximum number of tracks

class StoppingTracks;

class StoppingTracks : public art::EDProducer {
public:
  explicit StoppingTracks(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoppingTracks(StoppingTracks const &) = delete;
  StoppingTracks(StoppingTracks &&) = delete;
  StoppingTracks & operator = (StoppingTracks const &) = delete;
  StoppingTracks & operator = (StoppingTracks &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;
  /*
  bool beginRun(art::Run & r) override;
  bool beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  bool endRun(art::Run & r) override;
  bool endSubRun(art::SubRun & sr) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const  &fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;
  */

  bool isStoppingTrack( art::Ptr<recob::Track>  aTrack, art::FindManyP<anab::Calorimetry> fmcal );
 // THIS IS CHEATING, WE CAN PASS BETTER THAN THIS
private:

  // Declare member data here.
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fWC2TPCModuleLabel;
  std::string fWCTrackLabel;

  double fupstreamZPosition;
  
  double fParticleMass;
  int    fMinNofSpacePoints;
  double fLowLimitStop;
  double fUpLimitStop;
  bool   fCheckAssn;
  
  // === Storing the tracks Calorimetry Information
  int    trkhits[2];
  double trkke[2];
  double trkdedx[2][1000];
  double trkrr[2][1000];
  double trkpitchhit[2][1000];

  TH1F  *fTotNTrack;
  TH1F  *fStoppingTrack;
  TH1F  *fDifferenceTrack;
  TH1F  *fNTrackPassedCuts;
  
  TH1F *fdEdx;
  TH2F *fdEdxRange;
  TH1F *fFitPar;
};


StoppingTracks::StoppingTracks(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
  //std::string fTrackModuleLabel;
  //std::string fCalorimetryModuleLabel;

  produces<art::PtrVector<recob::Track> >();
  
}

void StoppingTracks::produce(art::Event & evt)
{
  std::unique_ptr<art::PtrVector<recob::Track> > TrackNonStoppingVector (new art::PtrVector<recob::Track> );
  (*TrackNonStoppingVector).clear();
		
  // #####################################
  // ### Getting the Track Information ###
  // #####################################
   
  //####TPC Track####
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
  
  // === Filling the tracklist from the tracklistHandle ===
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
  art::fill_ptr_vector(tracklist, trackListHandle);}
      
 // === Association between Calorimetry objects and Tracks ===
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  
  
  //std::cout<<"========================================="<<std::endl;
  std::cout<<"Run = "<<evt.run()<<", SubRun = "<<evt.subRun()<<", Evt = "<<evt.id().event()<<std::endl;
  // std::cout<<"========================================="<<std::endl;
  //std::set<int> tpcIndeces;
  //tpcIndeces.clear();
  
  if (fCheckAssn) 
    {
      std::cout<<"Checking Ass\n";
      
      // ###################################################                                                                                                                                
      // ### Getting the Wire Chamber Track  Information ###                                                                                                                                
      // ###################################################                                                                                                                                
      art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
      if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) 
         {evt.put(std::move(TrackNonStoppingVector)); return;}
      // === Association between WC Tracks and TPC Tracks ===                                                                                                                               
      art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
      if (!fWC2TPC.isValid())  
	{
	  evt.put(std::move(TrackNonStoppingVector));
	  return;
	}
      // === Loop on all the Assn WC-TPC tracks ===                                                          
      for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn )
	{
	  // === Get the TPC track ===
	  cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));
	  if (!trackWC2TPC.isValid())
	    {
	      evt.put(std::move(TrackNonStoppingVector));
	      return;
	    }
	  recob::Track const& aTrack(trackWC2TPC.ref());
	  // === Loop on all the TPC Tracks and check if the track ID is the same ===
	  // === THIS WAY IS NOT ELEGANT, BUT GETS THE JOBS DONE ===
	  for ( auto const& thisTrack : tracklist )
	    {
	      if (thisTrack->ID() == aTrack.ID())
		{
		  if (!isStoppingTrack(thisTrack,fmcal))  (*TrackNonStoppingVector).push_back(thisTrack);
		}
	    }
	}
    } else 
    {
      std::cout<<"Not checking Ass\n";
      // ### Looping over tracks and check every track if stopping###
      for ( auto const& thisTrack : tracklist )
	{ 
	  if (!isStoppingTrack(thisTrack,fmcal))  (*TrackNonStoppingVector).push_back(thisTrack);   
	}
      
    }
  
  evt.put(std::move(TrackNonStoppingVector));
}

void StoppingTracks::beginJob()
{
  // Declare checks histograms
  art::ServiceHandle<art::TFileService> tfs;

  fTotNTrack        = tfs->make<TH1F>("TotNTrack"        ,"TotNTrack       :Y:X"  ,30,-0.5,29.5);
  fStoppingTrack    = tfs->make<TH1F>("StoppingTrack"    ,"StoppingTrack   :Y:X"  ,30,-0.5,29.5);
  fDifferenceTrack  = tfs->make<TH1F>("DifferenceTrack"  ,"DifferenceTrack :Y:X"  ,30,-0.5,29.5);
  fNTrackPassedCuts = tfs->make<TH1F>("NTrackPassedCuts" ,"NTrackPassedCuts:Y:X"  ,30,-0.5,29.5);
  
  fdEdx = tfs->make<TH1F>("dEdx per each point", "dEdx per each point", 240,-10.,50.);
  fdEdxRange = tfs->make<TH2F>("dEdxVSresRange", "dEdxVSresRange for the tpc track selected", 100,0.,100., 100,0.,20.);
  fFitPar = tfs->make<TH1F>("expParLastHits fit", "expPar from lastHits fit", 100,-10.,10.);
  
}

void StoppingTracks::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel"              );
  fupstreamZPosition      = p.get< double      >("upstreamZPosition"     , 2.0   );
  fCalorimetryModuleLabel = p.get< std::string >("CalorimetryModuleLabel", "calo");
  fWC2TPCModuleLabel      = p.get< std::string >("WC2TPCModuleLabel"     , "WC2TPCtrk");
  fWCTrackLabel           = p.get< std::string >("fWCTrackLabel"         , "wctrack");
  fParticleMass           = p.get< double      >("ParticleMass"          , 139.57);//default particle: charged pion
  //Range for the p1 parameter of the fit dEdx vs RR (low RR) to select stopping particles 
  fLowLimitStop = p.get< double >("LowerLimitStoppingTrack", -0.43);
  fUpLimitStop  = p.get< double >("UpperLimitStoppingTrack", -0.35);
  
  fMinNofSpacePoints  = p.get< int >("fMinNofSpacePoints", 16);
  fCheckAssn          = p.get< bool   >("CheckAssn", false);
}

bool StoppingTracks::isStoppingTrack( art::Ptr<recob::Track> aTrack, art::FindManyP<anab::Calorimetry> fmcal )
{
  std::cout<<"FUNCTION IS STOPPING CALLED\n";
  /**
     In this function we decide if the given track is stopping or not
  */
  bool StoppingTrack = false;
  // If the calorimetry is not valid, you consider the particle stopping. 
  // You want the particle out from the interaction pool
  if (!fmcal.isValid()) return true; 
  
  // ########################################################## 
  // ### Looping over Calorimetry information for the track ###
  // ########################################################## 
  
  // ### Putting calo information for this track (i) into pointer vector ###
  std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(aTrack.key());
  
  // ### Looping over each calorimetry point (similar to SpacePoint) ###
  for (size_t j = 0; j<calos.size(); ++j)
    {
      // ### If we don't have calorimetry information for this plane skip ###
      if (!calos[j]->PlaneID().isValid) continue;
      
      // ### Grabbing this calorimetry points plane number (0 == induction, 1 == collection) ###
      int pl = calos[j]->PlaneID().Plane;
      
      // ### Skipping this point if the plane number doesn't make sense ###
      if (pl<0||pl>1) continue;
  
      // ### Recording the number of calorimetry points for this track in this plane ####
      trkhits[pl] = calos[j]->dEdx().size();
      
      // #### Recording the kinetic energy for this track in this plane ###
      trkke[pl] = calos[j]->KineticEnergy();
      
      // ###############################################
      // ### Looping over all the calorimetry points ###
      // ###############################################
      
      //if(pl == 1) std::cout << "Number of calo hits for this track in plane 1 " << calos[j]->dEdx().size() << std::endl;
      
      // double lastHitsdEdx[16]={0.};
      // double lastHitsRR[16]={0.};
      
      //double lastHitsdEdx[fMinNofSpacePoints]={0.};
      //double lastHitsRR[fMinNofSpacePoints]={0.};
      std::vector<double> lastHitsdEdx(fMinNofSpacePoints,0.);
      std::vector<double> lastHitsRR(fMinNofSpacePoints,0.);

      double ordereddEdx[1000]={0.};
      double orderedRR[1000]={0.};
      
      for (size_t k = 0; k<calos[j]->dEdx().size(); ++k)
	{
	  // ### If we go over 1000 points just skip them ###
	  if (k>=1000) continue;
	  
	  // ### Recording the dE/dX information for this calo point along the track in this plane ###
	  trkdedx[pl][k] = calos[j]->dEdx()[k];
	  
	  // ### Recording the residual range for this calo point along the track in this plane ###
	  trkrr[pl][k] = calos[j]->ResidualRange()[k];
	  
	  // ### Recording the pitch of this calo point along the track in this plane ###
	  trkpitchhit[pl][k] = calos[j]->TrkPitchVec()[k];
	  
	  //### Analyzing caloHits ONLY from collection plane - 1
	  if(pl == 1) {
	    size_t dimCalo = 0;
	    dimCalo = calos[j]->dEdx().size();
	    //Fill histos with calo info from collection plane for each track
	    fdEdx->Fill((calos[j]->dEdx()[k]));
	    fdEdxRange->Fill(calos[j]->ResidualRange()[k],(calos[j]->dEdx()[k]));
	    //In case the recorded CaloPoints are not ordered with decreasing RR: 
	    if(k < calos[j]->dEdx().size()-1){
	      if(calos[j]->ResidualRange()[k] > calos[j]->ResidualRange()[k+1]) {
		//If the previous caloHit RR is higher than the next, that's the starting point of the track	
		ordereddEdx[k]= calos[j]->dEdx()[k];
		orderedRR[k]= calos[j]->ResidualRange()[k];
	      }
	      else {
		ordereddEdx[dimCalo-k]= calos[j]->dEdx()[k];
		orderedRR[dimCalo-k]= calos[j]->ResidualRange()[k];	
	      }
	    }
	  }
	}//<---End calo points (k)
      
      if(pl == 1){
	//Actually to study the dEdx vs RR 
	//and to provide a good fit for distinguishing stopping particles, I take in account only tracks longer than around 8 cm
	if(calos[j]->dEdx().size() - fMinNofSpacePoints > 0){
	  size_t hj=0;
	  
	  hj=calos[j]->dEdx().size()-17;
	  int hjj=0;
	  while(hj > calos[j]->dEdx().size()-fMinNofSpacePoints && hj < calos[j]->dEdx().size()-1 ){
	    //lastHits are the one with lower RR
	    lastHitsdEdx[hjj]=ordereddEdx[hj];
	    lastHitsRR[hjj]=orderedRR[hj];
	    hj++;
	    hjj++;
	  }
	  
	  
	} else 
	  {
	    std::cout << "Too short track: "<< std::endl; 
	    return false;
	  }
      }
      
      int check=0;
      int h=0;
      //Doublechecking the lastHits vector has been filled
      while(h < fMinNofSpacePoints){
	if(lastHitsRR[h]==0.) {h++; check++;}
	else h++;
      }
      
      double p1=0.;
      
      if(check != fMinNofSpacePoints){
	TGraph *g1 = new TGraph(16,lastHitsRR.data(),lastHitsdEdx.data());
	TF1 *fitFcn = new TF1("fitFcn","[0]*pow(x,[1])",0.,10.);
	fitFcn->SetParameter(1,-0.4);
	g1->Fit("fitFcn");
	//In case of "invalid fit" we have to add a control <------
	p1 = fitFcn->GetParameter(1);
	//std::cout << "RR Fit Parameters "<< p1 << std::endl;
	fFitPar->Fill(p1);
      } 
      
      if(p1 < fUpLimitStop && p1 > fLowLimitStop) {
	//std::cout << "Possible stopping track "	<< std::endl;
	StoppingTrack = true;
      }
    }//<---End looping over calo points (j)	  
  
  return StoppingTrack;
}


DEFINE_ART_MODULE(StoppingTracks)

/*
bool StoppingTracks::isEscapingTrack( recob::Track aTrack )
{
 
  return false;
}
*/



