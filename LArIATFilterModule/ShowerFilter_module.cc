////////////////////////////////////////////////////////////////////////
// Class:       ShowerFilter
// Module Type: filter
// File:        ShowerFilter_module.cc
//
// Generated at Wed Feb 25 17:12:56 2016 by Irene Nutini using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

//This is a very *rough and preliminary* version for a Shower Filter

//It was written for the PionAnalysis  (reco and analysis in lariatsoft_v04_34_00)
//since there was still a 15% of electron induced showers contamination in our sample (especially for negative pions) after all the Selection Cuts and TPC/WC track matching  

// The filter looks at how many tracks are reconstructed for the event
// Then it counts how many of them are "short tracks"; I define as "short", tracks whose length is less than 5 cm (parameter that can be changed in the fcl file).
// We look for short tracks since in general with LineCluster and Trackers (no shower algorithms) the shower events are reconstructed as many short tracks
// If in the event there are less than 2 short tracks (parameter that can be changed in the fcl file) this is considered as a "track-like" event,
// otherwise it's considered as a "shower-like" event
// "shower-like" events are rejected

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
#include "LArIATDataProducts/WCTrack.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <TH2F.h>
#include "art_root_io/TFileService.h"
#include "cetlib/maybe_ref.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "LArIATDataProducts/WCTrack.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

const int kMaxTrack      = 1000; // unused  //maximum number of tracks

class ShowerFilter;

class ShowerFilter : public art::EDFilter {
public:
  explicit ShowerFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerFilter(ShowerFilter const &) = delete;
  ShowerFilter(ShowerFilter &&) = delete;
  ShowerFilter & operator = (ShowerFilter const &) = delete;
  ShowerFilter & operator = (ShowerFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p);
  bool CutCheck(std::vector<double> lengths, double lcut, double ncut);

private:
  int initialEvents = 0;
  int finalEvents = 0;
  // Declare member data here.
std::string fTrackModuleLabel;
std::string fWCTrackLabel;
std::string fWC2TPCModuleLabel;
int fnShortTks;
double fShortTkLength;
bool bData, bDoneWC2TPC;
bool bFilterOrSelect;
bool KeepEvent;
double fPassingRatioScale;
double fConeLength;
double fStartRadius;
double fEndRadius;
double fWC2TPCCircularXCenter;
double fWC2TPCCircularYCenter;
double fWC2TPCCircularRadius;
double fWC2TPCAlphaCut;
double backtrackcount=0;
double reversetrackcount=0; 
    int ntrks;
    TVector3 MatchedTrkStartDir;
    TVector3 MatchedTrkStartPos;
    TVector3 MatchedTrkEndPos;
    double trkLength[kMaxTrack]={0.};
  TH1F*  hnTracks;
    TH1F* fXface;
  TH1F* fYface;
  TH1F* fDistface;
  TH1F* fAlpha;
  TH1F* fRejectionModes;
  TH1F* fMomShowerPass;
  TH1F* fMomWCPass;
  TH1F* fMom;
  TH2F* fCheckPlot;
  TH1F* fCandidateTrackLength;
  TH2F* fXYface;
TVector3 TruthMom;
    //double n0=0;
//Histograms
};


ShowerFilter::ShowerFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool ShowerFilter::filter(art::Event & e)
{
  initialEvents++;
TruthMom.SetX(0.0);
TruthMom.SetY(0.0);
TruthMom.SetZ(0.0);
  // Implementation of required member function here.
   
  // #####################################
   // ### Getting the Track Information ###	
   // #####################################
   art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Track> > Tracklist; //<---Define tracklist as a pointer to recob::tracks	
   // === Filling the tracklist from the tracklistHandle ===	
   if (e.getByLabel(fTrackModuleLabel,trackListHandle))
      {art::fill_ptr_vector(Tracklist, trackListHandle);}
      
  //std::cout << "There are: " << tracklist.size() << " tracks produced for this event" << std::endl;   


  double xface = 0.0;
  double yface = 0.0;
  double mom=0.0;
  //double momo = 0.0;
  std::vector<double> WC_vec = {-1.0, -1.0, -1.0};
    if(bData && !bDoneWC2TPC)
    {
      std::cout<<"You've specified this event as data and no WC2TPC match has been made. You cannot run this filter on data without the WC2TPC Match code already run! Run WC2TPC Match then try again. Returning False!"<<std::endl;
      return false;
    }
  
    if(bData || bDoneWC2TPC) {
//    if(bVerbose) { std::cout << " Data event. " << std::endl; }

    // WC info
    art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
    std::vector<art::Ptr<ldp::WCTrack> > wctrack;
   
    if(e.getByLabel(fWCTrackLabel, wctrackHandle))
      art::fill_ptr_vector(wctrack, wctrackHandle);
    art::FindOneP<recob::Track> WCTrackAssn(wctrackHandle,   e, fWC2TPCModuleLabel);
    if(WCTrackAssn.isValid())
    {
      for (unsigned int indexAssn = 0; indexAssn < WCTrackAssn.size(); ++indexAssn ) 
      {
        // =========================                                                                                       
        // === Get the TPC track ===
	// =========================                                                                      
	cet::maybe_ref<recob::Track const> trackWC2TPC(*WCTrackAssn.at(indexAssn));
	if (!trackWC2TPC){return false;}
	if (trackWC2TPC)
	{
	  recob::Track MatchedTrack=*WCTrackAssn.at(indexAssn);
	  recob::Track::Vector_t tmpTrkStartDir=MatchedTrack.DirectionAtPoint(0);
    MatchedTrkStartDir.SetXYZ(tmpTrkStartDir.X(),tmpTrkStartDir.Y(),tmpTrkStartDir.Z());
// The Z component of the direction is negative sometimes, indicating the track is in reverse. If it's extremely backgoing (Z<-.95), flip the tracjectory to make it point forward. Also keep a count of how often this happens.
          if(MatchedTrkStartDir[2]<-0.95)
	  {
	    MatchedTrkStartDir=-1*MatchedTrkStartDir; ++backtrackcount;  
	  }
//also, I bet that means the track start and end position are reversed, so we'll make sure the track starts upstream by flipping the start and end point if the end point is upstream.
//Again keep a counter and hope it increases when the direction counter increases.
	
   MatchedTrkStartPos.SetXYZ(MatchedTrack.Vertex().X(), MatchedTrack.Vertex().Y(), MatchedTrack.Vertex().Z());
	 MatchedTrkEndPos  .SetXYZ(MatchedTrack.End().X()   , MatchedTrack.End().Y(),    MatchedTrack.End().Z());
	  if(MatchedTrkEndPos[2]<MatchedTrkStartPos[2]){MatchedTrkStartPos=MatchedTrkEndPos;  ++reversetrackcount;}
	  // std::cout<<"For this event, the Matched Track Starts at: ("<<MatchedTrkStartPos[0]<<", "<<MatchedTrkStartPos[1]<<", "<<MatchedTrkStartPos[2]<<")."<<std::endl;
	} //if trackWC2TPC
      } //indexAssn loop
    } //isValid
    //std::cout<<"Track with reversed direction: "<<backtrackcount<<". Tracks with reversed start/end points: "<<reversetrackcount<<"."<<std::endl;
    // requires single WC track
   // if(bVerbose) { std::cout << "wcsize = " << wctrack.size() << std::endl; }
    if(wctrack.size() != 1) { fRejectionModes->Fill(1.0); return false; }

    xface = wctrack[0]->XYFace(0);
    yface = wctrack[0]->XYFace(1);

    mom = wctrack[0]->Momentum();
    fMom->Fill(mom);

    WC_vec[0] = wctrack[0]->HitPosition(3, 0) - wctrack[0]->HitPosition(2, 0);
    WC_vec[1] = wctrack[0]->HitPosition(3, 1) - wctrack[0]->HitPosition(2, 1);
    WC_vec[2] = wctrack[0]->HitPosition(3, 2) - wctrack[0]->HitPosition(2, 2);

   // if(bVerbose) { std::cout << "Data WC projection (" << xface << ", " << yface << ")" << std::endl; }

  } 
  if(!bData && !bDoneWC2TPC){ 
  std::cout<<"This is an MC event without a WC2TPC Match already created. I'll try to make the WC2TPC match on the fly using the initial trajectory of the primary as a WCTrack stand-in."<<std::endl; 
    // ParticleInventoryService service ===
    //art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();

    std::vector<const simb::MCParticle* > geant_part;     
    for(size_t p = 0; p < plist.size(); ++p) { geant_part.push_back(plist.Particle(p)); }


    //fPassingEvents->Fill(0); // start of all 


    bool found = false;
    for(unsigned int i = 0; i < geant_part.size(); ++i ) {
      if(geant_part[i]->Process() != "primary") { continue; }

      // <x, y, z> = <sx, sy, sz> * t + <startx, starty, startz>
      // (z - startz) / pz = t
      // x = px * t + startx
      // y = px * t + srarty
      
      // intercept time

      double t = (0.0 - geant_part[i]->Vz()) / geant_part[i]->Pz(); 

      xface = geant_part[i]->Px() * t + geant_part[i]->Vx(); //declare
      yface = geant_part[i]->Py() * t + geant_part[i]->Vy(); //declare
    
      TVector3 wcVec;
      wcVec.SetX(geant_part[i]->Px());
      wcVec.SetY(geant_part[i]->Py());
      wcVec.SetZ(geant_part[i]->Pz());

      //momo = wcVec.Mag()*1000.; //declare
      //fMomo->Fill(momo); //declare

      WC_vec[0] = wcVec[0]; //declare
      WC_vec[1] = wcVec[1];
      WC_vec[2] = wcVec[2];
      TruthMom=wcVec;

      found = true;
    }

    if(!found) { 
      //if(bVerbose) { std::cout << "Could not find an MC WC-TPC match!" << std::endl; } 
      fRejectionModes->Fill(1.0); 
      return false; 
    }
}

if(!bData && !bDoneWC2TPC){ 
  
  int nMatched = 0.0;

  for(size_t i = 0; i < Tracklist.size(); i++) {

    //if(bVerbose) { 
     // std::cout << "Start of track: " << Tracklist[i]->Start().X() << " " << Tracklist[i]->Start().Y() << " " << Tracklist[i]->Start().Z() << std::endl;
     // std::cout << "End of track: "   << Tracklist[i]->End().X()   << " " << Tracklist[i]->End().Y()   << " " << Tracklist[i]->End().Z()   << std::endl;
    //}
    
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

    fXface->Fill(Tracklist[i]->Start().X() - xface); //declare
    fYface->Fill(Tracklist[i]->Start().Y() - yface);
    fXYface->Fill(Tracklist[i]->Start().X() - xface,Tracklist[i]->Start().Y() - yface);
    fDistface->Fill(TMath::Sqrt(TMath::Power(Tracklist[i]->Start().X() - xface, 2) + //declare
				TMath::Power(Tracklist[i]->Start().Y() - yface, 2)));
    fAlpha->Fill(acos(tpcz.Dot(wcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159)); //declare


    // Track spatial match using wc2tpc filer values dx=[-2,6]cm dy=[-3,6]cm

   
    // Alpha
    //std::cout << "alpha =" << acos(tpcz.Dot(wcz) / (wcz.Mag() * tpcz.Mag())) * (180.0 / 3.14159) << std::endl;

    //if(acos(wcz.Dot(tpcz) / (wcz.Mag() * tpcz.Mag())) >.2) { continue; }
    double XDist=Tracklist[i]->Start().X() - xface-fWC2TPCCircularXCenter;
    double YDist=Tracklist[i]->Start().Y() - yface-fWC2TPCCircularYCenter;
    
    if((XDist*XDist+YDist*YDist)>(fWC2TPCCircularRadius*fWC2TPCCircularRadius) || (acos(wcz.Dot(tpcz) / (wcz.Mag() * tpcz.Mag())) > fWC2TPCAlphaCut)){continue;}
    //std::cout << "space match " << std::endl;
    //closest_track_index = i;
    
    nMatched += 1;
    
    //Get the vectors necessary for the cone cut.
    MatchedTrkStartDir=Tracklist[i]->DirectionAtPoint<TVector3>(0);
    //recob::Track::Point_t tmp = Tracklist[i]->DirectionAtPoint<TVector3>(0);
    //MatchedTrkStartDir.SetXYZ( tmp.X(), tmp.Y(), tmp.Z() );

// The Z component of the direction is negative sometimes, indicating the track is in reverse. If it's extremely backgoing (Z<-.95), flip the tracjectory to make it point forward. Also keep a count of how often this happens.
    if(MatchedTrkStartDir[2]<-0.95)
    {
      MatchedTrkStartDir=-1*MatchedTrkStartDir; ++backtrackcount;  
    }
//also, I bet that means the track start and end position are reversed, so we'll make sure the track starts upstream by flipping the start and end point if the end point is upstream.
//Again keep a counter and hope it increases when the direction counter increases.
    MatchedTrkStartPos=Tracklist[i]->Vertex<TVector3>();
    MatchedTrkEndPos=Tracklist[i]->End<TVector3>();
    if(MatchedTrkEndPos[2]<MatchedTrkStartPos[2]){MatchedTrkStartPos=MatchedTrkEndPos;  ++reversetrackcount;}
  }
  
  if(nMatched != 1) {
    //if(bVerbose) { std::cout << "Found " <<nMatched<< " tracks. There can only be one!"<<std::endl; }
    fRejectionModes->Fill(nMatched);
    return false; 
  }
}
  fMomWCPass->Fill(mom);
   int shortTrack=0;
  // int longTrack=0;
   std::vector<double> trkLengthvect;
   std::vector<double> CandidateTracks;
   ntrks = Tracklist.size();
   for (size_t p = 0; p<Tracklist.size(); ++p){
     trkLengthvect.push_back(Tracklist[p]->Length()); //Get the length of the track   
//Time to use the cone method to find whether this track should matter or not.  First get the start and end position of this track, just in case the track is in reverse.

     TVector3 TrkStart=Tracklist[p]->Vertex<TVector3>();
     TVector3 TrkEnd=Tracklist[p]->End<TVector3>();
     
//Though the cone doesn't have to start at its geometric vertex , finding where that vertex would be is helpful for some trig. As the axis of the cone is along the wc2tpc track, the vertex is found projecting backwards along the wc2tpc track, using cone
//parameters to define how far back to project.

     double L_project=fStartRadius*fConeLength/(fEndRadius-fStartRadius);
     TVector3 ConeVertex=MatchedTrkStartPos-L_project*MatchedTrkStartDir;
//For the start and end pos of this track, we need the unit vector from the cone vertex to those points.
        
     TVector3 VertextoStart=(TrkStart-ConeVertex).Unit();
     TVector3 VertextoEnd=(TrkEnd-ConeVertex).Unit();
// The opening angle of this cone is a function of the cone parameters given in fcl. If this track starts/ends in the cone, then the arccos of the dot product of these unit vectors with the cone axis (MatchedTrackStartDir) must be less than that opening angle
     double ConeOpeningAngle=atan((fEndRadius-fStartRadius)/fConeLength); //in radians and should be positive (endradius>startradius)
     double StartAngleToAxis=acos(VertextoStart.Dot(MatchedTrkStartDir));
     
// If this track is the WC Matched Track, then the dot product above is 1, and arccos is nan. So we're gonna hygiene check that if the dot product is >.99999, we set the angle to 0;   
     if(VertextoStart.Dot(MatchedTrkStartDir)>.99999){StartAngleToAxis=0;}
     double EndAngleToAxis=acos(VertextoEnd.Dot(MatchedTrkStartDir));
//If the Track is perfectly straight and StartAngle was nan, then the end angle could be nan as well. So we're gonna be safe and do the same check.    
     if(VertextoEnd.Dot(MatchedTrkStartDir)>.99999){EndAngleToAxis=0; }
     //std::cout<<"Cone Opening Angle: "<<ConeOpeningAngle<<". Start track angle: "<<StartAngleToAxis<<". End track angle: "<<EndAngleToAxis<<std::endl;
     bool trackanglecontained=false;
     bool trackstartanglecontained=false;
     bool trackendanglecontained=false;
     
     //Check if the start or end of the track is contained. If either, then the track is a candidate to be included in short track shower cut. Using separate bools because I want to know which end (or both) to check for volume containment in a few lines.
     if(StartAngleToAxis >= 0 && StartAngleToAxis < ConeOpeningAngle ){trackstartanglecontained=true;}
     if(EndAngleToAxis >= 0 && EndAngleToAxis < ConeOpeningAngle){trackendanglecontained=true;}
     if(trackstartanglecontained || trackendanglecontained)
     {
       //std::cout<<"!!!!!!!!!!!Track Angle Contained!!!!!!!!"<<std::endl; 
       trackanglecontained=true;
     }
//We've checked the track is contained by angle, but we haven't checked whether the track starts/ends within the z boundaries of the cone, defined by the opening angle and the fStartRadius and fEndRadius. 
//Given we're measuring from the vertex of the cone, we have to add back in the distance from the WC2TPC matched Track from the vertex (L_project) to find the minimum distance along the axis. The max is fConeLength further.
     double zmin=L_project;
     double zmax=L_project+fConeLength;
     bool trackstartvolumecontained=false;         
     bool trackendvolumecontained=false; 
     bool trackvolumecontained=false; 
     if(trackstartanglecontained)
     {
       double trackdist=(TrkStart-ConeVertex).Mag()*cos(StartAngleToAxis);
       //std::cout<<"Start track dist: "<<trackdist<<std::endl;
       if (trackdist>=zmin && trackdist<=zmax)
       {
         //std::cout<<"Track Start Vol Contained"<<std::endl;
         trackstartvolumecontained=true;
       }
     }
     if(trackendanglecontained)
     {
       double trackdist=(TrkEnd-ConeVertex).Mag()*cos(EndAngleToAxis);
       //std::cout<<"End track dist: "<<trackdist<<std::endl;
       if (trackdist>=zmin && trackdist<=zmax)
       {
         //std::cout<<"Track End Vol Contained"<<std::endl; 
         trackendvolumecontained=true;
       }
     }     
     if(trackstartvolumecontained || trackendvolumecontained)
     {
       //std::cout<<"!!!!!!!!!!!Track Vol Contained!!!!!!!!"<<std::endl; 
       trackvolumecontained=true;
     }
//So if the track is within the opening angle of the cone and is within the bounds of the truncated cone we've defined, then the length of track should be included in the list of tracks to be used to decide if the event passes the shower cut.
      if(trackvolumecontained && trackanglecontained){CandidateTracks.push_back(Tracklist[p]->Length()); fCandidateTrackLength->Fill(Tracklist[p]->Length());}     
   }
   //std::cout<<"Starting with "<<Tracklist.size()<<" tracks, after the cone cut, we're down to "<<CandidateTracks.size()<<" tracks."<<std::endl;
   
   for(double lengthcut=1; lengthcut<11;++lengthcut)
   {
     for(double trkcut=1; trkcut<20; ++trkcut)
     {
       bool isgood=CutCheck(CandidateTracks,lengthcut,trkcut);
       if(isgood){fCheckPlot->Fill(lengthcut,trkcut,1./fPassingRatioScale);}
     }
   }  
   for (size_t p = 0; p<CandidateTracks.size(); ++p){

   if(CandidateTracks[p] < fShortTkLength) shortTrack++;
   
   }
   hnTracks->Fill(shortTrack);
   if(bFilterOrSelect)
   {
     if(shortTrack > fnShortTks) 
     {
       KeepEvent=false;
     }   
     else
     {
       if(bData){fMomShowerPass->Fill(mom);}
       if(!bData){fMomShowerPass->Fill(1000*TruthMom[2]);}
       KeepEvent=true;
   
     }
   }
   if(!bFilterOrSelect)
   {
     if(shortTrack > fnShortTks) 
     {
       if(bData){fMomShowerPass->Fill(mom);}
       if(!bData){fMomShowerPass->Fill(1000*TruthMom[2]);}
       KeepEvent=true;
     }   
     else
     {
       KeepEvent=false;   
     }
   }
   
   if (KeepEvent)   finalEvents++;
   return KeepEvent;       
}
bool ShowerFilter::CutCheck(std::vector<double> lengths, double lcut, double ncut)
{
  double ntemp=0;
  for(size_t ntrk=0; ntrk<lengths.size(); ++ntrk)
  {
    if(lengths[ntrk]<lcut){++ntemp;}    
  }
  if(ntemp>ncut){return false;}
  else{return true;}  
}
void ShowerFilter::beginJob()
{
  initialEvents = 0;
  finalEvents   = 0;
  // Implementation of optional member function here.
   art::ServiceHandle<art::TFileService> tfs;
   hnTracks=tfs->make<TH1F>("NShortTracks","NShortTracks",40,0,40);
  fXface = tfs->make<TH1F>("XFace","XFace",100, -20.0, 20.0);
  fYface = tfs->make<TH1F>("YFace","YFace",100, -20.0, 20.0);
  fDistface = tfs->make<TH1F>("DistFace","DistFace",100, 0.0, 50.0);
  fAlpha = tfs->make<TH1F>("Alpha","Alpha",100, 0.0, 50.0);
  fRejectionModes = tfs->make<TH1F>("RejectionModes","RejectionModes",10, 0.0, 10.0);
  fMomShowerPass=tfs->make<TH1F>("ShowerPassingMom","ShowerPassingMom",200,0,2000);
  fMomWCPass=tfs->make<TH1F>("WCTPCPassingMom","WCTPCPassingMom",200,0,2000);
  fMom=tfs->make<TH1F>("StartingMom","StartingMom",200,0,2000);
  fCheckPlot=tfs->make<TH2F>("NtrkLength","Number of Tracks Vs Length",12,0,12,20,0,20);
  fCandidateTrackLength=tfs->make<TH1F>("CandidateTrackLength","CandidateTrackLength",200,0,100);
  fXYface=tfs->make<TH2F>("WC2TPCCandidates","WC2TPCCandidates",100,-20,20,100,-20,20);
}

void ShowerFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here
  fTrackModuleLabel = p.get< std::string  >("TrackModuleLabel","pmtracktc");
  fWCTrackLabel = p.get< std::string  >("WCTrackLabel","wctrack");
  fShortTkLength = p.get< double >     ("ShortTkLength", 10.0);
  fnShortTks   = p.get< int >        ("nShortTks"  ,    5);
  bData=p.get<bool>("isData");
  fWC2TPCModuleLabel = p.get< std::string  >("WC2TPCModuleLabel", "WC2TPCtrk");
  fPassingRatioScale=p.get<double>("ScalingFactor",1);
  fConeLength=p.get<double>("ShowerConeLength",14);
  fStartRadius=p.get<double>("ShowerConeStartRadius",4);
  fEndRadius=p.get<double>("ShowerConeEndRadius",10);
  bFilterOrSelect=p.get<bool>("FilterShowers",true);  //If True, then remove events with too many short tracks. If False, keep only events with many short tracks.
  fWC2TPCCircularXCenter=p.get<double>("CircularCutXCenter",1.6);
  fWC2TPCCircularYCenter=p.get<double>("CircularCutYCenter",-.17);
  fWC2TPCCircularRadius=p.get<double>("CircularCutRadius",3.5);
  fWC2TPCAlphaCut=p.get<double>("AlphaCut",0.2);
  bDoneWC2TPC=p.get<bool>("DoneWC2TPCAlready");
}

void ShowerFilter::endJob()
{
  std::cout<<"&&&--&&&  Shower Filter Initial Events: "<<initialEvents<< " finalEvents " << finalEvents <<std::endl;
}
DEFINE_ART_MODULE(ShowerFilter)
