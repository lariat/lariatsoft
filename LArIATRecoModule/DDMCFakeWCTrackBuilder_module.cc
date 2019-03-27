////////////////////////////////////////////////////////////////////////
// Class:       DDMCFakeWCTrackBuilder
// Module Type: producer
// File:        DDMCFakeWCTrackBuilder_module.cc
//
// Generated at Fri Oct 16 14:58:18 2015 by Greg Pulliam using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#ifndef DDMCFAKEWCTRACKBUILDER_H
#define DDMCFAKEWCTRACKBUILDER_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include <map>
#include <iostream>
#include <fstream>
#include <math.h>       /* remainder */
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

#include <vector>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"


//ROOT Things
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATDataProducts/WCTrack.h"
#include "Utilities/DatabaseUtilityT1034.h"
#include "LArIATDataProducts/WCTrack.h"

#include <memory>
#include <utility>
#include <string>
#include <fstream>
namespace wct{
class DDMCFakeWCTrackBuilder;

class DDMCFakeWCTrackBuilder : public art::EDProducer {
public:
  explicit DDMCFakeWCTrackBuilder(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DDMCFakeWCTrackBuilder(DDMCFakeWCTrackBuilder const &) = delete;
  DDMCFakeWCTrackBuilder(DDMCFakeWCTrackBuilder &&) = delete;
  DDMCFakeWCTrackBuilder & operator = (DDMCFakeWCTrackBuilder const &) = delete;
  DDMCFakeWCTrackBuilder & operator = (DDMCFakeWCTrackBuilder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  //void beginJob(fhicl::ParameterSet const & p);
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
//  void endRun(art::Run & r) override;
//  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
//  void respondToCloseInputFile(art::FileBlock const & fb) override;
//  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
//  void respondToOpenInputFile(art::FileBlock const & fb) override;
//  void respondToOpenOutputFiles(art::FileBlock const & fb) override;					    
					    
  void plotTheTrackInformation( std::vector<double> reco_pz_list,
				std::vector<double> x_face_list,
				std::vector<double> y_face_list,
				std::vector<double> theta_list,
				std::vector<double> phi_list,
				std::vector<double> y_kink_list,
				std::vector<double> x_dist_list,
				std::vector<double> y_dist_list,
				std::vector<double> z_dist_list);
  void ResetTree();
				

private:

  // Declare member data here.

    //Offset ont he B field
    //float offset;

    std::string       fSlicerSourceLabel;


    //Hardware constants
    //int fNumber_wire_chambers;
    //int fNumber_wires_per_tdc;
    
    //int evtcounter = 0;
    //Histograms for plotting
    TH1F* fReco_Pz;
    TH1F* fY_Kink;
    TH1F* fX_Dist;
    TH1F* fY_Dist;
    TH1F* fZ_Dist;
    TH1F* fX_Face_Dist;
    TH1F* fY_Face_Dist;


  TH1F* ftrueX_Face;
  TH1F* ftrueY_Face;

  TH1F* fDeltaX_Face;
  TH1F* fDeltaY_Face;

    TH1F* fX_Face_Dist_Line;
    TH1F* fY_Face_Dist_Line;

    TH1F* fTheta_Dist;
    TH1F* fPhi_Dist;
    TH1F* fTrack_Type;
  //    std::vector<TH2F*> fRecodiff;
    TH1F* fWCDist;
    TTree* fTree;      
    //Misc
  bool fVerbose;
  // bool fPickyTracks;
  //  bool fCheckTracks;
    int WCx1,WCx2,WCx3,WCx4,WCy1,WCy2,WCy3,WCy4; //Arrays for the wires hit in each WC

  double minX =  0.0;
  double maxX = 47.0;
  double minY =-20.0;
  double maxY = 20.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 90.0;
  
  // ### Adding this to fake a WCTrack Object in MC ###
  float WC1Mult,WC2Mult,WC3Mult,WC4Mult;
  bool PickyTrackCheck;
  

};
  

  DDMCFakeWCTrackBuilder::DDMCFakeWCTrackBuilder(fhicl::ParameterSet const & p)
  {
    // Call appropriate produces<>() functions here.
    this->reconfigure(p);
    // Call appropriate produces<>() functions here.  
    produces<std::vector<ldp::WCTrack> >();
  }

void DDMCFakeWCTrackBuilder::produce(art::Event & e)
{
  //Reset the variable tree
  ResetTree();

  //Let's start with declaring the collection of WC I'm going to put on the event
  std::unique_ptr<std::vector<ldp::WCTrack> > WCTrackCol(new std::vector<ldp::WCTrack> );

  // These are a bunch of variables we need to fill
  // in order to contructed the WC Fake track. Some of them might be bogus

  std::vector<double> reco_pz_list;   // Put momentum here
  std::vector<double> reco_pz2M_list; // Put momentum here
  std::vector<double> y_kink_list;
  std::vector<double> x_dist_list;
  std::vector<double> y_dist_list;
  std::vector<double> z_dist_list;
  std::vector<double> x_face_list; // Projection on FF
  std::vector<double> y_face_list; // 
  std::vector<double> theta_list;  // Calculate your own from the MC
  std::vector<double> phi_list;    // Calculate your own from the MC
  std::vector<int> WC_vect;        // 1, 2, 3, 4 
  std::vector<float> hit_wire_vect; 
  float hit_position_vect[4][3];   // WC1 [0][xyz]
  int   WCMissed; 
  float residual;  

  //Bogus Filling (I'm not going to need them later)
  y_kink_list.push_back(-99999.);
  x_dist_list.push_back(-99999.);
  y_dist_list.push_back(-99999.);
  z_dist_list.push_back(-99999.);
  WC_vect.push_back(1);        // 1, 2, 3, 4 
  WC_vect.push_back(2);        // 1, 2, 3, 4 
  WC_vect.push_back(3);        // 1, 2, 3, 4 
  WC_vect.push_back(4);        // 1, 2, 3, 4 
  hit_wire_vect.push_back(-99999.);
  WCMissed = 0 ; 
  residual = 0.;  

  // Semi-Bogus filling (only position at WC4 needed)
  art::ServiceHandle<geo::Geometry> geom;
  //std::vector<geo::AuxDetGeo*> const & theAuxDetGeoVect = geom->AuxDetGeoVec();
  
  double centerOfDet[3] = {0,0,0};
  for( size_t iDet = 0; iDet < geom->NAuxDets() ; ++iDet ){
    geo::AuxDetGeo const& anAuxDetGeo = geom->AuxDet(iDet);
    std::string detName = anAuxDetGeo.Name();
    size_t wcnum = 999;
    if( detName == "volAuxDetSensitiveWC1") wcnum = 1;
    if( detName == "volAuxDetSensitiveWC2") wcnum = 2;
    if( detName == "volAuxDetSensitiveWC3") wcnum = 3;
    if( wcnum != 999 ){
      anAuxDetGeo.GetCenter(centerOfDet);
      hit_position_vect[wcnum-1][0] = centerOfDet[0] * CLHEP::cm;
      hit_position_vect[wcnum-1][1] = centerOfDet[1] * CLHEP::cm;
      hit_position_vect[wcnum-1][2] = centerOfDet[2] * CLHEP::cm;
    }
  }



  // Use the true MC info to construct the fake WC
  // Get the backtracker to recover true quantities
  //art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  
  TVector3 truePositionAtFF;
  TVector3 position;
  TVector3 momentum;
  double trueTheta;
  double truePhi;  
  double XProj;
  double YProj;
  
  
  
  
  sim::ParticleList::const_iterator itPart = plist.begin(),pend = plist.end();
  
  for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart)
     {
      const simb::MCParticle* pPart = (itPart++)->second;
      if (!pPart) {
      throw art::Exception(art::errors::LogicError)
      << "GEANT particle #" << iPart << " returned a null pointer";
      }
      
  std::string proc = pPart->Process();
  
  // Get the true particle and its process, skip whatever is not primary
  /*std::vector<const simb::MCParticle* > MCParticle;
  
  for(size_t p = 0; p < plist.size(); ++p) 
    {
    // ### Filling the vector with MC Particles ###
    MCParticle.push_back(plist.Particle(p));
    
    }
    
  for( unsigned int i = 0; i< MCParticle.size(); ++i)
      {  
      
      //auto mcPart = plist.Particle(p);
      std::string proc = MCParticle[i]->Process();*/
      if ( !(proc.find("primary") != std::string::npos) ) continue;
      if (fVerbose) std::cout<<"\n\nPrimary \n";
      // Get the True Trajectory point
      simb::MCTrajectory truetraj = pPart->Trajectory();
      auto firstTrjPtInTPC = truetraj.begin();
      
      if (fVerbose) std::cout<<"Trajectory size "<<truetraj.size()<<" firstTrjPtInTPC->first.Z() = "<<(firstTrjPtInTPC->first).Z()<<std::endl;
      if ( (firstTrjPtInTPC->first).Z() >  -50. ) continue; 
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  if (pos.Z() < minZ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    firstTrjPtInTPC = t;
	    break;
	  }
	}// End search for first point in TPC
      truePositionAtFF = (firstTrjPtInTPC->first).Vect();

      auto firstTrjPt      = truetraj.begin();
      position = (firstTrjPt->first).Vect();
      momentum = (firstTrjPt->second).Vect();
      
    }



  double Xline = momentum.X()*(-1.*position.Z()/momentum.Z()) + position.X();
  double Yline = momentum.Y()*(-1.*position.Z()/momentum.Z()) + position.Y();

  fX_Face_Dist_Line->Fill(Xline);
  fY_Face_Dist_Line->Fill(Yline);

  // Momentum direction
  truePhi   =  TMath::ATan2(momentum.Y(),momentum.X());
  trueTheta =  TMath::ACos(momentum.Z()/momentum.Mag());
 

  // Let's calculate the projected point at the TPC FF
  // We know the following: 
  // initial position (x1,y1,z1)
  // z final position z2 = 0 (aka TPC FF)
  // Angles between these 2 vectors
  // We use 
  // x2 - x1 = dx = d sin(theta)cos(phi)
  // y2 - y1 = dy = d sin(theta)sin(phi)
  // z2 - z1 = dz = d cos(theta)  ==> d = -z1/cos(theta)
  auto d = -1.*position.Z()/TMath::Cos(trueTheta) ;
  XProj  = position.X() + d * TMath::Sin(trueTheta) * TMath::Cos(truePhi) ;
  YProj  = position.Y() + d * TMath::Sin(trueTheta) * TMath::Sin(truePhi) ;
  
  if (truePositionAtFF.Z() > -10)
    {
      ftrueX_Face->Fill(truePositionAtFF.X());
      ftrueY_Face->Fill(truePositionAtFF.Y());

      fDeltaX_Face->Fill(truePositionAtFF.X() - XProj);
      fDeltaY_Face->Fill(truePositionAtFF.Y() - YProj);
    } 

  if (fVerbose) 
    {
      std::cout<<"X "<<position.X()<<" Y "<<position.Y() << " Z "<<position.Z()<<std::endl;
      std::cout<<"XProj "<<XProj<<" YProj "<<YProj<<" d "<<d <<std::endl;
      std::cout<<"trueTheta "<<trueTheta*57.2958<<" truePhi "<<truePhi*57.2958 <<std::endl;
    }
  //Fill the meaningful quantities
  reco_pz_list  .push_back(1000*momentum.Z());   // Put momentum here
  reco_pz2M_list.push_back(1000*momentum.Z());   // Put momentum here
  x_face_list.push_back(XProj); // Projection on FF
  y_face_list.push_back(YProj); // 
  theta_list .push_back(trueTheta);  // Calculate your own from the MC
  phi_list   .push_back(truePhi);    // Calculate your own from the MC
  hit_position_vect[3][0] = position.X();
  hit_position_vect[3][1] = position.Y();
  hit_position_vect[3][2] = position.Z();

  WCx1=hit_position_vect[0][0];
  WCx2=hit_position_vect[1][0];
  WCx3=hit_position_vect[2][0];
  WCx4=hit_position_vect[3][0];
  WCy1=hit_position_vect[0][1];
  WCy2=hit_position_vect[1][1];
  WCy3=hit_position_vect[2][1];
  WCy4=hit_position_vect[3][1];
  
  int iTrack = 0;  
			 
   // ### Going to set the cheated Monte Carlo Wire Chamber Track ###
   // ### to be a "Pick Track" since we are putting the position  ###
   // ###                        in by hand			  ###
   
   WC1Mult = 1;
   WC2Mult = 2;
   WC3Mult = 3;
   WC4Mult = 4;
   PickyTrackCheck = true;		 
   ldp::WCTrack the_track(   reco_pz_list[iTrack],
                             reco_pz2M_list[iTrack],
                             y_kink_list[iTrack],
                             x_dist_list[iTrack],
                             y_dist_list[iTrack],
                             z_dist_list[iTrack],
                             x_face_list[iTrack],
                             y_face_list[iTrack],
                             theta_list[iTrack],
                             phi_list[iTrack],
                             WC_vect,
                             hit_wire_vect,
                             hit_position_vect,
                             WCMissed,
                             residual,                             
                             WC1Mult,
                             WC2Mult,
                             WC3Mult,
                             WC4Mult,
                             WC1Mult,
                             WC2Mult,
                             WC3Mult,
                             WC4Mult,			     
                             PickyTrackCheck);

  (*WCTrackCol).push_back( the_track );
 
 
  //Plot the reconstructed momentum, y_kink, and delta X, Y, Z in histos
  plotTheTrackInformation(reco_pz_list,
			  x_face_list,
			  y_face_list,
			  theta_list,
			  phi_list,
			  y_kink_list,
			  x_dist_list,
			  y_dist_list,
			  z_dist_list);



  //Put objects into event (root file)
  e.put(std::move(WCTrackCol)); 
  fTree->Fill(); 
  if (fVerbose) std::cout<<"\n\n"; 
}
//==================================================================================================
void DDMCFakeWCTrackBuilder::ResetTree()
{
  
  WCx1=-99999;
  WCx2=-99999;
  WCx3=-99999;
  WCx4=-99999;
  WCy1=-99999;
  WCy2=-99999;
  WCy3=-99999;
  WCy4=-99999;
}

//===================================================================================
void DDMCFakeWCTrackBuilder::plotTheTrackInformation( std::vector<double> reco_pz_list,
						      std::vector<double> x_face_list,
						      std::vector<double> y_face_list,
						      std::vector<double> theta_list,
						      std::vector<double> phi_list,
						      std::vector<double> y_kink_list,
						      std::vector<double> x_dist_list,
						      std::vector<double> y_dist_list,
						      std::vector<double> z_dist_list)
{
  //Loop through the tracks and fill
  for( size_t iTrack = 0; iTrack < reco_pz_list.size(); ++iTrack ){
    fReco_Pz->Fill(reco_pz_list.at(iTrack));
    fY_Kink->Fill(y_kink_list.at(iTrack));
    fX_Dist->Fill(x_dist_list.at(iTrack));
    fY_Dist->Fill(y_dist_list.at(iTrack));
    fZ_Dist->Fill(z_dist_list.at(iTrack));
    fX_Face_Dist->Fill(x_face_list.at(iTrack));
    fY_Face_Dist->Fill(y_face_list.at(iTrack));
    fTheta_Dist->Fill(theta_list.at(iTrack));
    fPhi_Dist->Fill(phi_list.at(iTrack));
  }  
}

//=========================================================================================
void DDMCFakeWCTrackBuilder::beginJob()//fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fWCDist= tfs->make<TH1F>("WCCond","WC Conditions",7,0,7); 
  fTree=tfs->make<TTree>("WCTree","WCTree");
  fTree->Branch("WCx1",&WCx1,"WCx1/I");
  fTree->Branch("WCx2",&WCx2,"WCx2/I");
  fTree->Branch("WCx3",&WCx3,"WCx3/I");
  fTree->Branch("WCx4",&WCx4,"WCx4/I");
  fTree->Branch("WCy1",&WCy1,"WCy1/I");
  fTree->Branch("WCy2",&WCy2,"WCy2/I");
  fTree->Branch("WCy3",&WCy3,"WCy3/I");
  fTree->Branch("WCy4",&WCy4,"WCy4/I");

  
  //Hists that should probably stay for the production run.    
  fReco_Pz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum in XZ plane", 180, 0, 1800);
  fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",200,-10*3.1415926/180,10*3.141592654/180);
  fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",1200,-60,1260);
  fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",1200,-600,600);
  fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",1200,-60,1260);
  fX_Face_Dist = tfs->make<TH1F>("X_Face","X Location of Track's TPC Entry (mm)",1600,-200,1400);
  fY_Face_Dist = tfs->make<TH1F>("Y_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);
  
  fX_Face_Dist_Line = tfs->make<TH1F>("X_Face_Line","X Location of Track's TPC Entry (mm)",1600,-200,1400);
  fY_Face_Dist_Line = tfs->make<TH1F>("Y_Face_Line","Y Location of Track's TPC Entry (mm)",800,-400,400);

  ftrueX_Face = tfs->make<TH1F>("trueX_Face","X Location of Track's TPC Entry (mm)",1600,-200,1400);
  ftrueY_Face = tfs->make<TH1F>("trueY_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);


  fDeltaX_Face = tfs->make<TH1F>("DeltaX_Face","Delta X True - Proj FF (cm)",1200,-30,30);
  fDeltaY_Face = tfs->make<TH1F>("DeltaY_Face","Delta Y True = Proj FF (cm)",1200,-30,30);

  fTheta_Dist = tfs->make<TH1F>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",400,-.4,0.4);
  fPhi_Dist = tfs->make<TH1F>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",2000,-6.28318,6.28318);                   
  fReco_Pz->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
  fReco_Pz->GetYaxis()->SetTitle("Tracks per 10 MeV/c");
  fY_Kink->GetXaxis()->SetTitle("Reconstructed y_kink (radians)");
  fY_Kink->GetYaxis()->SetTitle("Tracks per 0.000872 radians");
  fX_Dist->GetXaxis()->SetTitle("X distance between US and DS track ends");
  fX_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
  fY_Dist->GetXaxis()->SetTitle("Y distance between US and DS track ends");
  fY_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
  fZ_Dist->GetXaxis()->SetTitle("Z distance between US and DS track ends");
  fZ_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
  fX_Face_Dist->GetXaxis()->SetTitle("X (mm)");
  fX_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
  fY_Face_Dist->GetXaxis()->SetTitle("Y (mm)");
  fY_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
  fTheta_Dist->GetXaxis()->SetTitle("Theta (radians)");
  fTheta_Dist->GetYaxis()->SetTitle("Tracks per .002 radians");
  fPhi_Dist->GetXaxis()->SetTitle("Phi (radians)");
  fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.00628 radians");
  
  fTrack_Type = tfs->make<TH1F>("TrackType","WCTrack conditions: 1=missHit,2=uniqueHits,3=lonelyHit,4=socialHits",4,0,4);
  fTrack_Type->GetYaxis()->SetTitle("# Events");
  fTrack_Type->GetXaxis()->SetTitle("Track Conditions");
}

void DDMCFakeWCTrackBuilder::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void DDMCFakeWCTrackBuilder::beginSubRun(art::SubRun & sr)
{
  /*
   std::cout<<"Module B: "<<fDDMCFakeWCTrackBuilderAlg.fMCMagneticField<<std::endl;
    // If the field override is not set, get the actual magnetic field value for the alg
    if( fDDMCFakeWCTrackBuilderAlg.fMCMagneticField == 0 ){ 
    	fDDMCFakeWCTrackBuilderAlg.loadXMLDatabaseTableForBField( sr.run(), sr.subRun() );}
  */
}

void DDMCFakeWCTrackBuilder::endJob()
{
  // Implementation of optional member function here.
}

// void DDMCFakeWCTrackBuilder::endRun(art::Run & r)
// {
//   // Implementation of optional member function here.
// }
// 
// void DDMCFakeWCTrackBuilder::endSubRun(art::SubRun & sr)
// {
//   // Implementation of optional member function here.
// }

void DDMCFakeWCTrackBuilder::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  //fNumber_wire_chambers = p.get<int>("NWC"); //4;  
  //fNumber_wires_per_tdc = p.get<int>("NWperTDC"); //64;
  fVerbose = p.get<bool>("Verbose", false);
  //fSlicerSourceLabel = p.get<std::string>("SourceLabel");
  //std::cout<<"Label WC: "<<fSlicerSourceLabel<<std::endl;
  //// fPickyTracks=p.get<bool>("PickyTracks");
  //fCheckTracks=p.get<bool>("CheckTracks");
  //offset = p.get<float>("BFieldOffset");
  
  
}
// 
// void DDMCFakeWCTrackBuilder::respondToCloseInputFile(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void DDMCFakeWCTrackBuilder::respondToCloseOutputFiles(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void DDMCFakeWCTrackBuilder::respondToOpenInputFile(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void DDMCFakeWCTrackBuilder::respondToOpenOutputFiles(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }

DEFINE_ART_MODULE(DDMCFakeWCTrackBuilder)
}//end namespace
#endif //DDMCFakeWCTrackBuilder_H
