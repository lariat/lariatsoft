////////////////////////////////////////////////////////////////////////
// Class:       ELossBeforeTPC
// Module Type: analyzer
// File:        ELossBeforeTPC_module.cc
//
// Generated at Tue June 29 11:21:46 2017 by Elena Gramellini
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################

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
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "LArIATDataProducts/WCTrack.h"

//#include "RawData/ExternalTrigger.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "LArIATDataProducts/WCTrack.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "Math/Functor.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "Fit/Fitter.h"

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

namespace lariat 
{
  class ELossBeforeTPC;
}

const int kMaxPt      = 500;
class lariat::ELossBeforeTPC : public art::EDAnalyzer 
{
public:
  explicit ELossBeforeTPC(fhicl::ParameterSet const & p);
  virtual ~ELossBeforeTPC();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const & p);
  double  distance(double, double, double, double, double, double );

private:

  TTree* fTree;
  int    run;
  int    subrun;
  int    eventN;

  double truePVtxX ;
  double truePVtxY ;
  double truePVtxZ ;
  double trueVtxX ;
  double trueVtxY ;
  double trueVtxZ ;
  double trueEndX ;
  double trueEndY ;
  double trueEndZ ;
  double initialKE;
  double finalKE ;


  
  int nPointsBeforeTPC = 0;
  double TruePosX[kMaxPt];
  double TruePosY[kMaxPt];
  double TruePosZ[kMaxPt];
  double TrueKE[kMaxPt];


  std::vector<std::string> G4Process; 
  const static int thisFreakingCount = 9;
  TH2D *hXZ; 
  TH2D *hYZ; 
  TH1D *hEDep[thisFreakingCount]; 
  TH1D *hEDep_tot; 
  
  double minX =  0.0;
  double maxX = 47.0;
  double minY =-20.0;
  double maxY = 20.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 90.0;

  /*
  double interestingZBoundaries[thisFreakingCount-1] = { -58.0, -52.765, -39.4, -20.,  
							 -8.0, -5.575, -5.565,  -4.995,  
							 -1.265, 0.};
  
  */
  double interestingZBoundaries[thisFreakingCount-1] = { -58.0, 
							 -52.765, 
							 -7.8,
							 -7.2, 
							 -5.6, 
							 -5.05,
							 -1.28, 
							 0.1   };
							 

  int evtsPreTPC  = 0;
  int evtsInTheMiddle = 0;
  int evtsPostTPC = 0;
  int throughgoing = 0;
  int interactingInTPC = 0;
};


lariat::ELossBeforeTPC::ELossBeforeTPC(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::ELossBeforeTPC::~ELossBeforeTPC()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::ELossBeforeTPC::reconfigure(fhicl::ParameterSet const & pset)
{
}

void lariat::ELossBeforeTPC::analyze(art::Event const & evt)
{

  bool verbose = false;
  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
  
  // Get the backtracker to recover true quantities
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  double E_Boundary[thisFreakingCount];
  double Z_Boundary[thisFreakingCount];
  double eDep[thisFreakingCount-1];
  for (int i = 0; i < thisFreakingCount; i++) 
    {
      E_Boundary[i] = 0.;
      Z_Boundary[i] = 0.;
      if (i < thisFreakingCount-1) eDep[i] = 0.; 
    }
  

  //Set a bogus branch quantities, to be overwritten  
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  eventN  = evt.event();
  
  truePVtxX = -999.;
  truePVtxY = -999.;
  truePVtxZ = -999.;
  trueVtxX = -999.;
  trueVtxY = -999.;
  trueVtxZ = -999.;
  trueEndX = -999.;
  trueEndY = -999.;
  trueEndZ = -999.;
  trueEndZ = -999.;
  initialKE= -999.;
  finalKE  = -999.;
  nPointsBeforeTPC = 0;
  for (int i = 0; i < 500; i++) {TruePosX[i] = -999.; TruePosY[i] = -999.; TruePosZ[i] = -999.; TrueKE[i] = -999.;}   

  // ### Looping over all the Geant4 particles from the BackTracker ###
  for(size_t p = 0; p < plist.size(); ++p) 
    {
      // Get the true particle and its process, skip whatever is not primary 
      auto mcPart = plist.Particle(p);
      std::string proc = mcPart->Process();
      if ( !(proc.find("primary") != std::string::npos) ) continue;

      // Get the True Trajectory point
      simb::MCTrajectory truetraj = mcPart->Trajectory();
      // Make Sure we get the beamline primary                                                                                                                       
      if ( ( (truetraj.begin())->first).Z() >  -50. ) continue;

      // Get the mass (in GeV)
      double mass = mcPart->Mass() ;
      if (verbose)  std::cout<<"mass "<<mass<<"\n";

      
      //Store the kinetic energy and momentum on z at WC4. Just for cross check 
      // Store information at WC4
      auto inTPCPoint  = truetraj.begin(); 
      auto Position0   = inTPCPoint->first;
      trueVtxX = Position0.X();
      trueVtxY = Position0.Y();
      trueVtxZ = Position0.Z();

      auto Momentum0   = inTPCPoint->second;
      truePVtxX = Momentum0.X();
      truePVtxY = Momentum0.Y();
      truePVtxZ = Momentum0.Z();      
      initialKE = 1000*(TMath::Sqrt(Momentum0.X()*Momentum0.X() + Momentum0.Y()*Momentum0.Y() + Momentum0.Z()*Momentum0.Z() + mass*mass ) - mass); //I want this in MeV
      
      E_Boundary[0] = initialKE;
      Z_Boundary[0] = trueVtxZ;
      //--------------------------------------------------------
      // Identify the first trajectory point inside the TPC
      // Loop From First TrajPoint --> First Point in TPC 
      // Stop when you get into the TPC
      int iCounter  = 0;
      int iCounter2 = 0;
      //auto oldThreshold = truetraj.begin();
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  hXZ->Fill(pos.Z(),pos.X());
	  hYZ->Fill(pos.Z(),pos.Y());

	  auto momentumT = t->second;
	  double thisKET = 1000*(TMath::Sqrt(momentumT.X()*momentumT.X() + 
					     momentumT.Y()*momentumT.Y() + 
					     momentumT.Z()*momentumT.Z() + 
					     mass*mass ) - mass); //I want this in MeV  
	  
	  ////std::cout<<iCounter2<<" "<<pos.X()<<" , "<<pos.Y()<<" , "<<pos.Z()<<" , "<< thisKET <<" \n ";
	  TruePosX[iCounter2] = pos.X();
	  TruePosY[iCounter2] = pos.Y();
	  TruePosZ[iCounter2] = pos.Z();
	  TrueKE  [iCounter2] = thisKET;
	  nPointsBeforeTPC++;
	  iCounter2++;


	  
	    if (pos.Z() >= interestingZBoundaries[iCounter]) 
	    {
	    auto MomentumThis = t->second;
	    double thisKE = 1000*(TMath::Sqrt(MomentumThis.X()*MomentumThis.X() + 
	    MomentumThis.Y()*MomentumThis.Y() + 
	    MomentumThis.Z()*MomentumThis.Z() + 
						mass*mass ) - mass); //I want this in MeV  
						iCounter++;
						E_Boundary[iCounter] = thisKE;
						Z_Boundary[iCounter] = pos.Z();
						//std::cout<<" "<<pos.Z()<<" "<<interestingZBoundaries[iCounter-1]<<" "<< E_Boundary[iCounter] <<"\n";  
						}
						
	  
	  if (pos.Z() < minZ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    inTPCPoint = t;
	    break; // WARNING WARNING: we need to make sure this is the first traj point in TPC!!! I checked, it is.
	  }
	}// End search for first point in TPC
      //std::cout<<"\n";
      // if the first point is not in the TPC, we're not interested in this event. 
      // Fill some histos for cross check
      
      evtsPreTPC++;    
      if (inTPCPoint == truetraj.begin()) continue;
      evtsPostTPC++;     
      
      auto PositionN   = inTPCPoint->first;
      trueEndX = PositionN.X();
      trueEndY = PositionN.Y();
      trueEndZ = PositionN.Z();

      auto MomentumN   = inTPCPoint->second;      
      finalKE = 1000*(TMath::Sqrt(MomentumN.X()*MomentumN.X() + MomentumN.Y()*MomentumN.Y() + MomentumN.Z()*MomentumN.Z() + mass*mass ) - mass); //I want this in MeV
      //      std::cout<<"Final E "<<finalKE<<"\n";
      E_Boundary[thisFreakingCount-1] = finalKE;
      Z_Boundary[thisFreakingCount-1] = trueEndZ;
      for (int i = 0; i < thisFreakingCount; i++ )
	if (verbose) std::cout<<i<<" "<<E_Boundary[i]<<" "<<Z_Boundary[i]<<"\n";
      
      //      std::cout<<"\n";

      for (int i = 0; i < thisFreakingCount-1; i++ ) 
	{ 
	  eDep[i] = E_Boundary[i] - E_Boundary[i+1];      
	  hEDep[i]->Fill( eDep[i] );
	}
      hEDep_tot->Fill(E_Boundary[0] - E_Boundary[thisFreakingCount-1]);
      
    }


  
  //for (int i = 0; i < 500; i++) std::cout<<TruePosX[i] <<" "<<TruePosX[i] <<" "<< TruePosX[i] <<" "<< TrueKE[i] <<" \n ";
  fTree->Fill();
  //  std::cout<<"\n";
  
  for (int i = 0; i < thisFreakingCount; i++) 
    {
      E_Boundary[i] = 0.;
      if (i < thisFreakingCount-1) eDep[i] = 0.; 
    } 

  run     = -999; 
  subrun  = -999; 
  eventN  = -999;

  truePVtxX = -999.;
  truePVtxY = -999.;
  truePVtxZ = -999.;
  trueVtxX = -999.;
  trueVtxY = -999.;
  trueVtxZ = -999.;
  trueEndX = -999.;
  trueEndY = -999.;
  trueEndZ = -999.;
  trueEndZ = -999.;
  initialKE= -999.;
  finalKE  = -999.;
  G4Process.clear();
  nPointsBeforeTPC = 0;
  for (int i = 0; i < 500; i++) {TruePosX[i] = -999.; TruePosY[i] = -999.; TruePosZ[i] = -999.; TrueKE[i] = -999.;}   
}

void lariat::ELossBeforeTPC::endJob()
{
  std::cout<<"-------------------------------------------"<<std::endl;
  std::cout<<"True Events pre-TPC .............. "<<evtsPreTPC<<std::endl;
  std::cout<<"True Events pre-TPC .............. "<<evtsInTheMiddle<<std::endl;
  std::cout<<"True Events post-TPC ............. "<<evtsPostTPC<<std::endl;
  std::cout<<"True Throughgoing    ............. "<<throughgoing<<std::endl;
  std::cout<<"True interactingInTPC ............ "<<interactingInTPC<<std::endl;
  std::cout<<"-------------------------------------------"<<std::endl;

 
}

double lariat::ELossBeforeTPC::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return d;
}
void lariat::ELossBeforeTPC::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hXZ     = tfs->make<TH2D>("hXZ"    , "hXZ"    ,  11000, -100,  10,  20,  24,  26); 
  hYZ     = tfs->make<TH2D>("hYZ"    , "hYZ"    ,  11000, -100,  10,  20,  -1,   1); 
  
  for (int i = 0; i < thisFreakingCount-1; i++ ) hEDep[i] = tfs->make<TH1D>(Form("hEDep_%.2f",interestingZBoundaries[i]) , Form("hEDep_Z_%.2f",interestingZBoundaries[i]),  400,   0, 20); 
  hEDep_tot    = tfs->make<TH1D>("hEDep_tot" ,"hEDep_Tot",  400,   0, 100); 
  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&run      ,"run/I");
  fTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  fTree->Branch("eventN"   ,&eventN   ,"eventN/I");
 
  fTree->Branch("truePVtxX" ,&truePVtxX ,"truePVtxX/D");
  fTree->Branch("truePVtxY" ,&truePVtxY ,"truePVtxY/D");
  fTree->Branch("truePVtxZ" ,&truePVtxZ ,"truePVtxZ/D");
  fTree->Branch("trueVtxX" ,&trueVtxX ,"trueVtxX/D");
  fTree->Branch("trueVtxY" ,&trueVtxY ,"trueVtxY/D");
  fTree->Branch("trueVtxZ" ,&trueVtxZ ,"trueVtxZ/D");
  fTree->Branch("trueEndX" ,&trueEndX ,"trueEndX/D");
  fTree->Branch("trueEndY" ,&trueEndY ,"trueEndY/D");
  fTree->Branch("trueEndZ" ,&trueEndZ ,"trueEndZ/D");
  fTree->Branch("finalKE"  ,&finalKE  ,"finalKE/D" );
  fTree->Branch("initialKE",&initialKE,"initialKE/D" );

  fTree->Branch("nPointsBeforeTPC",&nPointsBeforeTPC,"nPointsBeforeTPC/I" );
  fTree->Branch("TruePosX",TruePosX,"TruePosX[nPointsBeforeTPC]/D" );
  fTree->Branch("TruePosY",TruePosY,"TruePosY[nPointsBeforeTPC]/D" );
  fTree->Branch("TruePosZ",TruePosZ,"TruePosZ[nPointsBeforeTPC]/D" );
  fTree->Branch("TrueKE"  ,TrueKE  ,"TrueKE[nPointsBeforeTPC]/D" );

  fTree->Branch("G4Process",&G4Process);
}



DEFINE_ART_MODULE(lariat::ELossBeforeTPC)
