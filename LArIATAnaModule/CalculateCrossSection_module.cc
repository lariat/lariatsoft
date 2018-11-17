////////////////////////////////////////////////////////////////////////
// Class:       CalculateCrossSection
// Module Type: analyzer
// File:        CalculateCrossSection_module.cc
//
// This code proceeds as follows
//
// 1) Get reconstructed data products
// 2) Fill corresponding vectors
// 3) Manipulate those vectors so that they make sense for cross section results 
// 4) Store in TTree and Cross Section plots
// 
// To DO:
// [ ] If the track is inverted, do something
// [ ] Test
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
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"



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

namespace lariat 
{
  class CalculateCrossSection;
}

class lariat::CalculateCrossSection : public art::EDAnalyzer 
{
public:
  explicit CalculateCrossSection(fhicl::ParameterSet const & p);
  virtual ~CalculateCrossSection();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void clearVar();
  void reconfigure(fhicl::ParameterSet const & p);
  double distance(double, double, double, double, double, double );
  double energyLossCalculation();
  double energyLossCalculation(double, double, bool);
  bool isWithinActiveVolume(double, double, double);
  std::vector<double> apply_permutation(const std::vector<double>& vec,
					const std::vector<std::size_t>& p);
private:
  std::string fTrackModuleLabel; 
  std::string fCalorimetryModuleLabel;
  std::string fWCTrackLabel; 
  std::string fWC2TPCModuleLabel; 
  bool fisData;

  // Fiducial Volume  
  double minX =  1.0; //0.
  double maxX = 46.0; //47
  double minY =-15.0; //20
  double maxY = 15.0; //20
  double minZ =  1.0; // G10 at z=1.2
  double maxZ = 86.0;
  // dEdX maximum value for Cross Section hadron
  double hitDEDXThreshold = 40.;
  
  //Some counting
  int nevts           = 0;
  int matchedTracks   = 0;
  int validCaloOnColl = 0;

  //Cross Section Plots
  TH1D *hRecoMCInteractingKE; 
  TH1D *hRecoMCIncidentKE; 
  TH1D *hRecoMCCrossSection;
     
  // Ttree for cross check variables  
  TTree* fTree;
  int    run;
  int    subrun;
  int    eventN;
  bool   isTrackInteracting; // If this track is interacting put this to 1, if not, put this to 0
  double WCMom;
  double WC4X;
  double WC4Y;
  double mass;

  double EnLoss;
  double recoVtxX  ;
  double recoVtxY  ;
  double recoVtxZ  ;
  double recoEndX  ;
  double recoEndY  ;
  double recoEndZ  ;
  double recoKEFF  ;
  double recoLength;
  int    recoPoints = 241;
  double recoInteractingKE;
  int    nTracksTPCFF;  
  // Branch in TTree with multiple entries per event
  double recoIncidentKE[241]; 
  double recoDEDX[241] ;      
  double recoPitch[241];      
  double recoEnDep[241];      
  double recoResRange[241];   
  double recoZPosition[241]; 
  
};


lariat::CalculateCrossSection::CalculateCrossSection(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::CalculateCrossSection::~CalculateCrossSection()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::CalculateCrossSection::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel             = pset.get< std::string >("TrackModuleLabel"      , "pmtracktc");
  fCalorimetryModuleLabel       = pset.get< std::string >("CalorimetryModuleLabel", "calo"     );
  fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel"          , "wctrack"  );
  fWC2TPCModuleLabel      	= pset.get< std::string >("WC2TPCModuleLabel"     , "WC2TPCtrk");
  fisData      	= pset.get<  bool  >("isData"    , true);
  mass      	= pset.get< double >("fmass"     , 139.57018);
}

void lariat::CalculateCrossSection::analyze(art::Event const & evt)
{ 
  // Let's get a clean start
  clearVar();
 
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  eventN  = evt.event();
 
 // Let's take some Reco stuff...
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;   // Container of wc tracks  
  std::vector<art::Ptr<ldp::WCTrack> >     wctrack;         // Vector of wc tracks  
  art::Handle< std::vector<recob::Track> > trackListHandle; // Container of reco tracks
  std::vector<art::Ptr<recob::Track> >     tracklist;       // Vector of wc tracks  

  // Find which recontructed track is associated with the WC
  int matchedRecoTrkKey = -99999;
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;
  if(!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return;
  nevts++;
  
  // Fill the container of reco tracks
  art::fill_ptr_vector(wctrack, wctrackHandle);
  art::fill_ptr_vector(tracklist, trackListHandle);

  art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel); 
  if (fWC2TPC.isValid())   
    {
      for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn ) 
	{
	  cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn)); 
	  if (!trackWC2TPC) continue;   
	  recob::Track const& aTrack(trackWC2TPC.ref()); 
	  matchedRecoTrkKey = aTrack.ID(); // This checks out OK   
	}
    }// if there's an associated track
  // Now we know what ID identifies the reco track we're interested in.
 
  if (matchedRecoTrkKey == -99999) return;
  matchedTracks++;
  // Let's define the containers of the important information about the track for the XS
  // Energy deposited and other important vectors (the bulk of the XS study)
  std::vector<double> recoPitch_v      ;
  std::vector<double> recoDEDX_v       ;
  std::vector<double> recoEDep_v       ; // Energy deposition at each slice. Simply dEdX*Pitch
  std::vector<double> recoResR_v       ;
  std::vector<double> recoZPos_v       ;
  std::vector<double> recoIncidentKE_v ;  
  
  //Let's get the calorimetry of the matched track
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  art::Ptr<recob::Track> theRecoTrk; // this is our girl! 
  
  bool invertedTracking = false;   // Is the tracking inverted? If it is, we need to flip all the vectors later... 
  int  countIDmatch     = 0;       // Let's check that no ID is repeated

  //Loop on the reco tracks
  for (auto recoTrk : tracklist ) 
    {
      //we keep only the matched one
      if (recoTrk->ID() != matchedRecoTrkKey) continue;
      theRecoTrk = recoTrk; // here she is!!! Good girl! 
      countIDmatch++;
      if (countIDmatch > 1 ) throw cet::exception("Geometry")    /// <== ... throw exception
			       << "More than one matched tracks... somethings going on, may I smell... \n";
      
      // Let's make sure the track first point is actually it (directionality problem).
      auto realFirstValidPt = (theRecoTrk->TrajectoryPoint(theRecoTrk->FirstValidPoint())).position;
      auto realLastValidPt  = (theRecoTrk->TrajectoryPoint(theRecoTrk->LastValidPoint( ))).position;      
      //      if ( TMath::Abs(realFirstValidPt.Z() - realLastValidPt.Z()) <  0.2) continue; // This track was unusable anyhow
      if ( realFirstValidPt.Z() - realLastValidPt.Z() > 0)
	{
	  invertedTracking = true;
	  std::cout<<"point "<<realFirstValidPt.Z()<<" "<<realLastValidPt.Z()<<" "<<invertedTracking<<"\n";
	  auto bogusFirst  = realFirstValidPt;
	  realFirstValidPt = realLastValidPt;
	  realLastValidPt  = bogusFirst;
	}
      
      recoVtxX   = realFirstValidPt.X();         recoEndX   = realLastValidPt.X();
      recoVtxY   = realFirstValidPt.Y();         recoEndY   = realLastValidPt.Y();
      recoVtxZ   = realFirstValidPt.Z();         recoEndZ   = realLastValidPt.Z();
      
            // PileUp Check
      std::cout<<"\n\n------------------>"<<run<<" "<<subrun<<" "<<eventN<<"\n"<<"Track " << recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<" --- " << recoEndX<<" "<<recoEndY<<" "<<recoEndZ<<"\n";

      for (auto recoTrkPileUp : tracklist ) 
	{
	  if (recoTrkPileUp->ID() == matchedRecoTrkKey) continue; // Skip self track

	  auto zFirst = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->FirstValidPoint())).position).Z();
	  auto zLast  = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->LastValidPoint( ))).position).Z();
	 
	  auto minZPileUp = zFirst;
	  auto minXPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->FirstValidPoint())).position).X();
	  auto minYPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->FirstValidPoint())).position).Y();

	  auto maxZPileUp = zLast;
	  auto maxXPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->LastValidPoint())).position).X();
	  auto maxYPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->LastValidPoint())).position).Y();
	  
	  if (zLast < zFirst) 
	    {
	      minZPileUp = zLast;
	      minXPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->LastValidPoint())).position).X();
	      minYPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->LastValidPoint())).position).Y();

	      maxZPileUp = zFirst;
	      maxXPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->FirstValidPoint())).position).X();
	      maxYPileUp = ((recoTrkPileUp->TrajectoryPoint(recoTrkPileUp->FirstValidPoint())).position).Y();

	    }
	  
	  
	  double tracksDistance = TMath::Sqrt((recoEndX-minXPileUp)*(recoEndX-minXPileUp) + 
	  				      (recoEndY-minYPileUp)*(recoEndY-minYPileUp) + 
	  				      (recoEndZ-minZPileUp)*(recoEndZ-minZPileUp) );
	  
	  
	  
	  //	  if (minZPileUp < 14. && tracksDistance > 2.) {
	  if (minZPileUp < 14.) {
	    std::cout<<"          vX: "<<minXPileUp<<" vY: "<<minYPileUp<<" vZ: "<<minZPileUp<< " --  eX: "<<maxXPileUp<<" eY: "<<maxYPileUp<<" eZ: "<<maxZPileUp<<" tracksDistance: "<<tracksDistance<<"\n";
	    nTracksTPCFF++;}
	  minZPileUp = 0;
	} // End PileUp Check
      
      std::cout<<"nTracksTPCFF: "<<nTracksTPCFF<<"\n";
  
      size_t furtherstInZCaloPointIndex  = 0;
      double furtherstInZCaloPointZ      = 0.;
      // if the calorimetry for my reco track is valid, I'll fill the energy deposition at each slice
      if (fmcal.isValid()) {
	std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(theRecoTrk.key());
	
	for (size_t j = 0; j<calos.size(); ++j){ 
	  if (!calos[j]->PlaneID().isValid)   continue;
	  if  (calos[j]->PlaneID().Plane == 0) continue;  // Skip Induction Plane
	  validCaloOnColl++;
	  //double recoEnDeposited = 0.;
	  for (size_t k = 0; k<calos[j]->dEdx().size(); ++k)         
	    {
	      if (calos[j]->dEdx()[k] > hitDEDXThreshold ) continue;
	      //If the following happens, there's a mess in the calorimetry module, skip
	      if (calos[j]->XYZ()[k].Z() < 0 || calos[j]->XYZ()[k].Z() > 90. ) continue;
	      // Let's register the furtherst point in Z to determine if the track interacted or not. 
	      if (calos[j]->XYZ()[k].Z() > furtherstInZCaloPointZ) furtherstInZCaloPointIndex = k;
	      bool isThisPointInActiveVolume = isWithinActiveVolume(calos[j]->XYZ()[k].X(), calos[j]->XYZ()[k].Y(), calos[j]->XYZ()[k].Z());
	      if (!isThisPointInActiveVolume) continue; // I'm filling the track only with points in active volume
	      recoPitch_v.push_back(calos[j]->TrkPitchVec()[k]);
	      recoDEDX_v .push_back(calos[j]->dEdx()[k]);
	      recoEDep_v .push_back(calos[j]->dEdx()[k] * calos[j]->TrkPitchVec()[k]);
	      recoResR_v .push_back(calos[j]->ResidualRange()[k]);
	      recoZPos_v .push_back(calos[j]->XYZ()[k].Z());
	    } // Loop on calo points
	  
	  if (!furtherstInZCaloPointIndex) 
	    {//throw cet::exception("Geometry")    /// <== ... throw exception
					   //  << "Track doesn't go in? ... somethings going on, may I smell... \n";
	      recoPitch_v.clear();
	      recoPitch_v.clear();
	      recoDEDX_v .clear();
	      recoEDep_v .clear();
	      recoResR_v .clear();
	      recoZPos_v .clear();
	    }
	  // Let's check if this track is interacting
	  if (furtherstInZCaloPointIndex) isTrackInteracting =  isWithinActiveVolume(calos[j]->XYZ()[furtherstInZCaloPointIndex].X(), 
										     calos[j]->XYZ()[furtherstInZCaloPointIndex].Y(), 
										     calos[j]->XYZ()[furtherstInZCaloPointIndex].Z());
	  
	} // Loop on Collection Vs Induction
      }// Is calorimetry valid?
      
    }// Loop on tracks 



  WCMom = wctrack[0]->Momentum();
  WC4X  = wctrack[0]->HitPosition(3,0);
  WC4Y  = wctrack[0]->HitPosition(3,1);
  


  if (!recoDEDX_v.size()) return; // If there are no points, return
  if (countIDmatch != 1) return; // If for some reason I couldn't find the match track, I should return
  if (recoEndX < -100. || recoEndY < -100. || recoEndZ < -100.) throw cet::exception("Geometry")    /// <== ... throw exception
								  << "Very Weird End of track... somethings going on, may I smell... \n";
  if (WCMom < 0.) throw cet::exception("Geometry")    /// <== ... throw exception
		    << "WCMom < 0... somethings going on, may I smell... \n";
  
  if (mass < 0.) throw cet::exception("Geometry")    /// <== ... throw exception
		   << "Mass < 0... somethings going on, may I smell... \n";
  
  if (recoPitch_v.size() != recoDEDX_v.size()) throw cet::exception("Geometry")    /// <== ... throw exception
						 << "Different dimensions ... somethings going on, may I smell... \n";
  if (recoEDep_v.size() != recoResR_v.size()) throw cet::exception("Geometry")    /// <== ... throw exception
						<< "Different dimensions ... somethings going on, may I smell... \n";
  if (recoPitch_v.size() != recoZPos_v.size()) throw cet::exception("Geometry")    /// <== ... throw exception
						 << "Different dimensions ... somethings going on, may I smell... \n";
  if (recoPitch_v.size() != recoEDep_v.size()) throw cet::exception("Geometry")    /// <== ... throw exception
						 << "Different dimensions ... somethings going on, may I smell... \n";

  //If the z are inverted, flip the vectors
 if (invertedTracking)
    {
      //Is track inverted
      std::reverse(recoPitch_v.begin(),recoPitch_v.end());
      std::reverse(recoDEDX_v .begin(),recoDEDX_v .end());
      std::reverse(recoEDep_v .begin(),recoEDep_v .end());
      std::reverse(recoResR_v .begin(),recoResR_v .end());
      std::reverse(recoZPos_v .begin(),recoZPos_v .end());

    }

  // Let's check if this track is interacting
 //  isTrackInteracting = isWithinActiveVolume(recoEndX, recoEndY, recoEndZ); 
  //Kinetic energy at the Wire Chamber 4 for this mass hypothesis
  double WCKE             = TMath::Sqrt(WCMom*WCMom + mass*mass) - mass;  
  //How much energy is lost before getting to the TPC?
  double calculatedEnLoss = energyLossCalculation(); 
  if (fisData)
    {
      double WCTheta  = wctrack[0]->Theta();
      double WCPhi    = wctrack[0]->Phi();
      double tanThetaCosPhi = TMath::Tan(WCTheta) * TMath::Cos(WCPhi);
      double tanThetaSinPhi = TMath::Tan(WCTheta) * TMath::Sin(WCPhi);
      double den = TMath::Sqrt(1+tanThetaCosPhi*tanThetaCosPhi);
      double onTheFlyPz = WCMom/den;
      double onTheFlyPx = onTheFlyPz*tanThetaSinPhi;
      calculatedEnLoss = energyLossCalculation(WC4X,onTheFlyPx,fisData);
    }else
    {

      // Get the backtracker to recover true quantities
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      const sim::ParticleList& plist = pi_serv->ParticleList();

      for(size_t p = 0; p < plist.size(); ++p)
	{
	  auto mcPart = plist.Particle(p);
	  std::string proc = mcPart->Process();
	  if ( !(proc.find("primary") != std::string::npos) ) continue;
	  simb::MCTrajectory truetraj = mcPart->Trajectory();
	  auto inTPCPoint  = truetraj.begin();
	  auto Momentum0   = inTPCPoint->second;
	  double onTheFlyPx = 1000*Momentum0.X(); // the Momentum Needs To be in MeV
	  calculatedEnLoss = energyLossCalculation(WC4X,onTheFlyPx,fisData);
	}
    }

  //Initial energy at the TPC Front Face. 
  const double initialKE  = WCKE - calculatedEnLoss;
  recoIncidentKE_v.push_back(initialKE);

  // At this point, I assume the energy deposited vector are filled in the right direction
  double energyDepositedThusFar = 0.;
  for (size_t iDep = 0; iDep < recoEDep_v.size(); iDep++)
    {
      energyDepositedThusFar += recoEDep_v[iDep];
      recoIncidentKE_v.push_back(initialKE-energyDepositedThusFar);
    }

  //Fill the Incident and Interacting plots real quick
  for (auto inKE : recoIncidentKE_v) hRecoMCIncidentKE->Fill(inKE); 
  if (isTrackInteracting && recoIncidentKE_v.size())  hRecoMCInteractingKE  ->Fill(recoIncidentKE_v[recoIncidentKE_v.size()-1]); 



  EnLoss    = calculatedEnLoss;
  recoKEFF  = initialKE;
  recoInteractingKE = recoIncidentKE_v[recoIncidentKE_v.size()-1];
  for (size_t i = 0; i < recoIncidentKE_v.size(); i++) recoIncidentKE[i] = recoIncidentKE_v[i];
  double sumOfPitchs = 0;
  for (size_t i = 0; i < recoDEDX_v.size(); i++)
    {
      recoDEDX     [i] = recoDEDX_v [i];   
      recoPitch    [i] = recoPitch_v[i];   
      recoEnDep    [i] = recoEDep_v [i];   
      recoResRange [i] = recoResR_v [i]; 
      recoZPosition[i] = recoZPos_v [i];
      sumOfPitchs += recoPitch_v[i];
     }

  recoLength = sumOfPitchs;



  fTree->Fill();
  


  // Let's clean up after ourselves
  clearVar();  
  recoDEDX_v .clear();
  recoPitch_v.clear();
  recoEDep_v .clear();
  recoResR_v .clear();
  recoZPos_v .clear();
  recoIncidentKE_v.clear();
}

void lariat::CalculateCrossSection::clearVar()
{
  run       = -999 ; subrun     = -999 ;  eventN = -999;
  WCMom     = -999.;  
  WC4X      = -999.;  WC4Y      = -999.;
  recoVtxX  = -999.;  recoEndX  = -999.;
  recoVtxY  = -999.;  recoEndY  = -999.;
  recoVtxZ  = -999.;  recoEndZ  = -999.;
  recoKEFF  = -999.;  recoLength= -999.;
  recoInteractingKE   = -999.;
  nTracksTPCFF        = 1;
  isTrackInteracting  = 0; // If this track is interacting put this to 1, if not, put this to 0
  for (int i = 0; i < 240; i++){
    recoIncidentKE[i] = -999.;    recoDEDX[i]     = -999.;      recoPitch[i]     = -999.;      
    recoEnDep[i]      = -999.;    recoResRange[i] = -999.;      recoZPosition[i] = -999.;  
  }
}

void lariat::CalculateCrossSection::endJob()
{
  std::cout<<"\n\n################## Some checks ################# \n";
  std::cout<<"### All Events Up to here .............. "<< nevts            <<" ### \n";
  std::cout<<"### Match Tracks ....................... "<< matchedTracks   <<"  ### \n";
  std::cout<<"### Calo on Coll ....................... "<< validCaloOnColl <<"  ### \n";
  // ----------------------------------------------------------------
  // Create the cross section from the incident and interaction plots
  // ----------------------------------------------------------------
  //We shouldn't let a CPU handle large numbers if we don't have to
  // So, let's do the calculation on our own and just give the CPU what it needs
  // here is how you calculate the number density
  /*
    grams per kilo conversion = 1000  
    argon density rho         = 1396       kg/m^3
    molar mass                = 39.95      g/mol
    avogadro number           = 6.022e+23  number/mol
    number_density = rho*g_per_kg*avogadro/molar_mass = 0.21043 e+28 m^-3
  */
  // Then we use the number calculated for us, keeping the exponent on the side
  float number_density     = 2.1043083854818523; // * 10^28 m^-3
  float slab_width         = 0.45;//in cm = 4.5e-3 m
  float xsConversionFactor = 1./(number_density*slab_width); // *10^-26 m^2 ==> * 100 barn  
  
  
  // ===============================================================================================================
  //					MAKING THE CROSS-SECTION PLOT
  // ===============================================================================================================
  
  // ###################################################################
  // #### Looping over the exiting bins to extract the cross-section ###
  // ###################################################################
  for( int iBin = 1; iBin <= hRecoMCInteractingKE->GetNbinsX(); ++iBin )
    {
      // ### If an incident bin is equal to zero then skip that bin ###
      if( hRecoMCIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
      // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
      float conversion_to_barns = 100.;
      float crossSection = (hRecoMCInteractingKE->GetBinContent(iBin)/hRecoMCIncidentKE->GetBinContent(iBin)) * xsConversionFactor*conversion_to_barns;
      
      
      // ### Putting the value on the plot
      hRecoMCCrossSection->SetBinContent(iBin,crossSection);
      
      // ###########################################################
      // ### Calculating the error on the numerator of the ratio ###
      // ###########################################################
      float numError = pow(hRecoMCInteractingKE->GetBinContent(iBin),0.5);
      float num = hRecoMCInteractingKE->GetBinContent(iBin);
      
      // ### Putting in a protection against dividing by zero ###   
      if(num == 0){continue;}
      float term1 = numError/num;
      
      // #################################################
      // ### Calculating the error on the denomentator ###
      // #################################################
      float denomError = pow(hRecoMCIncidentKE->GetBinContent(iBin),0.5);
      float denom = hRecoMCIncidentKE->GetBinContent(iBin);
      if(denom == 0){continue;}
      
      // ### Putting in a protection against dividing by zero ###
      float term2 = denomError/denom;
      
      float totalError = crossSection * pow( ( (term1*term1) + (term2*term2) ),0.5) * xsConversionFactor;
      hRecoMCCrossSection->SetBinError(iBin,totalError);
      
    }//<---End iBin Loop

}

double lariat::CalculateCrossSection::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return d;
}

double lariat::CalculateCrossSection::energyLossCalculation()
{
  return 40.;
}

double lariat::CalculateCrossSection::energyLossCalculation(double x, double px, bool isData)
{
  // x in cm and px in MeV
  double discriminant = 0.0733*px + 1.3* x - 31; 
  if ( discriminant >  0 ) 
    {
      //particles going through the halo hole
      if (isData) return 17.5;
      else return 24.5;
    }
  else
    {
      //particles going through the halo paddle
      if (isData) return 25.5;
      else return 32.5;
    }
}

bool lariat::CalculateCrossSection::isWithinActiveVolume(double x, double y, double z)
{
  if (x < minX ) return false; 
  if (x > maxX ) return false;
  if (y < minY ) return false; 
  if (y > maxY ) return false;
  if (z < minZ ) return false; 
  if (z > maxZ ) return false;
  return true;
}

void lariat::CalculateCrossSection::beginJob()
{
  
  std::vector<double> vectorA;
  std::vector<double> vectorB;
  vectorA.push_back(1.);
  vectorA.push_back(3.);
  vectorA.push_back(2.);
  vectorB.push_back(14.);
  vectorB.push_back(1.);
  vectorB.push_back(2.);

  for (size_t is = 0; is < vectorA.size(); is++) std::cout<<vectorA[is]<<" "<<vectorB[is]<<"\n";

  //  auto p = sort_permutation(vectorA,
  //			    [](double const& a, double const& b){ a < b; });

  //vectorA = apply_permutation(vectorA, p);
  //vectorB = apply_permutation(vectorB, p);

  for (size_t is = 0; is < vectorA.size(); is++) std::cout<<vectorA[is]<<" "<<vectorB[is]<<"\n";
  
  art::ServiceHandle<art::TFileService> tfs;
  //----------> RecoMC Cross Section
  hRecoMCInteractingKE       = tfs->make<TH1D>("hRecoInteractingKE"      , "Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hRecoMCIncidentKE        = tfs->make<TH1D>("hRecoIncidentKE"         , "Incident Kinetic Energy [MeV]"   , 42, -100, 2000);
  hRecoMCCrossSection      = tfs->make<TH1D>("hRecoMCCrossSection"     , "Cross-Section [barn]"            , 42, -100, 2000);

  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&run      ,"run/I");
  fTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  fTree->Branch("eventN"   ,&eventN   ,"eventN/I");
  fTree->Branch("isTrackInteracting"  ,&isTrackInteracting   ,"isTrackInteracting/O");
  fTree->Branch("recoPoints"          ,&recoPoints           ,"recoPoints/I");
  fTree->Branch("nTracksTPCFF"       ,&nTracksTPCFF         ,"nTracksTPCFF/I");

  fTree->Branch("recoVtxX" ,&recoVtxX ,"recoVtxX/D");    fTree->Branch("recoEndX" ,&recoEndX ,"recoEndX/D");
  fTree->Branch("recoVtxY" ,&recoVtxY ,"recoVtxY/D");    fTree->Branch("recoEndY" ,&recoEndY ,"recoEndY/D");
  fTree->Branch("recoVtxZ" ,&recoVtxZ ,"recoVtxZ/D");    fTree->Branch("recoEndZ" ,&recoEndZ ,"recoEndZ/D");
  fTree->Branch("recoKEFF" ,&recoKEFF ,"recoKEFF/D");

  fTree->Branch("WCMom",&WCMom,"WCMom/D");
  fTree->Branch("WC4X" ,&WC4X ,"WC4X/D");
  fTree->Branch("WC4Y" ,&WC4Y ,"WC4Y/D");
  fTree->Branch("mass" ,&mass ,"mass/D");
  fTree->Branch("recoInteractingKE" ,&recoInteractingKE ,"recoInteractingKE/D");

  fTree->Branch("recoDEDX"      , recoDEDX      ,"recoDEDX[recoPoints]/D");
  fTree->Branch("recoPitch"     , recoPitch	,"recoPitch[recoPoints]/D");
  fTree->Branch("recoEnDep"     , recoEnDep	,"recoEnDep[recoPoints]/D");
  fTree->Branch("recoResRange"  , recoResRange	,"recoResRange[recoPoints]/D");
  fTree->Branch("recoZPosition" , recoZPosition	,"recoZPosition[recoPoints]/D");
  fTree->Branch("recoIncidentKE", recoIncidentKE,"recoIncidentKE[recoPoints]/D");

}


std::vector<double> lariat::CalculateCrossSection::apply_permutation(const std::vector<double>& vec,
								     const std::vector<std::size_t>& p)
{
  std::vector<double> sorted_vec(vec.size());
  std::transform(p.begin(), p.end(), sorted_vec.begin(), [&](std::size_t i){ return vec[i]; });
  return sorted_vec;
}



DEFINE_ART_MODULE(lariat::CalculateCrossSection)
