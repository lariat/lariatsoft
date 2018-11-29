////////////////////////////////////////////////////////////////////////
// Class:       IntTrackingEff
// Module Type: analyzer
// File:        IntTrackingEff_module.cc
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
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include <map>
#include <iostream>
#include <fstream>

namespace lariat 
{
  class TrackingEffWCMatch;
}

class lariat::TrackingEffWCMatch : public art::EDAnalyzer 
{
public:
  explicit TrackingEffWCMatch(fhicl::ParameterSet const & p);
  virtual ~TrackingEffWCMatch();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const & p);

private:

  std::string fTrackModuleLabel ;
  std::string fCalorimetryModuleLabel     ;
  std::string fWCTrackLabel     ;
  std::string fWC2TPCModuleLabel;

  double minX =  0.0;
  double maxX = 47.0;
  double minY =-20.0;
  double maxY = 20.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 91.0;

  int evtsPreTPC  = 0;
  int evtsInTheMiddle = 0;
  int evtsPostTPC = 0;

  // TTree section
  TTree* fTree;
  // Generic variables
  int    run;
  int    subrun;
  int    eventN;   
  // Truth Variables
  int    isPrimaryInTPC;
  std::vector<std::string> G4Process;
  double intAngle;
  double trueInteractionEnergy;  
  double initialKE;
  double trueVtxX ;
  double trueVtxY ;
  double trueVtxZ ;
  double trueEndX ;
  double trueEndY ;
  double trueEndZ ;
  double trueL;  
  double trueEnDeposited;  
  double trueKEFrontFace;
  // Reco Variables
  int    nReco;
  int    isTrackMatched;
  double recoVtxX ;
  double recoVtxY ;
  double recoVtxZ ;
  double recoEndX ;
  double recoEndY ;
  double recoEndZ ;
  double recoL; 
  double recoEnDeposited;  
  double recoInteractionEnergy;  
  double keWhen;
  double vtxDist;
  double endDist;
  int truePDGMatch = -99999;
};


lariat::TrackingEffWCMatch::TrackingEffWCMatch(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::TrackingEffWCMatch::~TrackingEffWCMatch()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::TrackingEffWCMatch::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel             = pset.get< std::string >("TrackModuleLabel"      , "pmtracktc");
  fCalorimetryModuleLabel       = pset.get< std::string >("CalorimetryModuleLabel", "calo"     );
  fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel"          , "wctrack"  );
  fWC2TPCModuleLabel      	= pset.get< std::string >("WC2TPCModuleLabel"     , "WC2TPCtrk");
}

void lariat::TrackingEffWCMatch::analyze(art::Event const & evt)
{

  bool verbose = false;

  if (verbose) std::cout<<"\n -------------------- new event\n";
  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
 
  // Get the backtracker to recover true quantities
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // ######################################
  // ### Making a vector of MCParticles ###
  // ######################################   
 
  //Set bogus branch quantities, to be overwritten
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  eventN  = evt.event(); 
  // Truth Variables
  isPrimaryInTPC        = 0;
  intAngle              = -999.;
  trueInteractionEnergy = -999.;  
  initialKE             = -999.;
  trueVtxX              = -999.;
  trueVtxY              = -999.;
  trueVtxZ              = -999.;
  trueEndX              = -999.;
  trueEndY              = -999.;
  trueEndZ              = -999.;
  trueL                 = -999.;  
  trueEnDeposited       = -999.;  
  //Reco Variables
  nReco                 = 0;
  isTrackMatched        = 0;
  recoVtxX              = -999.;
  recoVtxY              = -999.;
  recoVtxZ              = -999.;
  recoEndX              = -999.;
  recoEndY              = -999.;
  recoEndZ              = -999.;
  recoL                 = -999.; 
  recoEnDeposited       = -999.;  
  recoInteractionEnergy = -999.;  
  vtxDist = -999.;  
  endDist = -999.;  
  truePDGMatch = -99999;

  //----------------------------------------------> Truth and Only The Truth <--------------------------------------

  bool keepInteraction = false;
  // ### Looping over all the Geant4 particles from the BackTracker ###
  for(size_t p = 0; p < plist.size(); ++p) 
    {

      auto mcPart = plist.Particle(p);
      std::string proc = mcPart->Process();
      //Skip whatever is not primary
      if ( !(proc.find("primary") != std::string::npos) ) continue;
      // Get the True Trajectory point
      simb::MCTrajectory truetraj = mcPart->Trajectory();
      // Make Sure we get the beamline primary                                                                                                                       
      if ( ( (truetraj.begin())->first).Z() >  -50. ) continue;


      std::string interactionLabel = "";

      //Get simIDE associated to the primary
      geo::View_t view = geom->View(0); 
      auto simIDE_Prim  = bt->TrackIdToSimIDEs_Ps(mcPart->TrackId(),view);
      // Order them in order of increasing Z
      std::map<double, sim::IDE> orderedSimIDE;
      for (auto ide : simIDE_Prim ) orderedSimIDE[ide->z] = *ide;
      
      //--------------------------------------------------------
      // Identify the first trajectory point in TPC
      auto inTPCPoint  = truetraj.begin();
      auto Momentum0   = inTPCPoint->second; 
      double mass = mcPart->Mass() ;
      initialKE = 1000*(TMath::Sqrt(Momentum0.X()*Momentum0.X() + Momentum0.Y()*Momentum0.Y() + Momentum0.Z()*Momentum0.Z() + mass*mass ) - mass); //I want this in MeV

      // Loop From First TrajPoint --> First Point in TPC 
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  if (pos.Z() < minZ || pos.Z() > maxZ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    inTPCPoint = t;
	    break;
	  }
	}// End search for first point in TPC
      
      // if the first point is not in the TPC, we'll register some bogus info and move on.
      evtsPreTPC++;
      interactionLabel = "notArrived";
      if (inTPCPoint == truetraj.begin()) 
	{ 
	  auto vtx = (truetraj.begin())->first;
	  auto end = std::prev(truetraj.end())->first;
	  isPrimaryInTPC        = 0;
	  intAngle              = -999.;
	  trueInteractionEnergy = -999.;  
	  trueVtxX              = vtx.X();
	  trueVtxY              = vtx.Y();
	  trueVtxZ              = vtx.Z();
	  trueEndX              = end.X();;
	  trueEndY              = end.Y();;
	  trueEndZ              = end.Z();;
	  trueL                 = -999.;  
	  trueEnDeposited       = -999.;
	  truePDGMatch          = -99999; 
	  G4Process.push_back("notArrived"); 
	  break;
	}
      // Now we know the primary is in the TPC
      isPrimaryInTPC = 1;
      evtsPostTPC++;     
      
      // Identify the last interesting trajectory point in TPC
      auto finTPCPoint = std::prev(truetraj.end()); 
      // The last point is a bit more complicated:
      // if there's no interaction, then it is simply the last point in the TPC
      // if there's one or more interaction points, it's the first interaction point deemed interesting (no coulomb)
      // Take the interaction Map... check if there's something there
      auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
      if (thisTracjectoryProcessMap.size())
	{
	  for(auto const& couple: thisTracjectoryProcessMap) 
	    { 
	      // I'm not interested in the CoulombScattering, discart this case
	      if (verbose) std::cout<<(truetraj.KeyToProcess(couple.second))<<"\n";
	      if ((truetraj.KeyToProcess(couple.second)).find("CoulombScat")!= std::string::npos) continue;
	      
	      // Let's check if the interaction is in the the TPC
	      auto     interactionPos4D =  (truetraj.at(couple.first)).first ;	
	      if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
	      else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
	      else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
	      G4Process.push_back(truetraj.KeyToProcess(couple.second));
	      // If we made it here, then this is the first interesting interaction in the TPC
	      // Our job is done!!! Great!

	      finTPCPoint = truetraj.begin() + couple.first; 
	      keepInteraction = true;
	      
	      //--------------------- Int Angle ---------------------------
	      // Try to retreive the interaction angle
	      auto  prevInteractionPos4D = (truetraj.at(couple.first-1)).first ;
	      auto  prevInteractionPos3D = prevInteractionPos4D.Vect() ;
	      auto  interactionPos3D     = interactionPos4D.Vect() ;
	      auto  distanceBtwPoint     = interactionPos3D - prevInteractionPos3D;
	      
	      //Let's try to see if the next point exists
	      if (truetraj.size() > couple.first + 1) 
		{
		  // The particle doesn't die. No need to check for anything else.
		  auto nextInteractionPos4D =  (truetraj.at(couple.first+1)).first ;
		  auto nextInteractionPos3D =  nextInteractionPos4D.Vect() ;
		  auto distanceBtwPointNext =  nextInteractionPos3D - interactionPos3D;
		  intAngle = TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  );
		}else
		{ // The particle has come to an end. Let's check the daugthers.
		  if (mcPart->NumberDaughters() == 0 ) break;
		  double minAngle = 99999.;
		  for(size_t d = 0; d < plist.size(); ++d) 
		    {
		      auto mcDaught = plist.Particle(d);
		      if (mcDaught->Mother()  != 1 ) continue;
		      if (mcDaught->NumberTrajectoryPoints () < 2 ) continue ;
		      if (!(TMath::Abs(mcDaught->PdgCode())  == 13  ||
			    TMath::Abs(mcDaught->PdgCode())  == 11  ||
			    TMath::Abs(mcDaught->PdgCode())  == 211 ||
			    TMath::Abs(mcDaught->PdgCode())  == 321 || 
			    TMath::Abs(mcDaught->PdgCode())  == 2212) ) continue;	
		      auto daugthTraj = mcDaught->Trajectory();
		      if (daugthTraj.TotalLength() < 0.5 ) continue ;
		      
		      
		      auto daughtFirstPt  =  ((daugthTraj[0]).first).Vect() ;
		      auto daughtSecondPt =  ((daugthTraj[1]).first).Vect() ;
		      auto distanceBtwPointNext = daughtSecondPt - daughtFirstPt;
		      intAngle = TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  );
		      if ( minAngle > intAngle ) minAngle = intAngle;
		    }
		  intAngle = minAngle;
		}
	      
	      intAngle *= 57.2958;
	      break;
	    }// Loop on interaction points
	}// If there are G4 interactions
      

      // If I didn't find anything interesting in the intereaction map, let's loop back!
      if ( !keepInteraction )
	{ 
	  //Loop through the daugthers and look for decays/captures
	  bool isFirstDaugther = true;
	  for(size_t d = 0; d < plist.size(); ++d) 
	    {
	      auto mcDaught = plist.Particle(d);
	      if (mcDaught->Mother()  != 1 ) continue;
	      if ((mcDaught->Process()).find("astic")!= std::string::npos) continue;
	      simb::MCTrajectory trueDaugthTraj = mcDaught->Trajectory();
	      if (verbose) std::cout<<"Daughters "<<mcDaught->Process()<<" "
				    << (trueDaugthTraj.begin()->first).X()<<" , "<<(trueDaugthTraj.begin()->first).Y() <<" , "<<(trueDaugthTraj.begin()->first).Z() <<std::endl;	      
	     

	      if (trueDaugthTraj.begin()->first.Z() < minZ || trueDaugthTraj.begin()->first.Z() > maxZ) continue;
	      else if (trueDaugthTraj.begin()->first.X() <   minX || trueDaugthTraj.begin()->first.X() > maxX ) continue;
	      else if (trueDaugthTraj.begin()->first.Y() <   minY || trueDaugthTraj.begin()->first.Y() > maxY ) continue;
	      else {
		if (isFirstDaugther){G4Process.push_back(mcDaught->Process()); isFirstDaugther = false;}
		break;
	      }
	    }
	  for ( auto t = std::prev(truetraj.end()); t!= truetraj.begin(); t--)
	    {
	      auto pos = t->first;
	      
	      if (pos.Z() > maxZ) continue;
	      else if (pos.X() <   minX || pos.X() > maxX ) continue;
	      else if (pos.Y() <   minY || pos.Y() > maxY ) continue;
	      else {
		finTPCPoint = t;
		break;
	      }
	    }
	}
      
      if (!G4Process.size())
        {
          G4Process.push_back("throughgoing");
	}
      //--------------------------------------------------------
      // If the first and last point in the TPC are the same, 
      // it's impossible to calculate a track lenght inside the TPC
      // (less than 2 trjPoints in the TPC)
      // Discart this case
     
      if ( !std::distance(finTPCPoint,inTPCPoint) ) continue;

      //--------------------------------------------------------
      // LET'S RECAP: 
      // If we made it here, we identified the first and the last points inside the TPC
      // AND they are different. Awesome!
      
      // Calculate true lenght in TPC (loop only on points in TPC)
      double trueTPCLength = 0.;	  
      for ( auto iter = inTPCPoint; iter!= std::next(finTPCPoint); iter++)
	{  
	  // take the 4-vector position  for this point
	  auto posThis4D = iter ->first;  
	  // take the 3-vector position  for this point
	  auto posThis = posThis4D.Vect();	
	  //Let's take the previous point so that we can calculate distances
	  auto iter_pre  = iter; 
	  iter_pre--;
	  // take the 4-vector and 3-vector position for the previous point	      
	  auto posPrev4D = iter_pre  ->first;  
	  auto posPrev = posPrev4D.Vect();
	  
	  // Make the difference
	  auto distanceBtwPoint = posThis - posPrev;
	  // Calculate the distance
	  trueTPCLength += distanceBtwPoint.Mag();
	}

      trueVtxX = (inTPCPoint->first ).X();
      trueVtxY = (inTPCPoint->first ).Y();
      trueVtxZ = (inTPCPoint->first ).Z();
      trueEndX = (finTPCPoint->first).X();
      trueEndY = (finTPCPoint->first).Y();
      trueEndZ = (finTPCPoint->first).Z();
      trueL    = trueTPCLength;  


      // -------------------> The true energy part
      if (trueL < 0.47 ) continue;

      // Ok, now wer have the first and last point in the TPC and they are reasonably distant
      // Let's use them!
      
      // We want to chop up the points between the fist and list uniformely
      // and ordered by Z
      // Order them in order of increasing Z
      std::map<double, TVector3> orderedUniformTrjPts;
      // We want the first and last uniform point to coincide with the 
      // the first and last points we just found 
      auto positionVector0 = (inTPCPoint ->first).Vect(); 
      auto positionVector1 = (finTPCPoint->first).Vect(); 
      orderedUniformTrjPts[positionVector0.Z()] = positionVector0;
      orderedUniformTrjPts[positionVector1.Z()] = positionVector1;

      const double trackPitch = 0.47;
      // I do have space for at least one extra point, so let's put it there!
      // Calculate how many extra points I need to put between the new first point and the second TrajPoint
      int    nPts            = (int) (trueL/trackPitch);
      for (int iPt = 1; iPt <= nPts; iPt++ )
	{
	  auto newPoint = positionVector0 + iPt*(trackPitch/trueL) * (positionVector1 - positionVector0);
	  orderedUniformTrjPts[newPoint.Z()] = newPoint;
	}


      // If the distance between the last point and the second to last is less then 0.235
      // eliminate the second to last point
      
      auto lastPt         = (orderedUniformTrjPts.rbegin())->second;
      auto secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
      double lastDist = TMath::Sqrt((lastPt.X() - secondtoLastPt.X())*(lastPt.X() - secondtoLastPt.X()) +
				    (lastPt.Y() - secondtoLastPt.Y())*(lastPt.Y() - secondtoLastPt.Y()) +
				    (lastPt.Z() - secondtoLastPt.Z())*(lastPt.Z() - secondtoLastPt.Z())) ;
      
      
      if (lastDist < 0.235)
	{
	  orderedUniformTrjPts.erase((std::next(orderedUniformTrjPts.rbegin()))->first );
	} 

      
      // Calculate the initial kinetic energy
      auto initialMom =     inTPCPoint->second;
      trueKEFrontFace = 1000*(TMath::Sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y() + initialMom.Z()*initialMom.Z() + mass*mass ) - mass); 
      double kineticEnergy = trueKEFrontFace;

      auto old_it = orderedUniformTrjPts.begin();
      for (auto it = std::next(orderedUniformTrjPts.begin()); it != orderedUniformTrjPts.end(); it++, old_it++ )
	{
	  auto oldPos        = old_it->second;
	  auto currentPos    =     it->second;
	  
	  double uniformDist =  (currentPos - oldPos).Mag();
	  
	  //Calculate the energy deposited in this slice	  
	  auto old_iter = orderedSimIDE.begin();
	  double currentDepEnergy = 0.;
	  for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++)
	    {
	      auto currentIde = iter->second;
	      if ( currentIde.z < oldPos.Z()) continue;
	      if ( currentIde.z > currentPos.Z()) continue;
	      currentDepEnergy += currentIde.energy;
	    }// Determing which simIDE is within the current slice

	  // avoid overfilling super tiny energy depositions
	  if (currentDepEnergy/uniformDist < 0.1 ) continue;
	  //Calculate the current kinetic energy
	  kineticEnergy -= currentDepEnergy;
	  trueEnDeposited  = currentDepEnergy;

	}// Loop on OrderedPoints


      trueInteractionEnergy = kineticEnergy; 

      

    } // end geant particle list


  //----------------------------------------------> Reco and Only The Reco <--------------------------------------
  // Fake WC
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
  
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {art::fill_ptr_vector(wctrack, wctrackHandle);} 

  // Tracks
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
  // === Filling the tracklist from the tracklistHandle ===
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}


  // how many reco tracks exist out there?
  nReco = tracklist.size();
  
  // === Association between WC Tracks and TPC Tracks ===
  int TempTrackMatchedID = -1; 
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {
      art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
      if (fWC2TPC.isValid())
	{
	  // === Loop on all the Assn WC-TPC tracks === 
	  for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn ) 
	    {
	      // === Get the TPC track ===
	      cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));
	      if (!trackWC2TPC) continue;
	      recob::Track const& aTrack(trackWC2TPC.ref()); 
	      TempTrackMatchedID = aTrack.ID();
	    }//<----End indexAssn loop                                                                                                                       
	}//else std::cout<<"WC2TPC ass invalid"<<fWC2TPCModuleLabel<<std::endl;   
    } //else std::cout<<"Can't find WC "<<fWCTrackLabel<<std::endl;
  

  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  // Now TempTrackMatchedID contains the ID of the associated TPC track! YAY!
  // This is our guy!
  for(size_t iTrack=0; iTrack<tracklist.size();++iTrack)
    {
      //  for (auto recoTrk : tracklist )
      //{
      auto recoTrk = tracklist[iTrack];
      if (recoTrk->ID()  != TempTrackMatchedID) continue;
      // Let's make sure the track first point is actually it (directionality problem).
      auto realFirstValidPt = (recoTrk->TrajectoryPoint(recoTrk->FirstValidPoint())).position;
      auto realLastValidPt  = (recoTrk->TrajectoryPoint(recoTrk->LastValidPoint( ))).position;
      
      if ( realFirstValidPt.Z() - realLastValidPt.Z() > 0)
	{
	  auto bogusFirst  = realFirstValidPt;
	  realFirstValidPt = realLastValidPt;
	  realLastValidPt  = bogusFirst;
	}
      
      recoL    = recoTrk->Length();
      recoVtxX = (realFirstValidPt).X();
      recoVtxY = (realFirstValidPt).Y();
      recoVtxZ = (realFirstValidPt).Z();
      recoEndX = (realLastValidPt).X();
      recoEndY = (realLastValidPt).Y();
      recoEndZ = (realLastValidPt).Z();
      isTrackMatched        = 1;
      recoEnDeposited = 0.;

    
      vtxDist = TMath::Sqrt((recoVtxX - trueVtxX)*(recoVtxX - trueVtxX) +
			    (recoVtxY - trueVtxY)*(recoVtxY - trueVtxY) +
			    (recoVtxZ - trueVtxZ)*(recoVtxZ - trueVtxZ)) ;
    
      endDist = TMath::Sqrt((recoEndX - trueEndX)*(recoEndX - trueEndX) +
			    (recoEndY - trueEndY)*(recoEndY - trueEndY) +
			    (recoEndZ - trueEndZ)*(recoEndZ - trueEndZ)) ;

      if (fmcal.isValid())
	{
	  //	  std::cout<<"here\n";
	  // ### Putting calo information for this track (i) into pointer vector ###                                                                                                                           
	  std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(iTrack);
	  // ### Looping over each calorimetry point (similar to SpacePoint) ###                                                                                                                               
	  for (size_t j = 0; j<calos.size(); ++j)
	    {
	      // ### If we don't have calorimetry information for this plane skip ###                                                                                                                          
	      if (!calos[j]->PlaneID().isValid) continue;
	      // ### Skipping this point if the plane number doesn't make sense ###                                                                                                                            
	      if (calos[j]->PlaneID().Plane != 1) continue;
	      // ### Recording the number of calorimetry points for this track in this plane ####                                                                                                              
	      auto sizeDEdX = calos[j]->dEdx().size();
	      // #### Recording the kinetic energy for this track in this plane ###                                                                                                                            
	      keWhen = calos[j]->KineticEnergy();

	      // ###############################################                                                                                                                                               
	      // ### Looping over all the calorimetry points ###                                                                                            
	      // ###############################################                                                                                                                                               
	      for (size_t k = 0; k<sizeDEdX; ++k)  recoEnDeposited += calos[j]->dEdx()[k] * calos[j]->TrkPitchVec()[k];
	      
	    }//<---End looping over calo points (j)                           
	}
   

      recoInteractionEnergy = recoEnDeposited + 45.; // in MeV 
    }

  if (isTrackMatched){
    double minRecoTrueVtxDistance = 99999999.;
    //Find what is the true particle with the closest match
    for(size_t m = 0; m < plist.size(); ++m) 
      {
	auto matchedTruePart = plist.Particle(m);
	simb::MCTrajectory truetrajMatch = matchedTruePart->Trajectory();
	auto firstMatchTPCPoint =  truetrajMatch.begin();
	// Loop From First TrajPoint --> First Point in TPC 
	for ( auto t = truetrajMatch.begin(); t!= std::prev(truetrajMatch.end()); t++)
	  {
	    auto pos = t->first;
	    if (pos.Z() < minZ || pos.Z() > maxZ) continue;
	    else if (pos.X() < minX || pos.X() > maxX ) continue;
	    else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	    else {
	      firstMatchTPCPoint = t;
	      break;
	    }
	  }// End search for first point in TPC
	
	if (firstMatchTPCPoint == truetrajMatch.begin()) continue; // True particle not in TPC
	double distanceFromTrue =
	  ((firstMatchTPCPoint->first ).X() -  recoVtxX) * ((firstMatchTPCPoint->first ).X() -  recoVtxX) +
	  ((firstMatchTPCPoint->first ).Y() -  recoVtxY) * ((firstMatchTPCPoint->first ).Y() -  recoVtxY) +
	  ((firstMatchTPCPoint->first ).Z() -  recoVtxZ) * ((firstMatchTPCPoint->first ).Z() -  recoVtxZ) ;
	distanceFromTrue = TMath::Sqrt(distanceFromTrue);
	
	if (distanceFromTrue < minRecoTrueVtxDistance){
	  minRecoTrueVtxDistance = distanceFromTrue;
	  truePDGMatch = matchedTruePart->PdgCode();
	}
      }// End loop on particles
  }// End is Track Match

  bool debug = false;
  
  if (debug && isPrimaryInTPC && isTrackMatched && TMath::Abs(recoL-trueL) > 10. )
    {
      std::cout << std::fixed;
      std::cout << std::setprecision(2);
      std::cout<<"GBBB "<<recoL-trueL<<" "<< run<<" "<<subrun<<" "<<eventN<<"\n"
	       <<"GBBB TVtx: "<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<"\n"
	       <<"GBBB RVtx: "<<recoVtxX<<" "<<recoVtxY<<" "<<recoVtxZ<<"\n"
	     <<"GBBB TEnd: "<<trueEndX<<" "<<trueEndY<<" "<<trueEndZ<<"\n"
	     <<"GBBB REnd: "<<recoEndX<<" "<<recoEndY<<" "<<recoEndZ<<"\n\n\n";
    }
  fTree->Fill();
  
  
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  eventN  = evt.event();
  
  // Truth Variables
  isPrimaryInTPC        = 0;
  intAngle              = -999.;
  trueInteractionEnergy = -999.;  
  initialKE             = -999.;
  trueVtxX              = -999.;
  trueVtxY              = -999.;
  trueVtxZ              = -999.;
  trueEndX              = -999.;
  trueEndY              = -999.;
  trueEndZ              = -999.;
  trueL                 = -999.;  
  trueEnDeposited       = -999.;  
  //Reco Variables
  nReco                 = 0;
  isTrackMatched        = 0;
  recoVtxX              = -999.;
  recoVtxY              = -999.;
  recoVtxZ              = -999.;
  recoEndX              = -999.;
  recoEndY              = -999.;
  recoEndZ              = -999.;
  recoL                 = -999.; 
  recoEnDeposited       = -999.;  
  recoInteractionEnergy = -999.;  
  G4Process.clear();
  truePDGMatch         = -99999;
  
  
}

void lariat::TrackingEffWCMatch::endJob()
{
  std::cout<<"-------------------------------------------"<<std::endl;
  std::cout<<"True Events pre-TPC .............. "<<evtsPreTPC <<std::endl;
  std::cout<<"True Events post-TPC ............. "<<evtsPostTPC<<std::endl;
  std::cout<<"-------------------------------------------"<<std::endl;
}

void lariat::TrackingEffWCMatch::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&run      ,"run/I");
  fTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  fTree->Branch("eventN"   ,&eventN   ,"eventN/I");

  fTree->Branch("isPrimaryInTPC"    ,&isPrimaryInTPC ,"isPrimaryInTPC/I");
  fTree->Branch("isTrackMatched"    ,&isTrackMatched ,"isTrackMatched/I");

  fTree->Branch("recoInteractionEnergy"   ,&recoInteractionEnergy   ,"recoInteractionEnergy/D"); 
  fTree->Branch("trueInteractionEnergy"   ,&trueInteractionEnergy   ,"trueInteractionEnergy/D"); 
  fTree->Branch("trueEnDeposited"   ,&trueEnDeposited   ,"trueEnDeposited/D"); 
  fTree->Branch("recoEnDeposited"   ,&recoEnDeposited   ,"recoEnDeposited/D"); 
  fTree->Branch("trueKEFrontFace"   ,&trueKEFrontFace   ,"trueKEFrontFace/D"); 

  fTree->Branch("nReco"    ,&nReco    ,"nReco/I");
  fTree->Branch("initialKE",&initialKE,"initialKE/D"); 
  fTree->Branch("keWhen"   ,&keWhen   ,"keWhen/D"); 
  fTree->Branch("trueL"    ,&trueL    ,"trueL/D"); 
  fTree->Branch("recoL"    ,&recoL    ,"recoL/D"); 
  fTree->Branch("intAngle" ,&intAngle ,"intAngle/D");

  fTree->Branch("trueVtxX" ,&trueVtxX ,"trueVtxX/D");
  fTree->Branch("trueVtxY" ,&trueVtxY ,"trueVtxY/D");
  fTree->Branch("trueVtxZ" ,&trueVtxZ ,"trueVtxZ/D");
  fTree->Branch("trueEndX" ,&trueEndX ,"trueEndX/D");
  fTree->Branch("trueEndY" ,&trueEndY ,"trueEndY/D");
  fTree->Branch("trueEndZ" ,&trueEndZ ,"trueEndZ/D");

  fTree->Branch("recoVtxX" ,&recoVtxX ,"recoVtxX/D");
  fTree->Branch("recoVtxY" ,&recoVtxY ,"recoVtxY/D");
  fTree->Branch("recoVtxZ" ,&recoVtxZ ,"recoVtxZ/D");
  fTree->Branch("recoEndX" ,&recoEndX ,"recoEndX/D");
  fTree->Branch("recoEndY" ,&recoEndY ,"recoEndY/D");
  fTree->Branch("recoEndZ" ,&recoEndZ ,"recoEndZ/D");
  fTree->Branch("G4Process",&G4Process);

  fTree->Branch("vtxDist" ,&vtxDist ,"vtxDist/D");
  fTree->Branch("endDist" ,&endDist ,"endDist/D");
  fTree->Branch("truePDGMatch" ,&truePDGMatch ,"truePDGMatch/I");

}

DEFINE_ART_MODULE(lariat::TrackingEffWCMatch)

//  LocalWords:  GBBB
