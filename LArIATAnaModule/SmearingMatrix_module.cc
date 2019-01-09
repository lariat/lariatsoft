////////////////////////////////////////////////////////////////////////
// Class:       SmearingMatrix
// Module Type: analyzer
// File:        SmearingMatrix_module.cc
//
// Generated at Tue Oct 3 10:39:46 2018 by Elena Gramellini
// Final goal: produced the smearing matrix for N_Int and N_Inc
//
// Deliverables:
// [   ] ttree with: trueIntKE, recoIntKE, trueIncKE_v, trueIncKE_v
//
// N_Int
// [   ] Did I had a true interaction?
//       [   ] If yes, what was its energy?
// [   ] Did I reco an interaction? 
//       [   ] If yes, what was its energy?
//
// N_Inc
// [   ] Did I had a true incoming pion at this slice?
//       [   ] If yes, what was its energy?
// [   ] Did I reco an incoming pion at this slice?
//       [   ] If yes, what was its energy?
//
// To Do:
// [   ] define the tree
//       [  ] true interaction ke
//       [  ] reco interaction ke
//       [  ] true incident ke
//       [  ] reco incident ke
// [  ] reset variables
//
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


namespace lariat 
{
  class SmearingMatrix;
}

class lariat::SmearingMatrix : public art::EDAnalyzer 
{
public:
  explicit SmearingMatrix(fhicl::ParameterSet const & p);
  virtual ~SmearingMatrix();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const & p);
  double energyLossCalculation(double, double, bool);
  double distance(double, double, double, double, double, double );
  bool   isWithinActiveVolume(double, double, double);

private:
  void ResetVars();


  // True
  TH2D *hInteractionPointZX ;
  TH2D *hInteractionPointZY ;
  TH2D *hInitialPointZX ;
  TH2D *hInitialPointZY ;
  TH2D *hFinalPointZX ;
  TH2D *hFinalPointZY ;
  TH1D *hSimIDE_dEdX;
  TH1D *hSimIDE_dE;
  TH1D *hSimIDE_dX;
  TH1D *hDiffG4KE_SimIDE_El;
  TH1D *hDiffG4KE_SimIDE_Inel;
  TH1D *hDiffG4KE_SimIDE_Thro;
  // Reco
  TH2D *hRecoInteractionPointZX ;
  TH2D *hRecoInteractionPointZY ;
  TH2D *hRecoInitialPointZX ;
  TH2D *hRecoInitialPointZY ;
  TH2D *hRecoFinalPointZX ;
  TH2D *hRecoFinalPointZY ;
  TH1D *hReco_dEdX;
  TH1D *hReco_dE;
  TH1D *hReco_dX;




  TTree* tTree;
  int    run;
  int    subrun;
  int    event;
  int    nSlices  = 300;
  int    trueNPoints;
  int    recoNPoints;
  double trueLength;
  double trueTPCFFKE;
  double trueVtxX  = -999.;
  double trueVtxY  = -999.;
  double trueVtxZ  = -999.;
  double trueEndX  = -999.;
  double trueEndY  = -999.;
  double trueEndZ  = -999.;
  double finalKEG4 = -999.;

  double recoLenght;
  double recoTPCFFKE;
  double recoVtxX  = -999.;
  double recoVtxY  = -999.;
  double recoVtxZ  = -999.;
  double recoEndX  = -999.;
  double recoEndY  = -999.;
  double recoEndZ  = -999.;
  double trueIntKE ;
  double recoIntKE ;
  double trueIncKE_v[300] ;
  double recoIncKE_v[300] ;


  // True boundaries
  double minX =   0.;
  double maxX =  47.;
  double minY = -20.;
  double maxY =  20.;
  double minZ = -0.1;
  double maxZ =  90.;

  int evtsPreTPC  = 0;
  int evtsPostTPC = 0;  
  int interactingInTPC  = 0;
  int nEvtLongEnough    = 0;
  int nEvtMatchedTracks = 0;
  int matchedTracks     = 0;
  int validCaloOnColl   = 0;
  bool isTrackInteracting = true;

  bool fisData = false;

  //Fiducial Volume
  double fminX =  1.0; //0.
  double fmaxX = 46.0; //47
  double fminY =-15.0; //20
  double fmaxY = 15.0; //20
  double fminZ =  2.0; // G10 at z=1.2
  double fmaxZ = 86.0;


  double hitDEDXThreshold = 40.; 

  std::string fTrackModuleLabel; 
  std::string fCalorimetryModuleLabel;
  std::string fWCTrackLabel; 
  std::string fWC2TPCModuleLabel; 
  double rMass;

};


bool lariat::SmearingMatrix::isWithinActiveVolume(double x, double y, double z)
{
  if (x < minX ) return false; 
  if (x > maxX ) return false;
  if (y < minY ) return false; 
  if (y > maxY ) return false;
  if (z < minZ ) return false; 
  if (z > maxZ ) return false;
  return true;
}

lariat::SmearingMatrix::SmearingMatrix(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::SmearingMatrix::~SmearingMatrix()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::SmearingMatrix::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel             = pset.get< std::string >("TrackModuleLabel"      , "pmtrack");
  fCalorimetryModuleLabel       = pset.get< std::string >("CalorimetryModuleLabel", "calo"     );
  fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel"          , "wctrack"  );
  fWC2TPCModuleLabel      	= pset.get< std::string >("WC2TPCModuleLabel"     , "wctracktpctrackmatch");
  fisData                	= pset.get<  bool  >("isData"    , false);
  rMass                         = pset.get< double >("fmass"     , 139.57018);
}

double lariat::SmearingMatrix::energyLossCalculation(double x, double px, bool isData)
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

void lariat::SmearingMatrix::analyze(art::Event const & evt)
{
  bool verbose = false;
  ResetVars();      
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  event   = evt.event();

  //////////////////////////  
  //         Truth
  //////////////////////////  
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
 
  // Get the backtracker to recover true quantities
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  //Scope 1.: 
  // [] find an interaction that we can reconstruct
  // [] record its energy

  
  //bool keepInteraction = false;
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
      

      //Get simIDE associated to the primary
      geo::View_t view = geom->View(0); 
      auto simIDE_Prim  = bt->TrackIdToSimIDEs_Ps(mcPart->TrackId(),view);
      // Order them in order of increasing Z
      std::map<double, sim::IDE> orderedSimIDE;
      for (auto ide : simIDE_Prim ) orderedSimIDE[ide->z] = *ide;


      //Let's store the interaction type, so we can use it later to detemine if we want to keep this interaction
      std::string interactionLabel = "";
      
      // Get the mass (in GeV)
      double mass = mcPart->Mass() ;
      if (mass > 0.140 || mass < 0.130 ) {std::cout<<"process"<<proc<<" mass: "<<mass<<"Z ini "<<( (truetraj.begin())->first).Z()   <<" These are not pions \n"; return;} // 


      //Store the kinetic energy and momentum on z at WC4. Just for cross check 
      auto inTPCPoint  = truetraj.begin(); 
            
      //--------------------------------------------------------
      // Identify the first trajectory point inside the TPC
      // Loop From First TrajPoint --> First Point in TPC 
      // Stop when you get into the TPC
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  //std::cout<<pos.X() <<" "<<pos.Y() <<" "<<pos.Z() <<"\n";
	  if (pos.Z() < minZ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    inTPCPoint = t;
	    break;
	  }
	}// End search for first point in TPC
      

      
      // if the first point is not in the TPC, we're not interested in this event. 
      evtsPreTPC++;

      // If the pion didn't even arrive to the TPC, we don't care about it
      if (inTPCPoint == truetraj.begin()) continue;
      // Fill some histos for cross check
      evtsPostTPC++;     
      hInitialPointZX->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).X());
      hInitialPointZY->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).Y());
      
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
	      if ((truetraj.KeyToProcess(couple.second)).find("CoulombScat")!= std::string::npos) continue;
	      
	      // Let's check if the interaction is in the the TPC, 
	      // if it's not, I dont' care about it: continue!
	      auto     interactionPos4D =  (truetraj.at(couple.first)).first ;	
	      if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
	      else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
	      else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;

	      // Is it a interaction that I could see? 
	      double interactionAngle = 999999.; // This needs to be changed
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
		  interactionAngle = TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  );
		  //if (debug) std::cout<<"hopefully this is elastic. "<<truetraj.KeyToProcess(couple.second)<<"\n Angle: "<<interactionAngle<<"\n";
		}else
		{ // The particle has come to an end. Let's check the daugthers.
		  if (mcPart->NumberDaughters() == 0 ) break;
		  double maxAngle = 0.;
		  int numberOfTroubleDaugther = 0;
		  for(size_t d = 0; d < plist.size(); ++d) 
		    {
		      auto mcDaught = plist.Particle(d);
		      if (mcDaught->Mother()  != 1 ) continue;
		      auto daugthTraj = mcDaught->Trajectory();
		      
		      
		      if (mcDaught->NumberTrajectoryPoints () < 2 ) continue ;
		      if (!(TMath::Abs(mcDaught->PdgCode())  == 13  ||
			    TMath::Abs(mcDaught->PdgCode())  == 11  ||
			    TMath::Abs(mcDaught->PdgCode())  == 211 ||
			    TMath::Abs(mcDaught->PdgCode())  == 321 || 
			    TMath::Abs(mcDaught->PdgCode())  == 2212) ) continue;
		      if (daugthTraj.TotalLength() < 0.5 ) continue ; // I'm not going to reconstruct this daughter. No need to worry about her.
		      
		      numberOfTroubleDaugther++;	
		      
		      auto daughtFirstPt  =  ((daugthTraj[0]).first).Vect() ;
		      auto daughtSecondPt =  ((daugthTraj[1]).first).Vect() ;
		      auto distanceBtwPointNext = daughtSecondPt - daughtFirstPt;
		      interactionAngle = TMath::ACos(distanceBtwPointNext.Dot(distanceBtwPoint)/(distanceBtwPointNext.Mag()*distanceBtwPoint.Mag() )  );
		      if ( maxAngle < interactionAngle ) maxAngle = interactionAngle;
		    }
		  interactionAngle = maxAngle;
		  // If the track finishes without visible daugthers, we're likely to see it. Give it a big angle!
		  if (!numberOfTroubleDaugther) interactionAngle = 9999999.;//It's a huge angle: it's greater than the cut, so I don't remove this! But it also doesn't go into the angle plot
		  numberOfTroubleDaugther = 0;
		}
	    
	      
	      // If we made it here, then this is the first interesting interaction in the TPC
	      // Our job is done!!! Great! Store the interaction label and the iterator for the final point
	      interactionLabel = truetraj.KeyToProcess(couple.second);
	      finTPCPoint = truetraj.begin() + couple.first; 
	      //keepInteraction = true;
	      interactingInTPC++;
	      interactionAngle = 999999.;
	      break;
	    }// Loop on interaction points
	}// If there are G4 interactions      
      
      if (!interactionLabel.size()) 
	{
	  auto tempTPCEnd  = std::prev(truetraj.end()); 
	  for ( auto t = std::prev(truetraj.end()); t!= truetraj.begin(); t--)
	    {
	      auto pos = t->first;
	      //std::cout<<pos.X() <<" "<<pos.Y() <<" "<<pos.Z() <<"\n";
	      if (pos.Z() > maxZ) continue;
	      else if (pos.X() < minX || pos.X() > maxX ) continue;
	      else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	      else {
		tempTPCEnd = t;
		break;
	      }
	    }// End search for first point in TPC
	  finTPCPoint = tempTPCEnd;
	}
 
      if (verbose) std::cout<<interactionLabel<<" "<<(finTPCPoint->first).Z()<<" , "; 
      
      if (finTPCPoint == inTPCPoint) continue;
      auto posFin = finTPCPoint->first;
      auto posIni = inTPCPoint->first;
      //Let's record what the initial and final points are.
      auto totLength = distance(posFin.X(), posFin.Y(), posFin.Z(),posIni.X(), posIni.Y(), posIni.Z() );
      if (totLength < 2.0 ) continue;
      nEvtLongEnough++;

      hFinalPointZX->Fill(posFin.Z(), posFin.X());
      hFinalPointZY->Fill(posFin.Z(), posFin.Y());


      trueLength = totLength;
      trueVtxX   = posIni.X();
      trueVtxY   = posIni.Y();
      trueVtxZ   = posIni.Z();
      trueEndX   = posFin.X();
      trueEndY   = posFin.Y();
      trueEndZ   = posFin.Z();
      
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
      int    nPts            = (int) (totLength/trackPitch);
      for (int iPt = 1; iPt <= nPts; iPt++ )
	{
	  auto newPoint = positionVector0 + iPt*(trackPitch/totLength) * (positionVector1 - positionVector0);
	  orderedUniformTrjPts[newPoint.Z()] = newPoint;
	}


      // If the distance between the last point and the second to last is less then 0.235
      // eliminate the second to last point

      auto lastPt         = (orderedUniformTrjPts.rbegin())->second;
      auto secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
      double lastDist = distance(lastPt.X(),lastPt.Y(),lastPt.Z(),secondtoLastPt.X(),secondtoLastPt.Y(),secondtoLastPt.Z());

      if (lastDist < 0.235)
	{
	  orderedUniformTrjPts.erase((std::next(orderedUniformTrjPts.rbegin()))->first );
	} 

      
      // Ok, let's calculate the KE at the interaction point
      // and the true enrgy at each slice 
      // Calculate the initial kinetic energy
      auto initialMom =     inTPCPoint->second;
      double initialKE = 1000*(TMath::Sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y() + initialMom.Z()*initialMom.Z() + mass*mass ) - mass); 
      double kineticEnergy = initialKE;
      trueTPCFFKE = initialKE;

      trueNPoints = 0;
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
	  if (currentDepEnergy/uniformDist < 0.1 )   continue;

	  trueIncKE_v[trueNPoints] = kineticEnergy;	  
	  hSimIDE_dEdX->Fill(currentDepEnergy/uniformDist);
	  hSimIDE_dX  ->Fill(uniformDist);
	  hSimIDE_dE  ->Fill(currentDepEnergy);
	  
	  //Calculate the current kinetic energy
	  kineticEnergy -= currentDepEnergy;
	  trueNPoints++;
	}// Loop on ordered Points
      
      trueNPoints++;

      
      auto finalMomG4  =     finTPCPoint->second;
      finalKEG4 = 1000*(TMath::Sqrt(finalMomG4.X()*finalMomG4.X() + finalMomG4.Y()*finalMomG4.Y() + finalMomG4.Z()*finalMomG4.Z() + mass*mass ) - mass); 
      trueIntKE  = kineticEnergy; 


      //We're interested in Inelastic and Elastic Interacting with the last point
      if ( interactionLabel.find("Inelastic")!= std::string::npos ) {
	if (verbose) std::cout<< "Inelastic "<<finalKEG4<<" "<<trueIntKE<<"\n\n ";
	hDiffG4KE_SimIDE_Inel->Fill(finalKEG4 - trueIntKE);
	//trueIntKE = finalKEG4;
      }
      if ( interactionLabel.find("Elastic")  != std::string::npos ) {
	hDiffG4KE_SimIDE_El->Fill(finalKEG4 - trueIntKE);
	if (verbose) std::cout<< "Elastic "<<finalKEG4<<" "<<trueIntKE<<"\n\n ";
      }
      if (!interactionLabel.size()){ 
	interactionLabel = "throughgoing";
	hDiffG4KE_SimIDE_Thro->Fill(finalKEG4 - trueIntKE);
	//std::cout<<interactionLabel<<" "<<trueIntKE<<" "<<recoIntKE<<"\n";
	//trueIntKE = -999.;
	//if ( finalKEG4 - trueIntKE < -2 || finalKEG4 - trueIntKE > 2) std::cout<<(inTPCPoint->first).Z() <<" " << (finTPCPoint->first).Z() <<" "<<finalKEG4 - trueIntKE<<" "<< finalKEG4<<" "<< trueIntKE <<" \n";
	if (verbose) std::cout<< "Through "<<finalKEG4<<" "<<trueIntKE<<"\n\n ";
      }      
    }//MC Particles



  //////////////////////////  
  //         Reco
  //////////////////////////  
  
  
  // Let's take some Reco stuff...
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;   // Container of wc tracks  
  std::vector<art::Ptr<ldp::WCTrack> >     wctrack;         // Vector of wc tracks  
  art::Handle< std::vector<recob::Track> > trackListHandle; // Container of reco tracks
  std::vector<art::Ptr<recob::Track> >     tracklist;       // Vector of wc tracks  

  // Find which recontructed track is associated with the WC
  int matchedRecoTrkKey = -99999;
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;

  if(!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return;
  //nevts++;

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
	  //std::cout<<"point "<<realFirstValidPt.Z()<<" "<<realLastValidPt.Z()<<" "<<invertedTracking<<"\n";
	  auto bogusFirst  = realFirstValidPt;
	  realFirstValidPt = realLastValidPt;
	  realLastValidPt  = bogusFirst;
	}
      


      recoVtxX = realFirstValidPt.X();
      recoVtxY = realFirstValidPt.Y();
      recoVtxZ = realFirstValidPt.Z();
      recoEndX = realLastValidPt.X();
      recoEndY = realLastValidPt.Y();
      recoEndZ = realLastValidPt.Z();

      hRecoInitialPointZX->Fill(realFirstValidPt.Z(),realFirstValidPt.X());
      hRecoInitialPointZX->Fill(realFirstValidPt.Z(),realFirstValidPt.Y());
      hRecoFinalPointZX  ->Fill(realLastValidPt.Z(),realLastValidPt.X());
      hRecoFinalPointZY  ->Fill(realLastValidPt.Z(),realLastValidPt.Y());
  
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
	      
	      hReco_dEdX->Fill(calos[j]->dEdx()[k]);
	      hReco_dE  ->Fill(calos[j]->dEdx()[k] * calos[j]->TrkPitchVec()[k]);
	      hReco_dX  ->Fill(calos[j]->TrkPitchVec()[k]);

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


  double WCMom = wctrack[0]->Momentum();
  double WC4X  = wctrack[0]->HitPosition(3,0);


  if (!recoDEDX_v.size()) return; // If there are no points, return
  if (countIDmatch != 1) return; // If for some reason I couldn't find the match track, I should return
  if (recoEndX < -100. || recoEndY < -100. || recoEndZ < -100.) throw cet::exception("Geometry")    /// <== ... throw exception
								  << "Very Weird End of track... somethings going on, may I smell... \n";
  if (WCMom < 0.) throw cet::exception("Geometry")    /// <== ... throw exception
		    << "WCMom < 0... somethings going on, may I smell... \n";
  
  if (rMass < 0.) throw cet::exception("Geometry")    /// <== ... throw exception
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
  double WCKE             = TMath::Sqrt(WCMom*WCMom + rMass*rMass) - rMass;  
  //How much energy is lost before getting to the TPC?
  double calculatedEnLoss = 0; 
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
      for(size_t p = 0; p < plist.size(); ++p)
	{
	  auto mcPart = plist.Particle(p);
	  std::string proc = mcPart->Process();
	  if ( !(proc.find("primary") != std::string::npos) ) continue;
	  simb::MCTrajectory truetraj = mcPart->Trajectory();
	  // Make Sure we get the beamline primary                                                                                                                                                
	  if ( ( (truetraj.begin())->first).Z() >  -50. ) continue;
	  auto inTPCPoint  = truetraj.begin();
	  auto Momentum0   = inTPCPoint->second;
	  double onTheFlyPx = 1000*Momentum0.X(); // the Momentum Needs To be in MeV
	  calculatedEnLoss = energyLossCalculation(WC4X,onTheFlyPx,fisData);
	}
    }

  //Initial energy at the TPC Front Face. 
  const double initialKEReco  = WCKE - calculatedEnLoss;
  recoIncidentKE_v.push_back(initialKEReco);

  recoNPoints = (int)recoIncidentKE_v.size();
  // At this point, I assume the energy deposited vector are filled in the right direction
  double energyDepositedThusFar = 0.;
  for (size_t iDep = 0; iDep < recoEDep_v.size(); iDep++)
    {
      energyDepositedThusFar += recoEDep_v[iDep];
      recoIncidentKE_v.push_back(initialKEReco-energyDepositedThusFar);
    }

  //Fill the Incident and Interacting plots real quick
  for (size_t iKE = 0; iKE <   recoIncidentKE_v.size(); iKE++) recoIncKE_v[iKE] = recoIncidentKE_v[iKE];
  if (isTrackInteracting && recoIncidentKE_v.size())  recoIntKE =  recoIncidentKE_v[recoIncidentKE_v.size()-1]; 

  recoTPCFFKE = recoIncKE_v[0];  
  double sumOfPitchs = 0;
  for (size_t i = 0; i < recoDEDX_v.size(); i++) sumOfPitchs += recoPitch_v[i];
  recoLenght = sumOfPitchs;

  tTree->Fill();
}


void lariat::SmearingMatrix::endJob()
{
  std::cout<<"Evts        preTPC "<< evtsPreTPC  <<"\n";
  std::cout<<"Evts       postTPC "<< evtsPostTPC <<"\n";  
  std::cout<<"Evts   long enough "<< nEvtLongEnough <<"\n";  
  std::cout<<"Interacting In TPC "<< interactingInTPC <<"\n";
  std::cout<<"Evts with matches  "<< matchedTracks <<"\n";
}


void lariat::SmearingMatrix::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hInteractionPointZX   = tfs->make<TH2D>("hIntPointZX","hIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hInteractionPointZY   = tfs->make<TH2D>("hIntPointZY","hIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hInitialPointZX       = tfs->make<TH2D>("hInitialPointZX","hIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hInitialPointZY       = tfs->make<TH2D>("hInitialPointZY","hIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hFinalPointZX         = tfs->make<TH2D>("hFinalPointZX","hIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hFinalPointZY         = tfs->make<TH2D>("hFinalPointZY","hIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hSimIDE_dEdX          = tfs->make<TH1D>("hSimIDE_dEdX"       ,"dEdX SimIDE; dEdX [MeV/cm];", 200, 0,20);
  hSimIDE_dE            = tfs->make<TH1D>("hSimIDE_dE"         ,"dE SimIDE; dE [MeV/cm];", 200, 0,20);
  hSimIDE_dX            = tfs->make<TH1D>("hSimIDE_dX"         ,"dX SimIDE; dX [MeV/cm];", 200, 0,2);
  hDiffG4KE_SimIDE_El   = tfs->make<TH1D>("hDiffG4KE_SimIDE_El"  ,"Elastic    KE diff; G4 - SimIDE KE Elastic [MeV];", 200,-100,100);
  hDiffG4KE_SimIDE_Inel = tfs->make<TH1D>("hDiffG4KE_SimIDE_Inel","Inelast    KE diff; G4 - SimIDE KE InEla   [MeV];", 200,-1000,1000);
  hDiffG4KE_SimIDE_Thro = tfs->make<TH1D>("hDiffG4KE_SimIDE_Thro","Througoing KE diff; G4 - SimIDE KE InEla   [MeV];", 200,-100,100);

  hRecoInteractionPointZX   = tfs->make<TH2D>("hRecoIntPointZX","hRecoIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hRecoInteractionPointZY   = tfs->make<TH2D>("hRecoIntPointZY","hRecoIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hRecoInitialPointZX       = tfs->make<TH2D>("hRecoInitialPointZX","hRecoIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hRecoInitialPointZY       = tfs->make<TH2D>("hRecoInitialPointZY","hRecoIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hRecoFinalPointZX         = tfs->make<TH2D>("hRecoFinalPointZX","hIntPointZX;Z [cm]; X [cm]", 220,-10,100,120,-10,50); 
  hRecoFinalPointZY         = tfs->make<TH2D>("hRecoFinalPointZY","hIntPointZY;Z [cm]; Y [cm]", 220,-10,100,120,-30,30); 
  hReco_dEdX                = tfs->make<TH1D>("hReco_dEdX"       ,"dEdX SimIDE; dEdX [MeV/cm];", 200, 0,20);
  hReco_dE                  = tfs->make<TH1D>("hReco_dE"         ,"dE SimIDE; dE [MeV/cm];", 200, 0,20);
  hReco_dX                  = tfs->make<TH1D>("hReco_dX"         ,"dX SimIDE; dX [MeV/cm];", 200, 0,2);



  tTree = tfs->make<TTree>("smearTree","analysis tree");
  tTree->Branch("run"      ,&run      ,"run/I");
  tTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  tTree->Branch("event"    ,&event    ,"event/I");
  tTree->Branch("trueNPoints"    ,&trueNPoints    ,"trueNPoints/I");
  tTree->Branch("recoNPoints"    ,&recoNPoints    ,"recoNPoints/I");
  tTree->Branch("trueIntKE"  ,&trueIntKE   ,"trueIntKE/D"  );
  tTree->Branch("recoIntKE"  ,&recoIntKE   ,"recoIntKE/D"  );
  tTree->Branch("nSlices"    ,&nSlices     ,"nSlices/I");  
  tTree->Branch("trueIncKE"  , trueIncKE_v ,"trueIncKE[nSlices]/D"); 
  tTree->Branch("recoIncKE"  , recoIncKE_v ,"recoIncKE[nSlices]/D"); 

  tTree->Branch("trueLength"  ,&trueLength   ,"trueLength/D"  );
  tTree->Branch("trueTPCFFKE" ,&trueTPCFFKE  ,"trueTPCFFKE/D"  );
  tTree->Branch("trueVtxX"    ,&trueVtxX     ,"trueVtxX/D"  );
  tTree->Branch("trueVtxY"    ,&trueVtxY     ,"trueVtxY/D"  );
  tTree->Branch("trueVtxZ"    ,&trueVtxZ     ,"trueVtxZ/D"  );
  tTree->Branch("trueEndX"    ,&trueEndX     ,"trueEndX/D"  );
  tTree->Branch("trueEndY"    ,&trueEndY     ,"trueEndY/D"  );
  tTree->Branch("trueEndZ"    ,&trueEndZ     ,"trueEndZ/D"  );
  tTree->Branch("finalKEG4"   ,&finalKEG4    ,"finalKEG4/D"  );


  tTree->Branch("recoLenght"  ,&recoLenght   ,"recoLenght/D"  );
  tTree->Branch("recoTPCFFKE" ,&recoTPCFFKE  ,"recoTPCFFKE/D"  );
  tTree->Branch("recoVtxX"    ,&recoVtxX     ,"recoVtxX/D"  );
  tTree->Branch("recoVtxY"    ,&recoVtxY     ,"recoVtxY/D"  );
  tTree->Branch("recoVtxZ"    ,&recoVtxZ     ,"recoVtxZ/D"  );
  tTree->Branch("recoEndX"    ,&recoEndX     ,"recoEndX/D"  );
  tTree->Branch("recoEndY"    ,&recoEndY     ,"recoEndY/D"  );
  tTree->Branch("recoEndZ"    ,&recoEndZ     ,"recoEndZ/D"  );


}

double lariat::SmearingMatrix::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return d;
}

void lariat::SmearingMatrix::ResetVars()
{
  run      = 0;
  subrun   = 0;
  event    = 0 ;
  nSlices  = 300;
  trueNPoints = 0;
  recoNPoints = 0;
  trueLength = -999.;
  trueTPCFFKE= -999.;
  trueVtxX  = -999.;
  trueVtxY  = -999.;
  trueVtxZ  = -999.;
  trueEndX  = -999.;
  trueEndY  = -999.;
  trueEndZ  = -999.;
  finalKEG4 = -999.;

  recoLenght= -999.;
  recoTPCFFKE= -999.;
  recoVtxX  = -999.;
  recoVtxY  = -999.;
  recoVtxZ  = -999.;
  recoEndX  = -999.;
  recoEndY  = -999.;
  recoEndZ  = -999.;
  trueIntKE = -999.;
  recoIntKE = -999.;


  for (int i = 0; i < nSlices; i++) {trueIncKE_v[i] = -999.; recoIncKE_v[i] = -999.;}
}

DEFINE_ART_MODULE(lariat::SmearingMatrix)
