////////////////////////////////////////////////////////////////////////
// Class:       TrueXSAngleCut
// Module Type: analyzer
// File:        TrueXSAngleCut_module.cc
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
  class TrueXSAngleCut;
}

class lariat::TrueXSAngleCut : public art::EDAnalyzer 
{
public:
  explicit TrueXSAngleCut(fhicl::ParameterSet const & p);
  virtual ~TrueXSAngleCut();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const & p);
  double  distance(double, double, double, double, double, double );
  double  calculateAngle(double, double, double, double, double, double );

private:
  bool    debug = false;
  TH1D*   h_AngleEl;
  TH1D*   h_AngleInel;
  TH1D*   h_AngleEl_AfterCut;
  TH1D*   h_AngleInel_AfterCut;


  TH1D*   h_DE   ;
  TH1D*   h_DX   ;
  TH1D*   h_DEDX ;


  TH1D*   h_DEUniform   ;
  TH1D*   h_DXUniform   ;
  TH1D*   h_DEDXUniform ;


  TH1D*   h_DeltaE ;

  TH1D*   h_SimIDEDist ;

  TH1D*   h_UniformDistances ;

  TH1D *hInteractingKE; 
  TH1D *hInteractingKEEl; 
  TH1D *hInteractingKEElDep; 
  TH1D *hInteractingKEInel; 
  TH1D *hIncidentKE; 
  TH1D *hCrossSection;
  TH1D *hCrossSectionEl;
  TH1D *hCrossSectionInel;

  TH1D *hKEAtTPCFF; 
  TH1D *hInitialKE; 
  TH1D *hInitialPz; 


  TH2D *hXZ; 
  TH2D *hYZ; 
  TH2D *hXZPre; 
  TH2D *hYZPre; 

  TH2D *hdEVsdX; 
  TH2D *hdEVsKE; 

  
  TTree* fTree;
  int    run;
  int    subrun;
  int    eventN;

  double trueVtxX ;
  double trueVtxY ;
  double trueVtxZ ;
  double trueEndX ;
  double trueEndY ;
  double trueEndZ ;

  double finalKE ;
  std::vector<std::string> G4Process; 


  double angleCut = 0.08334; // MC Value: 0.08334, Data Value: 0.09248
  double minX =  0.0;
  double maxX = 47.0;
  double minY =-20.0;
  double maxY = 20.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 90.0;



  int evtsPreTPC  = 0;
  int evtsInTheMiddle = 0;
  int evtsPostTPC = 0;
  int throughgoing = 0;
  int interactingInTPC = 0;
};


lariat::TrueXSAngleCut::TrueXSAngleCut(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::TrueXSAngleCut::~TrueXSAngleCut()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::TrueXSAngleCut::reconfigure(fhicl::ParameterSet const & pset)
{
  angleCut = pset.get< double >("angle_Cut",  0.08334 ); // MC Value: 0.08334, Data Value: 0.09248 
  minX     = pset.get< double >("minX_Cut" ,  0.0 );
  maxX     = pset.get< double >("maxX_Cut" , 47.0 );
  minY     = pset.get< double >("minY_Cut" ,-20.0 );
  maxY     = pset.get< double >("maxY_Cut" , 20.0 );
  minZ     = pset.get< double >("minZ_Cut" ,  0.0 );
  maxZ     = pset.get< double >("maxZ_Cut" , 90.0 );
}

void lariat::TrueXSAngleCut::analyze(art::Event const & evt)
{
  if (debug) std::cout<<"\n";
  bool verbose = false;
  // std::cout<<"\n\n -------------------- new event\n";
  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
 
  // Get the backtracker to recover true quantities
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

 
  //Set a bogus branch quantities, to be overwritten  
  run     = evt.run(); 
  subrun  = evt.subRun(); 
  eventN  = evt.event();

  trueVtxX = -999.;
  trueVtxY = -999.;
  trueVtxZ = -999.;
  trueEndX = -999.;
  trueEndY = -999.;
  trueEndZ = -999.;
  

  bool keepInteraction = false;
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


      //Let's store the interaction type, so we can use it later in the Cross Section
      std::string interactionLabel = "";
      // Get the mass (in GeV)
      double mass = mcPart->Mass() ;
      if (verbose)  std::cout<<"mass "<<mass<<"\n";

      
      //Store the kinetic energy and momentum on z at WC4. Just for cross check 
      auto inTPCPoint  = truetraj.begin(); 
      auto Momentum0   = inTPCPoint->second;
      
      double KE0 = 1000*(TMath::Sqrt(Momentum0.X()*Momentum0.X() + Momentum0.Y()*Momentum0.Y() + Momentum0.Z()*Momentum0.Z() + mass*mass ) - mass); //I want this in MeV
      hInitialKE->Fill(KE0);
      hInitialPz->Fill(1000*Momentum0.Z());


      //--------------------------------------------------------
      // Identify the first trajectory point inside the TPC
      // Loop From First TrajPoint --> First Point in TPC 
      // Stop when you get into the TPC
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  if (pos.Z() < minZ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    inTPCPoint = t;
	    break;
	  }
	}// End search for first point in TPC
      
      // if the first point is not in the TPC, we're not interested in this event. 
      // Fill some histos for cross check
      evtsPreTPC++;
      hXZPre->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).X());
      hYZPre->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).Y());

      if (inTPCPoint == truetraj.begin()) continue;
      evtsPostTPC++;     
      hXZ   ->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).X());
      hYZ   ->Fill((inTPCPoint->first).Z(), (inTPCPoint->first).Y());
      
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
	      if (verbose) std::cout<<(truetraj.KeyToProcess(couple.second))<<" Position "<< ((truetraj.at(couple.first)).first).Z() <<"\n";
	      if ((truetraj.KeyToProcess(couple.second)).find("CoulombScat")!= std::string::npos) continue;
	      
	      // Let's check if the interaction is in the the TPC, 
	      // if it's not, I dont' care about it: continue!
	      auto     interactionPos4D =  (truetraj.at(couple.first)).first ;	
	      if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
	      else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
	      else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
	      if (debug) std::cout<<"Interaction in TPC, great! \n Let's check the angle\n";
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
		      if (debug) std::cout<< mcDaught->PdgCode()
					  <<" , length: "<< (mcDaught->Trajectory()).TotalLength() <<"\n";
				   //<<"First Point: "<< ((daugthTraj[0].first).Vect()).X() <<","<< ((daugthTraj[0].first).Vect()).Y()<<","<<((daugthTraj[0].first).Vect()).Z() <<"\n";

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
		  if (debug) {std::cout<<"Opefully this is inelastic. "<<truetraj.KeyToProcess(couple.second)<<"\n Angle: "<<interactionAngle<<" "<<numberOfTroubleDaugther<<"\n";
		    if (interactionAngle < 0.0001) std::cout<<"-------------------------------------------------> \n";
		  }
		  numberOfTroubleDaugther = 0;
		}
	      

	      if ((truetraj.KeyToProcess(couple.second)).find("Elastic")!= std::string::npos) h_AngleEl->Fill(interactionAngle);
	      else h_AngleInel->Fill(interactionAngle);	      
	      if (interactionAngle < angleCut )  continue;
	      
	      if ((truetraj.KeyToProcess(couple.second)).find("Elastic")!= std::string::npos) h_AngleEl_AfterCut->Fill(interactionAngle); 
	      else h_AngleInel_AfterCut->Fill(interactionAngle); 

	      
	      // If we made it here, then this is the first interesting interaction in the TPC
	      // Our job is done!!! Great! Store the interaction label and the iterator for the final point
	      interactionLabel = truetraj.KeyToProcess(couple.second);
	      finTPCPoint = truetraj.begin() + couple.first; 
	      keepInteraction = true;
	      interactingInTPC++;
	      interactionAngle = 999999.;
	      break;
	    }// Loop on interaction points
	}// If there are G4 interactions
      

      // If I didn't find anything interesting in the intereaction map, let's loop back!
      if ( !keepInteraction )
	{
	  //Loop on the daughters 
	  for(size_t d = 0; d < plist.size(); ++d) 
	    {
	      auto mcDaught = plist.Particle(d);
	      //We keep only the dauthers of the primary not coming from elastic or inelastic scattering
	      if (mcDaught->Mother()  != 1 ) continue;
	      if ((mcDaught->Process()).find("astic")!= std::string::npos) continue;
	      if ((mcDaught->Process()).find("CoulombScat")!= std::string::npos) continue;

	      //Is the daughter born inside the TPC? If yes, store the process which created it 
	      simb::MCTrajectory trueDaugthTraj = mcDaught->Trajectory();	      
	      if (trueDaugthTraj.begin()->first.Z() < minZ || trueDaugthTraj.begin()->first.Z() > maxZ) continue;
	      else if (trueDaugthTraj.begin()->first.X() <   minX || trueDaugthTraj.begin()->first.X() > maxX ) continue;
	      else if (trueDaugthTraj.begin()->first.Y() <   minY || trueDaugthTraj.begin()->first.Y() > maxY ) continue;
	      else {
		interactionLabel = mcDaught->Process();
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
		//std::cout<<"Daugthers: "<<interactionLabel<<" at ("<< pos.X()<<","<<pos.Y()<<","<<pos.Z() <<")\n";
		finTPCPoint = t;
		break;
	      }
	    }
	}      
 
      if (finTPCPoint == inTPCPoint) continue;
      auto posFin = finTPCPoint->first;
      auto posIni = inTPCPoint->first;
      //Let's record what the initial and final points are.
      trueVtxX = posIni.X();
      trueVtxY = posIni.Y();
      trueVtxZ = posIni.Z();
      trueEndX = posFin.X();
      trueEndY = posFin.Y();
      trueEndZ = posFin.Z();
      
      auto totLength = distance(posFin.X(), posFin.Y(), posFin.Z(),posIni.X(), posIni.Y(), posIni.Z() );
      if (totLength < 0.47 ) continue;

      // Ok, now wer have the first and last point in the TPC and they are reasonably distant
      // Let's use them!
      
      //But first, some stupid checks
      if (verbose)
	{
	  std::cout<<"True Vtx X: "<<posIni.X()<<" Y: "<<posIni.Y()<<" Z: "<<posIni.Z()<<"\n";
	  for (auto const& t : truetraj)
	    {
	      auto pos = t.first;
	      std::cout<<"------------------> "<<pos.X()<<" "<<pos.Y() <<" "<<pos.Z()<<" \n";
	    }
	  std::cout<<"True End X: "<<posFin.X()<<" Y: "<<posFin.Y()<<" Z: "<<posFin.Z()<<"-----> "<< interactionLabel<<"\n\n";
	  std::cout<<"\n\n\n";
	}


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

      //Some other stupid check
      if (verbose)
	{
	  lastPt         = (orderedUniformTrjPts.rbegin())->second;
	  secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
	  lastDist = distance(lastPt.X(),lastPt.Y(),lastPt.Z(),secondtoLastPt.X(),secondtoLastPt.Y(),secondtoLastPt.Z());
      
	  std::cout<<"True End X: "<<posFin.X()<<" Y: "<<posFin.Y()<<" Z: "<<posFin.Z()<<"-----> "<< interactionLabel<<"\n";
	  std::cout<<"True End X: "<<lastPt.X()<<" Y: "<<lastPt.Y()<<" Z: "<<lastPt.Z()<<"-----> "<< interactionLabel<<"\n";
	  std::cout<<"True End X: "<<secondtoLastPt.X()<<" Y: "<<secondtoLastPt.Y()<<" Z: "<<secondtoLastPt.Z()<<"-----> "<< interactionLabel<<"\n";
	  std::cout<<"lastDist "<<lastDist<<"\n\n\n";
	}



      /*

      //Let's try to create a more uniform trajpoint map
      auto old_t = inTPCPoint; 
      double residual   = 0.;
      double overflow   = 0.;
      const double trackPitch = 0.47;


      if (verbose) std::cout << "The distance is: " << std::distance(inTPCPoint,finTPCPoint) << '\n';
      if (std::distance(inTPCPoint,finTPCPoint) == 1)
	{
	
          double distance01       = distance(positionVector0.X(), positionVector0.Y(), positionVector0.Z(), positionVector1.X(), positionVector1.Y(), positionVector1.Z()) ;

	  // If the distance between trjpoints is smaller than the overflow, I don't have enough space for an extra point 
	  if (distance01 > trackPitch) 
	    {
	      // I do have space for at least one extra point, so let's put it there!
	      // Calculate how many extra points I need to put between the new first point and the second TrajPoint
	      int    nPts            = (int) (distance01/trackPitch);
	      for (int iPt = 1; iPt <= nPts; iPt++ )
		{
		  auto newPoint = positionVector0 + iPt*(trackPitch/distance01) * (positionVector1 - positionVector0);
		  orderedUniformTrjPts[newPoint.Z()] = newPoint;
		}
	    }
	}

      for ( auto t = std::next(inTPCPoint); t!= finTPCPoint ; t++, old_t++)
	{ 
	  // Create a line between these two points and divide it in segments of 4.7 mm along the 3-D trajectory
	  // Add the new points between the segments
	  auto positionVector0 = (old_t->first).Vect(); 
	  auto positionVector1 = (t->first).Vect(); 
	  if (verbose) std::cout<<"trjPtOld: "<<positionVector0.X()<<" , "<< positionVector0.Y()<<" , "<< positionVector0.Z()<<"\n";
	  if (verbose) std::cout<<"trjPtNow: "<<positionVector1.X()<<" , "<< positionVector1.Y()<<" , "<< positionVector1.Z()<<"\n";
	  double distance01       = distance(positionVector0.X(), positionVector0.Y(), positionVector0.Z(), positionVector1.X(), positionVector1.Y(), positionVector1.Z()) ;
	  double reducedDistance = distance01 - overflow;
	 
	  // If the distance between trjpoints is smaller than the overflow, I don't have enough space for an extra point 
	  if (distance01 < overflow)   continue;
	     
	  // I do have space for an extra point, so let's put it there!
	  TVector3 firstNewPoint = positionVector0 + (overflow/distance01) * (positionVector1 - positionVector0);
          orderedUniformTrjPts[firstNewPoint.Z()] = firstNewPoint;

	  // Calculate how many extra points I need to put between the new first point and the second TrajPoint
	  int    nPts            = (int) (reducedDistance/trackPitch);	  
	  if (verbose) std::cout<<"distance: "<<distance01<<" reducedDistance "<<reducedDistance<<" overflow "<<overflow<<" residual "<<residual
				<<" nPts "<<nPts<<" "<<reducedDistance/trackPitch<<" "<<reducedDistance - nPts*trackPitch <<"\n";
	  for (int iPt = 1; iPt <= nPts; iPt++ )
	    {
	      auto newPoint = firstNewPoint + iPt*(trackPitch/reducedDistance) * (positionVector1 - firstNewPoint);
	      orderedUniformTrjPts[newPoint.Z()] = newPoint;
	    }
	  
	  residual = reducedDistance - nPts*trackPitch; // new residual, for next point
	  overflow = trackPitch - residual;                   // new overflow, for next point
	}


      // Now we have a map of almost equally spaced trj-points
      // Let's calculate a bunch of things, like the dEdX and the XS!!!
      if (orderedUniformTrjPts.size() < 2 ) continue;
      */
      


      // Calculate the initial kinetic energy
      auto initialMom =     inTPCPoint->second;
      double initialKE = 1000*(TMath::Sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y() + initialMom.Z()*initialMom.Z() + mass*mass ) - mass); 
      hKEAtTPCFF->Fill(initialKE);
      double kineticEnergy = initialKE;

      auto old_it = orderedUniformTrjPts.begin();
      for (auto it = std::next(orderedUniformTrjPts.begin()); it != orderedUniformTrjPts.end(); it++, old_it++ )
	{
	  
	  if (verbose)  std::cout << it->first<<" : " << (it->second).Z() << std::endl ;

	  auto oldPos        = old_it->second;
	  auto currentPos    =     it->second;
	
	  double uniformDist =  (currentPos - oldPos).Mag();
	  h_UniformDistances->Fill(uniformDist);
	  
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

	  hdEVsdX->Fill(currentDepEnergy,(currentPos.Z()-oldPos.Z()) );
	  hdEVsKE->Fill(currentDepEnergy,kineticEnergy);
	  hIncidentKE->Fill(kineticEnergy);

	  /*
	  //Fill the Inelastic and Total Interacting with the last point
	  if (verbose) std::cout<<"InteractionLabel Size"<<interactionLabel.size()<<" "<<interactionLabel<<"  \n";
	  if (it == std::prev(orderedUniformTrjPts.end()) && interactionLabel.find("Inelastic")!= std::string::npos )
	    {
	      hInteractingKE->Fill(kineticEnergy);
	      hInteractingKEInel->Fill(kineticEnergy);
	    }
	  
	  //Fill the Elastic and Total Interacting with the last point
	  if (it == std::prev(orderedUniformTrjPts.end()) && interactionLabel.find("Elastic")!= std::string::npos )
	    {
	      h_DeltaE ->Fill(kineticEnergy -  1000*((finTPCPoint->second).E() - mass) );
	      hInteractingKEElDep->Fill(kineticEnergy);
	      hInteractingKE->Fill(kineticEnergy);
	    }
	  */
	  h_DEUniform->Fill(currentDepEnergy);
	  h_DXUniform->Fill(uniformDist);
	  h_DEDXUniform->Fill(currentDepEnergy/uniformDist);

	  
	}// Loop on OrderedPoints


      //Fill the Inelastic and Total Interacting with the last point
      
      if ( interactionLabel.find("Inelastic")!= std::string::npos )
	{
	  //std::cout<<"Interaction Label: "<<interactionLabel<<"\n";
	  hInteractingKE->Fill(kineticEnergy);
	  hInteractingKEInel->Fill(kineticEnergy);
	}
      
      //Fill the Elastic and Total Interacting with the last point
      if ( interactionLabel.find("Elastic")!= std::string::npos )
	{
	  //std::cout<<"Interaction Label: "<<interactionLabel<<"\n";
	  h_DeltaE ->Fill(kineticEnergy -  1000*((finTPCPoint->second).E() - mass) );
	  hInteractingKEElDep->Fill(kineticEnergy);
	  hInteractingKE->Fill(kineticEnergy);
	  auto MomentumF = finTPCPoint->second;
	  double KEF = 1000*(TMath::Sqrt(MomentumF.X()*MomentumF.X() + MomentumF.Y()*MomentumF.Y() + MomentumF.Z()*MomentumF.Z() + mass*mass ) - mass); //I want this in MeV
	  hInteractingKEEl->Fill(KEF);
	}

      finalKE = kineticEnergy;
      keepInteraction = false;
      if (!interactionLabel.size())
	{
	  throughgoing++;
	  G4Process.push_back("throughgoing");
	}else
	{
	  G4Process.push_back(interactionLabel);
	}
      
    }// MC Particle Loop
  
  

  fTree->Fill();
  
  run     = -999; 
  subrun  = -999; 
  eventN  = -999;
  
  trueVtxX = -999.;
  trueVtxY = -999.;
  trueVtxZ = -999.;
  trueEndX = -999.;
  trueEndY = -999.;
  trueEndZ = -999.;
  trueEndZ = -999.;
  finalKE  = -999.;
  G4Process.clear(); 
}

void lariat::TrueXSAngleCut::endJob()
{
  std::cout<<"-------------------------------------------"<<std::endl;
  std::cout<<"True Events pre-TPC .............. "<<evtsPreTPC<<std::endl;
  std::cout<<"True Events pre-TPC .............. "<<evtsInTheMiddle<<std::endl;
  std::cout<<"True Events post-TPC ............. "<<evtsPostTPC<<std::endl;
  std::cout<<"True Throughgoing    ............. "<<throughgoing<<std::endl;
  std::cout<<"True interactingInTPC ............ "<<interactingInTPC<<std::endl;
  std::cout<<"-------------------------------------------"<<std::endl;

  float rho = 1396; //kg/m^3
  float molar_mass = 39.95; //g/mol
  float g_per_kg = 1000; 
  float avogadro = 6.022e+23; //number/mol
  float number_density = rho*g_per_kg/molar_mass*avogadro;
  float slab_width = 0.0047;//in m


  // Calculate the Cross Section
  // ###################################################################
  // #### Looping over the exiting bins to extract the cross-section ###
  // ###################################################################
  for( int iBin = 1; iBin <= hInteractingKE->GetNbinsX(); ++iBin )
    {

      // ### If an incident bin is equal to zero then skip that bin ###
      if( hIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
   
      // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
      float TempCrossSection = (hInteractingKE->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);
   
      float elCrossSection   = ((hInteractingKEEl->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width) ) * (1/1e-28);
      float inelCrossSection = ((hInteractingKEInel->GetBinContent(iBin)/hIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width) ) * (1/1e-28);
      // ### Covert this into Barns ###
      float crossSection = TempCrossSection * (1/1e-28); 
	
      // ### Putting the value on the plot
      hCrossSection    ->SetBinContent(iBin,crossSection);
      hCrossSectionEl  ->SetBinContent(iBin,elCrossSection);
      hCrossSectionInel->SetBinContent(iBin,inelCrossSection);

      // ###########################################################
      // ### Calculating the error on the numerator of the ratio ###
      // ###########################################################

      float denomError = pow(hIncidentKE->GetBinContent(iBin),0.5);
      float denom = hIncidentKE->GetBinContent(iBin);
      if(denom == 0) continue; 
      float term2 = denomError/denom;
      
      float numError = pow(hInteractingKE->GetBinContent(iBin),0.5);
      float num = hInteractingKE->GetBinContent(iBin);

      float numErrorEl = pow(hInteractingKEEl->GetBinContent(iBin),0.5);
      float numEl = hInteractingKEEl->GetBinContent(iBin);

      float numErrorInel = pow(hInteractingKEInel->GetBinContent(iBin),0.5);
      float numInel = hInteractingKEInel->GetBinContent(iBin);

      // ### Putting in a protection against dividing by zero ###   
      if(num != 0){
	float term1 = numError/num;
	float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
	hCrossSection->SetBinError(iBin,totalError);
      }
      
      if(numEl != 0){
	float term1 = numErrorEl/numEl;
	float totalError = (elCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width)  *(1e26);
	hCrossSectionEl->SetBinError(iBin,totalError);
      }
      
      if(numInel != 0){
	float term1 = numErrorInel/numInel;
	float totalError = (inelCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width)  *(1e26);
	hCrossSectionInel->SetBinError(iBin,totalError);
      }
    }//<---End iBin Loop


}

double lariat::TrueXSAngleCut::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return d;
}

double lariat::TrueXSAngleCut::calculateAngle(double px1, double py1, double pz1, double px2, double py2, double pz2)
{
  double magP1   = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1 );
  double magP2   = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2 );
  double dotProd = px1*px2 + py1*py2 + pz2*pz2 ;
  double angle = -999999.;
  if (magP1*magP1 > 0.00000001) // this is what I call zero... not great, but will go for now
    {
      double ratio = dotProd/( magP1 * magP2);
      angle = TMath::ACos(ratio); 
    }

  return angle;
}


void lariat::TrueXSAngleCut::beginJob()
{
  
  art::ServiceHandle<art::TFileService> tfs;
  h_AngleEl     = tfs->make<TH1D>("h_AngleEl","h_AngleEl; Angle Elastic Scattering [rad]",1260, -6.30,6.30);   
  h_AngleInel   = tfs->make<TH1D>("h_AngleInel","h_AngleInel; Angle Inelastic Scattering [rad]",1260, -6.30,6.30);   
  h_AngleEl_AfterCut     = tfs->make<TH1D>("h_AngleEl_AfterCut","h_Angle; Angle Elastic Scattering [rad]",1260, -6.30,6.30);   
  h_AngleInel_AfterCut   = tfs->make<TH1D>("h_AngleInel_AfterCut","h_Angle; Angle Elastic Scattering [rad]",1260, -6.30,6.30);   


  h_DE   = tfs->make<TH1D>("h_DE","h_DE; Energy Deposited [MeV]",200, 0,100);   
  h_DX   = tfs->make<TH1D>("h_DX","h_DX; Distance between points  [cm]",400, 0,20);   
  h_DEDX = tfs->make<TH1D>("h_DEDEX","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   

  h_DEUniform   = tfs->make<TH1D>("h_DEUniform","h_DE; Energy Deposited [MeV]",200, 0,100);   
  h_DXUniform   = tfs->make<TH1D>("h_DXUniform","h_DX; Distance between points  [cm]",400, 0,20);   
  h_DEDXUniform = tfs->make<TH1D>("h_DEDEXUniform","h_DEDX; dE/dX [MeV/cm]",500, 0,50);   

  h_DeltaE = tfs->make<TH1D>("h_DeltaE","h_DeltaE; dEDep - TrjDE [MeV/cm]",500, -1000,1000);   
  h_SimIDEDist= tfs->make<TH1D>("h_SimIDEDist","h_SimIDEDist; h_SimIDEDist [cm]",1000, 0,10);   

  
  h_UniformDistances = tfs->make<TH1D>("h_UniformDistances","h_UniformDistances; Distance between uniform points  [cm]",500, 0,5);   

  hInitialPz     = tfs->make<TH1D>("hInitialPz"    , "Initial Pz [MeV/c]"    , 42, -100, 2000);
  hInitialKE     = tfs->make<TH1D>("hInitialKE"    , "Initial Kinetic Energy [MeV]"    , 42, -100, 2000);
  hKEAtTPCFF     = tfs->make<TH1D>("hKEAtTPCFF"    , "Kinetic Energy @ TPC FF [MeV]"   , 42, -100, 2000);


  hIncidentKE           = tfs->make<TH1D>("hIncidentKE"   , "Incident Kinetic Energy [MeV]"   , 42, -100, 2000); 
  hInteractingKE        = tfs->make<TH1D>("hInteractingKE", "Interacting Kinetic Energy [MeV]", 42, -100, 2000);
  hInteractingKEEl      = tfs->make<TH1D>("hInteractingKEEl", "Elastic Interacting Kinetic Energy [MeV]", 42, -100, 2000); 
  hInteractingKEElDep   = tfs->make<TH1D>("hInteractingKEElDep", "Dep Elastic Interacting Kinetic Energy [MeV]", 42, -100, 2000); 
  hInteractingKEInel    = tfs->make<TH1D>("hInteractingKEInel", "Inelastic Interacting Kinetic Energy [MeV]", 42, -100, 2000);


  hCrossSection     = tfs->make<TH1D>("hCrossSection"     , "Cross-Section [barn]"             , 42, -100, 2000);
  hCrossSectionEl   = tfs->make<TH1D>("hCrossSectionEl"   , "Elastic Cross-Section [barn]"     , 42, -100, 2000);
  hCrossSectionInel = tfs->make<TH1D>("hCrossSectionInel" , "Inelastic Cross-Section [barn]"   , 42, -100, 2000);


  hXZ     = tfs->make<TH2D>("hXZ"     , "hXZ"    , 110, -100, 10, 200, -100, 100);  
  hYZ     = tfs->make<TH2D>("hYZ"     , "hYZ"    , 110, -100, 10, 200, -100, 100); 
  hXZPre  = tfs->make<TH2D>("hXZPre"  , "hXZPre" , 110, -100, 10, 200, -100, 100); 
  hYZPre  = tfs->make<TH2D>("hYZPre"  , "hYZPre" , 110, -100, 10, 200, -100, 100); 


  hdEVsdX  = tfs->make<TH2D>("hdEVsdX"  , "hdEVsdX" , 504, -1, 50, 1100, -10, 100); 
  hdEVsKE  = tfs->make<TH2D>("hdEVsKE"  , "hdEVsKE" , 504, -1, 50, 220,  -10, 1000); 

  
  fTree = tfs->make<TTree>("effTree","analysis tree");
  fTree->Branch("run"      ,&run      ,"run/I");
  fTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  fTree->Branch("eventN"   ,&eventN   ,"eventN/I");
  

  fTree->Branch("trueVtxX" ,&trueVtxX ,"trueVtxX/D");
  fTree->Branch("trueVtxY" ,&trueVtxY ,"trueVtxY/D");
  fTree->Branch("trueVtxZ" ,&trueVtxZ ,"trueVtxZ/D");
  fTree->Branch("trueEndX" ,&trueEndX ,"trueEndX/D");
  fTree->Branch("trueEndY" ,&trueEndY ,"trueEndY/D");
  fTree->Branch("trueEndZ" ,&trueEndZ ,"trueEndZ/D");
  fTree->Branch("finalKE"  ,&finalKE  ,"finalKE/D" );
  fTree->Branch("G4Process",&G4Process);
}



DEFINE_ART_MODULE(lariat::TrueXSAngleCut)
