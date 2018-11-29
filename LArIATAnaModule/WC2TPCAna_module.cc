////////////////////////////////////////////////////////////////////////
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
  class WC2TPCAna;
}

class lariat::WC2TPCAna : public art::EDAnalyzer 
{
public:
  explicit WC2TPCAna(fhicl::ParameterSet const & p);
  virtual ~WC2TPCAna();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const & p);

private:
  double minX =  0.0;
  double maxX = 47.0;
  double minY =-20.0;
  double maxY = 20.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 90.0;
  
  int myCount1 = 0;
  int myCount2 = 0;
  
  // Declare member data here.
  
  std::string 	fTrackModuleLabel;
  std::string 	fWCTrackLabel; 		// The name of the producer that made tracks through the MWPCs

  art::ProductID WCProductId;

  // --- Histos for checking filter ---
  TTree *fTree;
  int    runN           = 0;
  int    subrunN        = 0;
  int    evtN           = 0;
  // Truth Variables
  int    primaryInTPC   = 0;
  double trueCTPCLength = -999.;
  // Variables for reco correctly matched
  double recoCTPCLength = -999.;
  double minDisTrueReco =  999.;
  double recoZC         =  999.;
  int    correctTrkNum  = -99;
  // Variables for everybody!
  double recoZ          =  999.;
  int    trackNumber    = -99;
  double DeltaX         = -999.;
  double DeltaY         = -999.;
  double alphaM         = -999.;

};


lariat::WC2TPCAna::WC2TPCAna(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::WC2TPCAna::~WC2TPCAna()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::WC2TPCAna::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel		= pset.get< std::string >("TrackModuleLabel");
  fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel","wctrack");
}

void lariat::WC2TPCAna::analyze(art::Event const & evt)
{

  myCount1 ++;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  runN           = evt.run();
  subrunN        = evt.subRun();
  evtN           = evt.event();
  DeltaX        = -999.;
  DeltaY        = -999.;
  alphaM        = -999.;
  trackNumber   = -99;
  correctTrkNum = -99;
  primaryInTPC  = 0;
  trueCTPCLength = 0.;


  
  // #################################################
  // ### (0) Is the primary in TPC? What's its length?
  // #################################################
  double trueVtxX = -999.;
  double trueVtxY = -999.;
  double trueVtxZ = -999.;

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

      //--------------------------------------------------------
      // Identify the first trajectory point in TPC
      auto inTPCPoint    = truetraj.begin();
      auto lastTPCPoint  = truetraj.begin();
      // Loop From First TrajPoint --> First Point in TPC
      bool firstPointTrue = true;
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++)
	{
	  auto pos = t->first;
	  if      (pos.Z() < minZ || pos.Z() > maxZ ) continue;
	  else if (pos.X() < minX || pos.X() > maxX ) continue;
	  else if (pos.Y() < minY || pos.Y() > maxY ) continue;
	  else {
	    primaryInTPC = 1;
	    lastTPCPoint = t;
	    if (firstPointTrue) {inTPCPoint = t; firstPointTrue = false;}
	  }
	}// End search for first point in TPC

      trueVtxX = (inTPCPoint->first).X();
      trueVtxY = (inTPCPoint->first).Y();
      trueVtxZ = (inTPCPoint->first).Z();

      for ( auto iter = inTPCPoint; iter!= lastTPCPoint; iter++)
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
	  trueCTPCLength += distanceBtwPoint.Mag();
	}

      if (trueCTPCLength < 0.47 ) primaryInTPC = 0; // The trajectory of this particle is impossible to reconstruct
    } // end geant particle list
  //------------------------> end of true stuff


    
  // #####################################
  // ### Getting the Track Information ###
  // #####################################
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
  
  // === Filling the tracklist from the tracklistHandle ===
  if (!evt.getByLabel(fTrackModuleLabel,trackListHandle))  return;
  art::fill_ptr_vector(tracklist, trackListHandle);


  
  // #########################################################
  // ### Grabbing associations for use later in the Filter ###
  // #########################################################  
  // === Association between SpacePoints and Tracks ===
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);


  
  // ###################################################
  // ### Getting the Wire Chamber Track  Information ###
  // ###################################################
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
      
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;
  art::fill_ptr_vector(wctrack, wctrackHandle);
  
  
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // 				Grabbing the wire chamber tracks
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  float wctrk_XFace[10] = {0.};
  float wctrk_YFace[10] = {0.};
  float wctrk_theta[10] = {0.};
  float wctrk_phi[10] = {0.};

  int   nwctrk = 0;
  
  // ##################################
  // ### Looping over the WC Tracks ###
  // ##################################
  for(size_t wcCount = 0; wcCount < wctrack.size(); wcCount++)
    {
      wctrk_XFace[nwctrk] = wctrack[nwctrk]->XYFace(0) ;//<---
      wctrk_YFace[nwctrk] = wctrack[nwctrk]->XYFace(1) ;//<---is this cm?
      wctrk_theta[nwctrk] = wctrack[nwctrk]->Theta();
      wctrk_phi[nwctrk]   = wctrack[nwctrk]->Phi();
      nwctrk++;
    }//<----End wcCount loop
  
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // 			   Grabbing the upstream most trajectory points
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  //double larStart[3];
  //double larEnd[3];
  //std::vector<double> trackStart;
  //std::vector<double> trackEnd;
  
  
  int ntrks = 0;
  
  // ###########################################################
  // ### Variables for storing the lowest Z trajectory point ###
  // ###    that is inside the active volume of the TPC      ###
  // ###########################################################
  float UpStreamTrjPointZ[100] = {999.};
  float UpStreamTrjPointY[100] = {999.};
  float UpStreamTrjPointX[100] = {999.}; //<---Initialized to a stupid high number on purpose
  
  float trackLength[100] = {999.}; //<---Initialized to a stupid high number on purpose
  
  float UpStreamTrjPointPHatY[100] = {999.};
  float UpStreamTrjPointPHatX[100] = {999.}; //<---Initialized to a stupid high number on purpose
  
  
  
  float FirstTrjPtZ = 999;
  float FirstTrjPtY = 999;
  float FirstTrjPtX = 999;
  
  float pHatX = 999;
  float pHatY = 999;
  
  float tpcTheta[100]= {0.};
  
  // ### Storing the trajectory points in a similar way to PionXS ###
  
  TVector3 z_hat(0,0,1);
  
  // #######################################
  // ### Looping over all the tpcTracks ###
  // ######################################
  for(size_t i=0; i<tracklist.size();++i)
    {
      // ### Clearing the vectors for each track ###
      //trackStart.clear();
      //trackEnd.clear();
      
      TVector3 p_hat_0;
      
      // ### Making temp variables to find the most upstream ###
      // ###         of all the trajectory points            ###
      FirstTrjPtZ = 999;
      FirstTrjPtY = 999;
      FirstTrjPtX = 999;
      
      pHatX = 999;
      pHatY = 999;
      trackLength[i] = tracklist[i]->Length();

      // ##############################################
      // ### Looping over all the trajectory points ###
      // ##############################################
      for( size_t iTrajPt = 0; iTrajPt < tracklist[i]->NumberTrajectoryPoints() ; ++iTrajPt )
	{
	  // ### Recording this point if it is the most upstream point ###
	  if(tracklist[i]->LocationAtPoint(iTrajPt).Z() < FirstTrjPtZ && 
	     tracklist[i]->LocationAtPoint(iTrajPt).Z() > 0.0 && tracklist[i]->LocationAtPoint(iTrajPt).Y() > -20.0 &&
	     tracklist[i]->LocationAtPoint(iTrajPt).Y() < 20.0 && tracklist[i]->LocationAtPoint(iTrajPt).X() > 0.0 &&
	     tracklist[i]->LocationAtPoint(iTrajPt).X() < 43.5)
	    {
	      
	      p_hat_0 = tracklist[i]->DirectionAtPoint(iTrajPt);
	      //Strange directionality convention - I'm reversing the direction vector
	      //if it's pointing in the negative Z direction
	      if( p_hat_0.Z() < 0 )
		{
		  p_hat_0.SetX(p_hat_0.X()*-1);
		  p_hat_0.SetY(p_hat_0.Y()*-1);
		  p_hat_0.SetZ(p_hat_0.Z()*-1);
		}
	      
	      // ###    Storing the temp variables of what is now      ###
	      // ### considered the most upstream point for this track ###
	      FirstTrjPtZ = tracklist[i]->LocationAtPoint(iTrajPt).Z();
	      FirstTrjPtY = tracklist[i]->LocationAtPoint(iTrajPt).Y();
	      FirstTrjPtX = tracklist[i]->LocationAtPoint(iTrajPt).X();
	      
	      pHatY = p_hat_0.Y();
	      pHatX = p_hat_0.X();
	    }//<---End only storing points if they are the lowest Z point
	}//<---End iTrajPt loop
      // ### Calculating the Theta for the TPC Track ###
      tpcTheta[i]=acos(z_hat.Dot(p_hat_0)/p_hat_0.Mag());
      
      // ###################################################
      // ### Saving for looping later the upstream point ###
      // ###       for every track in the list           ###
      // ###################################################
      //if(FirstTrjPtZ < 999)
      //{
      UpStreamTrjPointZ[i] = FirstTrjPtZ;  
      UpStreamTrjPointY[i] = FirstTrjPtY;
      UpStreamTrjPointX[i] = FirstTrjPtX;
      
      UpStreamTrjPointPHatY[i] = pHatY;
      UpStreamTrjPointPHatX[i] = pHatX;
      
      
      ntrks++;
      //}// <----End saving FirstPoint   
    }//<---End i loop
 
  
  // ############################################
  // ### Loop over all the eligible tpcTracks ###
  // ############################################
  
  //Fill the TTree even if there's no reco track in TPC
  if (ntrks == 0)
    {
      recoCTPCLength = -999.;
      minDisTrueReco =  999.;
      recoZC         =  999.;
      correctTrkNum  = -99;
      recoZ          =  999.;
      trackNumber    = -99;
      DeltaX         = -999.;
      DeltaY         = -999.;
      alphaM         = -999.;
      
      fTree->Fill(); // Fill for every possible WC-TPC couple	  
    }


  for(int iitrk = 0; iitrk < ntrks; iitrk++)
    {
      
      // ##################################
      // ### (1) Which one is the correct track?
      // ##################################
      if (primaryInTPC)
	{
	  //	  std::cout<<"LENGTH: "<< trackLength[iitrk] <<"\n";
	  if (trackLength[iitrk] > 2. ) {
	    //std::cout<<"primary: "<<primaryInTPC<<" "<<UpStreamTrjPointX[iitrk]<<" "<<UpStreamTrjPointY[iitrk]<<" "<<UpStreamTrjPointZ[iitrk]<<" "<<iitrk<<".................\n";
	    double dist = TMath::Sqrt( (trueVtxX - UpStreamTrjPointX[iitrk]) * (trueVtxX - UpStreamTrjPointX[iitrk]) +
				       (trueVtxY - UpStreamTrjPointY[iitrk]) * (trueVtxY - UpStreamTrjPointY[iitrk]) +
				       (trueVtxZ - UpStreamTrjPointZ[iitrk]) * (trueVtxZ - UpStreamTrjPointZ[iitrk]));
	    //std::cout<<"primary: "<<primaryInTPC<<" "<<correctTrkNum<<" "<<minDisTrueReco<<" "<<recoCTPCLength<<" "<<dist<<" "<< minDisTrueReco <<"................. 1\n";
	    if ( dist < minDisTrueReco) {
	      minDisTrueReco = dist; 
	      correctTrkNum  = iitrk; 
	      recoCTPCLength = trackLength[iitrk];
	      recoZC         = UpStreamTrjPointZ[iitrk];
	      //std::cout<<"primary: "<<primaryInTPC<<" "<<correctTrkNum<<" "<<minDisTrueReco<<" "<<recoCTPCLength<<" "<<dist<<"................. 2\n";
	    }
	  }
	}else 
	{
	  // if primary is not TPC there is no correct track
	  recoCTPCLength = -999.;
	  minDisTrueReco =  999.;
	  recoZC         =  999.;
	  correctTrkNum  = -99;
	}
    } // -----------------------> Correct Matched establish, if there's one :)


  // Let's see how we're doing not knowing the truth information....
  for(int aa = 0; aa < ntrks; aa++)
    {
      //p_hat_0 = tracklist[aa]->DirectionAtPoint(aa);
      float tpc_Theta=tpcTheta[aa];
      //std::cout<<"&&&& TPC Theta "<<tpcTheta[aa]<<" "<<std::endl;
      // ###########################################
      // ### Loop over all the eligible WCTracks ###
      // ###########################################
      for(int bb = 0; bb < nwctrk; bb++)
	{
	  // ### Grabbing the WCTrack Theta ###
	  float wcTheta = wctrk_theta[bb];
	  
	  // ### Grabbing the WCTRack Phi ###
	  float wcPhi = wctrk_phi[bb];
	  
	  // ### Using same convention as WCTrack to calculate phi ###
	  float phi = 0;
	  //Calculating phi (degeneracy elimination for the atan function)
	  //----------------------------------------------------------------------------------------------
	  if( UpStreamTrjPointPHatY[aa] > 0 && UpStreamTrjPointPHatX[aa] > 0 ){ 
	    
	    phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa]); 
	  }
	  else if( UpStreamTrjPointPHatY[aa] > 0 && UpStreamTrjPointPHatX[aa] < 0 ){ 
	    
	    phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+3.141592654; }
	  else if( UpStreamTrjPointPHatY[aa] < 0 && UpStreamTrjPointPHatX[aa] < 0 ){ 
	    
	    phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+3.141592654; }
	  else if( UpStreamTrjPointPHatY[aa] < 0 && UpStreamTrjPointPHatX[aa] > 0 ){ 
	    
	    phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+6.28318; }
	  else if( UpStreamTrjPointPHatY[aa] == 0 && UpStreamTrjPointPHatX[aa] == 0 ){ 
	    
	    phi = 0; }//defined by convention
	  else if( UpStreamTrjPointPHatY[aa] == 0 )
	    {
	      if( UpStreamTrjPointPHatX[aa] > 0 ){ 
		phi = 0; }
	      else{
		phi = 3.141592654;
	      }
	    }
	  else if( UpStreamTrjPointPHatX[aa] == 0 )
	    {
	      if( UpStreamTrjPointPHatY[aa] > 0 ){ phi = 3.141592654/2; }
	      else{ phi = 3.141592654*3/2; }
	    }
	  //----------------------------------------------------------------------------------------------
	  // ### Using TPC Phi ###
	  float tpcPhi = phi; 
	  
	  // #########################################################
	  // ### Define the unit vectors for the WC and TPC Tracks ###
	  // #########################################################
	  TVector3 theUnitVector_WCTrack;
	  TVector3 theUnitVector_TPCTrack;
	  
	  // === WCTrack Unit Vector ===
	  theUnitVector_WCTrack.SetX(sin(wcTheta)*cos(wcPhi));
	  theUnitVector_WCTrack.SetY(sin(wcTheta)*sin(wcPhi));
	  theUnitVector_WCTrack.SetZ(cos(wcTheta));
	  
	  // ### TPC Track Unit Vector ===
	  theUnitVector_TPCTrack.SetX(sin(tpc_Theta)*cos(tpcPhi));
	  theUnitVector_TPCTrack.SetY(sin(tpc_Theta)*sin(tpcPhi));
	  theUnitVector_TPCTrack.SetZ(cos(tpc_Theta));
	  
	  // ###########################################################
	  // ### Calculating the angle between WCTrack and TPC Track ###
	  // ###########################################################
	  float alpha = ( acos(theUnitVector_WCTrack.Dot(theUnitVector_TPCTrack)) )* (180.0/3.141592654);
	  
	  recoZ  = UpStreamTrjPointZ[aa]; 
	  DeltaX = UpStreamTrjPointX[aa] - (wctrk_XFace[bb]);
	  DeltaY = UpStreamTrjPointY[aa] - (wctrk_YFace[bb]);
	  alphaM = alpha;
	  trackNumber = aa;
	  fTree->Fill(); // Fill for every possible WC-TPC couple	  
	}//<---End bb loop
    }//<---End aa loop


  // Reset all variables to original values.
  runN           = 0;
  subrunN        = 0;
  evtN           = 0;
  //Truth Variables
  primaryInTPC   = 0;
  trueCTPCLength = -999.;
  //Variables for reco correctly matched
  recoCTPCLength = -999.;
  minDisTrueReco =  999.;
  recoZC         =  999.;
  correctTrkNum  = -99;
  //Variables for everybody!
  recoZ          =  999.;
  trackNumber    = -99;
  DeltaX         = -999.;
  DeltaY         = -999.;
  alphaM         = -999.;

}

void lariat::WC2TPCAna::endJob()
{
  std::cout<<"Counts: "<<myCount1<<" "<<  myCount2<<"\n";
}

void lariat::WC2TPCAna::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("everyMatchTree","analysis tree");
  fTree->Branch("runN"     ,&runN     ,"runN/I");
  fTree->Branch("subrunN"  ,&subrunN  ,"subrunN/I");
  fTree->Branch("evtN"     ,&evtN     ,"evtN/I");
  fTree->Branch("DeltaX"              ,&DeltaX               ,"DeltaX/D"           );
  fTree->Branch("DeltaY"              ,&DeltaY               ,"DeltaY/D"           );
  fTree->Branch("alphaM"              ,&alphaM               ,"alphaM/D"           );
  fTree->Branch("trackNumber"         ,&trackNumber          ,"trackNumber/I"      );
  fTree->Branch("correctTrkNum"       ,&correctTrkNum        ,"correctTrkNum/I"    );
  fTree->Branch("primaryInTPC"        ,&primaryInTPC         ,"primaryInTPC/I"     );
  fTree->Branch("trueTPCLength"       ,&trueCTPCLength       ,"trueTPCLength/D"    );
  fTree->Branch("recoCTPCLength"      ,&recoCTPCLength       ,"recoCTPCLength/D"   );
  fTree->Branch("recoZ"               ,&recoZ                ,"recoZ/D"            );
  fTree->Branch("recoZC"              ,&recoZC               ,"recoZC/D"           );
  fTree->Branch("minDisTrueReco"      ,&minDisTrueReco       ,"minDisTrueReco/D"   );

}

DEFINE_ART_MODULE(lariat::WC2TPCAna)

//  LocalWords:  ntrks
