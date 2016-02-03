////////////////////////////////////////////////////////////////////////
// Class:       WC2TPCTrackMatch
// Module Type: producer
// File:        WC2TPCTrackMatch_module.cc
//
// Generated at Wed Jan 27 15:06:47 2016 by Elena Gramellini using artmod
// from cetpkgsupport v1_08_06.
// Sort of detailed to do list...
// [ x ] Take all the wcTrack
// [ x ] Take all the tpcTrack
// [ x ] Store the ones that pass Jonathan's cuts as 
//       std::pair<int wcTrackKey, int tpcTrackKey> matchCandidate
// [ x ] check if a wcTrack or tpcTrack is used more than once
// [ x ] define a fig of merit to understand which match is the best
// [ x ] subtrack candidates from previous std::pair<int wcTrackKey, int tpcTrackKey> matchCandidate
// [   ] create association
// [ x ] store  association
////////////////////////////////////////////////////////////////////////


// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "LArIATDataProducts/WCTrack.h"

// ########################
// ### LArsoft Includes ###
// ########################
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"

// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

// #####################
// ### ROOT Includes ###
// #####################
#include <TH1F.h>
#include <TH2F.h>

class WC2TPCTrackMatch;

class WC2TPCTrackMatch : public art::EDProducer {
public:
  explicit WC2TPCTrackMatch(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WC2TPCTrackMatch(WC2TPCTrackMatch const &) = delete;
  WC2TPCTrackMatch(WC2TPCTrackMatch &&) = delete;
  WC2TPCTrackMatch & operator = (WC2TPCTrackMatch const &) = delete;
  WC2TPCTrackMatch & operator = (WC2TPCTrackMatch &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:
  
  // Declare member data here.
  
  std::string 	fTrackModuleLabel;
  std::string 	fWCTrackLabel; 		// The name of the producer that made tracks through the MWPCs
  double 	falpha;			// The angle between the WCTrack and TPC Track
  double	fDeltaXLow;		// The lower cut of the Delta X variable
  double 	fDeltaXHigh;		// The upper cut of the Delta X variable
  double	fDeltaYLow;		// The lower cut of the Delta Y variable
  double 	fDeltaYHigh;		// The upper cut of the Delta Y variable
  double	fMaxZPos;		// Maximum Z position of the TPC track we will consider
  int		fMaxMatchedTracks;	// The maximum number of tracks allowed to be matched per event
  
  // --- Histos for checking filter ---
  TH1F*		hDeltaX;
  TH1F*		hDeltaY;
  TH1F*		hAlpha;
  TH1F*		hnMatch;
  
};


// ---------------------- Parameter Setting ---------------------
WC2TPCTrackMatch::WC2TPCTrackMatch(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
  produces < art::Assns<ldp::WCTrack,recob::Track> >();
}


// ---------------------- Event Loop ---------------------------
void WC2TPCTrackMatch::produce(art::Event & evt)
{
  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
  
  // #####################################
  // ### Getting the Track Information ###
  // #####################################
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
  
  // === Filling the tracklist from the tracklistHandle ===
  if (!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return;
  art::fill_ptr_vector(tracklist, trackListHandle);

  // ###################################################
  // ### Getting the Wire Chamber Track  Information ###
  // ###################################################
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
   
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;
  art::fill_ptr_vector(wctrack, wctrackHandle);

  //##########################################################################
  // Make a std::unique_ptr<> for the thing you want to put into the event ###
  // because that handles the memory management for you                    ###
  //##########################################################################
  std::unique_ptr< art::Assns<ldp::WCTrack , recob::Track> > wcTpcTrackAssn(new art::Assns<ldp::WCTrack , recob::Track>);
  
  // #########################################################
  // ### Grabbing associations for use later in the Filter ###
  // #########################################################  
  // === Association between SpacePoints and Tracks ===
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  
  //std::cout<<"========================================="<<std::endl;
  //std::cout<<"Run = "<<evt.run()<<", SubRun = "<<evt.subRun()<<", Evt = "<<evt.id().event()<<std::endl;
  //std::cout<<"========================================="<<std::endl;
  
  
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // 				Grabbing the wire chamber tracks
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  float wctrk_XFace[10] = {0.};
  float wctrk_YFace[10] = {0.};
  float wctrk_theta[10] = {0.};
  float wctrk_phi[10] = {0.};
  int   wcTrackKey[10] = {0};
  int   nwctrk = 0;
  
  // ##################################
  // ### Looping over the WC Tracks ###
  // ##################################
  for(size_t wcCount = 0; wcCount < wctrack.size(); wcCount++)
    {
      wctrk_XFace[nwctrk] = wctrack[nwctrk]->XYFace(0) * 0.1;//<---Note: *0.1 to convert to cm
      wctrk_YFace[nwctrk] = wctrack[nwctrk]->XYFace(1) * 0.1;//<---Note: *0.1 to convert to cm
      wctrk_theta[nwctrk] = wctrack[nwctrk]->Theta();
      wctrk_phi[nwctrk]   = wctrack[nwctrk]->Phi();
      wcTrackKey[nwctrk]  = wctrack[nwctrk].key();
      nwctrk++;
    }//<----End wcCount loop
  
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  // 			   Grabbing the upstream most trajectory points
  // ---------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  
  
  int ntrks = 0;
  
  // ###########################################################
  // ### Variables for storing the lowest Z trajectory point ###
  // ###    that is inside the active volume of the TPC      ###
  // ###########################################################
  float UpStreamTrjPointZ[100] = {999.};
  float UpStreamTrjPointY[100] = {999.};
  float UpStreamTrjPointX[100] = {999.}; //<---Initialized to a stupid high number on purpose
  
  
  float UpStreamTrjPointPHatY[100] = {999.};
  float UpStreamTrjPointPHatX[100] = {999.}; //<---Initialized to a stupid high number on purpose
  
  int tpcTrackKey[100] = {0};


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
      trackStart.clear();
      trackEnd.clear();
      
      TVector3 p_hat_0;
      
      // ### Setting the track information into memory ###
      memset(larStart, 0, 3);
      memset(larEnd, 0, 3);
      tracklist[i]->Extent(trackStart,trackEnd); 
      tracklist[i]->Direction(larStart,larEnd);
      
      // ### Setting the number of Trajectory points for this track ###
      //nTrajPoints[i] = tracklist[ntrks]->NumberTrajectoryPoints();
      // ##############################################
      // ### Looping over all the trajectory points ###
      // ##############################################
      for( size_t iTrajPt = 0; iTrajPt < tracklist[i]->NumberTrajectoryPoints() ; ++iTrajPt )
	{
	  // ### Making temp variables to find the most upstream ###
	  // ###         of all the trajectory points            ###
	  FirstTrjPtZ = 999;
	  FirstTrjPtY = 999;
	  FirstTrjPtX = 999;
	  
	  pHatX = 999;
	  pHatY = 999;
	  
	  // ### Recording this point if it is the most upstream point ###
	  if(tracklist[i]->LocationAtPoint(iTrajPt).Z() < FirstTrjPtZ && 
	     tracklist[i]->LocationAtPoint(iTrajPt).Z() > 0.0 && tracklist[i]->LocationAtPoint(iTrajPt).Y() > -20.0 &&
	     tracklist[i]->LocationAtPoint(iTrajPt).Y() < 20.0 && tracklist[i]->LocationAtPoint(iTrajPt).X() > 0.0 &&
	     tracklist[i]->LocationAtPoint(iTrajPt).X() < 43.5)
	    {
	      
	      p_hat_0 = tracklist[i]->DirectionAtPoint(iTrajPt);
	      //Strange directionality convention - I'm reversing the direction vector
	      //if it's pointing in the negative X direction
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

      tpcTrackKey[i] = tracklist[i].key();
      
      ntrks++;
      //}// <----End saving FirstPoint   
    }//<---End i loop
  
  
  //---------------------------------------------------------------------------------------------------------------------------------
  // 						DOING THE MATCHING OF THE WCTRACK AND TPC TRACKS
  //---------------------------------------------------------------------------------------------------------------------------------
  
  // ### Initialize the counter for the number of matches ###
  int nWC_TPC_TrackMatch = 0;
  std::pair <int,int>  matchCandidate;
  std::map<float,std::pair <int,int>> mapCandidates; // the map takes the figure of merit as key and the couple as value 

  // ###################################################
  // ### Vectors for angles between TPC and WC Track ###
  // ###################################################
  
  
  
  // ############################################
  // ### Loop over all the eligible tpcTracks ###
  // ############################################
  for(int aa = 0; aa < ntrks; aa++)
    {
      float DeltaX_WC_TPC_Track = 999;
      float DeltaY_WC_TPC_Track = 999;
      //p_hat_0 = tracklist[aa]->DirectionAtPoint(aa);
      float tpc_Theta=tpcTheta[aa];
      
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
	  if( UpStreamTrjPointPHatY[aa] > 0 && UpStreamTrjPointPHatX[aa] > 0 ){ phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa]); }
	  else if( UpStreamTrjPointPHatY[aa] > 0 && UpStreamTrjPointPHatX[aa] < 0 ){ phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+3.141592654; }
	  else if( UpStreamTrjPointPHatY[aa] < 0 && UpStreamTrjPointPHatX[aa] < 0 ){ phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+3.141592654; }
	  else if( UpStreamTrjPointPHatY[aa] < 0 && UpStreamTrjPointPHatX[aa] > 0 ){ phi = atan(UpStreamTrjPointPHatY[aa]/UpStreamTrjPointPHatX[aa])+6.28318; }
	  else if( UpStreamTrjPointPHatY[aa] == 0 && UpStreamTrjPointPHatX[aa] == 0 ){ phi = 0; }//defined by convention
	  else if( UpStreamTrjPointPHatY[aa] == 0 )
	    {
	      if( UpStreamTrjPointPHatX[aa] > 0 ){ phi = 0; }
	      
	      else{ phi = 3.141592654; }
	      
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
	  
	  DeltaX_WC_TPC_Track = UpStreamTrjPointX[aa] - (wctrk_XFace[bb]);
	  DeltaY_WC_TPC_Track = UpStreamTrjPointY[aa] - (wctrk_YFace[bb]);
	  
	  
	  // ### Only counting the track if it passes the alpha, DeltaX and Delta Y ###
	  // ###       and is far enough upstream in the TPC in Z position          ###
	  if(alpha < falpha && 
	     DeltaX_WC_TPC_Track > fDeltaXLow && 
	     DeltaY_WC_TPC_Track > fDeltaYLow && 
	     DeltaX_WC_TPC_Track < fDeltaXHigh &&
             DeltaY_WC_TPC_Track < fDeltaYHigh && 
             UpStreamTrjPointZ[aa] < fMaxZPos){
	      
	      // ### Filling Histograms ###
	      hDeltaX->Fill(DeltaX_WC_TPC_Track);
	      hDeltaY->Fill(DeltaY_WC_TPC_Track);
	      hAlpha->Fill(alpha);
	      nWC_TPC_TrackMatch++;

	      // We store all possible wc-tpc track combination
	      // With their figure of merith as a key of the map.
	      std::pair <int,int>  matchCandidate (wcTrackKey[bb],tpcTrackKey[aa]);
	      // This figure of merith is a poor choice and needs to be revised
	      // We would do a better job with some MC tuning
	      // ROOT has an instrument that gives you the pdf of multiple cuts... ask Jeremy Hewes
	      float figOfMerit = DeltaX_WC_TPC_Track*2.5 + DeltaY_WC_TPC_Track*2.22 + alpha; 
	      mapCandidates[figOfMerit] = matchCandidate; 
	    }
	}//<---End bb loop
    }//<---End aa loop



  /////////////////////////////////////////////////
  //  Find the best match with then give metric  //
  ///////////////////////////////////////////////// 
  
  // With these sets we're going to check if we already used the wcTrack and the tpcTrack for the match
  std::set<int> tpcIndeces; 
  std::set<int>  wcIndeces;
  tpcIndeces.clear();
  wcIndeces.clear();

  // We loop on the keys of the map in reverse order:
  // from the highest to the lowest key: we select the 
  // best match first.
  int matchCounter = 0;
  for (std::map<float,std::pair<int, int>>::iterator it=mapCandidates.begin(); it!=mapCandidates.end(); ++it)
    {
      //Check if we used the wc and tpc indeces already
      bool wcB  =  wcIndeces.find(it->second.first ) !=  wcIndeces.end();
      bool tpcB = tpcIndeces.find(it->second.second) != tpcIndeces.end();
      if ( !wcB && !tpcB)
	{
	  //If we haven't, now we are
	  wcIndeces.insert(it->second.first) ;
	  tpcIndeces.insert(it->second.second);
	  // We create the assn here
	  util::CreateAssn(*this, evt, tracklist.at(it->second.second),wctrack.at(it->second.first), *wcTpcTrackAssn);
	  // Enable a max number of associations
	  ++matchCounter;
	  if(matchCounter >= fMaxMatchedTracks) break;
	}	
    }
  

  evt.put(std::move(wcTpcTrackAssn));
  
}//<----End Event Loop

// ---------------------- Begin Job ---------------------------
void WC2TPCTrackMatch::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  hDeltaX = tfs->make<TH1F>("hDeltaX", "#Delta X between TPC and WCtrk (TPC - WC)", 400, -20, 20);
  hDeltaY = tfs->make<TH1F>("hDeltaY", "#Delta Y between TPC and WCtrk (TPC - WC)", 400, -20, 20);
  hAlpha  = tfs->make<TH1F>("hAlpha", "#alpha (Angle between TPC and WC Track)", 110, -5, 50);
  hnMatch = tfs->make<TH1F>("hnMatch", "Number of matched TPC/WC Tracks", 20, -1, 9);
}

void WC2TPCTrackMatch::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::endJob()
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

// -------------------- FHICL Parameter Set --------------------
void WC2TPCTrackMatch::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel		= p.get< std::string >("TrackModuleLabel");
  fWCTrackLabel 		= p.get< std::string >("WCTrackLabel");
  falpha			= p.get< double >("alpha", 20.0);
  fDeltaXLow			= p.get< double >("DeltaXLow", -2.0);
  fDeltaXHigh			= p.get< double >("DeltaXHigh", 6.0);
  fDeltaYLow			= p.get< double >("DeltaYLow", -3.0);
  fDeltaYHigh			= p.get< double >("DeltaYHigh", 6.0);
  fMaxZPos			= p.get< double >("MaxZPos", 14.0);
  fMaxMatchedTracks		= p.get<  int   >("MaxMatchedTracks", 1);
}

void WC2TPCTrackMatch::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WC2TPCTrackMatch::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(WC2TPCTrackMatch)


//CreateAssn(WC2TPCTrackMatch&, art::Event&, const art::Ptr<ldp::WCTrack>&, const art::Ptr<recob::Track>&, art::Assns<ldp::WCTrack, recob::Track>&)
