////////////////////////////////////////////////////////////////////////
// Class:       AnaTreeT1034
// Module Type: analyzer
// File:        AnaTreeT1034_module.cc
//
// Generated at Tue Jul 14 11:21:46 2015 by Roberto Acciarri using artmod
// from cetpkgsupport v1_08_06.
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
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindOneP.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "SimpleTypesAndConstants/geo_types.h"
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RawData/ExternalTrigger.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "MCCheater/BackTracker.h"
#include "Simulation/SimChannel.h"
#include "Filters/ChannelFilter.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/ParticleID.h"
#include "RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "RawDataUtilities/TriggerDigitUtility.h"


// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxHits       = 20000; //maximum number of hits
const int kMaxTrackHits  = 1000;  //maximum number of space points
const int kMaxCluster    = 1000;  //maximum number of clusters
const int kMaxWCTracks   = 1000;   //maximum number of wire chamber tracks
const int kMaxTOF        = 100;   //maximum number of TOF objects

namespace lariat 
{
   class AnaTreeT1034;
}

class lariat::AnaTreeT1034 : public art::EDAnalyzer 
{
public:
   explicit AnaTreeT1034(fhicl::ParameterSet const & p);
   virtual ~AnaTreeT1034();

   // Required functions.
   void analyze(art::Event const & e) override;

   // Selected optional functions.
   void beginJob();
   void reconfigure(fhicl::ParameterSet const & p);

private:

   // === Function used to reset all the variables  ===
   void ResetVars();
  
   // === Storing information into TTree ====
   TTree* fTree;
   
   //=== Storing Run Information ===
   int run;			//<---Run Number
   int subrun;			//<---SubRun Number
   int event;			//<---Event Number
   double evttime;		//<---Event Time Stamp
   double efield[3];		//<---Electric Field 
   int t0;
   int trigtime[16];		//<---Trigger time
   
   // === Storing Track Information ===
   int ntracks_reco;		//<---Number of reconstructed tracks
   double trkvtxx[kMaxTrack];	//<---Track Start X Position (in cm)
   double trkvtxy[kMaxTrack];	//<---Track Start Y Position (in cm)
   double trkvtxz[kMaxTrack];	//<---Track Start Z Position (in cm)
   double trkendx[kMaxTrack];	//<---Track End X Position (in cm)
   double trkendy[kMaxTrack];	//<---Track End Y Position (in cm)
   double trkendz[kMaxTrack];	//<---Track End Z Position (in cm)
   
   double trkstartdcosx[kMaxTrack];	//<---Direction of the track in the x coordinate at its start point
   double trkstartdcosy[kMaxTrack];	//<---Direction of the track in the y coordinate at its start point
   double trkstartdcosz[kMaxTrack];	//<---Direction of the track in the z coordinate at its start point
   double trkenddcosx[kMaxTrack];	//<---Direction of the track in the x coordinate at its end point
   double trkenddcosy[kMaxTrack];	//<---Direction of the track in the y coordinate at its end point
   double trkenddcosz[kMaxTrack];	//<---Direction of the track in the z coordinate at its end point
   double trklength[kMaxTrack];		//<---Calculated length of the track
   double trkmomrange[kMaxTrack];	//<---Calculated track momentum from its length assuming a PID of 13
   double trkmommschi2[kMaxTrack];	//<---Calculated track momentum from multiple scattering using Chi2
   double trkmommsllhd[kMaxTrack];	//<---Calculated track momentum from multiple scattering
   
   
   // === Storing the tracks SpacePoints (Individual 3D points)
   int    ntrkhits[kMaxTrack];			//<---Number of SpacePoints associated with track	
   double trkx[kMaxTrack][kMaxTrackHits];	//<---X position of the spacepoint
   double trky[kMaxTrack][kMaxTrackHits];	//<---Y position of the spacepoint
   double trkz[kMaxTrack][kMaxTrackHits];	//<---Z position of the spacepoint
   double trkpitch[kMaxTrack][3];		//<---The track pitch in a given view (0 == Induction, 1 == Collection)
   
   // === Storing the tracks Calorimetry Information
   int    trkhits[kMaxTrack][2];
   double trkpida[kMaxTrack][2];
   double trkke[kMaxTrack][2];
   double trkdedx[kMaxTrack][2][1000];
   double trkrr[kMaxTrack][2][1000];
   double trkpitchhit[kMaxTrack][2][1000];
   
   // === Storing 2-d Hit information ===
   int    nhits;		//<---Number of 2-d hits in the event
   int    hit_plane[kMaxHits];	//<---Plane number of the hit
   int    hit_wire[kMaxHits];	//<---Wire number of the hit
   int    hit_channel[kMaxHits];//<---Channel number of the hit
   double hit_peakT[kMaxHits];	//<---Peak Time of the hit (in drift ticks)
   double hit_charge[kMaxHits];	//<---Number of ADC's assoicated with the hit (not in units of actual charge)
   double hit_ph[kMaxHits];	//<---Amplitude of the hit (usually the gaussian associated with the hit)
   double hit_tstart[kMaxHits];	//<---Start time of the hit (in drift ticks)
   double hit_tend[kMaxHits];	//<---End time of the hit (in drift ticks)
   int    hit_trkid[kMaxHits];	//<---The track ID associated with this hit (NEED TO CHECK IF THIS IS RIGHT)
   int    hit_pk[kMaxHits];	//<---This stores the raw ADC value (pedestal subtracted) for this hit
   int    hit_t[kMaxHits];	//<---Stores the time tick value for the current hit
   int    hit_ch[kMaxHits];	//<---Stores the hits "charge" in integrated raw ADC (pedestal subtracted)
   int    hit_fwhh[kMaxHits];	//<---Calculated Full width / half max value for the hit
   
   // === Storing 2-d Cluster Information ===
   int    nclus;
   double clustertwire[kMaxCluster];
   double clusterttick[kMaxCluster];
   double cluendwire[kMaxCluster];
   double cluendtick[kMaxCluster];
   int    cluplane[kMaxCluster];
   
   
   // === Storing Wire Chamber Track Information ===
   int nwctrks;
   double wctrk_XFaceCoor[kMaxWCTracks];	//<---The projected X position of the wctrack at the front face of the TPC
   double wctrk_YFaceCoor[kMaxWCTracks];	//<---The projected Y position of the wctrack at the front face of the TPC
   double wctrk_momentum[kMaxWCTracks];		//<---Reconstructed moomentum
   double wctrk_theta[kMaxWCTracks];		//<---angle of track w.r.t. z axis
   double wctrk_phi[kMaxWCTracks];		//<---angle of track w.r.t. x axis
   double wctrk_XDist[kMaxWCTracks];     	//<---X distance between upstream and downstream tracks
   double wctrk_YDist[kMaxWCTracks];     	//<---Y distance between upstream and downstream tracks
   double wctrk_ZDist[kMaxWCTracks];     	//<---Z distance between upstream and downstream tracks
   double XWireHist[kMaxWCTracks][1000];		//<---Coord in terms of wire number
   double YWireHist[kMaxWCTracks][1000];		//<---Coord in terms of wire number
   double XAxisHist[kMaxWCTracks][1000];		//<---coord in terms of units.
   double YAxisHist[kMaxWCTracks][1000];		//<---coord in terms of units.
   double Y_Kink[kMaxWCTracks];     		//<---angle in Y between upstream and downstream tracks.
   
   
   // === Storing Time of Flight information ===
   int ntof;
   double tofObject[kMaxTOF];		//<---The TOF calculated (in ns?) for this TOF object
   double tof_timestamp[kMaxTOF];	//<---Time Stamp for this TOF object
   
   
   // ==== NEED TO FIX THESE VARIABLES....FILLED WITH DUMMY VALUES FOR NOW ===
   
   
   
   double hit_rms[kMaxHits];
   double hit_nelec[kMaxHits];
   double hit_energy[kMaxHits];
   int    hit_trkkey[kMaxHits];
   int    hit_clukey[kMaxHits];
   
   
   std::string fTrigModuleLabel;
   std::string fClusterModuleLabel;
   std::string fHitsModuleLabel;
   std::string fTrackModuleLabel;
   std::string fCalorimetryModuleLabel;
   std::string fParticleIDModuleLabel;
   std::string fWCTrackLabel; 		// The name of the producer that made tracks through the MWPCs
   std::string fTOFModuleLabel;		// Name of the producer that made the TOF objects
};


lariat::AnaTreeT1034::AnaTreeT1034(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
   this->reconfigure(pset);
}

lariat::AnaTreeT1034::~AnaTreeT1034()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034::reconfigure(fhicl::ParameterSet const & pset)
{
   fTrigModuleLabel	 	= pset.get< std::string >("TriggerUtility");
   fHitsModuleLabel      	= pset.get< std::string >("HitsModuleLabel");
   fTrackModuleLabel		= pset.get< std::string >("TrackModuleLabel");
   fCalorimetryModuleLabel 	= pset.get< std::string >("CalorimetryModuleLabel");
   fParticleIDModuleLabel  	= pset.get< std::string >("ParticleIDModuleLabel");
   fClusterModuleLabel          = pset.get< std::string >("ClusterModuleLabel");
   fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel");
   
   return;
}

void lariat::AnaTreeT1034::analyze(art::Event const & evt)
{
   // #############################################
   // ### Reset variables before we get started ###
   // #############################################
   ResetVars();
   //std::cout<<"Check1"<<std::endl;

   // #######################################
   // ### Get potentially useful services ###
   // #######################################
   // === Geometry Service ===
   art::ServiceHandle<geo::Geometry> geom;
   // === Liquid Argon Properties Services ===
   art::ServiceHandle<util::LArProperties> larprop;
   // === Detector properties service ===
   art::ServiceHandle<util::DetectorProperties> detprop;
   // === BackTracker service ===
   art::ServiceHandle<cheat::BackTracker> bt;
  
   // === Run Number ===
   run = evt.run();
   // === Sub-Run Number ===
   subrun = evt.subRun();
   // === Event Number ===
   event = evt.id().event();
   
   std::cout<<std::endl;
   std::cout<<"========================================="<<std::endl;
   std::cout<<"Run = "<<run<<", SubRun = "<<subrun<<", Evt = "<<event<<std::endl;
   std::cout<<"========================================="<<std::endl;
   std::cout<<std::endl;
   
   // === Event Time ===
   art::Timestamp ts = evt.time();
   TTimeStamp tts(ts.timeHigh(), ts.timeLow());
   evttime = tts.AsDouble();
   
   // === Electric Field ===
   efield[0] = larprop->Efield(0);
   efield[1] = larprop->Efield(1);
   efield[2] = larprop->Efield(2);
   
   // === Trigger Offset ====
   t0 = detprop->TriggerOffset();
   
   // ###################################
   // ### Getting Trigger Information ###
   // ###################################
   art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
   std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
   
   if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
      {art::fill_ptr_vector(triglist, trigListHandle);}
   
   // #################################
   // ### Looping over the Triggers ###
   // #################################
   for (size_t i = 0; i<triglist.size(); ++i)
      {
      // === Recording the time of the various triggers ===
      trigtime[i] = triglist[i]->GetTrigTime();
      }
   
   // #####################################
   // ### Getting the Track Information ###
   // #####################################
   art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
   
   // === Filling the tracklist from the tracklistHandle ===
   if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
      {art::fill_ptr_vector(tracklist, trackListHandle);}
   
   // ###################################
   // ### Getting the Hit Information ###
   // ###################################
   art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define tracklist as a pointer to recob::tracks
   
   // === Filling the hitlist from the hitlistHandle ===
   if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      {art::fill_ptr_vector(hitlist, hitListHandle);}
   
   // ##########################################
   // ### Getting the 2D Cluster Information ###
   // ##########################################
   art::Handle< std::vector<recob::Cluster> > clusterListHandle; //<---Define clusterListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Cluster> > clusterlist; //<---Define cluster as a pointer to recob::cluster
   
   // === Filling the clusterlist from the clusterlistHandle ===
   if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      {art::fill_ptr_vector(clusterlist, clusterListHandle);}
      
   // ###################################################
   // ### Getting the Wire Chamber Track  Information ###
   // ###################################################
   art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
   std::vector<art::Ptr<ldp::WCTrack> > wctrack;
   
   if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
      {art::fill_ptr_vector(wctrack, wctrackHandle);}
      
   // ####################################################
   // ### Getting the Time of Flight (TOF) Information ###
   // ####################################################
   art::Handle< std::vector<ldp::TOF> > TOFColHandle;
   std::vector<art::Ptr<ldp::TOF> > tof;
   
   if(evt.getByLabel(fTOFModuleLabel,TOFColHandle))
      {art::fill_ptr_vector(tof, TOFColHandle);}
   
   // ##########################################################
   // ### Grabbing associations for use later in the AnaTool ###
   // ##########################################################
   //std::cout<<"Check2"<<std::endl;
   // === Associations between hits and raw digits ===
   art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
   // === Association between SpacePoints and Tracks ===
   art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
   // === Association between Tracks and 2d Hits ===
   art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
   // === Association between Calorimetry objects and Tracks ===
   art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
   // === Association between Particle ID objects (PID) and Tracks ===
   art::FindManyP<anab::ParticleID>  fmpid(trackListHandle, evt, fParticleIDModuleLabel);
   // ==== Association between Clusters and Hits ===
   art::FindManyP<recob::Cluster>     fmc(hitListHandle,   evt, fClusterModuleLabel);
  
   
   // ### Something to do with SimChannels...need to come back to ###
   std::vector<const sim::SimChannel*> fSimChannels;
   try
      {evt.getView("largeant", fSimChannels);}
   catch (art::Exception const&e){ }

   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   //							FILLING THE WIRE CHAMBER TRACK INFORMATION
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   nwctrks = wctrack.size();
   
   //int wct_count = 0;
   //std::cout<<nwctrks<<std::endl;
   // ########################################
   // ### Looping over Wire Chamber Tracks ###
   // ########################################
   for(size_t wct_count = 0; wct_count < wctrack.size(); wct_count++)
   //for(const auto& wctrack : (*wctrackHandle)) //trackHandle works somewhat like a pointer to a vector of tracks, so dereference the handle to loop over
    //the vector, then use each "track" as a ldp::WCTrack
      {
      
      //std::cout<<"wctrack[wct_count]->Momentum() = "<<wctrack[wct_count]->Momentum()<<std::endl;
      // ##############################################
      // ### Filling Wire Chamber Track information ###
      // ##############################################
      wctrk_XFaceCoor[wct_count] = wctrack[wct_count]->XYFace(0); //indices are 0 for x and 1 for y according to header for WCTrack
      wctrk_YFaceCoor[wct_count] = wctrack[wct_count]->XYFace(1); //indices are 0 for x and 1 for y according to header for WCTrack
      
      wctrk_momentum[wct_count] = wctrack[wct_count]->Momentum();
      wctrk_theta[wct_count] = wctrack[wct_count]->Theta();
      wctrk_phi[wct_count] = wctrack[wct_count]->Phi();
      wctrk_XDist[wct_count] = wctrack[wct_count]->DeltaDist(0);
      wctrk_YDist[wct_count] = wctrack[wct_count]->DeltaDist(1);
      wctrk_ZDist[wct_count] = wctrack[wct_count]->DeltaDist(2);
      Y_Kink[wct_count] = wctrack[wct_count]->YKink();
      
      // === Getting individual channel information ===
      for(size_t chIt = 0; 2*chIt+1 < wctrack[wct_count]->NHits(); ++chIt)
         {
	 XWireHist[wct_count][chIt] = wctrack[wct_count]->HitWire(2*chIt);
	 YWireHist[wct_count][chIt] = wctrack[wct_count]->HitWire(2*chIt+1);
	 
	 XAxisHist[wct_count][chIt] = wctrack[wct_count]->WC(2*chIt);
	 YAxisHist[wct_count][chIt] = wctrack[wct_count]->WC(2*chIt+1);
      
          }//<---End chIt loop
      
      }//<---end wctrack auto loop
      
	 
   
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   //							FILLING THE TIME OF FLIGHT INFORMATION
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   ntof = tof.size();
   // ################################
   // ### Looping over TOF objects ###
   // ################################
   for(size_t tof_count = 0; tof_count < tof.size(); tof_count++)
      {
      
      tofObject[tof_count] =  tof[tof_count]->SingleTOF(tof_count);
      tof_timestamp[tof_count] = tof[tof_count]->TimeStamp(tof_count);
      
      
      
      }//<---End tof_count loop
   
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   //							FILLING THE 2-D CLUSTER INFORMATION
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   
   // ### Saving the number of clusters in this event ###
   nclus = clusterlist.size();
   // ######################################################
   // ### Looping over the all the clusters in the event ###
   // ######################################################
   for (size_t i = 0; i<clusterlist.size(); ++i)
     {
     clustertwire[i] = clusterlist[i]->StartWire();
     clusterttick[i] = clusterlist[i]->StartTick();
     cluendwire[i] = clusterlist[i]->EndWire();
     cluendtick[i] = clusterlist[i]->EndTick();
     cluplane[i] = clusterlist[i]->Plane().Plane;
     }
     
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   //							FILLING THE 3-D TRACK INFORMATION
   // ----------------------------------------------------------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------------------------------------
   
   // ### Calling the track momentum calculator ###
   trkf::TrackMomentumCalculator trkm;
  
  // === Saving the number of tracks per event ===
  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  
  // ### Looping over tracks ###
  for(size_t i=0; i<tracklist.size();++i)
    {
    // ### Clearing the vectors for each track ###
    trackStart.clear();
    trackEnd.clear();
    
    // ### Setting the track information into memory ###
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    tracklist[i]->Extent(trackStart,trackEnd); 
    tracklist[i]->Direction(larStart,larEnd);
    
    // ### Recording the track vertex x, y, z location ###
    trkvtxx[i]        = trackStart[0];
    trkvtxy[i]        = trackStart[1];
    trkvtxz[i]        = trackStart[2];
    
    // ### Recording the track end point x, y, z location ###
    trkendx[i]        = trackEnd[0];
    trkendy[i]        = trackEnd[1];
    trkendz[i]        = trackEnd[2];
    
    // ### Recording the directional cosine at the start of the track ###
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    
    // ### Recording the directional cosine at the end of the track ###
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    
    // ### Recording the track length as calculated by the tracking module ###
    // ####                (each one may do it differently                 ###
    trklength[i]         = tracklist[i]->Length();
    
    // ### Calculating the track momentum using the momentum calculator ###
    trkmomrange[i]    = trkm.GetTrackMomentum(trklength[i],13);
    trkmommschi2[i]   = trkm.GetMomentumMultiScatterChi2(tracklist[i]);
    trkmommsllhd[i]   = trkm.GetMomentumMultiScatterLLHD(tracklist[i]);
    
    // ### Grabbing the SpacePoints associated with this track ###
    ntrkhits[i] = fmsp.at(i).size();
    std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
    
    // ########################################
    // ### Looping over all the SpacePoints ###
    // ########################################
    for (size_t j = 0; j<spts.size(); ++j)
      {
      // ### Recording the x, y, z location of every spacepoint in the track ###
      trkx[i][j] = spts[j]->XYZ()[0];
      trky[i][j] = spts[j]->XYZ()[1];
      trkz[i][j] = spts[j]->XYZ()[2];
      }//<----End SpacePoint loop (j)


    // ###########################################
    // ### Calculating the pitch for the track ###
    // ###########################################
    // === Looping over our two planes (0 == Induction, 1 == Collection) ===
    for (int j = 0; j<2; ++j)
      {
      // ### Putting in a protective try in case we can't set the pitch ###
      try
        {
	// ### If we are in the induction plane calculate the tracks pitch in that view ###
	if (j==0)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kU);
	// ### If we are in the collection plane calculate the tracks pitch in that view ###
	else if (j==1)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kV);
        }//<---End Try statement
      catch( cet::exception &e)
        {mf::LogWarning("AnaTree")<<"caught exeption "<<e<<"\n setting pitch to 0";
	trkpitch[i][j] = 0;
        }//<---End catch statement
    }// <---End looping over planes (j)
     

// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------
//							CALORIMERTY FROM THIS TRACK INFORMATION
// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------      
    
    // ########################################################## 
    // ### Looping over Calorimetry information for the track ###
    // ########################################################## 
    if (fmcal.isValid())
      {
      // ### Putting calo information for this track (i) into pointer vector ###
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(i);
      
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
	trkhits[i][pl] = calos[j]->dEdx().size();
	
	// #### Recording the kinetic energy for this track in this plane ###
	trkke[i][pl] = calos[j]->KineticEnergy();
	
	// ###############################################
	// ### Looping over all the calorimetry points ###
	// ###############################################
	for (size_t k = 0; k<calos[j]->dEdx().size(); ++k)
	  {
	  // ### If we go over 1000 points just skip them ###
	  if (k>=1000) continue;
	  
	  // ### Recording the dE/dX information for this calo point along the track in this plane ###
	  trkdedx[i][pl][k] = calos[j]->dEdx()[k];
	  
	  // ### Recording the residual range for this calo point along the track in this plane ###
	  trkrr[i][pl][k] = calos[j]->ResidualRange()[k];
	  
	  // ### Recording the pitch of this calo point along the track in this plane ###
	  trkpitchhit[i][pl][k] = calos[j]->TrkPitchVec()[k];
	  }//<---End calo points (k)
	  
       }//<---End looping over calo points (j)
       
     }//<---End checking Calo info is valid  
    

// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------
//						PARTICLE ID INFOR FROM THIS TRACK INFORMATION
// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------      
    
    
    // ################################################## 
    // ### Looping over PID information for the track ###
    // ################################################## 
    if (fmpid.isValid())
      {
      // ### Putting PID information for this track (i) into pointer vector ###
      std::vector<art::Ptr<anab::ParticleID> > pids = fmpid.at(i);
      for (size_t j = 0; j<pids.size(); ++j)
        {
	// ### Skip this PID info if not valid for this plane ###
	if (!pids[j]->PlaneID().isValid) continue;
	int pl = pids[j]->PlaneID().Plane;
	// ### Skipping this point if the plane number doesn't make sense ###
	if (pl<0||pl>1) continue;
	trkpida[i][pl] = pids[j]->PIDA();
        
	}//<---End PID loop (j)

      }//<---End checking PID info
      
  }//<---End track loop (i)

  nhits = hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){
    cet::maybe_ref<raw::RawDigit const> rdref(ford.at(i));
    unsigned int channel = hitlist[i]->Channel();
    geo::WireID wireid = hitlist[i]->WireID();
    hit_plane[i]   = wireid.Plane;
    hit_wire[i]    = wireid.Wire;
    hit_channel[i] = channel;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Integral();
    hit_ph[i]      = hitlist[i]->PeakAmplitude();
    hit_tstart[i]  = hitlist[i]->StartTick();
    hit_tend[i]    = hitlist[i]->EndTick();
    if (fmtk.isValid()){
      if (fmtk.at(i).size()!=0){
	hit_trkid[i] = fmtk.at(i)[0]->ID();
	hit_trkkey[i] = fmtk.at(i)[0].key();
      }
    }
    if (fmc.isValid()){
      if (fmc.at(i).size()!=0){
	hit_clukey[i] = fmc.at(i)[0].key();
      }
    }
    if (hit_plane[i]==1){//collection plane
      if( rdref.isValid() ){
	raw::RawDigit const& rd(rdref.ref());
	int dataSize = rd.Samples();
	short ped = rd.GetPedestal();
	std::vector<short> rawadc(dataSize);
	raw::Uncompress(rd.ADCs(),rawadc,rd.Compression());
	int t0 = hit_peakT[i] - 3*(hit_peakT[i]-hit_tstart[i]);
	if (t0<0) t0 = 0;
	int t1 = hit_peakT[i] + 3*(hit_peakT[i]-hit_tstart[i]);
	if (t1>=dataSize) t1 = dataSize-1;
	hit_pk[i] = -1;
	hit_t[i] = -1;
	for (int j = t0; j<=t1; ++j){
	  if (rawadc[j]-ped>hit_pk[i]){
	    hit_pk[i] = rawadc[j]-ped;
	    hit_t[i] = j;
	  }
	}
	hit_ch[i] = 0;
	hit_fwhh[i] = 0;
	double mean_t = 0;
	double mean_t2 = 0;
	for (int j = t0; j<=t1; ++j){
	  if (rawadc[j]-ped>=0.5*hit_pk[i]){
	    ++hit_fwhh[i];
	  }
	  if (rawadc[j]-ped>=0.1*hit_pk[i]){
	    hit_ch[i] += rawadc[j]-ped;
	    mean_t += j*(rawadc[j]-ped);
	    mean_t2 += j*j*(rawadc[j]-ped);
	  }
	}
	mean_t/=hit_ch[i];
	mean_t2/=hit_ch[i];
	hit_rms[i] = sqrt(mean_t2-mean_t*mean_t);
	if (!evt.isRealData()){
	  //	std::vector<sim::IDE> ides;	
	  //	bt->HitToSimIDEs(hitlist[i], ides);
	  hit_nelec[i] = 0;
	  hit_energy[i] = 0;
	  //	for (size_t j = 0; j<ides.size(); ++j){
	  //	  hit_nelec[i] += ides[j].numElectrons;
	  //	  hit_energy[i] += ides[j].energy;
	  //	}
	  const sim::SimChannel* chan = 0;
	  for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
	    if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
	  }
	  if (chan){
	    const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = chan->TDCIDEMap();
	    for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	      // loop over the vector of IDE objects.
	      
	      const std::vector<sim::IDE> idevec = (*mapitr).second;
	      
	      for(size_t iv = 0; iv < idevec.size(); ++iv){ 

		hit_nelec[i] += idevec[iv].numElectrons;
		hit_energy[i] += idevec[iv].energy;
	      }
	    }
	  }
	}
      //std::cout<<hit_wire[i]<<" "<<hit_peakT[i]<<" "<<hit_phlitude[i]<<" "<<hit_tend[i]-hit_tstart[i]<<" "<<hit_t[i]<<" "<<hit_pk[i]<<" "<<hit_ch[i]<<" "<<hit_fwhh[i]<<" "<<hit_rms[i]<<" "<<mean_t<<" "<<mean_t2<<std::endl;

      } // if cet::maybe_ref is valid
    }
  }
  fTree->Fill();



}


void lariat::AnaTreeT1034::beginJob()
{
   
   //std::cout<<"Check-1"<<std::endl;
   // Implementation of optional member function here.
   art::ServiceHandle<art::TFileService> tfs;
   fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("efield",efield,"efield[3]/D");
  fTree->Branch("t0",&t0,"t0/I");
  fTree->Branch("trigtime",trigtime,"trigtime[16]/I");
  fTree->Branch("nclus",&nclus,"nclus/I");
  fTree->Branch("clustertwire",clustertwire,"clustertwire[nclus]/D");
  fTree->Branch("clusterttick",clusterttick,"clusterttick[nclus]/D");
  fTree->Branch("cluendwire",cluendwire,"cluendwire[nclus]/D");
  fTree->Branch("cluendtick",cluendtick,"cluendtick[nclus]/D");
  fTree->Branch("cluplane",cluplane,"cluplane[nclus]/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trkvtxx",trkvtxx,"trkvtxx[ntracks_reco]/D");
  fTree->Branch("trkvtxy",trkvtxy,"trkvtxy[ntracks_reco]/D");
  fTree->Branch("trkvtxz",trkvtxz,"trkvtxz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("trklength",trklength,"trklength[ntracks_reco]/D");
  fTree->Branch("trkmomrange",trkmomrange,"trkmomrange[ntracks_reco]/D");
  fTree->Branch("trkmommschi2",trkmommschi2,"trkmommschi2[ntracks_reco]/D");
  fTree->Branch("trkmommsllhd",trkmommsllhd,"trkmommsllhd[ntracks_reco]/D");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][2]/D");
  fTree->Branch("trkhits",trkhits,"trkhits[ntracks_reco][2]/I");
  fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks_reco][2][1000]/D");
  fTree->Branch("trkrr",trkrr,"trkrr[ntracks_reco][2][1000]/D");
  fTree->Branch("trkpitchhit",trkpitchhit,"trkpitchhit[ntracks_reco][2][1000]/D");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco][2]/D");
  fTree->Branch("trkpida",trkpida,"trkpida[ntracks_reco][2]/D");
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");
  fTree->Branch("hit_tstart",hit_tstart,"hit_tstart[nhits]/D");
  fTree->Branch("hit_tend",hit_tend,"hit_tend[nhits]/D");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");
  fTree->Branch("hit_trkkey",hit_trkkey,"hit_trkkey[nhits]/I");
  fTree->Branch("hit_clukey",hit_clukey,"hit_clukey[nhits]/I");
  fTree->Branch("hit_pk",hit_pk,"hit_pk[nhits]/I");
  fTree->Branch("hit_t",hit_t,"hit_t[nhits]/I");
  fTree->Branch("hit_ch",hit_ch,"hit_ch[nhits]/I");
  fTree->Branch("hit_fwhh",hit_fwhh,"hit_fwhh[nhits]/I");
  fTree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/D");
  fTree->Branch("hit_nelec",hit_nelec,"hit_nelec[nhits]/D");
  fTree->Branch("hit_energy",hit_energy,"hit_energy[nhits]/D");
  
  fTree->Branch("nwctrks",&nwctrks,"nwctrks/I");
   fTree->Branch("wctrk_XFaceCoor",wctrk_XFaceCoor,"wctrk_XFaceCoor[nwctrks]/D");
   fTree->Branch("wctrk_YFaceCoor",wctrk_YFaceCoor,"wctrk_YFaceCoor[nwctrks]/D");
   fTree->Branch("wctrk_momentum",wctrk_momentum,"wctrk_momentum[nwctrks]/D");
   fTree->Branch("wctrk_theta",wctrk_theta,"wctrk_theta[nwctrks]/D");
   fTree->Branch("wctrk_phi",wctrk_phi,"wctrk_phi[nwctrks]/D");
   fTree->Branch("wctrk_XDist",wctrk_XDist,"wctrk_XDist[nwctrks]/D");
   fTree->Branch("wctrk_YDist",wctrk_YDist,"wctrk_YDist[nwctrks]/D");
   fTree->Branch("wctrk_ZDist",wctrk_ZDist,"wctrk_ZDist[nwctrks]/D");
   fTree->Branch("XWireHist",XWireHist,"XWireHist[nwctrks][1000]/D");
   fTree->Branch("YWireHist",YWireHist,"YWireHist[nwctrks][1000]/D");
   fTree->Branch("XAxisHist",XAxisHist,"XAxisHist[nwctrks][1000]/D");
   fTree->Branch("YAxisHist",YAxisHist,"YAxisHist[nwctrks][1000]/D");
   fTree->Branch("Y_Kink",Y_Kink,"Y_Kink[nwctrks]/D");
  
  fTree->Branch("ntof", &ntof, "ntof/I");
  fTree->Branch("tofObject", tofObject, "tofObject[ntof]/D");
  fTree->Branch("tof_timestamp", tof_timestamp, "tof_timestamp[ntof]/D"); 

   
}

void lariat::AnaTreeT1034::ResetVars()
{

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  for (int i = 0; i < 16; ++i){
     trigtime[i]=-99999;
  }
  nclus = -99999;
  for (int i = 0; i < kMaxCluster; ++i){
    clustertwire[i] = -99999;
    clusterttick[i] = -99999;
    cluendwire[i] = -99999;
    cluendtick[i] = -99999;
    cluplane[i] = -99999;
  }
  ntracks_reco = -99999;
  for (int i = 0; i < kMaxTrack; ++i){
    trkvtxx[i] = -99999;
    trkvtxy[i] = -99999;
    trkvtxz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    trklength[i] = -99999;
    trkmomrange[i] = -99999;
    trkmommschi2[i] = -99999;
    trkmommsllhd[i] = -99999;
    ntrkhits[i] = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trkx[i][j] = -99999;
      trky[i][j] = -99999;
      trkz[i][j] = -99999;
    }
    for (int j = 0; j<2; ++j){
      trkpitch[i][j] = -99999;
      trkhits[i][j] = -99999; 
      trkke[i][j] = -99999;
      trkpida[i][j] = -99999;
      for (int k = 0; k<1000; ++k){
	trkdedx[i][j][k] = -99999;
	trkrr[i][j][k] = -99999;
	trkpitchhit[i][j][k] = -99999;
      }
    }
  }
  nhits = -99999;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99999;
    hit_wire[i] = -99999;
    hit_channel[i] = -99999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
    hit_trkid[i] = -99999;
    hit_trkkey[i] = -99999;
    hit_clukey[i] = -99999;
    hit_tstart[i] = -99999;
    hit_tend[i] = -99999;
    hit_pk[i] = -99999;
    hit_t[i] = -99999;
    hit_ch[i] = -99999;
    hit_fwhh[i] = -99999;
    hit_rms[i] = -99999;
    hit_nelec[i] = -99999;
    hit_energy[i] = -99999;
  }
  
  nwctrks = -99999;
  for (int i = 0; i < kMaxWCTracks; i++)
  	{
   	wctrk_XFaceCoor[i] = -99999;	//<---The projected X position of the wctrack at the front face of the TPC
   	wctrk_YFaceCoor[i] = -99999;	//<---The projected Y position of the wctrack at the front face of the TPC
   	wctrk_momentum[i] = -99999;		//<---Reconstructed moomentum
   	wctrk_theta[i] = -99999;		//<---angle of track w.r.t. z axis
   	wctrk_phi[i] = -99999;		//<---angle of track w.r.t. x axis
   	wctrk_XDist[i] = -99999;     	//<---X distance between upstream and downstream tracks
   	wctrk_YDist[i] = -99999;     	//<---Y distance between upstream and downstream tracks
   	wctrk_ZDist[i] = -99999;
	
	for(int j = 0; j < 1000; j++)
		{
   		XWireHist[i][j] = -99999;		//<---Coord in terms of wire number
   		YWireHist[i][j] = -99999;		//<---Coord in terms of wire number
   		XAxisHist[i][j] = -99999;		//<---coord in terms of units.
   		YAxisHist[i][j] = -99999;		//<---coord in terms of units.
   		}//<---End j loop
   	Y_Kink[i] = -99999;
   
   	}//<---End I loop
  
  ntof = -99999;
  for (int i = 0; i < kMaxTOF; i++)
  	{
	tofObject[i] = -99999;
	tof_timestamp[i] = -99999;
	
	
	}//<---End i loop
  
}

DEFINE_ART_MODULE(lariat::AnaTreeT1034)
