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
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/RecoBaseArt/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
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

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTrack      = 1000;     //maximum number of tracks
const int kMaxHits       = 20000;    //maximum number of hits
const int kMaxTrackHits  = 1000;     //maximum number of space points
const int kMaxTrajHits   = 1000;     //maximum number of trajectory points
const int kMaxCluster    = 1000;     //maximum number of clusters
const int kMaxWCTracks   = 1000;     //maximum number of wire chamber tracks
const int kMaxTOF        = 100;      //maximum number of TOF objects
const int kMaxAG         = 100;      //maximum number of AG objects
const int kMaxPrimaryPart = 10;	    //maximum number of true primary particles
const int kMaxPrimaries  = 20000;    //maximum number of true particles tracked
const int kMaxShower     = 100;      //maximum number of Reconstructed showers
const int kMaxMCShower   = 1000;     //maximum number of MCShower Object
const int kMaxTruePrimaryPts = 5000; //maximum number of points in the true primary trajectory 
const int kMaxIDE = 5000; //maximum number of points in the true primary trajectory 

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
  //int trigtime[16];		//<---Trigger time
   
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
  int    trkWCtoTPCMath;		//<---Using an association to see if there was a match between WC and TPC
  					//    0 = match, 1 = no match
   
   
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
  double trkdqdx[kMaxTrack][2][1000];
  double trkrr[kMaxTrack][2][1000];
  double trkpitchhit[kMaxTrack][2][1000]; 
  double trkxyz[kMaxTrack][2][1000][3];

  // === Storing trajectory information for the track ===
  int nTrajPoint[kMaxTrack];			//<---Storing the number of trajectory points
  double pHat0_X[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the x-dir
  double pHat0_Y[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the y-dir
  double pHat0_Z[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the z-dir
  double trjPt_X[kMaxTrack][kMaxTrajHits];     //<---Storing the trajector point location in X
  double trjPt_Y[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Y
  double trjPt_Z[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Z


  std::vector<int>    InteractionPoint;         //<---Geant 4 Primary Trj Point Corresponding to the Interaction
  std::vector<int>    InteractionPointType;     //<---Geant 4 Primary Interaction Type

  // === Geant information for reconstruction track
  int trkg4id[kMaxHits];         //<---geant track id for the track

  int primarytrkkey;             //<---reco track index for primary particle
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
  int    hit_trkkey[kMaxHits]; //<---track index if hit is associated with a track
  float  hit_dQds[kMaxHits];   //<---hit dQ/ds
  float  hit_dEds[kMaxHits];   //<---hit dE/ds
  float  hit_ds[kMaxHits];     //<---hit ds
  float  hit_resrange[kMaxHits];//<---hit residual range
  float  hit_x[kMaxHits];        //<---hit x coordinate
  float  hit_y[kMaxHits];        //<---hit y coordinate
  float  hit_z[kMaxHits];        //<---hit z coordinate
   
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
  float  WC1xPos[kMaxWCTracks];
  float  WC1yPos[kMaxWCTracks];
  float  WC1zPos[kMaxWCTracks];
  float  WC2xPos[kMaxWCTracks];
  float  WC2yPos[kMaxWCTracks];
  float  WC2zPos[kMaxWCTracks];                 //<---The WC positions are relative to the lower front corner of the TPC
  float  WC3xPos[kMaxWCTracks];
  float  WC3yPos[kMaxWCTracks];
  float  WC3zPos[kMaxWCTracks];
  float  WC4xPos[kMaxWCTracks];
  float  WC4yPos[kMaxWCTracks];
  float  WC4zPos[kMaxWCTracks];
   
  // === Storing Time of Flight information ===
  int ntof;
  double tofObject[kMaxTOF];		//<---The TOF calculated (in ns?) for this TOF object
  double tof_timestamp[kMaxTOF];	//<---Time Stamp for this TOF object
   
  // === Storing Aerogel Counter Information ===

  int               nAG;
  long unsigned int HitTimeStamp1p10_1[kMaxAG]; //<---Pulse time stamp relative to (?) in units of (?)
  long unsigned int HitTimeStamp1p10_2[kMaxAG];
  long unsigned int HitTimeStamp1p06_1[kMaxAG];
  long unsigned int HitTimeStamp1p06_2[kMaxAG];

  float             HitPulseArea1p10_1[kMaxAG]; //<---Pulse area in uits of ns*mV(?) for given PMT
  float             HitPulseArea1p10_2[kMaxAG];
  float             HitPulseArea1p06_1[kMaxAG];
  float             HitPulseArea1p06_2[kMaxAG];

  bool              HitExist1p10_1[kMaxAG];     //<---Boolean of whether or not a pulse has been found for the given PMT
  bool              HitExist1p10_2[kMaxAG];
  bool              HitExist1p06_1[kMaxAG];
  bool              HitExist1p06_2[kMaxAG];
   
  // === Storing SimChannel Stuff ===
  int maxTrackIDE;
  double IDEEnergy[kMaxIDE]; 
  double IDEPos[kMaxIDE][3];

  // === Storing Geant4 MC Truth Information ===
  int no_primaries;				//<---Number of primary Geant4 particles in the event
  int geant_list_size;				//<---Number of Geant4 particles tracked
  int pdg[kMaxPrimaries];			//<---PDG Code number of this particle
  double StartPointx[kMaxPrimaries];		//<---X position that this Geant4 particle started at
  double StartPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle started at
  double StartPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle started at
  double Eng[kMaxPrimaries];			//<---Initial Energy of the particle
  double Px[kMaxPrimaries];			//<---Initial Px momentum of the particle
  double Py[kMaxPrimaries];			//<---Initial Py momentum of the particle
  double Pz[kMaxPrimaries];			//<---Initial Pz momentum of the particle

  double EndPointx[kMaxPrimaries];		//<---X position that this Geant4 particle ended at
  double EndPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle ended at
  double EndPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle ended at
  double EndEng[kMaxPrimaries];			//<---End Energy of the particle
  double EndPx[kMaxPrimaries];			//<---End Px momentum of the particle
  double EndPy[kMaxPrimaries];			//<---End Py momentum of the particle
  double EndPz[kMaxPrimaries];			//<---End Pz momentum of the particle

  int Process[kMaxPrimaries];	          	//<---Geant 4 process ID number
  // ### Recording the process as a integer ###
	  // 0 = primary
	  // 1 = PionMinusInelastic
	  // 2 = NeutronInelastic
	  // 3 = hadElastic
	  // 4 = nCapture
	  // 5 = CHIPSNuclearCaptureAtRest
	  // 6 = Decay
	  // 7 = KaonZeroLInelastic
	  // 8 = CoulombScat
	  // 9 = muMinusCaptureAtRest
	  //10 = ProtonInelastic
	  //11 = KaonPlusInelastic
	  //12 = hBertiniCaptureAtRest
  int NumberDaughters[kMaxPrimaries];		//<---Number of Daughters this particle has
  int TrackId[kMaxPrimaries];			//<---Geant4 TrackID number
  int Mother[kMaxPrimaries];			//<---TrackID of the mother of this particle
  int process_primary[kMaxPrimaries];		//<---Is this particle primary (primary = 1, non-primary = 0)
  std::vector<std::string> G4Process;         //<---The process which created this particle
  std::vector<std::string> G4FinalProcess;    //<---The last process which this particle went under

  // === Storing additionnal Geant4 MC Truth Information for the primary track only ===	   
  int NTrTrajPts[kMaxPrimaryPart];							 //<--Nb. of true points in the true primary trajectories
  double MidPosX[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--X position of a point in the true primary trajectory
  double MidPosY[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Y position of a point in the true primary trajectory  
  double MidPosZ[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Z position of a point in the true primary trajectory    
  double MidPx[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<- Px momentum of a point in the true primary trajectory
  double MidPy[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Py momentum of a point in the true primary trajectory
  double MidPz[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Pz momentum of a point in the true primary trajectory
 
  // ==== Storing MCShower MCTruth Information ===
   
  int     no_mcshowers;                         	//number of MC Showers in this event.
  int       mcshwr_origin[kMaxMCShower];            	//MC Shower origin information. 
  int       mcshwr_pdg[kMaxMCShower];            	//MC Shower particle PDG code.   
  int       mcshwr_TrackId[kMaxMCShower];        	//MC Shower particle G4 track ID.
  double    mcshwr_startX[kMaxMCShower];            	//MC Shower particle G4 startX 
  double    mcshwr_startY[kMaxMCShower];          	//MC Shower particle G4 startY 
  double    mcshwr_startZ[kMaxMCShower];          	//MC Shower particle G4 startZ
  double    mcshwr_endX[kMaxMCShower];            	//MC Shower particle G4 endX 
  double    mcshwr_endY[kMaxMCShower];            	//MC Shower particle G4 endY 
  double    mcshwr_endZ[kMaxMCShower];            	//MC Shower particle G4 endZ
  double    mcshwr_CombEngX[kMaxMCShower];            	//MC Shower Combined energy deposition information, Start Point X Position. 
  double    mcshwr_CombEngY[kMaxMCShower];            	//MC Shower Combined energy deposition information, Start Point Y Position.
  double    mcshwr_CombEngZ[kMaxMCShower];            	//MC Shower Combined energy deposition information, Start Point Z Position.
  double    mcshwr_CombEngPx[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  double    mcshwr_CombEngPy[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  double    mcshwr_CombEngPz[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  double    mcshwr_CombEngE[kMaxMCShower];            //MC Shower Combined energy deposition information, Energy
  double    mcshwr_dEdx[kMaxMCShower];           	//MC Shower dEdx, MeV/cm
  double    mcshwr_StartDirX[kMaxMCShower];      	//MC Shower Direction of begining of shower, X direction 
  double    mcshwr_StartDirY[kMaxMCShower];      	//MC Shower Direction of begining of shower, Y direction 
  double    mcshwr_StartDirZ[kMaxMCShower];      	//MC Shower Direction of begining of shower, Z direction 
  int       mcshwr_isEngDeposited[kMaxMCShower];  	//tells whether if this shower deposited energy in the detector or not.

   							//yes = 1; no =0;
  //MC Shower mother information
  int        mcshwr_Motherpdg[kMaxMCShower];       	//MC Shower's mother PDG code.
  int        mcshwr_MotherTrkId[kMaxMCShower];     	//MC Shower's mother G4 track ID.
  double     mcshwr_MotherstartX[kMaxMCShower];    	//MC Shower's mother  G4 startX .
  double     mcshwr_MotherstartY[kMaxMCShower];    	//MC Shower's mother  G4 startY .
  double     mcshwr_MotherstartZ[kMaxMCShower];    	//MC Shower's mother  G4 startZ .
  double     mcshwr_MotherendX[kMaxMCShower];          //MC Shower's mother  G4 endX   .
  double     mcshwr_MotherendY[kMaxMCShower];          //MC Shower's mother  G4 endY   .
  double     mcshwr_MotherendZ[kMaxMCShower];          //MC Shower's mother  G4 endZ   .
   
  //MC Shower ancestor information
  int        mcshwr_Ancestorpdg[kMaxMCShower];       	//MC Shower's ancestor PDG code.
  int        mcshwr_AncestorTrkId[kMaxMCShower];     	//MC Shower's ancestor G4 track ID.
  double     mcshwr_AncestorstartX[kMaxMCShower];   	//MC Shower's ancestor  G4 startX
  double     mcshwr_AncestorstartY[kMaxMCShower];   	//MC Shower's ancestor  G4 startY
  double     mcshwr_AncestorstartZ[kMaxMCShower];   	//MC Shower's ancestor  G4 startZ
  double     mcshwr_AncestorendX[kMaxMCShower];     	//MC Shower's ancestor  G4 endX 
  double     mcshwr_AncestorendY[kMaxMCShower];     	//MC Shower's ancestor  G4 endY
  double     mcshwr_AncestorendZ[kMaxMCShower];      	//MC Shower's ancestor  G4 endZ    
   
  // === Storing Shower Reco Information using ShowerReco3D ===

  int nshowers; ///number of showers per event
  int shwID[kMaxShower];//ID of the reco shower
  double CosStartShw[3][kMaxShower];
  double CosStartSigmaShw[3][kMaxShower];
  double CosStartXYZShw[3][kMaxShower];
  double CosStartXYZSigmaShw[3][kMaxShower];
  double TotalEShw[2][kMaxShower];/// total energy of the shower (under investigation...)
  //double TotalESigmaShw[2][kMaxShower];// not working
  double dEdxPerPlaneShw[2][kMaxShower];
  //double dEdxSigmaPerPlaneShw[2][kMaxShower];//not working
  double TotalMIPEShw[2][kMaxShower];
  //double TotalMIPESigmaShw[2][kMaxShower];//not working
  int BestPlaneShw[kMaxShower];	
  double LengthShw[kMaxShower];
   
   
   
  // ==== NEED TO FIX THESE VARIABLES....FILLED WITH DUMMY VALUES FOR NOW ===
  double hit_rms[kMaxHits];
  double hit_nelec[kMaxHits];
  double hit_energy[kMaxHits];
  //int    hit_trkkey[kMaxHits];
  int    hit_clukey[kMaxHits];
  
   
  //std::string fTrigModuleLabel;
  std::string fClusterModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fWCTrackLabel; 		// The name of the producer that made tracks through the MWPCs
  std::string fTOFModuleLabel;		// Name of the producer that made the TOF objects
  std::string fAGModuleLabel;         // Name of the producer that made the aerogel objects
  std::string fG4ModuleLabel;
  std::string fShowerModuleLabel;       // Producer that makes showers from clustering
  std::string fMCShowerModuleLabel;	// Producer name that makes MCShower Object
  std::string fSimChanModuleLabel; // Producer that makes SimIDE information
  std::string fWC2TPCModuleLabel;	// Producer which creates an association between WC and TPC Track

  calo::CalorimetryAlg fCalorimetryAlg;

};


lariat::AnaTreeT1034::AnaTreeT1034(fhicl::ParameterSet const & pset) 
  : EDAnalyzer(pset)
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(pset);
}

lariat::AnaTreeT1034::~AnaTreeT1034()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034::reconfigure(fhicl::ParameterSet const & pset)
{
  //fTrigModuleLabel	 	= pset.get< std::string >("TriggerUtility");
  fHitsModuleLabel      	= pset.get< std::string >("HitsModuleLabel");
  fTrackModuleLabel		= pset.get< std::string >("TrackModuleLabel");
  fCalorimetryModuleLabel 	= pset.get< std::string >("CalorimetryModuleLabel");
  fParticleIDModuleLabel  	= pset.get< std::string >("ParticleIDModuleLabel");
  fClusterModuleLabel          	= pset.get< std::string >("ClusterModuleLabel");
  fWCTrackLabel 		= pset.get< std::string >("WCTrackLabel");
  fTOFModuleLabel 		= pset.get< std::string >("TOFModuleLabel");
  fAGModuleLabel               	= pset.get< std::string >("AGModuleLabel");
  fG4ModuleLabel               	= pset.get< std::string >("G4ModuleLabel");
  fShowerModuleLabel           	= pset.get< std::string >("ShowerModuleLabel");
  fSimChanModuleLabel	        = pset.get< std::string >("SimChanModuleLabel");
  fMCShowerModuleLabel		= pset.get< std::string >("MCShowerModuleLabel");
  fWC2TPCModuleLabel      	= pset.get< std::string >("WC2TPCModuleLabel"     , "WC2TPCtrk");
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
  //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  // === Detector properties service ===
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // === BackTracker service ===
  art::ServiceHandle<cheat::BackTracker> bt;
  const sim::ParticleList& plist = bt->ParticleList();
   
   
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
  // Note: LArProperties::Efield() has moved to DetectorProperties/DetectorPropertiesService
  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);
   
  // === Trigger Offset ====
  t0 = detprop->TriggerOffset();
   
  /*
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
  */
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

  // ####################################################
  // ### Getting the Aerogel Information ###
  // ####################################################
  art::Handle< std::vector<ldp::AGCounter> > AGColHandle;
  std::vector<art::Ptr<ldp::AGCounter> > agc;

  if(evt.getByLabel(fAGModuleLabel,AGColHandle))
    {art::fill_ptr_vector(agc, AGColHandle);}
      
  // #####################################
  // ### Getting the Shower Information ###
  // #####################################
  art::Handle< std::vector<recob::Shower> > shwListHandle; 
  std::vector<art::Ptr<recob::Shower> > shwlist;
   
  // === Filling the shwlist from the shwlistHandle ===
  if (evt.getByLabel(fShowerModuleLabel,shwListHandle))
    {art::fill_ptr_vector(shwlist, shwListHandle);}

   
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
  // ==== Association between Clusters and Showers ===
  art::FindManyP<recob::Shower> fms (clusterListHandle, evt, fShowerModuleLabel);
  // ==== Association between Tracks and Hits
  art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
  
   
  // ### Something to do with SimChannels...need to come back to ###
  std::vector<const sim::SimChannel*> fSimChannels;
  try
    {evt.getView("largeant", fSimChannels);}
  catch (art::Exception const&e){ }
   
  // ###################################################################
  // ### Setting a boolian to only output MC info if this is MC-info ###
  // ###################################################################
  bool isdata = false;
  if (evt.isRealData())
    {isdata = true;}
	
  else isdata = false;
   
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE Sim Channel INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isdata) {

    art::Handle< std::vector<sim::SimChannel> > SimListHandle; 
    std::vector<art::Ptr<sim::SimChannel> > Simlist;    
    if(evt.getByLabel(fSimChanModuleLabel, SimListHandle))
      { art::fill_ptr_vector(Simlist, SimListHandle); }
  
    maxTrackIDE = 0;

    // Loop over the channels, the wires
    for(size_t nChan = 0; nChan < Simlist.size(); nChan++) {

      // Only getting one plane
      if(Simlist.at(nChan)->Channel() > 240) { break; }
      
      // Get the information of each wire
     const auto & wire = Simlist.at(nChan)->TDCIDEMap();

      // Looping over the IDEs in a wire, or looping over time
      //typedef std::map<unsigned short, std::vector< sim::IDE > >::iterator it_type;
      for(auto it = wire.begin(); it != wire.end(); it++) {

	// Looping over the IDEs in a given time tick
	for(size_t i = 0; i < it->second.size(); i++) {

	  // Only non-showering (nonnegative) primary track IDEs
	  if(it->second.at(i).trackID != 1) { continue; } 
 
	  IDEEnergy[maxTrackIDE] = it->second.at(i).energy; 
	    
	  IDEPos[maxTrackIDE][0] = it->second.at(i).x; 
	  IDEPos[maxTrackIDE][1] = it->second.at(i).y;
	  IDEPos[maxTrackIDE][2] = it->second.at(i).z;

	  maxTrackIDE += 1;
	  
	} // Loop over IDE
      } // Loop over Time
    } // Loop over Wire
  } // End of isdata


   
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE MCTruth Geant4 INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isdata)
    {
      art::Handle< std::vector<sim::MCShower> > mcshowerh;
      evt.getByLabel(fMCShowerModuleLabel, mcshowerh);
      
      // #######################################################
      // ### Check to make sure the MCShower Handle is valid ###
      // #######################################################
      if (mcshowerh.isValid())
	{
	  //std::cout<<mcshowerh->size()<<std::endl;
	  no_mcshowers = mcshowerh->size();
	  size_t shwr = 0;
	  // ############################################################
	  // ### Looping over all MC Showers (using uBooNEAna method) ###
	  // ############################################################
	  for(std::vector<sim::MCShower>::const_iterator imcshwr = mcshowerh->begin(); imcshwr != mcshowerh->end(); ++imcshwr)
	    {
	    
	      const sim::MCShower& mcshwr = *imcshwr;
	    
	      mcshwr_origin[shwr]          = mcshwr.Origin();
	      mcshwr_pdg[shwr]             = mcshwr.PdgCode();
	      mcshwr_TrackId[shwr]         = mcshwr.TrackID();
	      mcshwr_startX[shwr]          = mcshwr.Start().X();
	      mcshwr_startY[shwr]          = mcshwr.Start().Y();
	      mcshwr_startZ[shwr]          = mcshwr.Start().Z();
	      mcshwr_endX[shwr]            = mcshwr.End().X();
	      mcshwr_endY[shwr]            = mcshwr.End().Y(); 
	      mcshwr_endZ[shwr]            = mcshwr.End().Z();
	    
	      // ### Recording Energy inside the detector ###
	      if (mcshwr.DetProfile().E()!= 0)
		{
		  mcshwr_isEngDeposited[shwr] = 1;
		  mcshwr_CombEngX[shwr]        = mcshwr.DetProfile().X(); 
		  mcshwr_CombEngY[shwr]        = mcshwr.DetProfile().Y();
		  mcshwr_CombEngZ[shwr]        = mcshwr.DetProfile().Z();
		  mcshwr_CombEngPx[shwr]       = mcshwr.DetProfile().Px();
		  mcshwr_CombEngPy[shwr]       = mcshwr.DetProfile().Py();
		  mcshwr_CombEngPz[shwr]       = mcshwr.DetProfile().Pz(); 
		  mcshwr_CombEngE[shwr]        = mcshwr.DetProfile().E();
		  mcshwr_dEdx[shwr]            = mcshwr.dEdx();
		  mcshwr_StartDirX[shwr]       = mcshwr.StartDir().X();
		  mcshwr_StartDirY[shwr]       = mcshwr.StartDir().Y();
		  mcshwr_StartDirZ[shwr]       = mcshwr.StartDir().Z();
		}//<---End recording info inside the detector
	      else
		{mcshwr_isEngDeposited[shwr] = 0;}
	    
	      mcshwr_Motherpdg[shwr]       = mcshwr.MotherPdgCode();
	      mcshwr_MotherTrkId[shwr]     = mcshwr.MotherTrackID();
	      mcshwr_MotherstartX[shwr]    = mcshwr.MotherStart().X();
	      mcshwr_MotherstartY[shwr]    = mcshwr.MotherStart().Y();
	      mcshwr_MotherstartZ[shwr]    = mcshwr.MotherStart().Z();
	      mcshwr_MotherendX[shwr]      = mcshwr.MotherEnd().X();
	      mcshwr_MotherendY[shwr]      = mcshwr.MotherEnd().Y();
	      mcshwr_MotherendZ[shwr]      = mcshwr.MotherEnd().Z();
	      mcshwr_Ancestorpdg[shwr]     = mcshwr.AncestorPdgCode();
	      mcshwr_AncestorTrkId[shwr]   = mcshwr.AncestorTrackID();
	      mcshwr_AncestorstartX[shwr]  = mcshwr.AncestorStart().X();
	      mcshwr_AncestorstartY[shwr]  = mcshwr.AncestorStart().Y();
	      mcshwr_AncestorstartZ[shwr]  = mcshwr.AncestorStart().Z();
	      mcshwr_AncestorendX[shwr]    = mcshwr.AncestorEnd().X();
	      mcshwr_AncestorendY[shwr]    = mcshwr.AncestorEnd().Y();
	      mcshwr_AncestorendZ[shwr]    = mcshwr.AncestorEnd().Z();

      
	      shwr++;
            }//<---End imcshwr iterator loop
	    
	    
	}//<---Only going in if the handle is valid
    }//<---End only looking at MC information
   
   
   
   
   
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE MCTruth Geant4 INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isdata)
    {
      // ######################################
      // ### Making a vector of MCParticles ###
      // ######################################   
      std::vector<const simb::MCParticle* > geant_part;
      
      // ### Looping over all the Geant4 particles from the BackTracker ###
      for(size_t p = 0; p < plist.size(); ++p) 
		{
	  		// ### Filling the vector with MC Particles ###
	  		geant_part.push_back(plist.Particle(p)); 
		}
	
      //std::cout<<"No of geant part= "<<geant_part.size()<<std::endl;
      
      // ### Setting a string for primary ###
      std::string pri("primary");
      
      // ### Setting a string for PionMinusInelastic ###
      std::string PionMinusInelastic("pi-Inelastic");
      
      // ### Setting a string for NeutronInelastic ###
      std::string NeutronInelastic("neutronInelastic");
      
       // ### Setting a string for hadElastic ###
      std::string hadElastic("hadElastic");
      
      // ### Setting a string for nCapture ###
      std::string nCapture("nCapture");
      
      // This may not be called hBertiniCaptureAtRest ?
      // ### Setting a string for CHIPSNuclearCaptureAtRest ###
      std::string CHIPSNuclearCaptureAtRest("CHIPSNuclearCaptureAtRest");
      
      // ### Setting a string for Decay ###
      std::string Decay("Decay");
      
      // ### Setting a string for KaonZeroLInelastic ###
      std::string KaonZeroLInelastic("KaonZeroLInelastic");
      
      // ### Setting a string for CoulombScat ###
      std::string CoulombScat("CoulombScat");
      
      // ### Setting a string for muMinusCaptureAtRest ###
      std::string muMinusCaptureAtRest("muMinusCaptureAtRest");
      
      // ### Setting a string for ProtonInelastic ###
      std::string ProtonInelastic("protonInelastic");
      
      // ### Setting a string for Kaon Inelastic ###
      std::string KaonPlusInelastic("kaon+Inelastic");
      
      // ### Setting a string for BertiniCaptureAtRest
      std::string hBertiniCaptureAtRest("hBertiniCaptureAtRest");
      
      
      int primary=0;
      int geant_particle=0;
      
      // ############################################################
      // ### Determine the number of primary particles from geant ###
      // ############################################################
      for( unsigned int i = 0; i < geant_part.size(); ++i )
		{
	  		geant_particle++;
	  		// ### Counting the number of primary particles ###
	  		if(geant_part[i]->Process()==pri)
	    		{ primary++;}
		}//<---End i loop
	 
	
      // ### Saving the number of primary particles ###
      no_primaries=primary;
      // ### Saving the number of Geant4 particles ###
      geant_list_size=geant_particle;
       
      // ### Looping over all the Geant4 particles ###
      int iPrim = 0;
      for( unsigned int i = 0; i < geant_part.size(); ++i )
         {
   
      	 // ### If this particle is primary, set = 1 ###
	 if(geant_part[i]->Process()==pri)
	    {process_primary[i]=1;}
         // ### If this particle is not-primary, set = 0 ###
	 else
	    {process_primary[i]=0;}
          
	 // ### Recording the process as a integer ###
	 // 0 = primary
	 // 1 = PionMinusInelastic
   	 // 2 = NeutronInelastic
	 // 3 = hadElastic
	 // 4 = nCapture
	 // 5 = CHIPSNuclearCaptureAtRest
	 // 6 = Decay
	 // 7 = KaonZeroLInelastic
	 // 8 = CoulombScat
	 // 9 = muMinusCaptureAtRest
	 //10 = ProtonInelastic
	  
	 if(geant_part[i]->Process() == pri)
	    {Process[i] = 0;}
	     
	 if(geant_part[i]->Process() == PionMinusInelastic)
	    {Process[i] = 1;}
	     
	 if(geant_part[i]->Process() == NeutronInelastic)
	    {Process[i] = 2;}
	     
	 if(geant_part[i]->Process() == hadElastic)
	    {Process[i] = 3;}
	  
	 if(geant_part[i]->Process() == nCapture)
	    {Process[i] = 4;}
	     
	 if(geant_part[i]->Process() == CHIPSNuclearCaptureAtRest)
	    {Process[i] = 5;}
	  
	 if(geant_part[i]->Process() == Decay)
	    {Process[i] = 6;}
	  
	 if(geant_part[i]->Process() == KaonZeroLInelastic)
	    {Process[i] = 7;}
	     
	 if(geant_part[i]->Process() == CoulombScat)
	    {Process[i] = 8;}
	     
	 if(geant_part[i]->Process() == muMinusCaptureAtRest)
	    {Process[i] = 9;}
	     
	 if(geant_part[i]->Process() == ProtonInelastic)
	    {Process[i] = 10;}
			
	 if(geant_part[i]->Process() == KaonPlusInelastic)
	    {Process[i] = 11;}
			
	 if(geant_part[i]->Process() == hBertiniCaptureAtRest)
	    {Process[i] = 12;}
	     
	 //std::cout<<"Process = "<<geant_part[i]->Process()<<std::endl;		       

	 // ### Saving the particles mother TrackID ###
	 Mother[i]=geant_part[i]->Mother();
	 // ### Saving the particles TrackID ###
	 TrackId[i]=geant_part[i]->TrackId();
	 // ### Saving the PDG Code ###
	 pdg[i]=geant_part[i]->PdgCode();
	 // ### Saving the particles start and end Energy ###
	 Eng[i]=geant_part[i]->E();
	 EndEng[i]=geant_part[i]->EndE();
	  
	 // ### Saving the start and end Px, Py, Pz info ###
	 Px[i]=geant_part[i]->Px();
	 Py[i]=geant_part[i]->Py();
	 Pz[i]=geant_part[i]->Pz();
	 EndPx[i]=geant_part[i]->EndPx();
	 EndPy[i]=geant_part[i]->EndPy();
	 EndPz[i]=geant_part[i]->EndPz();
	  
	 // ### Saving the Start and End Point for this particle ###
	 StartPointx[i]=geant_part[i]->Vx();
	 StartPointy[i]=geant_part[i]->Vy();
	 StartPointz[i]=geant_part[i]->Vz();
	 EndPointx[i]=geant_part[i]->EndPosition()[0];
	 EndPointy[i]=geant_part[i]->EndPosition()[1];
	 EndPointz[i]=geant_part[i]->EndPosition()[2];

	 // ### Saving the processes for this particle ###
	 //std::cout<<"finding proc"<<std::endl;
	 G4Process.push_back( geant_part[i]->Process() );
	 G4FinalProcess.push_back( geant_part[i]->EndProcess() );
 	  
	 // ### Saving the number of Daughters for this particle ###
	 NumberDaughters[i]=geant_part[i]->NumberDaughters();

	 // ### Save intermediary information for the primary track
	 if(geant_part[i]->Process()==pri){
	 NTrTrajPts[i]=geant_part[i]->NumberTrajectoryPoints();
	 simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
				
	 int iPrimPt = 0;	
	 for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	    {
	    MidPosX[iPrim][iPrimPt] = truetraj.X(iPrimPt);
	    MidPosY[iPrim][iPrimPt] = truetraj.Y(iPrimPt);
	    MidPosZ[iPrim][iPrimPt] = truetraj.Z(iPrimPt);
	    MidPx[iPrim][iPrimPt] = truetraj.Px(iPrimPt);
	    MidPy[iPrim][iPrimPt] = truetraj.Py(iPrimPt);
	    MidPz[iPrim][iPrimPt] = truetraj.Pz(iPrimPt);
	    iPrimPt++;
	    }//<--End loop on true trajectory points
	 

	 	   // Yet an other scheme for interaction type
	             
	   // ### Recording the process as a integer ###
	   // 0 = NoInteractionNodaughters, thought going
	   // 1 = PionMinusInelastic
	   // 2 = NeutronInelastic
	   // 3 = hadElastic
	   // 4 = nCapture
	   // 5 = CHIPSNuclearCaptureAtRest
	   // 6 = Decay
	   // 7 = KaonZeroLInelastic
	   // 8 = CoulombScat
	   // 9 = muMinusCaptureAtRest
	   //10 = ProtonInelastic
	   //11 = Kaon+Inelastic
	   //12 = Kaon-Inelastic 
	   //13 = protonInelastic
	   //14 = Pi+Inelastic 


	   auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
	   
	   // Ok, if the size of the map is 0, all the action might happen at the end of the track
	   // So we check the daugthers:
	   //    - Case 1. There are daugthers:
	   //               * The interesting point is the last one
	   //               * The interaction type is the one that created the first daugther (this might be improved)
	   //    - Case 2. There are NO daugthers:
	   //              * We assign the interaction type to be 0: nothing happens, thought going particle
	   //              * The interesting point is the last one (might not be in the TPC)
	   if (!thisTracjectoryProcessMap.size())
	     {
	       int interestingPoint = (int) (NTrTrajPts[i] - 1);
	       InteractionPoint.push_back(interestingPoint);

	       if (NumberDaughters[i])
		 {		 
		   auto thePrimaryDaughterID = geant_part[i]-> Daughter(0); 
		   for( unsigned int iD = 0; iD < geant_part.size(); ++iD )
		     {
		       if (geant_part[iD]->TrackId() == thePrimaryDaughterID) 
			 {		      
			   if(geant_part[iD]->Process() == PionMinusInelastic)
			     {InteractionPointType.push_back(1);}
			   
			   if(geant_part[iD]->Process() == NeutronInelastic)
			     {InteractionPointType.push_back(2);}
			   
			   if(geant_part[iD]->Process() == hadElastic)
			     {InteractionPointType.push_back(3);}
			   
			   if(geant_part[iD]->Process() == nCapture)
			     {InteractionPointType.push_back(4);}
			   
			   if(geant_part[iD]->Process() == CHIPSNuclearCaptureAtRest)
			     {InteractionPointType.push_back(5);}
			   
			   if(geant_part[iD]->Process() == Decay)
			     {InteractionPointType.push_back(6);}
			   
			   if(geant_part[iD]->Process() == KaonZeroLInelastic)
			     {InteractionPointType.push_back(7);}
			   
			   if(geant_part[iD]->Process() == CoulombScat)
			     {InteractionPointType.push_back(8);}
			   
			   if(geant_part[iD]->Process() == muMinusCaptureAtRest)
			     {InteractionPointType.push_back(9);}
			   
			   if(geant_part[iD]->Process() == ProtonInelastic)
			     {InteractionPointType.push_back(10);}
			   
			   if(geant_part[iD]->Process() == KaonPlusInelastic)
			     {InteractionPointType.push_back(11);}
			   
			   if(geant_part[iD]->Process() == hBertiniCaptureAtRest)
			     {InteractionPointType.push_back(12);}
			 }
		     }
		 }else
		 {
		   InteractionPointType.push_back(0);              
		 }
	     }else
	     {
	       // The map is not zero: somthing interesting might happen in the middle of the track!!
	       for(auto const& couple: thisTracjectoryProcessMap) 
		 {
		   int interestingPoint = (int) couple.first;
		   InteractionPoint.push_back(interestingPoint);         	   
		   if ((truetraj.KeyToProcess(couple.second)).find("hadElastic")!= std::string::npos) InteractionPointType.push_back(3);           
		   if ((truetraj.KeyToProcess(couple.second)).find("pi-Inelastic")    != std::string::npos) InteractionPointType.push_back(1);           
		   if ((truetraj.KeyToProcess(couple.second)).find("pi+Inelastic")    != std::string::npos) InteractionPointType.push_back(14);           
		   if ((truetraj.KeyToProcess(couple.second)).find("kaon-Inelastic")  != std::string::npos) InteractionPointType.push_back(12);           
		   if ((truetraj.KeyToProcess(couple.second)).find("kaon+Inelastic")  != std::string::npos) InteractionPointType.push_back(11);           
		   if ((truetraj.KeyToProcess(couple.second)).find("protonInelastic") != std::string::npos) InteractionPointType.push_back(13);           
		   if ((truetraj.KeyToProcess(couple.second)).find("neutronInelastic")!= std::string::npos) InteractionPointType.push_back(2);           

		 }
	     }	
   
	 iPrim++;
	}//<--End if primary
     }//<--End loop on geant particles
	  
 }//<---End checking if this is MC 
   
   
 
   
   
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
      WC1xPos[wct_count] = wctrack[wct_count]->HitPosition(0,0);
      WC1yPos[wct_count] = wctrack[wct_count]->HitPosition(0,1);
      WC1zPos[wct_count] = wctrack[wct_count]->HitPosition(0,2);
      WC2xPos[wct_count] = wctrack[wct_count]->HitPosition(1,0);
      WC2yPos[wct_count] = wctrack[wct_count]->HitPosition(1,1);
      WC2zPos[wct_count] = wctrack[wct_count]->HitPosition(1,2);
      WC3xPos[wct_count] = wctrack[wct_count]->HitPosition(2,0);
      WC3yPos[wct_count] = wctrack[wct_count]->HitPosition(2,1);
      WC3zPos[wct_count] = wctrack[wct_count]->HitPosition(2,2);
      WC4xPos[wct_count] = wctrack[wct_count]->HitPosition(3,0);
      WC4yPos[wct_count] = wctrack[wct_count]->HitPosition(3,1);
      WC4zPos[wct_count] = wctrack[wct_count]->HitPosition(3,2);
//      std::cout<<"The WC4 x-position is: "<<wctrack[wct_count]->HitPosition(3,0)<<std::endl;
      
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
  size_t tof_counter = 0; // book-keeping
  for(size_t i = 0; i < tof.size(); i++) {

	  size_t number_tof = tof[i]->NTOF();

     for (size_t tof_idx = 0; tof_idx < number_tof; ++tof_idx) {
			tofObject[tof_counter] =  tof[i]->SingleTOF(tof_idx);
			tof_timestamp[tof_counter] = tof[i]->TimeStamp(tof_idx);
			++tof_counter;
     } // loop over TOF

  }//<---End tof_count loop

  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //                                                   FILLING THE AEROGEL COUNTER INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  //   ldp::AGCounter counter;
  //   std::cout<<"counter.GetNHits()"<<counter.GetNHits()<<std::endl;
  //   std::cout<<"counter.size()"<<counter.size()<<std::endl;
  //   std::cout<<"counter->size()"<<counter->size()<<std::endl;

  nAG = agc.size();
  // ################################
  // ### Looping over aerogel counter objects ###
  // ################################
  size_t agc_counter = 0; // book-keeping
  for(size_t i = 0; i < agc.size(); i++)
    {

      size_t number_agc = agc[i]->GetNHits();
    
      
     

      for (size_t agc_idx = 0; agc_idx < number_agc; ++agc_idx) {
	//        for (size_t agc_idx = 0; agc_idx < 1; ++agc_idx) {
	HitTimeStamp1p10_1[agc_counter]=agc[i]->GetHitTimeStamp1p10_1(agc_idx);
	HitTimeStamp1p10_2[agc_counter]=agc[i]->GetHitTimeStamp1p10_2(agc_idx);
	HitTimeStamp1p06_1[agc_counter]=agc[i]->GetHitTimeStamp1p06_1(agc_idx);
	HitTimeStamp1p06_2[agc_counter]=agc[i]->GetHitTimeStamp1p06_2(agc_idx);

	HitPulseArea1p10_1[agc_counter]=agc[i]->GetHitPulseArea1p10_1(agc_idx);
	HitPulseArea1p10_2[agc_counter]=agc[i]->GetHitPulseArea1p10_2(agc_idx);
	HitPulseArea1p06_1[agc_counter]=agc[i]->GetHitPulseArea1p06_1(agc_idx);
	HitPulseArea1p06_2[agc_counter]=agc[i]->GetHitPulseArea1p06_2(agc_idx);

	HitExist1p10_1[agc_counter]=agc[i]->GetHitExist1p10_1(agc_idx);
	HitExist1p10_2[agc_counter]=agc[i]->GetHitExist1p10_2(agc_idx);
	HitExist1p06_1[agc_counter]=agc[i]->GetHitExist1p06_1(agc_idx);
	HitExist1p06_2[agc_counter]=agc[i]->GetHitExist1p06_2(agc_idx);

	++agc_counter;
      } // loop over aerogel pulses

    }//<---End aerogel counters
     

  //        short unsigned int number_agc = agc[i]->GetNHits();

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
  //							FILLING THE SHOWER RECO INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  nshowers = shwlist.size();

  for (size_t i = 0; i<shwlist.size(); ++i)  // loop over showers
    {
      shwID[i] = shwlist[i]->ID();
      BestPlaneShw[i] = shwlist[i]->best_plane();
      LengthShw[i] = shwlist[i]->Length();

      for (size_t j = 0; j<3; ++j)
	{
	  CosStartShw[j][i] = shwlist[i]->Direction()[j];
	  // CosStartSigmaShw[j][i] = shwlist[i]->DirectionErr()[j];
	  CosStartXYZShw[j][i] = shwlist[i]->ShowerStart()[j];
	  //CosStartXYZSigmaShw[j][i] =  shwlist[i]->ShowerStartErr()[j];
	}

      for (int j = 0; j<2; ++j)/// looping over the 2 planes
	{
	  TotalEShw[j][i] = shwlist[i]->Energy()[j];
	  //TotalESigmaShw[j][i] = shwlist[i]->EnergyErr()[j];
	  dEdxPerPlaneShw[j][i] = shwlist[i]->dEdx()[j];
	  TotalMIPEShw[j][i] = shwlist[i]->MIPEnergy()[j];
	} 
    }    // end loop over showers


  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE 3-D TRACK INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
   
  // ### Calling the track momentum calculator ###
  trkf::TrackMomentumCalculator trkm;
  trkm.SetMinLength(10); //change the minimal track length requirement to 10 cm

  // === Saving the number of tracks per event ===
  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  
  // === Association between WC Tracks and TPC Tracks ===
  int TempTrackMatchedID = -1;                                                                                                                               
  /*art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
  
  if (fWC2TPC.isValid())
     {
     std::cout<<"Valid WC2TPC Track"<<std::endl;
     
     std::cout<<"fWC2TPC.size() = "<<fWC2TPC.size()<<std::endl;
     // === Loop on all the Assn WC-TPC tracks ===                                                          
     for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn )
         {
	 std::cout<<"indexAssn = "<<indexAssn<<std::endl;
	 // =========================
	 // === Get the TPC track ===
	 // =========================
	 auto fake = *fWC2TPC.at(indexAssn);
	 std::cout<<fake.ID()<<std::endl;
	 cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));
	 
	 recob::Track const& aTrack(trackWC2TPC.ref());
	 TempTrackMatchedID = aTrack.ID();
	 std::cout<<"TempTrackMatchedID = "<<TempTrackMatchedID<<std::endl;
	 
	 }//<----End indexAssn loop
     
     
     }//<---End checking that the WC2TPC is valid*/
  
  // ### Looping over tracks ###
  double maxtrackenergy = -1;
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
      
      // ### Stroing an integer for the match of a WC to TPC track ###
      int trackMatch = 1;
      if(TempTrackMatchedID == tracklist[i]->ID() )
         {
	 trackMatch = 0;
	 }//<---End match
	 
      // ### Setting the WC to TPC match ###
      trkWCtoTPCMath = trackMatch;
      
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
		trkpitch[i][j] = lar::util::TrackPitchInView(*tracklist[i], geo::kU);
	      // ### If we are in the collection plane calculate the tracks pitch in that view ###
	      else if (j==1)
		trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kV);
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
		  trkdqdx[i][pl][k] = calos[j]->dQdx()[k];
	  
		  // ### Recording the residual range for this calo point along the track in this plane ###
		  trkrr[i][pl][k] = calos[j]->ResidualRange()[k];
	  
		  // ### Recording the pitch of this calo point along the track in this plane ###
		  trkpitchhit[i][pl][k] = calos[j]->TrkPitchVec()[k];
                  trkxyz[i][pl][k][0] = calos[j]->XYZ()[k].X();
                  trkxyz[i][pl][k][1] = calos[j]->XYZ()[k].Y();
                  trkxyz[i][pl][k][2] = calos[j]->XYZ()[k].Z();
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

      // ----------------------------------------------------------------------------------------------------------------------------
      // ----------------------------------------------------------------------------------------------------------------------------
      //						RECORDING THE TRAJECTORY POINTS FOR THIS TRACK
      // ----------------------------------------------------------------------------------------------------------------------------
      // ---------------------------------------------------------------------------------------------------------------------------- 
   
      nTrajPoint[i] = tracklist[i]->NumberTrajectoryPoints();
   
      // ### Storing the trajectory points in a similar way to PionXS ###
      TVector3 p_hat_0;
      // ##############################################
      // ### Looping over all the trajectory points ###
      // ##############################################
      for( size_t iTrajPt = 0; iTrajPt < tracklist[i]->NumberTrajectoryPoints() ; ++iTrajPt )
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
      
      
	  pHat0_X[i][iTrajPt] = p_hat_0.X();
	  pHat0_Y[i][iTrajPt] = p_hat_0.Y();
	  pHat0_Z[i][iTrajPt] = p_hat_0.Z();
      
	  trjPt_X[i][iTrajPt] = tracklist[i]->LocationAtPoint(iTrajPt).X();
	  trjPt_Y[i][iTrajPt] = tracklist[i]->LocationAtPoint(iTrajPt).Y();
	  trjPt_Z[i][iTrajPt] = tracklist[i]->LocationAtPoint(iTrajPt).Z();

	}//<---End iTrajPt

      if (fmthm.isValid()){
	auto vhit = fmthm.at(i);
	auto vmeta = fmthm.data(i);
	for (size_t h = 0; h < vhit.size(); ++h){
	  if (vhit[h].key()<kMaxHits){
	    //hit_trkkey[vhit[h].key()] = tracklist[i].key();
	    if (vmeta[h]->Dx()){
	      hit_dQds[vhit[h].key()] = vhit[h]->Integral()*fCalorimetryAlg.LifetimeCorrection(vhit[h]->PeakTime())/vmeta[h]->Dx();
	      hit_dEds[vhit[h].key()] = fCalorimetryAlg.dEdx_AREA(vhit[h], vmeta[h]->Dx());
	    }
	    hit_ds[vhit[h].key()] = vmeta[h]->Dx();
	    hit_resrange[vhit[h].key()] = tracklist[i]->Length(vmeta[h]->Index());
	    hit_x[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).X();
	    hit_y[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Y();
	    hit_z[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Z();
	  }
	}//loop over all hits
      }//fmthm is valid   

      if (!isdata&&fmth.isValid()){
	// Find true track for each reconstructed track
	int TrackID = 0;
	std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);
      
	std::map<int,double> trkide;
	for(size_t h = 0; h < allHits.size(); ++h){
	  art::Ptr<recob::Hit> hit = allHits[h];
	  std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
	  for(size_t e = 0; e < TrackIDs.size(); ++e){
	    trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
	  }	    
	}
	// Work out which IDE despoited the most charge in the hit if there was more than one.
	double maxe = -1;
	double tote = 0;
	for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	  tote += ii->second;
	  if ((ii->second)>maxe){
	    maxe = ii->second;
	    TrackID = ii->first;
	  }
	  if ((ii->first) == 1){
	    if ((ii->second)>maxtrackenergy){
	      maxtrackenergy = ii->second;
	      primarytrkkey = tracklist[i].key();
	    }
	  }
	}
	// Now have trackID, so get PdG code and T0 etc.
	const simb::MCParticle *particle = bt->TrackIDToParticle(TrackID);
	if (particle){
	  trkg4id[i] = TrackID;
	}
      }   

    }//<---End track loop (i)

  nhits = hitlist.size();
  for (size_t i = 0; i<hitlist.size() && int(i)< kMaxHits; ++i){
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
	    const auto & tdcidemap = chan->TDCIDEMap();
	    for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	      // loop over the vector of IDE objects.
	      
	      const std::vector<sim::IDE> & idevec = (*mapitr).second;
	      
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
  //fTree->Branch("trigtime",trigtime,"trigtime[16]/I");
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
  fTree->Branch("trkWCtoTPCMath",trkWCtoTPCMath,"trkWCtoTPCMath/I");
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
  fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][2][1000]/D");
  fTree->Branch("trkrr",trkrr,"trkrr[ntracks_reco][2][1000]/D");
  fTree->Branch("trkpitchhit",trkpitchhit,"trkpitchhit[ntracks_reco][2][1000]/D");
  fTree->Branch("trkxyz",trkxyz,"trkxyz[ntracks_reco][2][1000][3]/D");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco][2]/D");
  fTree->Branch("trkpida",trkpida,"trkpida[ntracks_reco][2]/D");
  
  fTree->Branch("nTrajPoint", &nTrajPoint, "nTrajPoint[ntracks_reco]/I");
  fTree->Branch("pHat0_X", pHat0_X, "pHat0_X[ntracks_reco][1000]/D");
  fTree->Branch("pHat0_Y", pHat0_Y, "pHat0_Y[ntracks_reco][1000]/D");
  fTree->Branch("pHat0_Z", pHat0_Z, "pHat0_Z[ntracks_reco][1000]/D");
  fTree->Branch("trjPt_X", trjPt_X, "trjPt_X[ntracks_reco][1000]/D");
  fTree->Branch("trjPt_Y", trjPt_Y, "trjPt_Y[ntracks_reco][1000]/D");
  fTree->Branch("trjPt_Z", trjPt_Z, "trjPt_Z[ntracks_reco][1000]/D");
  fTree->Branch("trkg4id", trkg4id, "trkg4id[ntracks_reco]/I");
  fTree->Branch("primarytrkkey", primarytrkkey, "primarytrkkey/I");
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
  fTree->Branch("hit_dQds", hit_dQds, "hit_dQds[nhits]/F");
  fTree->Branch("hit_dEds", hit_dEds, "hit_dEds[nhits]/F");
  fTree->Branch("hit_ds", hit_ds, "hit_ds[nhits]/F");
  fTree->Branch("hit_resrange", hit_resrange, "hit_resrange[nhits]/F");
  fTree->Branch("hit_x", hit_x, "hit_x[nhits]/F");
  fTree->Branch("hit_y", hit_y, "hit_y[nhits]/F");
  fTree->Branch("hit_z", hit_z, "hit_z[nhits]/F");

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
  fTree->Branch("WC1xPos",WC1xPos,"WC1xPos[nwctrks]/F");
  fTree->Branch("WC1yPos",WC1yPos,"WC1yPos[nwctrks]/F");
  fTree->Branch("WC1zPos",WC1zPos,"WC1zPos[nwctrks]/F");
  fTree->Branch("WC2xPos",WC2xPos,"WC2xPos[nwctrks]/F");
  fTree->Branch("WC2yPos",WC2yPos,"WC2yPos[nwctrks]/F");
  fTree->Branch("WC2zPos",WC2zPos,"WC2zPos[nwctrks]/F");
  fTree->Branch("WC3xPos",WC3xPos,"WC3xPos[nwctrks]/F");
  fTree->Branch("WC3yPos",WC3yPos,"WC3yPos[nwctrks]/F");
  fTree->Branch("WC3zPos",WC3zPos,"WC3zPos[nwctrks]/F");
  fTree->Branch("WC4xPos",WC4xPos,"WC4xPos[nwctrks]/F");
  fTree->Branch("WC4yPos",WC4yPos,"WC4yPos[nwctrks]/F");
  fTree->Branch("WC4zPos",WC4zPos,"WC4zPos[nwctrks]/F");
  
  fTree->Branch("ntof", &ntof, "ntof/I");
  fTree->Branch("tofObject", tofObject, "tofObject[ntof]/D");
  fTree->Branch("tof_timestamp", tof_timestamp, "tof_timestamp[ntof]/D"); 

  fTree->Branch("nAG", &nAG, "nAG/I");
  fTree->Branch("HitTimeStamp1p10_1", HitTimeStamp1p10_1, "HitTimeStamp1p10_1[nAG]/D");
  fTree->Branch("HitTimeStamp1p10_2", HitTimeStamp1p10_2, "HitTimeStamp1p10_2[nAG]/D");
  fTree->Branch("HitTimeStamp1p06_1", HitTimeStamp1p06_1, "HitTimeStamp1p06_1[nAG]/D");
  fTree->Branch("HitTimeStamp1p06_2", HitTimeStamp1p06_2, "HitTimeStamp1p06_2[nAG]/D");
  fTree->Branch("HitPulseArea1p10_1", HitPulseArea1p10_1, "HitPulseArea1p10_1[nAG]/F");
  fTree->Branch("HitPulseArea1p10_2", HitPulseArea1p10_2, "HitPulseArea1p10_2[nAG]/F");
  fTree->Branch("HitPulseArea1p06_1", HitPulseArea1p06_1, "HitPulseArea1p06_1[nAG]/F");
  fTree->Branch("HitPulseArea1p06_2", HitPulseArea1p06_2, "HitPulseArea1p06_2[nAG]/F");
  fTree->Branch("HitExist1p10_1", HitExist1p10_1, "HitExist1p10_1[nAG]/O");
  fTree->Branch("HitExist1p10_2", HitExist1p10_2, "HitExist1p10_2[nAG]/O");
  fTree->Branch("HitExist1p06_1", HitExist1p06_1, "HitExist1p06_1[nAG]/O");
  fTree->Branch("HitExist1p06_2", HitExist1p06_2, "HitExist1p06_2[nAG]/O");

  fTree->Branch("maxTrackIDE", &maxTrackIDE, "maxTrackIDE/I");
  fTree->Branch("IDEEnergy", IDEEnergy, "IDEEnergy[maxTrackIDE]/D");
  fTree->Branch("IDEPos", IDEPos, "IDEPos[maxTrackIDE][3]/D");

  fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
  
  fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
  fTree->Branch("Eng",Eng,"Eng[geant_list_size]/D");
  fTree->Branch("Px",Px,"Px[geant_list_size]/D");
  fTree->Branch("Py",Py,"Py[geant_list_size]/D");
  fTree->Branch("Pz",Pz,"Pz[geant_list_size]/D");
  fTree->Branch("EndEng",EndEng,"EndEng[geant_list_size]/D");
  fTree->Branch("EndPx",EndPx,"EndPx[geant_list_size]/D");
  fTree->Branch("EndPy",EndPy,"EndPy[geant_list_size]/D");
  fTree->Branch("EndPz",EndPz,"EndPz[geant_list_size]/D");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/D");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/D");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/D");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/D");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/D");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/D");
  fTree->Branch("Process", Process, "Process[geant_list_size]/I");
  fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
  fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
  fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
  fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");
  fTree->Branch("G4Process",&G4Process);//,"G4Process[geant_list_size]");
  fTree->Branch("G4FinalProcess",&G4FinalProcess);//,"G4FinalProcess[geant_list_size]");  
  fTree->Branch("NTrTrajPts",NTrTrajPts,"NTrTrajPts[no_primaries]/I");
  fTree->Branch("MidPosX",MidPosX,"MidPosX[no_primaries][5000]/D");
  fTree->Branch("MidPosY",MidPosY,"MidPosY[no_primaries][5000]/D");
  fTree->Branch("MidPosZ",MidPosZ,"MidPosZ[no_primaries][5000]/D");
  fTree->Branch("MidPx",MidPx,"MidPx[no_primaries][5000]/D");
  fTree->Branch("MidPy",MidPy,"MidPy[no_primaries][5000]/D");
  fTree->Branch("MidPz",MidPz,"MidPz[no_primaries][5000]/D");
  fTree->Branch("InteractionPoint"         ,&InteractionPoint         );
  fTree->Branch("InteractionPointType"     ,&InteractionPointType     );


  fTree->Branch("no_mcshowers", &no_mcshowers, "no_mcshowers/I");
  fTree->Branch("mcshwr_origin", mcshwr_origin, "mcshwr_origin[no_mcshowers]/D");
  fTree->Branch("mcshwr_pdg", mcshwr_pdg, "mcshwr_pdg[no_mcshowers]/D");
  fTree->Branch("mcshwr_TrackId", mcshwr_TrackId, "mcshwr_TrackId[no_mcshowers]/I");
  fTree->Branch("mcshwr_startX", mcshwr_startX, "mcshwr_startX[no_mcshowers]/D");
  fTree->Branch("mcshwr_startY", mcshwr_startY, "mcshwr_startY[no_mcshowers]/D");
  fTree->Branch("mcshwr_startZ", mcshwr_startZ, "mcshwr_startZ[no_mcshowers]/D");
  fTree->Branch("mcshwr_endX", mcshwr_endX, "mcshwr_endX[no_mcshowers]/D");
  fTree->Branch("mcshwr_endY", mcshwr_endY, "mcshwr_endY[no_mcshowers]/D");
  fTree->Branch("mcshwr_endZ", mcshwr_endZ, "mcshwr_endZ[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngX", mcshwr_CombEngX, "mcshwr_CombEngX[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngY", mcshwr_CombEngY, "mcshwr_CombEngY[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngZ", mcshwr_CombEngZ, "mcshwr_CombEngZ[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngPx", mcshwr_CombEngPx, "mcshwr_CombEngPx[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngPy", mcshwr_CombEngPy, "mcshwr_CombEngPy[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngPz", mcshwr_CombEngPz, "mcshwr_CombEngPz[no_mcshowers]/D");
  fTree->Branch("mcshwr_CombEngE", mcshwr_CombEngE, "mcshwr_CombEngE[no_mcshowers]/D");
  fTree->Branch("mcshwr_dEdx", mcshwr_dEdx, "mcshwr_dEdx[no_mcshowers]/D");
  fTree->Branch("mcshwr_StartDirX", mcshwr_StartDirX, "mcshwr_StartDirX[no_mcshowers]/D");
  fTree->Branch("mcshwr_StartDirY", mcshwr_StartDirY, "mcshwr_StartDirY[no_mcshowers]/D");
  fTree->Branch("mcshwr_StartDirZ", mcshwr_StartDirZ, "mcshwr_StartDirZ[no_mcshowers]/D");
  fTree->Branch("mcshwr_isEngDeposited", mcshwr_isEngDeposited, "mcshwr_isEngDeposited[no_mcshowers]/I");
  fTree->Branch("mcshwr_Motherpdg", mcshwr_Motherpdg, "mcshwr_Motherpdg[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherTrkId", mcshwr_MotherTrkId, "mcshwr_MotherTrkId[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherstartX", mcshwr_MotherstartX, "mcshwr_MotherstartX[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherstartY", mcshwr_MotherstartY, "mcshwr_MotherstartY[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherstartZ", mcshwr_MotherstartZ, "mcshwr_MotherstartZ[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherendX", mcshwr_MotherendX, "mcshwr_MotherendX[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherendY", mcshwr_MotherendY, "mcshwr_MotherendY[no_mcshowers]/I");
  fTree->Branch("mcshwr_MotherendZ", mcshwr_MotherendZ, "mcshwr_MotherendZ[no_mcshowers]/I");
  
  fTree->Branch("mcshwr_Ancestorpdg", mcshwr_Ancestorpdg, "mcshwr_Ancestorpdg[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorTrkId", mcshwr_AncestorTrkId, "mcshwr_AncestorTrkId[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorstartX", mcshwr_AncestorstartX, "mcshwr_AncestorstartX[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorstartY", mcshwr_AncestorstartY, "mcshwr_AncestorstartY[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorstartZ", mcshwr_AncestorstartZ, "mcshwr_AncestorstartZ[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorendX", mcshwr_AncestorendX, "mcshwr_AncestorendX[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorendY", mcshwr_AncestorendY, "mcshwr_AncestorendY[no_mcshowers]/I");
  fTree->Branch("mcshwr_AncestorendZ", mcshwr_AncestorendZ, "mcshwr_AncestorendZ[no_mcshowers]/I");
      
  fTree->Branch("nshowers",&nshowers,"nshowers/I");
  fTree->Branch("shwID",shwID,"shwI[nshowers]/I");
  fTree->Branch("BestPlaneShw",BestPlaneShw,"BestPlaneShw[nshowers]/I");
  fTree->Branch("LengthShw",LengthShw,"LengthShw[nshowers]/D");
  fTree->Branch("CosStartShw",CosStartShw,"CosStartShw[3][1000]/D");
  // fTree->Branch("CosStartSigmaShw",CosStartSigmaShw,"CosStartSigmaShw[3][nshowers]/D");
  fTree->Branch("CosStartXYZShw",CosStartXYZShw,"CosStartXYZShw[3][1000]/D");
  //fTree->Branch("CosStartXYZSigmaShw",CosStartXYZSigmaShw,"CosStartXYZSigmaShw[3][nshowers]/D");
  fTree->Branch("TotalEShw",TotalEShw,"TotalEShw[2][1000]/D");
  //fTree->Branch("TotalESigmaShw",TotalESigmaShw,"TotalESigmaShw[2][nshowers]/D");
  fTree->Branch("dEdxPerPlaneShw",dEdxPerPlaneShw,"dEdxPerPlaneShw[2][1000]/D");
  //fTree->Branch("dEdxSigmaPerPlaneShw",dEdxSigmaPerPlaneShw,"dEdxSigmaPerPlaneShw[2][nshowers]/D");
  fTree->Branch("TotalMIPEShw",TotalMIPEShw,"TotalMIPEShw[2][1000]/D");
  //fTree->Branch("TotalMIPESigmaShw",TotalMIPESigmaShw,"TotalMIPESigmaShw[2][nshowers]/D");

}

void lariat::AnaTreeT1034::ResetVars()
{
  G4Process.clear();
  G4FinalProcess.clear();

  InteractionPoint.clear();
  InteractionPointType.clear();

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  //  for (int i = 0; i < 16; ++i){
  //     trigtime[i]=-99999;
  //  }
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
    trkWCtoTPCMath = -99999;
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
      pHat0_X[i][j] = -99999;
      pHat0_Y[i][j] = -99999;
      pHat0_Z[i][j] = -99999;
      trjPt_X[i][j] = -99999;
      trjPt_Y[i][j] = -99999;
      trjPt_Z[i][j] = -99999;
    }
    for (int j = 0; j<2; ++j){
      trkpitch[i][j] = -99999;
      trkhits[i][j] = -99999; 
      trkke[i][j] = -99999;
      trkpida[i][j] = -99999;
      for (int k = 0; k<1000; ++k){
			trkdedx[i][j][k] = -99999;
			trkdqdx[i][j][k] = -99999;
			trkrr[i][j][k] = -99999;
			trkpitchhit[i][j][k] = -99999;
			trkxyz[i][j][k][0] = -99999;
			trkxyz[i][j][k][1] = -99999;
			trkxyz[i][j][k][2] = -99999;
      }
    }
    trkg4id[i] = -9999;
  }
  primarytrkkey = -9999;
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
    hit_dQds[i] = -99999;
    hit_dEds[i] = -99999;
    hit_ds[i] = -99999;
    hit_resrange[i] = -99999;
    hit_x[i] = -99999;
    hit_y[i] = -99999;
    hit_z[i] = -99999;
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
      WC1xPos[i] = -99999;
      WC1yPos[i] = -99999;
      WC1zPos[i] = -99999;
      WC2xPos[i] = -99999;
      WC2yPos[i] = -99999;
      WC2zPos[i] = -99999;
      WC3xPos[i] = -99999;
      WC3yPos[i] = -99999;
      WC3zPos[i] = -99999;
      WC4xPos[i] = -99999;
      WC4yPos[i] = -99999;
      WC4zPos[i] = -99999;
   
    }//<---End I loop
  
  ntof = -99999;
  for (int i = 0; i < kMaxTOF; i++)
    {
      tofObject[i] = -99999;
      tof_timestamp[i] = -99999;
	
	
    }//<---End i loop
  
  nAG = -99999;
  for (int i = 0; i < kMaxAG; i++)
    {
      HitTimeStamp1p10_1[i] = -99999;
      HitTimeStamp1p10_2[i] = -99999;
      HitTimeStamp1p06_1[i] = -99999;
      HitTimeStamp1p06_2[i] = -99999;

      HitPulseArea1p10_1[i] = -99999;
      HitPulseArea1p10_2[i] = -99999;
      HitPulseArea1p06_1[i] = -99999;
      HitPulseArea1p06_2[i] = -99999;

      HitExist1p10_1[i] = -99999;
      HitExist1p10_2[i] = -99999;
      HitExist1p06_1[i] = -99999;
      HitExist1p06_2[i] = -99999;

    }//<---End i loop

  maxTrackIDE = -999;
  
  for(size_t i = 0; i < kMaxIDE; ++i) {
    IDEEnergy[i] = -999;

    IDEPos[i][0] = -999.9;
    IDEPos[i][1] = -999.9;
    IDEPos[i][2] = -999.9;

  } // End of maxTrackID loop

  no_primaries = -99999;
  geant_list_size=-999;
  for (int i = 0; i<kMaxPrimaries; ++i){
    pdg[i] = -99999;
    Eng[i] = -99999;
    Px[i] = -99999;
    Py[i] = -99999;
    Pz[i] = -99999;
    EndEng[i] = -99999;
    EndPx[i] = -99999;
    EndPy[i] = -99999;
    EndPz[i] = -99999;
    StartPointx[i] = -99999;
    StartPointy[i] = -99999;
    StartPointz[i] = -99999;
    EndPointx[i] = -99999;
    EndPointy[i] = -99999;
    EndPointz[i] = -99999;
    Process[i] = -99999;
    NumberDaughters[i] = -99999;
    Mother[i] = -99999;
    TrackId[i] = -99999;
    process_primary[i] = -99999;

  }

	for(int i = 0; i<kMaxPrimaryPart; i++){
		NTrTrajPts[i] = -99999;
		for(int j = 0; j<kMaxTruePrimaryPts; j++){	
			MidPosX[i][j] = -99999;
			MidPosY[i][j] = -99999;
			MidPosZ[i][j] = -99999;
   		MidPx[i][j] = -99999; 
   		MidPy[i][j] = -99999; 
   		MidPz[i][j] = -99999;
		}
	}

  nshowers = -99999;
  for (int i = 0; i<kMaxShower; ++i) 
    {
      shwID[i] = -99999;
      BestPlaneShw[i] = -99999;
      LengthShw[i] = -99999;
      for (int j = 0; j<3; ++j) 
	{
	  CosStartShw[j][i] = -99999;
	  //CosStartSigmaShw[j][i] = -99999;
	  CosStartXYZShw[j][i] = -99999;
	  // 	 CosStartXYZSigmaShw[j][i] = -99999;
	  // CosStartXYZSigmaShw[j][i] = -99999;
	}
      for (int j = 0; j<2; ++j) 
	{
	  TotalEShw[j][i] = -99999;
	  //TotalESigmaShw[j][i] = -99999;
	  //TotalESigmaShw[j][i] = -99999;
	  dEdxPerPlaneShw[j][i] = -99999;
	  //dEdxSigmaPerPlaneShw[j][i] = -99999;
	  TotalMIPEShw[j][i] = -99999;
	  //TotalMIPESigmaShw[j][i] = -99999;
	}
    }
  
}

DEFINE_ART_MODULE(lariat::AnaTreeT1034)
