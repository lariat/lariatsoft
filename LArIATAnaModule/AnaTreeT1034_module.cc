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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
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
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"

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
#include "LArIATRecoAlg/BeamlineMassAlg.h"
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

const int kMaxTrack      = 10000;     //maximum number of tracks
const int kMaxHits       = 20000;    //maximum number of hits
const int kMaxTrackHits  = 1000;      //maximum number of space points
const int kMaxTrajHits   = 1000;     //maximum number of trajectory points
const int kMaxCluster    = 1000;     //maximum number of clusters
const int kMaxWCTracks   = 1000;     //maximum number of wire chamber tracks
const int kMaxTOF        = 100;      //maximum number of TOF objects
const int kMaxAG         = 100;      //maximum number of AG objects
const int kMaxPrimaryPart = 50;     //maximum number of true primary particles
const int kMaxPrimaries  = 20000;    //maximum number of true particles tracked
const int kMaxShower     = 100;      //maximum number of Reconstructed showers
const int kMaxMCShower   = 1000;     //maximum number of MCShower Object
const int kMaxTruePrimaryPts = 5000; //maximum number of points in the true primary trajectory
const int kMaxIDE        = 5000; //maximum number of points in the true primary trajectory

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
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p);

private:
  // === Function used to reset all the variables  ===
  void ResetVars();
 
  // === Get path length of MCParticle in TPC ===
  float TrajLengthInTpcAV(const simb::MCParticle*);
  bool  IsPointInTpcAV(const simb::MCParticle*, int);
  bool  IsPointInTpcAV(TVector3&);

  // === Storing information into TTree ====
  TTree* fTree;
  
  // === Enable saving different types of data into tree ===
  bool  fSaveBeamlineInfo;
  bool  fSaveWireChamberHits;
  bool  fSaveGeantInfo;
  bool  fSaveGeantTrajectories;
  bool  fSaveMCShowerInfo;
  bool  fSaveAerogelInfo;
  bool  fSaveSimChannelInfo;
  bool  fSaveTrack3DSpacePoints;
  bool  fSaveTrackCalorimetry;
  bool  fSaveTrackTrajectories;

  // === Select mass range to save to tree
  std::vector<float> fSelectBeamlineMassRange;

  //=== Storing Run Information ===
  int run;			          //<---Run Number
  int subrun;			        //<---SubRun Number
  int event;			        //<---Event Number
  double evttime;		      //<---Event Time (UNIX time)
  //float trig_timestamp;   //<---Trigger timestmap within spill (0-60sec)
  float efield;		        //<---Electric Field
  float lifetime;         //<---Electron lifetime
  int t0;
  //int trigtime[16];		//<---Trigger time

  // === Storing Track Information ===
  int   ntracks_reco;		//<---Number of reconstructed tracks
  float trkvtxx[kMaxTrack];	//<---Track Start X Position (in cm)
  float trkvtxy[kMaxTrack];	//<---Track Start Y Position (in cm)
  float trkvtxz[kMaxTrack];	//<---Track Start Z Position (in cm)
  float trkendx[kMaxTrack];	//<---Track End X Position (in cm)
  float trkendy[kMaxTrack];	//<---Track End Y Position (in cm)
  float trkendz[kMaxTrack];	//<---Track End Z Position (in cm)
  float trkstartdcosx[kMaxTrack];	//<---Direction of the track in the x coordinate at its start point
  float trkstartdcosy[kMaxTrack];	//<---Direction of the track in the y coordinate at its start point
  float trkstartdcosz[kMaxTrack];	//<---Direction of the track in the z coordinate at its start point
  float trkenddcosx[kMaxTrack];	//<---Direction of the track in the x coordinate at its end point
  float trkenddcosy[kMaxTrack];	//<---Direction of the track in the y coordinate at its end point
  float trkenddcosz[kMaxTrack];	//<---Direction of the track in the z coordinate at its end point
  float trklength[kMaxTrack];		//<---Calculated length of the track
//  double trkmomrange[kMaxTrack];	//<---Calculated track momentum from its length assuming a PID of 13
//  double trkmommschi2[kMaxTrack];	//<---Calculated track momentum from multiple scattering using Chi2
//  double trkmommsllhd[kMaxTrack];	//<---Calculated track momentum from multiple scattering
  int    trkWCtoTPCMatch[kMaxTrack];	//<---Using an association to see if there was a match between WC and TPC
  					//    1 = match, 0 = no match
  
  // === Storing the tracks SpacePoints (Individual 3D points)
  int ntrkpts[kMaxTrack];			          //<---Number of SpacePoints associated with track
  float trkx[kMaxTrack][kMaxTrackHits];	//<---X position of the spacepoint
  float trky[kMaxTrack][kMaxTrackHits];	//<---Y position of the spacepoint
  float trkz[kMaxTrack][kMaxTrackHits];	//<---Z position of the spacepoint

  // === Storing the tracks Calorimetry Information
  int ntrkcalopts[kMaxTrack][2];      //<--- Number of calorimetry points (per plane) for dEdx, dQdx, rr
  float trkpida[kMaxTrack][2];        //<--- Track PIDA score
  float trkke[kMaxTrack][2];          //<--- Track kinetic energy deposited
  float trkdedx[kMaxTrack][2][kMaxTrackHits];  //<--- Track dE/dx
  float trkdqdx[kMaxTrack][2][kMaxTrackHits];  //<--- Track dQ/dx
  float trkrr[kMaxTrack][2][kMaxTrackHits];    //<--- Track residual range
  float trkpitch[kMaxTrack][2][kMaxTrackHits]; //<--- Track pitch
  float trkxyz[kMaxTrack][2][kMaxTrackHits][3];//<--- Track XYZ information from calorimetry

  // === Storing trajectory information for the track ===
  int nTrajPoint[kMaxTrack];			        //<---Storing the number of trajectory points
  float pHat0_X[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the x-dir
  float pHat0_Y[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the y-dir
  float pHat0_Z[kMaxTrack][kMaxTrajHits];	//<---Storing trajectory point in the z-dir
  float trjPt_X[kMaxTrack][kMaxTrajHits]; //<---Storing the trajector point location in X
  float trjPt_Y[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Y
  float trjPt_Z[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Z



  // === Geant information for reconstruction track
  int trkg4id[kMaxHits];         //<---geant track id for the track

  int primarytrkkey;             //<---reco track index for primary particle
  // === Storing 2-d Hit information ===
  int    nhits;		//<---Number of 2-d hits in the event
  int    hit_plane[kMaxHits];	//<---Plane number of the hit
  int    hit_wire[kMaxHits];	//<---Wire number of the hit
  int    hit_channel[kMaxHits];//<---Channel number of the hit
  int    hit_trkid[kMaxHits];	//<---The track ID associated with this hit
  int    hit_clusterid[kMaxHits];  //<---cluster ID associated with hit (NEED TO GET WORKING!)
  float  hit_peakT[kMaxHits];	//<---Peak Time of the hit (in drift ticks)
  float  hit_charge[kMaxHits];	//<---Number of ADC's assoicated with the hit (not in units of actual charge)
  float  hit_electrons[kMaxHits];	//<---Number of electrons in hit using CalAreaConstants
  float  hit_ph[kMaxHits];	//<---Amplitude of the hit (usually the gaussian associated with the hit)
  float  hit_rms[kMaxHits]; // <---Hit width (RMS) in time-ticks
  float  hit_tstart[kMaxHits];	//<---Start time of the hit (in drift ticks)
  float  hit_tend[kMaxHits];	//<---End time of the hit (in drift ticks)
  float  hit_driftT[kMaxHits];  //<--- Hit peak time, corrected for trigger/readout offsets
  float  hit_dQds[kMaxHits];   //<---hit dQ/ds
  float  hit_dEds[kMaxHits];   //<---hit dE/ds
  float  hit_ds[kMaxHits];     //<---hit ds
  float  hit_resrange[kMaxHits];//<---hit residual range
  float  hit_x[kMaxHits];        //<---hit x coordinate
  float  hit_y[kMaxHits];        //<---hit y coordinate
  float  hit_z[kMaxHits];        //<---hit z coordinate
  int    hit_g4id[kMaxHits];    //<--- G4 Track ID that made this hit
  float  hit_g4frac[kMaxHits];    //<--- Fraction of hit energy contributed by leading MCParticle (g4id)
  float  hit_g4nelec[kMaxHits];   //<--- Number of electrons collected at wire
  float  hit_g4energy[kMaxHits];  //<--- True deposited energy of hit
 
  // === Storing 2-d Cluster Information ===
  int    nclus;
  int   cluplane[kMaxCluster];
  float clustertwire[kMaxCluster];
  float clusterttick[kMaxCluster];
  float cluendwire[kMaxCluster];
  float cluendtick[kMaxCluster];


  // === Storing Wire Chamber Track Information ===
  int nwctrks;                          //<--- Number of wire chamber tracks reconstructed
  float wctrk_XFaceCoor[kMaxWCTracks];  //<---The projected X position of the wctrack at the front face of the TPC
  float wctrk_YFaceCoor[kMaxWCTracks];  //<---The projected Y position of the wctrack at the front face of the TPC
  float wctrk_momentum[kMaxWCTracks];   //<---Reconstructed momentum
  float wctrk_Px[kMaxWCTracks];         //<--- X-component of momentum
  float wctrk_Py[kMaxWCTracks];         //<--- Y-component
  float wctrk_Pz[kMaxWCTracks];         //<--- Z-component
  float wctrk_theta[kMaxWCTracks];      //<---angle of track w.r.t. z axis
  float wctrk_phi[kMaxWCTracks];        //<---angle of track w.r.t. x axis
  float wctrk_residual[kMaxWCTracks];   //How far off a straight line in YZ plane was the the WCTrack?
  int   wctrk_wcmissed[kMaxWCTracks];   //If this was made with 3 points, was it WC2 or WC3 missed?
  int   wctrk_picky[kMaxWCTracks];      //<--- track had *exactly* one X/Y hit in each WC
  float wctrk_XDist[kMaxWCTracks];      //<---X distance between upstream and downstream tracks
  float wctrk_YDist[kMaxWCTracks];      //<---Y distance between upstream and downstream tracks
  float wctrk_ZDist[kMaxWCTracks];      //<---Z distance between upstream and downstream tracks
  float wctrk_YKink[kMaxWCTracks];      //<---angle in Y between upstream and downstream tracks
  int   wctrk_WC1XMult[kMaxWCTracks];   //<--- number of hits found on X-axis wires in WC1
  int   wctrk_WC1YMult[kMaxWCTracks];   //<--- number of hits found on Y-axis wires in WC1
  int   wctrk_WC2XMult[kMaxWCTracks];   //<--- number of hits found on X-axis wires in WC2
  int   wctrk_WC2YMult[kMaxWCTracks];   //<--- number of hits found on Y-axis wires in WC2
  int   wctrk_WC3XMult[kMaxWCTracks];   //<--- number of hits found on X-axis wires in WC3
  int   wctrk_WC3YMult[kMaxWCTracks];   //<--- number of hits found on Y-axis wires in WC3
  int   wctrk_WC4XMult[kMaxWCTracks];   //<--- number of hits found on X-axis wires in WC4
  int   wctrk_WC4YMult[kMaxWCTracks];   //<--- number of hits found on Y-axis wires in WC4
  float XWireHist[kMaxWCTracks][1000];  //<---Coord in terms of wire number
  float YWireHist[kMaxWCTracks][1000];  //<---Coord in terms of wire number
  float XAxisHist[kMaxWCTracks][1000];  //<---coord in terms of units.
  float YAxisHist[kMaxWCTracks][1000];  //<---coord in terms of units.
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
  int   ntof;
  float tof[kMaxTOF];		//<---The TOF calculated (in ns?) for this TOF object
  float tof_timestamp[kMaxTOF];	//<---Time Stamp for this TOF object

  // === Storing calculated beamline mass ===
  float beamline_mass;	//<---Beamline Mass [MeV]

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
  float IDEEnergy[kMaxIDE];
  float IDEPos[kMaxIDE][3];

  // === Storing Geant4 MC Truth Information ===
  int no_primaries;				//<---Number of primary Geant4 particles in the event
  int geant_list_size;				//<---Number of Geant4 particles tracked
  int pdg[kMaxPrimaries];			//<---PDG Code number of this particle
  float Mass[kMaxPrimaries];			//<---Particle mass
  float StartPointx[kMaxPrimaries];		//<---X position that this Geant4 particle started at
  float StartPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle started at
  float StartPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle started at
  float Eng[kMaxPrimaries];			//<---Initial Energy of the particle
  float Px[kMaxPrimaries];			//<---Initial Px momentum of the particle
  float Py[kMaxPrimaries];			//<---Initial Py momentum of the particle
  float Pz[kMaxPrimaries];			//<---Initial Pz momentum of the particle
  float EndPointx[kMaxPrimaries];		//<---X position that this Geant4 particle ended at
  float EndPointy[kMaxPrimaries];		//<---Y position that this Geant4 particle ended at
  float EndPointz[kMaxPrimaries];		//<---Z position that this Geant4 particle ended at
  float EndEng[kMaxPrimaries];			//<---End Energy of the particle
  float EndPx[kMaxPrimaries];			//<---End Px momentum of the particle
  float EndPy[kMaxPrimaries];			//<---End Py momentum of the particle
  float EndPz[kMaxPrimaries];			//<---End Pz momentum of the particle
  float StartT[kMaxPrimaries];			//<---Start time of particle
  float EndT[kMaxPrimaries];                    //<---End time of particle
  float PathLenInTpcAV[kMaxPrimaries];          //<---Particle path length in TPC
  bool  StartInTpcAV[kMaxPrimaries];            //<---Is particle startpoint in TPC active vol?
  bool  EndInTpcAV[kMaxPrimaries];              //<---Is particle endpoint in TPC active vol?
  
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
  std::vector<int>    InteractionPoint;         //<---Geant 4 Primary Trj Point Corresponding to the Interaction
  std::vector<int>    InteractionPointType;     //<---Geant 4 Primary Interaction Type

  // === Storing additionnal Geant4 MC Truth Information for the primary track only ===
  int NTrTrajPts[kMaxPrimaryPart];                    //<--Nb. of true points in the true primary trajectories
  float MidPosX[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--X position of a point in the true primary trajectory
  float MidPosY[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Y position of a point in the true primary trajectory
  float MidPosZ[kMaxPrimaryPart][kMaxTruePrimaryPts];//<--Z position of a point in the true primary trajectory
  float MidPx[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<- Px momentum of a point in the true primary trajectory
  float MidPy[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Py momentum of a point in the true primary trajectory
  float MidPz[kMaxPrimaryPart][kMaxTruePrimaryPts];  //<--Pz momentum of a point in the true primary trajectory

  // ==== Storing MCShower MCTruth Information ===

  int     no_mcshowers;                                 //number of MC Showers in this event.
  int       mcshwr_origin[kMaxMCShower];                //MC Shower origin information.
  int       mcshwr_pdg[kMaxMCShower];                   //MC Shower particle PDG code.
  int       mcshwr_TrackId[kMaxMCShower];               //MC Shower particle G4 track ID.
  float    mcshwr_startX[kMaxMCShower];                //MC Shower particle G4 startX
  float    mcshwr_startY[kMaxMCShower];                //MC Shower particle G4 startY
  float    mcshwr_startZ[kMaxMCShower];                //MC Shower particle G4 startZ
  float    mcshwr_endX[kMaxMCShower];                  //MC Shower particle G4 endX
  float    mcshwr_endY[kMaxMCShower];                  //MC Shower particle G4 endY
  float    mcshwr_endZ[kMaxMCShower];                  //MC Shower particle G4 endZ
  float    mcshwr_CombEngX[kMaxMCShower];              //MC Shower Combined energy deposition information, Start Point X Position.
  float    mcshwr_CombEngY[kMaxMCShower];              //MC Shower Combined energy deposition information, Start Point Y Position.
  float    mcshwr_CombEngZ[kMaxMCShower];              //MC Shower Combined energy deposition information, Start Point Z Position.
  float    mcshwr_CombEngPx[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  float    mcshwr_CombEngPy[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  float    mcshwr_CombEngPz[kMaxMCShower];           //MC Shower Combined energy deposition information, Momentum X direction.
  float    mcshwr_CombEngE[kMaxMCShower];            //MC Shower Combined energy deposition information, Energy
  float    mcshwr_dEdx[kMaxMCShower];                  //MC Shower dEdx, MeV/cm
  float    mcshwr_StartDirX[kMaxMCShower];             //MC Shower Direction of begining of shower, X direction
  float    mcshwr_StartDirY[kMaxMCShower];             //MC Shower Direction of begining of shower, Y direction
  float    mcshwr_StartDirZ[kMaxMCShower];             //MC Shower Direction of begining of shower, Z direction
  int       mcshwr_isEngDeposited[kMaxMCShower];        //tells whether if this shower deposited energy in the detector or not.

                                                        //yes = 1; no =0;
  //MC Shower mother information
  int        mcshwr_Motherpdg[kMaxMCShower];            //MC Shower's mother PDG code.
  int        mcshwr_MotherTrkId[kMaxMCShower];          //MC Shower's mother G4 track ID.
  float     mcshwr_MotherstartX[kMaxMCShower];         //MC Shower's mother  G4 startX .
  float     mcshwr_MotherstartY[kMaxMCShower];         //MC Shower's mother  G4 startY .
  float     mcshwr_MotherstartZ[kMaxMCShower];         //MC Shower's mother  G4 startZ .
  float     mcshwr_MotherendX[kMaxMCShower];          //MC Shower's mother  G4 endX   .
  float     mcshwr_MotherendY[kMaxMCShower];          //MC Shower's mother  G4 endY   .
  float     mcshwr_MotherendZ[kMaxMCShower];          //MC Shower's mother  G4 endZ   .

  //MC Shower ancestor information
  int        mcshwr_Ancestorpdg[kMaxMCShower];          //MC Shower's ancestor PDG code.
  int        mcshwr_AncestorTrkId[kMaxMCShower];        //MC Shower's ancestor G4 track ID.
  float     mcshwr_AncestorstartX[kMaxMCShower];       //MC Shower's ancestor  G4 startX
  float     mcshwr_AncestorstartY[kMaxMCShower];       //MC Shower's ancestor  G4 startY
  float     mcshwr_AncestorstartZ[kMaxMCShower];       //MC Shower's ancestor  G4 startZ
  float     mcshwr_AncestorendX[kMaxMCShower];         //MC Shower's ancestor  G4 endX
  float     mcshwr_AncestorendY[kMaxMCShower];         //MC Shower's ancestor  G4 endY
  float     mcshwr_AncestorendZ[kMaxMCShower];         //MC Shower's ancestor  G4 endZ

  // === Storing Shower Reco Information using ShowerReco3D ===

  int nshowers; ///number of showers per event
  int shwID[kMaxShower];//ID of the reco shower
  float CosStartShw[3][kMaxShower];
  //float CosStartSigmaShw[3][kMaxShower]; // unused
  float CosStartXYZShw[3][kMaxShower];
  //float CosStartXYZSigmaShw[3][kMaxShower]; // unused
  float TotalEShw[2][kMaxShower];/// total energy of the shower (under investigation...)
  //float TotalESigmaShw[2][kMaxShower];// not working
  float dEdxPerPlaneShw[2][kMaxShower];
  //float dEdxSigmaPerPlaneShw[2][kMaxShower];//not working
  float TotalMIPEShw[2][kMaxShower];
  //float TotalMIPESigmaShw[2][kMaxShower];//not working
  int BestPlaneShw[kMaxShower];
  float LengthShw[kMaxShower];



  // ==== NEED TO FIX THESE VARIABLES....FILLED WITH DUMMY VALUES FOR NOW ===


  //std::string fTrigModuleLabel;
  std::string fClusterModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fWCTrackLabel;            // The name of the producer that made tracks through the MWPCs
  std::string fTOFModuleLabel;		// Name of the producer that made the TOF objects
  std::string fAGModuleLabel;         // Name of the producer that made the aerogel objects
  std::string fG4ModuleLabel;
  std::string fShowerModuleLabel;       // Producer that makes showers from clustering
  std::string fMCShowerModuleLabel;	// Producer name that makes MCShower Object
  std::string fSimChanModuleLabel; // Producer that makes SimIDE information
  std::string fWC2TPCModuleLabel;	// Producer which creates an association between WC and TPC Track

  calo::CalorimetryAlg fCaloAlg;
  BeamlineMassAlg      fBeamlineMassAlg;

};


lariat::AnaTreeT1034::AnaTreeT1034(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fCaloAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  , fBeamlineMassAlg(pset.get<fhicl::ParameterSet>("BeamlineMassAlg"))
{
  this->reconfigure(pset);
}

lariat::AnaTreeT1034::~AnaTreeT1034()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034::reconfigure(fhicl::ParameterSet const & pset)
{
  fSaveBeamlineInfo             = pset.get< bool >      ("SaveBeamlineInfo",      true);
  fSaveWireChamberHits          = pset.get< bool >      ("SaveWireChamberHits",   false);
  fSaveGeantInfo                = pset.get< bool >      ("SaveGeantInfo",         true);
  fSaveGeantTrajectories        = pset.get< bool >      ("SaveGeantTrajectories", false);
  fSaveSimChannelInfo           = pset.get< bool >      ("SaveSimChannelInfo",    true);
  fSaveMCShowerInfo             = pset.get< bool >      ("SaveMCShowerInfo",      false);
  fSaveAerogelInfo              = pset.get< bool >      ("SaveAerogelInfo",       false);
  fSaveTrackCalorimetry         = pset.get< bool >      ("SaveTrackCalorimetry",  true);
  fSaveTrackTrajectories        = pset.get< bool >      ("SaveTrackTrajectories", false);
  fSaveTrack3DSpacePoints       = pset.get< bool >      ("SaveTrack3DSpacePoints",false);
  fSelectBeamlineMassRange      = pset.get< std::vector<float> > ("SelectBeamlineMassRange",{-9,-9}); 
  fHitsModuleLabel              = pset.get< std::string >("HitsModuleLabel");
  fTrackModuleLabel             = pset.get< std::string >("TrackModuleLabel");
  fCalorimetryModuleLabel       = pset.get< std::string >("CalorimetryModuleLabel");
  fParticleIDModuleLabel        = pset.get< std::string >("ParticleIDModuleLabel");
  fClusterModuleLabel           = pset.get< std::string >("ClusterModuleLabel");
  fWCTrackLabel                 = pset.get< std::string >("WCTrackLabel");
  fTOFModuleLabel               = pset.get< std::string >("TOFModuleLabel");
  fAGModuleLabel                = pset.get< std::string >("AGModuleLabel");
  fG4ModuleLabel                = pset.get< std::string >("G4ModuleLabel");
  fShowerModuleLabel            = pset.get< std::string >("ShowerModuleLabel");
  fSimChanModuleLabel           = pset.get< std::string >("SimChanModuleLabel");
  fMCShowerModuleLabel          = pset.get< std::string >("MCShowerModuleLabel");
  fWC2TPCModuleLabel            = pset.get< std::string >("WC2TPCModuleLabel"     , "WC2TPCtrk");
  return;
}

void lariat::AnaTreeT1034::analyze(art::Event const & evt)
{
  // #############################################
  // ### Reset variables before we get started ###
  // #############################################
  ResetVars();

  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  //art::ServiceHandle<geo::Geometry> geom;
  // === Liquid Argon Properties Services ===
  //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  // === Detector properties service ===
  auto const* detprop   = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* detclock  = lar::providerFrom<detinfo::DetectorClocksService>();
  // === BackTracker service ===
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // Get run, subrun and event number
  run         = (int)evt.run();
  subrun      = (int)evt.subRun();
  event       = (int)evt.event();
  evttime     = (int)evt.getSubRun().beginTime().value();
  bool isData = (bool)evt.isRealData();
  // Get electric field and trigger offset
  efield      = detprop->Efield(0);
  lifetime    = detprop->ElectronLifetime();
  t0          = detprop->TriggerOffset();
 
  /*
  // Get the timestamp (within the spill cycle) from the opdetpulse
  // objects because I don't know where else this info is saved!
  art::Handle< std::vector< raw::OpDetPulse >> opdetHandle;
  evt.getByLabel("daq", "", opdetHandle);
  trig_timestamp = -1.;
  if( (size_t)opdetHandle->size() > 0 ){
    art::Ptr< raw::OpDetPulse > ThePulsePtr(opdetHandle,0); 
    trig_timestamp = ((float)(*ThePulsePtr).PMTFrame()*8.)/1.0e09;
  }
  */

  std::cout<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<"Run = "<<run<<", SubRun = "<<subrun<<", Evt = "<<event<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<std::endl;

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
  
  // ###################################
  // ### Getting the Hit Information ###
  // ###################################
  
  art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::Hit objects
  std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define tracklist as a pointer to recob::Hits
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}

  // #####################################
  // ### Getting the Track Information ###
  // #####################################
  
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}

  // ##########################################
  // ### Getting the 2D Cluster Information ###
  // ##########################################
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle; //<---Define clusterListHandle as a vector of recob::Cluster objects
  std::vector<art::Ptr<recob::Cluster> > clusterlist; //<---Define cluster as a pointer to recob::Clusters
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
  std::vector<art::Ptr<ldp::TOF> > toflist;
  if(evt.getByLabel(fTOFModuleLabel,TOFColHandle))
    {art::fill_ptr_vector(toflist, TOFColHandle);}

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
  if (evt.getByLabel(fShowerModuleLabel,shwListHandle))
    {art::fill_ptr_vector(shwlist, shwListHandle);}


  // ### Something to do with SimChannels...need to come back to ###
  std::vector<const sim::SimChannel*> fSimChannels;
  try
    {evt.getView("largeant", fSimChannels);}
  catch (art::Exception const&e){ }
  
  
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE Sim Channel INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  if(!isData) {

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
  } // End of isData



  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE MCShower Geant4 INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isData && fSaveMCShowerInfo )
    {
      art::Handle< std::vector<sim::MCShower> > mcshowerh;
      evt.getByLabel(fMCShowerModuleLabel, mcshowerh);

      // #######################################################
      // ### Check to make sure the MCShower Handle is valid ###
      // #######################################################
      if (mcshowerh.isValid())
        {

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
  if(!isData && fSaveGeantInfo )
    {
      // ######################################
      // ### Making a vector of MCParticles ###
      // ######################################
      std::vector<const simb::MCParticle* > geant_part;
    
      // ### Looping over all the Geant4 particles from the BackTracker ###
      for(size_t p = 0; p < plist.size(); ++p) 
        geant_part.push_back(plist.Particle(p)); 
      
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
	 if(geant_part[i]->Process()==pri)  process_primary[i]=1;
	 else                               process_primary[i]=0;
         
         int last         = geant_part[i]->NumberTrajectoryPoints()-1;
          
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
	  
	 if(geant_part[i]->Process() == pri)                Process[i] = 0;
	 if(geant_part[i]->Process() == PionMinusInelastic) Process[i] = 1;
	 if(geant_part[i]->Process() == NeutronInelastic)   Process[i] = 2;
	 if(geant_part[i]->Process() == hadElastic)         Process[i] = 3;
	 if(geant_part[i]->Process() == nCapture)           Process[i] = 4;
	 if(geant_part[i]->Process() == CHIPSNuclearCaptureAtRest)  Process[i] = 5;
	 if(geant_part[i]->Process() == Decay)              Process[i] = 6;
	 if(geant_part[i]->Process() == KaonZeroLInelastic) Process[i] = 7;
	 if(geant_part[i]->Process() == CoulombScat)        Process[i] = 8;
	 if(geant_part[i]->Process() == muMinusCaptureAtRest) Process[i] = 9;
	 if(geant_part[i]->Process() == ProtonInelastic)    Process[i] = 10;
	 if(geant_part[i]->Process() == KaonPlusInelastic)  Process[i] = 11;
	 if(geant_part[i]->Process() == hBertiniCaptureAtRest)  Process[i] = 12;
	     

	 // ### Saving the particles mother, trackID, pdg, and energy ###
	 Mother[i]    =geant_part[i]->Mother();
	 TrackId[i]   =geant_part[i]->TrackId();
	 pdg[i]       =geant_part[i]->PdgCode();
         Mass[i]      =geant_part[i]->Mass();
	 Eng[i]       =geant_part[i]->E();
	 EndEng[i]    =geant_part[i]->EndE();

	 // ### Saving the start and end Px, Py, Pz info ###
	 Px[i]        =geant_part[i]->Px();
	 Py[i]        =geant_part[i]->Py();
	 Pz[i]        =geant_part[i]->Pz();
	 EndPx[i]     =geant_part[i]->EndPx();
	 EndPy[i]     =geant_part[i]->EndPy();
	 EndPz[i]     =geant_part[i]->EndPz();
          
	 // ### Saving the Start and End Point for this particle ###
	 StartPointx[i] =geant_part[i]->Vx();
	 StartPointy[i] =geant_part[i]->Vy();
	 StartPointz[i] =geant_part[i]->Vz();
	 EndPointx[i]   =geant_part[i]->EndPosition()[0];
	 EndPointy[i]   =geant_part[i]->EndPosition()[1];
	 EndPointz[i]   =geant_part[i]->EndPosition()[2];

         // ### Saving the start and end time
         StartT[i]      =geant_part[i]->T();
         EndT[i]        =geant_part[i]->EndT();

         // ### Tpc AV
         PathLenInTpcAV[i]= TrajLengthInTpcAV(geant_part[i]);
         StartInTpcAV[i]  = IsPointInTpcAV(geant_part[i],0);
         EndInTpcAV[i]    = IsPointInTpcAV(geant_part[i],last);

	 // ### Saving the processes for this particle ###
	 G4Process.push_back(       geant_part[i]->Process()    );
	 G4FinalProcess.push_back(  geant_part[i]->EndProcess() );
 	  
	 // ### Saving the number of Daughters for this particle ###
	 NumberDaughters[i]=geant_part[i]->NumberDaughters();
        
         //std::cout<<"Saving particle "<<i<<"   PDG "<<pdg[i]<<" with energy "<<Eng[i]*1e3<<" MeV\n";
        
         // ### Correcting the "end energy"
         //
         // When particles interact in-flight, like in pi- inelastic interactions
         // for example, G4 sets the final kinetic energy to zero at its very
         // last step. This isn't helpful!! Instead, set the end energy to be the 
         // second-to-last trajectory point when this happens.
         int absPdg = fabs(pdg[i]);
         if(    geant_part[i]->NumberTrajectoryPoints() >= 2
            &&  ( absPdg==211 || absPdg == 13 || absPdg == 321 || absPdg == 2212) ){
         
            simb::MCTrajectory traj = geant_part[i]->Trajectory();
            int   Npts   = geant_part[i]->NumberTrajectoryPoints();
            int   last   = Npts-1;
            // if the final KE is ~ 0 (within 0.1 MeV) and the final
            // energy drop is significant (>10 MeV), then correct the final energy
            //float E_end  = traj.E(last);
            //float KE_end = E_end - geant_part[i]->Mass();
            //float dE_end = E_end - traj.E(last-1);
            //if( KE_end <= 0.0001 && fabs(dE_end) > 0.010 ) EndEng[i] = traj.E(last-1);
            EndEng[i] = traj.E(last-1);
         }
          
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
    
    
  if( fSaveBeamlineInfo ) {
    
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							FILLING THE WIRE CHAMBER TRACK INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    nwctrks = wctrack.size();
     
    // ########################################
    // ### Looping over Wire Chamber Tracks ###
    // ########################################
    for(size_t wct_count = 0; wct_count < wctrack.size(); wct_count++)
      //for(const auto& wctrack : (*wctrackHandle)) //trackHandle works somewhat like a pointer to a vector of tracks, so dereference the handle to loop over
      //the vector, then use each "track" as a ldp::WCTrack
      {
        
        
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
        wctrk_YKink[wct_count] = wctrack[wct_count]->YKink();
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
        wctrk_residual[wct_count]= wctrack[wct_count]->Residual();
        wctrk_wcmissed[wct_count]= wctrack[wct_count]->WCMissed();
        wctrk_picky[wct_count]= wctrack[wct_count]->IsPicky();
        wctrk_WC1XMult[wct_count]= wctrack[wct_count]->WC1XMult();
        wctrk_WC1YMult[wct_count]= wctrack[wct_count]->WC1YMult();
        wctrk_WC2XMult[wct_count]= wctrack[wct_count]->WC2XMult();
        wctrk_WC2YMult[wct_count]= wctrack[wct_count]->WC2YMult();
        wctrk_WC3XMult[wct_count]= wctrack[wct_count]->WC3XMult();
        wctrk_WC3YMult[wct_count]= wctrack[wct_count]->WC3YMult();
        wctrk_WC4XMult[wct_count]= wctrack[wct_count]->WC4XMult();
        wctrk_WC4YMult[wct_count]= wctrack[wct_count]->WC4YMult();
        
        // Calculate XYZ components of momentum vector using TVector3
        float p     = wctrack[wct_count]->Momentum();
        float theta = wctrack[wct_count]->Theta();
        float phi   = wctrack[wct_count]->Phi();
        TVector3 pvect(p,theta,phi);
        wctrk_Px[wct_count] = pvect.X();
        wctrk_Py[wct_count] = pvect.Y();
        wctrk_Pz[wct_count] = pvect.Z();
        
        /*
        // Compute wcP, wcPx and wcPy
        if (TMath::Cos(wctrk_theta[wct_count])) 
          { 
            float PT = wctrk_momentum[wct_count];
            float tanThetaCosPhi = TMath::Tan(wctrk_theta[wct_count]) * TMath::Cos(wctrk_phi[wct_count]);
            float tanThetaSinPhi = TMath::Tan(wctrk_theta[wct_count]) * TMath::Sin(wctrk_phi[wct_count]);
            float den = TMath::Sqrt(1+tanThetaCosPhi*tanThetaCosPhi);
            wctrk_momentum_z[wct_count] = PT/den;
            wctrk_momentum_y[wct_count] = wcPz[wct_count]*tanThetaSinPhi;
            wctrk_momentum_x[wct_count] = wcPz[wct_count]*tanThetaCosPhi;
            // Removing wcP -- should be same as "wctrk_momentum" above
            //wcP [wct_count] = wcPz[wct_count]* TMath::Sqrt(1+TMath::Tan(wctrk_theta[wct_count])*TMath::Tan(wctrk_theta[wct_count]));
          }
        */
        
        // === Getting individual channel information ===
        for(size_t chIt = 0; 2*chIt+1 < wctrack[wct_count]->NHits(); ++chIt)
          {
            if(float(chIt)!=wctrack[wct_count]->WCMissed()){  
              XWireHist[wct_count][chIt] = wctrack[wct_count]->HitWire(2*chIt);
              YWireHist[wct_count][chIt] = wctrack[wct_count]->HitWire(2*chIt+1);
              XAxisHist[wct_count][chIt] = wctrack[wct_count]->WC(2*chIt);
              YAxisHist[wct_count][chIt] = wctrack[wct_count]->WC(2*chIt+1);
            }//<-- end if chIt != WCMissed  
          }//<---End chIt loop
        
      }//<---end wctrack auto loop
  	 
     
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							FILLING THE TIME OF FLIGHT INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    ntof = toflist.size();
    
    // ################################
    // ### Looping over TOF objects ###
    // ################################
    size_t tof_counter = 0; // book-keeping
    for(size_t i = 0; i < toflist.size(); i++) {
      
      size_t number_tof = toflist[i]->NTOF();
  
      for (size_t tof_idx = 0; tof_idx < number_tof; ++tof_idx) {
        tof[tof_counter]       =  toflist[i]->SingleTOF(tof_idx);
        tof_timestamp[tof_counter]  = toflist[i]->TimeStamp(tof_idx)/1.0e9; // this is saved in "trig_timestamp" now
        ++tof_counter;
      } // loop over TOF
  
    }//<---End tof_count loop
    
    // Beamline Mass info
    if ( ntof == 1 && nwctrks == 1) 
      beamline_mass = fBeamlineMassAlg.GetMass(tof[0], wctrk_momentum[0]);
    
  }// Beamline Info


  
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //                                                   FILLING THE AEROGEL COUNTER INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  if( fSaveAerogelInfo ) {
     
    // ################################
    // ### Looping over aerogel counter objects ###
    // ################################
    nAG = agc.size();
    size_t agc_counter = 0; // book-keeping
    for(size_t i = 0; i < agc.size(); i++)
      {
        size_t number_agc = agc[i]->GetNHits();
        for (size_t agc_idx = 0; agc_idx < number_agc; ++agc_idx) {
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
  }


  // ##########################################################
  // ### Grabbing associations for use later in the AnaTool ###
  // ##########################################################
  try{
    // === Associations between hits and raw digits ===
    // === Association between SpacePoints and Tracks ===
    art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
    // === Association between Tracks and 2d Hits ===
    //art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
    // === Association between Calorimetry objects and Tracks ===
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    // === Association between Particle ID objects (PID) and Tracks ===
    art::FindManyP<anab::ParticleID>  fmpid(trackListHandle, evt, fParticleIDModuleLabel);
    // ==== Association between Clusters and Hits ===
    art::FindManyP<recob::Cluster>    fmcl(hitListHandle,   evt, fClusterModuleLabel);
    // ==== Association between Tracks and Hit Metadata
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);

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

      for (size_t j = 0; j<3; ++j){
        CosStartShw[j][i] = shwlist[i]->Direction()[j];
        CosStartXYZShw[j][i] = shwlist[i]->ShowerStart()[j];
        // CosStartSigmaShw[j][i] = shwlist[i]->DirectionErr()[j];
        //CosStartXYZSigmaShw[j][i] =  shwlist[i]->ShowerStartErr()[j];
      }

      for (int j = 0; j<2; ++j){ // looping over the 2 planes
        TotalEShw[j][i] = shwlist[i]->Energy()[j];
        //TotalESigmaShw[j][i] = shwlist[i]->EnergyErr()[j];
        dEdxPerPlaneShw[j][i] = shwlist[i]->dEdx()[j];
        TotalMIPEShw[j][i] = shwlist[i]->MIPEnergy()[j];
      }
    }// end loop over showers


    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							FILLING THE 3-D TRACK INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------

    // ### Calling the track momentum calculator ###
    // - change the minimal track length requirement to 10 cm
    trkf::TrackMomentumCalculator trkm{10.0};

    // === Saving the number of tracks per event ===
    ntracks_reco=tracklist.size();

    // === Association between WC Tracks and TPC Tracks ===
    int TempTrackMatchedID = -1;

    if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {
      art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
      if (fWC2TPC.isValid()){
        // === Loop on all the Assn WC-TPC tracks ===
        for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn ){
          // =========================
          // === Get the TPC track ===
          // =========================
          cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));
          if (!trackWC2TPC) continue;
          recob::Track const& aTrack(trackWC2TPC.ref());
          TempTrackMatchedID = aTrack.ID();
        }//<----End indexAssn loop
      }//<---End checking that the WC2TPC
    }


    // ### Looping over tracks ###
    double maxtrackenergy = -1;
    for(size_t i=0; i<tracklist.size();++i){
	
      //-----------------------------------------------------------------------
      // get track information
      //-----------------------------------------------------------------------
      //auto trackStartEnd = tracklist[i]->Extent();
      auto larStartEnd = tracklist[i]->Direction();
      // returns type std::pair<recob::Track::Point_t, recob::Track::Point_t>
      // returns type std::pair<recob::Track::Vector_t, recob::Track::Vector_t>
      
      auto const& startPt = tracklist[i]->Vertex();
      auto const& endPt   = tracklist[i]->End();

      // ### Storing an integer for the match of a WC to TPC track ###
      int trackMatch = 0;
      if(TempTrackMatchedID == tracklist[i]->ID() ){
        trackMatch = 1;
      }//<---End match
      
      // ### Setting the WC to TPC match ###
      trkWCtoTPCMatch[i] = trackMatch;
	
      // ### Recording the track vertex x, y, z location ###
      //trkvtxx[i]        = trackStartEnd.first.X();
      //trkvtxy[i]        = trackStartEnd.first.Y();
      //trkvtxz[i]        = trackStartEnd.first.Z();
      trkvtxx[i]        = startPt.X(); 
      trkvtxy[i]        = startPt.Y();
      trkvtxz[i]        = startPt.Z();
      
      // ### Recording the track end point x, y, z location ###
      //trkendx[i]        = trackStartEnd.second.X();
      //trkendy[i]        = trackStartEnd.second.Y();
      //trkendz[i]        = trackStartEnd.second.Z();
      trkendx[i]        = endPt.X(); 
      trkendy[i]        = endPt.Y();
      trkendz[i]        = endPt.Z();
      
      // ### Recording the directional cosine at the start of the track ###
      trkstartdcosx[i]  = larStartEnd.first.X();
      trkstartdcosy[i]  = larStartEnd.first.Y();
      trkstartdcosz[i]  = larStartEnd.first.Z();
      
      // ### Recording the directional cosine at the end of the track ###
      trkenddcosx[i]    = larStartEnd.second.X();
      trkenddcosy[i]    = larStartEnd.second.Y();
      trkenddcosz[i]    = larStartEnd.second.Z();
      
      // ### Recording the track length as calculated by the tracking module ###
      // ####                (each one may do it differently                 ###
      trklength[i]         = tracklist[i]->Length();
      
      // ### Calculating the track momentum using the momentum calculator ###
      //trkmomrange[i]    = trkm.GetTrackMomentum(trklength[i],13);
      //trkmommschi2[i]   = trkm.GetMomentumMultiScatterChi2(tracklist[i]);
      //trkmommsllhd[i]   = trkm.GetMomentumMultiScatterLLHD(tracklist[i]);
      
      // ### Grabbing the SpacePoints associated with this track ###
      ntrkpts[i] = fmsp.at(i).size();
      std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
      
      // ########################################
      // ### Looping over all the SpacePoints ###
      // ########################################
      for (size_t j = 0; j<spts.size(); ++j){
        // ### Recording the x, y, z location of every spacepoint in the track ###
        trkx[i][j] = spts[j]->XYZ()[0];
        trky[i][j] = spts[j]->XYZ()[1];
        trkz[i][j] = spts[j]->XYZ()[2];
      }//<----End SpacePoint loop (j)

      // ########################################################## 
      // ### Looping over Calorimetry information for the track ###
      // ########################################################## 
      if (fmcal.isValid()){
        // ### Putting calo information for this track (i) into pointer vector ###
        std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(i);
        
        // ### Looping over each calorimetry point (similar to SpacePoint) ###
        for (size_t j = 0; j<calos.size(); ++j){

          // ### If we don't have calorimetry information for this plane skip ###
          if (!calos[j]->PlaneID().isValid) continue;
          
          // ### Grabbing this calorimetry points plane number (0 == induction, 1 == collection) ###
          int pl = calos[j]->PlaneID().Plane;
          
          // ### Skipping this point if the plane number doesn't make sense ###
          if (pl<0||pl>1) continue;
          
          // ### Recording the number of calorimetry points for this track in this plane ####
          ntrkcalopts[i][pl] = calos[j]->dEdx().size();
          
          // #### Recording the kinetic energy for this track in this plane ###
          trkke[i][pl] = calos[j]->KineticEnergy();
          
          // ###############################################
          // ### Looping over all the calorimetry points ###
          // ###############################################
          for (size_t k = 0; k<calos[j]->dEdx().size(); ++k){
            
            // ### If we go over 1000 points just skip them ###
            if (k>=kMaxTrackHits) continue;
		    
            // ### Recording the dE/dX information for this calo point along the track in this plane ###
            trkdedx[i][pl][k] = calos[j]->dEdx()[k];
            trkdqdx[i][pl][k] = calos[j]->dQdx()[k];
            
            // ### Recording the residual range for this calo point along the track in this plane ###
            trkrr[i][pl][k] = calos[j]->ResidualRange()[k];
            
            // ### Recording the pitch of this calo point along the track in this plane ###
            trkpitch[i][pl][k] = calos[j]->TrkPitchVec()[k];
            trkxyz[i][pl][k][0] = calos[j]->XYZ()[k].X();
            trkxyz[i][pl][k][1] = calos[j]->XYZ()[k].Y();
            trkxyz[i][pl][k][2] = calos[j]->XYZ()[k].Z();
          }//<---End calo points (k)
		  
        }//<---End looping over calo objects (j)
	    
      }//<---End checking Calo info is valid 

      // ----------------------------------------------------------------------------------------------------------------------------
      // ----------------------------------------------------------------------------------------------------------------------------
      //						PARTICLE ID INFOR FROM THIS TRACK INFORMATION
      // ----------------------------------------------------------------------------------------------------------------------------
      // ----------------------------------------------------------------------------------------------------------------------------      
      
      // ################################################## 
      // ### Looping over PID information for the track ###
      // ################################################## 
      if (fmpid.isValid()){
        // ### Putting PID information for this track (i) into pointer vector ###
        std::vector<art::Ptr<anab::ParticleID> > pids = fmpid.at(i);
        for (size_t j = 0; j<pids.size(); ++j){
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
        recob::Track::Vector_t tmp = tracklist[i]->DirectionAtPoint(iTrajPt);
        p_hat_0.SetXYZ(tmp.X(),tmp.Y(),tmp.Z());
      
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

      // Find all associated hits and save information 
      // about the track they appear in...
      if (fmthm.isValid()){
        auto& vhit = fmthm.at(i);
        auto& vmeta = fmthm.data(i);
        for (size_t h = 0; h < vhit.size(); ++h){
          if (vhit[h].key()<kMaxHits){
            if (vmeta[h]->Dx()){
              hit_dQds[vhit[h].key()] = vhit[h]->Integral()*fCaloAlg.LifetimeCorrection(vhit[h]->PeakTime())/vmeta[h]->Dx();
              hit_dEds[vhit[h].key()] = fCaloAlg.dEdx_AREA(vhit[h], vmeta[h]->Dx());
            }
            // The hit "key" here doesn't correspond to the actual hitlist indices... 
            // need to loop through and find the match (ugh)
            int hit_index = -9;
            for(size_t k=0; k<hitlist.size(); k++){
              // consider only hits on same plane and wire
              if( hitlist[k]->WireID().Plane != vhit[h]->WireID().Plane ) continue;
              if( hitlist[k]->WireID().Wire  != vhit[h]->WireID().Wire ) continue;
              // if the time matches up exactly, we found our hit
              if( hitlist[k]->PeakTime() == vhit[h]->PeakTime() ){
                hit_index = k;
                break;
              }
            }
            if( hit_index < 0 ) continue;
            hit_trkid[hit_index] = i;
            hit_ds[hit_index] = vmeta[h]->Dx();
            hit_resrange[hit_index] = tracklist[i]->Length(vmeta[h]->Index());
            hit_x[hit_index] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).X();
            hit_y[hit_index] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Y();
            hit_z[hit_index] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Z();
          }
        }//loop over all hits
      }//fmthm is valid   

      // Find all IDEs for this track by going through the associated hits.
      if (!isData && fmth.isValid()){
        std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);
        // find true track for each reconstructed track
        int TrackID = 0;
        std::map<int,double> trkide;
        for(size_t h = 0; h < allHits.size(); ++h){
          art::Ptr<recob::Hit> hit = allHits[h];
          std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
          for(size_t e = 0; e < TrackIDs.size(); ++e) {
            trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
          }
        }
        // work out which IDE despoited the most charge in the hit if there was more than one.
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
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
        if (particle) trkg4id[i] = TrackID;
	    
      }//endif isMC and association is valid   

    }//<---End track loop (i)
    
    
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							SAVING WIRE HIT INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    

    //art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
    nhits = hitlist.size();
    for (size_t i = 0; i<hitlist.size() && int(i)< kMaxHits; ++i){
    // cet::maybe_ref<raw::RawDigit const> rdref(ford.at(i));
      unsigned int channel = hitlist[i]->Channel();
      geo::WireID wireid = hitlist[i]->WireID();
      hit_plane[i]    = wireid.Plane;
      hit_wire[i]     = wireid.Wire;
      hit_channel[i]  = channel;
      hit_peakT[i]    = hitlist[i]->PeakTime();
      hit_driftT[i]   = hitlist[i]->PeakTime()-detprop->GetXTicksOffset(wireid.Plane,0,0);
      hit_charge[i]   = hitlist[i]->Integral();
      hit_electrons[i]= fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),wireid.Plane); 
      hit_ph[i]       = hitlist[i]->PeakAmplitude();
      hit_tstart[i]   = hitlist[i]->StartTick();
      hit_tend[i]     = hitlist[i]->EndTick();
      hit_rms[i]      = hitlist[i]->RMS();
      
      // If this is MC, map hit to its G4 Track ID
      if(!evt.isRealData()){
        std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(hitlist[i]);
        if( trackIDEs.size() ) {
          // Loop through and find the leading TrackIDE, and keep
          // track of the total energy of ALL IDEs.
          float maxe = 0;
          float bestfrac = 0;
          int bestid = 0;
          int ne = 0;
          for(size_t j = 0; j < trackIDEs.size(); ++j){
            ne += trackIDEs[j].numElectrons;
            if( trackIDEs[j].energy > maxe ) {
              maxe = trackIDEs[j].energy;
              bestfrac = trackIDEs[j].energyFrac;
              bestid = trackIDEs[j].trackID;
            }
          }
          hit_g4id[i] = bestid;
          hit_g4frac[i] = bestfrac;
          hit_g4energy[i] = maxe;
          hit_g4nelec[i] = ne;
        }
      }

    }//end loop over hits

     
     

    
  } catch (art::Exception const&e){ }

  // ================================
  // Save to TTree if beamline mass
  // within selection range
  bool flag = true;
  if( fSelectBeamlineMassRange.size()==2 ) {
    float m1  = fSelectBeamlineMassRange[0];
    float m2  = fSelectBeamlineMassRange[1];
    if( m1 > 0 && beamline_mass < m1 ) flag = false;
    if( m2 > 0 && beamline_mass > m2 ) flag = false;
  }
  if( flag ) fTree->Fill();

}


void lariat::AnaTreeT1034::beginJob()
{
  
  

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  //fTree->Branch("trig_timestamp",&trig_timestamp,"trig_timestamp/F");
  fTree->Branch("efield",&efield,"efield/F");
  fTree->Branch("lifetime",&lifetime,"lifetime/F");
  fTree->Branch("t0",&t0,"t0/I");
  fTree->Branch("nclus",&nclus,"nclus/I");
  fTree->Branch("clustertwire",clustertwire,"clustertwire[nclus]/F");
  fTree->Branch("clusterttick",clusterttick,"clusterttick[nclus]/F");
  fTree->Branch("cluendwire",cluendwire,"cluendwire[nclus]/F");
  fTree->Branch("cluendtick",cluendtick,"cluendtick[nclus]/F");
  fTree->Branch("cluplane",cluplane,"cluplane[nclus]/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trkvtxx",trkvtxx,"trkvtxx[ntracks_reco]/F");
  fTree->Branch("trkvtxy",trkvtxy,"trkvtxy[ntracks_reco]/F");
  fTree->Branch("trkvtxz",trkvtxz,"trkvtxz[ntracks_reco]/F");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/F");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/F");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/F");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/F");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/F");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/F");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/F");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/F");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/F");
  fTree->Branch("trkWCtoTPCMatch",trkWCtoTPCMatch,"trkWCtoTPCMatch[ntracks_reco]/I");
  fTree->Branch("trklength",trklength,"trklength[ntracks_reco]/F");
  fTree->Branch("trkg4id", trkg4id, "trkg4id[ntracks_reco]/I");
  fTree->Branch("primarytrkkey", primarytrkkey, "primarytrkkey/I");
//  fTree->Branch("trkmomrange",trkmomrange,"trkmomrange[ntracks_reco]/F");
//  fTree->Branch("trkmommschi2",trkmommschi2,"trkmommschi2[ntracks_reco]/F");
//  fTree->Branch("trkmommsllhd",trkmommsllhd,"trkmommsllhd[ntracks_reco]/F");
  if( fSaveTrackCalorimetry ) {
    fTree->Branch("ntrkcalopts",ntrkcalopts,"ntrkcalopts[ntracks_reco][2]/I");
    fTree->Branch("trkpida",trkpida,"trkpida[ntracks_reco][2]/F");
    fTree->Branch("trkke",trkke,"trkke[ntracks_reco][2]/F");
    fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks_reco][2][1000]/F");
    fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks_reco][2][1000]/F");
    fTree->Branch("trkrr",trkrr,"trkrr[ntracks_reco][2][1000]/F");
    fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][2][1000]/F");
    fTree->Branch("trkxyz",trkxyz,"trkxyz[ntracks_reco][2][1000][3]/F");
  }
  if( fSaveTrack3DSpacePoints ) {
    fTree->Branch("ntrkpts",&ntrkpts,"ntrkpts[ntracks_reco]/I");
    fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/F");
    fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/F");
    fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/F");
  }
  if( fSaveTrackTrajectories ){
    fTree->Branch("nTrajPoint", &nTrajPoint, "nTrajPoint[ntracks_reco]/I");
    fTree->Branch("pHat0_X", pHat0_X, "pHat0_X[ntracks_reco][1000]/F");
    fTree->Branch("pHat0_Y", pHat0_Y, "pHat0_Y[ntracks_reco][1000]/F");
    fTree->Branch("pHat0_Z", pHat0_Z, "pHat0_Z[ntracks_reco][1000]/F");
    fTree->Branch("trjPt_X", trjPt_X, "trjPt_X[ntracks_reco][1000]/F");
    fTree->Branch("trjPt_Y", trjPt_Y, "trjPt_Y[ntracks_reco][1000]/F");
    fTree->Branch("trjPt_Z", trjPt_Z, "trjPt_Z[ntracks_reco][1000]/F");
  }
  
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F");
  fTree->Branch("hit_driftT",hit_driftT,"hit_driftT[nhits]/F");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
  fTree->Branch("hit_electrons",hit_electrons,"hit_electrons[nhits]/F");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/F");
  fTree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F");
  fTree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
  fTree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F");
  fTree->Branch("hit_g4nelec",hit_g4nelec,"hit_g4nelec[nhits]/F");
  fTree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F");
  fTree->Branch("hit_tstart",hit_tstart,"hit_tstart[nhits]/F");
  fTree->Branch("hit_tend",hit_tend,"hit_tend[nhits]/F");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");
  //fTree->Branch("hit_clusterid",hit_clusterid,"hit_clusterid[nhits]/I"); 
  //fTree->Branch("hit_pk",hit_pk,"hit_pk[nhits]/I");
  //fTree->Branch("hit_t",hit_t,"hit_t[nhits]/I");
  //fTree->Branch("hit_ch",hit_ch,"hit_ch[nhits]/I");
  //fTree->Branch("hit_fwhh",hit_fwhh,"hit_fwhh[nhits]/I");
  fTree->Branch("hit_dQds", hit_dQds, "hit_dQds[nhits]/F");
  fTree->Branch("hit_dEds", hit_dEds, "hit_dEds[nhits]/F");
  fTree->Branch("hit_ds", hit_ds, "hit_ds[nhits]/F");
  fTree->Branch("hit_resrange", hit_resrange, "hit_resrange[nhits]/F");
  fTree->Branch("hit_x", hit_x, "hit_x[nhits]/F");
  fTree->Branch("hit_y", hit_y, "hit_y[nhits]/F");
  fTree->Branch("hit_z", hit_z, "hit_z[nhits]/F");
 
  if( fSaveBeamlineInfo ) {
    fTree->Branch("beamline_mass", &beamline_mass, "beamline_mass/F");
    fTree->Branch("nwctrks",&nwctrks,"nwctrks/I");
    fTree->Branch("wctrk_XFaceCoor",wctrk_XFaceCoor,"wctrk_XFaceCoor[nwctrks]/F");
    fTree->Branch("wctrk_YFaceCoor",wctrk_YFaceCoor,"wctrk_YFaceCoor[nwctrks]/F");
    fTree->Branch("wctrk_theta",wctrk_theta,"wctrk_theta[nwctrks]/F");
    fTree->Branch("wctrk_phi",wctrk_phi,"wctrk_phi[nwctrks]/F");
    fTree->Branch("wctrk_momentum",wctrk_momentum,"wctrk_momentum[nwctrks]/F");
    fTree->Branch("wctrk_Px" ,wctrk_Px,"wctrk_Px[nwctrks]/F");
    fTree->Branch("wctrk_Py" ,wctrk_Py,"wctrk_Py[nwctrks]/F");
    fTree->Branch("wctrk_Pz" ,wctrk_Pz,"wctrk_Pz[nwctrks]/F");
    fTree->Branch("wctrk_residual", wctrk_residual,"wctrk_residual[nwctrks]/F");
    fTree->Branch("wctrk_wcmissed", wctrk_wcmissed,"wctrk_wcmissed[nwctrks]/I");
    fTree->Branch("wctrk_picky", wctrk_picky,"wctrk_picky[nwctrks]/I");
    fTree->Branch("wctrk_WC1XMult", wctrk_WC1XMult,"wctrk_WC1XMult[nwctrks]/I");
    fTree->Branch("wctrk_WC1YMult", wctrk_WC1YMult,"wctrk_WC1YMult[nwctrks]/I");
    fTree->Branch("wctrk_WC2XMult", wctrk_WC2XMult,"wctrk_WC2XMult[nwctrks]/I");
    fTree->Branch("wctrk_WC2YMult", wctrk_WC2YMult,"wctrk_WC2YMult[nwctrks]/I");
    fTree->Branch("wctrk_WC3XMult", wctrk_WC3XMult,"wctrk_WC3XMult[nwctrks]/I");
    fTree->Branch("wctrk_WC3YMult", wctrk_WC3YMult,"wctrk_WC3YMult[nwctrks]/I");
    fTree->Branch("wctrk_WC4XMult", wctrk_WC4XMult,"wctrk_WC4XMult[nwctrks]/I");
    fTree->Branch("wctrk_WC4YMult", wctrk_WC4YMult,"wctrk_WC4YMult[nwctrks]/I");
    fTree->Branch("wctrk_XDist",wctrk_XDist,"wctrk_XDist[nwctrks]/F");
    fTree->Branch("wctrk_YDist",wctrk_YDist,"wctrk_YDist[nwctrks]/F");
    fTree->Branch("wctrk_ZDist",wctrk_ZDist,"wctrk_ZDist[nwctrks]/F");
    fTree->Branch("wctrk_YKink",wctrk_YKink,"wctrk_YKink[nwctrks]/F");
    fTree->Branch("ntof", &ntof, "ntof/I");
    fTree->Branch("tof", tof, "tof[ntof]/F");
    fTree->Branch("tof_timestamp", tof_timestamp, "tof_timestamp[ntof]/F"); 
    if( fSaveWireChamberHits ) {
      fTree->Branch("XWireHist",XWireHist,"XWireHist[nwctrks][1000]/F");
      fTree->Branch("YWireHist",YWireHist,"YWireHist[nwctrks][1000]/F");
      fTree->Branch("XAxisHist",XAxisHist,"XAxisHist[nwctrks][1000]/F");
      fTree->Branch("YAxisHist",YAxisHist,"YAxisHist[nwctrks][1000]/F");
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
    }
  }

  if( fSaveAerogelInfo ) {
    fTree->Branch("nAG", &nAG, "nAG/I");
    fTree->Branch("HitTimeStamp1p10_1", HitTimeStamp1p10_1, "HitTimeStamp1p10_1[nAG]/F");
    fTree->Branch("HitTimeStamp1p10_2", HitTimeStamp1p10_2, "HitTimeStamp1p10_2[nAG]/F");
    fTree->Branch("HitTimeStamp1p06_1", HitTimeStamp1p06_1, "HitTimeStamp1p06_1[nAG]/F");
    fTree->Branch("HitTimeStamp1p06_2", HitTimeStamp1p06_2, "HitTimeStamp1p06_2[nAG]/F");
    fTree->Branch("HitPulseArea1p10_1", HitPulseArea1p10_1, "HitPulseArea1p10_1[nAG]/F");
    fTree->Branch("HitPulseArea1p10_2", HitPulseArea1p10_2, "HitPulseArea1p10_2[nAG]/F");
    fTree->Branch("HitPulseArea1p06_1", HitPulseArea1p06_1, "HitPulseArea1p06_1[nAG]/F");
    fTree->Branch("HitPulseArea1p06_2", HitPulseArea1p06_2, "HitPulseArea1p06_2[nAG]/F");
    fTree->Branch("HitExist1p10_1", HitExist1p10_1, "HitExist1p10_1[nAG]/O");
    fTree->Branch("HitExist1p10_2", HitExist1p10_2, "HitExist1p10_2[nAG]/O");
    fTree->Branch("HitExist1p06_1", HitExist1p06_1, "HitExist1p06_1[nAG]/O");
    fTree->Branch("HitExist1p06_2", HitExist1p06_2, "HitExist1p06_2[nAG]/O");
  }
    
  if( fSaveSimChannelInfo ) {
    fTree->Branch("maxTrackIDE", &maxTrackIDE, "maxTrackIDE/I");
    fTree->Branch("IDEEnergy", IDEEnergy, "IDEEnergy[maxTrackIDE]/F");
    fTree->Branch("IDEPos", IDEPos, "IDEPos[maxTrackIDE][3]/F");
  }
  
  if( fSaveGeantInfo ) {
    fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
    fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
    fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
    fTree->Branch("Mass",Mass,"Mass[geant_list_size]/F");
    fTree->Branch("Eng",Eng,"Eng[geant_list_size]/F");
    fTree->Branch("Px",Px,"Px[geant_list_size]/F");
    fTree->Branch("Py",Py,"Py[geant_list_size]/F");
    fTree->Branch("Pz",Pz,"Pz[geant_list_size]/F");
    fTree->Branch("EndEng",EndEng,"EndEng[geant_list_size]/F");
    fTree->Branch("EndPx",EndPx,"EndPx[geant_list_size]/F");
    fTree->Branch("EndPy",EndPy,"EndPy[geant_list_size]/F");
    fTree->Branch("EndPz",EndPz,"EndPz[geant_list_size]/F");
    fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
    fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
    fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
    fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
    fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
    fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
    fTree->Branch("StartT",StartT,"StartT[geant_list_size]/F");
    fTree->Branch("EndT",EndT,"EndT[geant_list_size]/F");
    fTree->Branch("PathLenInTpcAV",PathLenInTpcAV,"PathLenInTpcAV[geant_list_size]/F");
    fTree->Branch("StartInTpcAV",StartInTpcAV,"StartInTpcAV[geant_list_size]/O");
    fTree->Branch("EndInTpcAV",EndInTpcAV,"EndInTpcAV[geant_list_size]/O");
    fTree->Branch("Process", Process, "Process[geant_list_size]/I");
    fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
    fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
    fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
    fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");
    fTree->Branch("G4Process",&G4Process);//,"G4Process[geant_list_size]");
    fTree->Branch("G4FinalProcess",&G4FinalProcess);//,"G4FinalProcess[geant_list_size]");  
    if( fSaveGeantTrajectories ) { 
      fTree->Branch("NTrTrajPts",NTrTrajPts,"NTrTrajPts[no_primaries]/I");
      fTree->Branch("MidPosX",MidPosX,"MidPosX[no_primaries][5000]/F");
      fTree->Branch("MidPosY",MidPosY,"MidPosY[no_primaries][5000]/F");
      fTree->Branch("MidPosZ",MidPosZ,"MidPosZ[no_primaries][5000]/F");
      fTree->Branch("MidPx",MidPx,"MidPx[no_primaries][5000]/F");
      fTree->Branch("MidPy",MidPy,"MidPy[no_primaries][5000]/F");
      fTree->Branch("MidPz",MidPz,"MidPz[no_primaries][5000]/F");
      fTree->Branch("InteractionPoint"         ,&InteractionPoint         );
      fTree->Branch("InteractionPointType"     ,&InteractionPointType     );
    }
  }

  if( fSaveMCShowerInfo ) {
    fTree->Branch("no_mcshowers", &no_mcshowers, "no_mcshowers/I");
    fTree->Branch("mcshwr_origin", mcshwr_origin, "mcshwr_origin[no_mcshowers]/F");
    fTree->Branch("mcshwr_pdg", mcshwr_pdg, "mcshwr_pdg[no_mcshowers]/F");
    fTree->Branch("mcshwr_TrackId", mcshwr_TrackId, "mcshwr_TrackId[no_mcshowers]/I");
    fTree->Branch("mcshwr_startX", mcshwr_startX, "mcshwr_startX[no_mcshowers]/F");
    fTree->Branch("mcshwr_startY", mcshwr_startY, "mcshwr_startY[no_mcshowers]/F");
    fTree->Branch("mcshwr_startZ", mcshwr_startZ, "mcshwr_startZ[no_mcshowers]/F");
    fTree->Branch("mcshwr_endX", mcshwr_endX, "mcshwr_endX[no_mcshowers]/F");
    fTree->Branch("mcshwr_endY", mcshwr_endY, "mcshwr_endY[no_mcshowers]/F");
    fTree->Branch("mcshwr_endZ", mcshwr_endZ, "mcshwr_endZ[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngX", mcshwr_CombEngX, "mcshwr_CombEngX[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngY", mcshwr_CombEngY, "mcshwr_CombEngY[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngZ", mcshwr_CombEngZ, "mcshwr_CombEngZ[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngPx", mcshwr_CombEngPx, "mcshwr_CombEngPx[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngPy", mcshwr_CombEngPy, "mcshwr_CombEngPy[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngPz", mcshwr_CombEngPz, "mcshwr_CombEngPz[no_mcshowers]/F");
    fTree->Branch("mcshwr_CombEngE", mcshwr_CombEngE, "mcshwr_CombEngE[no_mcshowers]/F");
    fTree->Branch("mcshwr_dEdx", mcshwr_dEdx, "mcshwr_dEdx[no_mcshowers]/F");
    fTree->Branch("mcshwr_StartDirX", mcshwr_StartDirX, "mcshwr_StartDirX[no_mcshowers]/F");
    fTree->Branch("mcshwr_StartDirY", mcshwr_StartDirY, "mcshwr_StartDirY[no_mcshowers]/F");
    fTree->Branch("mcshwr_StartDirZ", mcshwr_StartDirZ, "mcshwr_StartDirZ[no_mcshowers]/F");
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
    fTree->Branch("LengthShw",LengthShw,"LengthShw[nshowers]/F");
    fTree->Branch("CosStartShw",CosStartShw,"CosStartShw[3][1000]/F");
    // fTree->Branch("CosStartSigmaShw",CosStartSigmaShw,"CosStartSigmaShw[3][nshowers]/F");
    fTree->Branch("CosStartXYZShw",CosStartXYZShw,"CosStartXYZShw[3][1000]/F");
    //fTree->Branch("CosStartXYZSigmaShw",CosStartXYZSigmaShw,"CosStartXYZSigmaShw[3][nshowers]/F");
    fTree->Branch("TotalEShw",TotalEShw,"TotalEShw[2][1000]/F");
    //fTree->Branch("TotalESigmaShw",TotalESigmaShw,"TotalESigmaShw[2][nshowers]/F");
    fTree->Branch("dEdxPerPlaneShw",dEdxPerPlaneShw,"dEdxPerPlaneShw[2][1000]/F");
    //fTree->Branch("dEdxSigmaPerPlaneShw",dEdxSigmaPerPlaneShw,"dEdxSigmaPerPlaneShw[2][nshowers]/F");
    fTree->Branch("TotalMIPEShw",TotalMIPEShw,"TotalMIPEShw[2][1000]/F");
    //fTree->Branch("TotalMIPESigmaShw",TotalMIPESigmaShw,"TotalMIPESigmaShw[2][nshowers]/F");
  }
      
}


bool lariat::AnaTreeT1034::IsPointInTpcAV(TVector3& v){
  // Get active volume boundary.
  art::ServiceHandle<geo::Geometry> geom;
  float   xmin = 0. +1e-8;
  float   xmax = 2.0 * geom->DetHalfWidth() - 1e-8;
  float   ymin = -geom->DetHalfHeight() + 1e-8;
  float   ymax = geom->DetHalfHeight()  - 1e-8;
  float   zmin = 0. +1e-8;
  float   zmax = geom->DetLength() - 1e-8;
  if(   (v.X()>=xmin && v.X()<=xmax)
    &&  (v.Y()>=ymin && v.Y()<=ymax)
    &&  (v.Z()>=zmin && v.Z()<=zmax) ) 
    return true;
  else 
    return false;
}

bool lariat::AnaTreeT1034::IsPointInTpcAV(const simb::MCParticle* part, int i){
  TVector3 p(part->Vx(i),part->Vy(i), part->Vz(i));
  return IsPointInTpcAV(p);
}


float lariat::AnaTreeT1034::TrajLengthInTpcAV(const simb::MCParticle* part){
  float L = 0.; 
  int nPts = part->NumberTrajectoryPoints();
  for(int i = 1; i < nPts; ++i) {
    // if both this and previous point are in AV, add to total traj length
    TVector3 pA(part->Vx(i-1),part->Vy(i-1),part->Vz(i-1));
    TVector3 pB(part->Vx(i),part->Vy(i),part->Vz(i));
    if( IsPointInTpcAV(pA) && IsPointInTpcAV(pB) ) {
      L += (pB-pA).Mag();
    }
  }
  return L;
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
  lifetime = -999;
  //trig_timestamp = -999;
  efield = -99999;
  t0 = -99999;
  nclus = 0;
  for (int i = 0; i < kMaxCluster; ++i){
    clustertwire[i] = -99999;
    clusterttick[i] = -99999;
    cluendwire[i] = -99999;
    cluendtick[i] = -99999;
    cluplane[i] = -99999;
  }
  ntracks_reco = 0;
  for (int i = 0; i < kMaxTrack; ++i){
    trkvtxx[i] = -99999;
    trkvtxy[i] = -99999;
    trkvtxz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkWCtoTPCMatch[i] = -99999;
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    trklength[i] = -99999;
    //trkmomrange[i] = -99999;
    //trkmommschi2[i] = -99999;
    //trkmommsllhd[i] = -99999;
    ntrkpts[i] = 0;
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
      ntrkcalopts[i][j] = 0; 
      trkke[i][j] = -99999;
      trkpida[i][j] = -99999;
      for (int k = 0; k<kMaxTrackHits; ++k){
        trkdedx[i][j][k] = -99999;
        trkdqdx[i][j][k] = -99999;
        trkrr[i][j][k] = -99999;
        trkpitch[i][j][k] = -99999;
        trkxyz[i][j][k][0] = -99999;
        trkxyz[i][j][k][1] = -99999;
        trkxyz[i][j][k][2] = -99999;
      }
    }
    trkg4id[i] = -9999;
  }
  primarytrkkey = -9999;
  nhits = 0;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -9;
    hit_wire[i] = -9;
    hit_channel[i] = -9;
    hit_peakT[i] = -99999;
    hit_driftT[i] = -99999;
    hit_charge[i] = -99999;
    hit_electrons[i] = -99999;
    hit_ph[i]     = -9999;
    hit_trkid[i]  = -9;
    hit_clusterid[i] = -999;
    hit_tstart[i] = -999;
    hit_tend[i] = -999;
    //hit_pk[i] = -99999;
    //hit_t[i] = -99999;
    //hit_ch[i] = -99999;
    //hit_fwhh[i] = -99999;
    hit_rms[i] = -999;
    hit_dQds[i] = -999;
    hit_dEds[i] = -999;
    hit_ds[i] = -999;
    hit_resrange[i] = -999;
    hit_x[i] = -999;
    hit_y[i] = -999;
    hit_z[i] = -999;
    hit_g4id[i] = -999;
    hit_g4frac[i] = -999;
    hit_g4nelec[i] = -999;
    hit_g4energy[i] = -999;
  }

  nwctrks = 0;

  for (int i = 0; i < kMaxWCTracks; i++)
    {
      wctrk_XFaceCoor[i] = -99999;	//<---The projected X position of the wctrack at the front face of the TPC
      wctrk_YFaceCoor[i] = -99999;	//<---The projected Y position of the wctrack at the front face of the TPC
      wctrk_momentum[i] = -99999;		//<---Reconstructed moomentum
      wctrk_Px[i] = -99999;
      wctrk_Py[i] = -99999;
      wctrk_Pz[i] = -99999;
      wctrk_theta[i] = -99999;		//<---angle of track w.r.t. z axis
      wctrk_phi[i] = -99999;		//<---angle of track w.r.t. x axis
      wctrk_XDist[i] = -99999;          //<---X distance between upstream and downstream tracks
      wctrk_YDist[i] = -99999;          //<---Y distance between upstream and downstream tracks
      wctrk_ZDist[i] = -99999;
      wctrk_YKink[i] = -99999;

      for(int j = 0; j < 1000; j++)
        {
          XWireHist[i][j] = -99999;		//<---Coord in terms of wire number
          YWireHist[i][j] = -99999;		//<---Coord in terms of wire number
          XAxisHist[i][j] = -99999;		//<---coord in terms of units.
          YAxisHist[i][j] = -99999;		//<---coord in terms of units.
        }//<---End j loop
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

  ntof = 0;
  beamline_mass = -99999;
  for (int i = 0; i < kMaxTOF; i++)
    {
      tof[i] = -99999;
      tof_timestamp[i] = -99999;
    }//<---End i loop

  nAG = 0;
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

  maxTrackIDE = 0;

  for(size_t i = 0; i < kMaxIDE; ++i) {
    IDEEnergy[i] = -999;

    IDEPos[i][0] = -999.9;
    IDEPos[i][1] = -999.9;
    IDEPos[i][2] = -999.9;

  } // End of maxTrackID loop

  no_primaries = 0;
  geant_list_size=0;
  for (int i = 0; i<kMaxPrimaries; ++i){
    pdg[i] = -99999;
    Mass[i]=-99999;
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
    StartT[i] = -99999;
    EndT[i]   = -99999;
    PathLenInTpcAV[i] = -9999;
    StartInTpcAV[i] = false;
    EndInTpcAV[i] = false;
    
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

  nshowers = 0;
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
          //     CosStartXYZSigmaShw[j][i] = -99999;
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
