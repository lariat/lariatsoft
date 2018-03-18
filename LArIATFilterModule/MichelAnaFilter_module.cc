////////////////////////////////////////////////////////////////////////
/// Class:       MichelAnaFilter
/// Module Type: filter
/// File:        MichelAnaFilter_module.cc
///
/// Input data products required:
///  - raw::OpDetPulses
///  - recob::Hits
///  - recob::Tracks
///  - anab::Calorimetry (not essential)
/// 
/// Output: 
///  - histogram file + TTree 
///
/// Primary analysis module for the Michel electron studies using data
/// from the light-based trigger on stopping cosmic muons.  It combines 
/// and replaces both MichelWfmReco_module and MichelMCAna_module.
///
/// This module takes in reconstructed events containing PMT data, hits, tracks,
/// and calorimetry objects.  It does PMT waveform reconstruction, looks 
/// for Michel-like events, and saves all relevant information to a TTree
/// to be analyzed by a standalone ROOT macro.  Michel electron charge-based
/// clustering and reconstruction is also performed, following the technique
/// from MicroBooNE.
///
/// Optical waveform processing algorithms from OpHitBuilderAlg are used
/// to find and integrate hits in the PMTs.
///
/// TODO:
///  [ ] Move Michel clustering routines to separate dedicated alg class.
///
/// Generated at Wed Jun 14 13:35:26 2017 by William Foreman using artmod
/// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

// LArIATSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"
#include "lardataobj/RawData/TriggerData.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom2.h"
#include "TVector2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TVirtualFFT.h"

// C++ includes
#include <iostream>
#include <memory>
#include <random>
#include <algorithm>

// Global constants and limits on tracks/hits
float RAD_TO_DEG = 180. / 3.1415926535;
const size_t kMaxOpHits = 20;
const size_t kMaxTracks = 50;

// michelcluster: pecial data structure to handle clustering/shower information
struct michelcluster {
  int               plane;
  int               seedHit;
  std::vector<int>  cluster;
  std::vector<int>  cluster_el;
  std::vector<int>  cluster_mu;
  std::vector<int>  shower;
  std::vector<bool> isInTrk;
  float             minCovAtBnd;
  int               muEndHit;
  float             elShowerFrac;
  TVector3          muEnd2D;
  TVector3          muEndDir2D;
  float             decayAngle2D;
  float             showerAngleMean;
  float             showerAngleRMS;
  int               extraHits;
  std::vector<int>   prof_hitKey;
  std::vector<float> prof_X;
  std::vector<float> prof_W;
  std::vector<float> prof_s;
  std::vector<float> prof_dQds;
  std::vector<float> prof_dQds_t;
  std::vector<float> prof_lin;
  int   bnd_i;
  float bnd_s;
  int maxQ_i;
  float maxQ_s;
  float totalCharge;
  float elShowerCharge;
  float elShowerEnergy;
  float elTrackCharge;
  float elTrackEnergy;
  float fracMuHitsLinear;
  float muAveLinearity;
  int   nPtsMuFit;
  float aveDriftTime;
  float aveX;
  
  michelcluster() {
    plane = -9;
    seedHit = -9;
    muEndHit = -9;
    elShowerFrac = -9.;
    decayAngle2D = -9;
    showerAngleMean = -9.;
    showerAngleRMS = -9.;
    muEnd2D.SetXYZ(-99., -99., -99.);
    muEndDir2D.SetXYZ(-99., -99., -99.);
    minCovAtBnd = 9.;
    bnd_i = -9;
    maxQ_i = -9;
    bnd_s = -999.;
    maxQ_s = -999.;
    extraHits = 0;
    fracMuHitsLinear = -9.;
    nPtsMuFit = -9;
    muAveLinearity = -9.;
    totalCharge = -999.;
    elShowerCharge = -999.;
    elShowerEnergy = -999.;
    elTrackCharge  = -999.;
    elTrackEnergy = -999.;
    aveDriftTime  = -9.;
    aveX          = -9.;
  }

};

//########################################################################################
class MichelAnaFilter;


//########################################################################################
class MichelAnaFilter : public art::EDFilter {

public:
  
  explicit MichelAnaFilter(fhicl::ParameterSet const & p);
  MichelAnaFilter(MichelAnaFilter const &) = delete;
  MichelAnaFilter(MichelAnaFilter &&) = delete;
  MichelAnaFilter & operator = (MichelAnaFilter const &) = delete;
  MichelAnaFilter & operator = (MichelAnaFilter &&) = delete;

  // Art functions
  bool  filter(art::Event & e) override;
  void  beginJob() override;
  void  endJob() override;
  void  reconfigure(fhicl::ParameterSet const & p) override;

  // Custom functions   
  void  ResetVariables();
  void  GetDetProperties();
  void  GetTruthInfo( art::Event &);
  void  ParticleTracker(const art::Event&, int, int,TH1D* h[], TH1D* he[]);
  bool  IsParticleDescendedFrom(const art::Handle< std::vector<simb::MCParticle>>&,int,int);
  bool  IsPointInFiducialVolume(TVector3, float, float, float);
  float CalcParticleDist(const simb::MCParticle& );
  void  PerfectWfmReco(raw::OpDetPulse & pulse);
  void  PropagatePhoton(float,float,float,float,float,bool,TH1D*,TH1D*);
  void  ReadPhotonLibrary(TVector3&,int,int,float&,float&,float&,float&);
  float CalcSigma(size_t,float,float);
  void  MuContamCorr(int);
  void  FindBestMatchedHit(int,float,int&,float&);
  void  FindBestMatchedHit(int,float,int&,float&,std::vector<int>&);
  void  Clustering(int,michelcluster&);
  void  Clustering(int,michelcluster&,float,std::vector<int>&);
  void  ProximityCluster(std::vector<float>&,std::vector<float>&,std::vector<float>&,std::vector<int>&,std::vector<int>&,int,int,float);
  void  CalcMichelShowerEnergy(michelcluster&);
  float CalcTruncatedMeanInProfile(std::vector<float>&,size_t,int,float);
  float CalcLocalLinearity(std::vector<float>&,std::vector<float>&,size_t,int);
  float GetVisibility(TVector3&,int);
  float GetVisibility(TVector3&,int,int);
  int   GetGlobalBin(const TH3D*,double,double,double);
  int   GetGlobalBin(const TH3D*,const TLorentzVector&);
  int   GetGlobalBin(const TH3D*,const TVector3&);
  float Integrate(const TH1D*,float,float); 
  float Integrate(const std::vector<float>&,int,int); 
  void  AddToAverageWfm(std::vector<float>&,short,short,TH1D*,int&,int);
  void  AddToAverageWfm(const TH1D*,short,short,short,TH1D*,int &,int);
  void  MakeClusteringGraphs(michelcluster&);
  void  MakeWfmGraphs(art::Event&);

private:
 
  // ................................................................................
  // Service objects 
  geo::GeometryCore     *fGeo;
  OpHitBuilderAlg       fOpHitBuilderAlg; 
  TTree*                fTree;
  
  // ................................................................................
  // Fhicl parameters 
  
  int                 fRandSeed;                // random seed for reproducibility
  
  bool                fFilter_PassAllEvents;    // pass all events through filter
  bool                fFilter_OpticalMode;      // only pass evts w/ optical Michel ID 
  bool                fFilter_Req3DShower;      // require 3D shower reco'd
  std::vector<float>  fFilter_DecayTime;        // acceptance range of decay times
  std::vector<int>    fFilter_NumTrackStopping; // acceptance range of num. stp trks
  std::vector<float>  fFilter_ElShowerEnergy;   // acceptance range of 2D shwr energy
  
  std::string         fHitsModule;          // producer of input recob::Hits
  std::string         fHitsInstance;        // instance name of input recob::Hits
  std::string         fTrackModule;         // producer of input recob::Tracks
  std::string         fTrackCalModule;      // producer of input anab::Calorimetry
  std::string         fSimProducerModule;   // producer of input sim::MCParticles
  std::string	      fLibraryFile;         // name of visibility library file
 
  bool                fUseCrossingMuons;    // reconstruct single crossing mu evts 
  bool                fLookAtTracks;        // attempt to look at recob::Track info
  bool                fReq1StpTrk;          // require 1 stp trk before doing optical reco
  float               fFiducialMarginX;     // fiducial margin X [cm]
  float               fFiducialMarginY;     // fidicual margin Y [cm]
  float               fFiducialMarginZ;     // fiducial margin Z [cm]
  float               fMuTrackLengthMin;    // min mu candidate track length [cm]
  
  std::vector<size_t> fSelectChannels;      // PMT channels to look at (0=HMM, 1=ETL)
  short               fTruncateWfm;         // truncate PMT wfms to this length
  short               fBaselineWindowLength;// pedestal region of wfm for baseline/RMS calculation
  std::vector<float>  fMaxWfmRMS;           // max wfm RMS (ADC)
  int                 fWfmSmoothingRange;   // local window for waveform smoothing (+/- val)
  std::string         fCorrectOvershootMode;// overshoot correciton mode ("int","amp") 
  bool                fMskBaselineSubtr;    // do masked baseline subtraction on wfms
  std::vector<float>  fGradHitThresh;       // gradient threshold for hit-finding (PMT-specific)
  std::vector<float>  fMinOpHitWidth;       // min ophit FWHM
  bool                fDeleteNarrowHits;    // delete/ignore hits that fail ophit width cut
  bool                fRequireHitMatching;  // req. hits match between PMTs
  std::vector<float>  fOpHitMatchThresh;    // op hits below this amp. are matched (PMT-specific)
  float               fMaxDiff_dT;          // max T diff btwn. matched ophits or decay times
  float               fGateDelay;           // min decay time allowed
  float               fMaxDecayTime;        // max decay time allowed
  short               fPromptWindow;        // prompt light int. window (default 100ns)
  short               fFullWindow;          // total light int. window (default 7us)
  
  std::vector<float>  fSmearSigT;           // integration time-dependent smearing on MC
  std::vector<float>  fSmearSigPE;          // photoelectron resolution smearing on MC
  std::vector<float>  fSmearFactor;         // final smear factor on MC total light
  std::vector<float>  fSmearFactorPrompt;   // final smear factor on MC prompt light (NOT USED)
  std::vector<float>  fSPE_ScaleFactor;     // scale PMT single photoelectron response
  std::vector<float>  fCollectionEff;       // PMT coll. efficiency for use in MC
  std::vector<float>  fQE_ScaleFactor;      // global QE scale factor (from optimization)
  std::vector<float>  fQE_ScaleFactor_Vis;  // QE scale factor for vis. light (NOT USED)
  std::vector<float>  fQE_ScaleFactor_VUV;  // QE scale factor for VUV light (NOT USED)
  float               fFastLightRatio;      // singlet-to-triplet ratio used in MC
  float               fLateLightTau;        // eff. late-light tau in MC (determines quenching)
  float               fTauT;                // "true" triplet dimer tau (influences quenching)
  float               fTauS;                // "true" singlet dimer tau
  bool                fApplyTrigEffCut;     // cut out MC events based on trigger eff. funct.
  std::vector<float>  fTrigEff_P;           // trigger eff. funct. param "P" (from optimization)
  std::vector<float>  fTrigEff_K;           // trigger eff. funct. param "K" (fixed to 8)
  
  bool                fMuContamCorr;        // correct for muon late-light contamination
  float               fMuContamCorr_dTmin;  // min dT required for mu contam. correction
  std::vector<float>  fMuContamCorr_EffTau; // eff. late-light tau for mu contam corr. (from DB!) 

  float               fMaxHitSeparation;    // max 2D dist (in WX space) btwn adjacent hits
  int                 fTruncMeanWindow;     // size local nbd. for trunc. mean in profile (+/-N)
  float               fTruncMeanP;          // frac of upper/lower tail to truncate
  int                 fMaxMatchDistance;    // dist (in hits) to find max Q/ds after trunc mean peak
  int                 fLocalLinearityWindow;// local nbd. in profile to calc local linearity (+/-N)
  float               fLinThresh;           // max local linearity allowed at mu-el boundary
  float               fLinTol;              // adjusts thresh based on RMS of mu track linearity
  float               fLinThreshMuFit;      // min local linearity for "good" mu track hits
  float               fShowerAcceptanceAngle;// opening angle of 2D electron shower cone
  
  std::vector<float>  fCalAreaConstants;    // scale hit integrals to charge (ADC -> e-)
  float               fWph;                 // energy to excitate/ionize Ar atom (19.5 eV)
  float               fWion;                // energy to produce e-ion pair (23.6 eV)
  float               fExcRatio;            // excitation ratio ( N_ex / N_i ) = 0.21
  float               fMvPerADC;            // scale ADC to mV
  float               fRecombFactor;        // recomb. survival frac for electron trk/shwr hits
  float               fRecombFactorTrack;   // recomb. survival frac for track-like el. hits
  float               fRecombFactorShower;  // recomb. survival frac for shwr-like el. hits

  float               fdTcut;               // min decay time for energy histograms 
  std::vector<float>  fPromptPECut;         // min prompt PE for dT histogram 
  int                 fMinClusterSize;      // min allowable cluster size for histos
  float               fMinElShowerFrac;     // min electron shower completeness frac for histos
  float               fMinFracMuHitsLinear; // min frac of good linearity mu trk hits for histos
  int                 fMinMuClusterSize;    // min mu cluster size for histos
  int                 fMinMuClusterHitsEndFit;  // min mu hits used in direction fit for histos
  int                 fMinElClusterSize;    // min el cluster size for histos
  int                 fMaxElClusterSize;    // max el cluster size for histos
  float               fMinMuLinearity;      // min average linearity of mu hits for histos
  int                 fMaxExtraHits;        // max num. of unclustered hits near mu end for histos
  float               fMuEndRadius;         // radius (WX space) around mu endpt to find extra hits
  float               fMaxDecayAngle2D;     // max 2D decay angle for histos
  float               fMinDecayAngle2D;     // min 2D decay angle for histos
  float               fMaxHitTimeDiff;      // max allowable time diff to match hits btwn planes
  int                 fMinNumPts3D;         // min num 3D pts for histos
  float               fMinFracHits3D;       // min frac of hits made into 3D points for histos

  int                 fBndOffset;           // force offset on cluster bnd. by +/- N hits (for 
                                            // estimatation of systematics)
  
  int                 fMaxSavedHitGraphs;   // save this many clustering displays
  int                 fMaxSavedWfmGraphs;   // save this many PMT waveform displays
  
  // ................................................................................
  // Event-level data that could be saved into a tree

  // Event identifiers
  bool                fIsMC;
  int                 fRunNumber;
  int                 fSubRunNumber;
  int                 fEventNumber;
  int                 fEventTime;
  float               fElectronLifetime;
  bool                fMichelOpticalID;
  
  // Track information
  int                 fNumTrack;
  int                 fNumTrackStopping;
  int                 fNumTrackCrossing;
  int                 fNumTrackContained;
  std::vector<int>    fMuTrackHits;
  int                 fMuTrackIndex;
  int                 fMuTrackID;
  float               fMuTrackLength;
  float               fMuTrackEnergy;
  float               fMuTrackZenithAngle;
  float               fMuTrackVertex_X;
  float               fMuTrackVertex_Y;
  float               fMuTrackVertex_Z;
  float               fMuTrackEnd_X;
  float               fMuTrackEnd_Y;
  float               fMuTrackEnd_Z;
  std::vector<double> fvMuTrackdEdx;
  std::vector<double> fvMuTrackResidualRange;
  bool                fMuTrackIsBraggPeaked;
  bool                fMuTrackIsCaloOrdered;
  TVector3            fMuTrackVertex;
  TVector3            fMuTrackEnd;
  int                 fCrossingMuTrackIndex;
  float               fCrossingMuVertex_X;
  float               fCrossingMuVertex_Y;
  float               fCrossingMuVertex_Z;
  float               fCrossingMuEnd_X;
  float               fCrossingMuEnd_Y;
  float               fCrossingMuEnd_Z;
  float               fCrossingMuLength;
  float               fCrossingMuPhotons;
  float               fCrossingMuPhotonsPrompt;
  float               fCrossingMuCharge;

  // Optical information
  float               fWfmRMS[2];
  float               fSPE[2];
  float               fSPE_err[2];
  int                 fNumOpHits0[2];
  int                 fNumOpHits[2];
  float               fdT[2];
  float               fAmplitude[2];
  float               fMuWidth[2];
  float               fWidth[2];
  
  // Hit times, trigger condition, saturation flag
  std::vector<short>  fvHitTimes0[2];
  std::vector<short>  fvHitTimes[2];
  std::vector<bool>   fvIsHitAtTrigger[2];
  std::vector<bool>   fvIsHitSaturated[2]; 

  // Raw hit integrals
  std::vector<float>  fvHitADC_100ns[2];
  std::vector<float>  fvHitADC_200ns[2];
  std::vector<float>  fvHitADC_300ns[2];
  std::vector<float>  fvHitADC_400ns[2];
  std::vector<float>  fvHitADC_500ns[2];
  std::vector<float>  fvHitADC_600ns[2];
  std::vector<float>  fvHitADC_700ns[2];
  std::vector<float>  fvHitADC_900ns[2];
  std::vector<float>  fvHitADC_1200ns[2];
  std::vector<float>  fvHitADC_1500ns[2];
  std::vector<float>  fvHitADC_1800ns[2];
  std::vector<float>  fvHitADC_7000ns[2];

  // Prepulse fit hit integrals, hit amplitude, width
  std::vector<float>  fvHitAmplitudes[2];
  std::vector<float>  fvHitWidth[2];
  std::vector<float>  fvHitADCpf_100ns[2];
  std::vector<float>  fvHitADCpf_7000ns[2];
  std::vector<short>  fvPrepulseX1[2];
  std::vector<float>  fvPrepulseZeroPoint[2];
  std::vector<float>  fvPrepulseBaseline[2];
  std::vector<float>  fvPrepulseRMS[2];
  std::vector<float>  fvPrepulseSlowNorm[2];
  std::vector<float>  fvPrepulseSlowTau[2];

  // Derived quantities when dT is defined (2 hits)
  float               fPE_promptRaw[2];
  float               fPE_totalRaw[2];
  float               fPE_prompt[2];
  float               fPE_total[2];
  float               fPromptFracRaw[2];
  float               fPromptFrac[2];
  float               fMuAmplitude[2];
  float               fMuPE_prompt[2];
  bool                fMuPulseSaturated[2];
  bool                fElPulseSaturated[2];
  float               fPE;
  float               fPEPrompt;
  float               fDecayTime;
    
  // 2D clustering data
  int                 fNumPlaneHits[2];
  float               fAveDriftTime;
  float               fAveX;
  float               fElOffsetT;
  int                 fElStartWire;
  int                 fMuEndWire;
  // ...collection plane
  float               fMinCovAtBnd;
  int                 fClusterSize;
  int                 fMuClusterSize;
  int                 fElClusterSize;
  int                 fExtraHits;
  float               fMuAveLinearity;
  float               fFracMuHitsLinear;
  int                 fMuClusterHitsEndFit;
  float               fDecayAngle2D;
  float               fElShowerAngleMean;
  float               fElShowerAngleRMS;
  int                 fElShowerSize;
  float               fElShowerFrac; 
  float               fTotalCharge;
  float               fElShowerCharge;
  float               fElTrackCharge;
  float               fElShowerEnergy;
  float               fElTrackEnergy;                 
  // ...induction plane
  float               fMinCovAtBnd_Pl0;
  int                 fClusterSize_Pl0;
  int                 fMuClusterSize_Pl0;
  int                 fElClusterSize_Pl0;
  int                 fExtraHits_Pl0;
  float               fMuAveLinearity_Pl0;
  float               fFracMuHitsLinear_Pl0;
  int                 fMuClusterHitsEndFit_Pl0;
  float               fDecayAngle2D_Pl0;
  float               fElShowerAngleMean_Pl0;
  float               fElShowerAngleRMS_Pl0;
  int                 fElShowerSize_Pl0;
  float               fElShowerFrac_Pl0; 
  float               fTotalCharge_Pl0;
  float               fElShowerCharge_Pl0;
  float               fElTrackCharge_Pl0;
  float               fElShowerEnergy_Pl0;
  float               fElTrackEnergy_Pl0;                 
  
  // 3D shower information
  float               fFracHits3D;
  int                 fNumPts3D;
  int                 fNumPts3DTrk;
  float               fElShowerVis;
  float               fElShowerVisCh[2];
//  TVector3            fElShowerCentroid;
  float               fElShowerCentroid_X;
  float               fElShowerCentroid_Y;
  float               fElShowerCentroid_Z;
  float               fElShowerPhotons;
  float               fElShowerEnergyQL;
//  TVector3            fMuEnd3D;
  float               fMuEnd3D_X;
  float               fMuEnd3D_Y;
  float               fMuEnd3D_Z;
  float               fMuEndHitTimeDiff;
  float               fElTrackdEdx;

  // Truth information
  int                 fTrue_MuCharge;
  float               fTrue_MuStopTime;
  TVector3            fTrue_MuTrackVertex;
  TVector3            fTrue_MuTrackEnd;
  float               fTrue_MuTrackVertex_X;
  float               fTrue_MuTrackVertex_Y;
  float               fTrue_MuTrackVertex_Z;
  float               fTrue_MuTrackEnd_X;
  float               fTrue_MuTrackEnd_Y;
  float               fTrue_MuTrackEnd_Z;
  float               fTrue_MuEnergyDep;
  int                 fTrue_NumBremmPhotons;
  float               fTrue_ElEnergy;
  float               fTrue_ElShowerPhotons;
  float               fTrue_ElTrackPhotons;
  float               fTrue_ElShowerPhotonsPrompt;
  float               fTrue_ElShowerPhotonsLate;
  float               fTrue_ElTrackEnergyDep;
  float               fTrue_ElTrackChargeDep;
  float               fTrue_ElShowerEnergyDep;
  float               fTrue_ElShowerChargeDep;
  float               fTrue_ElShowerVis;
  float               fTrue_ElShowerVisCh[2];
  TVector3            fTrue_ElMomentum;
  float               fTrue_ElMomentum_X;
  float               fTrue_ElMomentum_Y;
  float               fTrue_ElMomentum_Z;
  float               fTrue_ElAngle;
  bool                fTrue_IsElContained;
  float               fTrue_TotalChargeDep;
  float               fTrue_TotalEnergyDep;
  float               fTrue_dT;
  float               fTrue_Amplitude[2];
  float               fTrue_MuPE_prompt[2];
  float               fTrue_MuContam_prompt[2];
  float               fTrue_MuContam_total[2];
  float               fTrue_PE[2];
  float               fTrue_PE_prompt[2];
  float               fTrue_PE_total[2];
 
  // ..............................................................................
  // Histograms
  TH1D*               hEventCuts;
  TH1D*               hOpEventCuts;
  TH1D*               hWfmRMS[2];
  TH1D*               hBaselinePE_100ns[2];
  TH1D*               hBaselinePE_500ns[2];
  TH1D*               hBaselinePE_1000ns[2];
  TH1D*               hNumOpHits[2];
  TH2D*               hNumOpHitsCompare;
  TH2D*               hPECompare;
  TH1D*               hPE_promptRaw[2];
  TH1D*               hPE_promptRaw_dTcut[2];
  TH1D*               hPE_promptRaw_dTcut_shwr[2];
  TH1D*               hPE_totalRaw[2];
  TH1D*               hPE_totalRaw_dTcut[2];
  TH1D*               hPE_totalRaw_dTcut_shwr[2];
  TH1D*               hPE_prompt[2];
  TH1D*               hPE_prompt_dTcut[2];
  TH1D*               hPE_prompt_dTcut_shwr[2];
  TH1D*               hPE_total[2];
  TH1D*               hPE_total_dTcut[2];
  TH1D*               hPE_total_dTcut_shwr[2];
  TH2D*               hdT_vs_PE_promptRaw[2];
  TH2D*               hdT_vs_PE_prompt[2];
  TH2D*               hdT_vs_PE_totalRaw[2];
  TH2D*               hdT_vs_PE_total[2];
  TH2D*               hdT_vs_Amp[2];
  TH1D*               hPromptFracRaw[2];
  TH1D*               hPromptFrac[2];
  TH1D*               hPromptFracRaw_dTcut[2];
  TH1D*               hPromptFrac_dTcut[2];
  TH1D*               hMuTau[2];
  TH1D*               hAmplitude[2];
  TH1D*               hOpHitAmplitude[2];
  TH1D*               hOpHitAmplitude_unmatched[2];
  TH1D*               hHitPromptPE[2];
  TH1D*               hHitPromptPE_unmatched[2];
  TH1D*               hOpHitTime[2];
  TH1D*               hWidth[2];
  TH1D*               hWidth_AllHits[2];
  TH2D*               hAmplitudeVsWidth[2];
  TH2D*               hPromptPEVsWidth[2];
  TH2D*               hAmplitudeVsWidth_AllHits[2];
  TH2D*               hPromptPEVsWidth_AllHits[2];
  TH2D*               hPromptPEVsAmplitude[2];
  TH2D*               hMuPromptPEVsAmplitude[2];
  TH1D*               hdT_smallWidth[2];
  TH1D*               hdT[2];
  TH1D*               hDiff_dT;
  TH1D*               hDecayTime;
  TH1D*               hPE;
  TH1D*               hMudEdx;
  TH2D*               hMuResRangeVsdEdx;
  TH1D*               hExtraHits;
  TH1D*               hHitRMS[2];
  TH1D*               hHitX[2];
  TH1D*               hHitAmplitude[2];
  TH1D*               hHitCharge;
  TH1D*               hHitSummedADC;
  TH1D*               hHitIntegral;
  TH1D*               hTotalCharge;
  TH1D*               hClusterHitSeparation;
  TH1D*               hElShowerCharge;
  TH1D*               hElShowerPhotons;
  TH1D*               hElShowerEnergy;
  TH1D*               hElShowerEnergyQL;
  TH1D*               hElTrackdEdx;
  TH1D*               hQoverL;
  TH1D*               hQoverLCrsMuTrk;
  TH1D*               hClusterSize;
  TH1D*               hElShowerFrac;
  TH1D*               hDistMaxTdQdsToMaxQ;
  TH1D*               hFracMuHitsLinear;
  TH1D*               hMuAveLinearity;
  TH1D*               hMuClusterSize;
  TH1D*               hElClusterSize;
  TH1D*               hElShowerSize;
  TH1D*               hMuClusterHitsEndFit; 
  TH1D*               hDecayAngle2D;
  TH1D*               hRecomb;
  TH1D*               hRecombCrsMuTrk;
  TH1D*               hElShowerDepAngle2D;
  TH1D*               hElShowerAngleMean;
  TH1D*               hElShowerAngleRMS;
  TH1D*               hNumPts3D;
  TH1D*               hFracHits3D;
  TH1D*               hElOffsetT;
  TH2D*               hElClusterSizeCompare;
  TH2D*               hShowerSizeCompare;
  TH1D*               hHitTimeDiff;
  TH1D*               hMuEndHitTimeDiff;
  TH2D*               hMuEnd3D_ZX;
  TH2D*               hMuEnd3D_ZY;
  TH1D*               hCrsMuQPerCm;
  TH1D*               hCrsMuLPerCm;
  TH1D*               hElectronLifetime;
  TH1D*               hEffTau[2];
  TH1D*               hElectronLifetimeReco;
  TH1D*               hdQdxVsT;
  TH1D*               hdQdxVsT_N;
  TH1D*               hTrue_NumElectrons;
  TH1D*               hTrue_NumDRays;
  TH1D*               hTrue_DRayEnergy;
  TH1D*               hTrue_ElEnergy;
  TH1D*               hTrue_ElEnergyFree;
  TH1D*               hTrue_ElEnergyCapture;
  TH1D*               hTrue_ElTrackEnergyDep;
  TH1D*               hTrue_ElShowerEnergyDep; 
  TH1D*               hTrue_ElShowerChargeDep;
  TH1D*               hTrue_ElTrackChargeDep;
  TH1D*               hTrue_ElShowerPhotons;
  TH1D*               hTrue_PE[2];
  TH1D*               hTrue_PE_preTrig[2];
  TH1D*               hTrue_PE_prompt[2];
  TH1D*               hTrue_PE_total[2];
  TH1D*               hTrue_PE_total_dTcut_shwr[2];
  TH2D*               hTrue_PE_total_TrueVsReco[2];
  TH2D*               hTrue_MuLateLightContamination;
  TH2D*               hTrue_MuTrackEnd_ZX;
  TH2D*               hTrue_MuTrackEnd_ZY;
  TH1D*               hTrue_MuTrackVertexRes_X;
  TH1D*               hTrue_MuTrackVertexRes_Y;
  TH1D*               hTrue_MuTrackVertexRes_Z;
  TH1D*               hTrue_MuTrackEndRes_X;
  TH1D*               hTrue_MuTrackEndRes_Y;
  TH1D*               hTrue_MuTrackEndRes_Z;
  TH2D*               hTrue_EnergyDepVsRecoEnergy;
  TH1D*               hTrue_EnergyRes;
  TH1D*               hTrue_EnergyRes_Shwr3D;
  TH1D*               hTrue_EnergyResTrk;
  TH1D*               hTrue_EnergyResQL;
  TH2D*               hTrue_EnergyVsEnergyRes;
  TH1D*               hTrue_PERes[2];
  TH1D*               hTrue_VisRes;
  TH1D*               hTrue_PhotonRes;
  TH1D*               hTrue_TrajLength;
  TH1D*               hTrue_dEdx;
  TH1D*               hTrue_dEdx_ElTrk;
  TH1D*               hTrue_dEdx_ElShw;
  TH1D*               hTrue_dQdx;
  TH1D*               hTrue_dQdx_ElTrk;
  TH1D*               hTrue_dQdx_ElShw;
  TH2D*               hTrue_dEdx_vs_dQdx;
  TH2D*               hTrue_dEdx_vs_dQdx_ElTrk;
  TH2D*               hTrue_dEdx_vs_dQdx_ElShw;
  TH1D*               hTrue_RecombFactor;
  TH1D*               hTrue_RecombFactor_ElTrk;
  TH1D*               hTrue_RecombFactor_ElShw;
  TH1D*               hTrue_ElShowerDepVsDistFromMuTrackEnd;
  TH1D*               hTrue_VisTotal[2];
  TH1D*               hTrue_Vis[2][2];
  TH1D*               hTrue_X;
  TH1D*               hTrue_Y;
  TH1D*               hTrue_Z;
  TH1D*               hTrue_LightYield[2];
  TH1D*               hTrue_BremmPhotonLength;
  TH2D*               hTrue_ElEnergyVsNumBremmPhotons;
  TH1D*               hNumTrack;
  TH1D*               hNumTrackStopping;
  TH1D*               hNumTrackContained;
  TH1D*               hNumTrackCrossing;
  TH1D*               hTrackNode_X;
  TH1D*               hTrackNode_Y;
  TH1D*               hTrackNode_Z;
  TH2D*               hTrackNode_ZX;
  TH2D*               hTrackNode_ZY;
  TH2D*               hCrsTrackNode_ZX;
  TH2D*               hCrsTrackNode_ZY;
  TH1D*               hTrackLength;

  // ...........................................................
  // Average waveforms (muon, electron pulse)
  TH1D*               hAveWfm_mu[2];
  TH1D*               hAveWfm_el[2];
  int                 aveWfmCounter_mu[2];
  int                 aveWfmCounter_el[2];
  TH1D*               hAvePhelProfile_mu[2];
  int                 avePhelProfileCounter_mu[2];

  // ...........................................................
  // Constant, buffers, counters, misc... 
  std::vector<raw::OpDetPulse> fPulses;
  std::vector<short>    fIntegrationWindows;
  int                 fCachedRunNumber;
  char	              buffer[200], histName[200], histTitle[200];
  float               fEfield;
  float               fDriftVelocity; 
  float               fSamplingRate; 
  float               fXTicksOffset[2]; 
  int                 fElectronID;
  int                 fMuonID;
  float               fT0[2];
  float               fMuT0[2];
  size_t              numMichelWfms[2];
  size_t              numMuSaturated[2];
  size_t              numElSaturated[2];
  size_t              fNumOpticalID;
  size_t              fNumEvents;
  int                 fNumEvents_StpTrk;
  size_t              fNumEvents_Shwr;
  int                 fNumEvents_Shwr3D;
  int                 fNumEvents_Calibration;
  int                 fNumSavedHitGraphs;
  int                 fNumSavedWfmGraphs;
  TCanvas*            fTCanvas;
  TRandom2*           fRand;
  std::default_random_engine generator;
  TH1D*               hPMT_phelTimes[2];
  TH1D*               hPMT_phelTimes_electron[2];
  std::vector<float>  fPMT_wfm[2];
  std::vector<short>  fPMT_wfm_raw[2];
  std::vector<float>  fvbs[2];
 
  // ...........................................................
  // Vectors to store all the hit information
  std::vector<int>    fHitPlane;
  std::vector<int>    fHitWire;
  std::vector<float>  fHitT;
  std::vector<float>  fHitX;
  std::vector<float>  fHitW;
  std::vector<float>  fHitCharge;
  std::vector<float>  fHitIntegral;
  
  // ...........................................................
  // Photon visiblity and timing characterization libraries
  TFile* libraryFile;
  std::vector<std::vector< TH3D* >> libr_XYZ_visibility;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeMin;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeMean;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeRMS;
  
  // ...........................................................
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory truthDir          = tfs->mkdir("true");
  art::TFileDirectory graphDir_hits     = tfs->mkdir("graphs_hits");
  art::TFileDirectory graphDir_wfms     = tfs->mkdir("graphs_wfms");
  art::TFileDirectory diagDir           = tfs->mkdir("diagnostics");

};



//########################################################################################
MichelAnaFilter::MichelAnaFilter(fhicl::ParameterSet const & p)
: 
fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
{
  LOG_VERBATIM("MichelAnaFilter")
  <<"=======================================\n"
  <<"Configuring MichelAnaFilter..."; 

  // Read fhicl parameters
  this->reconfigure(p);
  
  // Reset counters
  fCachedRunNumber    = 0;
  fNumSavedHitGraphs  = 0;
  fNumSavedWfmGraphs  = 0;
  fNumOpticalID       = 0;
  fNumEvents          = 0;
  fNumEvents_StpTrk   = 0;
  fNumEvents_Shwr     = 0;
  fNumEvents_Shwr3D   = 0;
  fNumEvents_Calibration = 0;
  
  // Initialize data status ("real" data by default)
  fIsMC = false;
   
  // Random number generator
  fRand       = new TRandom2(fRandSeed);
  generator   .seed(fRandSeed);
    
  // Get a pointer to the geometry service provider
  fGeo = &*(art::ServiceHandle<geo::Geometry>());
  const geo::TPCGeo &tpc = fGeo->TPC(0);
  TVector3 frontFaceCenter = fGeo->GetTPCFrontFaceCenter(0,0);
  LOG_VERBATIM("MichelAnaFilter")
  <<"TPC front face center X: "<<frontFaceCenter.X()<<"\n"
  <<"TPC front face center Y: "<<frontFaceCenter.Y()<<"\n"
  <<"TPC front face center Z: "<<frontFaceCenter.Z()<<"\n"
  <<"TPC dimensions: dX = "<<fGeo->DetHalfWidth() * 2.<<"\n"
  <<"TPC dimensions: dY = "<<fGeo->DetHalfHeight() *2.<<"\n"
  <<"TPC dimensions: dZ = "<<fGeo->DetLength()<<"\n"
  <<"TPC min-max X : "<<tpc.MinX()<<" - "<<tpc.MaxX();
 
  // Initialize size of visiblity libraries
  libr_XYZ_visibility	  .resize(2);
  libr_XYZ_arrivalTimeMin .resize(2);
  libr_XYZ_arrivalTimeMean.resize(2);
  libr_XYZ_arrivalTimeRMS .resize(2);
  for( int i=0; i<2; i++){
    libr_XYZ_visibility[i].resize(2);
    libr_XYZ_arrivalTimeMin[i].resize(2);
    libr_XYZ_arrivalTimeMean[i].resize(2);
    libr_XYZ_arrivalTimeRMS[i].resize(2);
  }
 
  // Load the photon visibility and timing 3D libraries 
  if( fLibraryFile != "" ) {
    cet::search_path sp("FW_SEARCH_PATH");
    std::string fullname;
    sp.find_file(fLibraryFile.c_str(), fullname);
    if (fullname.empty()) {
      std::cout<<"VISIBILITY LIBRARY FILE NOT FOUND: "<<fLibraryFile.c_str()<<"\n";
      throw cet::exception("File not found");
    } else {
      LOG_VERBATIM("MichelAnaFilter")<<"Loading library file: "<<fLibraryFile;
      libraryFile                     = TFile::Open(fullname.c_str());
      libr_XYZ_visibility[0][0]	      = (TH3D*)libraryFile->Get("pmt0_XYZ_visibility_dir");
      libr_XYZ_visibility[0][1]	      = (TH3D*)libraryFile->Get("pmt0_XYZ_visibility_ref");
      libr_XYZ_arrivalTimeMin[0][0]   = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeMin_dir");
      libr_XYZ_arrivalTimeMin[0][1]   = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeMin_ref");
      libr_XYZ_arrivalTimeMean[0][0]  = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeMean_dir");
      libr_XYZ_arrivalTimeMean[0][1]  = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeMean_ref");
      libr_XYZ_arrivalTimeRMS[0][0]   = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeRMS_dir");
      libr_XYZ_arrivalTimeRMS[0][1]   = (TH3D*)libraryFile->Get("pmt0_XYZ_arrivalTimeRMS_ref");
      libr_XYZ_visibility[1][0]	      = (TH3D*)libraryFile->Get("pmt1_XYZ_visibility_dir");
      libr_XYZ_visibility[1][1]	      = (TH3D*)libraryFile->Get("pmt1_XYZ_visibility_ref");
      libr_XYZ_arrivalTimeMin[1][0]   = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeMin_dir");
      libr_XYZ_arrivalTimeMin[1][1]   = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeMin_ref");
      libr_XYZ_arrivalTimeMean[1][0]  = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeMean_dir");
      libr_XYZ_arrivalTimeMean[1][1]  = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeMean_ref");
      libr_XYZ_arrivalTimeRMS[1][0]   = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeRMS_dir");
      libr_XYZ_arrivalTimeRMS[1][1]   = (TH3D*)libraryFile->Get("pmt1_XYZ_arrivalTimeRMS_ref");
      libraryFile                     ->Close(); 
    }
  }
  
  // Reserve sizes for vector variables
  fMuContamCorr_EffTau  .reserve(2);
  fQE_ScaleFactor       .reserve(2);
  fPulses               .reserve(10);
  for(int i=0; i<2; i++) {
    fvHitTimes0[i]      .reserve(kMaxOpHits);
    fvHitTimes[i]       .reserve(kMaxOpHits);
    fvHitAmplitudes[i]  .reserve(kMaxOpHits);
    fvHitADC_100ns[i]   .reserve(kMaxOpHits);
    fvHitADC_200ns[i]   .reserve(kMaxOpHits);
    fvHitADC_300ns[i]   .reserve(kMaxOpHits);
    fvHitADC_400ns[i]   .reserve(kMaxOpHits);
    fvHitADC_500ns[i]   .reserve(kMaxOpHits);
    fvHitADC_600ns[i]   .reserve(kMaxOpHits);
    fvHitADC_700ns[i]   .reserve(kMaxOpHits);
    fvHitADC_900ns[i]   .reserve(kMaxOpHits);
    fvHitADC_1200ns[i]  .reserve(kMaxOpHits);
    fvHitADC_1500ns[i]  .reserve(kMaxOpHits);
    fvHitADC_1800ns[i]  .reserve(kMaxOpHits);
    fvHitADC_7000ns[i]  .reserve(kMaxOpHits);
    fvHitADCpf_100ns[i] .reserve(kMaxOpHits);
    fvHitADCpf_7000ns[i].reserve(kMaxOpHits);
    fvHitWidth[i]       .reserve(kMaxOpHits);
    fvPrepulseBaseline[i].reserve(kMaxOpHits);
    fvPrepulseRMS[i]    .reserve(kMaxOpHits);
    fvPrepulseSlowNorm[i].reserve(kMaxOpHits);
    fvPrepulseSlowTau[i].reserve(kMaxOpHits);
    fvPrepulseZeroPoint[i].reserve(kMaxOpHits);
    fvIsHitAtTrigger[i] .reserve(kMaxOpHits);
    fvIsHitSaturated[i] .reserve(kMaxOpHits);
    fvPrepulseX1[i]     .reserve(kMaxOpHits);
    
    // Initialize simulated arrival time histograms
    hPMT_phelTimes[i]           = new TH1D(Form("%i_phelTimes",i),Form("PE arrival times for PMT %i",i),28672,0.,28672.);
    hPMT_phelTimes_electron[i]  = new TH1D(Form("%i_phelTimes_electron",i),Form("PE arrival times for PMT %i",i),28672,0.,28672.);
    fPMT_wfm[i]                 .reserve(28672);
    fPMT_wfm_raw[i]             .reserve(28672);
    fvbs[i]                     .reserve(28672);
    
    // Initialize SPE
    fSPE[i]                     = 0.; 
    fSPE_err[i]                 = 0.;
    
    // PMT-specific counters
    numMichelWfms[i]            = 0;
    numMuSaturated[i]           = 0;
    numElSaturated[i]           = 0;
    aveWfmCounter_mu[i]         = 0;
    aveWfmCounter_el[i]         = 0;
    avePhelProfileCounter_mu[i] = 0;
    
  }
 
  // Define pulse integration windows
  fIntegrationWindows.resize(12);
  fIntegrationWindows[0]  = fPromptWindow;
  fIntegrationWindows[1]  = 200;
  fIntegrationWindows[2]  = 300;
  fIntegrationWindows[3]  = 400;
  fIntegrationWindows[4]  = 500;
  fIntegrationWindows[5]  = 600;
  fIntegrationWindows[6]  = 700;
  fIntegrationWindows[7]  = 900;
  fIntegrationWindows[8]  = 1200;
  fIntegrationWindows[9]  = 1500;
  fIntegrationWindows[10] = 1800;
  fIntegrationWindows[11] = fFullWindow;

  LOG_VERBATIM("MichelAnaFilter")
  <<"Configuration complete.\n"
  <<"=======================================";

}



//########################################################################################
void MichelAnaFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fRandSeed                 = p.get< int >                  ("RandSeed", 1989);
  fMuContamCorr_EffTau      = p.get< std::vector<float >>   ("MuContamCorr_EffTau",{0,0});
  fQE_ScaleFactor           = p.get< std::vector<float >>   ("QE_ScaleFactor",{1.,1.});
  fApplyTrigEffCut          = p.get< bool >                 ("ApplyTrigEffCut",false);
  fTrigEff_P                = p.get< std::vector< float >>  ("TrigEff_P",{80.,25.});
  fTrigEff_K                = p.get< std::vector< float >>  ("TrigEff_K",{8,8});
  fMuContamCorr             = p.get< bool >                 ("MuContamCorr",true);
  fMuContamCorr_dTmin       = p.get< float >                ("MuContamCorr_dTmin",2000);
  fMskBaselineSubtr         = p.get< bool >                 ("MskBaselineSubtr", false);
  fTruncateWfm              = p.get< short >                ("TruncateWfm",-1);
  fWfmSmoothingRange        = p.get< int  >                 ("WfmSmoothingRange",1);
  fCorrectOvershootMode     = p.get< std::string >          ("CorrectOvershootMode", "");
  fUseCrossingMuons         = p.get< bool >                 ("UseCrossingMuons",false);
  fLookAtTracks             = p.get< bool >                 ("LookAtTracks",true);
  fReq1StpTrk               = p.get< bool >                 ("Req1StpTrk",true);
  fFilter_PassAllEvents     = p.get< bool >                 ("Filter_PassAllEvents", false);
  fFilter_OpticalMode       = p.get< bool >                 ("Filter_OpticalMode", false); 
  fFilter_DecayTime         = p.get< std::vector<float >>   ("Filter_DecayTime",{0.,10000.});
  fFilter_Req3DShower       = p.get< bool >                 ("Filter_Req3DShower",false);
  fFilter_NumTrackStopping  = p.get< std::vector< int >>    ("Filter_NumTrackStopping",{1,1});
  fFilter_ElShowerEnergy    = p.get< std::vector< float >>  ("Filter_ElShowerEnergy",{-9999.,9999.});
  fSelectChannels           = p.get< std::vector< size_t >> ("SelectChannels", {1} );
  fPromptPECut              = p.get< std::vector< float >>  ("PromptPECut", {30.} );
  fMaxWfmRMS                = p.get< std::vector< float >>  ("MaxWfmRMS", {2.0} );
  fGradHitThresh            = p.get< std::vector< float >>  ("GradHitThresh",{-10.} );
  fMaxDiff_dT               = p.get< float >                ("MaxDiff_dT",15.);
  fdTcut                    = p.get< float >                ("dTcut",2000);
  fPromptWindow             = p.get< short >                ("PromptWindow",100);
  fFullWindow               = p.get< short >                ("FullWindow",7000);
  fBaselineWindowLength     = p.get< short >                ("BaselineWindowLength",1000); 
  fOpHitMatchThresh         = p.get< std::vector< float >>  ("OpHitMatchThresh",{300.,150.}); 
  fMinOpHitWidth            = p.get< std::vector< float >>  ("MinOpHitWidth",{8.,8.});
  fSPE_ScaleFactor          = p.get< std::vector< float >>  ("SPE_ScaleFactor",{1.,1.});
  fCollectionEff            = p.get< std::vector< float >>  ("CollectionEff",{1.,1.});
  fQE_ScaleFactor_Vis       = p.get< std::vector< float >>  ("QE_ScaleFactor_Vis",{1.,1.});
  fQE_ScaleFactor_VUV       = p.get< std::vector< float >>  ("QE_ScaleFactor_VUV",{1.,1.});
  fRequireHitMatching       = p.get< bool >                 ("RequireHitMatching",false);
  fDeleteNarrowHits         = p.get< bool >                 ("DeleteNarrowHits",false);
  fSmearSigT                = p.get< std::vector< float >>  ("SmearSigT",{0.02, 0.04 });
  fSmearSigPE               = p.get< std::vector< float >>  ("SmearSigPE",{0.50, 0.50 });
  fSmearFactor              = p.get< std::vector< float >>  ("SmearFactor",{0.,0.});
  fSmearFactorPrompt        = p.get< std::vector< float >>  ("SmearFactorPrompt",{0.,0.});
  fMvPerADC                 = p.get< float >                ("MvPerADC",0.195);
  fGateDelay                = p.get< float >                ("GateDelay",300);
  fMaxDecayTime             = p.get< float >                ("MaxDecayTime",7200);
  fTrackModule              = p.get< std::string >          ("TrackModule","pmtrack");
  fTrackCalModule           = p.get< std::string >          ("TrackCalModule","calo");
  fSimProducerModule        = p.get< std::string >          ("SimProducerModule","largeant");
  fHitsModule               = p.get< std::string >          ("HitsModule","gaushit");
  fHitsInstance               = p.get< std::string >          ("HitsInstance","");
  fLibraryFile              = p.get< std::string >          ("LibraryFile","");
  fFastLightRatio           = p.get< float >                ("FastLightRatio",0.30);
  fLateLightTau             = p.get< float >                ("LateLightTau",1300);
  fTauT                     = p.get< float >                ("TauT",1300);
  fTauS                     = p.get< float >                ("TauS",6.5);
  fCalAreaConstants         = p.get< std::vector< float >>  ("CalAreaConstants", {0.022, 0.050});
  fFiducialMarginX          = p.get< float >                ("FiducialMarginX",2.);
  fFiducialMarginY          = p.get< float >                ("FiducialMarginY",2.);
  fFiducialMarginZ          = p.get< float >                ("FiducialMarginZ",2.);
  fMuTrackLengthMin         = p.get< float >                ("MuTrackLengthMin",5.);
  fMaxHitSeparation         = p.get< float >                ("MaxHitSeparation",3.);
  fTruncMeanWindow          = p.get< int >                  ("TruncMeanWindow",5);
  fTruncMeanP               = p.get< float >                ("TruncMeanP",0.25);
  fMaxMatchDistance         = p.get< int >                  ("MaxMatchDistance",8);
  fLocalLinearityWindow     = p.get< int >                  ("LocalLinearityWindow",5);
  fLinTol                   = p.get< float >                ("LinTol",3.);
  fLinThresh                = p.get< float >                ("LinThresh",0.8);
  fMinClusterSize           = p.get< int >                  ("MinClusterSize",10);
  fMinElShowerFrac          = p.get< float >                ("MinElShowerFrac",0.);
  fMinElClusterSize         = p.get< int >                  ("MinElClusterSize",3);
  fMaxElClusterSize         = p.get< int >                  ("MaxElClusterSize",50);
  fMaxExtraHits             = p.get< int >                  ("MaxExtraHits",0);
  fMuEndRadius              = p.get< float >                ("MuEndRadius",3.);
  fMinDecayAngle2D          = p.get< float >                ("DecayAngle2DMin",10.);
  fMaxDecayAngle2D          = p.get< float >                ("DecayAngle2DMax",170.);
  fShowerAcceptanceAngle    = p.get< float >                ("ShowerAcceptanceAngle",60.);
  fWph                      = p.get< float >                ("Wph", 19.5 );
  fWion                     = p.get< float >                ("Wion", 23.6 );
  fExcRatio                 = p.get< float >                ("ExcRatio", 0.21 ); // Aprile 
  fRecombFactor             = p.get< float >                ("RecombFactor", 0.645);
  fRecombFactorTrack        = p.get< float >                ("RecombFactorTrack", 0.688);
  fRecombFactorShower       = p.get< float >                ("RecombFactorShower", 0.576);
  fMaxHitTimeDiff           = p.get< float >                ("MaxHitTimeDiff", 1.5 );
  fMinNumPts3D              = p.get< int >                  ("MinNumPts3D",2);
  fMinFracHits3D            = p.get< float >                ("MinFracHits3D",0.3);
  fLinThreshMuFit           = p.get< float >                ("LinThreshMuFit",0.8);
  fMinMuClusterHitsEndFit   = p.get< int >                  ("MinMuClusterHitsEndFit",5);
  fMinMuLinearity           = p.get< float >                ("MinMuLinearity",0.9);
  fMinMuClusterSize         = p.get< int >                  ("MinMuClusterSize",8);
  fMinFracMuHitsLinear      = p.get< float >                ("MinFracMuHitsLinear",0.5);  
  fBndOffset                = p.get< int >                  ("BndOffset",0);
  fMaxSavedHitGraphs        = p.get< int >                  ("MaxSavedHitGraphs",0);
  fMaxSavedWfmGraphs        = p.get< int >                  ("MaxSavedWfmGraphs",0);
}




//########################################################################################
void MichelAnaFilter::beginJob()
{
  
  // ====================================================================
  // Define TTree and set all branches
  fTree= tfs->make<TTree>("anatree","anatree");
  // ...event ID info, calibration values used
  fTree->Branch("RunNumber",                &fRunNumber,              "RunNumber/I");
  fTree->Branch("SubRunNumber",             &fSubRunNumber,           "SubRunNumber/I");
  fTree->Branch("EventNumber",              &fEventNumber,            "EventNumber/I");
  fTree->Branch("EventTime",                &fEventTime,              "EventTime/I");
  fTree->Branch("MichelOpticalID",          &fMichelOpticalID,        "MichelOpticalID/O");
  fTree->Branch("ElectronLifetime",         &fElectronLifetime,       "ElectronLifetime/F");
  fTree->Branch("MuContamCorr_EffTau",      &fMuContamCorr_EffTau);
  fTree->Branch("QE_ScaleFactor",           &fQE_ScaleFactor);
  //...track information
  fTree->Branch("NumTrack",                 &fNumTrack,               "NumTrack/I");
  fTree->Branch("NumTrackStopping",         &fNumTrackStopping,       "NumTrackStopping/I");
  fTree->Branch("MuTrackLength",            &fMuTrackLength,          "MuTrackLength/F");
  fTree->Branch("MuTrackEnergy",            &fMuTrackEnergy,          "MuTrackEnergy/F");
  fTree->Branch("MuTrackZenithAngle",       &fMuTrackZenithAngle,     "MuTrackZenithAngle/F");
  fTree->Branch("CrossingMuVertex_X",         &fCrossingMuVertex_X,       "CrossingMuVertex_X/F");
  fTree->Branch("CrossingMuVertex_Y",         &fCrossingMuVertex_Y,       "CrossingMuVertex_Y/F");
  fTree->Branch("CrossingMuVertex_Z",         &fCrossingMuVertex_Z,       "CrossingMuVertex_Z/F");
  fTree->Branch("CrossingMuEnd_X",         &fCrossingMuEnd_X,       "CrossingMuEnd_X/F");
  fTree->Branch("CrossingMuEnd_Y",         &fCrossingMuEnd_Y,       "CrossingMuEnd_Y/F");
  fTree->Branch("CrossingMuEnd_Z",         &fCrossingMuEnd_Z,       "CrossingMuEnd_Z/F");
  fTree->Branch("CrossingMuLength",         &fCrossingMuLength,       "CrossingMuLength/F");
  fTree->Branch("CrossingMuCharge",         &fCrossingMuCharge,       "CrossingMuCharge/F");
  fTree->Branch("CrossingMuPhotons",          &fCrossingMuPhotons,        "CrossingMuPhotons/F");
  fTree->Branch("CrossingMuPhotonsPrompt",  &fCrossingMuPhotonsPrompt,"CrossingMuPhotonsPrompt/F");
  // ...optical information
  fTree->Branch("NumOpHits0",               &fNumOpHits0,             "NumOpHits0[2]/I");
  fTree->Branch("NumOpHits",                &fNumOpHits,              "NumOpHits[2]/I");
  fTree->Branch("SPE",                      &fSPE,                    "SPE[2]/F");
  fTree->Branch("SPE_err",                  &fSPE_err,                "SPE_err[2]/F");
  fTree->Branch("dT",                       &fdT,                     "dT[2]/F");
  fTree->Branch("Amplitude",                &fAmplitude,              "Amplitude[2]/F");
  fTree->Branch("Width",                    &fWidth,                  "Width[2]/F");
  fTree->Branch("MuWidth",                  &fMuWidth,                "MuWidth[2]/F");
  fTree->Branch("PE_prompt",                &fPE_prompt,              "PE_prompt[2]/F");
  fTree->Branch("PE_total",                 &fPE_total,               "PE_total[2]/F");
  fTree->Branch("PE_promptRaw",             &fPE_promptRaw,           "PE_promptRaw[2]/F");
  fTree->Branch("PE_totalRaw",              &fPE_totalRaw,            "PE_totalRaw[2]/F");
  fTree->Branch("PromptFrac",               &fPromptFrac,             "PromptFrac[2]/F");
  fTree->Branch("MuAmplitude",              &fMuAmplitude,            "MuAmplitude[2]/F");
  fTree->Branch("MuPE_prompt",              &fMuPE_prompt,            "MuPE_prompt[2]/F");
  fTree->Branch("MuPulseSaturated",         &fMuPulseSaturated,       "MuPulseSaturated[2]/O");
  fTree->Branch("ElPulseSaturated",         &fElPulseSaturated,       "ElPulseSaturated[2]/O");
  fTree->Branch("DecayTime",                &fDecayTime,              "DecayTime/F");
  //...clustering/showering information
  fTree->Branch("NumPlaneHits",             &fNumPlaneHits,           "NumPlaneHits[2]/I");
  fTree->Branch("ElOffsetT",                &fElOffsetT,              "ElOffsetT/F");
  fTree->Branch("AveDriftTime",             &fAveDriftTime,           "AveDriftTime/F");
  fTree->Branch("AveX",                     &fAveX,                   "AveX/F");
  fTree->Branch("MinCovAtBnd",              &fMinCovAtBnd,            "MinCovAtBnd/F"); 
  fTree->Branch("ClusterSize",              &fClusterSize,            "ClusterSize/I"); 
  fTree->Branch("ElClusterSize",            &fElClusterSize,          "ElClusterSize/I");
  fTree->Branch("MuClusterSize",            &fMuClusterSize,          "MuClusterSize/I");
  fTree->Branch("ExtraHits",                &fExtraHits,              "ExtraHits/I");
  fTree->Branch("MuAveLinearity",           &fMuAveLinearity,         "MuAveLinearity/F");
  fTree->Branch("FracMuHitsLinear",         &fFracMuHitsLinear,       "FracMuHitsLinear/F");
  fTree->Branch("MuClusterHitsEndFit",      &fMuClusterHitsEndFit,    "MuClusterHitsEndFit/I");
  fTree->Branch("DecayAngle2D",             &fDecayAngle2D,           "DecayAngle2D/F");
  fTree->Branch("ElShowerAngleMean",        &fElShowerAngleMean,      "ElShowerAngleMean/F");
  fTree->Branch("ElShowerAngleRMS",         &fElShowerAngleRMS,       "ElShowerAngleRMS/F");
  fTree->Branch("ElShowerSize",             &fElShowerSize,           "ElShowerSize/I");
  fTree->Branch("ElShowerFrac",             &fElShowerFrac,           "ElShowerFrac/F");
  fTree->Branch("TotalCharge",              &fTotalCharge,            "TotalCharge/F");
  fTree->Branch("ElShowerCharge",           &fElShowerCharge,         "ElShowerCharge/F");
  fTree->Branch("ElShowerEnergy",           &fElShowerEnergy,         "ElShowerEnergy/F");
  fTree->Branch("ElTrackCharge",            &fElTrackCharge,          "ElTrackCharge/F");
  fTree->Branch("ElTrackEnergy",            &fElTrackEnergy,          "ElTrackEnergy/F");
  fTree->Branch("MinCovAtBnd_Pl0",          &fMinCovAtBnd_Pl0,        "MinCovAtBnd_Pl0/F"); 
  fTree->Branch("ClusterSize_Pl0",          &fClusterSize_Pl0,        "ClusterSize_Pl0/I"); 
  fTree->Branch("ElClusterSize_Pl0",        &fElClusterSize_Pl0,      "ElClusterSize_Pl0/I");
  fTree->Branch("MuClusterSize_Pl0",        &fMuClusterSize_Pl0,      "MuClusterSize_Pl0/I");
  fTree->Branch("ExtraHits_Pl0",            &fExtraHits_Pl0,          "ExtraHits_Pl0/I");
  fTree->Branch("MuAveLinearity_Pl0",       &fMuAveLinearity_Pl0,     "MuAveLinearity_Pl0/F");
  fTree->Branch("FracMuHitsLinear_Pl0",     &fFracMuHitsLinear_Pl0,   "FracMuHitsLinear_Pl0/F");
  fTree->Branch("MuClusterHitsEndFit_Pl0",  &fMuClusterHitsEndFit_Pl0,"MuClusterHitsEndFit_Pl0/I");
  fTree->Branch("DecayAngle2D_Pl0",         &fDecayAngle2D_Pl0,       "DecayAngle2D_Pl0/F");
  fTree->Branch("ElShowerAngleMean_Pl0",    &fElShowerAngleMean_Pl0,  "ElShowerAngleMean_Pl0/F");
  fTree->Branch("ElShowerAngleRMS_Pl0",     &fElShowerAngleRMS_Pl0,   "ElShowerAngleRMS_Pl0/F");
  fTree->Branch("ElShowerSize_Pl0",         &fElShowerSize_Pl0,       "ElShowerSize_Pl0/I");
  fTree->Branch("ElShowerFrac_Pl0",         &fElShowerFrac_Pl0,       "ElShowerFrac_Pl0/F");
  fTree->Branch("TotalCharge_Pl0",          &fTotalCharge_Pl0,        "TotalCharge_Pl0/F");
  fTree->Branch("ElShowerCharge_Pl0",       &fElShowerCharge_Pl0,     "ElShowerCharge_Pl0/F");
  fTree->Branch("ElShowerEnergy_Pl0",       &fElShowerEnergy_Pl0,     "ElShowerEnergy_Pl0/F");
  fTree->Branch("ElTrackCharge_Pl0",        &fElTrackCharge_Pl0,      "ElTrackCharge_Pl0/F");
  fTree->Branch("ElTrackEnergy_Pl0",        &fElTrackEnergy_Pl0,      "ElTrackEnergy_Pl0/F");
  fTree->Branch("MuEndHitTimeDiff",         &fMuEndHitTimeDiff,       "MuEndHitTimeDiff/F");
  fTree->Branch("ElShowerVis",              &fElShowerVis,            "ElShowerVis/F");
  fTree->Branch("ElShowerVisCh",            &fElShowerVisCh,          "ElShowerVisCh[2]/F");
  fTree->Branch("ElShowerPhotons",          &fElShowerPhotons,        "ElShowerPhotons/F");
  fTree->Branch("ElTrackdEdx",              &fElTrackdEdx,            "ElTrackdEdx/F");
  fTree->Branch("MuEnd3D_X",                &fMuEnd3D_X,              "MuEnd3D_X/F");
  fTree->Branch("MuEnd3D_Y",                &fMuEnd3D_Y,              "MuEnd3D_Y/F");
  fTree->Branch("MuEnd3D_Z",                &fMuEnd3D_Z,              "MuEnd3D_Z/F");
  fTree->Branch("ElShowerCentroid_X",       &fElShowerCentroid_X,     "ElShowerCentroid_X/F");
  fTree->Branch("ElShowerCentroid_Y",       &fElShowerCentroid_Y,     "ElShowerCentroid_Y/F");
  fTree->Branch("ElShowerCentroid_Z",       &fElShowerCentroid_Z,     "ElShowerCentroid_Z/F");
  fTree->Branch("FracHits3D",               &fFracHits3D,             "FracHits3D/F");
  fTree->Branch("NumPts3D",                 &fNumPts3D,               "NumPts3D/I");
  fTree->Branch("NumPts3DTrk",              &fNumPts3DTrk,            "NumPts3DTrk/I");
  //...truth information
  fTree->Branch("True_MuCharge",            &fTrue_MuCharge,          "True_MuCharge/I");
  fTree->Branch("True_MuTrackEnd_X",        &fTrue_MuTrackEnd_X,      "True_MuTrackEnd_X/F");
  fTree->Branch("True_MuTrackEnd_Y",        &fTrue_MuTrackEnd_Y,      "True_MuTrackEnd_Y/F");
  fTree->Branch("True_MuTrackEnd_Z",        &fTrue_MuTrackEnd_Z,      "True_MuTrackEnd_Z/F");
  fTree->Branch("True_ElAngle",             &fTrue_ElAngle,           "True_ElAngle/F");
  fTree->Branch("True_ElEnergy",            &fTrue_ElEnergy,          "True_ElEnergy/F");
  fTree->Branch("True_ElShowerEnergyDep",   &fTrue_ElShowerEnergyDep, "True_ElShowerEnergyDep/F");
  fTree->Branch("True_ElShowerChargeDep",   &fTrue_ElShowerChargeDep, "True_ElShowerChargeDep/F");
  fTree->Branch("True_ElTrackChargeDep",    &fTrue_ElTrackChargeDep,  "True_ElTrackChargeDep/F");
  fTree->Branch("True_ElShowerVis",       &fTrue_ElShowerVis,     "True_ElShowerVis/F");
  fTree->Branch("True_ElShowerVisCh",       &fTrue_ElShowerVisCh,     "True_ElShowerVisCh[2]/F");
  fTree->Branch("True_TotalEnergyDep",      &fTrue_TotalEnergyDep,    "True_TotalEnergyDep/F");
  fTree->Branch("True_IsElContained",       &fTrue_IsElContained,     "True_IsElContained/O");
  fTree->Branch("True_ElShowerPhotons",     &fTrue_ElShowerPhotons,   "True_ElShowerPhotons/F");
  fTree->Branch("True_ElTrackPhotons",      &fTrue_ElTrackPhotons,    "True_ElTrackPhotons/F");
  fTree->Branch("True_ElShowerPhotonsPrompt",&fTrue_ElShowerPhotonsPrompt,"True_ElShowerPhotonsPrompt/F");
  fTree->Branch("True_ElShowerPhotonsLate", &fTrue_ElShowerPhotonsLate,"True_ElShowerPhotonsLate/F");
  fTree->Branch("True_dT",                  &fTrue_dT,                "True_dT/F");
  fTree->Branch("True_PE_prompt",           &fTrue_PE_prompt,         "True_PE_prompt[2]/F");
  fTree->Branch("True_PE_total",            &fTrue_PE_total,          "True_PE_total[2]/F");
  fTree->Branch("True_MuPE_prompt",         &fTrue_MuPE_prompt,       "True_MuPE_prompt[2]/F");
 
  // ====================================================================
  // Light-based event cuts (won't really apply to MC until we have wfms simulated)
  hOpEventCuts=tfs->make<TH1D>("EventCutsLight","Event reduction (opt. reco + Michel ID)",9,0,9);
  hOpEventCuts->SetOption("HIST TEXT");
  hOpEventCuts->GetXaxis()->SetBinLabel(1,"Total evts");          // total events
  hOpEventCuts->GetXaxis()->SetBinLabel(2,"All opdets present");  // all opdets present
  hOpEventCuts->GetXaxis()->SetBinLabel(3,"Wfm RMS");             // wfms pass RMS cut
  hOpEventCuts->GetXaxis()->SetBinLabel(4,"2 hits");              // 2 hits in both PMTs
  hOpEventCuts->GetXaxis()->SetBinLabel(5,"#DeltaT match");       // dT match btwn. PMTs
  hOpEventCuts->GetXaxis()->SetBinLabel(6,"Hit widths");          // hit width cuts
  hOpEventCuts->GetXaxis()->SetBinLabel(7,"2nd hit time");        // 2nd hit at trigger sample
  hOpEventCuts->GetXaxis()->SetBinLabel(8,"2nd hit sat.");        // 2nd hit not saturated
  hOpEventCuts->GetXaxis()->SetBinLabel(9,"Min #DeltaT");         // dT > gate delay
  
  // ===================================================================
  // General event cuts
  hEventCuts=tfs->make<TH1D>("EventCuts","Event reduction",8,0,8);
  hEventCuts->SetOption("HIST TEXT");
  hEventCuts->GetXaxis()->SetBinLabel(1, "Total evts");       // total events
  hEventCuts->GetXaxis()->SetBinLabel(2, "1 stp trk");        // 1 stopping trk
  hEventCuts->GetXaxis()->SetBinLabel(3, "Optical ID");       // optical Michel ID
  hEventCuts->GetXaxis()->SetBinLabel(4, "Clstr bnd.");       // cluster boundary pt found
  hEventCuts->GetXaxis()->SetBinLabel(5, "Clstr size");       // cluster size cuts
  hEventCuts->GetXaxis()->SetBinLabel(6, "Extra hits");       // no extra hits around mu endpt
  hEventCuts->GetXaxis()->SetBinLabel(7, "#mu quality");      // mu linearity cuts
  hEventCuts->GetXaxis()->SetBinLabel(8, "2D decay angle");   // 2D decay angle range cut
 
  // ===================================================================
  // Electron lifetime histograms 
  hElectronLifetime       = diagDir.make<TH1D>("ElectronLifetime","Electron lifetime from database;Electron lifetime [#mus]",2500,-100.,2400.);
  hElectronLifetimeReco   = diagDir.make<TH1D>("ElectronLifetimeReco","Electron lifetime as measured from passing tracks;Electron lifetime [#mus]",100,0.,2000.);
  hdQdxVsT         = diagDir.make<TH1D>("dQdxVsT","dQ/dx vs. drift time for single crossing trks;Drift time [#mus];dQ/dx [ADC/cm]",30,10.,310.);
  hdQdxVsT_N       = diagDir.make<TH1D>("dQdxVsT_N",";Drift time [#mus];Number of entries",30,10.,310.);

  // ==================================================================
  // Effective tau (from DB)
  hEffTau[0]              = diagDir.make<TH1D>("0_EffTau","OpDet0: Effective late-light lifetime from DB table;#tau [ns]",70,0,1400.);
  hEffTau[1]              = diagDir.make<TH1D>("1_EffTau","OpDet1: Effective late-light lifetime from DB table;#tau [ns]",70,0,1400.);
 
  // ==================================================================
  // Light plots for the PMTs
  
  // Photoelectron and Michel energy ranges
  float x1_S[2]     = {   0.,     0.    };
  float x2_S[2]     = {   1500.,  500.  };
  int nbins_S[2]    = {   50,     50    };
  float x1_S100[2]    = {   0.,    0.    };
  float x2_S100[2]    = {   400.,  140.  };
  int nbins_S100[2]   = {   40,    35    };
  float x1_michelEnergy   = 0.;
  float x2_michelEnergy   = 120.;
  int Nbins_michelEnergy  = 120;
  
  hPE       = tfs->make<TH1D>("PE","Total Michel light for event;Michel light collected [PEs]",200,0,2000);
  
  // ----------------------------------------------------------------------
  // Loop over the PMTs and make PMT-specific histograms
  for(size_t i=0; i<fSelectChannels.size(); i++){
    size_t  ch     = fSelectChannels.at(i);
 
    // ------------------------------------------------------------
    // Average waveforms and binned photon arrival time profiles
    sprintf(histName,"%lu_AveWfm_mu",ch);
    sprintf(histTitle,"Averaged waveform of #mu pulse, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm_mu[ch]  = diagDir.make<TH1D>(histName,histTitle,4400,-400.,4000.);
    sprintf(histName,"%lu_AveWfm_el",ch);
    sprintf(histTitle,"Averaged waveform of Michel pulse, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm_el[ch]  = diagDir.make<TH1D>(histName,histTitle,7700,-700.,7000.);
    float bins[] = { 0, 100, 300, 500, 700, 900, 1200, 1500, 1800};
    int binnum = sizeof(bins)/sizeof(float) - 1;
    sprintf(histName,"%lu_AvePhelProfile_mu",ch);
    sprintf(histTitle,"Averaged photon arrival time profile of #mu pulse, opdet %lu;Time after #mu pulse [ns];PEs per ns",ch);
    hAvePhelProfile_mu[ch] = diagDir.make<TH1D>(histName,histTitle, binnum, bins);

    // ------------------------------------------------------------
    // PMT diagnostic plots
    sprintf(histName,"%lu_WfmRMS",ch);
    sprintf(histTitle,"Waveform RMS, opdet %lu;ADC",ch);
    hWfmRMS[ch]   = diagDir.make<TH1D>(histName,histTitle,100,0.,5.);
    sprintf(histName,"%lu_BaselinePE_100ns",ch);
    sprintf(histTitle,"Integrated baseline region (100ns), opdet %lu;PE",ch);
    hBaselinePE_100ns[ch] = diagDir.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_BaselinePE_500ns",ch);
    sprintf(histTitle,"Integrated baseline region (500ns), opdet %lu;PE",ch);
    hBaselinePE_500ns[ch] = diagDir.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_BaselinePE_1000ns",ch);
    sprintf(histTitle,"Integrated baseline region (1000ns), opdet %lu;PE",ch);
    hBaselinePE_1000ns[ch] = diagDir.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_NumOpHits",ch);
    sprintf(histTitle,"Opdet %lu;Number of optical hits",ch);
    hNumOpHits[ch]   = diagDir.make<TH1D>(histName,histTitle,20,0,20);
    sprintf(histName,"%lu_MuTau",ch);
    sprintf(histTitle,"Opdet %lu, Lifetime from #mu late-light fit (for contamination correction;ns",ch);
    hMuTau[ch]     = diagDir.make<TH1D>(histName,histTitle, 250, 400., 5400.);
   
    // --------------------------------------------------------------
    // True PMT scintillation visibilities 
    sprintf(histName,"%lu_True_Vis",ch);
    sprintf(histTitle,"Opdet %lu, total scintillation visiblity;Fractional visibility",ch);
    hTrue_VisTotal[ch]     = truthDir.make<TH1D>(histName,histTitle,200,-0.0001,0.0019);
    sprintf(histName,"%lu_True_Vis_dir",ch);
    sprintf(histTitle,"Opdet %lu, scintillation visiblity, direct light;Fractional visibility",ch);
    hTrue_Vis[ch][0]     = truthDir.make<TH1D>(histName,histTitle,200,-0.0001,0.0019);
    sprintf(histName,"%lu_True_Vis_ref",ch);
    sprintf(histTitle,"Opdet %lu, scintillation visiblity, reflected light;Fracitonal visibility",ch);
    hTrue_Vis[ch][1]     = truthDir.make<TH1D>(histName,histTitle,200,-0.0001,0.0019);
   
    // -------------------------------------------------------------
    // dT + decay time diagnostics 
    sprintf(histName,"%lu_OpHitTime",ch);
    sprintf(histTitle,"Opdet %lu;Time of op hits rel. to trigger [ns]",ch);
    hOpHitTime[ch]          = diagDir.make<TH1D>(histName,histTitle,200,-1000,1000);
    sprintf(histName,"%lu_dT",ch);
    sprintf(histTitle,"Opdet %lu;#DeltaT [ns]",ch);
    hdT[ch]          = tfs->make<TH1D>(histName,histTitle,750,0.,7500.);
    hdT_vs_Amp[ch]   = diagDir.make<TH2D>(Form("%lu_dT_vs_Amp",ch),Form("Opdet %lu;#DeltaT [ns];Amplitude [ADC]",ch),150,0.,7500.,102,-20,1000);
    hdT_vs_Amp[ch]  ->SetOption("colz");
    hdT_vs_PE_totalRaw[ch]  = diagDir.make<TH2D>(Form("%lu_dT_vs_PE_totalRaw",ch),Form("Opdet %lu (raw integral);#Delta T [ns];Total light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S[ch], x1_S[ch], x2_S[ch]);
    hdT_vs_PE_total[ch]     = diagDir.make<TH2D>(Form("%lu_dT_vs_PE_total",ch),Form("Opdet %lu;#Delta T [ns];Total light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S[ch], x1_S[ch], x2_S[ch]);
    hdT_vs_PE_promptRaw[ch]  = diagDir.make<TH2D>(Form("%lu_dT_vs_PE_promptRaw",ch),Form("Opdet %lu (raw integral);#Delta T [ns];Prompt light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S100[ch], x1_S100[ch], x2_S100[ch]);
    hdT_vs_PE_prompt[ch]     = diagDir.make<TH2D>(Form("%lu_dT_vs_PE_prompt",ch),Form("Opdet %lu;#Delta T [ns];Prompt light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S100[ch], x1_S100[ch], x2_S100[ch]);
    hdT_vs_PE_totalRaw[ch]  ->SetOption("colz");
    hdT_vs_PE_total[ch]     ->SetOption("colz");
    hdT_vs_PE_promptRaw[ch]  ->SetOption("colz");
    hdT_vs_PE_prompt[ch]     ->SetOption("colz");
    
    // -------------------------------------------------------------
    // General PMT histograms
    
    sprintf(histName,"%lu_PromptFracRaw",ch);
    sprintf(histTitle,"Opdet %lu (raw integral);Prompt fraction of Michel candidate pulse",ch);
    hPromptFracRaw[ch]     = diagDir.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFrac",ch);
    sprintf(histTitle,"Opdet %lu;Prompt fraction of Michel candidate pulse",ch);
    hPromptFrac[ch]     = diagDir.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFracRaw_dTcut",ch);
    sprintf(histTitle,"Opdet %lu (raw integral);Prompt fraction of Michel candidate pulse",ch);
    hPromptFracRaw_dTcut[ch]     = diagDir.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFrac_dTcut",ch);
    sprintf(histTitle,"Opdet %lu;Prompt fraction of Michel candidate pulse",ch);
    hPromptFrac_dTcut[ch]     = diagDir.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    
    sprintf(histName,"%lu_OpHitAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of ALL found hits [ADC]",ch);
    hOpHitAmplitude[ch]          = diagDir.make<TH1D>(histName,histTitle,120,0.,1200.);
    sprintf(histName,"%lu_OpHitAmplitude_unmatched",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of hits not matched to another PMT [ADC]",ch);
    hOpHitAmplitude_unmatched[ch]          = diagDir.make<TH1D>(histName,histTitle,120,0.,1200.);
    sprintf(histName,"%lu_Amplitude",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of Michel candidate pulse [ADC]",ch);
    hAmplitude[ch]          = diagDir.make<TH1D>(histName,histTitle,120,0.,1200.);
    
    sprintf(histName,"%lu_HitPromptPE",ch);
    sprintf(histTitle,"Opdet %lu;PromptPE of ALL found hits [pe]",ch);
    hHitPromptPE[ch]          = diagDir.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_HitPromptPE_unmatched",ch);
    sprintf(histTitle,"Opdet %lu;PromptPE of hits not matched to another PMT [pe]",ch);
    hHitPromptPE_unmatched[ch]          = diagDir.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    
    sprintf(histName,"%lu_Width",ch);
    sprintf(histTitle,"Opdet %lu;Width of Michel candidate pulse [ns]",ch);
    hWidth[ch]          = diagDir.make<TH1D>(histName,histTitle,120,0,30);
    sprintf(histName,"%lu_Width_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Width of all hit pulses [ns]",ch);
    hWidth_AllHits[ch]          = diagDir.make<TH1D>(histName,histTitle,120,0,30);
    
    sprintf(histName,"%lu_AmplitudeVsWidth",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of Michel candidate pulse [mV];Width [ns]",ch);
    hAmplitudeVsWidth[ch]          = diagDir.make<TH2D>(histName,histTitle,200,0.,200.,120,0,30);
    hAmplitudeVsWidth[ch]->SetOption("colz");
    
    sprintf(histName,"%lu_PromptPEVsWidth",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PEs];Pulse width [ns]",ch);
    hPromptPEVsWidth[ch]          = diagDir.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,30);
    hPromptPEVsWidth[ch]          ->SetOption("colz");
    
    sprintf(histName,"%lu_AmplitudeVsWidth_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude [ADC];Width [ns]",ch);
    hAmplitudeVsWidth_AllHits[ch]          = diagDir.make<TH2D>(histName,histTitle,120,0.,1200.,120,0,30);
    hAmplitudeVsWidth_AllHits[ch]->SetOption("colz");
    
    sprintf(histName,"%lu_PromptPEVsWidth_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light [PEs];Pulse width [ns]",ch);
    hPromptPEVsWidth_AllHits[ch]          = diagDir.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,30);
    hPromptPEVsWidth_AllHits[ch]          ->SetOption("colz");
    sprintf(histName,"%lu_PromptPEVsAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PEs];Pulse amplitude [ADC]",ch);
    hPromptPEVsAmplitude[ch]          = diagDir.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,600.);
    hPromptPEVsAmplitude[ch]          ->SetOption("colz");
    
    sprintf(histName,"%lu_MuPromptPEVsAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of muon candidate pulse [PEs];Pulse amplitude [ADC]",ch);
    hMuPromptPEVsAmplitude[ch]          = diagDir.make<TH2D>(histName,histTitle,100,0.,500.,120,0,1200.);
    hMuPromptPEVsAmplitude[ch]          ->SetOption("colz");
   
    sprintf(histName,"%lu_PE_totalRaw_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu (raw integral), #DeltaT > %3.1f#mus, Michel shwr reco'd;Total light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_totalRaw_dTcut_shwr[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch]);
    
    sprintf(histName,"%lu_PE_prompt",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PE]",ch);
    hPE_prompt[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_PE_prompt_dTcut",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus;Prompt light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_prompt_dTcut[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_PE_prompt_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd;Prompt light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_prompt_dTcut_shwr[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    
    sprintf(histName,"%lu_PE_total",ch);
    sprintf(histTitle,"Opdet %lu;Total light of Michel candidate pulse [PE]",ch);
    hPE_total[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch]);
    sprintf(histName,"%lu_PE_total_dTcut",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus;Total light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_total_dTcut[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch]);
    sprintf(histName,"%lu_PE_total_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd;Total light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_total_dTcut_shwr[ch]     = tfs->make<TH1D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch]);
   
    sprintf(histName,"%lu_PE_total_TrueVsReco",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd;True Light [PE];Reco Light [PE]",ch,fdTcut/1000);
    hTrue_PE_total_TrueVsReco[ch] = diagDir.make<TH2D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch],nbins_S[ch],x1_S[ch],x2_S[ch]);
    hTrue_PE_total_TrueVsReco[ch]->SetOption("colz");
    
    // -------------------------------------------------------------
    // "Truth" PMT histograms 
    sprintf(histName,"%lu_True_PE_preTrig",ch);
    sprintf(histTitle,"Opdet %lu: True number of detected photons (prior to trigger efficiency cut);True detected light [PE]",ch);
    hTrue_PE_preTrig[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_True_PE",ch);
    sprintf(histTitle,"Opdet %lu: True number of detected photons;True detected light [PE]",ch);
    hTrue_PE[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_TrueLY",ch);
    sprintf(histTitle,"True LY (opdet %lu);Light Yield [PE/MeV]",ch);
    hTrue_LightYield[ch]        = truthDir.make<TH1D>(histName,histTitle,200,0.,20.); 
    sprintf(histName,"%lu_True_PE_prompt",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PE]",ch);
    hTrue_PE_prompt[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_True_PE_total",ch);
    sprintf(histTitle,"Opdet %lu;Total light of Michel candidate pulse [PE]",ch);
    hTrue_PE_total[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_True_PE_total_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd (true);True detected light of Michel pulse [PE]",ch,fdTcut/1000);
    hTrue_PE_total_dTcut_shwr[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_True_PERes",ch);
    sprintf(histTitle,"Opdet %lu (#DeltaT > %3.1f#mus, shwr), Photoelectron resolution of Michel candidate pulse;(reco-true)/true",ch,fdTcut/1000);
    hTrue_PERes[ch]       = truthDir.make<TH1D>(histName,histTitle,150,-1.5,1.5);

  }
 
  // ===================================================================
  // Decay time histogram (after event cuts)
  hDecayTime              = tfs->make<TH1D>("DecayTime","Measured #mu decay time (prompt PE cut);#DeltaT [ns]",375,0.,7500.);
 
  // ==================================================================
  // Comparison plots between the two PMTs
  hNumOpHitsCompare   = diagDir.make<TH2D>("NumOpHitsCompare",";Num optical hits in PMT 0;Num optical hits in PMT 1",10,0.,10.,10,0.,10.);
  hNumOpHitsCompare   ->SetOption("colz");
  hPECompare          = diagDir.make<TH2D>("PECompare",";Total Michel candidate light in PMT 0;Total Michel candidate light in PMT 1",nbins_S[0],x1_S[0],x2_S[0],nbins_S[1],x1_S[1],x2_S[1]);
  hPECompare          ->SetOption("colz");
  hDiff_dT            = diagDir.make<TH1D>("Diff_dT","Decay time difference between PMTs;|#DeltaT_{1} - #DeltaT_{2}| [ns]",100,0.,100.);

  // ===================================================================
  // Wire charge-hit diagnostics  
  hHitX[0]              = diagDir.make<TH1D>("HitX_Ind","Calculated drift position from hit time (induction plane);X [cm]",300,-5.,55.);
  hHitX[1]              = diagDir.make<TH1D>("HitX_Col","Calculated drift position from hit time (collection plane);X [cm]",300,-5.,55.);
  hHitRMS[0]            = diagDir.make<TH1D>("HitRMS_Ind","Hit RMS (induction plane)",120,0.,30.);
  hHitRMS[1]            = diagDir.make<TH1D>("HitRMS_Col","Hit RMS (collection plane)",120,0.,30.);
  hHitAmplitude[0]      = diagDir.make<TH1D>("HitAmplitude_Ind","Hit amplitude (induction plane);ADC",150,0,300);
  hHitAmplitude[1]      = diagDir.make<TH1D>("HitAmplitude_Col","Hit amplitude (collection plane);ADC",150,0,300);
  hHitSummedADC            = diagDir.make<TH1D>("HitSummedADC","Hit summed ADC (collection plane);ADC",100,0, 50000);
  hHitIntegral            = diagDir.make<TH1D>("HitIntegral","Hit integral (collection plane);ADC",100,0, 50000);
  hHitCharge            = diagDir.make<TH1D>("HitCharge","Hit charge after lifetime correction (collection plane);Q (e-)",100,0, 1000000);
  hTotalCharge          = diagDir.make<TH1D>("TotalCharge","Total charge on collection plane;Charge Q [e-]",100,0,20000000);
  
  // ===================================================================
  // Reco track diagnostics
  hNumTrack             = diagDir.make<TH1D>("NumTrack","Number of tracks per event",30,0,30);
  hTrackNode_X          = diagDir.make<TH1D>("TrackNode_X","Distribution of reco track nodes;X [cm]", 600, -5., 55.);
  hTrackNode_Y          = diagDir.make<TH1D>("TrackNode_Y","Distribution of reco track nodes;Y [cm]", 500, -25., 25.);
  hTrackNode_Z          = diagDir.make<TH1D>("TrackNode_Z","Distribution of reco track nodes;Z [cm]", 1000, -5., 95.);
  hTrackNode_ZX         = diagDir.make<TH2D>("TrackNode_ZX","Reco track nodes (top-down view);Z [cm];X [cm]",200, -5., 95., 120, -5., 55.);
  hTrackNode_ZY         = diagDir.make<TH2D>("TrackNode_ZY","Reco track nodes (side view);Z [cm];Y [cm]",200, -5., 95., 100, -25., 25.);
  hTrackNode_ZX         ->SetOption("colz");
  hTrackNode_ZY         ->SetOption("colz");
  hCrsTrackNode_ZX         = diagDir.make<TH2D>("CrsTrackNode_ZX","Reco track nodes (top-down view);Z [cm];X [cm]",200, -5., 95., 120, -5., 55.);
  hCrsTrackNode_ZY         = diagDir.make<TH2D>("CrsTrackNode_ZY","Reco track nodes (side view);Z [cm];Y [cm]",200, -5., 95., 100, -25., 25.);
  hCrsTrackNode_ZX         ->SetOption("colz");
  hCrsTrackNode_ZY         ->SetOption("colz");
  hTrackLength          = diagDir.make<TH1D>("TrackLength","Track lengths;L [cm]",240,0.,120.);
  hNumTrackStopping     = diagDir.make<TH1D>("NumTrackStopping","Number of identified stopping tracks per event",20,0,20);
  hNumTrackCrossing     = diagDir.make<TH1D>("NumTrackCrossing","Number of identified passing tracks per event",20,0,20);
  hNumTrackContained    = diagDir.make<TH1D>("NumTrackContained","Number of identified contained tracks per event",20,0,20);
  hMudEdx               = diagDir.make<TH1D>("MudEdx","Reco dE/dx of tagged stopping muon track;dE/dx [MeV/cm]",100,0,10);
  hMuResRangeVsdEdx     = diagDir.make<TH2D>("MuResRangeVsdEdx","Residual range vs. dE/dx of tagged stopping muon track;Residual range [cm];dE/dx [MeV/cm]",100,0.,40.,100,0.,10.);
  hMuResRangeVsdEdx     ->SetOption("colz");
  
  // ======================================================================
  // 2D clustering and shower reco performance diagnostics  
  hClusterHitSeparation = diagDir.make<TH1D>("ClusterHitSeparation","Nearest neighbor hit distance in proximity clustering",50,0.,5.);
  hClusterSize          = diagDir.make<TH1D>("ClusterSize","Cluster size (num of hits)",300,0,300);
  hElShowerFrac        = diagDir.make<TH1D>("ElShowerFrac","Fraction of total non-muon hits included in electron shower",100,0,1.);
  hMuClusterSize        = diagDir.make<TH1D>("MuClusterSize","Number of #mu hits in cluster",200,0,200);
  hElClusterSize        = diagDir.make<TH1D>("ElClusterSize","Number of electron hits in cluster",200,0,200);
  hElShowerSize        = diagDir.make<TH1D>("ElShowerSize","Number of hits in shower (collection plane)",200,0,200);
  hExtraHits            = diagDir.make<TH1D>("ExtraHits","Extra (unclustered) hits near 2D #mu endpoint",20,0,20);
  hDistMaxTdQdsToMaxQ   = diagDir.make<TH1D>("DistMaxTdQdsToMaxQ","Distance btwn max trunc. mean Q/ds and max Q/ds;Hit_{Max dQ/ds} - Hit_{Max trunc. dQ/ds} [hits]",30,-15,15);
  hFracMuHitsLinear     = diagDir.make<TH1D>("FracMuHitsLinear","Fraction of #mu hits > covariance thresh",100,0,1.);
  hMuAveLinearity     = diagDir.make<TH1D>("MuAveLinearity","Average muon linearity",100,0,1.);
  hMuClusterHitsEndFit  = diagDir.make<TH1D>("MuClusterHitsEndFit","Number of #mu hits used for direction calcuation",30,0,30);
  hDecayAngle2D         = diagDir.make<TH1D>("DecayAngle2D","Michel decay angle rel. to terminal #mu direction;Angle [deg]",90,0.,180.);
  hElShowerDepAngle2D   = diagDir.make<TH1D>("ElShowerDepAngle2D","Opening angle of Michel shower depositions;Angle [deg]",90,0.,180.);
  hElShowerAngleMean    = diagDir.make<TH1D>("ElShowerAngleMean","Mean theta (rel to mu dir) of shower deposits;#theta [deg]",60,0.,60.);
  hElShowerAngleRMS    = diagDir.make<TH1D>("ElShowerAngleRMS","Angular RMS (rel to mu dir) of shower deposits;#sigma_{#theta} [deg]",60,0.,60.);
  hShowerSizeCompare   = diagDir.make<TH2D>("ShowerSizeCompare","Electron shower size on both planes;Collection Plane;Induction Plane",100,0,100,100,0,100);
  hShowerSizeCompare    ->SetOption("colz");
  hElClusterSizeCompare   = diagDir.make<TH2D>("ElClusterSizeCompare","Electron cluster size on both planes;Collection Plane;Induction Plane",100,0,100,100,0,100);
  hElClusterSizeCompare    ->SetOption("colz");
  hElOffsetT            = diagDir.make<TH1D>("ElOffsetT","Time difference between first electron / last #mu hit (same wire only);#DeltaT [#mus]",200,-20.,20);
  hHitTimeDiff          = diagDir.make<TH1D>("HitTimeDiff","Hit time differences between Michel shower hits in collection/induction planes;T(coll.) - T(ind.) [#mus]",200,-10.,10.);
  
  // ======================================================================
  // 3D shower reco performance diagnostics  
  hMuEndHitTimeDiff     = diagDir.make<TH1D>("MuEndHitTimeDiff","Hit time difference btwn candidate #mu end hits in collection/induction plane clusters;#DeltaT [#mus]",200,-5.,5.);
  hMuEnd3D_ZX         = diagDir.make<TH2D>("MuEnd3D_ZX","3D #mu endpoint (top-down view);Z [cm];X [cm]",500, -5., 95., 300, -5., 55.);
  hMuEnd3D_ZY         = diagDir.make<TH2D>("MuEnd3D_ZY","3D #mu endpoint (side view);Z [cm];Y [cm]",500, -5., 95., 250, -25., 25.);
  hMuEnd3D_ZX         ->SetOption("colz");
  hMuEnd3D_ZY         ->SetOption("colz");
  hFracHits3D            = diagDir.make<TH1D>("FracHits3D","Fraction of Michel shower hits made into 3D pts",120,0,1.2);
  hNumPts3D             = diagDir.make<TH1D>("NumPts3D","Number of 3D Michel shower hits",100,0,100);

  // ======================================================================
  // Michel electron charge and energy
  hElShowerCharge         = tfs->make<TH1D>("ElShowerCharge","Reco charge of Michel electron;Reconstructed #it{Q} [e-]",100,0,5000000);
  hElShowerPhotons        = tfs->make<TH1D>("ElShowerPhotons","Photons produced by Michel electron;Reconstructed #it{L} [#gamma]",100,0,5000000);
  hElShowerEnergy         = tfs->make<TH1D>("ElShowerEnergy","Reco Michel Electron Shower Energy;Reconstructed energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hElShowerEnergyQL       = tfs->make<TH1D>("ElShowerEnergyQL","Reco Michel Electron Shower Energy (Q+L);Reconstructed energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hElTrackdEdx            = tfs->make<TH1D>("ElTrackdEdx","Averaged dE/dx along 3D electron track;dE/dx [MeV/cm]",100,0.,10.);
  hQoverL                 = tfs->make<TH1D>("QoverL","Q/L for Michel electrons;Q/L",80,0.,4.);
  hQoverLCrsMuTrk         = tfs->make<TH1D>("QoverL_CrsMuTrk","Q/L for crossing muon tracks;Q/L",80,0.,4.);
  hRecomb                 = tfs->make<TH1D>("Recomb","Recombination (e- survival fraction) calculated from Q/L;Electron survival fraction #it{R}",140,-0.2,1.2);
  hRecombCrsMuTrk         = tfs->make<TH1D>("Recomb_CrsMuTrk","Recombination calculated from Q/L (crossing muons);Electron survival fraction #{R}",140,-0.2,1.2);
  hCrsMuQPerCm            = tfs->make<TH1D>("CrsMuQPerCm","Q/cm for crossing muons;Charge deposited per unit length [e-/cm]",140,0.,14e4);
  hCrsMuLPerCm            = tfs->make<TH1D>("CrsMuLPerCm","L/cm for crossing muons;Photons produced per unit length [#gamma/cm]",140,0.,14e4);
    
  // ======================================================================
  // Truth plots
  hTrue_X                 = truthDir.make<TH1D>("True_X","X position of energy depositions [cm]",600,-5.,55.);
  hTrue_Y                 = truthDir.make<TH1D>("True_Y","Y position of energy depositions [cm]",500,-25.,25.);
  hTrue_Z                 = truthDir.make<TH1D>("True_Z","Z position of energy depositions [cm]",1000,-5.,95.);
  hTrue_TrajLength        = truthDir.make<TH1D>("True_TrajLength","Particle trajectory length;L [cm]",500,0.,5.);
  hTrue_dQdx              = truthDir.make<TH1D>("True_dQdx","dQ/dx (avg. per particle) for Michel trk + shower depositions;dQ/dx [e-/cm]",80,0.,4e5);
  hTrue_dQdx_ElTrk        = truthDir.make<TH1D>("True_dQdx_ElTrk","dQ/dx (avg. per particle) for Michel trk depositions;dQ/dx [e-/cm]",80,0.,4e5);
  hTrue_dQdx_ElShw        = truthDir.make<TH1D>("True_dQdx_ElShw","dQ/dx (avg. per particle) for Michel shower depositions;dQ/dx [e-/cm]",80,0.,4e5);
  hTrue_dEdx              = truthDir.make<TH1D>("True_dEdx","dE/dx (avg. per particle) for Michel trk + shower depositions;dE/dx [MeV/cm]",60,0.,12.);
  hTrue_dEdx_ElTrk        = truthDir.make<TH1D>("True_dEdx_ElTrk","dE/dx (avg. per particle) for Michel trk depositions;dE/dx [MeV/cm]",60,0.,12.);
  hTrue_dEdx_ElShw        = truthDir.make<TH1D>("True_dEdx_ElShw","dE/dx (avg. per particle) for Michel shower depositions;dE/dx [MeV/cm]",60,0.,12.);
  hTrue_dEdx_vs_dQdx      = truthDir.make<TH2D>("True_dEdx_vs_dQdx","dE/dx vs. dQ/dx (avg. per particle) for Michel track+shower depositions;dE/dx [MeV/cm];dQ/dx [e-/cm]",
    80,0.,8.,
    125,0.,250e3);
  hTrue_dEdx_vs_dQdx_ElTrk      = truthDir.make<TH2D>("True_dEdx_vs_dQdx_ElTrk","dE/dx vs. dQ/dx (avg. per particle) for Michel track depositions;dE/dx [MeV/cm];dQ/dx [e-/cm]",
    80,0.,8.,
    125,0.,250e3);
  hTrue_dEdx_vs_dQdx_ElShw      = truthDir.make<TH2D>("True_dEdx_vs_dQdx_ElShw","dE/dx vs. dQ/dx (avg. per particle) for Michel shower depositions;dE/dx [MeV/cm];dQ/dx [e-/cm]",
    80,0.,8.,
    125,0.,250e3);
  hTrue_dEdx              ->SetOption("hist");
  hTrue_dQdx              ->SetOption("hist");
  hTrue_dEdx_vs_dQdx      ->SetOption("colz");
  hTrue_dEdx_vs_dQdx_ElTrk->SetOption("colz");
  hTrue_dEdx_vs_dQdx_ElShw->SetOption("colz");
  hTrue_RecombFactor      = truthDir.make<TH1D>("True_RecombFactor","Survival prob of ionization electrons for all Michel depositions;P",100,0.,1.);
  hTrue_RecombFactor_ElTrk      = truthDir.make<TH1D>("True_RecombFactor_ElTrk","Survival prob of ionization electrons for Michel tracks;P",100,0.,1.);
  hTrue_RecombFactor_ElShw      = truthDir.make<TH1D>("True_RecombFactor_ElShw","Survival prob of ionization electrons for Michel shower hoton depositions;P",100,0.,1.);
  hTrue_NumDRays        = truthDir.make<TH1D>("True_NumDRays","Number delta rays per event",50,0,50);
  hTrue_DRayEnergy      = truthDir.make<TH1D>("True_DRayEnergy","Energy deposited by delta rays off muon;Energy deposited [MeV]",30,0.,3.);
  hTrue_NumElectrons    = truthDir.make<TH1D>("True_NumElectrons","Number of electrons from IDE objects",400,0,4000);
  hTrue_ElEnergy        = truthDir.make<TH1D>("True_ElEnergy","Electron energy",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hTrue_ElEnergyFree    = truthDir.make<TH1D>("True_ElEnergyFree","True Michel Electron Energy (from #mu+)",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hTrue_ElEnergyCapture = truthDir.make<TH1D>("True_ElEnergyCapture","True Michel Electron Energy (from #mu-)",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hTrue_ElTrackEnergyDep     = truthDir.make<TH1D>("True_ElTrackEnergyDep","Michel track energy deposited;Energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hTrue_ElShowerEnergyDep     = truthDir.make<TH1D>("True_ElShowerEnergyDep","True Deposited Michel Electron Energy;Energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hTrue_ElTrackChargeDep     = truthDir.make<TH1D>("True_ElTrackChargeDep","True Deposited Michel Electron track Charge;Electrons [e-]",300,0,3000000);
  hTrue_ElShowerChargeDep     = truthDir.make<TH1D>("True_ElShowerChargeDep","True Deposited Michel Electron Charge;Electrons [e-]",300,0,3000000);
  hTrue_EnergyDepVsRecoEnergy = truthDir.make<TH2D>("True_EnergyDepVsRecoEnergy","True Michel energy deposited in LAr vs. reco energy;True deposited Michel electron energy [MeV];Reconstructed energy [MeV]",
    120,0,60,240,0,120);
  hTrue_EnergyDepVsRecoEnergy->SetOption("colz");
  hTrue_EnergyRes     = truthDir.make<TH1D>("True_EnergyRes","Michel electron deposited energy resolution;(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyResTrk     = truthDir.make<TH1D>("True_EnergyResTrk","Michel electron deposited energy resolution (electron track only);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyResQL     = truthDir.make<TH1D>("True_EnergyResQL","Michel electron deposited energy resolution (Q+L);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyRes_Shwr3D    = truthDir.make<TH1D>("True_EnergyRes_Shwr3d","Michel electron deposited energy resolution (3D shower locations found);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyVsEnergyRes     = truthDir.make<TH2D>("True_EnergyVsEnergyRes","Michel electron deposited energy resolution;True Deposited Energy [MeV];(reco-true)/true",60,0.,60.,150,-1.5,1.5);
  hTrue_EnergyVsEnergyRes     ->SetOption("colz");
  hTrue_ElShowerPhotons = truthDir.make<TH1D>("True_ElShowerPhotons","Number of scintillation photons produced",300,0,3000000);
  hTrue_MuLateLightContamination = truthDir.make<TH2D>("True_MuLateLightContamination","ETL PMT;#DeltaT [ns];Contamination Fraction", 75,0.,7500, 110, 0.,1.1);
  hTrue_MuLateLightContamination->SetOption("colz");
  hTrue_VisRes          = truthDir.make<TH1D>("True_VisRes","Resolution of reconstructed Michel shower vis.;(reco-true)/true",150,-1.5,1.5);
  hTrue_PhotonRes       = truthDir.make<TH1D>("True_PhotonRes","Photon resolution of Michel candidate;(reco-true)/true",150,-1.5,1.5);
  hTrue_MuTrackEnd_ZX   = truthDir.make<TH2D>("True_MuTrackEnd_ZX","True #mu endpoint;Z [cm];X [cm]",100,-5.,95, 57,-5.,52.);
  hTrue_MuTrackEnd_ZY   = truthDir.make<TH2D>("True_MuTrackEnd_ZY","True #mu endpoint;Z [cm];Y [cm]",100,-5.,95, 50,-25.,25.);
  hTrue_MuTrackEnd_ZX   ->SetOption("colz");
  hTrue_MuTrackEnd_ZY   ->SetOption("colz");
  hTrue_MuTrackVertexRes_X = truthDir.make<TH1D>("True_MuTrackVertexRes_X","Tagged #mu track startpoint resolution (X);Xreco - Xtrue [cm]",200,-2.,2);
  hTrue_MuTrackVertexRes_Y = truthDir.make<TH1D>("True_MuTrackVertexRes_Y","Tagged #mu track startpoint resolution (Y);Yreco - Ytrue [cm]",200,-2.,2);
  hTrue_MuTrackVertexRes_Z = truthDir.make<TH1D>("True_MuTrackVertexRes_Z","Tagged #mu track startpoint resolution (Z);Zreco - Ztrue [cm]",200,-2.,2);
  hTrue_MuTrackEndRes_X = truthDir.make<TH1D>("True_MuTrackEndRes_X","Tagged #mu track endpoint resolution (X);Xreco - Xtrue [cm]",200,-2.,2.);
  hTrue_MuTrackEndRes_Y = truthDir.make<TH1D>("True_MuTrackEndRes_Y","Tagged #mu track endpoint resolution (Y);Yreco - Ytrue [cm]",200,-2.,2);
  hTrue_MuTrackEndRes_Z = truthDir.make<TH1D>("True_MuTrackEndRes_Z","Tagged #mu track endpoint resolution (Z);Zreco - Ztrue [cm]",200,-2.,2);
  hTrue_ElShowerDepVsDistFromMuTrackEnd = truthDir.make<TH1D>("True_ElShowerDepVsDistFromMuTrackEnd","True distance of Michel shower depositions;Distance [cm]",100, 0., 100.);
  hTrue_BremmPhotonLength = truthDir.make<TH1D>("True_BremmPhotonLength","Bremm photon length;cm",100,0.,100);
  hTrue_ElEnergyVsNumBremmPhotons = truthDir.make<TH2D>("True_ElEnergyVsNumBremmPhotons",";True Michel electron energy [MeV];Number of Bremmstralung photons",60,0.,60., 30,0.,30.);
  hTrue_ElEnergyVsNumBremmPhotons ->SetOption("colz");

}





//########################################################################################
// Resets all member data variables/counters before the start of a new event
void MichelAnaFilter::ResetVariables()
{
  fMichelOpticalID          = false;
  fNumPlaneHits[0]          = 0;
  fNumPlaneHits[1]          = 0;
  fPulses.clear();
  fRunNumber                = -99;
  fSubRunNumber             = -99;
  fEventNumber              = -99;
  fEventTime                = -99;
  fDecayTime                = -999.;
  fPE                       = -999.;
  fPEPrompt                       = -999.;
  for(size_t i=0; i<2; i++) {
    fWfmRMS[i]              = -99.; 
    fNumOpHits0[i]          = -9; 
    fNumOpHits[i]           = -9; 
    fdT[i]                  = -999.;
    fWidth[i]               = -999.;
    fMuWidth[i]             = -999.;
    fAmplitude[i]           = -999.;
    fPE_total[i]            = -999.;
    fPE_totalRaw[i]         = -999.;
    fPE_prompt[i]           = -999.;
    fPE_promptRaw[i]        = -999.;
    fPromptFracRaw[i]       = -999.;
    fPromptFrac[i]          = -999.;
    fMuAmplitude[i]         = -999.;
    fMuPE_prompt[i]         = -999.;
    fMuPulseSaturated[i]    = false;
    fElPulseSaturated[i]    = false;
    fvHitTimes0[i]           .clear();   
    fvHitTimes[i]           .clear();   
    fvHitAmplitudes[i]      .clear();
    fvHitADC_100ns[i]      .clear();
    fvHitADC_200ns[i]      .clear();
    fvHitADC_300ns[i]      .clear();
    fvHitADC_400ns[i]      .clear();
    fvHitADC_500ns[i]      .clear();
    fvHitADC_600ns[i]      .clear();
    fvHitADC_700ns[i]      .clear();
    fvHitADC_900ns[i]      .clear();
    fvHitADC_1200ns[i]      .clear();
    fvHitADC_1500ns[i]      .clear();
    fvHitADC_1800ns[i]      .clear();
    fvHitADC_7000ns[i]      .clear();
    fvHitADC_7000ns[i]      .clear();
    fvHitADCpf_100ns[i]     .clear();
    fvHitADCpf_7000ns[i]    .clear();
    fvHitWidth[i]            .clear();
    fvPrepulseBaseline[i]   .clear();
    fvPrepulseZeroPoint[i]  .clear();
    fvPrepulseRMS[i]        .clear();
    fvPrepulseSlowNorm[i]   .clear();
    fvPrepulseSlowTau[i]    .clear();
    fvPrepulseX1[i]         .clear();
    fvIsHitAtTrigger[i]     .clear();
    fvIsHitSaturated[i]     .clear();
    fTrue_PE_prompt[i]       = -999.;
    fTrue_Amplitude[i]       = -999.;
    fTrue_MuPE_prompt[i]       = -999.;
    fTrue_MuContam_prompt[i] = -999.;
    fTrue_MuContam_total[i] = -999.;
    fTrue_PE_total[i]       = -999.;
    fTrue_PE[i]             = -999.;
    hPMT_phelTimes_electron[i]->Reset();
    hPMT_phelTimes[i]       ->Reset();
    fPMT_wfm[i]             .clear();
    fvbs[i]                 .clear();
    fT0[i]                  = 0;
    fMuT0[i]                = 0;
  }
  fNumTrack                = -9;
  fNumTrackStopping         = 0;
  fNumTrackCrossing          = 0;
  fNumTrackContained        = 0;
  fvMuTrackdEdx            .clear();
  fvMuTrackResidualRange   .clear();
  fMuTrackIsBraggPeaked     = false;
  fMuTrackIsCaloOrdered     = false;
  fMuTrackHits              .clear();
  fMuTrackID                = -999 ;
  fMuTrackIndex             = -999 ;
  fMuTrackLength            = -999.;
  fMuTrackEnergy            = -999.;
  fMuTrackZenithAngle            = -999.;
  fMuTrackVertex            .SetXYZ(-99.,-99.,-99.);
  fMuTrackVertex_X           = -99.;
  fMuTrackVertex_Y           = -99.;
  fMuTrackVertex_Z           = -99.;
  fMuTrackEnd_X           = -99.;
  fMuTrackEnd_Y           = -99.;
  fMuTrackEnd_Z           = -99.;
  fMuTrackEnd               .SetXYZ(-99.,-99.,-99.);

  // Mu calibration variables
  fCrossingMuVertex_X       = -999.;
  fCrossingMuVertex_Y       = -999.;
  fCrossingMuVertex_Z       = -999.;
  fCrossingMuEnd_X       = -999.;
  fCrossingMuEnd_Y       = -999.;
  fCrossingMuEnd_Z       = -999.;
  fCrossingMuTrackIndex     = -999 ;
  fCrossingMuLength        = -9.;
  fCrossingMuCharge        = -9999.;
  fCrossingMuPhotons         = -9999.;
  fCrossingMuPhotonsPrompt         = -9999.;

  // Shower/charge
  fHitPlane               .clear();
  fHitWire                .clear();
  fHitT                   .clear();
  fHitX                   .clear();
  fHitW                   .clear();
  fHitCharge              .clear();
  
  fElOffsetT            = -99.;
  fAveDriftTime         = -9.;
  fAveX                 = -9.;
  fMuEndWire            = -9;
  fElStartWire          = -9;
  
  fClusterSize          = -9;
  fMuClusterSize        = -9;
  fElClusterSize        = -9;
  fExtraHits            = -9;
  fMuAveLinearity       = -9.;
  fFracMuHitsLinear     = -9.;
  fMuClusterHitsEndFit  = -9;
  fDecayAngle2D         = -999;
  fElShowerAngleMean    = -9.;
  fElShowerAngleRMS     = -9.;
  fElShowerSize         = -9;
  fElShowerFrac         = -9.; 
  fTotalCharge          = -9.;
  fElShowerCharge       = -9.;
  fElTrackCharge        = -9.;
  fElShowerEnergy       = -9.;
  fElTrackEnergy        = -9.;                 
  
  fClusterSize_Pl0          = -9;
  fMuClusterSize_Pl0        = -9;
  fElClusterSize_Pl0        = -9;
  fExtraHits_Pl0            = -9;
  fMuAveLinearity_Pl0       = -9.;
  fFracMuHitsLinear_Pl0     = -9.;
  fMuClusterHitsEndFit_Pl0  = -9;
  fDecayAngle2D_Pl0         = -999;
  fElShowerAngleMean_Pl0    = -9.;
  fElShowerAngleRMS_Pl0     = -9.;
  fElShowerSize_Pl0         = -9;
  fElShowerFrac_Pl0         = -9.; 
  fTotalCharge_Pl0          = -9.;
  fElShowerCharge_Pl0       = -9.;
  fElTrackCharge_Pl0        = -9.;
  fElShowerEnergy_Pl0       = -9.;
  fElTrackEnergy_Pl0        = -9.;                 
  
  //fMuEnd3D                  .SetXYZ(-99., -99., -99.);
  fMuEndHitTimeDiff         = -999.;
  fElShowerPhotons          = -999.;
  fElTrackdEdx              = -999.;
  fFracHits3D               = -9.;
  fNumPts3D                 = -9;
  fNumPts3DTrk              = -9;
  fElShowerVis              = -0.09;
  fElShowerVisCh[0]         = -0.09;
  fElShowerVisCh[1]         = -0.09;
  //fElShowerCentroid         .SetXYZ(-99,-99,-99.);
  fElShowerCentroid_X       = -99.;
  fElShowerCentroid_Y       = -99.;
  fElShowerCentroid_Z       = -99.;
  fElShowerEnergyQL         = -999.;
  
  // Truth information
  fTrue_MuCharge          = 0;
  fTrue_MuStopTime        = -999.;
  fTrue_MuTrackEnd          .SetXYZ(-99.,-99.,-99.);
  fTrue_MuTrackEnd_X      = -99.;
  fTrue_MuTrackEnd_Y      = -99.;
  fTrue_MuTrackEnd_Z      = -99.;
  fTrue_MuTrackVertex_X      = -99.;
  fTrue_MuTrackVertex_Y      = -99.;
  fTrue_MuTrackVertex_Z      = -99.;
  fTrue_MuTrackVertex          .SetXYZ(-99.,-99.,-99.);
  fTrue_MuEnergyDep       = -999.;
  fTrue_ElMomentum        .SetXYZ(0.,0.,0.);
  fTrue_ElMomentum_X      = -9.;
  fTrue_ElMomentum_Y      = -9.;
  fTrue_ElMomentum_Z      = -9.;
  fTrue_ElAngle           = -999.;
  fTrue_ElEnergy          = -999.;
  fTrue_ElShowerVis       = -0.09;
  fTrue_ElShowerVisCh[0]       = -0.09;
  fTrue_ElShowerVisCh[1]       = -0.09;
  fTrue_ElShowerPhotons   = -999;
  fTrue_ElTrackPhotons   = -999;
  fTrue_ElShowerPhotonsPrompt   = -999;
  fTrue_ElShowerPhotonsLate   = -999;
  fTrue_ElTrackEnergyDep  = -999.;
  fTrue_ElTrackChargeDep  = -999.;
  fTrue_ElShowerChargeDep = -999.;
  fTrue_ElShowerEnergyDep = -999.;
  fTrue_IsElContained     = false;
  fTrue_TotalEnergyDep    = -999.;
  fTrue_TotalChargeDep    = -999.;
  fTrue_dT                = -999;
}



//########################################################################################
//
// The main MichelAnaFilter algorithm. This function takes in an event that has already gone
// through standard reconstruction, then carries out several major stages:
//  (1) Check reco tracks and identify events with stopping muons (for Michel clustering) 
//      or crossing muons (for calibration).
//  (2) Optical reconstruction of PMT waveforms, identification of Michel decay-like topology. 
//  (3) Charge clustering of reco hits to create Michel showers in 2D "wire-time(x)" space.
//      This procedure is adapted from MicroBooNE's Michel electron reconstruction, but with
//      several modifications and additions. Charge deposited by the Michel electron shower (Q) 
//      is computed as well as the total shower energy.
//  (4) Conversion of 2D electron clusters (on both planes) into a collection of 3D spacepoints
//      by matching hits in time between the planes. 
//  (5) Using MC scintillation visibility libraries, average fractional visibility of the 
//      3D Michel shower is calculated and used to scale up detected photoelectrons to determine
//      the number of VUV photons produced (L).  Michel shower energy using Q+L is computed.
//
//  This is an art filter function, so it can be used to pass/reject events to filter out 
//  a sample of Michel electron events.  This functionality is controlled by parameters with 
//  the "fFilter_" prefix.
//
//########################################################################################
bool MichelAnaFilter::filter(art::Event & e)
{
  // ===================================================================
  // Reset all the member data
  ResetVariables();
  fNumEvents++;
  
  // Get run, subrun and event number
  fRunNumber    = (int)e.run();
  fSubRunNumber = (int)e.subRun();
  fEventNumber  = (int)e.event();
  fEventTime    = (int)e.getSubRun().beginTime().value();
  fIsMC         = !(bool)e.isRealData(); 
  if( fRunNumber <= 0 ) return false;
   
  // Get detector properties
  GetDetProperties(); 
  
  LOG_VERBATIM("MichelAnaFilter")
  <<"\n"
  <<"Beginning filter: run "<<fRunNumber<<", subrun "<<fSubRunNumber<<", event "<<fEventNumber<<" (MC = "<<fIsMC<<")\n"
  <<"total evts processed: "<<fNumEvents
  <<", OpID: "<<fNumOpticalID
  <<", 2D showers: "<<fNumEvents_Shwr
  <<", 3D showers: "<<fNumEvents_Shwr3D
  <<", calib: "<<fNumEvents_Calibration; 
  LOG_VERBATIM("MichelAnaFilter")
  <<"Efield: "<<fEfield<<"  Drift vel.: "<<fDriftVelocity<<"  Electron lifetime: "<<fElectronLifetime
  <<"  SPE: "<<fSPE[0]<<","<<fSPE[1]<<"; eff. tau: "<<fMuContamCorr_EffTau[0]<<","<<fMuContamCorr_EffTau[1]<<"  (true TauT: "<<fTauT<<")";
  
  // Skip event if the electron lifetime isn't defined for this run
  hElectronLifetime ->Fill(fElectronLifetime); 
  hEffTau[0]        ->Fill(fMuContamCorr_EffTau[0]);
  hEffTau[1]        ->Fill(fMuContamCorr_EffTau[1]);
  if( fElectronLifetime <= 0 ) return false;
 
  
  // ==================================================================  
  // If MC, then get the Truth information. This function loops through
  // the MCParticle list and calculates truth-level quantities.
  if( fIsMC ) GetTruthInfo( e );


  // =================================================================
  // Trigger filter (not sure if this works... excluded for now)
  /*
  if( !fIsMC ) {
    // Implementation of required member function here.
    //Get the triggers from the event record
    art::Handle< std::vector<raw::Trigger> > triggerHandle;
    e.getByLabel("daq",triggerHandle);
    // Get Trigger info
    if (triggerHandle.isValid() && triggerHandle->size() > 0){
      uint32_t triggerBits = triggerHandle->at(0).TriggerBits();
      std::bitset<32> triggerBitSet(triggerBits);
      std::cout<<"TRIGGER BIT 13 = "<<triggerBitSet[13]<<"\n";
    }
  }
  */

 
  // ================================================================
  // Set booleans for various conditions used throughout the analysis
  bool promptPEcut          = true;
  bool isOneStoppingTrack   = false;
  bool isOneCrossingTrack   = false;
  bool isCalibrationEvent   = false;
  bool isClusterBndFound    = false;
  

  
  // *******************************************************************************
  // Save track info and look for single stopping tracks or crossing tracks
  // *******************************************************************************
  if( fLookAtTracks ) {
     
    // ================================================================
    // Get the tracks and their associated energy
    art::Handle< std::vector< recob::Track >> TrackHandle;
    e.getByLabel(fTrackModule,TrackHandle);
    // Association between Calorimetry objects and Tracks
    art::FindManyP<anab::Calorimetry> fmcal(TrackHandle, e, fTrackCalModule);
    fNumTrack = (int)TrackHandle->size();
    LOG_VERBATIM("MichelAnaFilter") << "Using track module "<<fTrackModule<<", number of tracks: "<<fNumTrack;
  
    // initialize the index of the identified stopping or crossing track 
    // as well as other variables (to prevent repeated loops thru tracks)
    int   crsTrkIndex     = -9;
    float crsTrkLength    = -9.; 
    int   stpTrkIndex     = -9;
    int   stpTrkID        = -9;
    float stpTrkLength    = -9.; 
    TVector3 stpTrkVertex;
    TVector3 stpTrkEnd;
    TVector3 crsTrkVertex;
    TVector3 crsTrkEnd; 
  
    
    // ================================================================
    // Loop over the tracks in the event
    for( int track_index = 0; track_index < std::min(fNumTrack,(int)kMaxTracks); track_index++){
      
      // ------------------------------------------------------------
      // Get the recob::Track object and record its endpoint/vertex, 
      // then save/plot a bunch of info about it.
      art::Ptr< recob::Track > TheTrackPtr(TrackHandle,track_index);
      recob::Track TheTrack = *TheTrackPtr;
      hTrackNode_X  ->Fill( TheTrack.Vertex().X() );
      hTrackNode_X  ->Fill( TheTrack.End().X() );
      hTrackNode_Y  ->Fill( TheTrack.Vertex().Y() );
      hTrackNode_Y  ->Fill( TheTrack.End().Y() );
      hTrackNode_Z  ->Fill( TheTrack.Vertex().Z() );
      hTrackNode_Z  ->Fill( TheTrack.End().Z() );
      hTrackNode_ZX ->Fill( TheTrack.End().Z(), TheTrack.End().X() );
      hTrackNode_ZX ->Fill( TheTrack.Vertex().Z(), TheTrack.Vertex().X() );
      hTrackNode_ZY ->Fill( TheTrack.End().Z(), TheTrack.End().Y() );
      hTrackNode_ZY ->Fill( TheTrack.Vertex().Z(), TheTrack.Vertex().Y() );
      hTrackLength  ->Fill( TheTrack.Length() );
      
      // ------------------------------------------------------------
      // Is track stopping, passing, or contained? 
      bool endIsInFid     = IsPointInFiducialVolume(TheTrack.End(),fFiducialMarginX,fFiducialMarginY,fFiducialMarginZ);
      bool vertexIsInFid  = IsPointInFiducialVolume(TheTrack.Vertex(),fFiducialMarginX,fFiducialMarginY,fFiducialMarginZ);
      bool isCrossingTrack  = (!endIsInFid && !vertexIsInFid);
      bool isContainedTrack = ( endIsInFid && vertexIsInFid);
      bool isStoppingTrack  = ( ((!endIsInFid && vertexIsInFid)||(endIsInFid && !vertexIsInFid)) 
                              &&(TheTrack.Length()>fMuTrackLengthMin));
      fNumTrackStopping   += isStoppingTrack; 
      fNumTrackCrossing   += isCrossingTrack;
      fNumTrackContained  += isContainedTrack;
      
      // ------------------------------------------------------------
      // Call the node that's outside the fiducial volume the "vertex".
      // If both are outside fiducial volume (ie, crossing track), then
      // whichever has higher 'Y' coordinate is the vertex.
      TVector3 vertex = TheTrack.Vertex();
      TVector3 end    = TheTrack.End();
      if( (isStoppingTrack && vertexIsInFid) || (isCrossingTrack && vertex.Y() < end.Y() ) ){
        vertex  = TheTrack.End();
        end     = TheTrack.Vertex();
      } 
      
      if( isStoppingTrack ) {
        stpTrkIndex = track_index;
        stpTrkID    = TheTrack.ID();
        stpTrkLength= TheTrack.Length();
        stpTrkVertex= vertex;
        stpTrkEnd   = end;
      }
      if( isCrossingTrack ) {
        crsTrkIndex   = track_index;
        crsTrkLength  = TheTrack.Length();
        crsTrkVertex  = vertex;
        crsTrkEnd     = end;
      }
      
      LOG_VERBATIM("MicheAnaFilter")
          <<"  track "<<track_index<<", N = "<<TheTrack.NumberTrajectoryPoints()<<"  length "<<TheTrack.Length()<<"  vertex("
          << TheTrack.Vertex().X() <<"," 
          << TheTrack.Vertex().Y() << "," 
          << TheTrack.Vertex().Z() << ")->InFiducial()="
          << vertexIsInFid 
          << "   end("
          << TheTrack.End().X() <<"," 
          << TheTrack.End().Y() << "," 
          << TheTrack.End().Z() << ")->InFiducial()="
          << endIsInFid;
    
    }//<-- end first pass loop through tracks
    LOG_VERBATIM("MichelAnaFilter") << "Done looping tracks: num stop/pass/contained "<<fNumTrackStopping<<" / "<<fNumTrackCrossing<<" / "<<fNumTrackContained;
    hNumTrack         ->Fill(fNumTrack);   
    hNumTrackStopping ->Fill(fNumTrackStopping);
    hNumTrackContained ->Fill(fNumTrackContained);
    hNumTrackCrossing ->Fill(fNumTrackCrossing);
  

    // ===================================================================
    // Identify events with 1 crossing track and no others for calibration
    if( fNumTrackCrossing == 1 && fNumTrack == 1 ) {
        isOneCrossingTrack    = true; 
        fCrossingMuTrackIndex = crsTrkIndex;
        fCrossingMuLength     = crsTrkLength; 
        fCrossingMuVertex_X   = crsTrkVertex.X();
        fCrossingMuVertex_Y   = crsTrkVertex.Y();
        fCrossingMuVertex_Z   = crsTrkVertex.Z();
        fCrossingMuEnd_X      = crsTrkEnd.X();
        fCrossingMuEnd_Y      = crsTrkEnd.Y();
        fCrossingMuEnd_Z      = crsTrkEnd.Z();
    }
    

    // ====================================================================
    // After all track info has been saved, identify a stopping muon-like
    // track in the event
    if( fNumTrackStopping == 1 && stpTrkIndex >= 0 ) {
      
      isOneStoppingTrack = true;
      fNumEvents_StpTrk++;
      
      // -------------------------------------------------------------  
      // Save info.
      fMuTrackIndex       = stpTrkIndex;
      fMuTrackID          = stpTrkID;
      fMuTrackVertex      = stpTrkVertex;
      fMuTrackEnd         = stpTrkEnd;
      fMuTrackVertex_X    = stpTrkVertex.X();
      fMuTrackVertex_Y    = stpTrkVertex.Y();
      fMuTrackVertex_Z    = stpTrkVertex.Z();
      fMuTrackEnd_X       = stpTrkEnd.X();
      fMuTrackEnd_Y       = stpTrkEnd.Y();
      fMuTrackEnd_Z       = stpTrkEnd.Z();
      fMuTrackLength      = stpTrkLength;
        
      // -------------------------------------------------------------  
      // Calculate zenith angle.
      TVector3 vert(0.,1.,0.);
      TVector3 tmp = (stpTrkVertex-stpTrkEnd);
      fMuTrackZenithAngle = tmp.Angle(vert)*RAD_TO_DEG;

      LOG_VERBATIM("MichelAnaFilter")
      <<"Tagged muon track, index = "<<fMuTrackIndex<<", stopping point = ("<<fMuTrackEnd_X<<","<<fMuTrackEnd_Y<<","<<fMuTrackEnd_Z<<")";
        
      // ------------------------------------------------------------
      // Fill truth plots if this is MC to compare to true end/vertex
      if( fTrue_MuTrackEnd_X > 0. ) {
        hTrue_MuTrackVertexRes_X->Fill( (fMuTrackVertex_X - fTrue_MuTrackVertex_X) );
        hTrue_MuTrackVertexRes_Y->Fill( (fMuTrackVertex_Y - fTrue_MuTrackVertex_Y));
        hTrue_MuTrackVertexRes_Z->Fill( (fMuTrackVertex_Z - fTrue_MuTrackVertex_Z));
        hTrue_MuTrackEndRes_X->Fill( (fMuTrackEnd_X - fTrue_MuTrackEnd_X));
        hTrue_MuTrackEndRes_Y->Fill( (fMuTrackEnd_Y - fTrue_MuTrackEnd_Y));
        hTrue_MuTrackEndRes_Z->Fill( (fMuTrackEnd_Z - fTrue_MuTrackEnd_Z));
      }
       
      // ------------------------------------------------------------ 
      // Look at the dEdx of this track
      if( fmcal.isValid() ){
          
        std::vector<art::Ptr<anab::Calorimetry> > calos_mu = fmcal.at(fMuTrackIndex);
        size_t  N_dEdx                = calos_mu[1]->dEdx().size();
        fMuTrackEnergy      = calos_mu[1]->KineticEnergy(); 

        // Save location, dEdx, and residual range vectors from the calorimetry object
        std::vector<TVector3> loc     = calos_mu[1]->XYZ();
        std::vector<double> dEdx      = calos_mu[1]->dEdx();
        std::vector<double> resRange  = calos_mu[1]->ResidualRange();
        fvMuTrackResidualRange = resRange;
        fvMuTrackdEdx = dEdx;
      
        // Find the max dEdx and res range where it peaks
        float dEdx_max    = -9999;
        float resRange_peak = 9999;
        for( size_t i=0; i<N_dEdx; i++ ){ 
          if( dEdx[i] > dEdx_max ) {
            dEdx_max = dEdx[i]; 
            resRange_peak = resRange[i];
          }
        }
        // The point corresponding to the identified muon endpoint (based on 
        // previous fiducial containment arguments) should have the smaller 
        // residual range.  Check if the residual range makes sense based
        // on what we expect.
        float distToEnd_firstPt = (loc[0]-fMuTrackEnd).Mag();
        float distToEnd_lastPt  = (loc[N_dEdx-1]-fMuTrackEnd).Mag();
        float resRange_firstPt  = resRange[0];
        float resRange_lastPt   = resRange[N_dEdx-1];
        bool scenario1 = ((distToEnd_firstPt > distToEnd_lastPt) && (resRange_firstPt > resRange_lastPt));
        bool scenario2 = ((distToEnd_firstPt < distToEnd_lastPt) && (resRange_firstPt < resRange_lastPt)); 
        if( scenario1 || scenario2 ) fMuTrackIsCaloOrdered = true;
        if( fMuTrackIsCaloOrdered && dEdx_max > 6.0 && resRange_peak < 1.0 ) 
          fMuTrackIsBraggPeaked = true;
      }//end if fmcal is valid
      
        
      // -------------------------------------------------------------
      // Finally, create vector of hit keys associated with this track
      art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(TrackHandle, e, fTrackModule);
      if (fmthm.isValid()){
        auto vhit = fmthm.at(fMuTrackID);
        for (size_t h = 0; h < vhit.size(); ++h) fMuTrackHits.push_back(vhit[h].key());
      }
    
    }//<--End if single stopping track
  }//<--end looking at tracks
  
  
  
  // **************************************************************************
  // Optical reconstruction: 
  //  - only perform IF there was either (a) exactly 1 stopping track, or
  //    (b) there was exactly 1 passing track and no other tracks in data.
  // **************************************************************************
  if( (!fReq1StpTrk)||(fReq1StpTrk && fMuTrackIndex >= 0)||(isOneCrossingTrack) ) { 
    
    // Get the PMT data, saving the OpDetPulses specified in fhicl
    if( !fIsMC ) {
      art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
      e.getByLabel("daq","",WaveformHandle);
      for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){
        art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
        raw::OpDetPulse thePulse = *ThePulsePtr;
        for( size_t i=0; i<fSelectChannels.size(); i++){ 
          if( fSelectChannels[i] == thePulse.OpChannel() ) fPulses.push_back(thePulse);
        }
      }
    } 

    // Skip event if we don't have the expected number of PMTs
    if( !fIsMC && fSelectChannels.size() != fPulses.size() ) return false;
  

    // ===============================================================
    // Loop over each PMT separately
    for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){

      size_t ch   = fPulses[ipmt].OpChannel(); 
      LOG_VERBATIM("MichelAnaFilter")<<"Analyzing optical channel "<<ch;
      
      // -----------------------------------------------------
      // If this is MC, do pseudo reco where truth quantities
      // are appropriated into reconstructed quantities
      if( fIsMC ) { 
        
        PerfectWfmReco(fPulses[ipmt]);
        std::vector<float> tmp(fPulses[ipmt].Waveform().begin(), fPulses[ipmt].Waveform().end());
        fPMT_wfm[ch]      = tmp;

      // -----------------------------------------------------
      // Otherwise, do full reconstruction on the PMT waveforms
      } else {
        
        // Set proper gradient hit threshold for hit-finding
        fOpHitBuilderAlg.Reset();
        fOpHitBuilderAlg.SetGradHitThresh( fGradHitThresh[ch] );
        fOpHitBuilderAlg.SetTau( fMuContamCorr_EffTau[ch] );
        
        // Save the PMT waveform
        if( fTruncateWfm == -1 ) fTruncateWfm = fPulses[ipmt].Waveform().size();
        std::vector<short> PMT_wfm_raw(fPulses[ipmt].Waveform().begin(), fPulses[ipmt].Waveform().begin() + fTruncateWfm); 
        std::vector<float> tmp(PMT_wfm_raw.begin(), PMT_wfm_raw.end());
        fPMT_wfm[ch] = tmp;
        fPMT_wfm_raw[ch] = PMT_wfm_raw;
        
        // Smooth out high-freq noise
        fOpHitBuilderAlg.SmoothOutVector(fPMT_wfm[ch],fWfmSmoothingRange);
        
        // Get waveform pedestal/RMS, and subtract off pedestal
        fOpHitBuilderAlg  .CalcBaselineAndRMS( fPMT_wfm[ch], 0, fBaselineWindowLength );
        fOpHitBuilderAlg  .SubtractBaseline( fPMT_wfm[ch] );
        fWfmRMS[ch]       = fOpHitBuilderAlg.GetRMS();

        // PMT overshoot correction
        std::vector<short> hitTimesTmp = fOpHitBuilderAlg.GetHits(fPMT_wfm[ch], (size_t)fPulses[ipmt].FirstSample() );
        if( fCorrectOvershootMode != "" && ch == 0 ) 
          fOpHitBuilderAlg.CorrectWfmOvershoot(fPMT_wfm[ch], hitTimesTmp, fCorrectOvershootMode );
        
        // Masked baseline subtraction
        fvbs[ch].resize(fPMT_wfm[ch].size());
        if( fMskBaselineSubtr ) {
          fOpHitBuilderAlg.MaskedBaselineSubtraction(fPMT_wfm[ch], fvbs[ch]);
          for(size_t i=0; i<fPMT_wfm[ch].size(); i++) fPMT_wfm[ch].at(i) = fPMT_wfm[ch].at(i) - fvbs[ch].at(i);
        }
        
        // Integrate 3 regions in the baseline to get a sense of the inherent spread
        hBaselinePE_100ns[ch]   -> Fill( Integrate(fPMT_wfm[ch],5,5+100) / fSPE[ch] );
        hBaselinePE_500ns[ch]   -> Fill( Integrate(fPMT_wfm[ch],5,5+500) / fSPE[ch] );
        hBaselinePE_1000ns[ch]  -> Fill( Integrate(fPMT_wfm[ch],5,5+1000) / fSPE[ch] );

        // Perform a "first-pass" hit-finding/filtering
        std::vector<short> hitTimes = fOpHitBuilderAlg.GetHits(fPMT_wfm[ch], (size_t)fPulses[ipmt].FirstSample() );
        for (size_t i=0; i< std::min(hitTimes.size(), kMaxOpHits); i++) fvHitTimes0[ch].push_back(hitTimes.at(i));
        LOG_VERBATIM("MichelAnaFilter") << "  found "<<hitTimes.size()<<" hits (thresh = "<<fGradHitThresh[ch]<<")";
        fNumOpHits0[ch] = fvHitTimes0[ch].size();
       
      } 
    }//<-- end first-pass loop over PMTs, doing waveform cleanup and hit-finding

   
    // ======================================================================
    // Now loop through PMT hits and do reconstruction on each found hit (data-only)
    if( !fIsMC ) { 
      for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
        size_t ch=fPulses.at(ipmt).OpChannel();
      
        fOpHitBuilderAlg.Reset();
        fOpHitBuilderAlg.SetTau( fMuContamCorr_EffTau[ch] );
        
        LOG_VERBATIM("MichelAnaFilter") << "Performing hit reconstruction on PMT "<<ch;
    
        // ---------------------------------------------------------------------
        // Loop over the hits, perform reconstruction over each
        std::vector<float> hit_info(fIntegrationWindows.size()+1); 
        std::vector<float> hit_info_raw(fIntegrationWindows.size()+1);
        for (size_t ihit=0; ihit< std::min(fvHitTimes0[ch].size(), kMaxOpHits); ihit++){
          
          short hitTime = fvHitTimes0[ch].at(ihit); 
          short prevHitTime = 0;
          if( fvHitTimes[ch].size() > 0 ) prevHitTime = fvHitTimes[ch].at(fvHitTimes[ch].size()-1);
          
          // --------------------------------------------------------------------
          // Get raw hit integrals; skip hits with raw prompt PE < 0
          hit_info_raw  = fOpHitBuilderAlg.GetHitInfo(fPMT_wfm[ch], hitTime, 0, fIntegrationWindows, false);
          if( hit_info_raw[1] < 0 ) continue;

          // ---------------------------------------------------------------------
          // Use prepulse fit to get hit's width, amplitude, and prefit integrals;
          hit_info      = fOpHitBuilderAlg.GetHitInfo(fPMT_wfm[ch], hitTime, prevHitTime, fIntegrationWindows, true);

          float hitAmplitude = hit_info[0]; // ADC
          float hitPromptPE = hit_info[1]/fSPE[ch]; // PE
          float pulseWidth =  fOpHitBuilderAlg.PulseWidth;
          
          LOG_VERBATIM("MichelAnaFilter") << "    "<<ihit<<"  T = "<<hitTime<<"   amp: "<<hitAmplitude<<"   width = "<<fOpHitBuilderAlg.PulseWidth;
          
          hWidth_AllHits[ch]          ->Fill( pulseWidth );
          hAmplitudeVsWidth_AllHits[ch]->Fill( hitAmplitude, pulseWidth );
          hPromptPEVsWidth_AllHits[ch]->Fill( hitPromptPE,  pulseWidth );
         
          
          //.......................................................................
          // Delete narrow hits
          if( fDeleteNarrowHits && pulseWidth < fMinOpHitWidth[ch] ) {
            LOG_VERBATIM("MichelAnaFilter")<<"    Deleting this hit...";
            for(short i=hitTime - (short)pulseWidth; i< hitTime+(short)pulseWidth; i++){
              float val = -1.*fOpHitBuilderAlg.fit_SlowNorm * exp( -1.*((float)i - fOpHitBuilderAlg.fit_ZeroPoint) / fOpHitBuilderAlg.fit_SlowTau );
              fPMT_wfm[ch].at(i) = val;
            }
            continue;
          }
     
          // .................................................................... 
          // This hit must either be matched to a hit on the other PMT(s) or 
          // its amplitude must be significant enough not to be ignored
          if( fRequireHitMatching && fPulses.size() > 1 ) {
            bool isMatched = true;
            if( hitAmplitude < fOpHitMatchThresh[ch] )  {
              LOG_VERBATIM("MichelAnaFilter") << "    Looking for matching hits on the other PMT(s)...";
              size_t nMatches = 0;
              for(size_t j=0; j<fPulses.size(); j++){
                size_t jch=fPulses.at(j).OpChannel();
                if( jch == ch ) continue;
                // loop through the hits found on this PMT
                for(size_t jhit=0; jhit<fvHitTimes0[jch].size(); jhit++){
                  float dT = fvHitTimes0[ch].at(ihit) - fvHitTimes0[jch].at(jhit);
                  if( fabs( dT ) < fMaxDiff_dT ) {
                    nMatches++;
                    LOG_VERBATIM("MichelAnaFilter")<<"      found match on PMT "<<jch<<" (dt = "<<dT<<")";
                    break;
                  }
                }
              }
              if( nMatches == fPulses.size()-1 ) {
                LOG_VERBATIM("MichelAnaFilter")<<"    found matches on all other PMTs, hit passes!";
              } else {
                isMatched = false;
                hOpHitAmplitude_unmatched[ch]->Fill(hitAmplitude);
                hHitPromptPE_unmatched[ch]->Fill(hitPromptPE);
              }
            }//endif hit > thresh ADC amplitude
            if( !isMatched) continue;
          }// endif > 1 PMTs

         
          // -----------------------------------------------------------------
          // Save hit time information
          fvHitTimes[ch]          .push_back(hitTime);
          fvIsHitAtTrigger[ch]    .push_back( fabs( hitTime - short(fPulses[ipmt].FirstSample()) ) <= fPMT_wfm[ch].size()*0.025 ); 
          hOpHitTime[ch]          ->Fill( hitTime - short(fPulses[ipmt].FirstSample()) );
          
          // -------------------------------------------------------------------
          // Check for saturation
          bool pulseSaturates = false;
          for( int j=0; j<100; j++){ if( fPMT_wfm_raw[ch].at(hitTime+j) == 0 ) pulseSaturates = true; }
          fvIsHitSaturated[ch]    .push_back( pulseSaturates ); 
          
          // ------------------------------------------------------------------
          // Save hit information (from pre-pulse fit reco)
          fvHitAmplitudes[ch]   .push_back(hit_info[0]);
          fvHitWidth[ch]        .push_back(fOpHitBuilderAlg.PulseWidth);
          fvHitADCpf_100ns[ch]  .push_back(hit_info[1]);
          fvHitADCpf_7000ns[ch] .push_back(hit_info[12]);
          fvPrepulseBaseline[ch].push_back(fOpHitBuilderAlg.prepulse_baseline);
          fvPrepulseRMS[ch]     .push_back(fOpHitBuilderAlg.prepulse_rms);
          fvPrepulseSlowNorm[ch].push_back(fOpHitBuilderAlg.fit_SlowNorm);
          fvPrepulseSlowTau[ch] .push_back(fOpHitBuilderAlg.fit_SlowTau); 
          fvPrepulseX1[ch]      .push_back(fOpHitBuilderAlg.prepulse_x1);     
          fvPrepulseZeroPoint[ch].push_back(fOpHitBuilderAlg.fit_ZeroPoint);
          
          // A couple histograms
          hOpHitAmplitude[ch]     ->Fill(hitAmplitude);
          hHitPromptPE[ch]     ->Fill(hitPromptPE);
          
          // Save the raw hit integrals
          fvHitADC_100ns[ch]    .push_back(hit_info_raw[1]);
          fvHitADC_200ns[ch]    .push_back(hit_info_raw[2]);
          fvHitADC_300ns[ch]    .push_back(hit_info_raw[3]);
          fvHitADC_400ns[ch]    .push_back(hit_info_raw[4]);
          fvHitADC_500ns[ch]    .push_back(hit_info_raw[5]);
          fvHitADC_600ns[ch]    .push_back(hit_info_raw[6]);
          fvHitADC_700ns[ch]    .push_back(hit_info_raw[7]);
          fvHitADC_900ns[ch]    .push_back(hit_info_raw[8]);
          fvHitADC_1200ns[ch]   .push_back(hit_info_raw[9]);
          fvHitADC_1500ns[ch]   .push_back(hit_info_raw[10]);
          fvHitADC_1800ns[ch]   .push_back(hit_info_raw[11]);
          fvHitADC_7000ns[ch]   .push_back(hit_info_raw[12]);
           
        }//<-- done looping over PMT hits
    
        LOG_VERBATIM("MichelAnaFilter")<<"    total hits saved: "<<fvHitTimes[ch].size();

      }//<-- end loop over PMTs

    }//<-- end hit reconstruction (DATA ONLY)
    LOG_VERBATIM("MichelAnaFilter")<<"Done reconstructing optical hits.";


    // ======================================================================
    // Redefine number of hits
    for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
      size_t ch   = fPulses[ipmt].OpChannel(); 
      fNumOpHits[ch]      = fvHitTimes[ch].size();
    }

    
    // ======================================================================
    // Now that all hits have been reconstructed, try and identify a Michel-like
    // decay topology
    if( !fReq1StpTrk || (fReq1StpTrk && fMuTrackIndex >= 0 ) ) {
      for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
        size_t ch   = fPulses[ipmt].OpChannel(); 

        // ------------------------------------------------------
        // Calc dT if there are two hits
        if( fNumOpHits[ch] == 2 ) {
          fdT[ch] = float( fvHitTimes[ch].at(1)-fvHitTimes[ch].at(0) );
          hdT[ch] ->Fill( fdT[ch] );
          hdT_vs_Amp[ch]->Fill( fdT[ch], fvHitAmplitudes[ch].at(1) );
          LOG_VERBATIM("MichelAnaFilter") << "Exactly 2 hits in PMT "<<ch<<"  dT = " << fdT[ch] << " ns.";
        
          // -------------------------------------------------------
          // If the SPE calibration for this PMT is known, calculate
          // photoelectron (PE) quantities
          if( fSPE[ch] > 0 ) { 
            
            // Define "raw" PE quantities
            fPE_promptRaw[ch] = fvHitADC_100ns[ch].at(1) / fSPE[ch];
            fPE_totalRaw[ch]  = fvHitADC_7000ns[ch].at(1) / fSPE[ch];
              
            // Muon late-light correction:
            // If dT > dTmin, do a fit to the muon late-light and subtract 
            // off the added PEs to the Michel pulse. Otherwise just use 
            // the integrals that used the pre-pulse fits from before.
            if( fMuContamCorr && fdT[ch] >= fMuContamCorr_dTmin ) {
              MuContamCorr(ch); //<-- this function defines fPE_prompt, fPE_total
            } else {
              fPE_prompt[ch]    = fvHitADCpf_100ns[ch].at(1) / fSPE[ch];
              fPE_total[ch]     = fvHitADCpf_7000ns[ch].at(1) / fSPE[ch];
            }//<-- end muon late-light correction
            
            //LOG_VERBATIM("MichelAnaFilter")
            //<<"  Prompt PE = "<<fPE_prompt[ch]<<"   Total PE = "<<fPE_total[ch];

          } //<-- endif SPE > 0

          // -----------------------------------------------------------------------
          // Apply smearing on MC to approximate resolution associated with detector 
          // effects and mu late light contamination correction.  "Natural" smearing 
          // has already been done on the simulated ADC integrals based on integration 
          // time, waveform RMS (from data), and SPE resolution.  Here we simply add 
          // an extra smearing contribution to the final total light integral based on 
          // optimized fits to data.
          if( fIsMC && fSmearFactor[ch] > 0 && fPE_total[ch] > 0 ) 
            fPE_total[ch]   += fRand->Gaus(0, fSmearFactor[ch]*sqrt(fPE_total[ch]));
          
          // Fill width vs. amplitude and PE histograms
          if( fdT[ch] >= fGateDelay ) {
            hWidth[ch]     ->Fill(fvHitWidth[ch].at(1));
            hAmplitudeVsWidth[ch]->Fill(fvHitAmplitudes[ch].at(1)*fMvPerADC,  fvHitWidth[ch].at(1));
            hPromptPEVsWidth[ch]->Fill(fPE_prompt[ch],                 fvHitWidth[ch].at(1));
          }

        } // endif nHits = 2

        // ---------------------------------------------------------------------- 
        // If dT was defined, and the second hit is the "triggered"
        // hit (in time), then we have muon and Michel candidates!
        // Also require BOTH hits to be above the width threshold.
        if(     fdT[ch] > 0 
            &&  isOneStoppingTrack
            &&  fvIsHitAtTrigger[ch].at(1) 
            &&  fvHitWidth[ch].at(0) >= fMinOpHitWidth[ch]
            &&  fvHitWidth[ch].at(1) >= fMinOpHitWidth[ch]  
            &&  fSPE[ch] > 0 ){
          
          fMuAmplitude[ch]  = fvHitAmplitudes[ch].at(0);
          fAmplitude[ch]    = fvHitAmplitudes[ch].at(1);
          fMuPE_prompt[ch]  = fvHitADC_100ns[ch].at(0) / fSPE[ch];
          fWidth[ch]        = fvHitWidth[ch].at(1);
          fMuWidth[ch]      = fvHitWidth[ch].at(0);
          fMuPulseSaturated[ch] = fvIsHitSaturated[ch].at(0);
          fElPulseSaturated[ch] = fvIsHitSaturated[ch].at(1);
          numMichelWfms[ch]++;
          numMuSaturated[ch] += (size_t)fMuPulseSaturated[ch];
          numElSaturated[ch] += (size_t)fElPulseSaturated[ch];
          LOG_VERBATIM("MichelAnaFilter")
          <<"  --> Michel candidate wfm found!  PE(prompt)= "<<fPE_prompt[ch]<<", PE(total)= "<<fPE_total[ch];
         
          // Replicate trigger efficiency in MC by scuplting the prompt PE dist
          // and cutting the event out based on some probability.
          if( fIsMC && fApplyTrigEffCut && fTrigEff_P[ch] > 0 ) {
            TF1 probCut("probCut","1./(1.+(x/[0])^[1] )");
            probCut.SetParameter(0,fTrigEff_P[ch]);
            probCut.SetParameter(1,fTrigEff_K[ch]);
            std::cout<<probCut.GetParameter(0)<<"   "<<probCut.GetParameter(1)<<"\n";
            float r = fRand->Rndm();
            float p = probCut.Eval(fPE_prompt[ch]);
            //std::cout<<"  prompt PE = "<<fPE_prompt[ch]<<"\n";
            //std::cout<<"  p         = "<<p<<"\n";
            //std::cout<<"  r         = "<<r<<"\n"; 
            if( r < p ) {
              LOG_VERBATIM("MichelAnaFilter")<<"MC event marked for removal via trigger eff. cut.";
              return false;
            }
          }
       
          // Define the prompt light fraction and fill some histograms 
          if( fPE_totalRaw[ch] > 0 ) fPromptFracRaw[ch] = fPE_promptRaw[ch] / fPE_totalRaw[ch];
          if( fPE_total[ch] > 0 ) {
            fPromptFrac[ch]   = fPE_prompt[ch] / fPE_total[ch];
            hTrue_MuLateLightContamination->Fill( fTrue_dT, fTrue_MuContam_total[ch] / fPE_totalRaw[ch] );
          }

          // Save average waveforms and PE arrival profiles
          if( fdT[ch] > 4050 ) {
              if( !fIsMC && fvHitTimes[ch].at(0) > 500 && fvHitTimes[ch].at(1)+7000 < (int)fPMT_wfm[ch].size() ) {
                AddToAverageWfm( fPMT_wfm[ch], fvHitTimes[ch].at(0) - 400, fvHitTimes[ch].at(0) + 4000, hAveWfm_mu[ch], aveWfmCounter_mu[ch], -1);
                AddToAverageWfm( fPMT_wfm[ch], fvHitTimes[ch].at(1) - 700, fvHitTimes[ch].at(1) + 7000, hAveWfm_el[ch], aveWfmCounter_el[ch], -1);
              } else
              if( fIsMC && fvHitTimes[ch].at(1)+7000 < (int)hPMT_phelTimes[ch]->GetNbinsX() ) {
                AddToAverageWfm( hPMT_phelTimes[ch], fvHitTimes[ch].at(0) - 400, fvHitTimes[ch].at(0) + 4000, fvHitTimes[ch].at(0), hAveWfm_mu[ch], aveWfmCounter_mu[ch], -1);
                AddToAverageWfm( hPMT_phelTimes[ch], fvHitTimes[ch].at(1) - 700, fvHitTimes[ch].at(1) + 7000, fvHitTimes[ch].at(1), hAveWfm_el[ch], aveWfmCounter_el[ch], -1);
              }
          }

        }// end Michel+muon candidate event

      }//<-- end loop over PMTs to search for Michel candidates
      

      // ---------------------------------------------------------------------------
      // Evaluate event reduction cuts on the light-based data and determine
      // the muon decay time if all optical detectors show consistent results
      hOpEventCuts->Fill(1-1);
      
      if( fSelectChannels.size() > 0 && fPulses.size() == fSelectChannels.size() ) {
        hOpEventCuts->Fill(2-1);

        // Check RMS of waveforms 
        bool goodWfms = true;
        for( auto & p : fPulses ) {
          if( fWfmRMS[p.OpChannel()] > fMaxWfmRMS[p.OpChannel()] ) goodWfms = false;
        }

        if( goodWfms ) {
          hOpEventCuts->Fill(3-1);

          // Check if all detectors saw exactly two hits, with consistent dT
          float maxdiff = -999.;
          float sum_dT = 0.;
          float sum_T = 0.;
          float sum_pe = 0.;
          size_t N = 0;
          for( auto & p : fPulses ) {
            if( fdT[p.OpChannel()] > 0. ) { 
              sum_T   += fvHitTimes[p.OpChannel()].at(1);
              sum_dT  += fdT[p.OpChannel()];
              sum_pe  += fPE_total[p.OpChannel()];
              N++; 
            }
            for( auto & p2 : fPulses ) {
              if( p.OpChannel() == p2.OpChannel()) continue; // avoid matching with same PMT
              float diff = fabs(fdT[p.OpChannel()] - fdT[p2.OpChannel()]);
              if( diff > maxdiff ) maxdiff = diff;
            }
          }
          if( maxdiff >= 0 ) hDiff_dT->Fill(maxdiff);

          // did all PMTs see 2 hits?
          if( N == fSelectChannels.size() ) {
            hOpEventCuts->Fill(4-1);
             
            // is time matching ok? (only applicable if N = NPmts > 1)
            if( (N==1) || (maxdiff <= fMaxDiff_dT) ) {

              float decayTime = sum_dT / N;
              hOpEventCuts->Fill(5-1);
        
              // Perform checks on 2nd pulse timing rel. to trigger,
              // saturation, and hit widths, prompt PE cut (for dT plots) 
              bool triggerCut = true;
              bool satCut     = true;
              bool widthCut   = true;
              for( auto & p : fPulses ) {
                size_t ch = p.OpChannel();
                if( !fIsMC ) { if( !fvIsHitAtTrigger[ch].at(1) ) triggerCut = false; }
                if( fElPulseSaturated[ch] ) satCut = false;
                if(    fMuWidth[ch] < fMinOpHitWidth[ch] 
                    || fWidth[ch]   < fMinOpHitWidth[ch]   ) widthCut = false;
                if( fPE_prompt[ch] < fPromptPECut[ch] ) promptPEcut = false;
              }   

                  
              // Check hit widths
              if( widthCut ) {
                hOpEventCuts->Fill(6-1);
                  
                if( triggerCut ) {
                  hOpEventCuts->Fill(7-1);
                
                  // Check for saturation of Michel candidate pulse
                  if( satCut ) {
                    hOpEventCuts->Fill(8-1);
              
                    if( decayTime >= fGateDelay && decayTime < fMaxDecayTime ) {
                      hOpEventCuts->Fill(9-1);
                      fMichelOpticalID = true;
                      fNumOpticalID++;
                      fDecayTime  = decayTime;
                      fPE         = sum_pe;
                      LOG_VERBATIM("MichelAnaFilter")<<"***OPTICAL MICHEL ID***";
                    }// dT > min
                  }// saturation cut
                }// trigger time cut
              }// width cut
            }// dTs match
          }//endif two hits
        }//endif rms cuts
      }//endif all wfms found
      
      
      // Make some light-related histograms
      if( fPulses.size() == 2 ) {
        hNumOpHitsCompare ->Fill(fNumOpHits[0],fNumOpHits[1]);
        if( fMichelOpticalID && fDecayTime > fdTcut) hPECompare ->Fill(fPE_total[0], fPE_total[1]);
      }
      if( fMichelOpticalID ) hPE->Fill(fPE);
      for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
        size_t ch   = fPulses[ipmt].OpChannel(); 
        hWfmRMS[ch]         ->Fill(fWfmRMS[ch]);
        hNumOpHits[ch]      ->Fill(fNumOpHits[ch]);
        if( fMichelOpticalID ) {
          hPE_prompt[ch]      ->Fill(fPE_prompt[ch]);
          hPE_total[ch]           ->Fill(fPE_total[ch]);
          hdT_vs_PE_total[ch]     ->Fill(fDecayTime, fPE_total[ch]);
          hdT_vs_PE_totalRaw[ch]  ->Fill(fDecayTime, fPE_totalRaw[ch]);
          hdT_vs_PE_prompt[ch]    ->Fill(fDecayTime, fPE_prompt[ch]);
          hdT_vs_PE_promptRaw[ch] ->Fill(fDecayTime, fPE_promptRaw[ch]);
          hTrue_PE[ch]        ->Fill(fTrue_PE[ch]);
          hTrue_PE_prompt[ch] ->Fill(fTrue_PE_prompt[ch]);
          hTrue_PE_total[ch]  ->Fill(fTrue_PE_total[ch]);
          hAmplitude[ch]      ->Fill(fAmplitude[ch]);
          hPromptPEVsAmplitude[ch]->Fill(fPE_prompt[ch],fAmplitude[ch]);
          hMuPromptPEVsAmplitude[ch]->Fill(fMuPE_prompt[ch],fMuAmplitude[ch]);
          hPromptFrac[ch]->Fill(fPromptFrac[ch]);
          hPromptFracRaw[ch]->Fill(fPromptFracRaw[ch]);
          if( fDecayTime > fdTcut ) {
            hPromptFrac_dTcut[ch]->Fill(fPromptFrac[ch]);
            hPromptFracRaw_dTcut[ch]->Fill(fPromptFracRaw[ch]);
            hPE_total_dTcut[ch]->Fill(fPE_total[ch]);
            hPE_prompt_dTcut[ch]->Fill(fPE_prompt[ch]);
          }
        }
      }//<-- endloop PMTs
      
      // Decay time histogram
      if( promptPEcut &&  fDecayTime >= fGateDelay && fDecayTime <= fMaxDecayTime ) hDecayTime->Fill(fDecayTime);

    }//<-- endif 1stptrk
    LOG_VERBATIM("MichelAnaFilter")<<"Done looping over PMTs for Michel optical ID.";



    // ==================================================================
    // If there's exactly 1 passing track, and 1 optical hit on each PMT,
    // then use this crossing muon event for calibration purposes
    if( fUseCrossingMuons && isOneCrossingTrack ) {
      
       
      // Check that both PMTs saw *one* unsaturated hit at nearly the same time
      float maxdiff = -999.;
      size_t N = 0;
      for( auto & p : fPulses ) {
        if( fNumOpHits[p.OpChannel()] == 1 ) {
          if( !fvIsHitSaturated[p.OpChannel()].at(0) ) {
            N++;
            for( auto & p2 : fPulses ) {
              if( p2.OpChannel() == p.OpChannel() ) continue;
              if( fNumOpHits[p2.OpChannel()] == 1 && !fvIsHitSaturated[p2.OpChannel()].at(0) ) {
                float diff = fabs( fvHitTimes[p.OpChannel()].at(0) - fvHitTimes[p2.OpChannel()].at(0) ); 
                if( diff > maxdiff ) maxdiff = diff;
              }
            }
          }
        }
      }
      if( N == fSelectChannels.size() ) {
        if( (N==1) || (maxdiff <= fMaxDiff_dT) ) {
          isCalibrationEvent = true;
          fNumEvents_Calibration++;
        }
      }

    }//<-- endif single PASSING track



  }//<-- endif: require single stopping or crossing track


  
  // **************************************************************************
  // Optical hit filtering mode: if turned on, then don't proceed to 
  // clustering stages and just pass/reject the event based on whether
  // a Michel electron was identified optically.
  // **************************************************************************
  if( fFilter_OpticalMode ) {
    if( fMichelOpticalID ) {
      LOG_VERBATIM("MichelAnaFilter")<<"Event PASSES optical filter.";
      return true;
    } else {
      return false;
    }
  }



  // **************************************************************************
  // IF this event is interesting to us based on the reconstruction up to now, 
  // save the wire hit information and apply time/offset and calibration scaling.
  // **************************************************************************
  if( (fMichelOpticalID && fMuTrackIndex >= 0) || ( isCalibrationEvent ) ){
   
    fTotalCharge      = 0;
    fTotalCharge_Pl0  = 0;

    // Get recob::Hit information
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (e.getByLabel(fHitsModule,fHitsInstance,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}
    
    // Assign a location to each hit, relative to projections on wire plane
    size_t  nWireHits = hitlist.size();
    fHitPlane     .reserve( nWireHits );
    fHitWire      .reserve( nWireHits );
    fHitT         .reserve( nWireHits );
    fHitX         .reserve( nWireHits );
    fHitW         .reserve( nWireHits );
    fHitCharge    .reserve( nWireHits );
    fHitIntegral  .reserve( nWireHits );
    for(size_t i=0; i<nWireHits; i++){
      float hitTime = fSamplingRate*(hitlist[i]->PeakTime()-fXTicksOffset[hitlist[i]->WireID().Plane]);
      fHitT         .push_back(hitTime); // us
      fHitX         .push_back(hitTime*fDriftVelocity);
      fHitWire      .push_back(hitlist[i]->Channel());
      fHitPlane     .push_back(hitlist[i]->WireID().Plane);
      if( hitlist[i]->Channel() < 240 ) fHitW.push_back(hitlist[i]->Channel()*0.4);
      else                              fHitW.push_back((hitlist[i]->Channel()-240)*0.4);
      fHitIntegral  .push_back( hitlist[i]->Integral() );
      float ltCorFac = 1.;
      if( hitTime >= 0 ) ltCorFac = exp(hitTime/fElectronLifetime);
      fHitCharge    .push_back( (hitlist[i]->Integral()*ltCorFac) / fCalAreaConstants[hitlist[i]->WireID().Plane]);
      hHitX[fHitPlane.at(i)]         ->Fill( fHitX.at(i) );
      hHitRMS[fHitPlane.at(i)]       ->Fill( hitlist[i]->RMS() );
      hHitAmplitude[fHitPlane.at(i)]  ->Fill( hitlist[i]->PeakAmplitude() );
      if( fHitPlane.at(i) == 1 ) {
        fTotalCharge  += fHitCharge.at(i);
        hHitIntegral  ->Fill(hitlist[i]->Integral());
        hHitCharge    ->Fill(fHitCharge.at(i));
        hHitSummedADC ->Fill(hitlist[i]->SummedADC());
      }
      if( fHitPlane.at(i) == 0 ) fTotalCharge_Pl0  += fHitCharge.at(i);
      fNumPlaneHits[ fHitPlane.at(i) ]++;
    }
    hTotalCharge->Fill( fTotalCharge );

    LOG_VERBATIM("MichelAnaFilter")
    <<"Total charge hits in this event: "<<fHitX.size()<<" ("<<fNumPlaneHits[0]<<" ind, "<<fNumPlaneHits[1]<<" coll)";
    
  }

  

  // **************************************************************************
  // If this is a calibraiton event, get track length, overall track charge
  // deposition (after lifetime corrections), and visibility-corrected light
  // **************************************************************************
  if( isCalibrationEvent ) {
    
    art::Handle< std::vector< recob::Track >> TrackHandle;
    e.getByLabel(fTrackModule,TrackHandle);
    art::Ptr< recob::Track > CrsTrackPtr(TrackHandle,fCrossingMuTrackIndex);
    recob::Track CrsTrack = *CrsTrackPtr;
    
    hCrsTrackNode_ZX->Fill( fCrossingMuVertex_Z, fCrossingMuVertex_X );
    hCrsTrackNode_ZX->Fill( fCrossingMuEnd_Z, fCrossingMuEnd_X );
    hCrsTrackNode_ZY->Fill( fCrossingMuVertex_Z, fCrossingMuVertex_Y );
    hCrsTrackNode_ZY->Fill( fCrossingMuEnd_Z, fCrossingMuEnd_Y );


    fCrossingMuCharge = fTotalCharge;
    /*
    // Get the total charge by adding up charge from all the associated hits
    // Track pitch
    float trkPitch = lar::util::TrackPitchInView(CrsTrack, geo::kV);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(TrackHandle, e, fTrackModule);
    std::vector<int> muTrackHits;
    if (fmthm.isValid()){
      auto vhit = fmthm.at(0);
      for (size_t h = 0; h < vhit.size(); ++h) muTrackHits.push_back(vhit[h].key());
    }
    fCrossingMuCharge = 0;
    for(size_t i=0; i<muTrackHits.size(); i++){
      if( fHitPlane.at( muTrackHits.at(i) ) == 1 ){
        fCrossingMuCharge += fHitCharge.at( muTrackHits.at(i) );
        hdQdxVsT    ->Fill( fHitT.at(muTrackHits.at(i)), fHitIntegral.at(muTrackHits.at(i))/trkPitch );
        hdQdxVsT_N  ->Fill( fHitT.at(muTrackHits.at(i)) );
      }
    }
    */

    // Find average optical visibility of the 3D traj. pts within the track
      float sum_pe    = 0;
      float sum_pe_pr = 0;
      float sum_vis   = 0.;
      for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
        size_t ch   = fPulses[ipmt].OpChannel(); 
        sum_pe    += fvHitADC_7000ns[ch].at(0) / fSPE[ch];
        sum_pe_pr += fvHitADC_100ns[ch].at(0) / fSPE[ch];
        for(size_t i=0; i<CrsTrack.NumberTrajectoryPoints(); i++){
          TVector3 point = CrsTrack.LocationAtPoint(i);
          sum_vis += GetVisibility( point , ch); 
        }
      }
      float mean_vis = sum_vis / float( CrsTrack.NumberTrajectoryPoints() );
      if( mean_vis > 0 && sum_pe > 0 ) {
        fCrossingMuPhotons        = sum_pe / mean_vis;
        fCrossingMuPhotonsPrompt  = fCrossingMuPhotons * (sum_pe_pr / sum_pe);
      }

      // Recomb. (exclude muons that pass through wire planes)
      if( fCrossingMuCharge > 0 && fCrossingMuPhotons > 0 && fCrossingMuLength > 20. 
        && fCrossingMuVertex_X > 2. && fCrossingMuEnd_X > 2. ){
        float QoverL = fCrossingMuCharge / fCrossingMuPhotons;
        float R = ( 1.+fExcRatio ) * QoverL / ( 1.+QoverL );
        hQoverLCrsMuTrk->Fill( QoverL );
        hRecombCrsMuTrk->Fill( R );
        hCrsMuQPerCm->Fill( fCrossingMuCharge / fCrossingMuLength );
        hCrsMuLPerCm->Fill( fCrossingMuPhotons / fCrossingMuLength );
      }


  }//<-- end calibration event information



  // *******************************************************************
  // Michel electron charge/shower reconstruction. If a stopping muon
  // track was identified AND a double-pulse was seein the optical 
  // waveforms, do Michel electron charge clustering in 2D/3D.
  // ******************************************************************
  
  // Define michelcluster objects for each plane
  michelcluster cl_pl0;
  michelcluster cl_pl1;

  if( fMuTrackIndex >= 0 && fMichelOpticalID ) {  

    // ===========================================================
    // Perform clustering on both wireplanes
    Clustering( 0, cl_pl0 );
    Clustering( 1, cl_pl1 );

    // ===========================================================
    // Save params to class member variables
    //...collection plane (default energy plane)
    fMinCovAtBnd        = cl_pl1.minCovAtBnd;
    fClusterSize        = cl_pl1.cluster.size();
    fMuClusterSize      = cl_pl1.cluster_mu.size();
    fElClusterSize      = cl_pl1.cluster_el.size();
    fFracMuHitsLinear   = cl_pl1.fracMuHitsLinear;
    fMuAveLinearity     = cl_pl1.muAveLinearity;
    fMuClusterHitsEndFit= cl_pl1.nPtsMuFit;
    fExtraHits          = cl_pl1.extraHits;
    fDecayAngle2D       = cl_pl1.decayAngle2D;
    fElShowerAngleMean    = cl_pl1.showerAngleMean;
    fElShowerAngleRMS     = cl_pl1.showerAngleRMS;
    fElShowerSize       = cl_pl1.shower.size();
    fElShowerFrac       = cl_pl1.elShowerFrac;
    //...induction plane
    fMinCovAtBnd_Pl0        = cl_pl0.minCovAtBnd;
    fClusterSize_Pl0        = cl_pl0.cluster.size();
    fMuClusterSize_Pl0      = cl_pl0.cluster_mu.size();
    fElClusterSize_Pl0      = cl_pl0.cluster_el.size();
    fFracMuHitsLinear_Pl0   = cl_pl0.fracMuHitsLinear;
    fMuAveLinearity_Pl0     = cl_pl0.muAveLinearity;
    fMuClusterHitsEndFit_Pl0= cl_pl0.nPtsMuFit;
    fExtraHits_Pl0          = cl_pl0.extraHits;
    fDecayAngle2D_Pl0       = cl_pl0.decayAngle2D;
    fElShowerAngleMean_Pl0    = cl_pl0.showerAngleMean;
    fElShowerAngleRMS_Pl0     = cl_pl0.showerAngleRMS;
    fElShowerSize_Pl0       = cl_pl0.shower.size();
    fElShowerFrac_Pl0       = cl_pl0.elShowerFrac;
    
    LOG_VERBATIM("MichelAnaFilter")
    <<"DONE CLUSTERING: \n"
    <<" coll. plane cluster size = "<<fClusterSize<<" (mu: "<<fMuClusterSize<<", el: "<<fElClusterSize<<")\n"
    <<" frac mu hits linear = "<<fFracMuHitsLinear<<"\n"
    <<" el shower completeness frac = "<<fElShowerFrac;
   
    // Average charge-weighted drift time and x-position of shower hits 
    fAveDriftTime       = cl_pl1.aveDriftTime;
    fAveX               = cl_pl1.aveX;
   
    // If a cluster boundary was found, set flag
    if( cl_pl1.bnd_i >= 0 ) isClusterBndFound = true;
    
    // ===========================================================
    // Measure the drift offset (X) of the first electron hit 
    // from the last muon hit
    if( isClusterBndFound && fMuClusterSize > 0 && fElClusterSize > 0 ){ 
      fElStartWire = fHitWire[ cl_pl1.cluster_el.at(0) ];
      fMuEndWire   = fHitWire[ cl_pl1.muEndHit ];
      fElOffsetT = fHitT[ cl_pl1.cluster_el.at(0) ] - fHitT[ cl_pl1.muEndHit ];
      if( fElStartWire == fMuEndWire ) 
        hElOffsetT ->Fill(fElOffsetT);
    }
  
    // ===========================================================
    // If cluster boundary (ie, muon end hit) found on both planes, then 
    // we can calculate the 3D muon endpoint position and get an idea of 
    // any systematic offsets in timing between the planes
    if( cl_pl1.muEndHit >= 0 && cl_pl0.muEndHit >= 0 ) {
      int hit0 = cl_pl0.muEndHit;
      int hit1 = cl_pl1.muEndHit;
      fMuEndHitTimeDiff = fHitT[hit1] - fHitT[hit0];
      hMuEndHitTimeDiff -> Fill( fMuEndHitTimeDiff );
      double y = -99.;
      double z = -99.;
      fGeo->ChannelsIntersect( fHitWire.at(hit0), fHitWire.at(hit1), y, z );
      //fMuEnd3D.SetXYZ(fHitX.at(hit1), y, z);
      fMuEnd3D_X = fHitX.at(hit1);
      fMuEnd3D_Y = y;
      fMuEnd3D_Z = z;
      hMuEnd3D_ZX->Fill( fMuEnd3D_Z, fMuEnd3D_X );
      hMuEnd3D_ZY->Fill( fMuEnd3D_Z, fMuEnd3D_Y );
    }
    
    
    // ===========================================================
    // Reconstruct an energy on both planes
    if( fDecayTime >= 0. && isOneStoppingTrack ) {
     
      // ...collection plane 
      if( fElShowerSize > 0 ) {
        CalcMichelShowerEnergy( cl_pl1 );
        fElShowerCharge = cl_pl1.elShowerCharge;
        fElShowerEnergy = cl_pl1.elShowerEnergy;
        fElTrackCharge  = cl_pl1.elTrackCharge;
        fElTrackEnergy  = cl_pl1.elTrackEnergy;
      }
      // ...induction plane
      if( fElShowerSize_Pl0 > 0 ) {
        CalcMichelShowerEnergy( cl_pl0 );
        fElShowerCharge_Pl0 = cl_pl0.elShowerCharge;
        fElShowerEnergy_Pl0 = cl_pl0.elShowerEnergy;
        fElTrackCharge_Pl0  = cl_pl0.elTrackCharge;
        fElTrackEnergy_Pl0  = cl_pl0.elTrackEnergy;
      }

    }//<-- end cut on dT>0, 1 stopping trk

    LOG_VERBATIM("MichelAnaFilter")
    <<"Calculated shower energy on coll. plane = "<<fElShowerEnergy; 
    
    
    // ===========================================================
    // Calculate energy resolution 
    float res = 0.;
    if( fIsMC && fElShowerEnergy > 0 ) {
      res = (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep;
      LOG_VERBATIM("MichelAnaFilter")<<"Energy resolution: "<<res;
    }

   
    
    // ===========================================================
    // Match hits between the two planes to form 3D points 
    // ===========================================================
    if( fElShowerSize > 0 && fElShowerSize_Pl0 > 0 ) {
        
        std::vector<int> shower1 = cl_pl1.shower;
        std::vector<int> shower0 = cl_pl0.shower;
        std::vector<bool> flag_1(cl_pl1.shower.size(), true);
        std::vector<bool> flag_0(cl_pl0.shower.size(), true);
        LOG_VERBATIM("MichelAnaFilter")
        <<"Beginning to match hits between two planes: "<<shower1.size()<<", "<<shower0.size()<<" hits each";
       
        // ---------------------------------------------------
        // Create array of hit time differences 
        std::vector< std::vector<float> > diff_array (shower1.size(), std::vector<float> (shower0.size(), 0));
        for(size_t i=0; i<shower1.size(); i++){
          for(size_t j=0; j<shower0.size(); j++){
            diff_array.at(i).at(j) = fHitT[shower1.at(i)] - fHitT[shower0.at(j)]; 
            hHitTimeDiff->Fill(diff_array.at(i).at(j));
          }
        }
         
        // ---------------------------------------------------
        // Create vectors to hold the 3D points (TVectors) and 
        // lists of hits from each plane associated with the points.
        std::vector<TVector3> Pts;
        std::vector<int>      Pts_hitKeys1;
        std::vector<int>      Pts_hitKeys0;

        // ---------------------------------------------------
        // Loop over all the hit differences and starting with the 
        // best-matched hits, calculate their XYZ point based on a
        // function in the geometry class. Do this until there are
        // no more matches to be made. Note that hits are required
        // to match by at least less than 2us.
        bool loopHits = true;
        while ( loopHits ) { 
          int min_i = -9;
          int min_j = -9;
          float minDiff = 99999.;
          for(size_t i=0; i<shower1.size(); i++){
            if( !flag_1.at(i) ) continue;
            for(size_t j=0; j<shower0.size(); j++){
             if( !flag_0.at(j) ) continue;
             if(       fabs(diff_array.at(i).at(j)) < fMaxHitTimeDiff
                   &&  fabs(diff_array.at(i).at(j)) < minDiff ) {
                 minDiff = diff_array.at(i).at(j);
                 min_i = i;
                 min_j = j;
               }
            }
          }
          if( min_i >= 0 ) {
            // Find y-z coordinate for this pair of hits
            double y = -99.;
            double z = -99.;
            int hit1 = shower1.at(min_i);
            int hit0 = shower0.at(min_j);
            fGeo->ChannelsIntersect( fHitWire.at(hit0), fHitWire.at(hit1), y, z );
            TVector3 pt( fHitX.at(hit1), y, z );
            Pts         .push_back(pt);
            Pts_hitKeys1 .push_back(hit1);
            Pts_hitKeys0 .push_back(hit0);
            flag_1.at(min_i) = false;
            flag_0.at(min_j) = false;
          } else {
            loopHits = false; 
          }
        }

       
        // ------------------------------------------------------
        // Count up the 3D points and calculate the fraction of points
        // on either plane that were successfully matched. 
        fNumPts3D   = Pts.size();
        fFracHits3D = float(Pts.size()) / float(std::min( fElShowerSize, fElShowerSize_Pl0 ));
        LOG_VERBATIM("MichelAnaFilter")
        <<"Found "<<Pts.size()<<" 3D points (frac = "<<fFracHits3D<<")";
      

        // ====================================================================
        // If we have some 3D points reconstructed, we can do more stuff...

        if( fNumPts3D > 0 && fElShowerEnergy > 0 ) { 
        
          // ----------------------------------------------------------------- 
          // Attempt to measure dE/dx along the electron portion of the track
          // (NOTE: This doesn't work so well yet!)
          
          // First, require that some number of 3D points are associated with
          // hits in the electron cluster (ie, the actual Michel el "track")
          // on each plane.
          std::vector<TVector3> PtsTrk;
          std::vector<int>      PtsTrk_hitKeys1;
          std::vector<int>      PtsTrk_i;
          for(size_t i=0; i<Pts.size(); i++){
            bool isInPl1Cluster = ( std::find(cl_pl1.cluster_el.begin(), cl_pl1.cluster_el.end(), Pts_hitKeys1.at(i) ) != cl_pl1.cluster_el.end() );
            bool isInPl0Cluster = ( std::find(cl_pl0.cluster_el.begin(), cl_pl0.cluster_el.end(), Pts_hitKeys0.at(i) ) != cl_pl0.cluster_el.end() );
            if( isInPl1Cluster && isInPl0Cluster ) {
              PtsTrk            .push_back( Pts.at(i) );
              PtsTrk_hitKeys1   .push_back( Pts_hitKeys1.at(i) );
            }
          }
          fNumPts3DTrk = PtsTrk.size();

          // If every hit in the 2D collection plane Michel cluster
          // was turned into a 3D point, we can calculate dE/dx
          if( fNumPts3DTrk == (int)cl_pl1.cluster_el.size() && fNumPts3DTrk >= 3 ) {
            
            // First sort the 3D electron track points based on their
            // sequence in the 2D cluster
            std::vector<TVector3> PtsTrk_sorted;
            std::vector<int>      PtsTrk_hitKeys1_sorted;
            for(size_t i=0; i<cl_pl1.cluster_el.size(); i++){
              for(int j=0; j<fNumPts3DTrk; j++) {
                if( cl_pl1.cluster_el.at(i) == PtsTrk_hitKeys1.at(j) ) {
                  PtsTrk_sorted.push_back( PtsTrk.at(j) );
                  PtsTrk_hitKeys1_sorted.push_back( PtsTrk_hitKeys1.at(j) );
                }
              }
            }
            PtsTrk          = PtsTrk_sorted;
            PtsTrk_hitKeys1 = PtsTrk_hitKeys1_sorted;
           
            // Find ave dE/dx of the electron track
            float dx_tot = 0.;
            for(int i=1; i<fNumPts3DTrk; i++) dx_tot += (PtsTrk.at(i) - PtsTrk.at(i-1)).Mag();
            if( dx_tot > 0. ) {
              fElTrackdEdx = cl_pl1.elTrackEnergy / dx_tot;
            }
             
        } //<-- endif fNum3DPtsTrk = cluster size >= 3 (electron dE/dx calculation)


        // -------------------------------------------------------------
        // Find average charge-weighted optical visibility of the 3D points
        fElShowerPhotons = 0;
        fElShowerVis = 0.;
        float sum_W = 0.;
        float sum_x = 0.;
        float sum_y = 0.;
        float sum_z = 0.;
        for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
          size_t ch   = fPulses[ipmt].OpChannel(); 
          fElShowerVisCh[ch] = 0.;
          float sum_W_ch = 0.;
          float sum_vis = 0.;
          for(size_t i=0; i<Pts.size(); i++){
            float W = fHitCharge.at( Pts_hitKeys1.at(i) );
            //float W = 1.;
            sum_W   += W;
            sum_W_ch += W;
            sum_vis += W* GetVisibility(Pts.at(i), ch);
            sum_x   += W* Pts.at(i).X();
            sum_y   += W* Pts.at(i).Y();
            sum_z   += W* Pts.at(i).Z();
          }
          fElShowerVisCh[ch]  =  sum_vis / sum_W_ch;
          fElShowerVis        += sum_vis / sum_W_ch;
          LOG_VERBATIM("MichelAnaFilter")<<"Ave vis for PMT "<<ch<<" = "<<sum_vis / sum_W;
        }//<-- end loop over PMTs
        if( sum_W > 0 ) {
          fElShowerCentroid_X = sum_x / sum_W;
          fElShowerCentroid_Y = sum_y / sum_W;
          fElShowerCentroid_Z = sum_z / sum_W;
          //fElShowerCentroid.SetXYZ(fElShowerCentroid_X, fElShowerCentroid_Y, fElShowerCentroid_Z);
        }
       

        // ---------------------------------------------------------
        // Calculate L and Q+L energy
        //   E = (N_i + N_ex) * W_ph
        //   N_i + N_ex = N_e + N_ph 
        if( fElShowerVis > 0 ) {
          fElShowerPhotons  = fPE / fElShowerVis;
          fElShowerEnergyQL = ( fElShowerPhotons + fElShowerCharge )* (fWph * 1e-6);
        }
        
      }//<-- endif Npts3D > 0 
    }//<-- endif shower size > 0 on both planes
  }//<-- endif tagged muon + Michel pulse found
   
   
   
    
    
  // *******************************************************************
  // Analysis of event reduction from these quality cuts, and filling histograms
  // *******************************************************************
  
  bool goodShower2D = false;
  bool goodShower3D = false;
  
  hEventCuts->Fill(1-1);
  
  // --------------------------------------------------------------------------- 
  // One stopping track
  if( fNumTrackStopping == 1 ) {
    hEventCuts->Fill(2-1); 
  
    // --------------------------------------------------------------------------- 
    // Optical ID
    if( fMichelOpticalID ) {
      hEventCuts->Fill(3-1);
      
      hShowerSizeCompare    ->Fill(fElShowerSize, fElShowerSize_Pl0);
      hElClusterSizeCompare ->Fill(fElClusterSize, fElClusterSize_Pl0);
      
      // --------------------------------------------------------------------------- 
      // Cluster bnd found
      if( isClusterBndFound ) {
        hEventCuts->Fill(4-1);    
        
        hClusterSize      ->Fill(fClusterSize);
        hMuClusterSize    ->Fill(fMuClusterSize);
        hElClusterSize    ->Fill(fElClusterSize);
        hFracMuHitsLinear ->Fill(fFracMuHitsLinear);
        hMuAveLinearity   ->Fill(fMuAveLinearity);
        hExtraHits        ->Fill(fExtraHits);
    
     
        // --------------------------------------------------------------------------- 
        // Cluster size cuts
        if(   fClusterSize >= fMinClusterSize
          &&  fMuClusterSize >= fMinMuClusterSize
          &&  fElClusterSize >= fMinElClusterSize
          &&  fElClusterSize <= fMaxElClusterSize  ){
          hEventCuts->Fill(5-1);
          
          // -------------------------------------------------------------------
          // extra hits
          if( fExtraHits <= fMaxExtraHits ){
            hEventCuts->Fill(6-1);

            // -----------------------------------------------------------------
            // muon linearity cuts
            if(     fFracMuHitsLinear > fMinFracMuHitsLinear
                &&  fMuAveLinearity > fMinMuLinearity ) {
        

              // ---------------------------------------------------------------
              // muon direction
              if( fMuClusterHitsEndFit >= fMinMuClusterHitsEndFit ) {
                hEventCuts->Fill(7-1);
                
                if( fElShowerFrac > 0 ) {
                  hElShowerSize->Fill( fElShowerSize );
                  hElShowerFrac->Fill( fElShowerFrac );
                  hDecayAngle2D->Fill( fDecayAngle2D ); 
                }
                
                // -------------------------------------------------------------
                // decay angle cut
                if( fDecayAngle2D > fMinDecayAngle2D && fDecayAngle2D < fMaxDecayAngle2D ) {

                  // -----------------------------------------------------------
                  // Shower frac cut
                  if( fElShowerFrac > fMinElShowerFrac && fElShowerEnergy > 0) {
                    hEventCuts->Fill(8-1);
                    
                    goodShower2D = true;
                    fNumEvents_Shwr++;
                    MakeClusteringGraphs(cl_pl1);

                    // ------------------------------
                    // 3D shower stuff
                    
                    // "good" shower on induction plane
                    if( fElShowerSize_Pl0 > 0 &&
                        fMuClusterSize_Pl0 >= fMinMuClusterSize &&
                        fElClusterSize_Pl0 >= fMinElClusterSize
                        //fMuClusterHitsEndFit_Pl0 >= fMinMuClusterHitsEndFit &&
                        //fMuAveLinearity_Pl0 > fMinMuLinearity &&
                        //fFracMuHitsLinear_Pl0 > fMinFracMuHitsLinear &&
                        //fDecayAngle2D_Pl0 > fMinDecayAngle2D && 
                        //fDecayAngle2D_Pl0 < fMaxDecayAngle2D 
                      ) {
                    
                      hNumPts3D   ->Fill(fNumPts3D);
                      hFracHits3D  -> Fill(fFracHits3D);  

                      // N pts 3D > min
                      if( fNumPts3D > fMinNumPts3D ) {
                        
                        // Frac hits
                        if( fFracHits3D > fMinFracHits3D ) {
                          
                          goodShower3D = true;        
                          fNumEvents_Shwr3D++;

                        }//frac hits

                      }// Npts3D
                      
                    }// good shower on induction 

                  }// shower frac
                } //decay angle   
              }// mu dir
            }//mu lin
          }//extra hits
        }// cluster size
      }//cluster bnd
    }// stop trk
  }// Michel optical ID



  // *******************************************************************
  // Make a few histograms (the rest will be made in a separate macro).
  // *******************************************************************
  for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
    size_t ch   = fPulses[ipmt].OpChannel(); 
    if( fDecayTime > fdTcut && goodShower2D ) {
          hPE_prompt_dTcut_shwr[ch]       ->Fill(fPE_prompt[ch]);
          hPE_total_dTcut_shwr[ch]        ->Fill(fPE_total[ch]);
          hPE_totalRaw_dTcut_shwr[ch]     ->Fill(fPE_totalRaw[ch]);
          hTrue_PE_total_dTcut_shwr[ch]   ->Fill(fTrue_PE_total[ch]);
          hTrue_PE_total_TrueVsReco[ch]   ->Fill(fTrue_PE_total[ch], fPE_total[ch]);
          if( fTrue_PE_total[ch] > 0 ) hTrue_PERes[ch]->Fill( (fPE_total[ch] - fTrue_PE_total[ch] ) / fTrue_PE_total[ch] );
    }
  }//<-- endloop PMTs
      
  if( fMuTrackIsCaloOrdered ) {
    for( size_t i=0; i<fvMuTrackdEdx.size(); i++){
      hMudEdx->Fill(fvMuTrackdEdx[i]);
      hMuResRangeVsdEdx->Fill(fvMuTrackResidualRange[i], fvMuTrackdEdx[i]);
    }
  }

  if( fIsMC ) {
    if( fMuonID >= 0 ) {
      hTrue_MuTrackEnd_ZX->Fill( fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_X);
      hTrue_MuTrackEnd_ZY->Fill( fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_Y);
    }
  }

  if( goodShower2D ) {
    hElShowerAngleMean  ->Fill( fElShowerAngleMean );
    hElShowerAngleRMS   ->Fill( fElShowerAngleRMS );
    hElShowerEnergy     ->Fill(fElShowerEnergy);
    hElShowerCharge     ->Fill(fElShowerCharge);
    
    if( goodShower3D ) {
      hElShowerPhotons    ->Fill(fElShowerPhotons);
      hElShowerEnergyQL   ->Fill(fElShowerEnergyQL);
      hElTrackdEdx        ->Fill(fElTrackdEdx);
      float QoverL = fElShowerCharge / fElShowerPhotons;
      float R = ( 1.+fExcRatio ) * QoverL / ( 1.+QoverL );
      hQoverL   ->Fill( QoverL );
      hRecomb   ->Fill( R ); 
    }

    if( fTrue_ElShowerEnergyDep > 0. ) {
      hTrue_ElEnergy                ->Fill(fTrue_ElEnergy);
      if( fTrue_MuCharge == 1 ) hTrue_ElEnergyFree->Fill(fTrue_ElEnergy);
      if( fTrue_MuCharge == -1 ) hTrue_ElEnergyCapture->Fill(fTrue_ElEnergy);
      hTrue_ElShowerChargeDep       ->Fill(fTrue_ElShowerChargeDep);
      hTrue_ElTrackChargeDep        ->Fill(fTrue_ElTrackChargeDep);
      hTrue_ElShowerEnergyDep       ->Fill(fTrue_ElShowerEnergyDep);
      hTrue_ElTrackEnergyDep        ->Fill(fTrue_ElTrackEnergyDep);
      hTrue_ElShowerPhotons         ->Fill(fTrue_ElShowerPhotons);
      hTrue_EnergyDepVsRecoEnergy   ->Fill(fTrue_ElShowerEnergyDep, fElShowerEnergy);
      hTrue_EnergyRes        ->Fill( (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
      hTrue_EnergyVsEnergyRes       ->Fill( fTrue_ElShowerEnergyDep, (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
      hTrue_EnergyResTrk     ->Fill( (fElTrackEnergy - fTrue_ElTrackEnergyDep) / fTrue_ElTrackEnergyDep );
      if( fElShowerEnergyQL > 0 && goodShower3D ) {
        hTrue_EnergyResQL ->Fill((fElShowerEnergyQL - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep);
        hTrue_EnergyRes_Shwr3D       ->Fill( (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
        hTrue_PhotonRes        ->Fill( (fElShowerPhotons - fTrue_ElShowerPhotons) / fTrue_ElShowerPhotons );
        hTrue_VisRes            ->Fill( (fElShowerVis - fTrue_ElShowerVis) / fTrue_ElShowerVis );
      }
    }

  }// end goodShower2D
  // Done making histograms 
  


  // *******************************************************
  // Fill the tree  
  // *******************************************************
  if( isCalibrationEvent || (fMuTrackIndex >= 0 && fMichelOpticalID) ) fTree->Fill();

  // Make waveform graphs
  if( fElShowerEnergy > 0 && fMichelOpticalID ) MakeWfmGraphs(e);

  LOG_VERBATIM("MichelAnaFilter")
  <<"Event summary: \n"
  <<"  isOneStoppingTrack?        "<<(int)isOneStoppingTrack<<"\n"
  <<"  OpticalMichelID?           "<<(int)fMichelOpticalID<<"\n"
  <<"  2D shower energy           "<<fElShowerEnergy<<"\n"
  <<"  3D points                  "<<fNumPts3D;


 
  // *******************************************************
  // Filter functionality
  // *******************************************************
  if( fFilter_PassAllEvents ) {
    return true;
  } else {
    return (
      fMichelOpticalID 
      && fDecayTime         >= fFilter_DecayTime[0]
      && fDecayTime         <= fFilter_DecayTime[1]
      && fNumTrackStopping  >= fFilter_NumTrackStopping[0] 
      && fNumTrackStopping  <= fFilter_NumTrackStopping[1] 
      && fElShowerEnergy    >= fFilter_ElShowerEnergy[0]
      && fElShowerEnergy    <= fFilter_ElShowerEnergy[1]
      );
  }
}



//########################################################################################
void MichelAnaFilter::endJob()
{
  TF1 fit("fit","[0]*exp(-x/[1])",300,7000);
  if( fNumEvents_Shwr == 0 ) {
    hTrue_EnergyRes->Fill(-9);
    hElShowerEnergy->Fill(-9); 
  }

  std::cout
  <<"============================================================\n"
  <<"            Ending MichelAnaFilter Job\n"
  <<"\n";

  std::cout
  <<"  ---------------------------------------------------------  \n"
  <<"                   Settings Used\n"
  <<"\n"
  <<"  GateDelay                              : "<<fGateDelay<<"\n"
  <<"  Max dT difference                      : "<<fMaxDiff_dT<<" ns\n";
  for(size_t i=0; i<fSelectChannels.size(); i++){
  size_t ch=fSelectChannels[i];
  std::cout
  <<"  OpDet"<<ch<<"\n"
  <<"    * Max RMS                            : "<<fMaxWfmRMS[ch]<<" ADC\n"
  <<"    * Gradient hit threshold             : "<<fGradHitThresh[ch]<<"\n"
  <<"    * Min pulse width                    : "<<fMinOpHitWidth[ch]<<" ADC\n"
  <<"    * Prompt PE cut                      : "<<fPromptPECut[ch]<<" PE\n";
  }
  std::cout
  <<"\n";

  std::cout
  <<"  MaxHitSeparation                       : "<<fMaxHitSeparation<<"\n"
  <<"  TruncMeanWindow                        : "<<fTruncMeanWindow<<"\n"
  <<"  TruncMeanP                             : "<<fTruncMeanP<<"\n"
  <<"  MaxMatchDistance                       : "<<fMaxMatchDistance<<"\n"
  <<"  LinThresh                       : "<<fLinThresh<<"\n"
  <<"  LinThreshMuFit                  : "<<fLinThreshMuFit<<"\n"
  <<"  MinClusterSize                         : "<<fMinClusterSize<<"\n"
  <<"  MinMuClusterSize                       : "<<fMinMuClusterSize<<"\n"
  <<"  MinElClusterSize                       : "<<fMinElClusterSize<<"\n"
  <<"  MinElShowerFrac                        : "<<fMinElShowerFrac<<"\n"
  <<"  MinFracMuHitsLinear                    : "<<fMinFracMuHitsLinear<<"\n"
  <<"\n"
  <<"  RFactor                                 : "<<fRecombFactor<<"\n"
  <<"  RFactorTrack                            : "<<fRecombFactorTrack<<"\n"
  <<"  RFactorShower                           : "<<fRecombFactorShower<<"\n"
  <<"\n";
  
  
  if( fIsMC ) {
  std::cout
  <<"  ---------------------------------------------------------  \n"
  <<"                    MC Truth Info                      \n"
  <<"\n"
  <<"  Ave energy deposited (true)            : "<<hTrue_ElShowerEnergyDep->GetMean()<<" MeV\n"
  <<"  Mean energy resolution                 : "<<hTrue_EnergyRes->GetMean()<<" (RMS = "<<hTrue_EnergyRes->GetRMS()<<")\n"
  <<"\n";
  }
  
  
  std::cout
  <<"  ---------------------------------------------------  \n"
  <<"                    Reco Summary                         \n"
  <<"\n"
  <<"  Total events                           : "<<fNumEvents<<"\n"
  <<"  Events w/ 1 stopping track             : "<<fNumEvents_StpTrk<<"\n"
  <<"  Events w/ optical Michel ID            : "<<fNumOpticalID<<"\n"
  <<"  Reconstructed Michel showers           : "<<fNumEvents_Shwr<<"\n"
  <<"  Events with showers on both planes     : "<<fNumEvents_Shwr3D<<"\n"
  <<"  Ave energy reconstructed               : "<<hElShowerEnergy->GetMean()<<" MeV\n" 
  <<"\n";

  // PMT diagnostic info
  for(size_t i=0; i<fSelectChannels.size(); i++){
  size_t ch=fSelectChannels[i];
  if( numMichelWfms[ch] == 0 ) {
      hPE_prompt[ch]->Fill(-9);
      hPE_total[ch]->Fill(-9);
      hPE_total_dTcut_shwr[ch]->Fill(-9);
      hPE_totalRaw_dTcut_shwr[ch]->Fill(-9);
  }
  std::cout
  <<"  OpDet"<<ch<<"\n"
  <<"    - Ave waveform RMS                   : "<<hWfmRMS[ch]->GetMean(1)<<" ADC\n"
  <<"    - Ave number of hits                 : "<<hNumOpHits[ch]->GetMean(1)<<" (RMS = "<<hNumOpHits[ch]->GetRMS(1)<<")\n"
  <<"    - Total Michel pulses ID'd           : "<<numMichelWfms[ch]<<"\n"
  <<"    - Muon pulse sat. rate               : "<<(float)numMuSaturated[ch]/numMichelWfms[ch]<<"\n"
  <<"    - Michel pulse sat. rate             : "<<(float)numElSaturated[ch]/numMichelWfms[ch]<<"\n"
  <<"    - Prompt frac (raw)                  : "<<hPromptFracRaw[ch]->GetMean()<<"\n"
  <<"    - Prompt frac (corrected)            : "<<hPromptFrac[ch]->GetMean()<<"\n"
  <<"    - Prompt frac (raw, dTcut)           : "<<hPromptFracRaw_dTcut[ch]->GetMean()<<"\n"
  <<"    - Prompt frac (corrected, dTcut)     : "<<hPromptFrac_dTcut[ch]->GetMean()<<"\n"
  <<"    - Prompt PE mean                     : "<<hPE_prompt[ch]->GetMean()<<"\n"
  <<"    - Raw total PE mean (dT>3us, shwr)   : "<<hPE_totalRaw_dTcut_shwr[ch]->GetMean()<<"\n"
  <<"    - Cor total PE mean (dT>3us, shwr)   : "<<hPE_total_dTcut_shwr[ch]->GetMean()<<"\n";
  if( fIsMC ) {
  std::cout
  <<"    - True total PE mean (dT>3us, shwr)  : "<<hTrue_PE_total_dTcut_shwr[ch]->GetMean()<<"\n";
  }
  std::cout
  <<"\n";
  
  std::cout
  <<"    - Num entries in ave muon wfm        : "<<aveWfmCounter_mu[ch]<<"\n"
  <<"    - Num entries in ave electron wfm    : "<<aveWfmCounter_el[ch]
  <<"\n\n";
  
  float t1 = 300;
  float t2 = 4000;
  std::cout<<"  Extracting tau from ave wfms, T=("<<t1<<","<<t2<<")ns\n";
  if( aveWfmCounter_mu[ch] > 0 ) {
    hAveWfm_mu[ch]->Scale(1./float(aveWfmCounter_mu[ch]));
    hAveWfm_mu[ch]->SetOption("hist");
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm_mu[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAveWfm_mu[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm_mu[ch]->Fit("fit","RNQ0");
    std::cout
    <<"    Mu wfm tau               :  "<<fit.GetParameter(1)<<" +/- "<<fit.GetParError(1)<<" ns\n";
  }

  if( avePhelProfileCounter_mu[ch] > 0 ) {
    hAvePhelProfile_mu[ch]->Scale(1./float(avePhelProfileCounter_mu[ch]));
    fit.SetRange(300,1800);
    fit.SetParameter(0, hAvePhelProfile_mu[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAvePhelProfile_mu[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAvePhelProfile_mu[ch]->Fit("fit","RNQ0");
  }

  if( aveWfmCounter_el[ch] > 0 ) {
    hAveWfm_el[ch]->Scale(1./float(aveWfmCounter_el[ch]));
    hAveWfm_el[ch]->SetOption("hist");

    // For the electron waveform, the pre-pulse baseline will likely
    // be already exponentially decreasing due to the late-light from 
    // the muon pulses. So fit/subtract this out first.
    fit.SetRange(hAveWfm_el[ch]->GetXaxis()->GetXmin(), -10);
    fit.SetParameter( 0, 0.1);
    fit.SetParLimits(0, 0, hAveWfm_el[ch]->GetMaximum() );
    fit.SetParameter( 1, 1200 );
    fit.SetParLimits(1, 500, 2000 );
    hAveWfm_el[ch]->Fit("fit","RNQ0");
    for(int i=1; i<=hAveWfm_el[ch]->GetNbinsX(); i++){
      hAveWfm_el[ch]->SetBinContent(i, hAveWfm_el[ch]->GetBinContent(i) - fit.Eval(hAveWfm_el[ch]->GetBinCenter(i)));
    }
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm_el[ch]->GetMaximum()/10.);
    fit.SetParLimits(0, 0, hAveWfm_el[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm_el[ch]->Fit("fit","RNQ0");
    std::cout
    <<"    Electron wfm tau         :  "<<fit.GetParameter(1)<<" +/- "<<fit.GetParError(1)<<" ns\n";
  }
  std::cout<<"\n";

  }
 
  // -------------------------------------------
  // Estimate e lifetime
  if( hdQdxVsT->GetEntries() > 0 ) {
  hdQdxVsT->Divide(hdQdxVsT_N);
  fit.SetRange(10.,310.);
  fit.SetParameter(0, hdQdxVsT->GetMaximum() );
  fit.SetParLimits(0, 0., hdQdxVsT->GetMaximum() * 3. );
  fit.SetParameter(1, 800.);
  fit.SetParLimits(1, 50., 3000.);
  hdQdxVsT->Fit("fit","RQN0");
  hElectronLifetimeReco->Fill(fit.GetParameter(1) );
  std::cout
  <<"  Reco'd electron lifetime               : "<<fit.GetParameter(1)<<" +/- "<<fit.GetParError(1)<<" us\n\n";
  }
     
  std::cout
  <<"  Saved PMT waveforms                    : "<<fNumSavedWfmGraphs<<" / "<<fMaxSavedWfmGraphs<<"\n"
  <<"  Saved hit-clustering graphs            : "<<fNumSavedHitGraphs<<" / "<<fMaxSavedHitGraphs<<"\n"
  <<"\n";  
  
  std::cout
  <<"=======================================================\n";

}



//########################################################################################
// This function performs an "idealized" PMT waveform reconstruction over truth-level
// optical information, since we don't have a full PMT electronics response simulated.
void MichelAnaFilter::PerfectWfmReco(raw::OpDetPulse &pulse) 
{
  size_t ch   = pulse.OpChannel(); 

  fWfmRMS[ch] = 0.;

  if( fMuT0[ch] > 0.  ) fvHitTimes[ch].push_back( fMuT0[ch] );
  if( fT0[ch] > 0.    ) fvHitTimes[ch].push_back( fT0[ch]   );
  fNumOpHits[ch] = fvHitTimes[ch].size();
  
  for(int i=0; i<fNumOpHits[ch]; i++){
        
        float T = fvHitTimes[ch].at(i) - 1;
        float PE_100 = Integrate(hPMT_phelTimes[ch],T,T+100);
        float PE_200 = Integrate(hPMT_phelTimes[ch],T,T+200);
        float PE_300 = Integrate(hPMT_phelTimes[ch],T,T+300);
        float PE_400 = Integrate(hPMT_phelTimes[ch],T,T+400);
        float PE_500 = Integrate(hPMT_phelTimes[ch],T,T+500);
        float PE_600 = Integrate(hPMT_phelTimes[ch],T,T+600);
        float PE_700 = Integrate(hPMT_phelTimes[ch],T,T+700);
        float PE_900 = Integrate(hPMT_phelTimes[ch],T,T+900);
        float PE_1200 = Integrate(hPMT_phelTimes[ch],T,T+1200);
        float PE_1500 = Integrate(hPMT_phelTimes[ch],T,T+1500);
        float PE_1800 = Integrate(hPMT_phelTimes[ch],T,T+1800);
        float PE_7000 = Integrate(hPMT_phelTimes[ch],T,T+7000);

        if( i==1) fvIsHitAtTrigger[ch]    .push_back( true );
        else      fvIsHitAtTrigger[ch]    .push_back( false );
        fvIsHitSaturated[ch]    .push_back( false );

          float tmp = 0;
          tmp = PE_100;
          float pe_100 =          tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_200-PE_100;  
          float pe_200 = pe_100 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_300-PE_200;  
          float pe_300 = pe_200 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_400-PE_300;  
          float pe_400 = pe_300 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_500-PE_400;  
          float pe_500 = pe_400 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_600-PE_500;  
          float pe_600 = pe_500 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_700-PE_600;  
          float pe_700 = pe_600 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 100. ));
          tmp = PE_900-PE_700;  
          float pe_900 = pe_700 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 200. ));
          tmp = PE_1200-PE_900;  
          float pe_1200 = pe_900 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 300. ));
          tmp = PE_1500-PE_1200;  
          float pe_1500 = pe_1200 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 300. ));
          tmp = PE_1800-PE_1500;  
          float pe_1800 = pe_1500 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 300. ));
          tmp = PE_7000-PE_1800;  
          float pe_7000 = pe_1800 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, 5200. ));

          fvHitADC_100ns[ch]      .push_back( pe_100    *fSPE[ch] );
          fvHitADC_200ns[ch]      .push_back( pe_200    *fSPE[ch] );
          fvHitADC_300ns[ch]      .push_back( pe_300    *fSPE[ch] );
          fvHitADC_400ns[ch]      .push_back( pe_400    *fSPE[ch] );
          fvHitADC_500ns[ch]      .push_back( pe_500    *fSPE[ch] );
          fvHitADC_600ns[ch]      .push_back( pe_600    *fSPE[ch] );
          fvHitADC_700ns[ch]      .push_back( pe_700    *fSPE[ch] );
          fvHitADC_900ns[ch]      .push_back( pe_900    *fSPE[ch] );
          fvHitADC_1200ns[ch]     .push_back( pe_1200    *fSPE[ch] );
          fvHitADC_1500ns[ch]     .push_back( pe_1500    *fSPE[ch] );
          fvHitADC_1800ns[ch]     .push_back( pe_1800   *fSPE[ch] );
          fvHitADC_7000ns[ch]     .push_back( pe_7000   *fSPE[ch] );
        
        fvPrepulseBaseline[ch]  .push_back(0);
        fvPrepulseRMS[ch]       .push_back(0);
        fvPrepulseSlowNorm[ch]  .push_back(0);
        fvPrepulseSlowTau[ch]   .push_back(0); 
        fvPrepulseX1[ch]        .push_back(0);  
        fvPrepulseZeroPoint[ch] .push_back(0);
        
        // The prefit integrals are only really sensical in data;
        // in MC, we will just assign them as truth for now
        if( i==0 ) {
          fvHitADCpf_100ns[ch]    .push_back( fTrue_MuPE_prompt[ch] * fSPE[ch] );
          fvHitADCpf_7000ns[ch]   .push_back( PE_7000 * fSPE[ch] );
        } else
        if( i==1 ) {
          fvHitADCpf_100ns[ch]    .push_back( fTrue_PE_prompt[ch] * fSPE[ch] );
          fvHitADCpf_7000ns[ch]   .push_back( fTrue_PE_total[ch] * fSPE[ch]);
        }
          
        fvHitWidth[ch]          .push_back(fMinOpHitWidth[ch]);
        fvHitAmplitudes[ch]     .push_back(0.); // convert PEs to amplitude
        
  }

}



//########################################################################################
// Given a PMT channel, interated photoelectron count, and integration time, calculate the 
// smearing sigma based on noise information provided in the fcl (extracted from data).
float MichelAnaFilter::CalcSigma(size_t ch, float PE, float T ) {
  float sigPE = 0.;
  float sigT  = 0.;
  if( fSmearSigPE[ch] > 0 ) sigPE = fSmearSigPE[ch] *sqrt(PE);
  if( fSmearSigT[ch]  > 0 ) sigT  = fSmearSigT[ch]  *sqrt(T);
  return sqrt( pow(sigPE,2) + pow(sigT,2) );
}



//########################################################################################
// Function for determining if a point [cm] is inside or outside predefined fiducial volume 
// (fiducial margins [cm] = fx, fy, fz)
/*
bool MichelAnaFilter::IsPointInFiducialVolume(
  TVector3 p, 
  float dx1 = 0., float dx2 = 0.,
  float dy1 = 0., float dy2 = 0.,
  float dz1 = 0., float dz2 = 0. )
{
  if(    (p.X() >= 0. + dx1) && (p.X() <= 2.*fGeo->DetHalfWidth() - dx2)
      && (p.Y() >= -1.*fGeo->DetHalfHeight() + dy1) && (p.Y() <= fGeo->DetHalfHeight() - dy2)
      && (p.Z() >= 0. + dz1) && (p.Z() <= fGeo->DetLength() - dz2) )
  {
    return true;
  } else {
    return false;
  }
}
*/

//########################################################################################
// Function for determining if a point [cm] is inside or outside predefined fiducial volume 
// (fiducial margins [cm] = fx, fy, fz)
bool MichelAnaFilter::IsPointInFiducialVolume(TVector3 p, float dx = 0., float dy = 0., float dz = 0.)
{
  if(    (p.X() >= 0. + dx) && (p.X() <= 2.*fGeo->DetHalfWidth() - dx)
      && (p.Y() >= -1.*fGeo->DetHalfHeight() + dy) && (p.Y() <= fGeo->DetHalfHeight() - dy)
      && (p.Z() >= 0. + dz) && (p.Z() <= fGeo->DetLength() - dz) )
  {
    return true;
  } else {
    return false;
  }
}



//########################################################################################
// Retrieve all detector properties (TPC and PMTs)
void MichelAnaFilter::GetDetProperties(){

  if( fRunNumber != fCachedRunNumber ) {
    fCachedRunNumber = fRunNumber;
    
    // ================================================   
    // Update detector properties
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fEfield             = detprop->Efield(0);           // kV/cm
    fElectronLifetime   = detprop->ElectronLifetime();  // microseconds
    fDriftVelocity      = detprop->DriftVelocity();
    fXTicksOffset[0]    = detprop->GetXTicksOffset(0,0,0);
    fXTicksOffset[1]    = detprop->GetXTicksOffset(1,0,0);
    fSamplingRate       = detprop->SamplingRate()*1e-3;

    // ================================================   
    // Get SPE values for PMTs
    if( fIsMC ) {
      fSPE[0]     = 50.;
      fSPE_err[0] = 25.;
      fSPE[1]     = 50.;
      fSPE_err[1] = 25.;
    } else {
    cet::search_path sp("FW_SEARCH_PATH");
    for(size_t i = 0; i < 2; i++){
      sprintf(buffer,"LArIATPhotodetectorSER_ch%lu_Michels.txt",i);
      std::string fullname;
      sp.find_file(buffer, fullname);
      if (fullname.empty()) {
        throw cet::exception("File not found");
      } else {
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        size_t lineNum = 0; int run1 = 0.; int run2 = 0.;
        while (std::getline(inFile,line)) {
          // skip first line (column headers)
          lineNum++;
          if(lineNum == 1 ) continue;
          float spe, speErr, width, tau;
          std::istringstream ss(line);
          ss >> run1 >> run2 >> spe >> speErr >> width >> tau;
          if( fRunNumber >= run1 && fRunNumber <= run2 ){
            fSPE[i]	          = spe*fSPE_ScaleFactor[i];
            fSPE_err[i]           = speErr*fSPE_ScaleFactor[i];
            if( tau > 0 ) {
              fMuContamCorr_EffTau[i]   = tau;            
            }
            break;
          }
        }
      }
    } 
         
    }// end getting SPE
  
  }// endif new run

}



//########################################################################################
// This function should 
//   - fill all member truth variables
//   - create artificial OpDetPulses (fPulses) to be 
//     passed on to later stages for quasi-waveform reco
void MichelAnaFilter::GetTruthInfo(art::Event &e){

  LOG_VERBATIM("MichelAnaFilter")<<"GetTruthInfo";
  
  // Get MCParticle information, and exit if none is found
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  e.getByLabel(fSimProducerModule,particleHandle);
  
  // Reset electron/muon ID
  fElectronID = -9;
  fMuonID     = -9;

  LOG_VERBATIM("MichelAnaFilter")<<"This event has "<<particleHandle->size()<<" MCParticles";

  // ********************************************************
  // Loop through MC particles
  // *******************************************************
  for( auto const& particle : (*particleHandle) )
  {
    // Save useful info
    int	  last	      = particle.NumberTrajectoryPoints() - 1;
    size_t nTrajPts   = particle.NumberTrajectoryPoints();
    TVector3 locStart = particle.Position(0).Vect();
    TVector3 locEnd   = particle.Position(last).Vect();
    TVector3 momStart = particle.Momentum(0).Vect();
    TVector3 momEnd   = particle.Momentum(last).Vect();
    bool startsInVol  = IsPointInFiducialVolume(locStart,0.,0.,0.);
    bool endsInVol    = IsPointInFiducialVolume(locEnd,0.,0.,0.);
    double dL         = (locStart-locEnd).Mag();
    double P0         = particle.P(0)*1000.;
    double E0         = 1000.*( particle.E(0) - particle.Mass() );
    double Ef         = 1000.*( particle.E(last) - particle.Mass() );
    double dE         = (E0-Ef);
    bool isFromEl     = IsParticleDescendedFrom(particleHandle, particle.TrackId(), fElectronID);
   
    // Printout particle list
    if(0){
    printf("  %3i PDG: %5i (%i->%i) Npts=%3lu dL=%6.1f P0=%10.6f Etot=%10.6f E=%10.6f dE=%10.6f  T=%11.6f-%11.6f   moth=%3i   %30s    fromE? %i  Ndaughters=%i\n",
      particle.TrackId(),
      particle.PdgCode(),
      (int)startsInVol,
      (int)endsInVol,
      nTrajPts,
      dL,
      P0,
      particle.E(0)*1000.,
      E0,
      dE,
      particle.T(0),
      particle.T(last),
      particle.Mother(),
      particle.Process().c_str(),
      isFromEl,
      particle.NumberDaughters()
      );
    }
    
    // =======================================================
    // Look for muon
    if( particle.Process() == "primary" && abs(particle.PdgCode()) == 13 ) {

      // save muon track ID, charge, and stop time	
      fMuonID           = particle.TrackId();
      fTrue_MuCharge    = -1.*(particle.PdgCode()/abs(particle.PdgCode()));
      fTrue_MuStopTime  = particle.T(last);
      
      // save muon endpoint
      fTrue_MuTrackEnd        = particle.Position(last).Vect();
      fTrue_MuTrackEnd_X      = fTrue_MuTrackEnd.X();
      fTrue_MuTrackEnd_Y      = fTrue_MuTrackEnd.Y();
      fTrue_MuTrackEnd_Z      = fTrue_MuTrackEnd.Z();

      // find where muon first enters TPC
      for(size_t i=0; i<particle.NumberTrajectoryPoints(); i++) {
        if( IsPointInFiducialVolume( particle.Position(i).Vect(), 0, 0, 0) ) {
          fTrue_MuTrackVertex   = particle.Position(i).Vect();
          fTrue_MuTrackVertex_X = fTrue_MuTrackVertex.X();
          fTrue_MuTrackVertex_Y = fTrue_MuTrackVertex.Y();
          fTrue_MuTrackVertex_Z = fTrue_MuTrackVertex.Z();
          break;
        }
      }

      LOG_VERBATIM("MichelAnaFilter") 
      << "  --> Found primary muon!  Charge = "<<fTrue_MuCharge
      <<", ("<<fTrue_MuTrackVertex_X<<","<<fTrue_MuTrackVertex_Y<<","<<fTrue_MuTrackVertex_Z<<") "
      <<" --> ("<<fTrue_MuTrackEnd_X<<","<<fTrue_MuTrackEnd_Y<<","<<fTrue_MuTrackEnd_Z<<")"
      <<"  L = "<< (fTrue_MuTrackVertex - fTrue_MuTrackEnd).Mag();
        
    }//<-- endif muon


    // =======================================================
    // Look for Michel electron
    if( abs(particle.PdgCode()) == 11 && particle.Mother() == fMuonID ) {
      
      // since the Geant4 "Process" param is useless in case of mu-, 
      // make some clever cuts to ensure we really have a Michel:
      if( -1.*particle.PdgCode()/abs(particle.PdgCode()) == fTrue_MuCharge && (particle.T(0) - fTrue_MuStopTime) > 5.0 ){
        
        // save electron track ID, decay time, energy, momentum, and angle
        fElectronID     = particle.TrackId();
        fTrue_dT        = particle.T(0) - fTrue_MuStopTime;
        fTrue_ElEnergy  = particle.E(0)*1000.;
        fTrue_ElMomentum = particle.Momentum(0).Vect();
        fTrue_ElMomentum_X  = fTrue_ElMomentum.X();
        fTrue_ElMomentum_Y  = fTrue_ElMomentum.Y();
        fTrue_ElMomentum_Z  = fTrue_ElMomentum.Z();
        fTrue_ElAngle   = fTrue_ElMomentum.Angle(fTrue_MuTrackEnd - fTrue_MuTrackVertex);
          
        // is the electron fully contained?
        fTrue_IsElContained = true;
        for(size_t i=0; i<particle.NumberTrajectoryPoints(); i++) {
          if( !IsPointInFiducialVolume(particle.Position(i).Vect(), 0., 0., 0.) ) {
            fTrue_IsElContained = false;
            break;
          }
        }
        
        LOG_VERBATIM("MichelAnaFilter")
        << "  --> Found Michel electron!  dT = "<<fTrue_dT<<", KE = "<<fTrue_ElEnergy<<", angle (rel. to muon) = "<<fTrue_ElAngle*RAD_TO_DEG<<" deg";

      }//<-- endif Michel electron
    }//<-- endif electron
 

  }//<-- end loop over MCParticles
  
  
  // ******************************************************************
  // If the michel electron was found, then simulate the scintillation signals
  // ******************************************************************
  if( fMuonID > 0 && fElectronID > 0 ) {

    // This stage loops through the IDEs from SimChannels and 
    // calculates the scintillation light produced
    ParticleTracker(e, fMuonID, fElectronID, hPMT_phelTimes, hPMT_phelTimes_electron); 

    // Loop over PMTs
    for( auto & i : fSelectChannels ) {

      // Approximate the "hit time" T0 as the time of arrival of the first photon
      fMuT0[i]                = (float)hPMT_phelTimes[i]->FindFirstBinAbove(0,1);
      fT0[i]                  = (float)hPMT_phelTimes_electron[i]->FindFirstBinAbove(0,1);
      
      // Save prompt, total PEs by integrating the photon arrival 
      // time histograms produced in "ParticleTracker"
      float muX1 = fMuT0[i]-1;
      float elX1 = fT0[i]-1;
      fTrue_PE_prompt[i]      = Integrate(hPMT_phelTimes_electron[i],elX1,fT0[i]+fPromptWindow);
      fTrue_PE_total[i]       = Integrate(hPMT_phelTimes_electron[i],elX1,fT0[i]+fFullWindow);
      fTrue_MuPE_prompt[i]    = Integrate(hPMT_phelTimes[i],         muX1,fMuT0[i]+fPromptWindow);
      fTrue_MuContam_prompt[i]= Integrate(hPMT_phelTimes[i],elX1,fT0[i]+fPromptWindow) - fTrue_PE_prompt[i];
      fTrue_MuContam_total[i] = Integrate(hPMT_phelTimes[i],elX1,fT0[i]+fFullWindow) - fTrue_PE_total[i];
      hTrue_LightYield[i]     ->Fill( hPMT_phelTimes[i]->GetEntries() / fTrue_TotalEnergyDep );
      hTrue_PE_preTrig[i]     ->Fill( fTrue_PE[i] );
    }

    LOG_VERBATIM("MichelAnaFilter")
    <<"  Total charge deposited in event: "<<fTrue_TotalChargeDep<<" electrons, total energy = "<<fTrue_TotalEnergyDep<<" MeV\n"
    <<"  Total energy deposited by electron: "<<fTrue_ElShowerEnergyDep<<"\n"
    <<"    #pe on HMM (prompt/total) = "<<fTrue_PE_prompt[0]<<" / "<<fTrue_PE_total[0]<<"\n"
    <<"    #pe on ETL (prompt/total) = "<<fTrue_PE_prompt[1]<<" / "<<fTrue_PE_total[1];
  }
  
  
  // ******************************************************************
  // Turn the histograms into waveforms and then create and fill fPulses  
  // *****************************************************************
  for( size_t ii=0; ii<fSelectChannels.size(); ii++){ 
    std::vector<short> wfm_tmp(10,0);
    fPulses.push_back(raw::OpDetPulse ( 
      static_cast <unsigned short> (fSelectChannels.at(ii)), 
      wfm_tmp, 
      static_cast <unsigned int> (0), 
      static_cast <unsigned int> (8600))
      );
  }
  
}



//########################################################################################
// Boolean function which loops through the lineaage of an MCParticle and determines if it
// is descended from a given mother particle.
bool MichelAnaFilter::IsParticleDescendedFrom(const art::Handle< std::vector<simb::MCParticle>>& particleHandle, int particleID, int motherID){
  if( particleID < 0 || motherID < 0 ) return false;
  if( motherID == 0 || motherID == particleID ) return true;
  bool isDescended = false;
  std::vector<int> lineage;
  lineage.reserve(200);
  lineage.push_back(motherID);
  for( auto const& particle : (*particleHandle) )
  {
    // update lineage
    for(size_t i=0; i<lineage.size(); i++){
      if( particle.Mother() == lineage[i] ) {
	lineage.push_back(particle.TrackId());
	if( particle.TrackId() == particleID ) isDescended = true;
      }
    }
    if(particle.TrackId() == particleID) return isDescended;
  }
  return false;
}



//########################################################################################
// Calculate total distance traversed by a particle
float MichelAnaFilter::CalcParticleDist( const simb::MCParticle& p){
  float L = 0; 
  if( p.NumberTrajectoryPoints() > 1 ) {
    for(size_t i=1; i<p.NumberTrajectoryPoints(); i++){
      if( IsPointInFiducialVolume(p.Position(i).Vect())&&(IsPointInFiducialVolume(p.Position(i-1).Vect()))){
        L += (p.Position(i).Vect()-p.Position(i-1).Vect()).Mag();
      }
    }
  }
  return L;
}



//########################################################################################
// This function loops all the MCParticles in the event, calculates energy depositions, and
// fills the photon arrival time histograms for the PMTs. This requires that the particle ID
// of the muon and Michel electron have already been identified. 
void MichelAnaFilter::ParticleTracker(const art::Event& e, int muID, int elID, TH1D* hPhelTimes[], TH1D* hPhelTimes_electron[]){

  LOG_VERBATIM("MichelAnaFilter")<<"Particle tracker...";
 
  // Get IDEs
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > simList;
  if(e.getByLabel(fSimProducerModule,SimListHandle)) { 
    art::fill_ptr_vector(simList, SimListHandle); 
  }
  
  // Get MCParticle information
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  e.getByLabel(fSimProducerModule,particleHandle);

  size_t numDRays       = 0;
  fTrue_NumBremmPhotons = 0;
  fTrue_TotalChargeDep  =0.;
  fTrue_TotalEnergyDep  =0.;
  
  // Keep track of total charge (electrons escaping recomb.) 
  // for Michel electrons and their shower products, and
  // total number of electron-ion pairs. These are used to 
  // calculate recombination factors.
  int   Ni_trk  = 0;
  int   Ni_shwr = 0;
  float Q_trk   = 0.;
  float Q_shwr  = 0.;
  
  if( elID > 0 ) {
    fTrue_ElShowerVis     = 0.; // total frac. vis
    fTrue_ElShowerVisCh[0] = 0.; // total frac. vis: ETL
    fTrue_ElShowerVisCh[1] = 0.; // total frac. vis: HMM
    fTrue_ElTrackPhotons = 0;
    fTrue_ElShowerPhotons   = 0;
    fTrue_ElShowerPhotonsPrompt   = 0;
    fTrue_ElShowerPhotonsLate   = 0;
    fTrue_PE[0] = 0;
    fTrue_PE[1] = 0;
    fTrue_ElTrackChargeDep  = 0.;
    fTrue_ElTrackEnergyDep  = 0.;
    fTrue_ElShowerChargeDep = 0.;
    fTrue_ElShowerEnergyDep = 0.;
  }

  // Loop through particle list, and for each particle with the common 
  // ancestor, determine energy loss deposited INTO LAr at each step.
  for( auto const& particle : (*particleHandle) )
  {
    
    // Is this particle descended from electron, or is it the muon?
    bool isFromElectron = IsParticleDescendedFrom(particleHandle,particle.TrackId(),elID);
    bool isElectron     = (particle.TrackId() == elID);
    bool isMuon         = (particle.TrackId() == muID);
    
    // Bremm photon plots (not sure if this is right...)
    if( isFromElectron && particle.PdgCode() == 22 ) {
      fTrue_NumBremmPhotons++;
      hTrue_BremmPhotonLength->Fill( (particle.EndPosition().Vect() - particle.Position(0).Vect()).Mag() );
    }
    
    // Count delta rays: electrons NOT originating from the Michel electron,
    // which are born prior to the muon's stop time
    if( !isFromElectron && !isMuon && particle.T(0) <= fTrue_MuStopTime && abs(particle.PdgCode()) == 11) numDRays++;

    // Skip photons and unrelated particles
    if( particle.PdgCode() == 22 ) continue;
    int Npts = particle.NumberTrajectoryPoints();
    if( Npts < 1 ) continue;

    // Particle start time
    float particleTime = particle.T(0);
    

    // Particle's deposited charge and energy
    float particleDepEnergy = 0.;
    float particleDepCharge = 0.;
    
    // Loop over track's IDEs and convert each to scintillation 
    for(size_t nChan = 0; nChan < simList.size(); nChan++) {
      if(simList.at(nChan)->Channel() > 240) break;
      // Get the information of each wire
      const auto & wire = simList.at(nChan)->TDCIDEMap();
      // Looping over the IDEs in a wire, or looping over time
      for(auto it = wire.begin(); it != wire.end(); it++) {
        // Looping over the IDEs in a given time tick
        for(size_t i = 0; i < it->second.size(); i++) {
          if(it->second.at(i).trackID != particle.TrackId()) continue;

          // get location and energy deposited
          float dE   = it->second.at(i).energy;
          float x    = it->second.at(i).x;
          float y    = it->second.at(i).y;
          float z    = it->second.at(i).z;
          TVector3   loc(x,y,z);
          hTrue_X   ->Fill(loc.X());
          hTrue_Y   ->Fill(loc.Y());
          hTrue_Z   ->Fill(loc.Z());
          
          // if we're beyond the cathode, skip
          if( loc.X() > fGeo->DetHalfWidth()*2. ) continue;

          // correct numElectrons for drift attenuation
          float T = x / fDriftVelocity;
          int numElectrons = it->second.at(i).numElectrons * TMath::Exp( T / fElectronLifetime );
          hTrue_NumElectrons->Fill(it->second.at(i).numElectrons);
          
          // Use deposited energy to back-calculate scintillation.
          // In NEST, they use W_ph = 19.5 eV, divide up the energy deposited into 
          // excitons/ionization, etc.  Stochastic variation using Resolution Scale 
          // of 0.107 is applied to calculate the total quanta produced.  For nuclear 
          // recoils, there is a further reduction due to the "L" factor as per Lindhard
          // theory ("YieldFactor") but for now we ignore this.
          // http://nusoft.fnal.gov/larsoft/doxsvn/html/NestAlg_8cxx_source.html
          float meanNumQuanta = dE / ( fWph * 1e-6 );
          int totalQuanta = meanNumQuanta;
          //float sigma = sqrt( 0.107 * meanNumQuanta );
          //totalQuanta = (int)fRand->Gaus( meanNumQuanta, sigma );
          if( totalQuanta <= 0 ) continue;
          int numExc      = totalQuanta*(fExcRatio / (1. + fExcRatio));
          int numIon      = totalQuanta - numExc;
          int numPhotons  = totalQuanta - numElectrons;
          if( numPhotons <= 0 ) numPhotons = 0; 

          // add up the energy
          fTrue_TotalChargeDep        += numElectrons; 
          fTrue_TotalEnergyDep        += dE;
          particleDepCharge           += numElectrons;
          particleDepEnergy           += dE;
          
          // Keep track of total charge, light, and quanta for the different components
          // of the Michel shower (electron ionization track and displaced shower products)
          if( isElectron ){
            Ni_trk  += numIon; 
            Q_trk   += numElectrons;
            fTrue_ElTrackPhotons      += numPhotons;
            fTrue_ElTrackEnergyDep    += dE;
          }
          if( isFromElectron || isElectron ){
            fTrue_ElShowerEnergyDep   += dE;
            fTrue_ElShowerPhotons     += numPhotons;
            hTrue_ElShowerDepVsDistFromMuTrackEnd->Fill( (loc-fTrue_MuTrackEnd).Mag() );
          }
          if( !isElectron && isFromElectron ) {
            Ni_shwr += numIon;
            Q_shwr  += numElectrons;
          }
          if( !isFromElectron && !isMuon && particle.T(0) <= fTrue_MuStopTime ) {
            hTrue_DRayEnergy          ->Fill(dE);
          }

          // Divide photons into prompt/delayed
          std::binomial_distribution<int> distPhotonsFast(numPhotons, fFastLightRatio/(1.+fFastLightRatio));
          int numPhotons_fast = distPhotonsFast(generator);
          int numPhotons_late = numPhotons - numPhotons_fast;
          if( isFromElectron || isElectron ) {
            fTrue_ElShowerPhotonsPrompt += numPhotons_fast;
            fTrue_ElShowerPhotonsLate   += numPhotons_late;
          }
          
          // Calculate quenching factor based on given "measured" late-light time const,
          // and quench the prompt and late light accordingly.
          double tauFast = fTauS; // theory value
          double tauLate = fTauT; // theory value
          double quenchFacFast = 1.;
          double quenchFacLate = 1.;
          if( fLateLightTau > 0 && fLateLightTau < tauLate ) {
            double tauQuench = pow( 1./fLateLightTau - 1./tauLate , -1.);
            double tauLateEff = fLateLightTau;
            double tauFastEff = pow( 1./tauQuench + 1/tauFast , -1);
            quenchFacFast = tauFastEff/tauFast;
            quenchFacLate = tauLateEff/tauLate;
            std::binomial_distribution<int> fastQuenched( numPhotons_fast, quenchFacFast );
            numPhotons_fast = fastQuenched(generator);
            std::binomial_distribution<int> lateQuenched( numPhotons_late, quenchFacLate );
            numPhotons_late = lateQuenched(generator);
            tauFast = tauFastEff;
            tauLate = tauLateEff;
          }
          // update num photons after quenching
          numPhotons = numPhotons_fast + numPhotons_late;

          // Loop over the PMTs and for each, run through the list of photons
          // (after quenching) and draw a random variable to determine detection.
          // For now, we will assume each PMT's visibility is low enough such
          // that probabilities are effectively independent, so this is valid.
          int Ndet[2] = {0,0};
          for(size_t ipmt=0; ipmt < fSelectChannels.size(); ipmt++){
            size_t ch = fSelectChannels.at(ipmt);
            
            // Read library information for direct/visible
            float vis[2]    = {0,0};
            float tmin[2]   = {0,0};
            float tmean[2]  = {0,0};
            float trms[2]   = {0.0};
            ReadPhotonLibrary(loc, ch, 0, vis[0], tmin[0], tmean[0], trms[0]);
            ReadPhotonLibrary(loc, ch, 1, vis[1], tmin[1], tmean[1], trms[1]);
            float total_vis_ch = vis[0]+vis[1];
            hTrue_VisTotal[ch]->Fill(total_vis_ch);
            hTrue_Vis[ch][0]  ->Fill(vis[0]);
            hTrue_Vis[ch][1]  ->Fill(vis[1]);
            

            // ------------------------------- 
            // Determine number of phots detected by this PMT
            std::binomial_distribution<int> nDetDistFast( numPhotons_fast, total_vis_ch );
            std::binomial_distribution<int> nDetDistLate( numPhotons_late, total_vis_ch );
            int Ndet_fast = nDetDistFast(generator); 
            int Ndet_late = nDetDistLate(generator); 
            Ndet[ch] = Ndet_fast + Ndet_late;
            if( isFromElectron || isElectron ){ 
              fTrue_PE[ch] += Ndet[ch];
              // weighted vis
              fTrue_ElShowerVisCh[ch] += total_vis_ch * float(numPhotons);
              fTrue_ElShowerVis      +=  total_vis_ch * float(numPhotons); 
            }
            
            // Loop over those photons 
            for(int iph = 0; iph < Ndet[ch]; iph++){
              
              // direct or reflected light?
              int isvis = 0;
              if( fRand->Rndm() > vis[0]/total_vis_ch ) isvis = 1;
              
              // assign a time constant
              float tauLAr = 0.;
              if( iph < Ndet_fast)  tauLAr = tauFast;
              else                  tauLAr = tauLate;
                              
              // ok, so we now have a detected photon and we know everything
              // we need to know about it, so assign arrival times.     
              PropagatePhoton( particleTime,  tauLAr, tmin[isvis], tmean[isvis], trms[isvis], isFromElectron, hPhelTimes[ch], hPhelTimes_electron[ch]);
              
            }// end loop over photons
          }// end loop over PMTs
        }//<-- endloop over IDEs
      }//<-- endloop over time ticks
    }//<-- loop over channels

    // ----------------------------------------------------------------------------------
    // If this particle was Michel electron or one of its shower products, 
    // record total dE/dx and dQ/dx into histograms
    if( (isElectron || isFromElectron) ) {
      
      // traj. length (if fully contained in TPC)
      float particleTrajLength = 0.;
      for( size_t i=0; i<particle.NumberTrajectoryPoints(); i++)
      {
        if( i > 0 ) particleTrajLength += (particle.Position(i).Vect()-particle.Position(i-1).Vect()).Mag();
        if( !IsPointInFiducialVolume( particle.Position(i).Vect(), 0., 0., 0.)){
          particleTrajLength  = 0.;
          break;
        }
      }

      // now calculate dE/dx
      if( particleTrajLength > 0 ) {
        hTrue_TrajLength    ->Fill( particleTrajLength ); 
        
        if( particleTrajLength > 0.1 && particleDepEnergy > 0. ) {
          float dEdx = particleDepEnergy / particleTrajLength;
          float dQdx = particleDepCharge / particleTrajLength;
          hTrue_dEdx          ->Fill( dEdx, particleDepEnergy );
          hTrue_dQdx          ->Fill( dQdx, particleDepCharge );
          hTrue_dEdx_vs_dQdx  ->Fill( dEdx, dQdx);
          if( isElectron ) {
            hTrue_dEdx_ElTrk  ->Fill( dEdx ); 
            hTrue_dQdx_ElTrk  ->Fill( dQdx );
            hTrue_dEdx_vs_dQdx_ElTrk  ->Fill( dEdx, dQdx);
          }
          if( !isElectron && isFromElectron ) {
            hTrue_dEdx_ElShw  ->Fill( dEdx );
            hTrue_dQdx_ElShw  ->Fill( dQdx );
            hTrue_dEdx_vs_dQdx_ElShw  ->Fill( dEdx, dQdx);
          }
        }// endif L > 0.1 && dE > 0
      }// endif > 0
    }// endif electron shower product

  }//<-- end loop over particles
 
  hTrue_NumDRays->Fill(numDRays); 
  hTrue_ElEnergyVsNumBremmPhotons->Fill(fTrue_ElEnergy, fTrue_NumBremmPhotons);
  
  // ----------------------------------------------------------
  // Normalize the weighted true visibilities
  fTrue_ElShowerVisCh[0]  /= fTrue_ElShowerPhotons;
  fTrue_ElShowerVisCh[1]  /= fTrue_ElShowerPhotons;
  fTrue_ElShowerVis       /= fTrue_ElShowerPhotons;

  // -----------------------------------------------------------    
  // Record the survival yield, or, ratio of surviving ionizaiton electrons
  // out of the total number of quanta produced by the energy deposit. The mean
  // from these histograms will be the scale factors we use for charge deposits
  // in data when calculating Q-based energy. (NOTE: this is for Michel electrons
  // and their shower products only.)
  fTrue_ElTrackChargeDep    = Q_trk;
  fTrue_ElShowerChargeDep   = Q_shwr;
  fTrue_ElShowerChargeDep   = Q_trk + Q_shwr;
  hTrue_RecombFactor        ->Fill( (Q_trk + Q_shwr) / float(Ni_trk + Ni_shwr) );
  hTrue_RecombFactor_ElTrk  ->Fill( Q_trk / float(Ni_trk) );
  hTrue_RecombFactor_ElShw  ->Fill( Q_shwr / float(Ni_shwr) );

  LOG_VERBATIM("MichelAnaFilter")<<"Done with particle tracker.";

}



//########################################################################################
// Fast-propagate a photon to a PMT, given some visibility, propagation time info
void MichelAnaFilter::PropagatePhoton(float t0, float tlar, float tmin, float tmean, float trms, bool isElectron, TH1D* h, TH1D* h_el){
  
  if( tmean < tmin ) return;
  double time = t0;
    
  // add the LAr time component (23% fast, 77% slow for MIPs)
  //   fast light: 6 +/- 2 ns
  //   late light: 1590 +/- 100 ns 
  //   ( J Chem Phys vol 91 (1989) 1469 E Morikawa et al )
  time  += fRand->Exp(tlar);
    
  // add physical propagation time
  double dt = 0.;
  while( dt <= tmin ) { 
    dt = 0.; 
    dt = fRand->Gaus(tmean,trms);
  }
  time += dt;
      
  // add TPB timing based on results from E. Segreto (arxiv:1411.4524)
  //   60% of emission has (5 +/- 5)ns lifetime
  //   30% of emission has (49 +/- 1)ns lifetime
  //   8% of emission has (3550 +/- 500)ns lifetime
  //   2% of emission has (309 +/- 10)ns lifetime
  double r = fRand->Rndm();
      if      ( 0.    <= r && r < 0.60 )  { time += fRand->Exp(5.); }
      else if ( 0.60  <= r && r < 0.90 )  { time += fRand->Exp(49.);  }
      else if ( 0.90  <= r && r < 0.98 )  { time += fRand->Exp(3550); }
      else			          { time += fRand->Exp(309.); }
    
  h ->Fill(time);
  if(isElectron) h_el ->Fill(time);

}



//########################################################################################
// Given a vector of hits and some reference 'X' coordinate, find the hit that is best
// matched to this reference. 
void MichelAnaFilter::FindBestMatchedHit( int plane, float x, int &bestMatchedHit, float &minDist){
  std::vector<int> hits;
  FindBestMatchedHit(plane, x, bestMatchedHit, minDist, hits);
}
void MichelAnaFilter::FindBestMatchedHit( int plane, float x, int &bestMatchedHit, float &minDist, std::vector<int>& hits){
  bestMatchedHit = -9;
  minDist  = 9999.;
  int N = fHitX.size();
  for(int i=0; i<N; i++ ){
    if( fHitPlane.at(i) != plane ) continue;
    // If the vector of hits provided is > 0, then require
    // that the hit be one of them
    if( hits.size() > 0 && (std::find(hits.begin(), hits.end(), i) == hits.end() )) {
      continue;
    }
    float dist  = fabs( fHitX.at(i) - x );
    if( dist < minDist ){
      minDist = dist;
      bestMatchedHit = i;
    }
  }
}



//########################################################################################
// The primary all-in-one clustering routine which takes in some empty michelcluster data
// object and performs a proximity clustering procedure. A 2D shower is reconstructed if the
// mu-electron boundary is found. This function requires that the candidate stopping mu track
// was already identified and the class member variables fMuTrackVertex,fMuTrackHits, and 
// fMuTrackIndex are defined.

void MichelAnaFilter::Clustering( int plane, michelcluster& cl) { 
  if( fMuTrackIndex < 0 ) return;
  Clustering( plane, cl, float(fMuTrackVertex.X()), fMuTrackHits);
}

void MichelAnaFilter::Clustering( int plane, michelcluster& cl, float muTrackX, std::vector<int> &muTrackHits) { 

  // Must have already done 3D tracking and tagged a candidate muon 
  // track for determination of seed, and also calculated X/W positions for all hits.
  if( fHitX.size() == 0 ) return;

  // Assign plane
  cl.plane = plane;
  LOG_VERBATIM("MichelAnaFilter")
  <<"Begin clustering on plane "<<plane<<"\n"
  <<"Looking for hit best matched to mu track startX = "<<muTrackX<<" cm .... "<<muTrackHits.size()<<" total hits associated with mu trk";
 
  // -------------------------------------------------------------
  // Find hit that is best matched to the muon's vertex X coordinate
  int muStartHit = -9;
  float minDistStart = 9999.;
  FindBestMatchedHit( plane, muTrackX, muStartHit, minDistStart, muTrackHits );
  LOG_VERBATIM("MichelAnaFilter")
  <<"Mu start hit dist = "<<minDistStart<<", key = "<<muStartHit;
  if( muStartHit < 0 ) return; 
  
  // --------------------------------------------------------------- 
  // Perform clustering of the muon + electron hits, using the hit
  // best matched to the reco vertex of our stopping track as the seed.
  std::vector<int> cluster;
  int seedHit = muStartHit;
  if( minDistStart  < 5.0 ) ProximityCluster(fHitX, fHitW, fHitCharge, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  
  // Starting again at the seedpoint, try clustering again (excluding hits 
  // previously clustered) to see if the initial clustering was incomplete.
  std::vector<int> cluster2 = cluster;
  ProximityCluster(fHitX, fHitW, fHitCharge, fHitPlane, cluster2, plane, seedHit, fMaxHitSeparation);

  // If cluster2 is same size as cluster, then no new hits were added in the
  // reclustering and we should just use the original cluster.
  // If cluster2 is larger, it means the original cluster was incomplete and 
  // we need to redefine the cluster using new endpoint
  if( cluster2.size() > cluster.size() ) {
    seedHit = cluster2.at(cluster2.size()-1);
    cluster.clear(); 
    ProximityCluster(fHitX, fHitW, fHitCharge, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  } 
 
  /* 
  // One final check to ensure the endpoint closest to the mu start X is
  // the one used as the seed    
  float x_end   = fHitX.at(cluster[cluster.size()-1]);
  float x_start = fHitX.at(cluster[0]);
  if( fabs(x_end - fMuTrackVertex.X()) < fabs(x_start - fMuTrackVertex.X()) ) {
    seedHit = cluster.at(cluster.size()-1);
    cluster.clear(); 
    ProximityCluster(fHitX, fHitW, fHitCharge, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  }
  */

  // Assign seed hit and cluster
  cl.seedHit = seedHit; 
  cl.cluster = cluster;

  // Assign 2D start-point for electron 
  TVector3 elStart2D, elDir2D;
  elStart2D .SetXYZ(-99.,-99.,-99.); 
  elDir2D   .SetXYZ(-99.,-99.,-99.); 
 
  // Find total charge of cluster
  cl.totalCharge        = 0.; 
  for(size_t i=0; i<cluster.size(); i++) cl.totalCharge += fHitCharge[cluster.at(i)];
  
  // -------------------------------------------------------------------
  // if cluster meets min required size, compute profiles
  if( (int)cluster.size() >= fMinClusterSize ) {
    
    // -------------------------------------------------------------------
    // Create dQ, dQ/ds, X, W profiles
    float s = 0.;
    size_t             nprof = cluster.size()-1;
    std::vector<float> profile_ds(nprof, 0.);
    std::vector<float> profile_s(nprof, 0.);
    std::vector<float> profile_dQ(nprof, 0.);
    std::vector<float> profile_dQds(nprof, 0.);
    std::vector<float> profile_X(nprof,0.);
    std::vector<float> profile_W(nprof,0.);
    std::vector<int>   profile_hitKey(nprof, -9);
    for(size_t i=1; i<cluster.size(); i++){
      TVector3 loc(fHitW[cluster.at(i)], fHitX[cluster.at(i)], 0.);
      TVector3 loc_prev( fHitW[cluster.at(i-1)], fHitX[cluster.at(i-1)],0.);
      float ds = (loc - loc_prev).Mag();
      s += ds;
      profile_hitKey.at(i-1)= cluster.at(i);
      profile_dQ.at(i-1)    = fHitCharge[cluster.at(i)];
      profile_X.at(i-1)     = fHitX[cluster.at(i)];
      profile_W.at(i-1)     = fHitW[cluster.at(i)];
      profile_dQds.at(i-1)  = fHitCharge[cluster.at(i)] / ds;
      profile_s.at(i-1)     = s;
      profile_ds.at(i-1)    = ds;
    }

    // -------------------------------------------------------------------
    // Make truncated mean version of dQ, dQds to smooth it out
    std::vector<float> profile_dQds_t(nprof, 0.);
    for( size_t i=0; i<nprof; i++)
      profile_dQds_t.at(i) = CalcTruncatedMeanInProfile(profile_dQds, i, fTruncMeanWindow, fTruncMeanP);
    
    // -------------------------------------------------------------------
    // Make local linearity profile
    std::vector<float> profile_linearity(nprof, 0.);
    for(size_t i=0; i<nprof; i++)
      profile_linearity.at(i) = CalcLocalLinearity(profile_X, profile_W, i, fLocalLinearityWindow);


    // ======================================================================
    // Determine boundary based on dQ, linearity profiles
    //  - find max in truncated dQ
    //  - search neighborhood for largest dQ hit value
    //  - require this coincide with dip in linearity (<0.80) 
   
    // First locate the truncated dQ/ds maximum, ignoring the edges
    int bnd1 = -99, bnd2 = -99;
    float dQmax = -9., dQtmax = -9.;
    size_t k1 = 2;
    size_t k2 = (nprof-1)-2;
    for(size_t i=k1; i<=k2; i++){
      if( profile_dQds_t.at(i) > dQtmax ) {
        dQtmax = profile_dQds_t.at(i);
        bnd1 = i;
      }
    }

    LOG_VERBATIM("MichelAnaFilter") 
    <<"Determining boundary point in dQds, dQds_t_max = "<<dQtmax<<" at index "<<bnd1;

    // ----------------------------------------------------------------------
    // The truncated mean peak will precede the charge drop-off boundary, 
    // so only look forward
    if( bnd1 > 0 ) { 
      for(size_t i=0; i< nprof; i++){
        int dist = int(i) - bnd1;
        if( ( dist >=0 ) && ( dist <= fMaxMatchDistance ) && profile_dQds.at(i) > dQmax ){
          dQmax = profile_dQds.at(i);
          bnd2 = i;
        }
      }
    }
    
    LOG_VERBATIM("MichelAnaFilter") 
    <<"dQds_max = "<<dQmax<<" at index "<<bnd2;

    hDistMaxTdQdsToMaxQ->Fill(bnd2 - bnd1);
    if( bnd2 >= 0 ) {
      cl.maxQ_i = bnd2;
      cl.maxQ_s = profile_s.at(bnd2);
    }
    LOG_VERBATIM("MichelAnaFilter") 
    <<"maxQ_s = "<<cl.maxQ_s;
   
    // ----------------------------------
    // Calculate the mean/RMS of the linearity of the hits preceding this one
    float linearityThresh = fLinThresh;
    if( fLinTol > 0 && bnd1+1 >= 10 ) {
      float linMean = -999.; 
      float linRMS  = 999.;
      int   N = 0;
      float sum = 0.;
      float sum_rms = 0.;
      for(int i=0; i<bnd1; i++) { sum += profile_linearity.at(i); N++;}
      linMean = sum / N;
      for(int i=0; i<bnd1; i++) sum_rms += pow( profile_linearity.at(i)-linMean, 2);
      linRMS = sqrt( sum_rms / N );
      float th = linMean - fLinTol * linRMS;
      if( th > fLinThresh ) linearityThresh = th;
    }
    
    LOG_VERBATIM("MichelAnaFilter")
    <<"Using linearity threshold: "<<linearityThresh;

     // ----------------------------------
    // Require that the candidate boundary hit resides in a region of 
    // low linearity.
    float minCov = 9;
    for(int i = std::max(int(0), bnd2-1); i<= std::min(int(nprof-1), bnd2+1); i++){
      if( profile_linearity.at(i) < minCov ) minCov = profile_linearity.at(i);
    }
    cl.minCovAtBnd = minCov;
    if(  minCov > linearityThresh ) bnd2 = -9;

    LOG_VERBATIM("MichelAnaFilter")<<"Lin threshold = "<<linearityThresh<<", minCov = "<<cl.minCovAtBnd<<" ... was boundary found?";

    if( bnd2 >= 0 ) {
      LOG_VERBATIM("MichelAnaFilter")<<"YES!  bnd2 = "<<bnd2;
      if( fBndOffset != 0 ){
        LOG_VERBATIM("MichelAnaFilter")<<"Offsetting bnd by "<<fBndOffset<<" hits";
        bnd2 += fBndOffset;
      }
      cl.bnd_i = bnd2;
      cl.bnd_s = profile_s.at(bnd2);
      cl.muEnd2D.SetXYZ(profile_W.at(bnd2), profile_X.at(bnd2), 0.);
      cl.muEndHit = profile_hitKey.at(bnd2); 
    }
   
    // Assign profiles to michelcluster object
    cl.prof_hitKey  = profile_hitKey;
    cl.prof_X       = profile_X;
    cl.prof_W       = profile_W;
    cl.prof_s       = profile_s;
    cl.prof_dQds    = profile_dQds;
    cl.prof_dQds_t  = profile_dQds_t;
    cl.prof_lin     = profile_linearity;
       
     
    // ======================================================================
    // Quality check calculations:
    //  - cluster sizes (total and mu/el subclusters)
    //  - linearity of muon track
    //  - electron-tagged hit cluster should not be co-linear w/ muon
    //  - check for unclustered ("extra") hits around mu endpoint
    //  - Bragg peak area (really necessary?)
    
    if( cl.bnd_i >= 0 ) {
     
      // -----------------------------------------------------------------
      // Cluster together the "muon" and "electron" portions of the cluster
      TGraph g_mu, g_el;
      TVector3 muEndDir_Pt1(-9.,-9.,0.);
      TVector3 muEndDir_Pt2 = cl.muEnd2D;
      cl.cluster_mu.push_back(seedHit);
      for(size_t i=0; i< nprof ; i++){
        if( (int)i <= bnd2 ) {
          cl.cluster_mu.push_back(profile_hitKey.at(i));
          if( profile_linearity.at(i) > fLinThreshMuFit && bnd2-(int)i < 20 ) {
            if( muEndDir_Pt1.X() == -9. ) muEndDir_Pt1.SetXYZ(profile_W.at(i), profile_X.at(i), 0.);
            muEndDir_Pt2.SetXYZ(profile_W.at(i), profile_X.at(i), 0.);
            g_mu.SetPoint(g_mu.GetN(),profile_W.at(i), profile_X.at(i));
          }
        } else { 
          cl.cluster_el.push_back(profile_hitKey.at(i));
          g_el.SetPoint(g_el.GetN(),profile_W.at(i), profile_X.at(i));
        }
      }
     
      int muHits = cl.cluster_mu.size();
      int elHits = cl.cluster_el.size();
      cl.nPtsMuFit = g_mu.GetN();
      hMuClusterHitsEndFit->Fill(cl.nPtsMuFit);
     
      // If we have > 0 electron track hits in the cluster, define the electron
      // start-point as the first hit after the boundary point 
      if( elHits > 0 ) elStart2D.SetXYZ( fHitW.at(cl.cluster_el[0]), fHitX.at(cl.cluster_el[0]), 0.);

      // ------------------------------------------------------------------
      // Count up hits in muon segment of cluster and find fraction of hits
      // with linearity above threshold
      int numMuHitsGood = 0;
      float sum_lin = 0;
      int N = 0;
      for(int i=0; i<=bnd2; i++){
        sum_lin += profile_linearity.at(i);
        N++;
        if( profile_linearity.at(i) > fLinThreshMuFit ) numMuHitsGood++;  
      }
      if( N > 0 ) cl.muAveLinearity = sum_lin / float(N);
      if( muHits > 0 ) cl.fracMuHitsLinear = (float)numMuHitsGood / muHits;
      
      LOG_VERBATIM("MichelAnaFilter")
      <<"Linearity check... good mu hits = "<<numMuHitsGood<<" (frac = "<<cl.fracMuHitsLinear<<")";
    
      // --------------------------------------------------------------------
      // Count extra unclustered hits around muon endpoint and electron startpoint
      cl.extraHits   = 0;
      for(size_t i=0; i<fHitX.size(); i++){
        if( fHitPlane.at(i) != plane ) continue; 
        if( std::find(cluster.begin(), cluster.end(), i ) != cluster.end() ) continue;
        TVector3 loc( fHitW.at(i), fHitX.at(i), 0. );
        float distMu = (loc - cl.muEnd2D).Mag();
        float distEl = (loc - elStart2D).Mag();
        if( distMu <= fMuEndRadius || distEl <= fMuEndRadius ) cl.extraHits++;
      }

      // ------------------------------------------------------------------
      // Find the 2D direction of mu/el subclusters
      if( cl.nPtsMuFit  >= fMinMuClusterHitsEndFit ){
        
        // First compute a "simple" (ie, dumb)  muon vector by drawing
        // a vector btwn first and last hits to be used in the fit
        TVector3 muEndDir_simple = muEndDir_Pt2 - muEndDir_Pt1;
        
        // Do linear fit to the good muon hits and get the slope 
        TF1 muTr("muTr","[0]*x+[1]",0.,100.);
        muTr.SetParameter(0, 1.);
        muTr.SetParameter(1, 1.);
        g_mu.Fit("muTr","Q");
        float m = muTr.GetParameter(0);
        
        // Convert slope into a TVector
        cl.muEndDir2D.SetXYZ( 1., m, 0. );
        cl.muEndDir2D.SetMag(muEndDir_simple.Mag());

        // Use the simple vector to get a rough idea of the angle,
        // and reverse the slope-calculated vector if necessary
        float ang = muEndDir_simple.Angle(cl.muEndDir2D) * RAD_TO_DEG;
        if( ang > 90 ) {
          cl.muEndDir2D.SetX( -1.* cl.muEndDir2D.X() );
          cl.muEndDir2D.SetY( -1.* cl.muEndDir2D.Y() );
        }

        // Now calculate the direction of the electron cluster, but this time just
        // draw a vector from the starting point to the charge-weighted centroid
        if( elHits >= 2 ) {
          float sumW = 0., sumX = 0., sumWeight = 0.;
          for(auto hit : cl.cluster_el ){
            sumX += fHitX.at(hit)*fHitCharge.at(hit);
            sumW += fHitW.at(hit)*fHitCharge.at(hit);
            sumWeight += fHitCharge.at(hit);
          }
          float meanX = sumX / sumWeight;
          float meanW = sumW / sumWeight;
          TVector3 elDir_Pt1 = elStart2D;
          TVector3 elDir_Pt2(meanW, meanX, 0.);
          elDir2D = elDir_Pt2 - elDir_Pt1;
          cl.decayAngle2D = cl.muEndDir2D.Angle(elDir2D) * RAD_TO_DEG;
          LOG_VERBATIM("MichelAnaFilter")
          <<"Angle between muon and electron clusters: "<<cl.decayAngle2D;
        }

      }//<-- endif enough muon points for directional fit

    }//<-- endif a boundary point was determined
  
  }//<-- endif cluster size cut
    

  // =====================================================================
  // Electron shower clustering: if a decay angle was calculated above, then
  // proceed to cluster all hits along the reconstructed direction in 2D.
  
  if( cl.decayAngle2D > 0 ){
    float sum_t = 0.; 
    float sum_x = 0.;
    float sum_w = 0.;
  
    // histogram to keep track of anglular spread of shower
    // (charge-weighted)
    TH1D hAng("hAng","ang",180,0,180);
  
    for(size_t i=0; i<fHitX.size(); i++){
       
        if( fHitPlane.at(i) != plane ) continue;
        
        // Skip unphysical hit times, as well as any hits with small X
        if( fHitT.at(i) < 0 || fHitX.at(i) < 1. ) continue;
        
        // Skip muon hits
        if( std::find(cl.cluster_mu.begin(), cl.cluster_mu.end(), i ) != cl.cluster_mu.end() ) continue;

        bool isInElCluster = false;
        if( std::find(cl.cluster_el.begin(), cl.cluster_el.end(), i ) != cl.cluster_el.end() ) isInElCluster = true;
        
        // Find angle of this hit relative to Michel direction (2D)
        TVector3 loc(fHitW.at(i), fHitX.at(i), 0.);
        TVector3 hv = loc - elStart2D; 
        float ang = hv.Angle(elDir2D) * RAD_TO_DEG;
        
        if( cl.plane == 1 ) hElShowerDepAngle2D->Fill(ang);
        
        // If this part of the original electron-portion of total cluster, or if 
        // it lies within 2D acceptance angle cone, add it to the shower.
        if( isInElCluster || ang <= fShowerAcceptanceAngle ) {
          cl.shower     .push_back( i );
          cl.isInTrk    .push_back(isInElCluster);
          hAng          .Fill(ang,fHitCharge.at(i));
          sum_t         += fHitT.at(i) * fHitCharge.at(i);
          sum_x         += fHitX.at(i) * fHitCharge.at(i);
          sum_w         += fHitCharge.at(i);
        }//<-- end if in 2D cone
        
    }//<-- end loop over hits

    if( sum_w > 0. ) {
      cl.aveDriftTime  = sum_t / sum_w;
      cl.aveX          = sum_x / sum_w;
    }
  
    cl.elShowerFrac = cl.shower.size() / float( fNumPlaneHits[plane] - (int)cl.cluster_mu.size());
    cl.showerAngleMean     = hAng.GetMean();
    cl.showerAngleRMS     = hAng.GetRMS();
  
  } // end Michel shower 2D clustering (if decay angle found)

}



//########################################################################################
void MichelAnaFilter::CalcMichelShowerEnergy( michelcluster& cl ) { 

  cl.elShowerCharge     = 0.;
  cl.elShowerEnergy     = 0.;
  cl.elTrackEnergy      = 0.;
  cl.elTrackCharge      = 0.;

  for(size_t i=0; i<cl.shower.size(); i++){
    
    bool    isInTrack   = cl.isInTrk.at(i);
    size_t  hitID       = cl.shower.at(i);

    float rfac = fRecombFactorShower;
    if( isInTrack ) rfac = fRecombFactorTrack;
          
    // Add to Michel charge
    cl.elShowerCharge   += fHitCharge.at(hitID);
    if( isInTrack ) 
      cl.elTrackCharge  += fHitCharge.at(hitID);
         
    // Scale charge to energy
    //   dE = N_i * W_ion
    //   N_i = dQ/R
    //   W_ion = W_ph * (1+alpha)
    float dE = ( fHitCharge.at(hitID) / rfac ) * (fWion * 1e-6);
    cl.elShowerEnergy   += dE;
    if( isInTrack ) 
      cl.elTrackEnergy  += dE;

  }//<-- end loop over hits
 
  LOG_VERBATIM("MichelAnaFilter")<<"Plane "<<cl.plane
  <<": Michel electron charge= "<<cl.elShowerCharge
  <<", energy= "<<cl.elShowerEnergy;

}



//########################################################################################
// A simple clustering algorithm that starts at some "seed" hit and clusters hits based on
// proximity in 2D 'WX' space.
void MichelAnaFilter::ProximityCluster(
std::vector<float>& hitX, std::vector<float>& hitW, std::vector<float>& hitCharge, std::vector<int>& hitPlane, 
std::vector<int>& cluster, int plane, int seedHit, float distThresh)
{
  LOG_VERBATIM("MichelAnaFilter")
  <<"Beginning proximity-based clustering; input cluster has "<<cluster.size()<<" hits to start with.\n"
  <<"  seed hit = "<<hitW[seedHit]<<","<<hitX[seedHit]<<"\n"
  <<"  looking at plane "<<plane;
 
  // Count all hits that are in the plane of choice.
  size_t hitsInPlane = 0;
  for(size_t i=0; i<hitX.size(); i++){
    if( hitPlane.at(i) == plane ) hitsInPlane++;
  }
  
  // Check if the seed hit is already in cluster; if not, add it
  if( std::find(cluster.begin(), cluster.end(), seedHit ) == cluster.end() ) 
    cluster.push_back(seedHit);

  std::vector<int> cluster1 = cluster;
  std::vector<int> cluster2 = cluster;
  bool clustering;
  bool firstPass;

  // ----------------------------------------------------------------
  // Round one of clustering
  clustering = true;
  firstPass  = true;
  while(clustering){
    float minDist = 9999.;
    int   minKey  = -9.;
    size_t prevHit = cluster1[cluster1.size()-1];
    if( firstPass ) { prevHit = seedHit; firstPass = false; }
    // Loop over all hits and update
    //std::cout<<"  beginning loop over hits (cluster1 size = "<<cluster1.size()<<")\n";
    for(size_t i=0; i<hitX.size(); i++){
      // Skip hits not in plane, or already in cluster
      if( hitPlane.at(i) != plane ) continue;
      if( std::find(cluster1.begin(), cluster1.end(), i ) != cluster1.end() ) continue;
      TVector3 hit_loc(hitX[i], hitW[i], 0.);
      TVector3 hit_loc_prev(hitX[prevHit], hitW[prevHit], 0.);
      TVector3 dir = hit_loc - hit_loc_prev;
      float dist = dir.Mag();
      if( dist < minDist ) { minDist = dist; minKey = i; }
    }
    if( minKey >= 0 ) hClusterHitSeparation->Fill(minDist);
    // If we found a closest hit, add it to cluster
    if( minKey >= 0 && minDist < distThresh ){
      cluster1.push_back( minKey );
    // Otherwise, stop clustering
    } else {
      clustering = false;
    }
  }
  
  // ----------------------------------------------------------------
  // Now repeat, but this time excluding the hits that were clustered
  // together in the first round
  clustering = true;
  firstPass  = true;
  while(clustering){
    float minDist = 9999.;
    int   minKey  = -9.;
    size_t prevHit = cluster2[cluster2.size()-1];
    if( firstPass ) { prevHit = seedHit; firstPass = false; }
    // Loop over all hits and update
    //std::cout<<"  beginning loop over hits (cluster2 size = "<<cluster2.size()<<")\n";
    for(size_t i=0; i<hitX.size(); i++){
      // Skip hits not in plane, or already in cluster
      if( hitPlane.at(i) != plane ) continue;
      if( std::find(cluster1.begin(), cluster1.end(), i ) != cluster1.end() ) continue;
      if( std::find(cluster2.begin(), cluster2.end(), i ) != cluster2.end() ) continue;
      TVector3 hit_loc(hitX[i], hitW[i], 0.);
      TVector3 hit_loc_prev(hitX[prevHit], hitW[prevHit], 0.);
      TVector3 dir = hit_loc - hit_loc_prev;
      float dist = dir.Mag();
      if( dist < minDist ) { minDist = dist; minKey = i; }
    }
    if( minKey >= 0 ) hClusterHitSeparation->Fill(minDist);
    // If we found a closest hit, add it to cluster
    if( minKey >= 0 && minDist < distThresh ){
      cluster2.push_back( minKey );
    // Otherwise, stop clustering
    } else {
      clustering = false;
    }
  }

  // --------------------------------------
  // Pick the larger cluster
  if( cluster1.size() > cluster2.size() ) { cluster = cluster1; } 
  else {                                    cluster = cluster2; }
 
  LOG_VERBATIM("MichelAnaFilter")
  <<"  total hits in plane "<<plane<<": "<<hitsInPlane<<", cluster size: "<<cluster.size();
}



//########################################################################################
// Calculate the truncated mean in the neighborhood of a sample within a clustered hit profile.
float MichelAnaFilter::CalcTruncatedMeanInProfile(std::vector<float>& v, size_t index, int nb, float p){
  
  //???
  // Resize neighborhood depending on how close to edge we are
  //nb = std::min( nb, (int)index );
  //nb = std::min( nb, (int)(v.size()-index-1));
  //???
  
  // Restrict neighborhood (nb) at edges
  if( nb > 1 && (index < 2 || index >= v.size() - 2 )) nb = 1;

  // Make vector of neighborhood values
  int k1 = std::max(int(0), int(index - nb));
  int k2 = std::min(int(v.size()-1), int(index + nb));
  std::vector<float> neighborhood;
  for(int i=k1; i<=k2; i++) {
    neighborhood.push_back(v.at(i));
    // if we're at an edge, add extra weight to the edge point
    if( i == 0 || i == (int)v.size() - 1 ) neighborhood.push_back(v.at(i));
  }
 
  return fOpHitBuilderAlg.CalcTruncatedMean(neighborhood, p);
}



//########################################################################################
// Calculate the local linearity within a neighborhood of a sample within a clustered hit profile.
float MichelAnaFilter::CalcLocalLinearity(std::vector<float>& vx, std::vector<float>& vy, size_t index, int nb){

  // Don't even bother if vectors don't match up in size,
  // or if we don't have at least 3 points
  if( vx.size() < 3 || vx.size() != vy.size() ) return 0.;
  
  // Restrict neighborhood (nb) at edges
  if( nb > 1 && (index < 2 || index >= vx.size() - 2 )) nb = 1;
  
  //???
  // Resize neighborhood depending on how close to edge we are
  //nb = std::min( nb, (int)index );
  //nb = std::min( nb, (int)(vx.size()-index-1));
  //???

  // Define index range
  int k1  = std::max(int(0), int(index - nb));
  int k2  = std::min(int(vx.size()-1), int(index + nb));

  // Average X and Y
  int N   = 0;
  float sumX = 0., sumY = 0.;
  for(int i=k1; i<=k2; i++){
    sumX += vx.at(i);
    sumY += vy.at(i);
    N++;
  }
  float aveX = sumX / N;
  float aveY = sumY / N;

  // Covariance, std dev
  float cov = 0.;
  float sum_dX2 = 0.;
  float sum_dY2 = 0.;
  for(int i=k1; i<=k2; i++){
    cov     += ((vx.at(i) - aveX) * (vy.at(i) - aveY)) / N;
    sum_dX2 += pow(vx.at(i) - aveX, 2);
    sum_dY2 += pow(vy.at(i) - aveY, 2);
  }
  float stDevX = sqrt( sum_dX2 / N );
  float stDevY = sqrt( sum_dY2 / N );
  if( stDevX*stDevY != 0 ) 
    return pow( fabs(cov) / (stDevX * stDevY), 2);
  else 
    return 1.;
}



//########################################################################################
// This function performs a fit to the (binned) muon late-light and extrapolates this fit
// forward to estimate how much of this light spills over into the Michel electron integration
// window. It then subtracts this off of the integrated PE for the Michel.
//  TODO: Make number of bins/bin widths variable
void MichelAnaFilter::MuContamCorr(int ch) {
  
  LOG_VERBATIM("MichelAnaFilter")
  <<"Correcting for muon late-light...";
  
  const int N = 7;
  float binW[N] = { 200., 200., 200., 200., 300.,   300., 300.}; 
  float x[N]    = { 200., 400., 600., 800., 1050., 1350., 1650.};
  float y[N];
  float dx[N], dy[N];
  
  // Check that there are two hits, separated by at least (500ns?)
  if( fvHitTimes[ch].size() != 2 || fdT[ch] < x[N-1] + 50 || fSPE[ch] <= 0 || fPE_totalRaw[ch] < 0) return;

  float adc100 = fvHitADC_100ns[ch].at(0);
  float adc300 = fvHitADC_300ns[ch].at(0);
  float adc500 = fvHitADC_500ns[ch].at(0);
  float adc700 = fvHitADC_700ns[ch].at(0);
  float adc900 = fvHitADC_900ns[ch].at(0);
  float adc1200 = fvHitADC_1200ns[ch].at(0);
  float adc1500 = fvHitADC_1500ns[ch].at(0);
  float adc1800 = fvHitADC_1800ns[ch].at(0);

  y[0] = ((adc300-adc100) / fSPE[ch])/binW[0];
  y[1] = ((adc500-adc300) / fSPE[ch])/binW[1];
  y[2] = ((adc700-adc500) / fSPE[ch])/binW[2];
  y[3] = ((adc900-adc700) / fSPE[ch])/binW[3];
  y[4] = ((adc1200-adc900) / fSPE[ch])/binW[4];
  y[5] = ((adc1500-adc1200) / fSPE[ch])/binW[5];
  y[6] = ((adc1800-adc1500) / fSPE[ch])/binW[6];

  for(int i=0; i<N; i++){
    dx[i] = binW[i]/2.;
    if( y[i] >= 0 ) dy[i] = sqrt(y[i]*binW[i]) / binW[i];
    else            dy[i] = 0.;
    hAvePhelProfile_mu[ch]->Fill(x[i],y[i]);
    avePhelProfileCounter_mu[ch]++;
  }
  /*
  std::cout
  <<"  ADC integrals\n"
  <<"  100ns: "<<adc100<<"\n"
  <<"  300ns: "<<adc300<<"\n"
  <<"  500ns: "<<adc500<<"\n"
  <<"  700ns: "<<adc700<<"\n"
  <<"  900ns: "<<adc900<<"\n"
  <<"  1200ns: "<<adc1200<<"\n"
  <<"  1500ns: "<<adc1500<<"\n"
  <<"  1800ns: "<<adc1800<<"\n";

  std::cout
  <<"  PE bins (PEs per ns):\n"
  <<"  100-300: "<<y[0]<<"\n"
  <<"  300-500: "<<y[1]<<"\n"
  <<"  500-700: "<<y[2]<<"\n"
  <<"  700-900: "<<y[3]<<"\n"
  <<"  900-1200: "<<y[4]<<"\n"
  <<"  1200-1500: "<<y[5]<<"\n"
  <<"  1500-1800: "<<y[6]<<"\n";
  */

  TF1 fit("fit","[0]*exp(-x/[1])",300.,1800.);
  TGraphErrors gr(N,x,y,dx,dy);
  fit.SetParameter(0, y[0] );  
  fit.SetParameter(1, 1500. );
  fit.SetParLimits(1, 0., 5000. );
  if( fMuContamCorr_EffTau[ch] > 0 ) 
    fit.FixParameter(1, fMuContamCorr_EffTau[ch] );
  gr.Fit("fit","RN0");
  float tau = fit.GetParameter(1);
  float Chi2 = fit.GetChisquare() / fit.GetNDF(); 
  LOG_VERBATIM("MichelAnaFilter")<<"  tau = "<<tau<<", red chi2 = "<<Chi2;
  hMuTau[ch]->Fill(tau);
 
  // Save info to be plotted later with waveforms
  fvPrepulseX1[ch].at(1)        = fvHitTimes[ch].at(0) + 300;
  fvPrepulseSlowNorm[ch].at(1)  = fit.GetParameter(0)*fSPE[ch]*-1.;
  fvPrepulseZeroPoint[ch].at(1) = fvHitTimes[ch].at(0);
  fvPrepulseSlowTau[ch].at(1)   = fit.GetParameter(1);
  
  int k1 = fdT[ch]-10;
  int k2a = k1 + fPromptWindow + 10;
  int k2 = k1 + fFullWindow;
  float sum_100 = 0.;
  float sum_total = 0.;
  for(int i = k1; i < k2; i++){
    sum_total += fit.Eval(i);
    if( i < k2a ) sum_100 += fit.Eval(i); 
  }

  // Correct off the contamination
  fPE_prompt[ch] = fPE_promptRaw[ch] - sum_100;
  fPE_total[ch]  = fPE_totalRaw[ch] - sum_total;
  
  LOG_VERBATIM("MichelAnaFilter")
  <<"  pe100 contam: "<<sum_100<<", total contam: "<<sum_total<<"\n"
  <<"  raw prompt: "<<fPE_promptRaw[ch]<<", total: "<<fPE_totalRaw[ch]<<"\n"
  <<"  corrected prompt: "<<fPE_prompt[ch]<<", total: "<<fPE_total[ch];


}



//########################################################################################
int MichelAnaFilter::GetGlobalBin(const TH3D* h, double x, double y, double z){
  int xbin = h->GetXaxis()->FindBin(x);
  int ybin = h->GetYaxis()->FindBin(y);
  int zbin = h->GetZaxis()->FindBin(z);
  return h->GetBin(xbin,ybin,zbin);
}
int MichelAnaFilter::GetGlobalBin(const TH3D* h, const TLorentzVector& v){
  return GetGlobalBin(h, v.X(), v.Y(), v.Z() );
}
int MichelAnaFilter::GetGlobalBin(const TH3D* h, const TVector3& v){
  return GetGlobalBin(h, v.X(), v.Y(), v.Z() );
}



//########################################################################################
float MichelAnaFilter::GetVisibility(TVector3 &loc, int ch){
  return GetVisibility(loc, ch, 0) + GetVisibility(loc, ch, 1);
}



//########################################################################################
float MichelAnaFilter::GetVisibility(TVector3 &loc, int ch, int isVis ){
  float tmin = -999; // -999 = don't fetch time info
  float tmean = 0.;
  float trms = 0.;
  float vis = 0.;
  ReadPhotonLibrary(loc, ch, isVis, vis, tmin, tmean, trms);
  return vis;
}



//########################################################################################
void MichelAnaFilter::ReadPhotonLibrary(TVector3 &loc, int ch, int isVis, float& vis, float& tmin, float& tmean, float& trms){
  
  int   gbin  = GetGlobalBin(libr_XYZ_visibility[ch][isVis],loc);

  vis   =   fCollectionEff[ch]
            * fQE_ScaleFactor[ch]
            * libr_XYZ_visibility[ch][isVis]->GetBinContent(gbin);
  
  if( isVis == 0 ) vis *= fQE_ScaleFactor_VUV[ch];
  if( isVis == 1 ) vis *= fQE_ScaleFactor_Vis[ch];

  // Check for bad visibility
  if( vis <= 0. ) vis = 0.;

  // Get timing information for this voxel
  if( tmin != -999 ) {
    tmin  =   libr_XYZ_arrivalTimeMin[ch][isVis]->GetBinContent(gbin);  
    tmean =   libr_XYZ_arrivalTimeMean[ch][isVis]->GetBinContent(gbin);  
    trms  =   libr_XYZ_arrivalTimeRMS[ch][isVis]->GetBinContent(gbin);
  
    // Check for bad time values, and set to some reasonable defaults
    // if necessary
    if( tmin <= 0 || tmean <= 0 || trms <= 0 || tmean < tmin ) {
      tmin  = 2.;
      tmean = 5.;
      trms  = 3.;
    }
  }

}



//########################################################################################
float MichelAnaFilter::Integrate(const TH1D* h, float x1, float x2){
  int Nbins = h->GetNbinsX();
  double integral = 0.;
  for(int i=1; i<Nbins; i++ ){
    double x = h->GetBinLowEdge(i);
    if( x >= x1 and x <= x2 ) integral += h->GetBinContent(i);
  }
  return integral;
}



//########################################################################################
float MichelAnaFilter::Integrate(const std::vector<float> &vec, int x1, int x2){
  float integral = 0.;
  for(int i=x1; i<x2; i++ ) integral += vec.at(i);
  return integral;
}



//########################################################################################
void MichelAnaFilter::AddToAverageWfm( std::vector<float> &wfm, short t1, short t2, TH1D* h, int &counter, int polarity){
  int nbins = h->GetNbinsX();
  short T = t2 - t1;
  if( nbins ==  T ) {
    float xmin = h->GetXaxis()->GetXmin();
    for(short i=0; i<T; i++) h->Fill(xmin+i, polarity*wfm.at(t1+i) );
    counter++;
  } else {
    LOG_VERBATIM("MichelAnaFilter")<<"MichelAnaFilter::AddToAverageWfm ERROR: Number of bins ("<<nbins<<") does not match range given ("<<T<<")";
  }
}

//########################################################################################
void MichelAnaFilter::AddToAverageWfm( const TH1D* hin, short t1, short t2, short thit, TH1D* h, int &counter, int polarity){
  LOG_VERBATIM("MichelAnaFilter")<<"Adding to average MC waveform: "<<t1<<", "<<t2<<", thit = "<<thit;
  int nin = hin->GetNbinsX();
  double tmin = h->GetXaxis()->GetXmin();
  double tmax = h->GetXaxis()->GetXmax();
  for( short i=0; i<nin; i++){
    double t = hin->GetBinLowEdge(i) - double(thit);
    if  ( t < tmin ) { 
      continue; 
    } else {
      double n = hin->GetBinContent(i);
      h->Fill(t, n);
    }
    if( t > tmax ) break;
  }
  counter++;
}



//########################################################################################
// Plots PMT waveforms 
void MichelAnaFilter::MakeWfmGraphs(art::Event& e){

    for(size_t i=0; i<fPulses.size(); i++){
      size_t ch=fPulses.at(i).OpChannel();

      // -----------------------------------------------------------------------------------
      if( !fIsMC && fNumSavedWfmGraphs < fMaxSavedWfmGraphs ) {
        
        short x1 = 0;
        short x2 = 8500+fIntegrationWindows[fIntegrationWindows.size()-1];
        if( fNumOpHits[ch] > 0 ) x1 = std::max(0,fvHitTimes[ch].at(0) - 1000);

        fNumSavedWfmGraphs++;

        TGraph g_wfm;
        g_wfm.SetMarkerColor(kBlack);
        g_wfm.SetLineColor(kBlack);
        TGraph g_bs;
        g_bs.SetMarkerColor(kRed);
        g_bs.SetLineColor(kRed);
        TGraph g_grad;
        g_grad.SetMarkerColor(kViolet+2);
        g_grad.SetLineColor(kViolet+2);
        std::vector<float> g = fOpHitBuilderAlg.MakeGradient(fPMT_wfm[ch]);
        for(size_t i=0; i<fPMT_wfm[ch].size(); i++) {
          g_wfm .SetPoint(  g_wfm.GetN(),  i, -1.*fPMT_wfm[ch].at(i));
          g_grad.SetPoint(  g_grad.GetN(), i, -1.*g.at(i) );
          if( fvbs[ch].at(i) > -999. ) g_bs.SetPoint( g_bs.GetN(), i, -1.*fvbs[ch].at(i) );
        }
         
        sprintf(histName, "wfm_ch%lu_%i_r%i_sr%i_e%i", ch, fNumSavedWfmGraphs, e.run(), e.subRun(), e.id().event());
        fTCanvas     = graphDir_wfms.make<TCanvas>(histName,"c1",700,600);
        TPad c1("c1","c1",0.0,0.0,1.0,1.0);
        c1.Draw(); 
        c1.Divide(1,2);
          
        c1.cd(1);

        TMultiGraph mg;
        if(g_wfm.GetN() > 0 )    mg.Add(&g_wfm,"APL");
        mg.Draw("a");
        mg.SetTitle(histName);
        mg.GetYaxis()->SetTitle("Inverted PMT Signal [ADC]");
        mg.GetXaxis()->SetTitle("Time Tick [ns]");
        mg.GetXaxis()->SetLimits(x1,x2);
        mg.GetYaxis()->CenterTitle();
        mg.GetXaxis()->SetTitleSize(0.05);
        mg.GetYaxis()->SetTitleSize(0.05);
        mg.GetXaxis()->SetLabelSize(0.05);
        mg.GetYaxis()->SetLabelSize(0.05);
        mg.GetYaxis()->SetTitleOffset(0.8);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        gPad->Update();
        /*
        TLine hitLA[20];
        for(int i=0; i<std::min(20,fNumOpHits[ch]); i++){
          hitLA[i].SetX1( fvHitTimes[ch].at(i) );
          hitLA[i].SetX2( fvHitTimes[ch].at(i) );
          hitLA[i].SetY1( gPad->GetUymin() );
          hitLA[i].SetY2( gPad->GetUymax() );
          hitLA[i].SetLineColor(kRed);
          hitLA[i].SetLineWidth(1);
          hitLA[i].SetLineStyle(7);   
          hitLA[i].Draw();
        }
        g_wfm.Draw();
        */
        // Draw pre-pulse fit
        TF1 pfit("pfit","[0] + [1]*exp(-(x-[2])/[3])",0.,30000.);
        pfit.SetLineColor(kRed);
        pfit.SetLineStyle(1);
        pfit.SetLineWidth(1);
        TPaveText pt(0.6,0.7,0.9,0.9);
        pt.SetTextAlign(33);
        pt.SetBorderSize(0);
        pt.SetFillColor(0);
          //short T_mu = fvHitTimes[ch].at(0);
          short T_el = fvHitTimes[ch].at(1);
          pfit.SetRange(fvPrepulseX1[ch].at(1),T_el+fIntegrationWindows[fIntegrationWindows.size()-1]);
          pfit.SetParameter(0,0); // baseline already subtracted
          pfit.SetParameter(1,-1.*fvPrepulseSlowNorm[ch].at(1));
          pfit.SetParameter(2,fvPrepulseZeroPoint[ch].at(1));
          pfit.SetParameter(3,fvPrepulseSlowTau[ch].at(1));
          pfit.Draw("same");
          sprintf(buffer,"#DeltaT =  %f", fdT[ch]); pt.AddText(buffer);
          sprintf(buffer,"Prompt = %f PE", fPE_prompt[ch]); pt.AddText(buffer);
          sprintf(buffer,"Full = %f PE", fPE_total[ch]); pt.AddText(buffer);
          pt.Draw();

        c1.cd(2);
        g_grad.Draw("AL"); 
        g_grad.GetYaxis()->SetTitle("Signal Gradient");
        g_grad.GetXaxis()->SetTitle("Time Tick [ns]");
        g_grad.GetXaxis()->SetLimits(x1,x2);
        g_grad.GetYaxis()->CenterTitle();
        g_grad.GetXaxis()->SetTitleSize(0.05);
        g_grad.GetYaxis()->SetTitleSize(0.05);
        g_grad.GetXaxis()->SetLabelSize(0.05);
        g_grad.GetYaxis()->SetLabelSize(0.05);
        g_grad.GetYaxis()->SetTitleOffset(0.8);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        gPad->Update();
        TLine hitLB[20];
        for(int i=0; i<std::min(20,fNumOpHits[ch]); i++){
          hitLB[i].SetX1( fvHitTimes[ch].at(i) );
          hitLB[i].SetX2( fvHitTimes[ch].at(i) );
          hitLB[i].SetY1( gPad->GetUymin() );
          hitLB[i].SetY2( gPad->GetUymax() );
          hitLB[i].SetLineColor(kRed);
          hitLB[i].SetLineWidth(1);
          hitLB[i].SetLineStyle(7);   
          hitLB[i].Draw();
        }
        g_grad.Draw();

        
        
        fTCanvas->Update();
        fTCanvas->Write(histName);
        
      }//<-- done drawing graphs
    } // Done looop PMT
    // -----------------------------------------------------------
}

//########################################################################################
void MichelAnaFilter::MakeClusteringGraphs( michelcluster& cl ) {
  
  if( fNumSavedHitGraphs >= fMaxSavedHitGraphs ) return;
  
  int plane = cl.plane;
  if ( plane < 0 || plane > 1 ) return;
 
  size_t nprof = cl.prof_dQds.size();
  if( nprof == 0 ) return;

  float elShowerEnergy = cl.elShowerEnergy;

  // Add all hits to a graph
  TGraph        g_hits;
  g_hits.SetMarkerColor(kBlack);
  g_hits.SetMarkerStyle(20);
  g_hits.SetMarkerSize(0.9);
  for(size_t i=0; i<fHitX.size(); i++)
    if( fHitPlane.at(i) == plane ) g_hits.SetPoint(g_hits.GetN(),fHitW.at(i),fHitX.at(i));
  
  // Mark the seed hit
  TGraph g_seed;
  g_seed.SetMarkerColor(kGreen+2);
  g_seed.SetMarkerStyle(29);
  g_seed.SetMarkerSize(2.0);
  if( cl.seedHit >= 0 ) g_seed.SetPoint(0,fHitW.at(cl.seedHit), fHitX.at(cl.seedHit));
    
  // Add clustered hits to graph
  TGraph g_cluster;
  g_cluster.SetMarkerColor(kBlue);
  g_cluster.SetLineColor(kBlue);
  g_cluster.SetMarkerStyle(20);
  g_cluster.SetMarkerSize(0.5);
  for(size_t i=0; i<cl.cluster.size(); i++) g_cluster.SetPoint(g_cluster.GetN(),fHitW.at(cl.cluster[i]),fHitX.at(cl.cluster[i]));
   
  // Mark hits that were included in the Michel shower
  TGraph g_shower;
  g_shower.SetMarkerColor(kOrange+7);
  g_shower.SetMarkerStyle(4);
  g_shower.SetMarkerSize(1.4);
  for(size_t i=0; i<cl.shower.size(); i++) g_shower.SetPoint(g_shower.GetN(),fHitW.at(cl.shower[i]),fHitX.at(cl.shower[i]));
    
  // Make dQ and truncated mean dQ graphs 
  TGraph g_dQt;
  TGraph g_dQ;
    g_dQ.SetMarkerColor(kBlack);
    g_dQ.SetLineColor(kBlack);
    g_dQ.SetLineStyle(2);
    g_dQ.SetMarkerStyle(20);
    g_dQ.SetMarkerSize(0.4);
    g_dQt.SetMarkerColor(kBlue-4);
    g_dQt.SetLineColor(kBlue-4);
    g_dQt.SetMarkerStyle(20);
    g_dQt.SetMarkerSize(0.8);
    g_dQt.SetLineWidth(2);
    for(size_t i=0; i<nprof; i++){
      g_dQ.SetPoint(g_dQ.GetN(), cl.prof_s.at(i), cl.prof_dQds.at(i) );
      g_dQt.SetPoint(g_dQt.GetN(),cl.prof_s.at(i),cl.prof_dQds_t.at(i));
    }
    
  TGraph g_lin;
  g_lin.SetMarkerColor(kBlue-4);
  g_lin.SetLineColor(kBlue-4);
  g_lin.SetMarkerStyle(20);
  g_lin.SetMarkerSize(0.8);
  g_lin.SetLineWidth(2);
  for(size_t i=0; i<nprof; i++)
      g_lin.SetPoint(g_lin.GetN(),cl.prof_s.at(i),cl.prof_lin.at(i));
  
  TGraph  g_bnd;
  g_bnd.SetMarkerColor(kRed);
  g_bnd.SetMarkerStyle(29);
  g_bnd.SetMarkerSize(2.0);
  //if( cl.bnd_s >= 0 ) g_bnd.SetPoint(0, cl.prof_W.at(cl.bnd_i), cl.prof_X.at(cl.bnd_i) );
  if( cl.maxQ_i >= 0 ) g_bnd.SetPoint(0, cl.prof_W.at(cl.maxQ_i), cl.prof_X.at(cl.maxQ_i) );
        
  sprintf(buffer, "clstr%i_r%i_sr%i_e%i_pl%d_%3.0f_MeV", fNumSavedHitGraphs, fRunNumber, fSubRunNumber, fEventNumber, plane, elShowerEnergy);
  if( fIsMC && fTrue_ElShowerEnergyDep > 0. ) 
    sprintf(buffer, "clstr%i_r%i_sr%i_e%i_pl%d_Reco-%3.0fMeV_True-%3.0fMeV", fNumSavedHitGraphs, fRunNumber, fSubRunNumber, fEventNumber, plane, elShowerEnergy, fTrue_ElShowerEnergyDep);
  fTCanvas     = graphDir_hits.make<TCanvas>(buffer,"c2",700,700);
  TPad c1("c1","c1",0.0,0.60,1.0,1.0);
  TPad c2("c2","c2",0.0,0.0,1.0,0.60);
  c1.Draw();
  c2.Draw();
  c2.Divide(1,2);
        
  c1.cd();
  TMultiGraph g;
  if(g_shower.GetN() > 0 )  g.Add(&g_shower,"P");
  if(g_hits.GetN() > 0 )    g.Add(&g_hits,"AP");
  if(g_cluster.GetN() > 0)  g.Add(&g_cluster,"PL");
  if(g_seed.GetN() > 0)  g.Add(&g_seed,"P");
  if(g_bnd.GetN() > 0 )     g.Add(&g_bnd,"P");
  g.Draw("a");
  g.SetTitle(buffer);
  g.GetXaxis()->SetTitle("Wire coordinate [cm]");
  g.GetYaxis()->SetTitle("Drift coordinate [cm]");
  g.GetXaxis()->SetTitleSize(0.05);
  g.GetYaxis()->SetTitleSize(0.05);
  g.GetXaxis()->SetLabelSize(0.05);
  g.GetYaxis()->SetLabelSize(0.05);
  g.GetYaxis()->SetTitleOffset(0.8);
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);
        
  double ts = 0.07;
  double ls = 0.06;

  c2.cd(1);
  TMultiGraph g2;
  g2.Add(&g_dQ,"APL");
  g2.Add(&g_dQt,"APL");
  g2.Draw("a");
  g2.GetYaxis()->SetTitle("Charge density #it{q/ds} [e-/cm]");
  g2.GetYaxis()->CenterTitle();
  g2.GetXaxis()->SetTitleSize(ts);
  g2.GetXaxis()->SetLabelSize(ls);
  g2.GetYaxis()->SetTitleSize(ts);
  g2.GetYaxis()->SetLabelSize(ls);
  g2.GetYaxis()->SetTitleOffset(0.7);
  gPad->SetBottomMargin(0.01);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->Update();
  TLine maxQLineA(cl.maxQ_s,gPad->GetUymin(),cl.maxQ_s,gPad->GetUymax());
  maxQLineA.SetLineColor(kRed);
  maxQLineA.SetLineWidth(1);
  maxQLineA.SetLineStyle(2);
  TLine bndLineA(cl.bnd_s,gPad->GetUymin(),cl.bnd_s,gPad->GetUymax());
  bndLineA.SetLineColor(kRed);
  bndLineA.SetLineWidth(2);
  bndLineA.SetLineStyle(7);
  if( cl.maxQ_s == cl.bnd_s ) maxQLineA.SetLineWidth(0);
  if( elShowerEnergy >= 0. ) bndLineA.SetLineStyle(1);
  maxQLineA.Draw();
  bndLineA.Draw();

  c2.cd(2);
  TMultiGraph g3;
  g3.Add(&g_lin,"APL");
  g3.Draw("a");
  g3.GetXaxis()->SetTitle("Projected 2D Distance [cm]");
  g3.GetYaxis()->SetTitle("Local Linearity #chi^{2}");
  g3.GetYaxis()->CenterTitle();
  g3.GetXaxis()->SetTitleSize(ts);
  g3.GetXaxis()->SetLabelSize(ls);
  g3.GetYaxis()->SetTitleSize(ts);
  g3.GetYaxis()->SetLabelSize(ls);
  g3.GetYaxis()->SetTitleOffset(0.7);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.15);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->Update();
  TLine maxQLineB(cl.maxQ_s,gPad->GetUymin(),cl.maxQ_s,gPad->GetUymax());
  maxQLineB.SetLineColor(kRed);
  maxQLineB.SetLineWidth(1);
  maxQLineB.SetLineStyle(2);
  TLine bndLineB(cl.bnd_s,gPad->GetUymin(),cl.bnd_s,gPad->GetUymax());
  bndLineB.SetLineColor(kRed);
  bndLineB.SetLineWidth(2);
  bndLineB.SetLineStyle(7);
  if( cl.maxQ_s == cl.bnd_s) maxQLineB.SetLineWidth(0);
  if( elShowerEnergy >= 0. ) bndLineB.SetLineStyle(1);
  maxQLineB.Draw();
  bndLineB.Draw();
  
  fTCanvas->Update();
  fTCanvas->Write(buffer);
  
  fNumSavedHitGraphs++;
   
}


DEFINE_ART_MODULE(MichelAnaFilter)
