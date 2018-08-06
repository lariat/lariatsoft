////////////////////////////////////////////////////////////////////////
/// Class:       MichelAna
/// Module Type: filter
/// File:        MichelAna_module.cc
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
#include "lardata/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcorealg/CoreUtils/ProviderUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
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
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

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

// MichelCluster: special data structure to handle clustering/shower information
struct MichelCluster {
  int               plane;
  int               seedHit;
  std::vector<int>  cluster;
  std::vector<int>  cluster_el;
  std::vector<int>  cluster_mu;
  std::vector<int>  shower;
  std::vector<bool> isInTrk;
  float             minCovAtBnd;
  float             truncMeanSlope;
  int               muEndHit;
  float             elShowerFrac;
  TVector3          muEnd2D;
  TVector3          muEndDir2D;
  float             muEnd2D_X;
  float             muEnd2D_W;
  float             elDir2D_X;
  float             elDir2D_W;
  float             decayAngle2D;
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
  float muCharge;
  float elShowerCharge;
  float elShowerEnergy;
  float elTrackCharge;
  float elTrackEnergy;
  float fracMuHitsLinear;
  float muAveLinearity;
  int   nPtsMuFit;
  float aveDriftTime;
  float aveX;
  
  // initialize all quantities in the
  // creation of new cluster object  
  MichelCluster() {
    plane = -9;
    seedHit = -9;
    muEndHit = -9;
    muEnd2D_X = -9.;
    muEnd2D_W = -9.;
    elDir2D_X = 0.;
    elDir2D_W = 0.;
    elShowerFrac = -9.;
    decayAngle2D = -9;
    muEnd2D.SetXYZ(-99., -99., -99.);
    muEndDir2D.SetXYZ(-99., -99., -99.);
    minCovAtBnd = 9.;
    truncMeanSlope = 9.;
    bnd_i = -9;
    maxQ_i = -9;
    bnd_s = -999.;
    maxQ_s = -999.;
    extraHits = 0;
    fracMuHitsLinear = -9.;
    nPtsMuFit = -9;
    muAveLinearity = -9.;
    totalCharge = -999.;
    muCharge       = -999.;
    elShowerCharge = -999.;
    elShowerEnergy = -999.;
    elTrackCharge  = -999.;
    elTrackEnergy = -999.;
    aveDriftTime  = -99.;
    aveX          = -99.;
  }

};

//########################################################################################
class MichelAna;


//########################################################################################
class MichelAna : public art::EDFilter 
{

  public:
  
  explicit MichelAna(fhicl::ParameterSet const & p);
  MichelAna(MichelAna const &) = delete;
  MichelAna(MichelAna &&) = delete;
  MichelAna & operator = (MichelAna const &) = delete;
  MichelAna & operator = (MichelAna &&) = delete;

  // Art functions
  bool  filter(art::Event & e) override;
  void  beginJob() override;
  void  endJob() override;
  void  reconfigure(fhicl::ParameterSet const & p);

  // Custom functions   
  void  ResetVariables();
  void  GetDetProperties();
  void  GetTruthInfo( art::Event &);
  void  ParticleTracker(const art::Event&, int, int,TH1D* h[], TH1D* hm[], TH1D* he[]);
  bool  IsParticleDescendedFrom(const art::Handle< std::vector<simb::MCParticle>>&,int,int);
  float CalcParticleHeritageDisplacement(const art::Handle< std::vector<simb::MCParticle>>&,const simb::MCParticle&, int );
  bool  IsPointInFiducialVolume(TVector3, float, float, float);
  bool  IsPointInFiducialVolume(TVector3, float, float, float, float, float, float);
  bool  IsParticleContained(const simb::MCParticle& );
  float CalcParticleDist(const simb::MCParticle& );
  void  PerfectWfmReco(raw::OpDetPulse & pulse);
  void  PropagatePhoton(float,float,float,float,float,bool,bool,TH1D*,TH1D*,TH1D*);
  void  ReadPhotonLibrary(TVector3&,int,int,float&,float&,float&,float&);
  float CalcSigma(size_t,float,float);
  void  MuContamCorr(int);
  float TauConversion(float);
  float CorrectForQuenching(float, float, float);
  float PeIntegrationCorrection(float, float, float, float);
  void  FindBestMatchedHit(int,float,int&,float&);
  void  FindBestMatchedHit(int,float,int&,float&,std::vector<int>&);
  void  Clustering(int,MichelCluster&);
  void  Clustering(int,MichelCluster&,float);
  void  Clustering(int,MichelCluster&,float,std::vector<int>&,float);
  void  ProximityCluster(std::vector<float>&,std::vector<float>&,std::vector<int>&,std::vector<int>&,int,int,float);
//  void  SuperSonicCluster(std::vector<float>&,std::vector<float>&,std::vector<int>&,std::vector<int>&,int,int,float);
  void  CalcMichelShowerEnergy(MichelCluster&);
  float CalcTruncatedMeanInProfile(std::vector<float>&,size_t,int,float);
  float CalcLocalLinearity(std::vector<float>&,std::vector<float>&,size_t,int);
  float GetVisibility(TVector3,int);
  float GetVisibility(TVector3,int,int);
  int   GetGlobalBin(const TH3D*,double,double,double);
  int   GetGlobalBin(const TH3D*,const TLorentzVector&);
  int   GetGlobalBin(const TH3D*,const TVector3&);
  float Integrate(const TH1D*,float,float); 
  float Integrate(const std::vector<float>&,int,int); 
  void  AddToAverageWfm(std::vector<float>&,short,short,TH1D*,int&,int);
  void  AddToAverageWfm(const TH1D*,short,short,short,TH1D*,int &,int);
  void  MakeClusteringGraphs(MichelCluster&);
  void  Display3DEvent(std::vector<TVector3>&, std::vector<TVector3>&); 
  void  MakeWfmGraphs(art::Event&);
  float T_from_X(float,int);
  
  float ElectronsFromADCArea( float, int ); 
  float dQdx_from_dEdx(float);

private:
 
  // ................................................................................
  // Algs and service objects
  geo::GeometryCore     *fGeo;
  OpHitBuilderAlg       fOpHitBuilderAlg; 
  TTree*                fTree;
  TTree*                fTreeCrsMu;

  calo::CalorimetryAlg  fCaloAlg; // CAN'T GET THIS TO WORK?!?
    
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
  float               fEdgeMarginX;         // fiducial margin X [cm]
  float               fEdgeMarginY;         // fidicual margin Y [cm]
  float               fEdgeMarginZ;         // fiducial margin Z [cm]
  float               fMuTrackLengthMin;    // min mu candidate track length [cm]
  float               fMuTrackZenithAngleMax; // max mu candidate zenith angle [deg]
  float               fMuResRangeOffset;    // offset on residual range [cm] 
  
  std::vector<bool>   fUsePmtForCalo;
  bool                fMakeAveWfms;         // fill average waveforms 
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
  
  float               fElectronLifetimeErr; // error on lifetime for MC
  std::vector<float>  fSmearSigT;           // integration time-dependent smearing on MC
  std::vector<float>  fSmearSigPE;          // photoelectron resolution smearing on MC
  std::vector<float>  fSmearFactor;         // final smear factor on MC total light
  std::vector<float>  fSmearFactorPrompt;   // final smear factor on MC prompt light (NOT USED)
  std::vector<float>  fCollectionEff;       // PMT coll. efficiency for use in MC
  std::vector<float>  fQE_ScaleFactor;      // global QE scale factor (from optimization)
  std::vector<float>  fQE_ScaleFactor_Vis;  // QE scale factor for vis. light (NOT USED)
  std::vector<float>  fQE_ScaleFactor_VUV;  // QE scale factor for VUV light (NOT USED)
  float               fSingletToTripletRatio;      // singlet-to-triplet ratio used in MC
  float               fLateLightTau;        // eff. late-light tau in MC (determines quenching)
  float               fTauT;                // "true" triplet dimer tau (influences quenching)
  float               fTauS;                // "true" singlet dimer tau
  bool                fApplyTrigEffCut;     // cut out MC events based on trigger eff. funct.
  std::vector<float>  fTrigEff_P;           // trigger eff. funct. param "P" (from optimization)
  std::vector<float>  fTrigEff_K;           // trigger eff. funct. param "K" (fixed to 8)
 
  bool                fMuContamCorr;        // correct for muon late-light contamination
  float               fMuContamCorr_dTmin;  // min dT required for mu contam. correction
  std::vector<float>  fMuContamCorr_EffTau; // eff. late-light tau for mu contam corr. (from DB!)
  bool                fCorrectQuenching;    // apply quenching corrections to late light

  bool                fUseMaxDropAsBnd;
  float               fMaxHitSeparation;    // max 2D dist (in WX space) btwn adjacent hits
  int                 fTruncMeanWindow;     // size local nbd. for trunc. mean in profile (+/-N)
  float               fTruncMeanP;          // frac of upper/lower tail to truncate
  int                 fMaxMatchDistance;    // dist (in hits) to find max Q/ds after trunc mean peak
  int                 fLocalLinearityWindow;// local nbd. in profile to calc local linearity (+/-N)
  float               fLinThresh;           // max local linearity allowed at mu-el boundary
  float               fLinTol;              // adjusts thresh based on RMS of mu track linearity
  float               fLinThreshMuFit;      // min local linearity for "good" mu track hits
  float               fShowerAcceptanceAngle;// opening angle of 2D electron shower cone
  std::vector<float>  fMinHitAmp;           // min allowable hit amplitude to include in shower
  std::vector<float>  fMinHitCharge;        // min allowable hit charge to include in shower
  
//  std::vector<float>  fCalAreaConstants;    // scale hit integrals to charge (ADC -> e-)
  float               fWph;                 // energy to excitate/ionize Ar atom (19.5 eV)
  float               fWion;                // energy to produce e-ion pair (23.6 eV)
  float               fExcRatio;            // excitation ratio ( N_ex / N_i ) = 0.21
  bool                fUseUniformRecomb;
  float               fRecombFactor;        // recomb. survival frac for electron trk/shwr hits
  float               fRecombFactorTrack;   // recomb. survival frac for track-like el. hits
  float               fRecombFactorShower;  // recomb. survival frac for shwr-like el. hits
//  float               fChargeSmearFactor;

  float               fdTcut;               // min decay time for energy histograms 
  std::vector<float>  fPromptPECut;         // min prompt PE for dT histogram 
  int                 fMinClusterSize;      // min allowable cluster size for histos
  float               fMinElShowerFrac;     // min electron shower completeness frac for histos
  float               fMinFracMuHitsLinear; // min frac of good linearity mu trk hits for histos
  int                 fMinMuClusterSize;    // min mu cluster size for histos
  int                 fMinMuClusterHitsEndFit;  // min mu hits used in direction fit for histos
  float               fMinTruncMeanSlope;
  int                 fMinElClusterSize;    // min el cluster size for histos
  int                 fMaxElClusterSize;    // max el cluster size for histos
  float               fMinMuLinearity;      // min average linearity of mu hits for histos
  int                 fMaxExtraHits;        // max num. of unclustered hits near mu end for histos
  float               fMaxDecayAngle2D;     // max 2D decay angle for histos
  float               fMinDecayAngle2D;     // min 2D decay angle for histos
  float               fMaxHitTimeDiffFac;   // max allowable time diff to match hits btwn planes (xhitRMS)
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
  float               fElectronLifetimeFromDB;
  float               fElectronLifetime;
  bool                fMichelOpticalID;
  bool                fSingleOpticalHit;
  
  // Track information
  int                 fNumTrack;
  int                 fNumTrackStopping;
  int                 fNumTrackCrossing;
  int                 fNumTrackContained;
  std::vector<int>    fMuTrackHits;
  std::vector<int>    fBgTrackHits;
  int                 fMuTrackIndex;
  int                 fMuTrackID;
  float               fMuTrackLength;
  float               fMuTrackEnergy;
  float               fMuTrackZenithAngle;
  float               fMuTrackInclineAngle;
  TVector3            fMuTrackEndDir;
  float               fMuTrackVertex_X;
  float               fMuTrackVertex_Y;
  float               fMuTrackVertex_Z;
  float               fMuTrackEnd_X;
  float               fMuTrackEnd_Y;
  float               fMuTrackEnd_Z;
  std::vector<TVector3> fvMuTrkXYZ[2];
  std::vector<double> fvMuTrkdEdx[2];
  std::vector<double> fvMuTrkdQdx[2];
  std::vector<double> fvMuTrkdADCdx[2];
  std::vector<double> fvMuTrkResRange[2];
  std::vector<double> fvMuTrkPitch[2];
//  std::vector<double> fvMuClusterdEdx;
//  std::vector<double> fvMuClusterdQdx;
//  std::vector<double> fvMuClusterdADCdx;
//  std::vector<double> fvMuClusterResRange;
//  std::vector<double> fvMuClusterdEdx_Pl0;
//  std::vector<double> fvMuClusterdQdx_Pl0;
//  std::vector<double> fvMuClusterdADCdx_Pl0;
//  std::vector<double> fvMuClusterResRange_Pl0;
  bool                fMuTrackIsCaloOrdered;
  TVector3            fMuTrackVertex;
  TVector3            fMuTrackEnd;
//  float               fMuTrackPitch[2];
  
  // Crossing muon information (calibration)
  TVector3            fCrsMuVertex;
  TVector3            fCrsMuEnd;
  int                 fCrsMuTrackIndex;
  float               fCrsMuVertex_X;
  float               fCrsMuVertex_Y;
  float               fCrsMuVertex_Z;
  float               fCrsMuEnd_X;
  float               fCrsMuEnd_Y;
  float               fCrsMuEnd_Z;
  float               fCrsMuEnergy;
  float               fCrsMuLength;
  float               fCrsMuPhotons;
  float               fCrsMuPhotonsPrompt;
  float               fCrsMuCharge;
  float               fCrsMuCharge_Pl0;
  float               fCrsMuIntegral;
  float               fCrsMuIntegral_Pl0;
  float               fCrsMuPitch[2];
  float               fCrsMuHitFrac;
  float               fCrsMuPE_total[2];
  float               fCrsMuPE_prompt[2];
  float               fTrue_CrsMuCharge;
  float               fTrue_CrsMuPhotons;

  // Optical information
  float               fTriggerTime;
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
  std::vector<float>  fvHitADC_total[2];

  // Prepulse fit hit integrals, hit amplitude, width
  std::vector<float>  fvHitAmplitudes[2];
  std::vector<float>  fvHitWidth[2];
  std::vector<float>  fvHitADCpf_100ns[2];
  std::vector<float>  fvHitADCpf_total[2];
  std::vector<short>  fvPrepulseX1[2];
  std::vector<float>  fvPrepulseZeroPoint[2];
  std::vector<float>  fvPrepulseBaseline[2];
  std::vector<float>  fvPrepulseRMS[2];
  std::vector<float>  fvPrepulseSlowNorm[2];
  std::vector<float>  fvPrepulseSlowTau[2];

  // Derived quantities when dT is defined (2 hits)
  float               fMuContam_prompt[2];
  float               fMuContam_total[2];
  float               fPE_promptRaw[2];
  float               fPE_totalRaw[2];
  float               fPE_prompt[2];
  float               fPE_total[2];
  float               fPE_total_qc[2];
  float               fPromptFracRaw[2];
  float               fPromptFrac[2];
  float               fMuAmplitude[2];
  float               fMuPE_prompt[2];
  float               fMuPE_total[2];
  bool                fMuPulseSaturated[2];
  bool                fElPulseSaturated[2];
  float               fPE_comb;
  float               fPE_comb_qc;
  float               fPE_prompt_comb;
  float               fDecayTime;
    
  // 2D clustering data
  int                 fNumPlaneHits[2];
  float               fTotalCharge[2];
  float               fTotalIntegral[2];
  // ...collection plane
  float               fFracDropoffAtBnd;
  float               fTruncMeanSlope;
  float               fCovAtBnd;
  int                 fClusterSize;
  int                 fMuClusterSize;
  int                 fElClusterSize;
  int                 fExtraHits;
  float               fMuAveLinearity;
  float               fFracMuHitsLinear;
  int                 fMuClusterHitsEndFit;
  float               fDecayAngle2D;
  float               fMuEnd2D_W;
  float               fMuEnd2D_X;
  float               fElDir2D_W;
  float               fElDir2D_X;
  int                 fElShowerSize;
  float               fElShowerFrac; 
  float               fElShowerCharge;
  float               fElTrackCharge;
  float               fElShowerEnergy;
  float               fElTrackEnergy;                 
  float               fMuCharge;
  // ...induction plane
  float               fFracDropoffAtBnd_Pl0;
  float               fTruncMeanSlope_Pl0;
  float               fCovAtBnd_Pl0;
  int                 fClusterSize_Pl0;
  int                 fMuClusterSize_Pl0;
  int                 fElClusterSize_Pl0;
  int                 fExtraHits_Pl0;
  float               fMuAveLinearity_Pl0;
  float               fFracMuHitsLinear_Pl0;
  int                 fMuClusterHitsEndFit_Pl0;
  float               fDecayAngle2D_Pl0;
  float               fMuEnd2D_W_Pl0;
  float               fMuEnd2D_X_Pl0;
  float               fElDir2D_W_Pl0;
  float               fElDir2D_X_Pl0;
  int                 fElShowerSize_Pl0;
  float               fElShowerFrac_Pl0; 
  float               fElShowerCharge_Pl0;
  float               fElTrackCharge_Pl0;
  float               fElShowerEnergy_Pl0;
  float               fElTrackEnergy_Pl0;                 
  
  float               fProjDist2DToWires;
  
  // 3D shower information
  float               fFracHits3D;
  int                 fNumPts3D;
  float               fElShowerVis;
  float               fElShowerVisCh[2];
  TVector3            fElShowerCentroid;
  float               fElShowerCentroid_X;
  float               fElShowerCentroid_Y;
  float               fElShowerCentroid_Z;
  float               fElShowerPhotons;
  float               fElShowerPhotonsPrompt;
  float               fElShowerEnergyQL;
  TVector3            fMuEnd3D;
  float               fMuEnd3D_X;
  float               fMuEnd3D_Y;
  float               fMuEnd3D_Z;
  float               fMuEndHitTimeDiff;
  TVector3            fElDir3D;
  float               fElDir3D_X;
  float               fElDir3D_Y;
  float               fElDir3D_Z;
  float               fProjDist3DToWires;
  float               fProjPtOnWires_Y;
  float               fProjPtOnWires_Z;

  // Truth information
  bool                fTrue_IsMuStopping;
  int                 fTrue_MuSign;
  float               fTrue_MuEnergy;
  float               fTrue_MuStopTime;
  float               fTrue_MuTrackLength;
  TVector3            fTrue_MuTrackVertex;
  TVector3            fTrue_MuTrackEnd;
  float               fTrue_MuTrackVertex_X;
  float               fTrue_MuTrackVertex_Y;
  float               fTrue_MuTrackVertex_Z;
  float               fTrue_MuTrackEnd_X;
  float               fTrue_MuTrackEnd_Y;
  float               fTrue_MuTrackEnd_Z;
  float               fTrue_MuEnergyDep;
  float               fTrue_MuCharge;
  float               fTrue_MuPhotons;
  int                 fTrue_NumBremmPhotons;
  float               fTrue_ElEnergy;
  float               fTrue_ElShowerCentroid_X;
  float               fTrue_ElShowerCentroid_Y;
  float               fTrue_ElShowerCentroid_Z;
  float               fTrue_ElShowerPhotons;
  float               fTrue_ElTrackPhotons;
  float               fTrue_ElShowerPhotonsPrompt;
  float               fTrue_ElShowerPhotonsLate;
  float               fTrue_ElTrackEnergyDep;
  float               fTrue_ElTrackChargeDep;
  float               fTrue_ElShowerEnergyDep;
  float               fTrue_ElShowerChargeDep;
  float               fTrue_ElShowerTotalTrajLength;
  float               fTrue_ElShowerIons;
  float               fTrue_ElShowerVis;
  float               fTrue_ElShowerVisCh[2];
  TVector3            fTrue_ElMomentum;
  float               fTrue_ElMomentum_X;
  float               fTrue_ElMomentum_Y;
  float               fTrue_ElMomentum_Z;
  float               fTrue_ProjDist3DToWires;
  float               fTrue_ElAngle;
  bool                fTrue_IsElTrkContained;
  bool                fTrue_IsElShwrContained;
  float               fTrue_ContainmentFrac;
  float               fTrue_TotalChargeDep;
  float               fTrue_TotalEnergyDep;
  float               fTrue_dT;
  float               fTrue_Amplitude[2];
  float               fTrue_MuPE_prompt[2];
  float               fTrue_MuPE_total[2];
  float               fTrue_MuContam_prompt[2];
  float               fTrue_MuContam_total[2];
  float               fTrue_MuPE;
  float               fTrue_PE0;
  float               fTrue_PE;
  float               fTrue_PEPrompt;
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
  TH1D*               hContamFrac[2];
  TH1D*               hContamFrac_dTcut[2];
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
  TH1D*               hMuResRangePitch[2];
//  TH2D*               hMuResRangeVsdEdxHyp[2];
  TH2D*               hMuResRangeVsdEdx[2];
  TH2D*               hMuResRangeVsdQdx[2];
  TH2D*               hMuResRangeVsdADCdx[2];
//  TH2D*               hMuClusterResRangeVsdEdx[2];
//  TH2D*               hMuClusterResRangeVsdADCdx[2];
//  TH2D*               hMuClusterResRangeVsdQdx[2];
  //TH2D*               hMu_dEdxHyp_vs_dQdx;
  TH2D*               hTrue_MuResRangeVsdEdx;
  TH2D*               hTrue_MuResRangeVsdEdx_HR;
  TH1D*               hExtraHits;
  TH1D*               hNumHits[2];
  TH1D*               hHitRMS[2];
  TH1D*               hHitRMS_ElShower[2];
  TH1D*               hHitT[2];
  TH1D*               hHitX[2];
  TH1D*               hHitAmp[2];
  TH1D*               hHitAmp_ElShower[2];
  TH2D*               hHitAmp_vs_HitRMS[2];
  TH2D*               hHitAmp_vs_HitIntegral[2];
  TH2D*               hHitIntegral_vs_HitRMS[2];
  TH1D*               hHitCharge[2];
  TH1D*               hHitCharge_ElShower[2];
  TH1D*               hHitIntegral[2];
  TH1D*               hHitIntegral_ElShower[2];
  TH1D*               hTotalCharge;
  TH1D*               hElTrackFrac;
  TH1D*               hProjDist3DToWires;
  TH1D*               hProjDist3DToWiresRes;
  TH2D*               hProjPtOnWires;
  TH1D*               hClusterHitSeparation;
  TH1D*               hFracDropoffAtBnd;
  TH1D*               hElShowerVis;
  TH1D*               hElShowerCharge;
  TH1D*               hElShowerCharge3D;
  TH1D*               hElTrackCharge;
  TH1D*               hElShowerPhotons;
  TH1D*               hElShowerEnergy;
  TH1D*               hElShowerEnergyQL;
  TH2D*               hElShowerEnergyVsPromptPE[2];
//  TH1D*               hElTrackdEdx;
  TH1D*               hQoverL;
  TH1D*               hQoverLCrsMuTrk;
  TH1D*               hTruncMeanSlope; 
  TH1D*               hTruncMeanSlopeChi2; 
  TH1D*               hClusterSize;
  TH1D*               hCovAtBnd;
  TH2D*               hElHitChargeRatio;
  TH1D*               hElShowerFrac;
  TH1D*               hLeftoverHits;
  TH1D*               hDistMaxTdQdsToMaxQ;
  TH1D*               hDistMaxTdQdsToMaxDrop;
  TH1D*               hFracMuHitsLinear;
  TH1D*               hMuAveLinearity;
  TH1D*               hMuClusterSize;
  TH1D*               hElClusterSize;
  TH1D*               hElShowerSize;
  TH1D*               hMuClusterHitsEndFit; 
  TH1D*               hDecayAngle2D;
  TH1D*               hRecomb;
  TH1D*               hRecombCrsMuTrk;
  TH1D*               hElShowerDepDist;
  TH1D*               hElShowerDepAngle2D;
  TH1D*               hNumPts3D;
  TH1D*               hFracHits3D;
  TH2D*               hElClusterSizeCompare;
  TH2D*               hShowerSizeCompare;
  TH1D*               hHitTimeDiff;
  TH1D*               hHitXDiff;
  TH1D*               hMuEndHitTimeDiff;
  TH1D*               hMuEndDiff3D;
  TH1D*               hMuEnd3D_X;
  TH1D*               hMuEnd3D_Y;
  TH1D*               hMuEnd3D_Z;
  TH2D*               hMuEnd3D_ZX;
  TH2D*               hMuEnd3D_ZY;
  TH1D*               hElectronLifetimeFromDB;
  TH1D*               hElectronLifetime;
  TH1D*               hEffTau[2];
  TH1D*               hElShowerCentroid_X;
  TH1D*               hElShowerCentroid_Y;
  TH1D*               hElShowerCentroid_Z;
  TH2D*               hElShowerCentroid_ZX;
  TH2D*               hElShowerCentroid_ZY;
 
  TH1D*               hMuTrackZenithAngle;
  TH1D*               hMuTrackInclineAngle;
  TH1D*               hMuTrackInclineAngle_dEdx;
  
  TH1D*               hTrue_ResQ;
  TH1D*               hTrue_ResQ_Mu;
  TH1D*               hTrue_ResE_CrsMu;
  TH1D*               hTrue_ResQ_CrsMu[2];
  TH1D*               hTrue_ResQ_CrsMuTrk[2];
  TH1D*               hTrue_ResQ_Shwr3D;
  TH1D*               hTrue_ResL;
  TH1D*               hTrue_ResL_CrsMu;
  TH1D*               hTrue_ResPE;
  TH1D*               hTrue_ResPEPrompt;
  TH1D*               hTrue_ResPE0;
  TH1D*               hTrue_ResVis;
  TH1D*               hTrue_EnergyRes;
  TH1D*               hTrue_EnergyRes_Shwr3D;
  TH1D*               hTrue_EnergyResTrk;
  TH1D*               hTrue_EnergyResTrk_Shwr3D;
  TH1D*               hTrue_EnergyResQL;
  TH2D*               hTrue_EnergyVsEnergyRes;
  TH1D*               hTrue_ResPEch[2];

  TH1D*               hTrue_DriftStochasticity;
 
  TH1D*               hTrue_MuEnergy; 
  TH1D*               hTrue_MuAveMomentum;
  TH1D*               hTrue_ElEnergy;
  TH1D*               hTrue_ElEnergyFree;
  TH1D*               hTrue_ElEnergyCapture;
  TH1D*               hTrue_ElTrackEnergyDep;
  TH1D*               hTrue_ElShowerEnergyDep; 
  TH1D*               hTrue_ElShowerChargeDep;
  TH1D*               hTrue_ElTrackChargeDep;
  TH1D*               hTrue_ElShowerPhotons;
  TH1D*               hTrue_PE_preTrig[2];
  TH1D*               hTrue_PE_prompt[2];
  TH1D*               hTrue_PE_total[2];
  TH1D*               hTrue_PE_total_dTcut_shwr[2];
  TH2D*               hTrue_PE_total_TrueVsReco[2];
  TH2D*               hTrue_MuLateLightContamination;
  TH2D*               hTrue_dT_vs_MuContamRes[2];
  TH2D*               hTrue_MuTrackEnd_ZX;
  TH2D*               hTrue_MuTrackEnd_ZY;
  TH1D*               hTrue_MuTrkVtxRes_X;
  TH1D*               hTrue_MuTrkVtxRes_Y;
  TH1D*               hTrue_MuTrkVtxRes_Z;
  TH1D*               hTrue_MuClsEndRes;
  TH1D*               hTrue_MuClsEndRes_X;
  TH1D*               hTrue_MuClsEndRes_Y;
  TH1D*               hTrue_MuClsEndRes_Z;
  TH1D*               hTrue_MuClsEndResDir;
  TH1D*               hTrue_MuTrkEndRes;
  TH1D*               hTrue_MuTrkEndRes_X;
  TH1D*               hTrue_MuTrkEndRes_Y;
  TH1D*               hTrue_MuTrkEndRes_Z;
  TH1D*               hTrue_MuTrkEndResDir;
  TH1D*               hTrue_MuTrkEndResTrans;
  TH1D*               hTrue_MuTrkEndResDir_cut;
  TH1D*               hTrue_MuTrkEndResTrans_cut;
  TH2D*               hTrue_EnergyDepVsRecoEnergy;
  TH1D*               hTrue_PhotonRes;
  TH1D*               hTrue_dEdx;
  TH1D*               hTrue_dEdx_ElTrk;
  TH1D*               hTrue_dEdx_ElShw;
  TH1D*               hTrue_dQdx;
  TH1D*               hTrue_dEdx_CrsMu;
  TH1D*               hTrue_dEdx_CrsMu_4mm;
  TH1D*               hTrue_dEdx_CrsMu_1cm;
  TH1D*               hTrue_dEdx_CrsMu_2cm;
  TH1D*               hTrue_dEdx_CrsMu_5cm;
  TH1D*               hTrue_dEdx_CrsMu_10cm;
  TH1D*               hTrue_dQdx_CrsMu;
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
  TH1D*               hTrue_IDE_X;
  TH1D*               hTrue_IDE_Y;
  TH1D*               hTrue_IDE_Z;
  TH1D*               hTrue_IDE_NumElectrons;
  TH1D*               hTrue_LightYield[2];
  TH1D*               hTrue_PartDispFromElectron;
  TH1D*               hTrue_BremmPhotonLength;
  TH2D*               hTrue_ElEnergyVsNumBremmPhotons;
  TH1D*               hNumTrack;
  TH1D*               hNumTrackStopping;
  TH1D*               hNumTrackContained;
  TH1D*               hNumTrackCrossing;
  TH1D*               hTrackNode_X;
  TH1D*               hTrackNode_X_MuTrk;
  TH1D*               hTrackNode_Y;
  TH1D*               hTrackNode_Z;
  TH2D*               hTrackNode_ZX;
  TH2D*               hTrackNode_ZY;
  TH1D*               hTrackZenithAngle;
  TH1D*               hTrackLength;
  TH2D*               hTrackLength_vs_ZenithAngle;
 
  TH1D*               hCrsMuLength; 
  TH1D*               hCrsMuHitFrac; 
  TH2D*               hCrsTrackNode_ZX;
  TH2D*               hCrsTrackNode_ZY;
  TH1D*               hCrsTrackdT_ACP;
  TH1D*               hCrsMuPitch[2];
  TH1D*               hCrsMuADCPerCm;
  TH1D*               hCrsMuADCPerCm_Pl0;
  TH1D*               hCrsMuQPerCm;
  TH1D*               hCrsMuQPerCm_Pl0;
  TH1D*               hTrue_CrsMuQPerCm;
  TH1D*               hTrue_CrsMuLPerCm;
  TH1D*               hTrue_CrsMuLength;
  TH1D*               hCrsMuLPerCm;
  TH1D*               hCrsMu_dADCdx[2];
  TH1D*               hCrsMu_dQdx[2];
  TH1D*               hCrsMu_dEdx[2];
  
//  TH2D*               hHitIntegral_vs_SummedADC_ElShower;
//  TH2D*               hHitIntegral_vs_SummedADC_MuTrk;
//  TH2D*               hHitIntegral_vs_DiffSummedADC_ElShower;
//  TH2D*               hHitIntegral_vs_DiffSummedADC_ElTrk;
//  TH2D*               hHitIntegral_vs_DiffSummedADC_MuTrk;
//  TH2D*               hHitAmp_vs_DiffSummedADC_ElShower;
//  TH2D*               hHitAmp_vs_DiffSummedADC_ElTrk;
//  TH2D*               hHitAmp_vs_DiffSummedADC_MuTrk;
//  TH2D*               hTrackInclineAngle_vs_DiffSummedADC;
  
  /* 
  TH1D* h_tmp;
  TH1D* h_tmp_n;
  TH1D* h_tmp_HR;
  TH1D* h_tmp_HRn;
  */
  
  // Angle bins
  float               fAngleBins[4];
  float               fAngleBins_X1[4];
  TH2D*               hMuResRangeVsdEdx_AngleBins[4];
  TH2D*               hMuResRangeVsdQdx_AngleBins[4];
  TH2D*               hMuResRangeVsdADCdx_AngleBins[4];

  // ...........................................................
  // Average waveforms (muon, electron pulse)
  TH1D*               hAveWfm[2];
  TH1D*               hAveWfm_mu[2];
  TH1D*               hAveWfm_el[2];
  TH1D*               hAveWfm_crsmu[2];
  int                 aveWfmCounter[2];
  int                 aveWfmCounter_mu[2];
  int                 aveWfmCounter_el[2];
  int                 aveWfmCounter_crsmu[2];
  TH1D*               hAvePhelProfile_mu[2];
  int                 avePhelProfileCounter_mu[2];

  // ...........................................................
  // Constant, buffers, counters, misc... 
  std::vector<raw::OpDetPulse> fPulses;
  std::vector<short>    fIntegrationWindows;
  int                 fCachedRunNumber;
  char	              buffer[200], histName[200], histTitle[200];
  float               fEfield;
  float               fDriftVelocity[3]; 
  float               fSamplingRate; 
  float               fXTicksOffset[2]; 
  float               fTriggerOffset;
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
  int                 fNumSaved3DGraphs;
  int                 fNumSavedWfmGraphs;
  TCanvas*            fTCanvas;
  TRandom2*           fRand;
  std::default_random_engine generator;
  TH1D*               hPMT_phelTimes[2];
  TH1D*               hPMT_phelTimes_electron[2];
  TH1D*               hPMT_phelTimes_muon[2];
  std::vector<float>  fPMT_wfm[2];
  std::vector<short>  fPMT_wfm_raw[2];
  std::vector<float>  fvbs[2];

  // ..........................................................
  // Vectors to store optical calibration information to be read
  // from a file at the start of the job
  std::vector<int>    fCalibRunStart[2];
  std::vector<int>    fCalibRunEnd[2];
  std::vector<float>  fCalibSPE[2];
  std::vector<float>  fCalibSPEerr[2];
  std::vector<float>  fCalibTau[2]; 
 
  // ...........................................................
  // Vectors to store all the hit information
  std::vector<int>    fHitPlane;
  std::vector<int>    fHitWire;
  std::vector<float>  fHitT;
  std::vector<float>  fHitX;
  std::vector<float>  fHitW;
  std::vector<float>  fHitIntegral;
  std::vector<float>  fHitSummedADC;
  std::vector<float>  fHitCharge;
  std::vector<float>  fHitRMS;
  std::vector<float>  fHitAmplitude;
  
  // ...........................................................
  // Photon visiblity and timing characterization libraries
  TFile* libraryFile;
  std::vector<std::vector< TH3D* >> libr_XYZ_visibility;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeMin;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeMean;
  std::vector<std::vector< TH3D* >> libr_XYZ_arrivalTimeRMS;
 
  // ...........................................................
  // True muon endpoint
  TH2D* hTrueMuEndpt_ZX;
  TH3D* hTrueMuEndpt;
  
  // ...........................................................
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory truthDir          = tfs->mkdir("true");
  art::TFileDirectory graphDir_hits     = tfs->mkdir("graphs_hits");
  art::TFileDirectory graphDir_wfms     = tfs->mkdir("graphs_wfms");
  art::TFileDirectory diagDir           = tfs->mkdir("diagnostics");
  art::TFileDirectory diagDirOp         = diagDir.mkdir("optical");
  art::TFileDirectory diagDirChHit      = diagDir.mkdir("charge_hits");
  art::TFileDirectory crsMuDir          = tfs->mkdir("crsmu");
  art::TFileDirectory dEdxDir           = tfs->mkdir("dEdx_profiles");

};



//########################################################################################
MichelAna::MichelAna(fhicl::ParameterSet const & p)
 : 
  fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
  ,fCaloAlg(p.get<fhicl::ParameterSet>("CaloAlg"))
{
  LOG_VERBATIM("MichelAna")
  <<"=======================================\n"
  <<"Configuring MichelAna..."; 
  
  // Read fhicl parameters
  this->reconfigure(p);

  if( fUseUniformRecomb ) {
    fRecombFactorTrack = fRecombFactor;
    fRecombFactorShower = fRecombFactor; 
  }

  fWion = 1000./util::kGeVToElectrons; // 23.6 eV, converted to MeV
  //fWph  = 19.5e-6;
 
  /* 
  h_tmp = new TH1D("h_tmp","h_tmp",
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmax());
  h_tmp_n = new TH1D("h_tmp_n","h_tmp_n",
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmax());
  h_tmp_HR = new TH1D("h_tmp_HR","h_tmp",
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmax());
  h_tmp_HRn = new TH1D("h_tmp_HRn","h_tmp_n",
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmax());
  */

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
  LOG_VERBATIM("MichelAna")
  <<"TPC front face center X: "<<frontFaceCenter.X()<<"\n"
  <<"TPC front face center Y: "<<frontFaceCenter.Y()<<"\n"
  <<"TPC front face center Z: "<<frontFaceCenter.Z()<<"\n"
  <<"TPC dimensions: dX = "<<fGeo->DetHalfWidth() * 2.<<"\n"
  <<"TPC dimensions: dY = "<<fGeo->DetHalfHeight() *2.<<"\n"
  <<"TPC dimensions: dZ = "<<fGeo->DetLength()<<"\n"
  <<"TPC min-max X : "<<tpc.MinX()<<" - "<<tpc.MaxX()<<"\n"
  <<"Plane 0 location: "<<tpc.PlaneLocation(0)[0]<<", "<<tpc.PlaneLocation(0)[1]<<", \n"
  <<"Plane 1 location: "<<tpc.PlaneLocation(1)[0]<<", "<<tpc.PlaneLocation(1)[1]<<", \n";
//  <<"Plane 2 location: "<<tpc.PlaneLocation(2)[0]<<", "<<tpc.PlaneLocation(2)[1]<<", \n";
 
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
      LOG_VERBATIM("MichelAna")<<"Loading library file: "<<fLibraryFile;
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

  // Define true muon endpt histogram
  hTrueMuEndpt = truthDir.make<TH3D>("TrueMuEndpt","",24,0.,47.5,20,-20.,20.,45,0.,90.);
  hTrueMuEndpt_ZX = truthDir.make<TH2D>("TrueMuEndpt_ZX","",45,0.,90.,24,0.,47.5);
  
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
    fvHitADC_total[i]   .reserve(kMaxOpHits);
    fvHitADCpf_100ns[i] .reserve(kMaxOpHits);
    fvHitADCpf_total[i] .reserve(kMaxOpHits);
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
    hPMT_phelTimes_muon[i]  = new TH1D(Form("%i_phelTimes_muon",i),Form("PE arrival times for PMT %i",i),28672,0.,28672.);
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
    aveWfmCounter[i]         = 0;
    aveWfmCounter_mu[i]         = 0;
    aveWfmCounter_el[i]         = 0;
    aveWfmCounter_crsmu[i]         = 0;
    avePhelProfileCounter_mu[i] = 0;
  
    fMinHitAmp            .reserve(2);
    fMinHitCharge            .reserve(2);
    
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
    
  // Make vector of SPE values
  cet::search_path sp("FW_SEARCH_PATH");
  for(size_t i = 0; i < 2; i++){
    sprintf(buffer,"LArIATPhotodetectorSER_ch%lu_Michels.txt",i);
    std::string fullname;
    sp.find_file(buffer, fullname);
    if (fullname.empty()) {
      throw cet::exception("PMT calibration file not found!");
    } else {
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        size_t lineNum = 0; int run1 = 0.; int run2 = 0.;
        while (std::getline(inFile,line)) {
          if(lineNum == 1 ) continue;
          float spe, speErr, width, tau;
          std::istringstream ss(line);
          ss >> run1 >> run2 >> spe >> speErr >> width >> tau;
          fCalibRunStart[i].push_back(run1);
          fCalibRunEnd[i].push_back(run2);
          fCalibSPE[i].push_back(spe);
          fCalibSPEerr[i].push_back(speErr);
          fCalibTau[i].push_back(tau);
        }
    }
  }

  // Angle bins (copying ArgoNeuT proton/deuteron study)
  fAngleBins[3] = 40.;
  fAngleBins[2] = 50.;
  fAngleBins[1] = 60.;
  fAngleBins[0] = 80.;
  fAngleBins_X1[3] = 20.;
  fAngleBins_X1[2] = 47.;
  fAngleBins_X1[1] = 55.;
  fAngleBins_X1[0] = 70.;

  LOG_VERBATIM("MichelAna")
  <<"Configuration complete.\n"
  <<"=======================================";

}



//########################################################################################
void MichelAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRandSeed                 = p.get< int >                  ("RandSeed", 1989);
  fUsePmtForCalo            = p.get< std::vector<bool> >    ("UsePmtForCalo",{true,true});
  fMinHitAmp                = p.get< std::vector<float >>   ("MinHitAmp",{0,0});
  fMinHitCharge                = p.get< std::vector<float >>   ("MinHitCharge",{0,0});
  fMuContamCorr_EffTau      = p.get< std::vector<float >>   ("MuContamCorr_EffTau",{0,0});
  fQE_ScaleFactor           = p.get< std::vector<float >>   ("QE_ScaleFactor",{1.,1.});
  fApplyTrigEffCut          = p.get< bool >                 ("ApplyTrigEffCut",false);
  fTrigEff_P                = p.get< std::vector< float >>  ("TrigEff_P",{80.,25.});
  fTrigEff_K                = p.get< std::vector< float >>  ("TrigEff_K",{8,8});
  fCorrectQuenching         = p.get< bool >                 ("CorrectQuenching",false);
  fWph                      = p.get< float >                ("Wph",19.5e-6);
  fMakeAveWfms              = p.get< bool >                 ("MakeAveWfms",true);
  fMuContamCorr             = p.get< bool >                 ("MuContamCorr",true);
  fMuContamCorr_dTmin       = p.get< float >                ("MuContamCorr_dTmin",2000);
  fMskBaselineSubtr         = p.get< bool >                 ("MskBaselineSubtr", false);
  fTruncateWfm              = p.get< short >                ("TruncateWfm",-1);
  fWfmSmoothingRange        = p.get< int  >                 ("WfmSmoothingRange",1);
  fCorrectOvershootMode     = p.get< std::string >          ("CorrectOvershootMode", "");
  fUseCrossingMuons         = p.get< bool >                 ("UseCrossingMuons",true);
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
  fCollectionEff            = p.get< std::vector< float >>  ("CollectionEff",{1.,1.});
  fQE_ScaleFactor_Vis       = p.get< std::vector< float >>  ("QE_ScaleFactor_Vis",{1.,1.});
  fQE_ScaleFactor_VUV       = p.get< std::vector< float >>  ("QE_ScaleFactor_VUV",{1.,1.});
  fRequireHitMatching       = p.get< bool >                 ("RequireHitMatching",false);
  fDeleteNarrowHits         = p.get< bool >                 ("DeleteNarrowHits",false);
  fElectronLifetimeErr      = p.get< float >                ("ElectronLifetimeErr",0.);
  fSmearSigT                = p.get< std::vector< float >>  ("SmearSigT",{0.02, 0.04 });
  fSmearSigPE               = p.get< std::vector< float >>  ("SmearSigPE",{0.00, 0.00 });
  fSmearFactor              = p.get< std::vector< float >>  ("SmearFactor",{0.,0.});
  fSmearFactorPrompt        = p.get< std::vector< float >>  ("SmearFactorPrompt",{0.,0.});
  fGateDelay                = p.get< float >                ("GateDelay",300);
  fMaxDecayTime             = p.get< float >                ("MaxDecayTime",7200);
  fTrackModule              = p.get< std::string >          ("TrackModule","pmtrack");
  fTrackCalModule           = p.get< std::string >          ("TrackCalModule","calo");
  fSimProducerModule        = p.get< std::string >          ("SimProducerModule","largeant");
  fHitsModule               = p.get< std::string >          ("HitsModule","gaushit");
  fHitsInstance               = p.get< std::string >          ("HitsInstance","");
  fLibraryFile              = p.get< std::string >          ("LibraryFile","");
  fSingletToTripletRatio           = p.get< float >                ("SingletToTripletRatio",0.30);
  fLateLightTau             = p.get< float >                ("LateLightTau",1200);
  fTauT                     = p.get< float >                ("TauT",1300);
  fTauS                     = p.get< float >                ("TauS",6.5);
  //fCalAreaConstants         = p.get< std::vector< float >>  ("CalAreaConstants", {0.022, 0.050});
  fEdgeMarginX              = p.get< float >                ("EdgeMarginX",1.);
  fEdgeMarginY              = p.get< float >                ("EdgeMarginY",1.);
  fEdgeMarginZ              = p.get< float >                ("EdgeMarginZ",1.);
  fMuTrackLengthMin         = p.get< float >                ("MuTrackLengthMin",5.);
  fMuTrackZenithAngleMax    = p.get< float >                ("MuTrackZenithAngleMax",180.);
  fMuResRangeOffset         = p.get< float >                ("MuResRangeOffset",0.10);
  fUseMaxDropAsBnd          = p.get< bool >                 ("UseMaxDropAsBnd",false);
  fMaxHitSeparation         = p.get< float >                ("MaxHitSeparation",2.5);
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
  fMinDecayAngle2D          = p.get< float >                ("DecayAngle2DMin",10.);
  fMaxDecayAngle2D          = p.get< float >                ("DecayAngle2DMax",170.);
  fShowerAcceptanceAngle    = p.get< float >                ("ShowerAcceptanceAngle",60.);
  fExcRatio                 = p.get< float >                ("ExcRatio", 0.21 ); // Aprile 
  fUseUniformRecomb         = p.get< bool >             ("UseUniformRecomb", false);
  fRecombFactor             = p.get< float >                ("RecombFactor", 0.645);
  fRecombFactorTrack        = p.get< float >                ("RecombFactorTrack", 0.688);
  fRecombFactorShower       = p.get< float >                ("RecombFactorShower", 0.576);
//  fChargeSmearFactor        = p.get< float >                ("ChargeSmearFactor",0.);
  fMaxHitTimeDiffFac        = p.get< float >             ("MaxHitTimeDiffFac", 1.0 );
  fMinNumPts3D              = p.get< int >                  ("MinNumPts3D",2);
  fMinFracHits3D            = p.get< float >                ("MinFracHits3D",0.3);
  fLinThreshMuFit           = p.get< float >                ("LinThreshMuFit",0.8);
  fMinTruncMeanSlope        = p.get< float >                ("MinTruncMeanSlope",1000.);
  fMinMuClusterHitsEndFit   = p.get< int >                  ("MinMuClusterHitsEndFit",5);
  fMinMuLinearity           = p.get< float >                ("MinMuLinearity",0.9);
  fMinMuClusterSize         = p.get< int >                  ("MinMuClusterSize",8);
  fMinFracMuHitsLinear      = p.get< float >                ("MinFracMuHitsLinear",0.5);  
  fBndOffset                = p.get< int >                  ("BndOffset",0);
  fMaxSavedHitGraphs        = p.get< int >                  ("MaxSavedHitGraphs",0);
  fMaxSavedWfmGraphs        = p.get< int >                  ("MaxSavedWfmGraphs",0);
}




//########################################################################################
void MichelAna::beginJob()
{
  
  // ====================================================================
  // Define TTree and set all branches
  fTree= tfs->make<TTree>("anatree","anatree");
  // ...event ID info, calibration values used
  fTree->Branch("RunNumber",                &fRunNumber,              "RunNumber/I");
  fTree->Branch("SubRunNumber",             &fSubRunNumber,           "SubRunNumber/I");
//  fTree->Branch("EventNumber",              &fEventNumber,            "EventNumber/I");
  fTree->Branch("EventTime",                &fEventTime,              "EventTime/I");
  fTree->Branch("MichelOpticalID",          &fMichelOpticalID,        "MichelOpticalID/O");
  fTree->Branch("ElectronLifetime",         &fElectronLifetime,       "ElectronLifetime/F");
  fTree->Branch("ElectronLifetimeFromDB",         &fElectronLifetimeFromDB,       "ElectronLifetimeFromDB/F");
  fTree->Branch("MuContamCorr_EffTau",      &fMuContamCorr_EffTau);
  fTree->Branch("QE_ScaleFactor",           &fQE_ScaleFactor);
  //...track information
  fTree->Branch("NumTrack",                 &fNumTrack,               "NumTrack/I");
  fTree->Branch("NumTrackStopping",         &fNumTrackStopping,       "NumTrackStopping/I");
  fTree->Branch("NumTrackContained",         &fNumTrackContained,       "NumTrackContained/I");
  fTree->Branch("NumTrackCrossing",         &fNumTrackCrossing,       "NumTrackCrossing/I");
  fTree->Branch("MuTrackLength",            &fMuTrackLength,          "MuTrackLength/F");
  fTree->Branch("MuTrackEnergy",            &fMuTrackEnergy,          "MuTrackEnergy/F");
  fTree->Branch("MuTrackZenithAngle",       &fMuTrackZenithAngle,     "MuTrackZenithAngle/F");
  fTree->Branch("MuTrackInclineAngle",       &fMuTrackInclineAngle,     "MuTrackInclineAngle/F");
  fTree->Branch("MuTrackEnd_X",             &fMuTrackEnd_X,           "MuTrackEnd_X/F");
  fTree->Branch("MuTrackEnd_Y",             &fMuTrackEnd_Y,           "MuTrackEnd_Y/F");
  fTree->Branch("MuTrackEnd_Z",             &fMuTrackEnd_Z,           "MuTrackEnd_Z/F");
//  fTree->Branch("MuTrackPitch",          &fMuTrackPitch,        "MuTrackPitch[2]/F");
  // ...optical information
//  fTree->Branch("NumOpHits0",               &fNumOpHits0,             "NumOpHits0[2]/I");
  fTree->Branch("NumOpHits",                &fNumOpHits,              "NumOpHits[2]/I");
  fTree->Branch("SPE",                      &fSPE,                    "SPE[2]/F");
  //fTree->Branch("SPE_err",                  &fSPE_err,                "SPE_err[2]/F");
  fTree->Branch("dT",                       &fdT,                     "dT[2]/F");
  fTree->Branch("Amplitude",                &fAmplitude,              "Amplitude[2]/F");
  fTree->Branch("Width",                    &fWidth,                  "Width[2]/F");
  fTree->Branch("MuWidth",                  &fMuWidth,                "MuWidth[2]/F");
  fTree->Branch("PE_comb",                &fPE_comb,              "PE_comb/F");
  fTree->Branch("PE_comb_qc",                &fPE_comb_qc,              "PE_comb_qc/F");
  fTree->Branch("MuContam_prompt",                &fMuContam_prompt,              "MuContam_prompt[2]/F");
  fTree->Branch("MuContam_total",                &fMuContam_total,              "MuContam_total[2]/F");
  fTree->Branch("PE_prompt",                &fPE_prompt,              "PE_prompt[2]/F");
  fTree->Branch("PE_total",                 &fPE_total,               "PE_total[2]/F");
  fTree->Branch("PE_total_qc",                 &fPE_total_qc,               "PE_total_qc[2]/F");
  fTree->Branch("PE_promptRaw",             &fPE_promptRaw,           "PE_promptRaw[2]/F");
  fTree->Branch("PE_totalRaw",              &fPE_totalRaw,            "PE_totalRaw[2]/F");
  fTree->Branch("PromptFrac",               &fPromptFrac,             "PromptFrac[2]/F");
  fTree->Branch("MuAmplitude",              &fMuAmplitude,            "MuAmplitude[2]/F");
  fTree->Branch("MuPE_prompt",              &fMuPE_prompt,            "MuPE_prompt[2]/F");
  fTree->Branch("MuPE_total",              &fMuPE_total,            "MuPE_total[2]/F");
  fTree->Branch("MuPulseSaturated",         &fMuPulseSaturated,       "MuPulseSaturated[2]/O");
  fTree->Branch("ElPulseSaturated",         &fElPulseSaturated,       "ElPulseSaturated[2]/O");
  fTree->Branch("DecayTime",                &fDecayTime,              "DecayTime/F");
  //...clustering/showering information
  fTree->Branch("NumPlaneHits",             &fNumPlaneHits,           "NumPlaneHits[2]/I");
  fTree->Branch("FracDropoffAtBnd",              &fFracDropoffAtBnd,            "FracDropoffAtBnd/F"); 
  fTree->Branch("CovAtBnd",              &fCovAtBnd,            "CovAtBnd/F"); 
  fTree->Branch("TruncMeanSlope",              &fTruncMeanSlope,            "TruncMeanSlope/F"); 
  fTree->Branch("ClusterSize",              &fClusterSize,            "ClusterSize/I"); 
  fTree->Branch("MuEnd2D_X",            &fMuEnd2D_X,          "MuEnd2D_X/F");
  fTree->Branch("MuEnd2D_W",            &fMuEnd2D_W,          "MuEnd2D_W/F");
  fTree->Branch("ElDir2D_X",            &fElDir2D_X,          "ElDir2D_X/F");
  fTree->Branch("ElDir2D_W",            &fElDir2D_W,          "ElDir2D_W/F");
  fTree->Branch("ElClusterSize",            &fElClusterSize,          "ElClusterSize/I");
  fTree->Branch("MuClusterSize",            &fMuClusterSize,          "MuClusterSize/I");
  fTree->Branch("ExtraHits",                &fExtraHits,              "ExtraHits/I");
  fTree->Branch("MuAveLinearity",           &fMuAveLinearity,         "MuAveLinearity/F");
  fTree->Branch("FracMuHitsLinear",         &fFracMuHitsLinear,       "FracMuHitsLinear/F");
  fTree->Branch("MuClusterHitsEndFit",      &fMuClusterHitsEndFit,    "MuClusterHitsEndFit/I");
  fTree->Branch("DecayAngle2D",             &fDecayAngle2D,           "DecayAngle2D/F");
  fTree->Branch("ElShowerSize",             &fElShowerSize,           "ElShowerSize/I");
  fTree->Branch("ElShowerFrac",             &fElShowerFrac,           "ElShowerFrac/F");
  fTree->Branch("ProjDist3DToWires",        &fProjDist3DToWires,           "ProjDist3DToWires/F");
  fTree->Branch("ProjPtOnWires_Y",          &fProjPtOnWires_Y,        "ProjPtOnWires_Y/F");
  fTree->Branch("ProjPtOnWires_Z",          &fProjPtOnWires_Z,        "ProjPtOnWires_Z/F");
  fTree->Branch("TotalCharge",              &fTotalCharge,            "TotalCharge[2]/F");
  fTree->Branch("TotalIntegral",            &fTotalIntegral,       "TotalIntegral[2]/F");
  fTree->Branch("ElShowerCharge",           &fElShowerCharge,         "ElShowerCharge/F");
  fTree->Branch("ElShowerEnergy",           &fElShowerEnergy,         "ElShowerEnergy/F");
  fTree->Branch("ElTrackCharge",            &fElTrackCharge,          "ElTrackCharge/F");
  fTree->Branch("ElTrackEnergy",            &fElTrackEnergy,          "ElTrackEnergy/F");
  fTree->Branch("MuCharge",           &fMuCharge,         "MuCharge/F");
  fTree->Branch("FracDropoffAtBnd_Pl0",              &fFracDropoffAtBnd_Pl0,            "FracDropoffAtBnd_Pl0/F"); 
  fTree->Branch("CovAtBnd_Pl0",          &fCovAtBnd_Pl0,        "CovAtBnd_Pl0/F"); 
  fTree->Branch("TruncMeanSlope_Pl0",              &fTruncMeanSlope_Pl0,            "TruncMeanSlope_Pl0/F"); 
  fTree->Branch("ClusterSize_Pl0",          &fClusterSize_Pl0,        "ClusterSize_Pl0/I"); 
  fTree->Branch("MuEnd2D_X_Pl0",            &fMuEnd2D_X_Pl0,          "MuEnd2D_X_Pl0/F");
  fTree->Branch("MuEnd2D_W_Pl0",            &fMuEnd2D_W_Pl0,          "MuEnd2D_W_Pl0/F");
  fTree->Branch("ElDir2D_X_Pl0",            &fElDir2D_X_Pl0,          "ElDir2D_X_Pl0/F");
  fTree->Branch("ElDir2D_W_Pl0",            &fElDir2D_W_Pl0,          "ElDir2D_W_Pl0/F");
  fTree->Branch("ElClusterSize_Pl0",        &fElClusterSize_Pl0,      "ElClusterSize_Pl0/I");
  fTree->Branch("MuClusterSize_Pl0",        &fMuClusterSize_Pl0,      "MuClusterSize_Pl0/I");
  fTree->Branch("ExtraHits_Pl0",            &fExtraHits_Pl0,          "ExtraHits_Pl0/I");
  fTree->Branch("MuAveLinearity_Pl0",       &fMuAveLinearity_Pl0,     "MuAveLinearity_Pl0/F");
  fTree->Branch("FracMuHitsLinear_Pl0",     &fFracMuHitsLinear_Pl0,   "FracMuHitsLinear_Pl0/F");
  fTree->Branch("MuClusterHitsEndFit_Pl0",  &fMuClusterHitsEndFit_Pl0,"MuClusterHitsEndFit_Pl0/I");
  fTree->Branch("DecayAngle2D_Pl0",         &fDecayAngle2D_Pl0,       "DecayAngle2D_Pl0/F");
  fTree->Branch("ElShowerSize_Pl0",         &fElShowerSize_Pl0,       "ElShowerSize_Pl0/I");
  fTree->Branch("ElShowerFrac_Pl0",         &fElShowerFrac_Pl0,       "ElShowerFrac_Pl0/F");
  fTree->Branch("ElShowerCharge_Pl0",       &fElShowerCharge_Pl0,     "ElShowerCharge_Pl0/F");
  fTree->Branch("ElShowerEnergy_Pl0",       &fElShowerEnergy_Pl0,     "ElShowerEnergy_Pl0/F");
  fTree->Branch("ElTrackCharge_Pl0",        &fElTrackCharge_Pl0,      "ElTrackCharge_Pl0/F");
  fTree->Branch("ElTrackEnergy_Pl0",        &fElTrackEnergy_Pl0,      "ElTrackEnergy_Pl0/F");
  fTree->Branch("MuEndHitTimeDiff",         &fMuEndHitTimeDiff,       "MuEndHitTimeDiff/F");
  fTree->Branch("ElShowerVis",              &fElShowerVis,            "ElShowerVis/F");
  fTree->Branch("ElShowerVisCh",            &fElShowerVisCh,          "ElShowerVisCh[2]/F");
  fTree->Branch("ElShowerPhotons",          &fElShowerPhotons,        "ElShowerPhotons/F");
  fTree->Branch("ElShowerPhotonsPrompt",          &fElShowerPhotonsPrompt,        "ElShowerPhotonsPrompt/F");
//  fTree->Branch("ElTrackdEdx",              &fElTrackdEdx,            "ElTrackdEdx/F");
  fTree->Branch("MuEnd3D_X",                &fMuEnd3D_X,              "MuEnd3D_X/F");
  fTree->Branch("MuEnd3D_Y",                &fMuEnd3D_Y,              "MuEnd3D_Y/F");
  fTree->Branch("MuEnd3D_Z",                &fMuEnd3D_Z,              "MuEnd3D_Z/F");
  fTree->Branch("ElDir3D_X",                &fElDir3D_X,              "ElDir3D_X/F");
  fTree->Branch("ElDir3D_Y",                &fElDir3D_Y,              "ElDir3D_Y/F");
  fTree->Branch("ElDir3D_Z",                &fElDir3D_Z,              "ElDir3D_Z/F");
  fTree->Branch("ElShowerCentroid_X",       &fElShowerCentroid_X,     "ElShowerCentroid_X/F");
  fTree->Branch("ElShowerCentroid_Y",       &fElShowerCentroid_Y,     "ElShowerCentroid_Y/F");
  fTree->Branch("ElShowerCentroid_Z",       &fElShowerCentroid_Z,     "ElShowerCentroid_Z/F");
  fTree->Branch("FracHits3D",               &fFracHits3D,             "FracHits3D/F");
  fTree->Branch("NumPts3D",                 &fNumPts3D,               "NumPts3D/I");
  //...truth information
  fTree->Branch("True_IsMuStopping",      &fTrue_IsMuStopping,        "True_IsMuStopping/O");
  fTree->Branch("True_MuSign",            &fTrue_MuSign,          "True_MuSign/I");
  fTree->Branch("True_MuCharge",          &fTrue_MuCharge,            "True_MuCharge/F");
  fTree->Branch("True_MuPhotons",          &fTrue_MuPhotons,            "True_MuPhotons/F");
  fTree->Branch("True_MuEnergyDep",          &fTrue_MuEnergyDep,            "True_MuEnergyDep/F");
  fTree->Branch("True_MuTrackEnd_X",        &fTrue_MuTrackEnd_X,      "True_MuTrackEnd_X/F");
  fTree->Branch("True_MuTrackEnd_Y",        &fTrue_MuTrackEnd_Y,      "True_MuTrackEnd_Y/F");
  fTree->Branch("True_MuTrackEnd_Z",        &fTrue_MuTrackEnd_Z,      "True_MuTrackEnd_Z/F");
  fTree->Branch("True_ProjDist3DToWires",   &fTrue_ProjDist3DToWires, "True_ProjDist3DToWires/F");
  fTree->Branch("True_MuPE",                &fTrue_MuPE,              "True_MuPE/F");
  fTree->Branch("True_MuPE_prompt",         &fTrue_MuPE_prompt,       "True_MuPE_prompt[2]/F");
  fTree->Branch("True_MuPE_total",         &fTrue_MuPE_total,         "True_MuPE_total[2]/F");
  fTree->Branch("True_MuContam_prompt",     &fTrue_MuContam_prompt,   "MuContam_prompt[2]/F");
  fTree->Branch("True_MuContam_total",      &fTrue_MuContam_total,    "MuContam_total[2]/F");
  fTree->Branch("True_ElAngle",             &fTrue_ElAngle,           "True_ElAngle/F");
  fTree->Branch("True_ElEnergy",            &fTrue_ElEnergy,          "True_ElEnergy/F");
  fTree->Branch("True_ElShowerEnergyDep",   &fTrue_ElShowerEnergyDep, "True_ElShowerEnergyDep/F");
  fTree->Branch("True_ElShowerChargeDep",   &fTrue_ElShowerChargeDep, "True_ElShowerChargeDep/F");
  fTree->Branch("True_ElShowerTotalTrajLength",   &fTrue_ElShowerTotalTrajLength, "True_ElShowerTotalTrajLength/F");
  fTree->Branch("True_ElShowerIons",   &fTrue_ElShowerIons, "True_ElShowerIons/F");
  fTree->Branch("True_ElTrackChargeDep",    &fTrue_ElTrackChargeDep,  "True_ElTrackChargeDep/F");
  fTree->Branch("True_ElTrackEnergyDep",    &fTrue_ElTrackEnergyDep,  "True_ElTrackEnergyDep/F");
  fTree->Branch("True_ElShowerVis",       &fTrue_ElShowerVis,     "True_ElShowerVis/F");
  fTree->Branch("True_ElShowerVisCh",       &fTrue_ElShowerVisCh,     "True_ElShowerVisCh[2]/F");
  fTree->Branch("True_TotalEnergyDep",      &fTrue_TotalEnergyDep,    "True_TotalEnergyDep/F");
  fTree->Branch("True_IsElTrkContained",       &fTrue_IsElTrkContained,     "True_IsElTrkContained/O");
  fTree->Branch("True_IsElShwrContained",       &fTrue_IsElShwrContained,     "True_IsElShwrContained/O");
  fTree->Branch("True_ContainmentFrac",     &fTrue_ContainmentFrac,   "True_ContainmentFrac/F");
  fTree->Branch("True_ElShowerCentroid_X",  &fTrue_ElShowerCentroid_X,"True_ElShowerCentroid_X/F");
  fTree->Branch("True_ElShowerCentroid_Y",  &fTrue_ElShowerCentroid_Y,"True_ElShowerCentroid_Y/F");
  fTree->Branch("True_ElShowerCentroid_Z",  &fTrue_ElShowerCentroid_Z,"True_ElShowerCentroid_Z/F");
  fTree->Branch("True_ElShowerPhotons",     &fTrue_ElShowerPhotons,   "True_ElShowerPhotons/F");
  fTree->Branch("True_ElTrackPhotons",      &fTrue_ElTrackPhotons,    "True_ElTrackPhotons/F");
  fTree->Branch("True_ElShowerPhotonsPrompt",&fTrue_ElShowerPhotonsPrompt,"True_ElShowerPhotonsPrompt/F");
  fTree->Branch("True_ElShowerPhotonsLate", &fTrue_ElShowerPhotonsLate,"True_ElShowerPhotonsLate/F");
  fTree->Branch("True_dT",                  &fTrue_dT,                "True_dT/F");
  fTree->Branch("True_PE",           &fTrue_PE,         "True_PE/F");
  fTree->Branch("True_PE_prompt",           &fTrue_PE_prompt,         "True_PE_prompt[2]/F");
  fTree->Branch("True_PE_total",            &fTrue_PE_total,          "True_PE_total[2]/F");

  // TTree dedicated to crossing muon calibration
  fTreeCrsMu = crsMuDir.make<TTree>("anatree_crsmu","anatree_crsmu");
  fTreeCrsMu->Branch("RunNumber",                &fRunNumber,              "RunNumber/I");
  fTreeCrsMu->Branch("SubRunNumber",             &fSubRunNumber,           "SubRunNumber/I");
  fTreeCrsMu->Branch("EventNumber",              &fEventNumber,            "EventNumber/I");
  fTreeCrsMu->Branch("EventTime",                &fEventTime,              "EventTime/I");
  fTreeCrsMu->Branch("MichelOpticalID",          &fMichelOpticalID,        "MichelOpticalID/O");
  fTreeCrsMu->Branch("ElectronLifetime",         &fElectronLifetime,       "ElectronLifetime/F");
  fTreeCrsMu->Branch("ElectronLifetimeFromDB",         &fElectronLifetimeFromDB,       "ElectronLifetimeFromDB/F");
  fTreeCrsMu->Branch("MuContamCorr_EffTau",      &fMuContamCorr_EffTau);
  fTreeCrsMu->Branch("QE_ScaleFactor",           &fQE_ScaleFactor);
  fTreeCrsMu->Branch("CrsMuVertex_X",         &fCrsMuVertex_X,       "CrsMuVertex_X/F");
  fTreeCrsMu->Branch("CrsMuVertex_Y",         &fCrsMuVertex_Y,       "CrsMuVertex_Y/F");
  fTreeCrsMu->Branch("CrsMuVertex_Z",         &fCrsMuVertex_Z,       "CrsMuVertex_Z/F");
  fTreeCrsMu->Branch("CrsMuEnd_X",         &fCrsMuEnd_X,       "CrsMuEnd_X/F");
  fTreeCrsMu->Branch("CrsMuEnd_Y",         &fCrsMuEnd_Y,       "CrsMuEnd_Y/F");
  fTreeCrsMu->Branch("CrsMuEnd_Z",         &fCrsMuEnd_Z,       "CrsMuEnd_Z/F");
  fTreeCrsMu->Branch("CrsMuLength",         &fCrsMuLength,       "CrsMuLength/F");
  fTreeCrsMu->Branch("CrsMuCharge",         &fCrsMuCharge,       "CrsMuCharge/F");
  fTreeCrsMu->Branch("CrsMuCharge_Pl0",         &fCrsMuCharge_Pl0,       "CrsMuCharge_Pl0/F");
  fTreeCrsMu->Branch("CrsMuIntegral",         &fCrsMuIntegral,       "CrsMuIntegral/F");
  fTreeCrsMu->Branch("CrsMuIntegral_Pl0",         &fCrsMuIntegral_Pl0,       "CrsMuIntegral_Pl0/F");
  fTreeCrsMu->Branch("CrsMuPE_total",              &fCrsMuPE_total,            "CrsMuPE_total[2]/F");
  fTreeCrsMu->Branch("CrsMuPE_prompt",              &fCrsMuPE_prompt,            "CrsMuPE_prompt[2]/F");
  fTreeCrsMu->Branch("CrsMuPhotons",          &fCrsMuPhotons,        "CrsMuPhotons/F");
  fTreeCrsMu->Branch("CrsMuPhotonsPrompt",  &fCrsMuPhotonsPrompt,"CrsMuPhotonsPrompt/F");
  fTreeCrsMu->Branch("CrsMuPitch",          &fCrsMuPitch,        "CrsMuPitch[2]/F");
  fTreeCrsMu->Branch("CrsMuHitFrac",          &fCrsMuHitFrac,        "CrsMuHitFrac/F");
  fTreeCrsMu->Branch("CrsMuEnergy",          &fCrsMuEnergy,        "CrsMuEnergy/F");
  fTreeCrsMu->Branch("True_CrsMuCharge",         &fTrue_CrsMuCharge,       "True_CrsMuCharge/F");
  fTreeCrsMu->Branch("True_CrsMuPhotons",         &fTrue_CrsMuPhotons,       "True_CrsMuPhotons/F");
 
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
  hElectronLifetimeFromDB       = diagDir.make<TH1D>("ElectronLifetimeFromDB","Electron lifetime from database;Electron lifetime [#mus]",2500,-100.,2400.);
  hElectronLifetime       = diagDir.make<TH1D>("ElectronLifetime","Electron lifetime (after MC smearing);Electron lifetime [#mus]",2500,-100.,2400.);

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
  int Nbins_michelEnergy  = 60;
  
  hPE       = tfs->make<TH1D>("PE","Total Michel light for event;Michel light collected [PEs]",200,0,2000);
  
  // ----------------------------------------------------------------------
  // Loop over the PMTs and make PMT-specific histograms
  for(size_t i=0; i<fSelectChannels.size(); i++){
    size_t  ch     = fSelectChannels.at(i);
 
    // ------------------------------------------------------------
    // Average waveforms and binned photon arrival time profiles
    sprintf(histName,"%lu_AveWfm",ch);
    sprintf(histTitle,"Averaged scintillation waveform, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm[ch]  = diagDirOp.make<TH1D>(histName,histTitle,7700,-700.,7000.);
    sprintf(histName,"%lu_AveWfm_mu",ch);
    sprintf(histTitle,"Averaged waveform of #mu pulse, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm_mu[ch]  = diagDirOp.make<TH1D>(histName,histTitle,2200,-200.,2000.);
    sprintf(histName,"%lu_AveWfm_el",ch);
    sprintf(histTitle,"Averaged waveform of Michel pulse, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm_el[ch]  = diagDirOp.make<TH1D>(histName,histTitle,7700,-700.,7000.);
    sprintf(histName,"%lu_AveWfm_crsmu",ch);
    sprintf(histTitle,"Averaged waveform of crossing muon pulse, opdet %lu;Time Tick [ns];Averaged Signal [ADC]",ch);
    hAveWfm_crsmu[ch]  = diagDirOp.make<TH1D>(histName,histTitle,7700,-700.,7000.);
    float bins[] = { 0, 400, 500, 700, 900, 1200, 1500, 1800};
    int binnum = sizeof(bins)/sizeof(float) - 1;
    sprintf(histName,"%lu_AvePhelProfile_mu",ch);
    sprintf(histTitle,"Averaged photon arrival time profile of #mu pulse, opdet %lu;Time after #mu pulse [ns];PEs per ns",ch);
    hAvePhelProfile_mu[ch] = diagDirOp.make<TH1D>(histName,histTitle, binnum, bins);

    // ------------------------------------------------------------
    // PMT diagnostic plots
    sprintf(histName,"%lu_WfmRMS",ch);
    sprintf(histTitle,"Waveform RMS, opdet %lu;ADC",ch);
    hWfmRMS[ch]   = diagDirOp.make<TH1D>(histName,histTitle,100,0.,5.);
    sprintf(histName,"%lu_BaselinePE_100ns",ch);
    sprintf(histTitle,"Integrated baseline region (100ns), opdet %lu;PE",ch);
    hBaselinePE_100ns[ch] = diagDirOp.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_BaselinePE_500ns",ch);
    sprintf(histTitle,"Integrated baseline region (500ns), opdet %lu;PE",ch);
    hBaselinePE_500ns[ch] = diagDirOp.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_BaselinePE_1000ns",ch);
    sprintf(histTitle,"Integrated baseline region (1000ns), opdet %lu;PE",ch);
    hBaselinePE_1000ns[ch] = diagDirOp.make<TH1D>(histName,histTitle,200,-8.,8.);
    sprintf(histName,"%lu_NumOpHits",ch);
    sprintf(histTitle,"Opdet %lu;Number of optical hits",ch);
    hNumOpHits[ch]   = diagDirOp.make<TH1D>(histName,histTitle,20,0,20);
    sprintf(histName,"%lu_MuTau",ch);
    sprintf(histTitle,"Opdet %lu, Lifetime from #mu late-light fit (for contamination correction);ns",ch);
    hMuTau[ch]     = diagDirOp.make<TH1D>(histName,histTitle, 250, 400., 5400.);
    hContamFrac[ch] = diagDirOp.make<TH1D>(Form("%lu_ContamFrac",ch),Form("Opdet %lu, fraction of estimated contamination from #mu late-light;Fraction",ch),120.,0.,1.2);
    hContamFrac_dTcut[ch] = diagDirOp.make<TH1D>(Form("%lu_ContamFrac_dTcut",ch),Form("Opdet %lu, fraction estimated contamination from #mu late-light (dT > %4.0f ns);Fraction",ch,fdTcut),120.,0.,1.2);
   
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
    hOpHitTime[ch]          = diagDirOp.make<TH1D>(histName,histTitle,200,-1000,1000);
    sprintf(histName,"%lu_dT",ch);
    sprintf(histTitle,"Opdet %lu;#DeltaT [ns]",ch);
    hdT[ch]          = diagDirOp.make<TH1D>(histName,histTitle,750,0.,7500.);
    hdT_vs_Amp[ch]   = diagDirOp.make<TH2D>(Form("%lu_dT_vs_Amp",ch),Form("Opdet %lu;#DeltaT [ns];Amplitude [ADC]",ch),150,0.,7500.,102,-20,1000);
    hdT_vs_Amp[ch]  ->SetOption("colz");
    hdT_vs_PE_totalRaw[ch]  = diagDirOp.make<TH2D>(Form("%lu_dT_vs_PE_totalRaw",ch),Form("Opdet %lu (raw integral);#Delta T [ns];Total light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S[ch], x1_S[ch], x2_S[ch]);
    hdT_vs_PE_total[ch]     = diagDirOp.make<TH2D>(Form("%lu_dT_vs_PE_total",ch),Form("Opdet %lu;#Delta T [ns];Total light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S[ch], x1_S[ch], x2_S[ch]);
    hdT_vs_PE_promptRaw[ch]  = diagDirOp.make<TH2D>(Form("%lu_dT_vs_PE_promptRaw",ch),Form("Opdet %lu (raw integral);#Delta T [ns];Prompt light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S100[ch], x1_S100[ch], x2_S100[ch]);
    hdT_vs_PE_prompt[ch]     = diagDirOp.make<TH2D>(Form("%lu_dT_vs_PE_prompt",ch),Form("Opdet %lu;#Delta T [ns];Prompt light of Michel candidate pulse [pe]",ch),
      75,0.,7500.,nbins_S100[ch], x1_S100[ch], x2_S100[ch]);
    hdT_vs_PE_totalRaw[ch]  ->SetOption("colz");
    hdT_vs_PE_total[ch]     ->SetOption("colz");
    hdT_vs_PE_promptRaw[ch]  ->SetOption("colz");
    hdT_vs_PE_prompt[ch]     ->SetOption("colz");
    
    // -------------------------------------------------------------
    // General PMT histograms
    
    sprintf(histName,"%lu_PromptFracRaw",ch);
    sprintf(histTitle,"Opdet %lu (raw integral);Prompt fraction of Michel candidate pulse",ch);
    hPromptFracRaw[ch]     = diagDirOp.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFrac",ch);
    sprintf(histTitle,"Opdet %lu;Prompt fraction of Michel candidate pulse",ch);
    hPromptFrac[ch]     = diagDirOp.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFracRaw_dTcut",ch);
    sprintf(histTitle,"Opdet %lu (raw integral);Prompt fraction of Michel candidate pulse",ch);
    hPromptFracRaw_dTcut[ch]     = diagDirOp.make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    sprintf(histName,"%lu_PromptFrac_dTcut",ch);
    sprintf(histTitle,"Opdet %lu;Prompt fraction of Michel candidate pulse",ch);
    hPromptFrac_dTcut[ch]     = tfs->make<TH1D>(histName,histTitle, 140, -0.2, 1.2);
    
    
    sprintf(histName,"%lu_OpHitAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of ALL found hits [ADC]",ch);
    hOpHitAmplitude[ch]          = diagDirOp.make<TH1D>(histName,histTitle,120,0.,1200.);
    sprintf(histName,"%lu_OpHitAmplitude_unmatched",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of hits not matched to another PMT [ADC]",ch);
    hOpHitAmplitude_unmatched[ch]          = diagDirOp.make<TH1D>(histName,histTitle,120,0.,1200.);
    sprintf(histName,"%lu_Amplitude",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of Michel candidate pulse [ADC]",ch);
    hAmplitude[ch]          = diagDirOp.make<TH1D>(histName,histTitle,120,0.,1200.);
    
    sprintf(histName,"%lu_HitPromptPE",ch);
    sprintf(histTitle,"Opdet %lu;PromptPE of ALL found hits [pe]",ch);
    hHitPromptPE[ch]          = diagDirOp.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_HitPromptPE_unmatched",ch);
    sprintf(histTitle,"Opdet %lu;PromptPE of hits not matched to another PMT [pe]",ch);
    hHitPromptPE_unmatched[ch]          = diagDirOp.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    
    sprintf(histName,"%lu_Width",ch);
    sprintf(histTitle,"Opdet %lu;Width of Michel candidate pulse [ns]",ch);
    hWidth[ch]          = diagDirOp.make<TH1D>(histName,histTitle,120,0,30);
    sprintf(histName,"%lu_Width_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Width of all hit pulses [ns]",ch);
    hWidth_AllHits[ch]          = diagDirOp.make<TH1D>(histName,histTitle,120,0,30);
    
    sprintf(histName,"%lu_AmplitudeVsWidth",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude of Michel candidate pulse [mV];Width [ns]",ch);
    hAmplitudeVsWidth[ch]          = diagDirOp.make<TH2D>(histName,histTitle,200,0.,200.,120,0,30);
    hAmplitudeVsWidth[ch]->SetOption("colz");
    
    sprintf(histName,"%lu_PromptPEVsWidth",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PEs];Pulse width [ns]",ch);
    hPromptPEVsWidth[ch]          = diagDirOp.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,30);
    hPromptPEVsWidth[ch]          ->SetOption("colz");
    
    sprintf(histName,"%lu_AmplitudeVsWidth_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Amplitude [ADC];Width [ns]",ch);
    hAmplitudeVsWidth_AllHits[ch]          = diagDirOp.make<TH2D>(histName,histTitle,120,0.,1200.,120,0,30);
    hAmplitudeVsWidth_AllHits[ch]->SetOption("colz");
    
    sprintf(histName,"%lu_PromptPEVsWidth_AllHits",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light [PEs];Pulse width [ns]",ch);
    hPromptPEVsWidth_AllHits[ch]          = diagDirOp.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,30);
    hPromptPEVsWidth_AllHits[ch]          ->SetOption("colz");
    sprintf(histName,"%lu_PromptPEVsAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PEs];Pulse amplitude [ADC]",ch);
    hPromptPEVsAmplitude[ch]          = diagDirOp.make<TH2D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch],120,0,600.);
    hPromptPEVsAmplitude[ch]          ->SetOption("colz");
    
    sprintf(histName,"%lu_MuPromptPEVsAmplitude",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of muon candidate pulse [PEs];Pulse amplitude [ADC]",ch);
    hMuPromptPEVsAmplitude[ch]          = diagDirOp.make<TH2D>(histName,histTitle,100,0.,500.,120,0,1200.);
    hMuPromptPEVsAmplitude[ch]          ->SetOption("colz");
   
    sprintf(histName,"%lu_PE_totalRaw_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu (raw integral), #DeltaT > %3.1f#mus, Michel shwr reco'd;Total light of Michel candidate pulse [PE]",ch,fdTcut/1000);
    hPE_totalRaw_dTcut_shwr[ch]     = diagDirOp.make<TH1D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch]);
    
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
  
    hElShowerEnergyVsPromptPE[ch]=tfs->make<TH2D>(Form("%lu_ElShowerEnergyVsPromptPE",ch),";Reconstructed Q-only energy [MeV];Prompt light [pe]",
      Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy,
      nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
      hElShowerEnergyVsPromptPE[ch]->SetOption("colz");
  
    sprintf(histName,"%lu_PE_total_TrueVsReco",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd;True Light [PE];Reco Light [PE]",ch,fdTcut/1000);
    hTrue_PE_total_TrueVsReco[ch] = diagDirOp.make<TH2D>(histName,histTitle,nbins_S[ch],x1_S[ch],x2_S[ch],nbins_S[ch],x1_S[ch],x2_S[ch]);
    hTrue_PE_total_TrueVsReco[ch]->SetOption("colz");
    
    // -------------------------------------------------------------
    // "Truth" PMT histograms 
    sprintf(histName,"%lu_True_PE_preTrig",ch);
    sprintf(histTitle,"Opdet %lu: True number of detected photons (prior to trigger efficiency cut);True detected light [PE]",ch);
    hTrue_PE_preTrig[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_TrueLY",ch);
    sprintf(histTitle,"True LY (opdet %lu);Light Yield [PE/MeV]",ch);
    hTrue_LightYield[ch]        = truthDir.make<TH1D>(histName,histTitle,300,0.,30.); 
    sprintf(histName,"%lu_True_PE_prompt",ch);
    sprintf(histTitle,"Opdet %lu;Prompt light of Michel candidate pulse [PE]",ch);
    hTrue_PE_prompt[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S100[ch],x1_S100[ch],x2_S100[ch]);
    sprintf(histName,"%lu_True_PE_total",ch);
    sprintf(histTitle,"Opdet %lu;Total light of Michel candidate pulse [PE]",ch);
    hTrue_PE_total[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_True_PE_total_dTcut_shwr",ch);
    sprintf(histTitle,"Opdet %lu, #DeltaT > %3.1f#mus, Michel shwr reco'd (true);True detected light of Michel pulse [PE]",ch,fdTcut/1000);
    hTrue_PE_total_dTcut_shwr[ch]     = truthDir.make<TH1D>(histName,histTitle,nbins_S[ch], x1_S[ch], x2_S[ch]);
    sprintf(histName,"%lu_True_ResPEch",ch);
    sprintf(histTitle,"Opdet %lu (#DeltaT > %3.1f#mus, shwr), Photoelectron resolution of Michel candidate pulse;(reco-true)/true",ch,fdTcut/1000);
    hTrue_ResPEch[ch]       = truthDir.make<TH1D>(histName,histTitle,150,-1.5,1.5);
    
    hTrue_dT_vs_MuContamRes[ch] = truthDir.make<TH2D>(Form("%lu_True_dT_vs_MuContamRes",ch),";#DeltaT [ns];(reco-true)/true", 75,0.,7500, 100, -1.0,1.0);
    hTrue_dT_vs_MuContamRes[ch]->SetOption("colz");

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
  for(int i=0; i<2; i++){
    hNumHits[i]                 = diagDirChHit.make<TH1D>(Form("pl%d_NumHits",i),          Form("Number of wire hits in event (plane %d)",i),200,0,1000);
    hHitT[i]                    = diagDirChHit.make<TH1D>(Form("pl%d_HitT",i),             Form("Calculated drift hit (plane %d);X [cm]",i),400,-50.,350.);
    hHitX[i]                    = diagDirChHit.make<TH1D>(Form("pl%d_HitX",i),             Form("Calculated drift position from hit time (plane %d);X [cm]",i),325,-10.,55.);
    hHitRMS[i]                  = diagDirChHit.make<TH1D>(Form("pl%d_HitRMS",i),           Form("Hit RMS, ALL hits (plane %d)",i),120,0.,30.);
    hHitAmp[i]                  = diagDirChHit.make<TH1D>(Form("pl%d_HitAmp",i),           Form("Hit amplitude, ALL hits (plane %d);ADC",i),300,0,300);
    hHitRMS_ElShower[i]         = diagDirChHit.make<TH1D>(Form("pl%d_HitRMS_ElShower",i),  Form("Hit RMS (plane %d)",i),120,0.,30.);
    hHitAmp_ElShower[i]         = diagDirChHit.make<TH1D>(Form("pl%d_HitAmp_ElShower",i),  Form("Hit amplitude (plane %d);ADC",i),300,0,300);
    hHitIntegral[i]             = diagDirChHit.make<TH1D>(Form("pl%d_HitIntegral",i),      Form("Hit integral (plane %d);ADC",i),150,0, 15000);
    hHitIntegral_ElShower[i]    = diagDirChHit.make<TH1D>(Form("pl%d_HitIntegral_ElShower",i),Form("Hit integral (plane %d);ADC",i),150,0, 15000);
    hHitAmp_vs_HitRMS[i]        = diagDirChHit.make<TH2D>(Form("pl%d_HitAmp_vs_HitRMS",i), Form("Plane %d;Hit peak amplitude [ADC];Hit width [ADC]",i),250,0.,250.,120,0.,30.);
      hHitAmp_vs_HitRMS[i]      ->SetOption("colz"); 
    hHitAmp_vs_HitIntegral[i]   = diagDirChHit.make<TH2D>(Form("pl%d_HitAmp_vs_HitIntegral",i),Form("Plane %d;Hit amplitude [ADC];Hit integral [ADC]",i),250,0.,250.,250,0.,5000);
      hHitAmp_vs_HitIntegral[i]  ->SetOption("colz");
    hHitIntegral_vs_HitRMS[i]   = diagDirChHit.make<TH2D>(Form("pl%d_HitIntegral_vs_HitRMS",i),Form("Plane %d;Hit integral [ADC];Hit width [ADC]",i),250,0.,5000.,120,0.,30.);
      hHitIntegral_vs_HitRMS[i] ->SetOption("colz"); 
    hHitCharge[i]               = diagDirChHit.make<TH1D>(Form("pl%d_HitCharge",i),        Form("Hit charge after lifetime correction (plane %d);Q (e-)",i),150,0, 300e3);
    hHitCharge_ElShower[i]               = diagDirChHit.make<TH1D>(Form("pl%d_HitCharge_ElShower",i),        Form("Hit charge after lifetime correction (plane %d);Q (e-)",i),150,0, 300e3);
  }
  hTotalCharge          = diagDir.make<TH1D>("TotalCharge","Total charge on collection plane;Charge Q [e-]",100,0,20000000);
 
  /* 
  hHitIntegral_vs_SummedADC_ElShower = diagDir.make<TH2D>("HitIntegral_vs_SummedADC_ElShower","Electron shower +track hits;Hit integral [ADC];Summed ADC [ADC]",
    220,-1000,10000,
    220,-1000,10000);
    hHitIntegral_vs_SummedADC_ElShower->SetOption("colz");
  
  hHitIntegral_vs_SummedADC_MuTrk = diagDir.make<TH2D>("HitIntegral_vs_SummedADC_MuTrk","Muon track hits;Hit integral [ADC];Summed ADC [ADC]",
    220,-1000,10000,
    220,-1000,10000);
    hHitIntegral_vs_SummedADC_MuTrk->SetOption("colz");
  
  hHitIntegral_vs_DiffSummedADC_ElShower = diagDir.make<TH2D>("HitIntegral_vs_DiffSummedADC_ElShower","Displaced electron shower hits;Hit integral [ADC];%#Delta (integral - sum ADC)",
    200,0.,10000,
    300,-0.50,1.00);
    hHitIntegral_vs_DiffSummedADC_ElShower->SetOption("colz");
  hHitIntegral_vs_DiffSummedADC_ElTrk = diagDir.make<TH2D>("HitIntegral_vs_DiffSummedADC_ElTrk","Electron track hits;Hit integral [ADC];%#Delta (integral - sum ADC)",
    200,0.,10000,
    300,-0.50,1.0);
    hHitIntegral_vs_DiffSummedADC_ElTrk->SetOption("colz");
  
  hHitIntegral_vs_DiffSummedADC_MuTrk = diagDir.make<TH2D>("HitIntegral_vs_DiffSummedADC_MuTrk","Muon track hits;Hit integral [ADC];%#Delta (integral - sum ADC)",
    200,0.,10000,
    300,-0.50,1.00);
    hHitIntegral_vs_DiffSummedADC_MuTrk->SetOption("colz");
  
  hHitAmp_vs_DiffSummedADC_ElShower = diagDir.make<TH2D>("HitAmp_vs_DiffSummedADC_ElShower","Displaced electron shower hits;Hit amplitude [ADC];%#Delta (integral - sum ADC)",
    250,0.,250.,
    300,-0.50,1.00);
    hHitAmp_vs_DiffSummedADC_ElShower->SetOption("colz");
  
  hHitAmp_vs_DiffSummedADC_ElTrk = diagDir.make<TH2D>("HitAmp_vs_DiffSummedADC_ElTrk","Electron track hits;Hit amplitude [ADC];%#Delta (integral - sum ADC)",
    250,0.,250.,
    300,-0.50,1.00);
    hHitAmp_vs_DiffSummedADC_ElTrk->SetOption("colz");
  
  hHitAmp_vs_DiffSummedADC_MuTrk = diagDir.make<TH2D>("HitAmp_vs_DiffSummedADC_MuTrk","Muon track hits;Hit amplitude [ADC];%#Delta (integral - sum ADC)",
    250,0.,250,
    300,-0.50,1.00);
    hHitAmp_vs_DiffSummedADC_MuTrk->SetOption("colz");
  
  hTrackInclineAngle_vs_DiffSummedADC = diagDir.make<TH2D>("TrackInclineAngle_vs_DiffSummedADC",";Track inclination angle [deg];%#Delta (integral - sum ADC)",
    100,0.,90.,
    150,-0.15,0.15);
  hTrackInclineAngle_vs_DiffSummedADC->SetOption("colz");
 */
  
  // ===================================================================
  // Reco track diagnostics
  hNumTrack             = diagDir.make<TH1D>("NumTrack","Number of tracks per event",30,0,30);
  hTrackNode_X          = diagDir.make<TH1D>("TrackNode_X","Distribution of reco track nodes;X [cm]", 600, -5., 55.);
  hTrackNode_X_MuTrk    = diagDir.make<TH1D>("TrackNode_X_MuTrk","Distribution of reco track vertex, tagged mu trk (corrected for decay time);X [cm]", 600, -5., 55.);
  hTrackNode_Y          = diagDir.make<TH1D>("TrackNode_Y","Distribution of reco track nodes;Y [cm]", 500, -25., 25.);
  hTrackNode_Z          = diagDir.make<TH1D>("TrackNode_Z","Distribution of reco track nodes;Z [cm]", 1000, -5., 95.);
  hTrackNode_ZX         = diagDir.make<TH2D>("TrackNode_ZX","Reco track nodes (top-down view);Z [cm];X [cm]",200, -5., 95., 120, -5., 55.);
  hTrackNode_ZY         = diagDir.make<TH2D>("TrackNode_ZY","Reco track nodes (side view);Z [cm];Y [cm]",200, -5., 95., 100, -25., 25.);
  hTrackNode_ZX         ->SetOption("colz");
  hTrackNode_ZY         ->SetOption("colz");
  hTrackZenithAngle          = diagDir.make<TH1D>("TrackZenithAngle","Track zenith angle;Zenith angle [deg]",180.,0.,180.);
  hTrackLength          = diagDir.make<TH1D>("TrackLength","Track lengths;L [cm]",240,0.,120.);
  hTrackLength_vs_ZenithAngle          = diagDir.make<TH2D>("TrackLength_vs_ZenithAngle","Track length vs. zenith angle;L [cm];Zenith angle [deg]",180,0.,90.,180.,0.,180.);
  hTrackLength_vs_ZenithAngle          ->SetOption("colz");
  hNumTrackStopping     = diagDir.make<TH1D>("NumTrackStopping","Number of identified stopping tracks per event",20,0,20);
  hNumTrackCrossing     = diagDir.make<TH1D>("NumTrackCrossing","Number of identified passing tracks per event",20,0,20);
  hNumTrackContained    = diagDir.make<TH1D>("NumTrackContained","Number of identified contained tracks per event",20,0,20);
  hMuTrackZenithAngle   = diagDir.make<TH1D>("MuTrackZenithAngle","Zenith angle of stopping #mu track;Zenith angle [deg]",180.,0.,180.);
  hMuTrackInclineAngle   = diagDir.make<TH1D>("MuTrackInclineAngle","Inclination angle of stopping #mu track (rel. to E-field);Incline angle [deg]",180.,0.,180.);
  hMuTrackInclineAngle_dEdx   = dEdxDir.make<TH1D>("MuTrackInclineAngle_dEdx","Inclination angle of stopping #mu track (rel. to E-field);Incline angle [deg]",180.,0.,180.);
 
  float r_max = 30;
  float r_bins = 75; 
  float dEdx_max = 25.;
  int dEdx_bins = 250;
  float dADCdx_max = 25000.;
  int dADCdx_bins = 250;
  float dQdx_max = 400e3;
  int dQdx_bins = 200;

  hMuResRangePitch[1]       = dEdxDir.make<TH1D>("MuResRangePitch","Pitch of hits within tracks used in residual range study (coll plane);Pitch [cm]",160,0.,8.);
  hMuResRangePitch[0]       = dEdxDir.make<TH1D>("MuResRangePitch_Pl0","Pitch of hits within tracks used in residual range study (ind plane);Pitch [cm]",160,0.,8.);

  //hMuResRangeVsdEdxHyp[1]     = dEdxDir.make<TH2D>("MuResRangeVsdEdxHyp","(dE/dx)_{hyp} vs. residual range of tagged stopping #mu track (coll plane);Residual range [cm];(dE/dx)_{hyp} [MeV/cm]",
  //  r_bins,0.,r_max,  dEdx_bins*10,0.,dEdx_max);
  hMuResRangeVsdEdx[1]     = dEdxDir.make<TH2D>("MuResRangeVsdEdx","Residual range vs. dE/dx of tagged stopping #mu track (coll plane);Residual range [cm];dE/dx [MeV/cm]",
    r_bins,0.,r_max,  dEdx_bins,0.,dEdx_max);
  hMuResRangeVsdQdx[1]     = dEdxDir.make<TH2D>("MuResRangeVsdQdx","Residual range vs. dQ/dx of tagged stopping #mu track (coll plane);Residual range [cm];dQ/dx [e-/cm]",
    r_bins,0.,r_max,  dQdx_bins,0.,dQdx_max);
  hMuResRangeVsdADCdx[1]     = dEdxDir.make<TH2D>("MuResRangeVsdADCdx","Residual range vs. dQ/dx of tagged stopping #mu track (coll plane);Residual range [cm];dQ/dx [ADC/cm]",
    r_bins,0.,r_max,  dADCdx_bins,0.,dADCdx_max);
//  hMuResRangeVsdEdxHyp[0]     = dEdxDir.make<TH2D>("MuResRangeVsdEdxHyp_Pl0","(dE/dx)_{hyp} vs. residual range of tagged stopping #mu track (ind plane);Residual range [cm];(dE/dx)_{hyp} [MeV/cm]",
//    r_bins,0.,r_max,  dEdx_bins*10,0.,dEdx_max);
  hMuResRangeVsdEdx[0]     = dEdxDir.make<TH2D>("MuResRangeVsdEdx_Pl0","Residual range vs. dE/dx of tagged stopping #mu track (ind plane);Residual range [cm];dE/dx [MeV/cm]",
    r_bins,0.,r_max,  dEdx_bins,0.,dEdx_max);
  hMuResRangeVsdQdx[0]     = dEdxDir.make<TH2D>("MuResRangeVsdQdx_Pl0","Residual range vs. dQ/dx of tagged stopping #mu track (ind plane);Residual range [cm];dQ/dx [e-/cm]",
    r_bins,0.,r_max,  dQdx_bins,0.,dQdx_max);
  hMuResRangeVsdADCdx[0]     = dEdxDir.make<TH2D>("MuResRangeVsdADCdx_Pl0","Residual range vs. dQ/dx of tagged stopping #mu track (ind plane);Residual range [cm];dQ/dx [ADC/cm]",
    r_bins,0.,r_max,  dADCdx_bins,0.,dADCdx_max);
  /*hMuClusterResRangeVsdEdx[1]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdEdx","Residual range vs. dE/dx of tagged #mu cluster (coll. plane);Residual range [cm];dE/dx [MeV/cm]",
    r_bins,0.,r_max,  dEdx_bins,0.,dEdx_max);
  hMuClusterResRangeVsdQdx[1]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdQdx","Residual range vs. dQ/dx of tagged #mu cluster (coll. plane);Residual range [cm];dQ/dx [e-/cm]",
    r_bins,0.,r_max,  dQdx_bins,0.,dQdx_max);
  hMuClusterResRangeVsdADCdx[1]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdADCdx","Residual range vs. dQ/dx (ADC) of tagged #mu cluster (coll. plane);Residual range [cm];dQ/dx [ADC/cm]",
    r_bins,0.,r_max,  dADCdx_bins,0.,dADCdx_max);
  hMuClusterResRangeVsdEdx[0]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdEdx_Pl0","Residual range vs. dE/dx of tagged #mu cluster (ind. plane);Residual range [cm];dE/dx [MeV/cm]",
    r_bins,0.,r_max,  dEdx_bins,0.,dEdx_max);
  hMuClusterResRangeVsdQdx[0]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdQdx_Pl0","Residual range vs. dQ/dx of tagged #mu cluster (ind. plane);Residual range [cm];dQ/dx [e-/cm]",
    r_bins,0.,r_max,  dQdx_bins,0.,dQdx_max);
  hMuClusterResRangeVsdADCdx[0]     = dEdxDir.make<TH2D>("MuClusterResRangeVsdADCdx_Pl0","Residual range vs. dQ/dx (ADC) of tagged #mu cluster (ind. plane);Residual range [cm];dQ/dx [ADC/cm]",
    r_bins,0.,r_max,  dADCdx_bins,0.,dADCdx_max);
  */
  //hMu_dEdxHyp_vs_dQdx             = dEdxDir.make<TH2D>("Mu_dEdxHyp_vs_dQdx",";(dE/dx)_{hyp} [MeV/cm];dQ/dx [e-/cm]",  150,0.,15., 200,0,4e5);
  //hMu_dEdxHyp_vs_dQdx             ->SetOption("colz");
  
  // Angle bins
  for(size_t i=0; i<4; i++){
    hMuResRangeVsdEdx_AngleBins[i] = dEdxDir.make<TH2D>(
      Form("MuResRangeVsdEdx_AngleBin_%2.0f",fAngleBins[i]),"Residual range vs. dE/dx of tagged stopping #mu track (coll plane);Residual range [cm];dE/dx [MeV/cm]",
      r_bins,0.,r_max,  dEdx_bins,0.,dEdx_max);
      hMuResRangeVsdEdx_AngleBins[i]->SetOption("colz");
    hMuResRangeVsdQdx_AngleBins[i] = dEdxDir.make<TH2D>(
      Form("MuResRangeVsdQdx_AngleBin_%2.0f",fAngleBins[i]),"Residual range vs. dQ/dx of tagged stopping #mu track (coll plane);Residual range [cm];dQ/dx [e-/cm]",
      r_bins,0.,r_max,  dQdx_bins,0.,dQdx_max);
      hMuResRangeVsdQdx_AngleBins[i]->SetOption("colz");
    hMuResRangeVsdADCdx_AngleBins[i] = dEdxDir.make<TH2D>(
      Form("MuResRangeVsdADCdx_AngleBin_%2.0f",fAngleBins[i]),"Residual range vs. dADC/dx of tagged stopping #mu track (coll plane);Residual range [cm];dADC/dx [ADC/cm]",
      r_bins,0.,r_max,  dADCdx_bins,0.,dADCdx_max);
      hMuResRangeVsdADCdx_AngleBins[i]->SetOption("colz");
  }
  
  
  
  //hMuResRangeVsdEdxHyp[0]    ->SetOption("colz");
  hMuResRangeVsdEdx[0]     ->SetOption("colz");
  hMuResRangeVsdADCdx[0]     ->SetOption("colz");
  hMuResRangeVsdQdx[0]     ->SetOption("colz");
//  hMuResRangeVsdEdxHyp[1]    ->SetOption("colz");
  hMuResRangeVsdEdx[1]     ->SetOption("colz");
  hMuResRangeVsdADCdx[1]     ->SetOption("colz");
  hMuResRangeVsdQdx[1]     ->SetOption("colz");
  /*
  hMuClusterResRangeVsdADCdx[0]     ->SetOption("colz");
  hMuClusterResRangeVsdEdx[1]     ->SetOption("colz");
  hMuClusterResRangeVsdQdx[1]     ->SetOption("colz");
  hMuClusterResRangeVsdADCdx[1]     ->SetOption("colz");
  hMuClusterResRangeVsdQdx[0]     ->SetOption("colz");
  hMuClusterResRangeVsdEdx[0]     ->SetOption("colz");
  */


  hTrue_MuResRangeVsdEdx     = dEdxDir.make<TH2D>("True_MuResRangeVsdEdx","Residual range vs. dE/dx of tagged stopping muon track (dL = 0.4cm);Residual range [cm];dE/dx [MeV/cm]",
    r_bins,0.,r_max,dEdx_bins*2,0.,dEdx_max);
    hTrue_MuResRangeVsdEdx     ->SetOption("colz");
  hTrue_MuResRangeVsdEdx_HR   = dEdxDir.make<TH2D>("True_MuResRangeVsdEdx_HR","Residual range vs. dE/dx of tagged stopping muon track (dL = 0.1cm);Residual range [cm];dE/dx [MeV/cm]",
    r_bins*4,0.,r_max,dEdx_bins*2,0.,dEdx_max);
    hTrue_MuResRangeVsdEdx_HR   ->SetOption("colz");
  
  // ======================================================================
  // 2D clustering and shower reco performance diagnostics  
  hClusterHitSeparation = diagDir.make<TH1D>("ClusterHitSeparation","Nearest neighbor hit distance in proximity clustering",50,0.,5.);
  hTruncMeanSlope       = diagDir.make<TH1D>("TruncMeanSlope","Slope of truncated mean dQ/ds at Bragg peak;e-/cm",225,-5000.,40000);
  hTruncMeanSlopeChi2   = diagDir.make<TH1D>("TruncMeanSlopeChi2","Chi2 of slope",100,0.,20.);
  hFracDropoffAtBnd     = diagDir.make<TH1D>("FracDropoffAtBnd","Fractional change in charge at cluster boundary",140,-1.2,1.2);
  hCovAtBnd             = diagDir.make<TH1D>("CovAtBnd","Minimum covariance (linearity) in neighborhood of bnd",100,0.,1.);
  hClusterSize          = diagDir.make<TH1D>("ClusterSize","Cluster size (num of hits)",300,0,300);
  hElHitChargeRatio     = diagDir.make<TH2D>("ElHitChargeRatio",";Distance from muon boundary [hits];Q/Q_{bnd}",30,0,30,100,0.,2.0);
    hElHitChargeRatio   ->SetOption("colz");
  hLeftoverHits         = diagDir.make<TH1D>("LeftoverHits","Number of non-muon hits not included in electron shower",60,0,60);
  hElShowerFrac        = diagDir.make<TH1D>("ElShowerFrac","Fraction of total non-muon hits included in electron shower",120,0,1.2);
  hElTrackFrac          = diagDir.make<TH1D>("ElTrackFrac","Charge fraction of contiguous trk-like hits in shower",110,0.,1.1);
  hMuClusterSize        = diagDir.make<TH1D>("MuClusterSize","Number of #mu hits in cluster",200,0,200);
  hElClusterSize        = diagDir.make<TH1D>("ElClusterSize","Number of electron hits in cluster",200,0,200);
  hElShowerSize        = diagDir.make<TH1D>("ElShowerSize","Number of hits in shower (collection plane)",200,0,200);
  hExtraHits            = diagDir.make<TH1D>("ExtraHits","Extra (unclustered) hits in close proximity to all clustered hits",50,0,50);
  hDistMaxTdQdsToMaxQ   = diagDir.make<TH1D>("DistMaxTdQdsToMaxQ","Distance btwn max trunc. mean Q/ds and max Q/ds;Hit_{Max dQ/ds} - Hit_{Max trunc. dQ/ds} [hits]",30,-15,15);
  hDistMaxTdQdsToMaxDrop   = diagDir.make<TH1D>("DistMaxTdQdsToMaxDrop","Distance btwn max trunc. mean Q/ds and max Q/ds drop;Hit_{Max dQ/ds} - Hit_{Max trunc. dQ/ds} [hits]",30,-15,15);
  hFracMuHitsLinear     = diagDir.make<TH1D>("FracMuHitsLinear","Fraction of #mu hits > covariance thresh",120,0,1.2);
  hMuAveLinearity     = diagDir.make<TH1D>("MuAveLinearity","Average muon linearity",110,0,1.2);
  hMuClusterHitsEndFit  = diagDir.make<TH1D>("MuClusterHitsEndFit","Number of #mu hits used for direction calcuation",30,0,30);
  hDecayAngle2D         = diagDir.make<TH1D>("DecayAngle2D","Michel decay angle rel. to terminal #mu direction;Angle [deg]",90,0.,180.);
  hElShowerDepDist   = diagDir.make<TH1D>("ElShowerDepDist","Distance of Michel shower depositions;Distance in W-X space [cm]",100,0.,40.);
  hElShowerDepAngle2D   = diagDir.make<TH1D>("ElShowerDepAngle2D","Opening angle of Michel shower depositions;Angle [deg]",90,0.,180.);
  hShowerSizeCompare   = diagDir.make<TH2D>("ShowerSizeCompare","Electron shower size on both planes;Collection Plane;Induction Plane",100,0,100,100,0,100);
  hShowerSizeCompare    ->SetOption("colz");
  hElClusterSizeCompare   = diagDir.make<TH2D>("ElClusterSizeCompare","Electron cluster size on both planes;Collection Plane;Induction Plane",100,0,100,100,0,100);
  hElClusterSizeCompare    ->SetOption("colz");
  hHitTimeDiff          = diagDir.make<TH1D>("HitTimeDiff","Hit time differences between Michel shower hits in collection/induction planes;T(coll.) - T(ind.) [#mus]",200,-10.,10.);
  hHitXDiff          = diagDir.make<TH1D>("HitXDiff","Hit X differences between Michel shower hits in collection/induction planes;X(coll.) - X(ind.) [cm]",200,-10.,10.);
  
  // ======================================================================
  // 3D shower reco performance diagnostics  
  hMuEndHitTimeDiff     = diagDir.make<TH1D>("MuEndHitTimeDiff","Hit time difference btwn candidate #mu end hits in collection/induction plane clusters;T_{coll} - T_{ind} [#mus]",300,-15.,15.);
  hMuEndDiff3D          = diagDir.make<TH1D>("MuEndDiff3D","Spatial separation between #mu endpts from 3D track vs. clustering;#Delta d [cm]",200,0.,10.);
  hMuEnd3D_X            = diagDir.make<TH1D>("MuEnd3D_X","3D #mu endpoint;X [cm]",11, -5., 50.);
  hMuEnd3D_Y            = diagDir.make<TH1D>("MuEnd3D_Y","3D #mu endpoint;Y [cm]",10, -25., 25.);
  hMuEnd3D_Z            = diagDir.make<TH1D>("MuEnd3D_Z","3D #mu endpoint;Z [cm]",20, -5., 95.);
  hMuEnd3D_ZX         = diagDir.make<TH2D>("MuEnd3D_ZX","3D #mu endpoint (top-down view);Z [cm];X [cm]",100, -5., 95., 60, -5., 55.);
  hMuEnd3D_ZY         = diagDir.make<TH2D>("MuEnd3D_ZY","3D #mu endpoint (side view);Z [cm];Y [cm]",100, -5., 95.,40, -25., 25.);
  hMuEnd3D_ZX         ->SetOption("colz");
  hMuEnd3D_ZY         ->SetOption("colz");
  hElShowerCentroid_X            = diagDir.make<TH1D>("ElShowerCentroid_X","Shower centroid;X [cm]",11, -5., 50.);
  hElShowerCentroid_Y            = diagDir.make<TH1D>("ElShowerCentroid_Y","Shower centroid;Y [cm]",10, -25., 25.);
  hElShowerCentroid_Z            = diagDir.make<TH1D>("ElShowerCentroid_Z","Shower centroid;Z [cm]",20, -5., 95.);
  hElShowerCentroid_ZX         = diagDir.make<TH2D>("ElShowerCentroid_ZX","Shower centroid (top-down view);Z [cm];X [cm]",100, -5., 95., 60, -5., 55.);
  hElShowerCentroid_ZY         = diagDir.make<TH2D>("ElShowerCentroid_ZY","Shower centroid (side view);Z [cm];Y [cm]",100, -5., 95.,40, -25., 25.);
  hElShowerCentroid_ZX         ->SetOption("colz");
  hElShowerCentroid_ZY         ->SetOption("colz");
  hFracHits3D            = diagDir.make<TH1D>("FracHits3D","Fraction of Michel shower hits made into 3D pts",120,0,1.2);
  hNumPts3D             = diagDir.make<TH1D>("NumPts3D","Number of 3D Michel shower hits",100,0,100);
  hProjDist3DToWires    = diagDir.make<TH1D>("ProjDist3DToWires","Projected 3D distance of shower to wires;L [cm]",300,0.,300.);
  hProjDist3DToWiresRes = diagDir.make<TH1D>("ProjDist3DToWiresRes","L - L_true / L_true",100,-1.0,1.0);
  hProjPtOnWires        = diagDir.make<TH2D>("ProjPtOnWires","Projected intersection point of Michel direction and wireplanes;Z [cm];Y [cm]",75,-30., 120.,35,-35.,35.);
    hProjPtOnWires      ->SetOption("colz");

  // ======================================================================
  // Michel electron charge and energy
  hElShowerCharge         = tfs->make<TH1D>("ElShowerCharge","Reco charge of Michel electron;Reconstructed Q [e-]",150,0,3000000);
  hElShowerCharge3D         = tfs->make<TH1D>("ElShowerCharge3D","Reco charge of Michel electron;Reconstructed Q [e-]",150,0,3000000);
  hElShowerVis            = tfs->make<TH1D>("ElShowerVis","Visibility of reco'd Michel shower;Fractional photon visibility",100,0.,0.002);
  hElShowerPhotons        = tfs->make<TH1D>("ElShowerPhotons","Photons produced by Michel electron;Reconstructed L [#gamma]",150,0,3000000);
  hElShowerEnergy         = tfs->make<TH1D>("ElShowerEnergy","Reco Michel Electron Shower Energy;Reconstructed energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hElShowerEnergyQL       = tfs->make<TH1D>("ElShowerEnergyQL","Reco Michel Electron Shower Energy (Q+L);Reconstructed energy [MeV]",Nbins_michelEnergy, x1_michelEnergy, x2_michelEnergy);
  hElTrackCharge         = tfs->make<TH1D>("ElTrackCharge","Reco charge of Michel electron (track-like components);Reconstructed Q [e-]",150,0,3000000);
//  hElTrackdEdx            = tfs->make<TH1D>("ElTrackdEdx","Averaged dE/dx along 3D electron track;dE/dx [MeV/cm]",100,0.,10.);
  hQoverL                 = tfs->make<TH1D>("QoverL","Q/L for Michel electrons;Q/L",80,0.,4.);
  hRecomb                 = tfs->make<TH1D>("Recomb","Recombination (e- survival fraction) calculated from Q/L;Electron survival fraction, R",140,-0.2,1.2);

  // ======================================================================
  // Crossing muon calibration plots 
  hTrue_ResE_CrsMu        = crsMuDir.make<TH1D>("True_ResE_CrsMu","Energy resolution of crossing mu track;(E_{reco}-E_{true})/E_{true}", 320,-0.8,0.8);
  hTrue_ResQ_CrsMu[0]        = crsMuDir.make<TH1D>("True_ResQ_CrsMu_Pl0","Resolution of Q for crossing #mu (ind. plane);(Q_{reco}-Q_{true}) / Q_{true}",320,-0.8,0.8);
  hTrue_ResQ_CrsMu[1]        = crsMuDir.make<TH1D>("True_ResQ_CrsMu","Resolution of Q for crossing #mu (coll. plane);(Q_{reco}-Q_{true}) / Q_{true}",320,-0.8,0.8);
  hTrue_ResQ_CrsMuTrk[0]        = crsMuDir.make<TH1D>("True_ResQ_CrsMuTrk_Pl0","Resolution of Q for crossing #mu (ind. plane);(Q_{reco}-Q_{true}) / Q_{true}",320,-0.8,0.8);
  hTrue_ResQ_CrsMuTrk[1]        = crsMuDir.make<TH1D>("True_ResQ_CrsMuTrk","Resolution of Q for crossing #mu (coll. plane);(Q_{reco}-Q_{true}) / Q_{true}",320,-0.8,0.8);
  hTrue_ResL_CrsMu        = crsMuDir.make<TH1D>("True_ResL_CrsMu","Photon resolution of crossing #mu;(reco-true)/true",320,-0.8,0.8);
  hCrsTrackdT_ACP         = crsMuDir.make<TH1D>("CrsTrackdT_ACP","Time difference between nodes for anode-cathode piercing tracks;#DeltaT [#mus]",100,270,370);
  hCrsTrackNode_ZX        = crsMuDir.make<TH2D>("CrsTrackNode_ZX","Reco track nodes (top-down view);Z [cm];X [cm]",200, -5., 95., 120, -5., 55.);
  hCrsTrackNode_ZY        = crsMuDir.make<TH2D>("CrsTrackNode_ZY","Reco track nodes (side view);Z [cm];Y [cm]",200, -5., 95., 100, -25., 25.);
  hCrsTrackNode_ZX        ->SetOption("colz");
  hCrsTrackNode_ZY        ->SetOption("colz");
  hCrsMuLength            = crsMuDir.make<TH1D>("CrsMuLength","Crossing muon length;L [cm]",100,0.,100.);
  hTrue_CrsMuLength       = crsMuDir.make<TH1D>("True_CrsMuLength","Crossing muon length;L [cm]",100,0.,100.);
  hCrsMuHitFrac           = crsMuDir.make<TH1D>("CrsMuHitFrac","Fraction of collection plane hits in #mu track",101,0.,1.01);
  hCrsMuPitch[0]             = crsMuDir.make<TH1D>("CrsMuPitch_Pl0","Crossing #mu dL (induction plane)",200.,0.,10.);
  hCrsMuPitch[1]             = crsMuDir.make<TH1D>("CrsMuPitch","Crossing #mu dL",200.,0.,10.);
  hQoverLCrsMuTrk         = crsMuDir.make<TH1D>("QoverL_CrsMuTrk","Q/L for crossing muon tracks;Q/L",80,0.,4.);
  hRecombCrsMuTrk         = crsMuDir.make<TH1D>("Recomb_CrsMuTrk","Recombination calculated from Q/L (crossing muons);Electron survival fraction, R",140,-0.2,1.2);
  hCrsMuADCPerCm          = crsMuDir.make<TH1D>("CrsMuADCPerCm","ADC/cm for crossing #mu (coll. plane);Charge deposited per unit length [ADC/cm]",150,0.,15000);
  hCrsMuADCPerCm_Pl0      = crsMuDir.make<TH1D>("CrsMuADCPerCm_Pl0","ADC/cm for crossing #mu (ind. plane);Charge deposited per unit length [ADC/cm]",150,0.,15000);
  hCrsMuQPerCm            = crsMuDir.make<TH1D>("CrsMuQPerCm","Q/cm for crossing #mu (coll. plane);Charge deposited per unit length [e-/cm]",140,0.,14e4);
  hCrsMuQPerCm_Pl0        = crsMuDir.make<TH1D>("CrsMuQPerCm_Pl0","Q/cm for crossing #mu (ind. plane);Charge deposited per unit length [e-/cm]",140,0.,14e4);
  hCrsMuLPerCm            = crsMuDir.make<TH1D>("CrsMuLPerCm","L/cm for crossing muons;Photons produced per unit length [#gamma/cm]",140,0.,14e4);
  hTrue_CrsMuQPerCm       = crsMuDir.make<TH1D>("True_CrsMuQPerCm","True Q/cm for crossing #mu;Charge deposited per unit length [e-/cm]",140,0.,14e4);
  hTrue_CrsMuLPerCm       = crsMuDir.make<TH1D>("True_CrsMuLPerCm","True L/cm for crossing #mu;Photons produced per unit length [#gamma/cm]",140,0.,14e4);
  hCrsMu_dADCdx[1]             = crsMuDir.make<TH1D>("CrsMu_dADCdx","Hit-by-hit ADC/dx for crossing #mu tracks (col. plane);ADC/dx [e-/cm]",200,0.,20000);
  hCrsMu_dQdx[1]             = crsMuDir.make<TH1D>("CrsMu_dQdx","Hit-by-hit dQ/dx for crossing #mu tracks (col. plane);dQ/dx [e-/cm]",100,0.,2e5);
  hCrsMu_dEdx[1]             = crsMuDir.make<TH1D>("CrsMu_dEdx","Hit-by-hit dE/dx for crossing #mu tracks (col. plane);dE/dx [MeV/cm]",80,0.,8.);
  hCrsMu_dADCdx[0]             = crsMuDir.make<TH1D>("CrsMu_dADCdx_Pl0","Hit-by-hit ADC/dx for crossing #mu tracks (ind. plane);ADC/dx [e-/cm]",200,0.,20000);
  hCrsMu_dQdx[0]             = crsMuDir.make<TH1D>("CrsMu_dQdx_Pl0","Hit-by-hit dQ/dx for crossing #mu tracks (ind. plane);dQ/dx [e-/cm]",100,0.,2e5);
  hCrsMu_dEdx[0]             = crsMuDir.make<TH1D>("CrsMu_dEdx_Pl0","Hit-by-hit dE/dx for crossing #mu tracks (ind. plane);dE/dx [MeV/cm]",80,0.,8.);
   
  // ======================================================================
  // Truth plots
  hTrue_ResQ_Mu           = truthDir.make<TH1D>("True_ResQ_Mu","Resolution of Q for stopping #mu;(Q_{reco}-Q_{true}) / Q_{true}",200,-1.0,1.0);
  hTrue_ResQ              = truthDir.make<TH1D>("True_ResQ","Resolution of Q;(Q_{reco}-Q_{true}) / Q_{true}",200,-1.0,1.0);
  hTrue_ResQ_Shwr3D              = truthDir.make<TH1D>("True_ResQ_Shwr3D","Resolution of Q;(Q_{reco}-Q_{true}) / Q_{true}",200,-1.0,1.0);
  hTrue_ResL              = truthDir.make<TH1D>("True_ResL","Photon resolution of Michel candidate;(reco-true)/true",200,-1.,1.);
  hTrue_ResVis            = truthDir.make<TH1D>("True_ResVis","Resolution of reconstructed Michel shower vis.;(reco-true)/true",200,-1.0,1.0);
  hTrue_ResPE             = truthDir.make<TH1D>("True_ResPE","PE resolution of Michel candidate;(reco-true)/true",200,-1.0,1.0);
  hTrue_ResPEPrompt       = truthDir.make<TH1D>("True_ResPEPrompt","Prompt PE resolution of Michel candidate;(reco-true)/true",200,-1.0,1.0);
  hTrue_ResPE0            = truthDir.make<TH1D>("True_ResPE0","Smearing due to limited integration window;(S_{true}-S_{true}^{0})/S_{true}^{0}",125,-0.20,0.05);
  hTrue_EnergyRes     = truthDir.make<TH1D>("True_EnergyRes","Michel electron deposited energy resolution;(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyResTrk     = truthDir.make<TH1D>("True_EnergyResTrk","Michel electron deposited energy resolution (electron track only);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyResQL     = truthDir.make<TH1D>("True_EnergyResQL","Michel electron deposited energy resolution (Q+L);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyRes_Shwr3D    = truthDir.make<TH1D>("True_EnergyRes_Shwr3d","Michel electron deposited energy resolution (3D shower locations found);(reco-true)/true",150,-1.5,1.5);
  hTrue_EnergyResTrk_Shwr3D  = truthDir.make<TH1D>("True_EnergyResTrk_Shwr3d","Michel electron deposited energy resolution (electron track only);(reco-true)/true",150,-1.5,1.5);
  hTrue_IDE_X                 = truthDir.make<TH1D>("True_IDE_X","X position of energy depositions [cm]",600,-5.,55.);
  hTrue_IDE_Y                 = truthDir.make<TH1D>("True_IDE_Y","Y position of energy depositions [cm]",500,-25.,25.);
  hTrue_IDE_Z                 = truthDir.make<TH1D>("True_IDE_Z","Z position of energy depositions [cm]",1000,-5.,95.);
  hTrue_IDE_NumElectrons = truthDir.make<TH1D>("True_IDE_NumElectrons","Number electrons in IDEs (drift attenuation applied)",2500,0.,2500);
  hTrue_dQdx              = truthDir.make<TH1D>("True_dQdx","dQ/dx (avg. per particle) for Michel trk + shower depositions;dQ/dx [e-/cm]",300,0.,15e4);
  hTrue_dQdx_CrsMu           = crsMuDir.make<TH1D>("True_dQdx_CrsMu","Ave dQ/dx for muon track;dQ/dx [e-/cm]",300,0.,15e4);
  hTrue_dQdx_ElTrk        = truthDir.make<TH1D>("True_dQdx_ElTrk","dQ/dx (avg. per particle) for Michel trk depositions;dQ/dx [e-/cm]",300,0.,15e4);
  hTrue_dQdx_ElShw        = truthDir.make<TH1D>("True_dQdx_ElShw","dQ/dx (avg. per particle) for Michel shower depositions;dQ/dx [e-/cm]",300,0.,15e4);
  hTrue_dEdx              = truthDir.make<TH1D>("True_dEdx","dE/dx (avg. per particle) for Michel trk + shower depositions;dE/dx [MeV/cm]",240,0.,12.);
  hTrue_dEdx_CrsMu           = crsMuDir.make<TH1D>("True_dEdx_CrsMu","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_CrsMu_4mm           = crsMuDir.make<TH1D>("True_dEdx_CrsMu_4mm","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_CrsMu_1cm           = crsMuDir.make<TH1D>("True_dEdx_CrsMu_1cm","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_CrsMu_2cm           = crsMuDir.make<TH1D>("True_dEdx_CrsMu_2cm","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_CrsMu_5cm           = crsMuDir.make<TH1D>("True_dEdx_CrsMu_5cm","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_CrsMu_10cm           = crsMuDir.make<TH1D>("True_dEdx_CrsMu_10cm","Ave dE/dx for muon track;dE/dx [MeV/cm]",500,0.,5.);
  hTrue_dEdx_ElTrk        = truthDir.make<TH1D>("True_dEdx_ElTrk","dE/dx (avg. per particle) for Michel trk depositions;dE/dx [MeV/cm]",240,0.,12.);
  hTrue_dEdx_ElShw        = truthDir.make<TH1D>("True_dEdx_ElShw","dE/dx (avg. per particle) for Michel shower depositions;dE/dx [MeV/cm]",240,0.,12.);
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
  hTrue_RecombFactor      = truthDir.make<TH1D>("True_RecombFactor","Survival prob of ionization electrons for all Michel depositions;P",500,0.,1.);
    hTrue_RecombFactor    ->SetOption("hist");
  hTrue_RecombFactor_ElTrk      = truthDir.make<TH1D>("True_RecombFactor_ElTrk","Survival prob of ionization electrons for Michel tracks;P",500,0.,1.);
    hTrue_RecombFactor_ElTrk    ->SetOption("hist");
  hTrue_RecombFactor_ElShw      = truthDir.make<TH1D>("True_RecombFactor_ElShw","Survival prob of ionization electrons for Michel shower photon depositions;P",500,0.,1.);
    hTrue_RecombFactor_ElShw    ->SetOption("hist");
  //hTrue_NumElectrons    = truthDir.make<TH1D>("True_NumElectrons","Number of electrons from IDE objects",40,0,4000);
  hTrue_DriftStochasticity = truthDir.make<TH1D>("True_DriftStochasticity","Fractional smearing on hit charge due to drift stochasticity",100,0.0,0.05);
  hTrue_MuEnergy        = truthDir.make<TH1D>("True_MuEnergy","Initial muon kinetic energy;Muon energy [MeV]",1000.,0.,10e3);
  hTrue_MuAveMomentum        = truthDir.make<TH1D>("True_MuAveMomentum","Average momentum of muon in TPC;Momentum [MeV/c]",1000,0.,10e3);
  hTrue_ElEnergy        = truthDir.make<TH1D>("True_ElEnergy","Electron energy",30, 0., 60.);
  hTrue_ElEnergyFree    = truthDir.make<TH1D>("True_ElEnergyFree","True Michel Electron Energy (from #mu+)",30,0.,60.);
  hTrue_ElEnergyCapture = truthDir.make<TH1D>("True_ElEnergyCapture","True Michel Electron Energy (from #mu-)",30,0.,60.);
  hTrue_ElTrackEnergyDep     = truthDir.make<TH1D>("True_ElTrackEnergyDep","Michel track energy deposited;Energy [MeV]",30,0.,60.);
  hTrue_ElShowerEnergyDep     = truthDir.make<TH1D>("True_ElShowerEnergyDep","True Deposited Michel Electron Energy;Energy [MeV]",30,0.,60.);
  hTrue_ElTrackChargeDep     = truthDir.make<TH1D>("True_ElTrackChargeDep","True Deposited Michel Electron track Charge;Electrons [e-]",300,0,3000000);
  hTrue_ElShowerChargeDep     = truthDir.make<TH1D>("True_ElShowerChargeDep","True Deposited Michel Electron Charge;Electrons [e-]",300,0,3000000);
  hTrue_EnergyDepVsRecoEnergy = truthDir.make<TH2D>("True_EnergyDepVsRecoEnergy","True Michel energy deposited in LAr vs. reco energy;True deposited Michel electron energy [MeV];Reconstructed energy [MeV]",
    120,0,60,240,0,120);
  hTrue_EnergyDepVsRecoEnergy->SetOption("colz");
  hTrue_EnergyVsEnergyRes     = truthDir.make<TH2D>("True_EnergyVsEnergyRes","Michel electron deposited energy resolution;True Deposited Energy [MeV];(reco-true)/true",60,0.,60.,150,-1.5,1.5);
  hTrue_EnergyVsEnergyRes     ->SetOption("colz");
  hTrue_ElShowerPhotons = truthDir.make<TH1D>("True_ElShowerPhotons","Number of scintillation photons produced",300,0,3000000);
  hTrue_MuLateLightContamination = truthDir.make<TH2D>("True_MuLateLightContamination","ETL PMT;#DeltaT [ns];Contamination Fraction", 75,0.,7500, 110, 0.,1.1);
  hTrue_MuLateLightContamination->SetOption("colz");
  hTrue_MuTrackEnd_ZX   = truthDir.make<TH2D>("True_MuTrackEnd_ZX","True #mu endpoint;Z [cm];X [cm]",100,-5.,95, 57,-5.,52.);
  hTrue_MuTrackEnd_ZY   = truthDir.make<TH2D>("True_MuTrackEnd_ZY","True #mu endpoint;Z [cm];Y [cm]",100,-5.,95, 50,-25.,25.);
  hTrue_MuTrackEnd_ZX   ->SetOption("colz");
  hTrue_MuTrackEnd_ZY   ->SetOption("colz");
  hTrue_MuTrkVtxRes_X = truthDir.make<TH1D>("True_MuTrkVtxRes_X","Tagged #mu track startpoint resolution (X);Xreco - Xtrue [cm]",200,-2.,2);
  hTrue_MuTrkVtxRes_Y = truthDir.make<TH1D>("True_MuTrkVtxRes_Y","Tagged #mu track startpoint resolution (Y);Yreco - Ytrue [cm]",200,-2.,2);
  hTrue_MuTrkVtxRes_Z = truthDir.make<TH1D>("True_MuTrkVtxRes_Z","Tagged #mu track startpoint resolution (Z);Zreco - Ztrue [cm]",200,-2.,2);
  hTrue_MuClsEndRes        = truthDir.make<TH1D>("True_MuClsEndRes",";dL [cm]",100,0.,5.);
  hTrue_MuClsEndRes_X        = truthDir.make<TH1D>("True_MuClsEndRes_X",";dL [cm]",100,-2.5,2.5);
  hTrue_MuClsEndRes_Y        = truthDir.make<TH1D>("True_MuClsEndRes_Y",";dL [cm]",100,-2.5,2.5);
  hTrue_MuClsEndRes_Z        = truthDir.make<TH1D>("True_MuClsEndRes_Z",";dL [cm]",100,-2.5,2.5);
  hTrue_MuClsEndResDir      = truthDir.make<TH1D>("True_MuClsEndResDir","Resolution along direction of #mu;dL [cm]",100,-5.,5.);
  hTrue_MuTrkEndRes        = truthDir.make<TH1D>("True_MuTrkEndRes",";dL [cm]",100,0.,5.);
  hTrue_MuTrkEndRes_X         = truthDir.make<TH1D>("True_MuTrkEndRes_X",";dL [cm]",100,-2.5,2.5);
  hTrue_MuTrkEndRes_Y         = truthDir.make<TH1D>("True_MuTrkEndRes_Y",";dL [cm]",100,-2.5,2.5);
  hTrue_MuTrkEndRes_Z         = truthDir.make<TH1D>("True_MuTrkEndRes_Z",";dL [cm]",100,-2.5,2.5);
  hTrue_MuTrkEndResDir        = truthDir.make<TH1D>("True_MuTrkEndResDir","#mu endpt separation (reco-true) projected along #mu reco dir;#Delta L_{dir} [cm]",160,-4.,4.);
  hTrue_MuTrkEndResTrans      = truthDir.make<TH1D>("True_MuTrkEndResTrans","#mu enpt separation (reco-true) transverse to #mu reco dir;#Delta L_{T} [cm]",80,0.,4.);
  hTrue_MuTrkEndResDir_cut    = truthDir.make<TH1D>("True_MuTrkEndResDir_cut","#mu endpt separation (reco-true) projected along #mu reco dir (drift time match cut);#Delta L_{dir} [cm]",160,-4.,4.);
  hTrue_MuTrkEndResTrans_cut  = truthDir.make<TH1D>("True_MuTrkEndResTrans_cut","#mu enpt separation (reco-true) transverse to #mu reco dir (drift time match cut);#Delta L_{T} [cm]",80,0.,4.);
  hTrue_ElShowerDepVsDistFromMuTrackEnd = truthDir.make<TH1D>("True_ElShowerDepVsDistFromMuTrackEnd","True distance of Michel shower depositions;Distance [cm]",100, 0., 100.);
  hTrue_PartDispFromElectron  = truthDir.make<TH1D>("True_PartDispFromElectron",";cm",500,0.,50.);
  hTrue_BremmPhotonLength = truthDir.make<TH1D>("True_BremmPhotonLength","Bremm photon length;cm",100,0.,100);
  hTrue_ElEnergyVsNumBremmPhotons = truthDir.make<TH2D>("True_ElEnergyVsNumBremmPhotons",";True Michel electron energy [MeV];Number of Bremmstralung photons",60,0.,60., 30,0.,30.);
  hTrue_ElEnergyVsNumBremmPhotons ->SetOption("colz");

}





//########################################################################################
// Resets all member data variables/counters before the start of a new event
void MichelAna::ResetVariables()
{
  fMichelOpticalID          = false;
  fSingleOpticalHit         = false;
  fPulses.clear();
  fRunNumber                = -99;
  fSubRunNumber             = -99;
  fEventNumber              = -99;
  fEventTime                = -99;
  fDecayTime                = -999.;
  fPE_comb                  = -999.;
  fPE_comb_qc               = -999.;
  fPE_prompt_comb           = -999.;
  fTriggerTime              = -99.;
  for(size_t i=0; i<2; i++) {
    fWfmRMS[i]              = -99.; 
    fNumOpHits0[i]          = -9; 
    fNumOpHits[i]           = -9; 
    fdT[i]                  = -999.;
    fWidth[i]               = -999.;
    fMuWidth[i]             = -999.;
    fAmplitude[i]           = -999.;
    fPE_total[i]            = -999.;
    fPE_total_qc[i]            = -999.;
    fPE_totalRaw[i]         = -999.;
    fPE_prompt[i]           = -999.;
    fPE_promptRaw[i]        = -999.;
    fPromptFracRaw[i]       = -999.;
    fPromptFrac[i]          = -999.;
    fMuAmplitude[i]         = -999.;
    fMuContam_prompt[i]     = -999.;
    fMuContam_total[i]     = -999.;
    fMuPE_prompt[i]         = -999.;
    fMuPE_total[i]          = -999.;
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
    fvHitADC_total[i]      .clear();
    fvHitADCpf_100ns[i]     .clear();
    fvHitADCpf_total[i]    .clear();
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
    fTrue_MuPE_total[i]       = -999.;
    fTrue_MuContam_prompt[i] = -999.;
    fTrue_MuContam_total[i] = -999.;
    fTrue_PE_total[i]       = -999.;
    hPMT_phelTimes_electron[i]->Reset();
    hPMT_phelTimes_muon[i]->Reset();
    hPMT_phelTimes[i]       ->Reset();
    fPMT_wfm[i]             .clear();
    fPMT_wfm_raw[i]         .clear();
    fvbs[i]                 .clear();
    fT0[i]                  = 0;
    fMuT0[i]                = 0;
  }
  fNumTrack                = -9;
  fNumTrackStopping         = 0;
  fNumTrackCrossing          = 0;
  fNumTrackContained        = 0;
  
  //fvMuClusterdEdx            .clear();
  //fvMuClusterdQdx            .clear();
  //fvMuClusterdADCdx            .clear();
  //fvMuClusterdEdx_Pl0            .clear();
  //fvMuClusterdQdx_Pl0            .clear();
  //fvMuClusterdADCdx_Pl0            .clear();
  //fvMuClusterResRange   .clear();
  //fvMuClusterdEdx_Pl0            .clear();
  //fvMuClusterdQdx_Pl0            .clear();
  //fvMuClusterdADCdx_Pl0            .clear();
  //fvMuClusterResRange_Pl0   .clear();
  
  fMuTrackIsCaloOrdered     = false;
  fMuTrackHits              .clear();
//  fBgTrackHits              .clear();
  fMuTrackID                = -999 ;
  fMuTrackIndex             = -999 ;
  fMuTrackLength            = -999.;
  fMuTrackEnergy            = -999.;
  fMuTrackZenithAngle            = -999.;
  fMuTrackInclineAngle            = -999.;
  fMuTrackVertex            .SetXYZ(-99.,-99.,-99.);
  fMuTrackEnd               .SetXYZ(-99.,-99.,-99.);
  fMuTrackEndDir            .SetXYZ(0., 0., 0.);
  fMuTrackVertex_X           = -99.;
  fMuTrackVertex_Y           = -99.;
  fMuTrackVertex_Z           = -99.;
  fMuTrackEnd_X           = -99.;
  fMuTrackEnd_Y           = -99.;
  fMuTrackEnd_Z           = -99.;
  //fMuTrackPitch[0]        = -9.;
  //fMuTrackPitch[1]        = -9.;

  // Mu calibration variables
  fCrsMuVertex         .SetXYZ(-99.,-99.,-99.);
  fCrsMuEnd         .SetXYZ(-99.,-99.,-99.);
  fCrsMuVertex_X       = -999.;
  fCrsMuVertex_Y       = -999.;
  fCrsMuVertex_Z       = -999.;
  fCrsMuEnd_X       = -999.;
  fCrsMuEnd_Y       = -999.;
  fCrsMuEnd_Z       = -999.;
  fCrsMuHitFrac      = -9.;
  fCrsMuPitch[0]       = -9.;
  fCrsMuPitch[1]       = -9.;
  fCrsMuTrackIndex     = -999 ;
  fCrsMuLength        = -9.;
  fCrsMuEnergy        = -9.;
  fCrsMuCharge        = -9999.;
  fCrsMuCharge_Pl0        = -9999.;
  fCrsMuIntegral        = -9999.;
  fCrsMuIntegral_Pl0        = -9999.;
  fCrsMuPE_total[0]           = -999.; 
  fCrsMuPE_total[1]           = -999.; 
  fCrsMuPE_prompt[0]           = -999.; 
  fCrsMuPE_prompt[1]           = -999.; 
  fCrsMuPhotons         = -9999.;
  fCrsMuPhotonsPrompt         = -9999.;
  fTrue_CrsMuCharge    = -9.;
  fTrue_CrsMuPhotons   = -9.;

  // Shower/charge
  fHitPlane               .clear();
  fHitWire                .clear();
  fHitT                   .clear();
  fHitX                   .clear();
  fHitW                   .clear();
  fHitCharge              .clear();
  fHitIntegral            .clear();
  fHitSummedADC           .clear();
  fHitRMS                 .clear();
  fHitAmplitude           .clear();
  
  fFracDropoffAtBnd     = -9;
  fFracDropoffAtBnd_Pl0 = -9;
  fTruncMeanSlope     = -9;
  fTruncMeanSlope_Pl0 = -9;
  
  fClusterSize          = -9;
  fMuEnd2D_X            = -9.;
  fMuEnd2D_W            = -9.;
  fElDir2D_X            = 0.;
  fElDir2D_W            = 0.;
  fMuClusterSize        = -9;
  fElClusterSize        = -9;
  fExtraHits            = -9;
  fMuAveLinearity       = -9.;
  fFracMuHitsLinear     = -9.;
  fMuClusterHitsEndFit  = -9;
  fDecayAngle2D         = -999;
  fElShowerSize         = -9;
  fElShowerFrac         = -9.;
  for(size_t i=0; i<2; i++){
    fNumPlaneHits[i]          = 0;
    fTotalIntegral[i]     = -9.;
    fTotalCharge[i]       = -9.;
    fvMuTrkXYZ[i]              .clear();
    fvMuTrkdEdx[i]             .clear();
    fvMuTrkdQdx[i]             .clear();
    fvMuTrkdADCdx[i]           .clear();
    fvMuTrkResRange[i]         .clear();
    fvMuTrkPitch[i]         .clear();
  }
  fMuCharge             = -9.;
  fElShowerCharge       = -9.;
  fElTrackCharge        = -9.;
  fElShowerEnergy       = -9.;
  fElTrackEnergy        = -9.;                 
  
  fClusterSize_Pl0          = -9;
  fMuEnd2D_X_Pl0            = -9.;
  fMuEnd2D_W_Pl0            = -9.;
  fElDir2D_X_Pl0            = 0.;
  fElDir2D_W_Pl0            = 0.;
  fMuClusterSize_Pl0        = -9;
  fElClusterSize_Pl0        = -9;
  fExtraHits_Pl0            = -9;
  fMuAveLinearity_Pl0       = -9.;
  fFracMuHitsLinear_Pl0     = -9.;
  fMuClusterHitsEndFit_Pl0  = -9;
  fDecayAngle2D_Pl0         = -999;
  fElShowerSize_Pl0         = -9;
  fElShowerFrac_Pl0         = -9.; 
  fElShowerCharge_Pl0       = -9.;
  fElTrackCharge_Pl0        = -9.;
  fElShowerEnergy_Pl0       = -9.;
  fElTrackEnergy_Pl0        = -9.;                 
  
  fProjDist2DToWires        = -999.;
  fProjDist3DToWires        = -999.;
  fProjPtOnWires_Y          = -999.;
  fProjPtOnWires_Z          = -999.;
  
  fMuEnd3D                  .SetXYZ(-999., -999., -999.);
  fMuEnd3D_X                = fMuEnd3D.X();
  fMuEnd3D_Y                = fMuEnd3D.Y();
  fMuEnd3D_Z                = fMuEnd3D.Z();
  fElDir3D                  .SetXYZ(0.,0.,0.);
  fElDir3D_X                = fElDir3D.X();
  fElDir3D_Y                = fElDir3D.Y();
  fElDir3D_Z                = fElDir3D.Z();
  fMuEndHitTimeDiff         = -999.;
  fElShowerPhotons          = -999.;
  fElShowerPhotonsPrompt          = -999.;
  //fElTrackdEdx              = -999.;
  fFracHits3D               = -9.;
  fNumPts3D                 = -9;
  fElShowerVis              = -0.09;
  fElShowerVisCh[0]         = -0.09;
  fElShowerVisCh[1]         = -0.09;
  fElShowerCentroid         .SetXYZ(-99,-99,-99.);
  fElShowerCentroid_X       = -99.;
  fElShowerCentroid_Y       = -99.;
  fElShowerCentroid_Z       = -99.;
  fElShowerEnergyQL         = -999.;
  
  // Truth information
  fTrue_PE0                  = -999.;
  fTrue_PE                  = -999.;
  fTrue_PEPrompt                  = -999.;
  fTrue_MuPE              = -999.;
  fTrue_MuSign          = 0;
  fTrue_MuEnergy          = -999.;
  fTrue_IsMuStopping          = false;
  fTrue_MuStopTime        = -999.;
  fTrue_MuTrackLength          = -9.;
  fTrue_MuTrackEnd          .SetXYZ(-99.,-99.,-99.);
  fTrue_MuTrackEnd_X      = -99.;
  fTrue_MuTrackEnd_Y      = -99.;
  fTrue_MuTrackEnd_Z      = -99.;
  fTrue_MuTrackVertex_X      = -99.;
  fTrue_MuTrackVertex_Y      = -99.;
  fTrue_MuTrackVertex_Z      = -99.;
  fTrue_MuTrackVertex          .SetXYZ(-99.,-99.,-99.);
  fTrue_MuEnergyDep       = -999.;
  fTrue_MuCharge       = -999.;
  fTrue_MuPhotons       = -999.;
  fTrue_ProjDist3DToWires = -999.;
  fTrue_ElMomentum        .SetXYZ(0.,0.,0.);
  fTrue_ElMomentum_X      = -9.;
  fTrue_ElMomentum_Y      = -9.;
  fTrue_ElMomentum_Z      = -9.;
  fTrue_ElAngle           = -999.;
  fTrue_ElEnergy          = -999.;
  fTrue_ElShowerVis       = -0.09;
  fTrue_ElShowerVisCh[0]       = -0.09;
  fTrue_ElShowerVisCh[1]       = -0.09;
  fTrue_ElShowerCentroid_X = -999.;
  fTrue_ElShowerCentroid_Y = -999.;
  fTrue_ElShowerCentroid_Z = -999.;
  fTrue_ElShowerPhotons   = -999;
  fTrue_ElTrackPhotons   = -999;
  fTrue_ElShowerPhotonsPrompt   = -999;
  fTrue_ElShowerPhotonsLate   = -999;
  fTrue_ElTrackEnergyDep  = -999.;
  fTrue_ElTrackChargeDep  = -999.;
  fTrue_ElShowerChargeDep = -999.;
  fTrue_ElShowerTotalTrajLength = -999.;
  fTrue_ElShowerIons = -999.;
  fTrue_ElShowerEnergyDep = -999.;
  fTrue_IsElTrkContained     = false;
  fTrue_IsElShwrContained     = false;
  fTrue_ContainmentFrac = -9.;
  fTrue_TotalEnergyDep    = -999.;
  fTrue_TotalChargeDep    = -999.;
  fTrue_dT                = -999;
}



//########################################################################################
//
// The main MichelAna algorithm. This function takes in an event that has already gone
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
bool MichelAna::filter(art::Event & e)
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
  
  // ================================================   
  // Get detector properties
  auto const* detprop     = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fEfield                 = detprop->Efield(0);   
  fElectronLifetimeFromDB = detprop->ElectronLifetime();  // microseconds
  fElectronLifetime       = fElectronLifetimeFromDB;
  fDriftVelocity[0]       = detprop->DriftVelocity(fEfield,detprop->Temperature());
  fDriftVelocity[1]       = detprop->DriftVelocity(detprop->Efield(1),detprop->Temperature());
  fDriftVelocity[2]       = detprop->DriftVelocity(detprop->Efield(2),detprop->Temperature());
  fTriggerOffset          = detprop->TriggerOffset();
  fXTicksOffset[0]        = detprop->GetXTicksOffset(0,0,0);
  fXTicksOffset[1]        = detprop->GetXTicksOffset(1,0,0);
  fSamplingRate           = detprop->SamplingRate()*1e-3;
 
  /* 
  std::cout<<"Trigger offset "<<detprop->TriggerOffset()<<"\n"; 
  std::cout<<"XTicksOffset (ind) "<<fXTicksOffset[0]<<"\n";
  std::cout<<"XTicksOffset (col) "<<fXTicksOffset[1]<<"\n";
  */

  // If this is MC, smear the "measured" electron lifetime
  // to simulate uncertainty estimated in data
  if( fIsMC && fElectronLifetimeErr > 0 ) 
    fElectronLifetime += fRand->Gaus(0., fElectronLifetimeFromDB*fElectronLifetimeErr);
  
  // Get SPE info
  GetDetProperties(); 
 
  LOG_VERBATIM("MichelAna")
  <<"\n"
  <<"Beginning filter: run "<<fRunNumber<<", subrun "<<fSubRunNumber<<", event "<<fEventNumber<<" (MC = "<<fIsMC<<")\n"
  <<"total evts processed: "<<fNumEvents
  <<", 1 stp trk: "<<fNumEvents_StpTrk
  <<", OpID: "<<fNumOpticalID
  <<", 2D showers: "<<fNumEvents_Shwr
  <<", 3D showers: "<<fNumEvents_Shwr3D
  <<", calib: "<<fNumEvents_Calibration; 
  LOG_VERBATIM("MichelAna")
  <<"Efield: "<<fEfield<<"  Drift vel.: "<<fDriftVelocity[0]<<"  Electron lifetime: "<<fElectronLifetimeFromDB<<"  XTicks offset[1]: "<<fXTicksOffset[1]<<"\n"
  <<"  SPE: "<<fSPE[0]<<","<<fSPE[1]<<"; eff. tau: "<<fMuContamCorr_EffTau[0]<<","<<fMuContamCorr_EffTau[1]<<"  (true TauT: "<<fTauT<<")";
  
  // Skip event if the electron lifetime isn't defined for this run
  hElectronLifetime       ->Fill(fElectronLifetime); 
  hElectronLifetimeFromDB ->Fill(fElectronLifetimeFromDB); 
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
    
    // Associations between tracks and hits
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(TrackHandle, e, fTrackModule);

    fNumTrack = (int)TrackHandle->size();
    LOG_VERBATIM("MichelAna") << "Using track module "<<fTrackModule<<", number of tracks: "<<fNumTrack;
  
    // initialize the index of the identified stopping or crossing track 
    // as well as other variables (to prevent repeated loops thru tracks)
    int   crsTrkIndex     = -9;
    float crsTrkLength    = -9.;
    float crsTrkIncline   = -9.;
    int   stpTrkIndex     = -9;
    int   stpTrkID        = -9;
    float stpTrkLength    = -9.;
    float stpTrkZenith    = -9.;
    float stpTrkIncline   = -9.;
    TVector3 stpTrkEndDir;
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
      // Where are the track nodes rel. to fiducial volume
      bool isCrossingTrack = false;
      bool isStoppingTrack = false;
      bool isContainedTrack = false;
      bool endIsAtEdge     = !IsPointInFiducialVolume(TheTrack.End(),fEdgeMarginX,fEdgeMarginY,fEdgeMarginZ);
      bool vertexIsAtEdge  = !IsPointInFiducialVolume(TheTrack.Vertex(),fEdgeMarginX,fEdgeMarginY,fEdgeMarginZ);

      // ------------------------------------------------------------
      // First categorize the trivial cases: crossing and contained trks
      if      ( endIsAtEdge && vertexIsAtEdge   ) { isCrossingTrack = true;}
      else if ( !endIsAtEdge && !vertexIsAtEdge  ) { isContainedTrack = true;}
      
      // ------------------------------------------------------------
      // Call the node that's outside the fiducial volume the "vertex".
      // If both are outside fiducial volume (ie, crossing track), then
      // whichever has higher 'Y' coordinate is the vertex.
      TVector3 vertex = TheTrack.Vertex();
      TVector3 end    = TheTrack.End();
      bool trkIsReversed = false;
      if( !isContainedTrack ) {
        // swap vertex/end if either of these is met:
        if(     (isCrossingTrack  && vertex.Y() < end.Y())
            ||  (!isCrossingTrack && endIsAtEdge)  ){
          vertex  = TheTrack.End();
          end     = TheTrack.Vertex();
          trkIsReversed = true;
          bool ee = endIsAtEdge;
          bool ve = vertexIsAtEdge;
          endIsAtEdge = ve;
          vertexIsAtEdge = ee;
        }
      }
      

      // --------------------------------------------------------------
      // Calculate zenith angle and angle of inclination rel. to E-field
      float zenithAngle = -9.;
      float inclineAngle = -9.;
      TVector3 vert(0.,-1.,0.);
      TVector3 dir=(end-vertex);
      if( dir.Mag() ) {
        dir.SetMag(1);
        zenithAngle = vert.Angle(dir);
        inclineAngle  = asin( sqrt( dir.Y()*dir.Y() + dir.Z()*dir.Z() ));
        hTrackZenithAngle->Fill( zenithAngle*RAD_TO_DEG ); 
        hTrackLength_vs_ZenithAngle->Fill( TheTrack.Length(), zenithAngle*RAD_TO_DEG );
      }
      
      // -------------------------------------------------------------
      // Now to determine if this is a candidate "stopping mu track". 
      // We've already defined the vertex as the track node that is 
      // at the edge of the TPC.
      if(   !isCrossingTrack && !isContainedTrack //&& endIsInFid  
        &&  vertex.Y() > -1.*fGeo->DetHalfHeight()+fEdgeMarginY
        &&  TheTrack.Length() > fMuTrackLengthMin
        //&&  zenithAngle*RAD_TO_DEG < fMuTrackZenithAngleMax 
        ){
        isStoppingTrack = true;
      } 

      fNumTrackStopping   += isStoppingTrack; 
      fNumTrackCrossing   += isCrossingTrack;
      fNumTrackContained  += isContainedTrack;
      
      LOG_VERBATIM("MicheAnaFilter")
          <<"  track "<<track_index<<", N = "<<TheTrack.NumberTrajectoryPoints()<<"  length "<<TheTrack.Length()<<"  vertex("
          << vertex.X() <<"," 
          << vertex.Y() << "," 
          << vertex.Z() << ")->InFiducial()="
          << !vertexIsAtEdge
          << "   end("
          << end.X() <<"," 
          << end.Y() << "," 
          << end.Z() << ")->InFiducial()="
          << !endIsAtEdge;
      
      if( isStoppingTrack ) {
        stpTrkIndex = track_index;
        stpTrkID    = TheTrack.ID();
        stpTrkLength= TheTrack.Length();
        stpTrkVertex= vertex;
        stpTrkEnd   = end;
        stpTrkZenith = zenithAngle;
        stpTrkIncline = inclineAngle;
        stpTrkEndDir  = TheTrack.EndDirection();
        if( trkIsReversed ) stpTrkEndDir = -1.*TheTrack.VertexDirection();
      }
      if( isCrossingTrack ) {
        crsTrkIndex   = track_index;
        crsTrkLength  = TheTrack.Length();
        crsTrkVertex  = vertex;
        crsTrkEnd     = end;
        crsTrkIncline = inclineAngle;
      }
     
      // ------------------------------------------------------------- 
      // Create a collection of hits associated with any "background"
      // tracks, which include
      //  a) crossing tracks
      //  b) other tracks with L > some threshold that have small zenith angle (?)
      //if( isCrossingTrack && fmthm.isValid() ){
      //  auto vhit = fmthm.at(TheTrack.ID());
      //  for (size_t h = 0; h < vhit.size(); ++h) fBgTrackHits.push_back(vhit[h].key());
      //}
      
    
    }//<-- end first pass loop through tracks
    LOG_VERBATIM("MichelAna") << "Done looping tracks: num stop/pass/contained "<<fNumTrackStopping<<" / "<<fNumTrackCrossing<<" / "<<fNumTrackContained;
    hNumTrack         ->Fill(fNumTrack);   
    hNumTrackStopping ->Fill(fNumTrackStopping);
    hNumTrackContained ->Fill(fNumTrackContained);
    hNumTrackCrossing ->Fill(fNumTrackCrossing);
  

    // ===================================================================
    // Identify events with 1 crossing track and no others for calibration
    if( fNumTrackCrossing == 1 && fNumTrack == 1 ) {
        isOneCrossingTrack = true;
        fCrsMuTrackIndex = crsTrkIndex;
        fCrsMuLength     = crsTrkLength;
        fCrsMuVertex_X   = crsTrkVertex.X();
        fCrsMuVertex_Y   = crsTrkVertex.Y();
        fCrsMuVertex_Z   = crsTrkVertex.Z();
        fCrsMuEnd_X      = crsTrkEnd.X();
        fCrsMuEnd_Y      = crsTrkEnd.Y();
        fCrsMuEnd_Z      = crsTrkEnd.Z();
        fCrsMuVertex     = crsTrkVertex;
        fCrsMuEnd        = crsTrkEnd;
        fMuTrackInclineAngle= crsTrkIncline;
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
      fMuTrackZenithAngle = stpTrkZenith;
      fMuTrackInclineAngle= stpTrkIncline;
      fMuTrackEndDir      = stpTrkEndDir;
      
      hMuTrackZenithAngle->Fill( fMuTrackZenithAngle*RAD_TO_DEG);
      hMuTrackInclineAngle->Fill( fMuTrackInclineAngle*RAD_TO_DEG);

      LOG_VERBATIM("MichelAna")
      <<"Tagged muon track, index = "<<fMuTrackIndex<<", stopping point = ("<<fMuTrackEnd_X<<","<<fMuTrackEnd_Y<<","<<fMuTrackEnd_Z<<")";
    
      // Track pitch on both planes
      /*
      TVector3 dir = (fMuTrackEnd-fMuTrackVertex);
      dir.SetMag(1.);
      geo::View_t views[2]={ geo::kU, geo::kV };
      for(size_t i=0; i<2; i++){
        float wireVertAng = fGeo->WireAngleToVertical(views[i], 0, 0) - 0.5*TMath::Pi();
        float wirePitch = fGeo->WirePitch(views[i], 0, 0);
        float cosgamma  = fabs( sin(wireVertAng)*dir.Y() + cos(wireVertAng)*dir.Z() );
        if( cosgamma > 0 ) {
          fMuTrackPitch[i] = wirePitch / cosgamma;
        }
      }
      */
        
      // ------------------------------------------------------------
      // Fill truth plots if this is MC to compare to true end/vertex
      if( fTrue_MuTrackEnd_X > 0. ) {
        hTrue_MuTrkVtxRes_X->Fill( (fMuTrackVertex_X - fTrue_MuTrackVertex_X) );
        hTrue_MuTrkVtxRes_Y->Fill( (fMuTrackVertex_Y - fTrue_MuTrackVertex_Y));
        hTrue_MuTrkVtxRes_Z->Fill( (fMuTrackVertex_Z - fTrue_MuTrackVertex_Z));
      }
       
      // ------------------------------------------------------------ 
      // Look at the dEdx of this track
      // Association between Calorimetry objects and Tracks
      art::FindManyP<anab::Calorimetry> fmcal(TrackHandle, e, fTrackCalModule);
      
      if( fmcal.isValid() ){
        
        std::vector<art::Ptr<anab::Calorimetry> > calos_mu = fmcal.at(fMuTrackIndex);
        
        size_t N;
        size_t plane = 0;
       
        for(size_t i=0; i<calos_mu.size(); i++){
          plane = calos_mu[i]->PlaneID().Plane;
          N     = calos_mu[i]->dEdx().size();
          if( plane == 1 ) fMuTrackEnergy = calos_mu[i]->KineticEnergy();
          if( plane >= 0 && plane < 2 && N > 0 ) {
            fvMuTrkResRange[plane]  = calos_mu[i]->ResidualRange();
            fvMuTrkdEdx[plane]      = calos_mu[i]->dEdx();
            fvMuTrkXYZ[plane]       = calos_mu[i]->XYZ();
            fvMuTrkPitch[plane]     = calos_mu[i]->TrkPitchVec();
            for(size_t j=0; j<N; j++){
              float T = calos_mu[i]->XYZ().at(j).X() / fDriftVelocity[0];
              float eLifetimeCorr = exp( T / fElectronLifetime );            
              float dADCdx = eLifetimeCorr*calos_mu[i]->dQdx().at(j);
              float dQdx_e = fCaloAlg.ElectronsFromADCArea(dADCdx,plane);
              fvMuTrkdADCdx[plane] .push_back( dADCdx );
              fvMuTrkdQdx[plane]   .push_back( dQdx_e);
            }
          }
        }
        
        // The point corresponding to the identified muon endpoint (based on 
        // previous fiducial containment arguments) should have the smaller 
        // residual range.  Check if the residual range makes sense based
        // on what we expect. Check on collection plane (ordering should be 
        // the same on other planes).
        N = fvMuTrkXYZ[1].size();
        if( N > 0 ) {
          float distToEnd_firstPt = (fvMuTrkXYZ[1].at(0)    -fMuTrackEnd).Mag();
          float distToEnd_lastPt  = (fvMuTrkXYZ[1].at(N-1)  -fMuTrackEnd).Mag();
          float resRange_firstPt  = fvMuTrkResRange[1].at(0);
          float resRange_lastPt   = fvMuTrkResRange[1].at(N-1);
          // scenario 1 (normal order):   last pt is closer to end, residual range starts large and gets smaller
          // scenario 2 (inverted order): first pt is closer to end, but residual starts small and gets larger so it's fine
          bool scenario1 = ((distToEnd_firstPt > distToEnd_lastPt) && (resRange_firstPt > resRange_lastPt));
          bool scenario2 = ((distToEnd_firstPt < distToEnd_lastPt) && (resRange_firstPt < resRange_lastPt)); 
          if( scenario1 || scenario2 ) fMuTrackIsCaloOrdered = true;
        }

      }//end if fmcal is valid
        
      // -------------------------------------------------------------
      // Finally, create vector of hit keys associated with this track
      if( fmthm.isValid() ) {
        auto vhit = fmthm.at(fMuTrackID);
        for (size_t h = 0; h < vhit.size(); ++h) {
          int hit_index = vhit[h].key();
          fMuTrackHits.push_back(hit_index);
        }
      }

    }//<--End if single stopping track
  }//<--end looking at tracks
  
 
  
  // **************************************************************************
  // Optical reconstruction: 
  //  - only perform IF there was either (a) exactly 1 stopping track, or
  //    (b) there was exactly 1 passing track and no other tracks in data.
  // **************************************************************************
  if( !fLookAtTracks || !fReq1StpTrk || (fReq1StpTrk && fMuTrackIndex >= 0) || (isOneCrossingTrack) ) { 
    
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
      LOG_VERBATIM("MichelAna")<<"Analyzing optical channel "<<ch;
     
      if( fTriggerTime < 0 ) fTriggerTime = fPulses[ipmt].FirstSample();
      
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
//        fOpHitBuilderAlg.SetTau( TauConversion(fMuContamCorr_EffTau[ch]) );
        
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
        if( fCorrectOvershootMode != "" && ch == 0 ){ 
          std::vector<short> hitTimesTmp = fOpHitBuilderAlg.GetHits(fPMT_wfm[ch], (size_t)fPulses[ipmt].FirstSample() );
          fOpHitBuilderAlg.CorrectWfmOvershoot(fPMT_wfm[ch], hitTimesTmp, fCorrectOvershootMode );
        }
        
        // Masked baseline subtraction
        if( fMskBaselineSubtr ) {
//          std::cout<<"Doing masked baseline subtraction for ch "<<ch<<"\n";
          fvbs[ch].resize(fPMT_wfm[ch].size());
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
        LOG_VERBATIM("MichelAna") << "  found "<<hitTimes.size()<<" hits (thresh = "<<fGradHitThresh[ch]<<")";
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
//        fOpHitBuilderAlg.SetTau( TauConversion(fMuContamCorr_EffTau[ch]) );
        
        LOG_VERBATIM("MichelAna") << "Performing hit reconstruction on PMT "<<ch<<" (#hits = "<<fvHitTimes0[ch].size()<<")";
    
        // ---------------------------------------------------------------------
        // Loop over the hits, perform reconstruction over each
        std::vector<float> hit_info(fIntegrationWindows.size()+1); 
        std::vector<float> hit_info_raw(fIntegrationWindows.size()+1);
        for (size_t ihit=0; ihit< std::min(fvHitTimes0[ch].size(), kMaxOpHits); ihit++){
          
          short hitTime = fvHitTimes0[ch].at(ihit); 
          short prevHitTime = 0;
          if( fvHitTimes[ch].size() > 0 ) prevHitTime = fvHitTimes[ch].at(fvHitTimes[ch].size()-1);
         
          // --------------------------------------------------------------------
          // Get raw hit integrals; skip hits with raw prompt PE < 0 or where amplitude undefined
          hit_info_raw  = fOpHitBuilderAlg.GetHitInfo(fPMT_wfm[ch], hitTime, 0, fIntegrationWindows, false);
          if( hit_info_raw[1] < 0 ) continue;

          // ---------------------------------------------------------------------
          // Use prepulse fit to get hit's width, amplitude, and prefit integrals;
          hit_info      = fOpHitBuilderAlg.GetHitInfo(fPMT_wfm[ch], hitTime, prevHitTime, fIntegrationWindows, true);

          float hitAmplitude = hit_info[0]; // ADC
          float hitPromptPE = hit_info[1]/fSPE[ch]; // PE
          float pulseWidth =  fOpHitBuilderAlg.PulseWidth;
          
          LOG_VERBATIM("MichelAna") << "    "<<ihit<<"  T = "<<hitTime<<"   amp: "<<hitAmplitude<<"   width = "<<fOpHitBuilderAlg.PulseWidth;
          
          hWidth_AllHits[ch]          ->Fill( pulseWidth );
          hAmplitudeVsWidth_AllHits[ch]->Fill( hitAmplitude, pulseWidth );
          hPromptPEVsWidth_AllHits[ch]->Fill( hitPromptPE,  pulseWidth );
         
          if( hitAmplitude < 0. ) continue;

          //.......................................................................
          // Delete narrow hits
          if( fDeleteNarrowHits && pulseWidth < fMinOpHitWidth[ch] ) {
            LOG_VERBATIM("MichelAna")<<"    Skipping this hit...  promptPE (prepulse fit) = "<<hitPromptPE;
            for(short i=hitTime - (short)pulseWidth; i< hitTime+3.*(short)pulseWidth; i++){
              float val = fOpHitBuilderAlg.fit_SlowNorm * exp( -1.*((float)i - fOpHitBuilderAlg.fit_ZeroPoint) / fOpHitBuilderAlg.fit_SlowTau );
              //float val = 0.;
              //std::cout<<"      "<<i<<"  "<<fPMT_wfm[ch].at(i)<<" --> "<<val<<"\n";
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
              LOG_VERBATIM("MichelAna") << "    Looking for matching hits on the other PMT(s)...";
              size_t nMatches = 0;
              for(size_t j=0; j<fPulses.size(); j++){
                size_t jch=fPulses.at(j).OpChannel();
                if( jch == ch ) continue;
                // loop through the hits found on this PMT
                for(size_t jhit=0; jhit<fvHitTimes0[jch].size(); jhit++){
                  float dT = fvHitTimes0[ch].at(ihit) - fvHitTimes0[jch].at(jhit);
                  if( fabs( dT ) < fMaxDiff_dT ) {
                    nMatches++;
                    LOG_VERBATIM("MichelAna")<<"      found match on PMT "<<jch<<" (dt = "<<dT<<")";
                    break;
                  }
                }
              }
              if( nMatches == fPulses.size()-1 ) {
                LOG_VERBATIM("MichelAna")<<"    found matches on all other PMTs, hit passes!";
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
          fvHitADCpf_total[ch]  .push_back(hit_info[12]);
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
          fvHitADC_total[ch]    .push_back(hit_info_raw[12]);
           
        }//<-- done looping over PMT hits
    
        LOG_VERBATIM("MichelAna")<<"    total hits saved: "<<fvHitTimes[ch].size();

      }//<-- end loop over PMTs

    }//<-- end hit reconstruction (DATA ONLY)
    LOG_VERBATIM("MichelAna")<<"Done reconstructing optical hits.";


    // ======================================================================
    // Redefine number of hits
    for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
      size_t ch   = fPulses[ipmt].OpChannel(); 
      fNumOpHits[ch]      = fvHitTimes[ch].size();
    }

    
    // ======================================================================
    // Now that all hits have been reconstructed, try and identify a Michel-like
    // decay topology
    if( !fReq1StpTrk || !fLookAtTracks || (fReq1StpTrk && fMuTrackIndex >= 0 ) ) {
      for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
        size_t ch   = fPulses[ipmt].OpChannel(); 

        // ------------------------------------------------------
        // Calc dT if there are two hits
        if( fNumOpHits[ch] == 2 ) {
          fdT[ch] = float( fvHitTimes[ch].at(1)-fvHitTimes[ch].at(0) );
          hdT[ch] ->Fill( fdT[ch] );
          hdT_vs_Amp[ch]->Fill( fdT[ch], fvHitAmplitudes[ch].at(1) );
          LOG_VERBATIM("MichelAna") << "Exactly 2 hits in PMT "<<ch<<"  dT = " << fdT[ch] << " ns.";
        
          // -------------------------------------------------------
          // If the SPE calibration for this PMT is known, calculate
          // photoelectron (PE) quantities
          if( fSPE[ch] > 0 ) { 
            
            // Define "raw" PE quantities
            fPE_promptRaw[ch] = fvHitADC_100ns[ch].at(1) / fSPE[ch];
            fPE_totalRaw[ch]  = fvHitADC_total[ch].at(1) / fSPE[ch];

            // Initialize prompt/total light with raw quantities
            fPE_prompt[ch]    = fPE_promptRaw[ch];
            fPE_total[ch]     = fPE_totalRaw[ch];
              
            // Muon late-light correction:
            // If dT > dTmin, do a fit to the muon late-light and subtract 
            // off the added PEs to the Michel pulse. Otherwise just use 
            // the integrals that used the pre-pulse fits from before.
            if( fMuContamCorr ) MuContamCorr(ch); //<-- this function defines fPE_prompt, fPE_total
          
            if( fIsMC ) {
              hTrue_dT_vs_MuContamRes[ch]->Fill( fdT[ch], (fMuContam_total[ch]-fTrue_MuContam_total[ch])/fTrue_MuContam_total[ch] );
            }
           
            //LOG_VERBATIM("MichelAna")
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
            hAmplitudeVsWidth[ch]->Fill(fvHitAmplitudes[ch].at(1)*0.195,  fvHitWidth[ch].at(1));
            hPromptPEVsWidth[ch]->Fill(fPE_prompt[ch],                 fvHitWidth[ch].at(1));
          }
            
          // For L, we want to correct for expected quenching. We use
          // the effective late-light lifetime from the ETL. We don't bother
          // correcting the prompt light, since this quenching is negligible.
          if( fCorrectQuenching && fPE_total[ch] > 0 ) 
            fPE_total_qc[ch] = CorrectForQuenching(fPE_total[ch], fPE_prompt[ch], fMuContamCorr_EffTau[ch]); 
          else 
            fPE_total_qc[ch] = fPE_total[ch];

        } // endif nHits = 2

        // ---------------------------------------------------------------------- 
        // If dT was defined, and the second hit is the "triggered"
        // hit (in time), then we have muon and Michel candidates!
        // Also require BOTH hits to be above the width threshold.
        if(     fdT[ch] > 0 
            //&&  ((!fReq1StpTrk)||(isOneStoppingTrack))
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
          LOG_VERBATIM("MichelAna")
          <<"  --> Michel candidate wfm found!  PE(prompt)= "<<fPE_prompt[ch]<<", PE(total)= "<<fPE_total[ch];
         
          // Replicate trigger efficiency in MC by scuplting the prompt PE dist
          // and cutting the event out based on some probability.
          if( fIsMC && fApplyTrigEffCut && fTrigEff_P[ch] > 0 ) {
            TF1 probCut("probCut","1./(1.+(x/[0])^[1] )");
            probCut.SetParameter(0,fTrigEff_P[ch]);
            probCut.SetParameter(1,fTrigEff_K[ch]);
            //std::cout<<probCut.GetParameter(0)<<"   "<<probCut.GetParameter(1)<<"\n";
            float r = fRand->Rndm();
            float p = probCut.Eval(fPE_prompt[ch]);
            //std::cout<<"  prompt PE = "<<fPE_prompt[ch]<<"\n";
            //std::cout<<"  p         = "<<p<<"\n";
            //std::cout<<"  r         = "<<r<<"\n"; 
            if( r < p ) {
              LOG_VERBATIM("MichelAna")<<"MC event marked for removal via trigger eff. cut.";
              return false;
            }
          }
       
          // Define the prompt light fraction and fill some histograms 
          if( fPE_totalRaw[ch] > 0 ) fPromptFracRaw[ch] = fPE_promptRaw[ch] / fPE_totalRaw[ch];
          if( fPE_total_qc[ch] > 0 ) {
            fPromptFrac[ch]   = fPE_prompt[ch] / fPE_total_qc[ch];
            hTrue_MuLateLightContamination->Fill( fTrue_dT, fTrue_MuContam_total[ch] / fPE_totalRaw[ch] );
          }

          // Save average waveforms and PE arrival profiles
          if( fMakeAveWfms ) {
            
            if( fdT[ch] > 2000 && !fvIsHitSaturated[ch].at(0) ) {
              if( !fIsMC && fvHitTimes[ch].at(0) > 500 && fvHitTimes[ch].at(1)+7000 < (int)fPMT_wfm[ch].size() ) {
                AddToAverageWfm( fPMT_wfm[ch], fvHitTimes[ch].at(0) - 200, fvHitTimes[ch].at(0) + 2000, hAveWfm_mu[ch], aveWfmCounter_mu[ch], -1);
                AddToAverageWfm( fPMT_wfm[ch], fvHitTimes[ch].at(1) - 700, fvHitTimes[ch].at(1) + 7000, hAveWfm_el[ch], aveWfmCounter_el[ch], -1);
              } 
              if( fIsMC && fvHitTimes[ch].at(1)+7000 < (int)hPMT_phelTimes[ch]->GetNbinsX() ) {
                AddToAverageWfm( hPMT_phelTimes[ch], fvHitTimes[ch].at(0) - 200, fvHitTimes[ch].at(0) + 2000, fvHitTimes[ch].at(0), hAveWfm_mu[ch], aveWfmCounter_mu[ch], -1);
                AddToAverageWfm( hPMT_phelTimes[ch], fvHitTimes[ch].at(1) - 700, fvHitTimes[ch].at(1) + 7000, fvHitTimes[ch].at(1), hAveWfm_el[ch], aveWfmCounter_el[ch], -1);
              }
            } 
            else if( fIsMC && fvHitTimes[ch].at(0)+7000 < (int)hPMT_phelTimes[ch]->GetNbinsX() ) {
              AddToAverageWfm( hPMT_phelTimes_muon[ch], fvHitTimes[ch].at(0) - 700, fvHitTimes[ch].at(0) + 7000, fvHitTimes[ch].at(0), hAveWfm[ch], aveWfmCounter[ch], -1);
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
//          std::cout<<"ch "<<p.OpChannel()<<" has RMS "<<fWfmRMS[p.OpChannel()]<<"\n";
          if( fWfmRMS[p.OpChannel()] > fMaxWfmRMS[p.OpChannel()] ) goodWfms = false;
        }


        if( goodWfms ) {
          hOpEventCuts->Fill(3-1);

          // Check if all detectors saw exactly two hits, with consistent dT
          float maxdiff = -999.;
          float sum_dT = 0.;
          float sum_T = 0.;
          float sum_pe = 0.;
          float sum_pe_qc = 0.;
          float sum_pe_pr = 0.;
          size_t N = 0;
          for( auto & p : fPulses ) {
            if( fdT[p.OpChannel()] > 0. && fPE_prompt[p.OpChannel()] > 0 && fPE_total[p.OpChannel()] > 0 ) { 
              sum_T       += fvHitTimes[p.OpChannel()].at(1);
              sum_dT      += fdT[p.OpChannel()];
              N++;
              if( fUsePmtForCalo[p.OpChannel()] ) {
                sum_pe      += fPE_total[p.OpChannel()];
                sum_pe_qc   += fPE_total_qc[p.OpChannel()];
                sum_pe_pr   += fPE_prompt[p.OpChannel()];
              }
            }
            for( auto & p2 : fPulses ) {
              if( p.OpChannel() == p2.OpChannel()) continue; // avoid matching with same PMT
              float diff = fabs(fdT[p.OpChannel()] - fdT[p2.OpChannel()]);
              if( diff > maxdiff ) maxdiff = diff;
            }
          }
          if( sum_pe > 0 && maxdiff >= 0 ) hDiff_dT->Fill(maxdiff);
          

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
                      fPE_comb         = sum_pe;
                      fPE_comb_qc      = sum_pe_qc;
                      fPE_prompt_comb   = sum_pe_pr;
                      LOG_VERBATIM("MichelAna")<<"***OPTICAL MICHEL ID***";
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
      if( fMichelOpticalID ) {
        hPE         ->Fill(fPE_comb);
      }
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
          hTrue_PE_prompt[ch] ->Fill(fTrue_PE_prompt[ch]);
          hTrue_PE_total[ch]  ->Fill(fTrue_PE_total[ch]);
          hAmplitude[ch]      ->Fill(fAmplitude[ch]);
          hPromptPEVsAmplitude[ch]->Fill(fPE_prompt[ch],fAmplitude[ch]);
          hMuPromptPEVsAmplitude[ch]->Fill(fMuPE_prompt[ch],fMuAmplitude[ch]);
          hPromptFrac[ch]->Fill(fPromptFrac[ch]);
          hPromptFracRaw[ch]->Fill(fPromptFracRaw[ch]);
          hContamFrac[ch]->Fill( fMuContam_total[ch] / fPE_totalRaw[ch] );
          if( fDecayTime > fdTcut ) {
            if( fPE_prompt[ch] > fPromptPECut[ch] ) { 
              hPromptFrac_dTcut[ch]->Fill(fPromptFrac[ch]);
              hPromptFracRaw_dTcut[ch]->Fill(fPromptFracRaw[ch]);
            }
            hContamFrac_dTcut[ch]->Fill( fMuContam_total[ch] / fPE_totalRaw[ch] );
            hPE_total_dTcut[ch]->Fill(fPE_total[ch]);
            hPE_prompt_dTcut[ch]->Fill(fPE_prompt[ch]);
          }
        }
      }//<-- endloop PMTs
      
      // Decay time histogram
      if( promptPEcut &&  fDecayTime >= fGateDelay && fDecayTime <= fMaxDecayTime ) hDecayTime->Fill(fDecayTime);

    }//<-- endif 1stptrk
    LOG_VERBATIM("MichelAna")<<"Done looping over PMTs for Michel optical ID.";
    

    // =================================================================
    // Is there exactly 1 hit?
    // Check that both PMTs saw *one* unsaturated hit at nearly the same time
    //std::cout<<"Checking that exactly 1 hit found...\n";
    float maxdiff = -999.;
    size_t N = 0;
    for( auto & p : fPulses ) {
      LOG_VERBATIM("MichelAna")
      <<"Ch "<<p.OpChannel()<<"   Nhits "<<fNumOpHits[p.OpChannel()];
      if( fNumOpHits[p.OpChannel()] == 1 ) {
        bool isThereSomethingThere = ( fIsMC || fvHitADC_100ns[p.OpChannel()].at(0) > 100);
        LOG_VERBATIM("MichelAna")
          <<"Ch "<<p.OpChannel()<<", is there a hit? "<<isThereSomethingThere<<"   "
          <<"sat? "<<fvIsHitSaturated[p.OpChannel()].at(0);

        if( !fvIsHitSaturated[p.OpChannel()].at(0) && isThereSomethingThere && fvIsHitAtTrigger[p.OpChannel()].at(0) ) {
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
      if( (N==1) || (maxdiff > 0. && maxdiff <= fMaxDiff_dT) ) {
        fSingleOpticalHit = true;
      }
    }

    // ==================================================================
    // If there's exactly 1 passing track, and 1 optical hit on each PMT,
    // then use this crossing muon event for calibration purposes
    if( fUseCrossingMuons && isOneCrossingTrack && fSingleOpticalHit ) {
      
      LOG_VERBATIM("MichelAna")<<"Exactly 1 stopping track found -- can we use it for calibration?";
        
      isCalibrationEvent = true;
      fNumEvents_Calibration++;
      
      // Save ave wfms
      if( isCalibrationEvent ) {
        for( auto & p : fPulses ) {
          int ch = p.OpChannel();
          if( fMakeAveWfms && fvHitTimes[ch].at(0)+7000 < 18000 ){
            if( !fIsMC )  AddToAverageWfm( fPMT_wfm[ch],        fvHitTimes[ch].at(0) - 700, fvHitTimes[ch].at(0) + 7000,                        hAveWfm_crsmu[ch], aveWfmCounter_crsmu[ch], -1);
            else          AddToAverageWfm( hPMT_phelTimes[ch],  fvHitTimes[ch].at(0) - 700, fvHitTimes[ch].at(0) + 7000, fvHitTimes[ch].at(0),  hAveWfm_crsmu[ch], aveWfmCounter_crsmu[ch], -1);
          }
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
    if( fMichelOpticalID || fSingleOpticalHit ) {
      LOG_VERBATIM("MichelAna")<<"Event PASSES optical filter.";
      return true;
    } else {
      return false;
    }
  }

    
    
    
    // First fill some diagnostic histograms going over all hits
    // Get recob::Hit information (all hits, before filter)
    art::Handle< std::vector<recob::Hit> > hitListHandleAll;
    std::vector<art::Ptr<recob::Hit> > hitlistAll;
    if (e.getByLabel(fHitsModule,"allhits",hitListHandleAll))
    {art::fill_ptr_vector(hitlistAll, hitListHandleAll);}
    for(size_t i=0; i<hitlistAll.size(); i++){
      int hitPlane = hitlistAll[i]->WireID().Plane;
      hHitAmp_vs_HitRMS[hitPlane]->Fill( hitlistAll[i]->PeakAmplitude(), hitlistAll[i]->RMS() );
      hHitAmp_vs_HitIntegral[hitPlane]->Fill( hitlistAll[i]->PeakAmplitude(), hitlistAll[i]->Integral() );
      hHitIntegral_vs_HitRMS[hitPlane]->Fill( hitlistAll[i]->Integral(), hitlistAll[i]->RMS() );
      hHitAmp_vs_HitIntegral[hitPlane]->Fill( hitlistAll[i]->PeakAmplitude(), hitlistAll[i]->Integral() );
      hHitRMS[hitPlane]             ->Fill( hitlistAll[i]->RMS() );
      hHitAmp[hitPlane]       ->Fill( hitlistAll[i]->PeakAmplitude() );
    }
    

  // **************************************************************************
  // IF this event is interesting to us based on the reconstruction up to now, 
  // save the wire hit information and apply time/offset and calibration scaling.
  // **************************************************************************
  int nHits[2]={0};
  if( (fMichelOpticalID && fMuTrackIndex >= 0) || ( isCalibrationEvent ) ){

    for(size_t i=0; i<2; i++){ 
      fTotalIntegral[i]     = 0.; 
      fTotalCharge[i]       = 0.;
    }
  
    // Get recob::Hit information
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (e.getByLabel(fHitsModule,fHitsInstance,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}
    size_t  nWireHits = hitlist.size();

    // Assign a location to each hit, relative to projections on wire plane
    fHitRMS                 .reserve( nWireHits );
    fHitAmplitude                 .reserve( nWireHits );
    fHitPlane                 .reserve( nWireHits );
    fHitWire                  .reserve( nWireHits );
    fHitT                     .reserve( nWireHits );
    fHitX                     .reserve( nWireHits );
    fHitW                     .reserve( nWireHits );
    fHitIntegral              .reserve( nWireHits );
    fHitSummedADC             .reserve( nWireHits );
    fHitCharge                .reserve( nWireHits );
    for(size_t i=0; i<nWireHits; i++){
      int hitPlane = hitlist[i]->WireID().Plane;
      float hitTime0  = fSamplingRate*(hitlist[i]->PeakTime()-fTriggerOffset ); 
      float hitTime   = fSamplingRate*(hitlist[i]->PeakTime()-fXTicksOffset[hitPlane]);
      float ltCorFac = 1.;
      if( hitTime0 >= 0 ) { ltCorFac = exp(hitTime0/fElectronLifetime);}
      fHitPlane     .push_back(hitPlane);
      fHitT         .push_back(hitTime); // us
      fHitX         .push_back(hitTime*fDriftVelocity[0]);
      fHitWire      .push_back(hitlist[i]->Channel());
      if( hitlist[i]->Channel() < 240 ) fHitW.push_back(hitlist[i]->Channel()*0.4);
      else                              fHitW.push_back((hitlist[i]->Channel()-240)*0.4);
      fHitIntegral                  .push_back( hitlist[i]->Integral() );
      fHitSummedADC                 .push_back( hitlist[i]->SummedADC() );
      fHitCharge                    .push_back( fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(), hitPlane) * ltCorFac );
      fHitRMS                       .push_back( (hitlist[i]->RMS() ) );
      fHitAmplitude                 .push_back( (hitlist[i]->PeakAmplitude() ) );
      // Stochasticy in drift attenuation is NOT simulated. So,
      // add it back in for the MC.
      if( fIsMC ) {
        float prob  = 1./exp(hitTime0/fElectronLifetimeFromDB);
        float mean  = fHitCharge.at(i)*prob;
        float sigma = sqrt(mean*(1.-prob));
        float fac   = (fRand->Gaus(mean,sigma)*ltCorFac)/fHitCharge.at(i);
        fHitCharge.at(i)    *= fac;
        fHitIntegral.at(i)  *= fac;
        hTrue_DriftStochasticity->Fill(fabs(1-fac));
      }
      fTotalIntegral[hitPlane]      += (float)hitlist[i]->Integral()*ltCorFac;
      fTotalCharge[hitPlane]        += fHitCharge.at(i);
      fNumPlaneHits[hitPlane]       ++;
      hHitCharge[hitPlane]    ->Fill(fHitCharge.at(i));
      hHitT[hitPlane]               ->Fill( fHitT.at(i) );
      hHitX[hitPlane]               ->Fill( fHitX.at(i) );
      hHitIntegral[hitPlane]  ->Fill(hitlist[i]->Integral());
      nHits[hitPlane]++;
//      if( hitPlane == 1 ) {
//        if( std::find(fMuTrackHits.begin(), fMuTrackHits.end(), i ) != fMuTrackHits.end() ) {
//          float integral = hitlist[i]->Integral();
//          float sumadc    = hitlist[i]->SummedADC();
//          hHitIntegral_vs_DiffSummedADC_MuTrk->Fill( integral, (integral-sumadc)/sumadc );
//          hHitAmp_vs_DiffSummedADC_MuTrk     ->Fill( hitlist[i]->PeakAmplitude(), (integral-sumadc)/sumadc );
//          hHitIntegral_vs_SummedADC_MuTrk     ->Fill( integral, sumadc );
//        }
//      }
    }
    hTotalCharge->Fill(fTotalCharge[1]);
    hNumHits[0]->Fill(nHits[0]);
    hNumHits[1]->Fill(nHits[1]);
    

    LOG_VERBATIM("MichelAna")
    <<"Total charge hits in this event: "<<fHitX.size()<<" ("<<fNumPlaneHits[0]<<" ind, "<<fNumPlaneHits[1]<<" coll)";
    
  }

  

  // **************************************************************************
  // If this is a calibraiton event, get track length, overall track charge
  // deposition (after lifetime corrections), and visibility-corrected light
  // **************************************************************************
  if( isCalibrationEvent ) {
      
    
    // Get the 3D track object
    art::Handle< std::vector< recob::Track >> TrackHandle;
    e.getByLabel(fTrackModule,TrackHandle);
    art::Ptr< recob::Track > CrsTrackPtr(TrackHandle,fCrsMuTrackIndex);
    recob::Track CrsTrack = *CrsTrackPtr;
      
    // Get calo object
    art::FindManyP<anab::Calorimetry> fmcal(TrackHandle, e, fTrackCalModule);
    
    // Plot end/vertex onto histogram maps
    hCrsTrackNode_ZX->Fill( fCrsMuVertex_Z, fCrsMuVertex_X );
    hCrsTrackNode_ZX->Fill( fCrsMuEnd_Z, fCrsMuEnd_X );
    hCrsTrackNode_ZY->Fill( fCrsMuVertex_Z, fCrsMuVertex_Y );
    hCrsTrackNode_ZY->Fill( fCrsMuEnd_Z, fCrsMuEnd_Y );
   
    // If this is ACP track, record dT for drift time check
    bool isACP = false;
    if( (fCrsMuVertex_X < 2. && fCrsMuEnd_X > 45.5 ) || (fCrsMuEnd_X < 2. && fCrsMuVertex_X > 45.5 ) )
      isACP=true;
    
    // True charge/light (for MC, use total charge deposited
    // since it's the only particle simulated)
    fTrue_CrsMuCharge       = fTrue_TotalChargeDep;
    fTrue_CrsMuPhotons      = fTrue_MuPhotons;
    
    // Get average track pitch in each plane
    TVector3 dir = (fCrsMuEnd - fCrsMuVertex);
    dir.SetMag(1.);
    geo::View_t views[2]={ geo::kU, geo::kV };
    for(size_t i=0; i<2; i++){
      float wireVertAng = fGeo->WireAngleToVertical(views[i], 0, 0) - 0.5*TMath::Pi();
      float wirePitch = fGeo->WirePitch(views[i], 0, 0);
      float cosgamma  = fabs( sin(wireVertAng)*dir.Y() + cos(wireVertAng)*dir.Z() );
      if( cosgamma > 0 ) 
        fCrsMuPitch[i] = wirePitch / cosgamma;
    }
    
    // Calculate hit fraction by counting up hits on collection plane
    float muTrkQ = 0.;
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(TrackHandle, e, fTrackModule);
    int Nhits = 0;
    float minTime = 9999;
    float maxTime = -9999;
    if (fmthm.isValid()){
      auto vhit = fmthm.at(0);
      for (size_t h = 0; h < vhit.size(); ++h) {
        if( fHitPlane.at(vhit[h].key()) == 1 ) {
          float T = fHitT.at(vhit[h].key());
          if( T < minTime ) minTime = T;
          if( T > maxTime ) maxTime = T;
          muTrkQ += fHitCharge.at(vhit[h].key());
          Nhits++;
        }
      }
    }
    if( fNumPlaneHits[1] > 0 ) fCrsMuHitFrac = float(Nhits) / fNumPlaneHits[1]; 
    hCrsMuHitFrac->Fill( fCrsMuHitFrac );
    if( isACP ) hCrsTrackdT_ACP->Fill( maxTime - minTime );

    // From the calo object, get energy information. Also use
    // XYZ points to calculate the total vis.
    float mean_vis   = 0.;
    float muEnergy = 0.;
    if( fmcal.isValid() ) {
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(fCrsMuTrackIndex);
      for(size_t i=0; i<calos.size(); i++){
        int plane = calos[i]->PlaneID().Plane;
        int N     = calos[i]->dEdx().size();
        if( plane == 1 && N > 1 ) {
          //fCrsMuEnergy = calos[i]->KineticEnergy();
          for(int j=0; j<N; j++){
          
            // visibility for this point
            for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++)
              mean_vis += GetVisibility( calos[i]->XYZ().at(j), fPulses[ipmt].OpChannel() )/float(N); 
        
            // dE/dx for this point
            float T = calos[i]->XYZ().at(j).X() / fDriftVelocity[0];
            float eLifetimeCorr = exp( T / fElectronLifetime );            
            float dADCdx = eLifetimeCorr*calos[i]->dQdx().at(j);
            float dQdx_e = fCaloAlg.ElectronsFromADCArea(dADCdx,plane);
            float dEdx   = detprop->ModBoxCorrection( dQdx_e );
            float pitch  = calos[i]->TrkPitchVec().at(j);
            muEnergy += dEdx*pitch;
            hCrsMuPitch[plane]->Fill( pitch );
            if( pitch < 1.2 ) {
              hCrsMu_dQdx[plane]    ->Fill( dQdx_e ); 
              hCrsMu_dADCdx[plane]  ->Fill( dADCdx ); 
              hCrsMu_dEdx[plane]    ->Fill( detprop->ModBoxCorrection( dQdx_e ) );
            }

          }
        }
      }
    }
    fCrsMuEnergy = muEnergy;
        
    // Assign total reco'd charge
    fCrsMuIntegral         = fTotalIntegral[1];
    fCrsMuCharge           = fTotalCharge[1];
    fCrsMuIntegral_Pl0     = fTotalIntegral[0];
    fCrsMuCharge_Pl0       = fTotalCharge[0];
    
    // Assign total reco'd light and photons
    float sum_pe = 0;
    float sum_pe_pr = 0;
    for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
      size_t ch   = fPulses[ipmt].OpChannel();
      fCrsMuPE_total[ch]  = fvHitADC_total[ch].at(0) / fSPE[ch];
      fCrsMuPE_prompt[ch] = fvHitADC_100ns[ch].at(0) / fSPE[ch];
      if( fUsePmtForCalo[ch] ) {
        sum_pe    += fCrsMuPE_total[ch]; 
        sum_pe_pr += fCrsMuPE_prompt[ch];
      }
    }
    if( mean_vis > 0 && sum_pe > 0 ) {
      fCrsMuPhotons        = sum_pe / mean_vis;
      fCrsMuPhotonsPrompt  = sum_pe_pr / mean_vis;
      if( fCorrectQuenching )
        fCrsMuPhotons = CorrectForQuenching( fCrsMuPhotons, fCrsMuPhotonsPrompt, fMuContamCorr_EffTau[1]);
      hTrue_ResL_CrsMu          ->Fill( (fCrsMuPhotons-fTrue_CrsMuPhotons)/fTrue_CrsMuPhotons);
    }
    
    
    // Recomb and calibration plots (exclude muons that pass through wire planes)
    if( fCrsMuCharge > 0 && fCrsMuPhotons > 0 && fCrsMuLength > 10. 
      && fCrsMuVertex_X > 2. && fCrsMuEnd_X > 2. 
      && fCrsMuHitFrac > 0.95 ){
      float QoverL = fCrsMuCharge / fCrsMuPhotons;
      float R = ( 1.+fExcRatio ) / (1. + 1./QoverL );
      hQoverLCrsMuTrk->Fill( QoverL );
      hRecombCrsMuTrk->Fill( R );
      if( fCrsMuPitch[0] < 1.2 ) {
        hCrsMuADCPerCm_Pl0->Fill( fCrsMuIntegral_Pl0 / fCrsMuLength );
        hCrsMuQPerCm_Pl0  ->Fill( fCrsMuCharge_Pl0 / fCrsMuLength );
      }
      if( fCrsMuPitch[1] < 1.2 ) {
        hCrsMuADCPerCm    ->Fill( fCrsMuIntegral / fCrsMuLength );
        hCrsMuQPerCm      ->Fill( fCrsMuCharge / fCrsMuLength );
      }
      hCrsMuLPerCm      ->Fill( fCrsMuPhotons / fCrsMuLength );
    }
   
    hCrsMuLength  ->Fill( fCrsMuLength );
    hTrue_CrsMuLength->Fill( fTrue_MuTrackLength );
    hTrue_CrsMuQPerCm         ->Fill( fTrue_CrsMuCharge / fCrsMuLength );
    hTrue_CrsMuLPerCm         ->Fill( fTrue_CrsMuPhotons / fCrsMuLength );
    hTrue_ResQ_CrsMu[1]          ->Fill( (fCrsMuCharge - fTrue_CrsMuCharge) / fTrue_CrsMuCharge );
    hTrue_ResQ_CrsMu[0]          ->Fill( (fCrsMuCharge_Pl0 - fTrue_CrsMuCharge) / fTrue_CrsMuCharge );
    if( fCrsMuHitFrac > 0.95 ) {
      hTrue_ResE_CrsMu        ->Fill( ( fCrsMuEnergy - fTrue_MuEnergyDep ) / fTrue_MuEnergyDep );
      hTrue_ResQ_CrsMuTrk[1]          ->Fill( (muTrkQ - fTrue_CrsMuCharge) / fTrue_CrsMuCharge);
    }
//    hTrue_ResQ_CrsMuTrk[1]          ->Fill( ( muTrkQ - fTrue_CrsMuCharge ) / fTrue_CrsMuCharge );
//    hTrue_ResQ_CrsMuTrk[0]          ->Fill( (fCrsMuCharge_Pl0 - fTrue_CrsMuCharge_Trk) / fTrue_CrsMuCharge_Trk );
//    hTrue_ResQ_CrsMu[1]          ->Fill( (fCrsMuCharge - fTrue_TotalChargeDep) / fTrue_TotalChargeDep );
//    hTrue_ResQ_CrsMu[0]          ->Fill( (fCrsMuCharge_Pl0 - fTrue_TotalChargeDep) / fTrue_TotalChargeDep );
    

  }//<-- end calibration event information






  // *******************************************************************
  // Michel electron charge/shower reconstruction. If a stopping muon
  // track was identified AND a double-pulse was seein the optical 
  // waveforms, do Michel electron charge clustering in 2D/3D.
  // ******************************************************************
  
  // Define MichelCluster objects for each plane
  MichelCluster cl_pl0;
  MichelCluster cl_pl1;

  if( fMuTrackIndex >= 0 && fMichelOpticalID ) {  

    // ==========================================================
    // Fill histogram with track endpoints in drift direction, 
    // corrected for the measured muon decay time offset
    if( !fIsMC ) {
      hTrackNode_X_MuTrk->Fill( fMuTrackVertex_X + fDriftVelocity[0]*(fDecayTime/1000.) );
    }

    // ===========================================================
    // Perform clustering on both wireplanes
    Clustering( 0, cl_pl0 );
    Clustering( 1, cl_pl1 );

    // It's expected that event topology will not always be resolvable
    // on both planes. So if the muon-electron boundary is identified
    // successfully on one plane but not the other, then repeat clustering
    // on the failed plane without the requirements on linearity. Instead,
    // only require that the boundary hit be close in drift time ('x') to
    // the found boundary on the other plane.
    bool recluster = false;
    if( cl_pl1.bnd_i >= 0 && cl_pl0.bnd_i < 0 ) {
      MichelCluster newcluster;
      cl_pl0 = newcluster;
      Clustering(0, cl_pl0, cl_pl1.muEnd2D_X );
      recluster = true;
    } else 
    if( cl_pl1.bnd_i < 0 && cl_pl0.bnd_i >= 0 ) {
      MichelCluster newcluster;
      cl_pl1 = newcluster;
      Clustering(1, cl_pl1, cl_pl0.muEnd2D_X );
      recluster = true;
    }

    // ===========================================================
    // Save params to class member variables
    //...collection plane (default energy plane)
    fCovAtBnd               = cl_pl1.minCovAtBnd;
    fTruncMeanSlope         = cl_pl1.truncMeanSlope;
    fMuEnd2D_X              = cl_pl1.muEnd2D_X; 
    fMuEnd2D_W              = cl_pl1.muEnd2D_W;
    fElDir2D_X              = cl_pl1.elDir2D_X; 
    fElDir2D_W              = cl_pl1.elDir2D_W; 
    fClusterSize            = cl_pl1.cluster.size();
    fMuClusterSize          = cl_pl1.cluster_mu.size();
    fElClusterSize          = cl_pl1.cluster_el.size();
    fFracMuHitsLinear       = cl_pl1.fracMuHitsLinear;
    fMuAveLinearity         = cl_pl1.muAveLinearity;
    fMuClusterHitsEndFit    = cl_pl1.nPtsMuFit;
    fExtraHits              = cl_pl1.extraHits;
    fDecayAngle2D           = cl_pl1.decayAngle2D;
    fElShowerSize           = cl_pl1.shower.size();
    fElShowerFrac           = cl_pl1.elShowerFrac;
    //...induction plane
    fCovAtBnd_Pl0           = cl_pl0.minCovAtBnd;
    fTruncMeanSlope_Pl0     = cl_pl0.truncMeanSlope;
    fMuEnd2D_X_Pl0          = cl_pl0.muEnd2D_X; 
    fMuEnd2D_W_Pl0          = cl_pl0.muEnd2D_W; 
    fElDir2D_X_Pl0          = cl_pl0.elDir2D_X; 
    fElDir2D_W_Pl0          = cl_pl0.elDir2D_W; 
    fClusterSize_Pl0        = cl_pl0.cluster.size();
    fMuClusterSize_Pl0      = cl_pl0.cluster_mu.size();
    fElClusterSize_Pl0      = cl_pl0.cluster_el.size();
    fFracMuHitsLinear_Pl0   = cl_pl0.fracMuHitsLinear;
    fMuAveLinearity_Pl0     = cl_pl0.muAveLinearity;
    fMuClusterHitsEndFit_Pl0= cl_pl0.nPtsMuFit;
    fExtraHits_Pl0          = cl_pl0.extraHits;
    fDecayAngle2D_Pl0       = cl_pl0.decayAngle2D;
    fElShowerSize_Pl0       = cl_pl0.shower.size();
    fElShowerFrac_Pl0       = cl_pl0.elShowerFrac;
    
    hCovAtBnd           ->Fill( fCovAtBnd );
  
    
    LOG_VERBATIM("MichelAna")
    <<"DONE CLUSTERING: \n"
    <<" coll. plane cluster size = "<<fClusterSize<<" (mu: "<<fMuClusterSize<<", el: "<<fElClusterSize<<")\n"
    <<" frac mu hits linear = "<<fFracMuHitsLinear<<"\n"
    <<" el shower completeness frac = "<<fElShowerFrac;
   
    // Find projected distance to wireplanes
    if( fElDir2D_X < 0. ) {
      float dw = (fElDir2D_W / fElDir2D_X)*fMuEnd2D_X;
      fProjDist2DToWires = sqrt( pow(fMuEnd2D_X,2) + pow(dw,2));
    }
   
    // If a cluster boundary was found, set flag
    if( cl_pl1.bnd_i >= 0 ) isClusterBndFound = true;
      
    // Record fractional drop-off for each plane
    if( cl_pl1.muEndHit >= 0 && cl_pl1.cluster_el.size() > 0 ) {
      fFracDropoffAtBnd = 
        -1.*(fHitCharge[ cl_pl1.muEndHit ] - fHitCharge[ cl_pl1.cluster_el.at(0) ])/fHitCharge[cl_pl1.muEndHit];
    }
    if( cl_pl0.muEndHit >= 0 && cl_pl0.cluster_el.size() > 0 ) {
      fFracDropoffAtBnd_Pl0 = 
        -1.*(fHitCharge[ cl_pl0.muEndHit ] - fHitCharge[ cl_pl0.cluster_el.at(0) ])/fHitCharge[cl_pl0.muEndHit];
    }
    hFracDropoffAtBnd->Fill( fFracDropoffAtBnd);

    // ===========================================================
    // If cluster boundary (ie, muon end hit) found on both planes, then 
    // we can calculate the 3D muon endpoint position and get an idea of 
    // any systematic offsets in timing between the planes
    if( cl_pl1.muEndHit >= 0 && cl_pl0.muEndHit >= 0 ) {
      int hit0 = cl_pl0.muEndHit;
      int hit1 = cl_pl1.muEndHit;
      fMuEndHitTimeDiff = fHitT[hit1] - fHitT[hit0];
      if( !recluster ) hMuEndHitTimeDiff -> Fill( fMuEndHitTimeDiff );
      double dt = fMaxHitTimeDiffFac * fSamplingRate * sqrt(pow(fHitRMS[hit1],2) + pow(fHitRMS[hit0],2) );
      if( fabs(fMuEndHitTimeDiff) < dt ){
        double y = -99.;
        double z = -99.;
        fGeo->ChannelsIntersect( fHitWire.at(hit0), fHitWire.at(hit1), y, z );
        fMuEnd3D.SetXYZ((fHitX.at(hit1)+fHitX.at(hit0))/2., y, z);
        fMuEnd3D_X = fMuEnd3D.X();
        fMuEnd3D_Y = y;
        fMuEnd3D_Z = z;
          
        if( fIsMC ) {
          TVector3 delta          = (fMuEnd3D - fTrue_MuTrackEnd);
          hTrue_MuClsEndRes       ->Fill( delta.Mag() );
          hTrue_MuClsEndRes_X     ->Fill( delta.X() );
          hTrue_MuClsEndRes_Y     ->Fill( delta.Y() );
          hTrue_MuClsEndRes_Z     ->Fill( delta.Z() );
          if( fMuTrackEndDir.Mag() > 0. ) {
            hTrue_MuClsEndResDir   ->Fill( delta.Dot(fMuTrackEndDir) );
          }
        }

      }//endif mu hit time diff < cut
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
        fMuCharge       = cl_pl1.muCharge;
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

  /*
    // ===========================================================
    // Calculate energy resolution 
    float res = 0.;
    if( fIsMC && fElShowerEnergy > 0 ) {
      res = (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep;
      LOG_VERBATIM("MichelAna")<<"Energy resolution: "<<res;
    }
    */
   
   
    
    // ===========================================================
    // Match hits between the two planes to form 3D points 
    // ===========================================================
    if( fElShowerSize > 0 && fElShowerSize_Pl0 > 0 ) {
      
      std::vector<int> shower1 = cl_pl1.shower;
      std::vector<int> shower0 = cl_pl0.shower;
      std::vector<bool> flag_1(cl_pl1.shower.size(), true);
      std::vector<bool> flag_0(cl_pl0.shower.size(), true);
      LOG_VERBATIM("MichelAna")
      <<"Beginning to match hits between two planes: "<<shower1.size()<<", "<<shower0.size()<<" hits each";
      
      // ---------------------------------------------------
      // Create array of hit time differences 
      std::vector< std::vector<float> > diff_array (shower1.size(), std::vector<float> (shower0.size(), 0));
      std::vector< std::vector<float> > dt_array (shower1.size(), std::vector<float> (shower0.size(), 0));
      for(size_t i=0; i<shower1.size(); i++){
        for(size_t j=0; j<shower0.size(); j++){
          diff_array.at(i).at(j)  = fHitT[shower1.at(i)] - fHitT[shower0.at(j)]; 
          dt_array.at(i).at(j)    = fSamplingRate* sqrt(pow(fHitRMS[shower1.at(i)],2) + pow(fHitRMS[shower0.at(j)],2) );
          hHitTimeDiff->Fill(diff_array.at(i).at(j));
          hHitXDiff   ->Fill(fHitX[shower1.at(i)] - fHitX[shower0.at(j)]);
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
           if(       fabs(diff_array.at(i).at(j)) < fMaxHitTimeDiffFac*dt_array.at(i).at(j)
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
      LOG_VERBATIM("MichelAna")
      <<"Found "<<Pts.size()<<" 3D points (frac = "<<fFracHits3D<<")";
     
      // -------------------------------------------------------
      // Make 3D projections of muon and electron
      

      // ====================================================================
      // If we have some 3D points reconstructed, we can do more stuff...

      if( fNumPts3D > 0 && fElShowerEnergy > 0 ) { 
    
        // -------------------------------------------------------------
        // Find average charge-weighted optical visibility of the 3D points
        fElShowerPhotons = 0;
        fElShowerPhotonsPrompt = 0;
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
          if( fUsePmtForCalo[ch] ) fElShowerVis        += fElShowerVisCh[ch];
          LOG_VERBATIM("MichelAna")<<"Ave vis for PMT "<<ch<<" = "<<sum_vis / sum_W;
        }//<-- end loop over PMTs
        if( sum_W > 0 ) {
          fElShowerCentroid_X = sum_x / sum_W;
          fElShowerCentroid_Y = sum_y / sum_W;
          fElShowerCentroid_Z = sum_z / sum_W;
          fElShowerCentroid.SetXYZ(fElShowerCentroid_X, fElShowerCentroid_Y, fElShowerCentroid_Z);
        }
      

        // ---------------------------------------------------------
        // Calculate L and Q+L energy
        //   E = (N_i + N_ex) * W_ph
        //   N_i + N_ex = N_e + N_ph 
        if( fElShowerVis > 0 ) {
          fElShowerPhotons        = fPE_comb_qc / fElShowerVis;
          fElShowerPhotonsPrompt  = fPE_prompt_comb / fElShowerVis;
          fElShowerEnergyQL       = ( fElShowerPhotons + fElShowerCharge )*fWph;
        }
     
        // ----------------------------------------------------------
        // Determine direction in 3D for the electron shower
        if( fElShowerCentroid_X > 0. ) {
          // use 3D muon endpoint if it was defined. otherwise, use the 3D
          // stopping track endpoint
          float x = -9.;
          float y = -999.;
          float z = -999.;
          if( fMuEnd3D.X() > 0. ) {
            fElDir3D = (fElShowerCentroid - fMuEnd3D);
            x = fMuEnd3D.X();
            y = fMuEnd3D.Y();
            z = fMuEnd3D.Z();
          } else 
          if( fMuTrackEnd.X() > 0. ) {
            fElDir3D = (fElShowerCentroid - fMuTrackEnd);
            x = fMuTrackEnd.X();
            y = fMuTrackEnd.Y();
            z = fMuTrackEnd.Z();
          }
          if( fElDir3D.Mag() > 0 ) {
            fElDir3D.SetMag(1.);
            if( fElDir3D.X() < 0. ) {
              TVector3 dirToWires(-1., 0., 0.);
              float cos_theta = cos(fElDir3D.Angle(dirToWires));
              if( cos_theta > 0. && x > 0. ) {
                fProjDist3DToWires  = x /  cos_theta;
                fProjPtOnWires_Y    = y + fabs(x/fElDir3D.X()) * fElDir3D.Y();
                fProjPtOnWires_Z    = z + fabs(x/fElDir3D.X()) * fElDir3D.Z();
              }
            }
          }
        
        }
      
      }//<-- endif Npts3D > 0 
    }//<-- endif shower size > 0 on both planes
  }//<-- endif tagged muon + Michel pulse found
   
   
   
    
    
  // *******************************************************************
  // Analysis of event reduction from these quality cuts, and filling histograms
  // *******************************************************************
  
  bool goodShower2D     = false;
  bool goodShower3D     = false;
  
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

            if( fTruncMeanSlope > fMinTruncMeanSlope ) {
            
              // -----------------------------------------------------------------
              // muon linearity cuts
              if(     fFracMuHitsLinear > fMinFracMuHitsLinear
                  &&  fMuAveLinearity > fMinMuLinearity ) {
        

                // ---------------------------------------------------------------
                // muon direction
                if( fMuClusterHitsEndFit >= fMinMuClusterHitsEndFit ) {
                  hEventCuts->Fill(7-1);
                  
                  hDecayAngle2D->Fill( fDecayAngle2D ); 
                  
                  // -------------------------------------------------------------
                  // decay angle cut
                  if( fDecayAngle2D > fMinDecayAngle2D && fDecayAngle2D < fMaxDecayAngle2D ) {
                  
                    if( fElShowerFrac > 0 ) {
                      hElShowerSize->Fill( fElShowerSize );
                      hElShowerFrac->Fill( fElShowerFrac );
                      hLeftoverHits->Fill( fElShowerSize*(1./fElShowerFrac - 1.) );
                      hElTrackFrac->Fill( fElTrackCharge / fElShowerCharge );
                    }

                    // -----------------------------------------------------------
                    // Shower frac cut
                    if( fElShowerFrac > fMinElShowerFrac && fElShowerEnergy > 0) {
                      hEventCuts->Fill(8-1);
                      
                      LOG_VERBATIM("MichelAna")
                      <<"Good shower!";
                      
                      goodShower2D = true;
                      fNumEvents_Shwr++;
                     
                      float trkRes = (fElTrackEnergy - fTrue_ElTrackEnergyDep)/fTrue_ElTrackEnergyDep; 
                      //float res = (fElShowerEnergy-fTrue_ElShowerEnergyDep)/fTrue_ElShowerEnergyDep;
                        
                      if( fIsMC && (trkRes < -0.5 || fElTrackEnergy < 4.) ) 
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
                      
                        LOG_VERBATIM("MichelAna")
                        <<"Good 3D shower!";
                             
                        hNumPts3D   ->Fill(fNumPts3D);
                        hFracHits3D  -> Fill(fFracHits3D);  

                        // N pts 3D > min
                        if( fNumPts3D > fMinNumPts3D ) {
                         
                          // Frac hits
                          if( fFracHits3D > fMinFracHits3D ) {
                            
                            goodShower3D = true;        
                            fNumEvents_Shwr3D++;

                            hProjDist3DToWires->Fill( fProjDist3DToWires );
                            hProjPtOnWires->Fill( fProjPtOnWires_Z, fProjPtOnWires_Y );
                            if( fIsMC && fTrue_ProjDist3DToWires > 0. ) 
                              hProjDist3DToWiresRes->Fill( (fProjDist3DToWires-fTrue_ProjDist3DToWires)/fTrue_ProjDist3DToWires );

                          }//frac hits

                        }// Npts3D
                        
                      }// good shower on induction 

                    }// shower frac
                  } //decay angle   
                }// mu dir
              }//mu lin
            } // trunc mean slope
          }//extra hits
        }// cluster size
      }//cluster bnd
    }// op ID
  }// stp trk

  // *******************************************************************
  // Make a few histograms (the rest will be made in a separate macro).
  // *******************************************************************
  for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){
    size_t ch   = fPulses[ipmt].OpChannel(); 
    if( fDecayTime > fdTcut && goodShower2D ) {
          hElShowerEnergyVsPromptPE[ch]->Fill( fElShowerEnergy, fPE_prompt[ch] );
          hPE_prompt_dTcut_shwr[ch]       ->Fill(fPE_prompt[ch]);
          hPE_total_dTcut_shwr[ch]        ->Fill(fPE_total[ch]);
          hPE_totalRaw_dTcut_shwr[ch]     ->Fill(fPE_totalRaw[ch]);
          hTrue_PE_total_dTcut_shwr[ch]   ->Fill(fTrue_PE_total[ch]);
          hTrue_PE_total_TrueVsReco[ch]   ->Fill(fTrue_PE_total[ch], fPE_total[ch]);
          if( fTrue_PE_total[ch] > 0 ) hTrue_ResPEch[ch]->Fill( (fPE_total[ch] - fTrue_PE_total[ch] ) / fTrue_PE_total[ch] );
    }
  }//<-- endloop PMTs
 
  
  // --------------------------------------------------------------------
  // If the identified stopping mu track's residual range was found to be
  // in the right order (relative to the likely endpoint) and its endpoint
  // X coordinate is consistent with the endpoint found in the clustering
  // stage (2D), then fill residual range histograms.
  if( fDecayTime > 300. && fMuTrackIsCaloOrdered ) {
     
    TVector3 delta;
    float delta_p = -99.;
    float delta_t = -99.;

    // Which angle bin does this track fall in?
    int   angleBin = -9;
    for(size_t i=0; i<4; i++){
      if( fMuTrackInclineAngle*RAD_TO_DEG >= fAngleBins_X1[i] ) {
        angleBin = i;
        break;
      }
    }

    // Look at accuracy of track endpoint vs. true endpoint
    if( fIsMC ) {
        delta = (fMuTrackEnd-fTrue_MuTrackEnd);
        hTrue_MuTrkEndRes   ->Fill( delta.Mag() );
        hTrue_MuTrkEndRes_X ->Fill( delta.X()   );
        hTrue_MuTrkEndRes_Y ->Fill( delta.Y()   );
        hTrue_MuTrkEndRes_Z ->Fill( delta.Z()   );
        if( fMuTrackEndDir.Mag() > 0. ){
          delta_p = delta.Dot(fMuTrackEndDir);
          delta_t = sqrt( pow(delta.Mag(),2) - pow(delta.Dot(fMuTrackEndDir),2));
          hTrue_MuTrkEndResDir  ->Fill( delta_p );
          hTrue_MuTrkEndResTrans->Fill( delta_t );
        }
    }

    // Only select very well-reconstructed 3D tracks by requiring
    // endpoint to be close to that found in clustering earlier
    float diff = (fMuTrackEnd-fMuEnd3D).Mag();
    hMuEndDiff3D->Fill( diff );
    if( diff >= 0. && diff < 0.2 ) {
     
      // Keep track of incline angles of selected tracks
      hMuTrackInclineAngle_dEdx->Fill( fMuTrackInclineAngle*RAD_TO_DEG );
      
      if( delta_t >= 0. ) {
        hTrue_MuTrkEndResDir_cut  ->Fill( delta_p );
        hTrue_MuTrkEndResTrans_cut->Fill( delta_t );
      }

      // Now fill dQ/dx and dE/dx information
      for(size_t ipl=0; ipl<2; ipl++){
        if( fvMuTrkResRange[ipl].size() < 5 ) continue;
        for( size_t i=0; i<fvMuTrkdEdx[ipl].size(); i++){
          hMuResRangePitch[ipl]       ->Fill( fvMuTrkPitch[ipl].at(i) );
          float cf = 1.0;
          if( !fIsMC ) cf = exp((fDecayTime/1000.)/fElectronLifetime);
          fvMuTrkdQdx[ipl].at(i)    *= cf;
          fvMuTrkdADCdx[ipl].at(i)  *= cf;
          fvMuTrkdEdx[ipl].at(i)    = detprop->ModBoxCorrection( fvMuTrkdQdx[ipl].at(i) );
          if( fvMuTrkPitch[ipl].at(i) < 1.2 ) {
            float rr = fvMuTrkResRange[ipl].at(i) + fMuResRangeOffset;
            if( fMuTrackInclineAngle*RAD_TO_DEG > 20. ) {
              hMuResRangeVsdEdx[ipl]    ->Fill(rr, fvMuTrkdEdx[ipl].at(i));
              hMuResRangeVsdQdx[ipl]    ->Fill(rr, fvMuTrkdQdx[ipl].at(i));
              hMuResRangeVsdADCdx[ipl]  ->Fill(rr, fvMuTrkdADCdx[ipl].at(i));
            }
            if( angleBin >= 0 && ipl == 1 ) {
              hMuResRangeVsdEdx_AngleBins[angleBin]    ->Fill(rr, fvMuTrkdEdx[ipl].at(i));
              hMuResRangeVsdQdx_AngleBins[angleBin]    ->Fill(rr, fvMuTrkdQdx[ipl].at(i));
              hMuResRangeVsdADCdx_AngleBins[angleBin]  ->Fill(rr, fvMuTrkdADCdx[ipl].at(i));
            } // angle bin plots
          } // trkPitch cut
        } // loop over dEdx points
      } // loop over planes
    
    }//endif endpoints match in time (x)
  }//endif muon track cal points are ordered w.r.t. res range

  /* 
  // Similarly, fill dQ/dx vs. res range info for the cluster too
  if( fvMuClusterResRange.size() > 0 ) {
    for(size_t i=0; i<fvMuClusterdEdx.size(); i++){
      hMuClusterResRangeVsdADCdx[1]->Fill(fvMuClusterResRange[i], fvMuClusterdADCdx[i]);
      hMuClusterResRangeVsdQdx[1]->Fill(fvMuClusterResRange[i], fvMuClusterdQdx[i]);
      hMuClusterResRangeVsdEdx[1]->Fill(fvMuClusterResRange[i], fvMuClusterdEdx[i]);
    }
  }
  if( fvMuClusterResRange_Pl0.size() > 0 ) {
    for(size_t i=0; i<fvMuClusterdEdx_Pl0.size(); i++){
      hMuClusterResRangeVsdADCdx[0]->Fill(fvMuClusterResRange_Pl0[i], fvMuClusterdADCdx_Pl0[i]);
      hMuClusterResRangeVsdQdx[0]->Fill(fvMuClusterResRange_Pl0[i], fvMuClusterdQdx_Pl0[i]);
      hMuClusterResRangeVsdEdx[0]->Fill(fvMuClusterResRange_Pl0[i], fvMuClusterdEdx_Pl0[i]);
    }
  }
  */

  if( fIsMC ) {
    if( fMuonID >= 0 ) {
      hTrue_MuTrackEnd_ZX->Fill( fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_X);
      hTrue_MuTrackEnd_ZY->Fill( fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_Y);
      hTrue_RecombFactor  ->Fill( fTrue_ElShowerChargeDep / fTrue_ElShowerIons, fTrue_ElShowerIons );
    }
  }

  if( goodShower2D ) {
    hElShowerEnergy     ->Fill(fElShowerEnergy);
    hElShowerCharge     ->Fill(fElShowerCharge);
    hElTrackCharge     ->Fill(fElTrackCharge);
    if( fMuEnd3D_X > 0. ) {
      hMuEnd3D_X->Fill( fMuEnd3D_X);
      hMuEnd3D_Y->Fill( fMuEnd3D_Y);
      hMuEnd3D_Z->Fill( fMuEnd3D_Z);
      hMuEnd3D_ZX->Fill( fMuEnd3D_Z, fMuEnd3D_X );
      hMuEnd3D_ZY->Fill( fMuEnd3D_Z, fMuEnd3D_Y );
    }
    if( goodShower3D && fDecayTime > fdTcut ) {
      hElShowerVis        ->Fill( fElShowerVis );
      hElShowerCentroid_X   ->Fill( fElShowerCentroid_X);
      hElShowerCentroid_Y   ->Fill( fElShowerCentroid_Y);
      hElShowerCentroid_Z   ->Fill( fElShowerCentroid_Z);
      hElShowerCentroid_ZX  ->Fill( fElShowerCentroid_Z, fElShowerCentroid_X);
      hElShowerCentroid_ZY  ->Fill( fElShowerCentroid_Z, fElShowerCentroid_Y);
      hElShowerCharge3D     ->Fill(fElShowerCharge);
      hElShowerPhotons    ->Fill(fElShowerPhotons);
      hElShowerEnergyQL   ->Fill(fElShowerEnergyQL);
      float R = ( 1.+fExcRatio ) / (1. + fElShowerPhotons/fElShowerCharge);
      hQoverL   ->Fill( fElShowerCharge/fElShowerPhotons );
      hRecomb   ->Fill( R ); 
    }

    if( fTrue_ElShowerEnergyDep > 0. ) {
      
      float trkRes = (fElTrackEnergy - fTrue_ElTrackEnergyDep) / fTrue_ElTrackEnergyDep;

      hTrue_ElEnergy                ->Fill(fTrue_ElEnergy);
      if( fTrue_MuSign == 1 ) hTrue_ElEnergyFree->Fill(fTrue_ElEnergy);
      if( fTrue_MuSign == -1 ) hTrue_ElEnergyCapture->Fill(fTrue_ElEnergy);
      hTrue_ElShowerChargeDep       ->Fill(fTrue_ElShowerChargeDep);
      hTrue_ElTrackChargeDep        ->Fill(fTrue_ElTrackChargeDep);
      hTrue_ElShowerEnergyDep       ->Fill(fTrue_ElShowerEnergyDep);
      hTrue_ElTrackEnergyDep        ->Fill(fTrue_ElTrackEnergyDep);
      hTrue_ElShowerPhotons         ->Fill(fTrue_ElShowerPhotons);
      hTrue_EnergyDepVsRecoEnergy   ->Fill(fTrue_ElShowerEnergyDep, fElShowerEnergy);
      hTrue_EnergyRes        ->Fill( (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
      hTrue_EnergyVsEnergyRes       ->Fill( fTrue_ElShowerEnergyDep, (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
      std::cout<<"Track energy: "<<fElTrackEnergy<<", true= "<< fTrue_ElTrackEnergyDep<<"\n";
      std::cout<<"trkRes = "<<trkRes<<"\n";
      hTrue_EnergyResTrk     ->Fill( trkRes );
      if( trkRes < -0.6 ) {
        std::cout<<"LOW TRACK RESOLUTION: "<<trkRes<<"\n";
      }
      hTrue_ResQ          ->Fill( (fElShowerCharge - fTrue_ElShowerChargeDep) / fTrue_ElShowerChargeDep);
      hTrue_ResQ_Mu          ->Fill( ( fMuCharge - fTrue_MuCharge ) / fTrue_MuCharge );
      if( fIsMC && fElShowerEnergyQL > 0 && goodShower3D && fDecayTime > fdTcut ) {
        hTrue_ResQ_Shwr3D          ->Fill( (fElShowerCharge - fTrue_ElShowerChargeDep) / fTrue_ElShowerChargeDep);
        hTrue_EnergyResQL ->Fill((fElShowerEnergyQL - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep);
        hTrue_EnergyRes_Shwr3D       ->Fill( (fElShowerEnergy - fTrue_ElShowerEnergyDep) / fTrue_ElShowerEnergyDep );
        hTrue_EnergyResTrk_Shwr3D     ->Fill( trkRes );
        hTrue_ResL        ->Fill( (fElShowerPhotons - fTrue_ElShowerPhotons) / fTrue_ElShowerPhotons );
        hTrue_ResVis      ->Fill( (fElShowerVis - fTrue_ElShowerVis) / fTrue_ElShowerVis );
        hTrue_ResPE       ->Fill( (fPE_comb - fTrue_PE) / fTrue_PE );
        hTrue_ResPEPrompt       ->Fill( (fPE_prompt_comb - fTrue_PEPrompt) / fTrue_PEPrompt );
        hTrue_ResPE0       ->Fill( (fTrue_PE - fTrue_PE0) / fTrue_PE0 );
      }
    }

  }// end goodShower2D
  // Done making histograms 
 


  // *******************************************************
  // Fill the tree  
  // *******************************************************
  if( fMuTrackIndex >= 0 && fMichelOpticalID && isClusterBndFound ) fTree->Fill();
  if( isCalibrationEvent ) fTreeCrsMu->Fill();
  //if( isCalibrationEvent || (fMuTrackIndex >= 0 && fMichelOpticalID && isClusterBndFound) ) fTree->Fill();

  // Make waveform graphs
  if( fElShowerEnergy > 0 && fMichelOpticalID ) MakeWfmGraphs(e);

  LOG_VERBATIM("MichelAna")
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

//#######################################################################################
float MichelAna::CorrectForQuenching(float L_tot, float L_prompt, float tau){
  float quenchCorrFactor = 1.0;
  float estimatedTau    = TauConversion(fMuContamCorr_EffTau[1]); // formula to account for TPB convolution effects (use ETL)

  // first correct for limited integration window 
  //L_tot = PeIntegrationCorrection(L_prompt, L_tot, fFullWindow, TauConversion(tau));

  // scale up light to compensate for quenching
  if( estimatedTau < fTauT ) quenchCorrFactor = (fTauT*exp(-100./fTauT)) / (estimatedTau*exp(-100./estimatedTau));
  float lateLight = (L_tot - L_prompt)*quenchCorrFactor;
  return L_prompt + lateLight;
}

//########################################################################################
float MichelAna::PeIntegrationCorrection(float pe_pr, float pe_total, float window, float tau){
  float late_light = pe_total - pe_pr;
  float totalExp = exp(-100./tau);
  float totalCol = exp(-100./tau) - exp(-window/tau);
  late_light *= totalExp/totalCol;
  return pe_pr + late_light;
}

//#####################################
float MichelAna::TauConversion( float tau_eff ) {
//  return 1.06*tau_eff - 264.;
  // values computed 5/24/2018
  return 0.965*tau_eff - 41.7;
}

//########################################################################################
void MichelAna::endJob()
{
  if( fNumEvents_Shwr == 0 ) {
    hTrue_EnergyRes->Fill(-9);
    hElShowerEnergy->Fill(-9); 
  }
    

  std::cout
  <<"============================================================\n"
  <<"            Ending MichelAna Job\n"
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
  <<"  LinThresh                              : "<<fLinThresh<<"\n"
  <<"  LinThreshMuFit                         : "<<fLinThreshMuFit<<"\n"
  <<"  MinClusterSize                         : "<<fMinClusterSize<<"\n"
  <<"  MinMuClusterSize                       : "<<fMinMuClusterSize<<"\n"
  <<"  MinElClusterSize                       : "<<fMinElClusterSize<<"\n"
  <<"  MinElShowerFrac                        : "<<fMinElShowerFrac<<"\n"
  <<"  MinFracMuHitsLinear                    : "<<fMinFracMuHitsLinear<<"\n"
  <<"\n"
  <<"  RFactor                                : "<<fRecombFactor<<"\n"
  <<"  RFactorTrack                           : "<<fRecombFactorTrack<<"\n"
  <<"  RFactorShower                          : "<<fRecombFactorShower<<"\n"
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
  <<"  Ave Michel charge reconstructed        : "<<hElShowerCharge->GetMean()<<" e-\n" 
  <<"     - track-like component              : "<<hElTrackCharge->GetMean()<<" e-\n" 
  <<"     - shower-like component             : "<<hElShowerCharge->GetMean()-hElTrackCharge->GetMean()<<" e-\n" 
  <<"  Ave Michel charge (3D cut)             : "<<hElShowerCharge3D->GetMean()<<" e-\n" 
  <<"  Ave energy reconstructed (Q)           : "<<hElShowerEnergy->GetMean()<<" MeV\n" 
  <<"  Ave energy reconstructed (Q+L)         : "<<hElShowerEnergyQL->GetMean()<<" MeV\n" 
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
  <<"    - Num entries in ave 'true' wfm      : "<<aveWfmCounter[ch]<<"\n"
  <<"    - Num entries in ave muon wfm        : "<<aveWfmCounter_mu[ch]<<"\n"
  <<"    - Num entries in ave electron wfm    : "<<aveWfmCounter_el[ch]<<"\n"
  <<"    - Num entries in ave crs mu wfm      : "<<aveWfmCounter_crsmu[ch]<<"\n"
  <<"\n\n";
  
  
  
  // =====================================================================
  // Average waveform analysis and fits
  TF1 fit("fit","[0]*exp(-x/[1])",100,7000);
  TF1 fit2("fit2","[0]*(exp(-x/[1])+[2]*exp(-x/3550))",100,7000);
  float t1 = 400;
  float t2 = 2000;
  float a, tau, a2, tau2, b2;
  float da, dtau, da2, dtau2, db2;
  float int0, int100, int2000, int7000;

  // --------------------------------------------------------------------
  // MC "truth" waveform
  if( aveWfmCounter[ch] > 0 ) {
    hAveWfm[ch]->Scale(1./float(aveWfmCounter[ch]));
    hAveWfm[ch]->SetOption("hist");
    // single exponential fit (400 ns to 2 us)
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAveWfm[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm[ch]->Fit("fit","RNQ0");
    a     = fit.GetParameter(0);
    da    = fit.GetParError(0);
    tau   = fit.GetParameter(1);
    dtau  = fit.GetParError(1); 
    // Now fit with full range adding long TPB component (400ns-7us)
    fit2.SetRange(t1,7000);
    fit2.SetParameter(0,fit.GetParameter(0));
    fit2.SetParameter(1,fit.GetParameter(1));
    fit2.SetParameter(2,0.05);
    fit2.SetParLimits(2,0.,1.0);
    hAveWfm[ch]->Fit("fit2","RNQ0");
    a2    = fit2.GetParameter(0);
    da2   = fit2.GetParError(0);
    tau2  = fit2.GetParameter(1);
    dtau2 = fit2.GetParError(1); 
    b2    = fit2.GetParameter(2);
    db2   = fit2.GetParError(2);
    // Figure out ratios of prompt light to total
    int0    = hAveWfm[ch]->GetCumulative()->GetBinContent(hAveWfm[ch]->GetCumulative()->GetXaxis()->FindBin(-10.));
    int100  = hAveWfm[ch]->GetCumulative()->GetBinContent(hAveWfm[ch]->GetCumulative()->GetXaxis()->FindBin(100.)) - int0;
    int2000 = hAveWfm[ch]->GetCumulative()->GetBinContent(hAveWfm[ch]->GetCumulative()->GetXaxis()->FindBin(2000.)) - int0;
    int7000 = hAveWfm[ch]->GetCumulative()->GetBinContent(hAveWfm[ch]->GetCumulative()->GetXaxis()->FindBin(7000.-1.)) - int0;
    std::cout
    <<"  -------------------------------------------\n"
    <<"  Ave MC waveform:\n"
    <<"   Exp fit, 400ns-2us:\n"
    <<"     - a     = "<<a<<" +/- "<<da<<"\n"
    <<"     - tau   = "<<tau<<" +/- "<<dtau<<" ns\n"
    <<"   Exp + long TPB fit, 400ns-7us:\n"
    <<"     - a     = "<<a2<<" +/- "<<da2<<"\n"
    <<"     - b     = "<<b2<<" +/- "<<db2<<"\n"
    <<"     - tau   = "<<tau2<<" +/- "<<dtau2<<" ns\n"
    <<"   Prompt integral & ratios:\n"
    <<"     - 100ns      = "<<int100<<"\n" 
    <<"     - 100ns/2us  = "<<int100/int2000<<"\n"
    <<"     - 100ns/7us  = "<<int100/int7000<<"\n";
  }
  
  // --------------------------------------------------------------------
  // Muon pulse waveform
  if( aveWfmCounter_mu[ch] > 0 ) {
    hAveWfm_mu[ch]->Scale(1./float(aveWfmCounter_mu[ch]));
    hAveWfm_mu[ch]->SetOption("hist");
    // single exponential fit (400 ns to 2 us)
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm_mu[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAveWfm_mu[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm_mu[ch]->Fit("fit","RNQ0");
    a     = fit.GetParameter(0);
    da    = fit.GetParError(0);
    tau   = fit.GetParameter(1);
    dtau  = fit.GetParError(1); 
    // Now fit with full range adding long TPB component (400ns-7us)
    fit2.SetRange(t1,t2);
    fit2.SetParameter(0,fit.GetParameter(0));
    fit2.SetParameter(1,fit.GetParameter(1));
    fit2.SetParameter(2,0.05);
    fit2.SetParLimits(2,0.,1.0);
    hAveWfm_mu[ch]->Fit("fit2","RNQ0");
    a2    = fit2.GetParameter(0);
    da2   = fit2.GetParError(0);
    tau2  = fit2.GetParameter(1);
    dtau2 = fit2.GetParError(1); 
    b2    = fit2.GetParameter(2);
    db2   = fit2.GetParError(2);
    // Figure out ratios of prompt light to total
    int0    = hAveWfm_mu[ch]->GetCumulative()->GetBinContent(hAveWfm_mu[ch]->GetCumulative()->GetXaxis()->FindBin(-10.));
    int100  = hAveWfm_mu[ch]->GetCumulative()->GetBinContent(hAveWfm_mu[ch]->GetCumulative()->GetXaxis()->FindBin(100.)) - int0;
    int2000 = hAveWfm_mu[ch]->GetCumulative()->GetBinContent(hAveWfm_mu[ch]->GetCumulative()->GetXaxis()->FindBin(2000.-1.)) - int0;
    std::cout
    <<"  -------------------------------------------\n"
    <<"  Muon pulse waveform:\n"
    <<"   Exp fit, 400ns-2us:\n"
    <<"     - a     = "<<a<<" +/- "<<da<<"\n"
    <<"     - tau   = "<<tau<<" +/- "<<dtau<<" ns\n"
    <<"   Exp + long TPB fit, 400ns-2us:\n"
    <<"     - a     = "<<a2<<" +/- "<<da2<<"\n"
    <<"     - b     = "<<b2<<" +/- "<<db2<<"\n"
    <<"     - tau   = "<<tau2<<" +/- "<<dtau2<<" ns\n"
    <<"   Prompt integral & ratios:\n"
    <<"     - 100ns      = "<<int100<<"\n" 
    <<"     - 100ns/2us  = "<<int100/int2000<<"\n";
  }
  
  // --------------------------------------------------------------------
  // Electron waveform
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
    // single exponential fit (400 ns to 2 us)
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm_el[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAveWfm_el[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm_el[ch]->Fit("fit","RNQ0");
    a     = fit.GetParameter(0);
    da    = fit.GetParError(0);
    tau   = fit.GetParameter(1);
    dtau  = fit.GetParError(1); 
    // Now fit with full range adding long TPB component (400ns-7us)
    fit2.SetRange(t1,7000);
    fit2.SetParameter(0,fit.GetParameter(0));
    fit2.SetParameter(1,fit.GetParameter(1));
    fit2.SetParameter(2,0.05);
    fit2.SetParLimits(2,0.,1.0);
    hAveWfm_el[ch]->Fit("fit2","RNQ0");
    a2    = fit2.GetParameter(0);
    da2   = fit2.GetParError(0);
    tau2  = fit2.GetParameter(1);
    dtau2 = fit2.GetParError(1); 
    b2    = fit2.GetParameter(2);
    db2   = fit2.GetParError(2);
    // Figure out ratios of prompt light to total
    int0    = hAveWfm_el[ch]->GetCumulative()->GetBinContent(hAveWfm_el[ch]->GetCumulative()->GetXaxis()->FindBin(-10.));
    int100  = hAveWfm_el[ch]->GetCumulative()->GetBinContent(hAveWfm_el[ch]->GetCumulative()->GetXaxis()->FindBin(100.)) - int0;
    int2000 = hAveWfm_el[ch]->GetCumulative()->GetBinContent(hAveWfm_el[ch]->GetCumulative()->GetXaxis()->FindBin(2000.)) - int0;
    int7000 = hAveWfm_el[ch]->GetCumulative()->GetBinContent(hAveWfm_el[ch]->GetCumulative()->GetXaxis()->FindBin(7000.-1.)) - int0;
    std::cout
    <<"  -------------------------------------------\n"
    <<"  Electron pulse waveform:\n"
    <<"   Exp fit, 400ns-2us:\n"
    <<"     - a     = "<<a<<" +/- "<<da<<"\n"
    <<"     - tau   = "<<tau<<" +/- "<<dtau<<" ns\n"
    <<"   Exp + long TPB fit, 400ns-7us:\n"
    <<"     - a     = "<<a2<<" +/- "<<da2<<"\n"
    <<"     - b     = "<<b2<<" +/- "<<db2<<"\n"
    <<"     - tau   = "<<tau2<<" +/- "<<dtau2<<" ns\n"
    <<"   Prompt integral & ratios:\n"
    <<"     - 100ns      = "<<int100<<"\n" 
    <<"     - 100ns/2us  = "<<int100/int2000<<"\n"
    <<"     - 100ns/7us  = "<<int100/int7000<<"\n";
  }
  
  // --------------------------------------------------------------------
  // Crossing muon waveform
  if( aveWfmCounter_crsmu[ch] > 0 ) {
    hAveWfm_crsmu[ch]->Scale(1./float(aveWfmCounter_crsmu[ch]));
    hAveWfm_crsmu[ch]->SetOption("hist");
    // single exponential fit (400 ns to 2 us)
    fit.SetRange(t1,t2);
    fit.SetParameter(0, hAveWfm_crsmu[ch]->GetMaximum()/10. );
    fit.SetParLimits(0, 0, hAveWfm_crsmu[ch]->GetMaximum() );
    fit.SetParameter(1, 1200 );
    fit.SetParLimits(1, 500, 2000);
    hAveWfm_crsmu[ch]->Fit("fit","RNQ0");
    a     = fit.GetParameter(0);
    da    = fit.GetParError(0);
    tau   = fit.GetParameter(1);
    dtau  = fit.GetParError(1); 
    // Now fit with full range adding long TPB component (400ns-7us)
    fit2.SetRange(t1,7000);
    fit2.SetParameter(0,fit.GetParameter(0));
    fit2.SetParameter(1,fit.GetParameter(1));
    fit2.SetParameter(2,0.05);
    fit2.SetParLimits(2,0.,1.0);
    hAveWfm_crsmu[ch]->Fit("fit2","RNQ0");
    a2    = fit2.GetParameter(0);
    da2   = fit2.GetParError(0);
    tau2  = fit2.GetParameter(1);
    dtau2 = fit2.GetParError(1); 
    b2    = fit2.GetParameter(2);
    db2   = fit2.GetParError(2);
    // Figure out ratios of prompt light to total
    int0    = hAveWfm_crsmu[ch]->GetCumulative()->GetBinContent(hAveWfm_crsmu[ch]->GetCumulative()->GetXaxis()->FindBin(-10.));
    int100  = hAveWfm_crsmu[ch]->GetCumulative()->GetBinContent(hAveWfm_crsmu[ch]->GetCumulative()->GetXaxis()->FindBin(100.)) - int0;
    int2000 = hAveWfm_crsmu[ch]->GetCumulative()->GetBinContent(hAveWfm_crsmu[ch]->GetCumulative()->GetXaxis()->FindBin(2000.)) - int0;
    int7000 = hAveWfm_crsmu[ch]->GetCumulative()->GetBinContent(hAveWfm_crsmu[ch]->GetCumulative()->GetXaxis()->FindBin(7000.-1.)) - int0;
    std::cout
    <<"  -------------------------------------------\n"
    <<"  Crossing muon waveform:\n"
    <<"   Exp fit, 400ns-2us:\n"
    <<"     - a     = "<<a<<" +/- "<<da<<"\n"
    <<"     - tau   = "<<tau<<" +/- "<<dtau<<" ns\n"
    <<"   Exp + long TPB fit, 400ns-7us:\n"
    <<"     - a     = "<<a2<<" +/- "<<da2<<"\n"
    <<"     - b     = "<<b2<<" +/- "<<db2<<"\n"
    <<"     - tau   = "<<tau2<<" +/- "<<dtau2<<" ns\n"
    <<"   Prompt integral & ratios:\n"
    <<"     - 100ns      = "<<int100<<"\n" 
    <<"     - 100ns/2us  = "<<int100/int2000<<"\n"
    <<"     - 100ns/7us  = "<<int100/int7000<<"\n";
  }
  
  std::cout<<"\n";
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
void MichelAna::PerfectWfmReco(raw::OpDetPulse &pulse) 
{
  size_t ch   = pulse.OpChannel(); 

  fWfmRMS[ch] = 0.;

  if( fMuT0[ch] > 0.  ) fvHitTimes[ch].push_back( fMuT0[ch] );
  if( fT0[ch] > 0.    ) fvHitTimes[ch].push_back( fT0[ch]   );
  fNumOpHits[ch] = fvHitTimes[ch].size();
  
  for(int i=0; i<fNumOpHits[ch]; i++){
        
        float T = fvHitTimes[ch].at(i) - 5;
        if( T < 0 ) T = 0;
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
        float PE_total = Integrate(hPMT_phelTimes[ch],T,T+fFullWindow);

        fvIsHitAtTrigger[ch]          .push_back(true);
//        if( i==1) fvIsHitAtTrigger[ch]    .push_back( true );
//        else      fvIsHitAtTrigger[ch]    .push_back( false );
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
          tmp = PE_total-PE_1800;  
          float pe_total = pe_1800 + tmp + fRand->Gaus(0, CalcSigma( ch,tmp, fFullWindow-1800. ));

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
          fvHitADC_total[ch]      .push_back( pe_total   *fSPE[ch] );
        
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
          fvHitADCpf_total[ch]   .push_back( PE_total * fSPE[ch] );
        } else
        if( i==1 ) {
          fvHitADCpf_100ns[ch]    .push_back( fTrue_PE_prompt[ch] * fSPE[ch] );
          fvHitADCpf_total[ch]   .push_back( fTrue_PE_total[ch] * fSPE[ch]);
        }
          
        fvHitWidth[ch]          .push_back(fMinOpHitWidth[ch]);
        fvHitAmplitudes[ch]     .push_back(0.); // convert PEs to amplitude
        
  }

}



//########################################################################################
// Given a PMT channel, interated photoelectron count, and integration time, calculate the 
// smearing sigma based on noise information provided in the fcl (extracted from data).
float MichelAna::CalcSigma(size_t ch, float PE, float T ) {
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
bool MichelAna::IsPointInFiducialVolume(
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
bool MichelAna::IsPointInFiducialVolume(TVector3 p, float dx = 0., float dy = 0., float dz = 0.)
{
  return IsPointInFiducialVolume(p, dx, dx, dy, dy, dz, dz);
}
bool MichelAna::IsPointInFiducialVolume(TVector3 p, 
  float dx1, float dx2, 
  float dy1, float dy2,
  float dz1, float dz2)
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

//#######################################################################################
// Function for determining if an MCParticle is fully contained in TPC
bool MichelAna::IsParticleContained(const simb::MCParticle& p){
  for( size_t i=0; i<p.NumberTrajectoryPoints(); i++){
    if( !IsPointInFiducialVolume( p.Position(i).Vect(),0.,0.,0.)) 
      return false;
  }
  return true;
}


//########################################################################################
// Retrieve all detector properties (TPC and PMTs)
void MichelAna::GetDetProperties( ){

  if( fRunNumber != fCachedRunNumber ) {
    fCachedRunNumber = fRunNumber;
    
    // ================================================   
    // Get SPE values for PMTs
    if( fIsMC ) {
      fSPE[0]     = 50.;
      fSPE_err[0] = 25.;
      fSPE[1]     = 50.;
      fSPE_err[1] = 25.;
    } else {
      for(size_t i = 0; i < 2; i++){
        for(size_t j=0; j<fCalibRunStart[i].size(); j++){
          if( fRunNumber >= fCalibRunStart[i].at(j) && fRunNumber <= fCalibRunEnd[i].at(j) ){
            fSPE[i]     = fCalibSPE[i].at(j);
            fSPE_err[i] = fCalibSPEerr[i].at(j);
            if( fCalibTau[i].at(j) > 0 ) fMuContamCorr_EffTau[i] = fCalibTau[i].at(j);
            break; 
          }
        }
      }
    }
  
  }// endif new run

}



//########################################################################################
// This function should 
//   - fill all member truth variables
//   - create artificial OpDetPulses (fPulses) to be 
//     passed on to later stages for quasi-waveform reco
void MichelAna::GetTruthInfo(art::Event &e){

  LOG_VERBATIM("MichelAna")<<"GetTruthInfo";
  
  // Get MCParticle information, and exit if none is found
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  e.getByLabel(fSimProducerModule,particleHandle);
  
  // Reset electron/muon ID
  fElectronID = -9;
  fMuonID     = -9;

  LOG_VERBATIM("MichelAna")<<"This event has "<<particleHandle->size()<<" MCParticles";
  if( particleHandle->size() > 1000 ) return;

  // ********************************************************
  // Loop through MC particles
  // *******************************************************
  for( auto const& particle : (*particleHandle) )
  {
    size_t nTrajPts   = particle.NumberTrajectoryPoints();
    int	  last	      = (int)nTrajPts-1;

    if(0){
      TVector3 locStart = particle.Position(0).Vect();
      TVector3 locEnd   = particle.Position(last).Vect();
  //    TVector3 momStart = particle.Momentum(0).Vect();
  //    TVector3 momEnd   = particle.Momentum(last).Vect();
      bool startsInVol  = IsPointInFiducialVolume(locStart,0.,0.,0.);
      bool endsInVol    = IsPointInFiducialVolume(locEnd,0.,0.,0.);
      double dL         = (locStart-locEnd).Mag();
      double P0         = particle.P(0)*1000.;
      double E0         = 1000.*( particle.E(0) - particle.Mass() );
      double Ef         = 1000.*( particle.E(last) - particle.Mass() );
      double dE         = (E0-Ef);
      bool isElectron     = ( particle.TrackId() == fElectronID || (abs(particle.PdgCode())==11 && particle.Mother() == fElectronID) ); ;
      bool isFromElectron = IsParticleDescendedFrom(particleHandle,particle.TrackId(),fElectronID);
      bool isElectronShwr = (isFromElectron || isElectron );
      float displacement = CalcParticleHeritageDisplacement(particleHandle,particle,fElectronID);
   
      // Printout particle list
      printf("  %3i PDG: %5i (%i->%i) Npts=%3lu dL=%7.3f P0=%10.6f Etot=%10.6f E=%10.5f dE=%10.5f  T=%11.4f-%11.4f   moth=%3i   %10s  trk? %i  shwr? %i  displ=%7.3f  Ndaught=%i\n",
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
        isElectron,
        isElectronShwr,
        displacement,
        particle.NumberDaughters()
        );
    }
    
    // =======================================================
    // Look for muon
    if( particle.Process() == "primary" && abs(particle.PdgCode()) == 13 ) {

      // save muon track ID, charge, and stop time, and whether muon is
      // crossing or not
      fMuonID           = particle.TrackId();
      fTrue_MuSign    = -1.*(particle.PdgCode()/abs(particle.PdgCode()));
      fTrue_MuStopTime  = particle.T(last);
      fTrue_MuEnergy    = 1000.*( particle.E(0) - particle.Mass() );
      hTrue_MuEnergy    ->Fill( fTrue_MuEnergy );
      
      // save muon endpoint
      fTrue_MuTrackEnd        = particle.Position(last).Vect();
      fTrue_MuTrackEnd_X      = fTrue_MuTrackEnd.X();
      fTrue_MuTrackEnd_Y      = fTrue_MuTrackEnd.Y();
      fTrue_MuTrackEnd_Z      = fTrue_MuTrackEnd.Z();
      
      // find where muon first enters TPC
      //size_t iVertex = 0;
      bool inTPC = false;
      bool isInFid = false;
      float P0mu = 0.;
      float P1mu = 0.;
      float dLmu = 0.;
      for(size_t i=1; i<nTrajPts; i++) {
        isInFid = IsPointInFiducialVolume( particle.Position(i).Vect(), -0.1, -0.1, -0.1 );
        if( !inTPC && isInFid ) {
          inTPC = true;
          P0mu = particle.P(i);
          //iVertex = i;
          fTrue_MuTrackVertex   = particle.Position(i).Vect();
          fTrue_MuTrackVertex_X = fTrue_MuTrackVertex.X();
          fTrue_MuTrackVertex_Y = fTrue_MuTrackVertex.Y();
          fTrue_MuTrackVertex_Z = fTrue_MuTrackVertex.Z();
        } 
        else if( inTPC && isInFid ) { 
          float dl = (particle.Position(i).Vect()-particle.Position(i-1).Vect() ).Mag();
          dLmu += dl;
        }
        // if the muon track ends (stops) or it exits the TPC,
        else if( (inTPC && !isInFid) || i == nTrajPts-1 ) {
          inTPC = false;
          P1mu = particle.P(i);
          break;
        } 
      }
      if( dLmu > 0 ) fTrue_MuTrackLength = dLmu;
      hTrue_MuAveMomentum->Fill( 0.5*(P0mu+P1mu)*1000. );
      
      LOG_VERBATIM("MichelAna") 
      << "  --> Found primary muon!  Charge = "<<fTrue_MuSign
      <<", ("<<fTrue_MuTrackVertex_X<<","<<fTrue_MuTrackVertex_Y<<","<<fTrue_MuTrackVertex_Z<<") "
      <<" --> ("<<fTrue_MuTrackEnd_X<<","<<fTrue_MuTrackEnd_Y<<","<<fTrue_MuTrackEnd_Z<<")"
      <<"  L = "<<dLmu
      <<"  Ldev = "<< (fTrue_MuTrackVertex - fTrue_MuTrackEnd).Mag() / fTrue_MuTrackLength;
     
      
      // starting from point where muon first enters TPC,
      // record res. range vs. dE/dx info
      if( IsPointInFiducialVolume( fTrue_MuTrackEnd, 0, 0, 0)){
        
        fTrue_IsMuStopping = true;
       
        /* 
        //std::cout<<"Walking backward from mu endpoint...\n";
        // walk backward from endpoint
        float resRange = 0.;
        for(size_t i=particle.NumberTrajectoryPoints()-1; i>iVertex; i--){
          float ke1 = particle.E(i-1) - particle.Mass();
          float ke2 = particle.E(i) - particle.Mass();
          float dE = (ke1-ke2)*1000.; // convert GeV --> MeV
          float dL = (particle.Position(i).Vect()-particle.Position(i-1).Vect()).Mag();
          //std::cout <<"  i="<<i<<"  T= "<<particle.T(i)<<"   p2 = "<<particle.P(i)<<"   ke1= "<<ke1<<"   ke2="<<ke2<<"   dE = "<<dE<<"   dL = "<<dL<<"   dE/dx= "<<dE/dL<<"   resRange = "<<resRange+dL/2.<<"\n";
          if( dL > 0.1 ) hTrue_MuResRangeVsdEdx->Fill( resRange + dL/2., dE/dL );
          resRange += dL;
        }
        */

      } else {
        fTrue_IsMuStopping = false;
      }

        
    }//<-- endif muon


    // =======================================================
    // Look for Michel electron
    if( abs(particle.PdgCode()) == 11 && particle.Mother() == fMuonID ) {
      
      // since the Geant4 "Process" param is useless in case of mu-, 
      // make some clever cuts to ensure we really have a Michel:
      if( -1.*particle.PdgCode()/abs(particle.PdgCode()) == fTrue_MuSign && (particle.T(0) - fTrue_MuStopTime) > 5.0 ){
        
        // save electron track ID, decay time, energy, momentum, and angle
        fElectronID     = particle.TrackId();
        fTrue_dT        = particle.T(0) - fTrue_MuStopTime;
        fTrue_ElEnergy  = (particle.E(0)-particle.Mass())*1000.;
        fTrue_ElMomentum = particle.Momentum(0).Vect();
        fTrue_ElMomentum_X  = fTrue_ElMomentum.X();
        fTrue_ElMomentum_Y  = fTrue_ElMomentum.Y();
        fTrue_ElMomentum_Z  = fTrue_ElMomentum.Z();
        fTrue_ElAngle   = fTrue_ElMomentum.Angle(fTrue_MuTrackEnd - fTrue_MuTrackVertex);
        if( fTrue_ElMomentum_X < 0. ) {
          TVector3 dirToWires(-1., 0., 0.);
          float cos_theta = cos(fTrue_ElMomentum.Angle(dirToWires));
          if( cos_theta > 0. ) fTrue_ProjDist3DToWires = fTrue_MuTrackEnd.X() / cos_theta;
        }
          
        LOG_VERBATIM("MichelAna")
        << "  --> Found Michel electron!  dT = "<<fTrue_dT<<", KE = "<<fTrue_ElEnergy<<", angle (rel. to muon) = "<<fTrue_ElAngle*RAD_TO_DEG<<" deg";

      }//<-- endif Michel electron
    }//<-- endif electron
 

  }//<-- end loop over MCParticles
  
  if( fMuonID >= 0  && fElectronID >= 0 ){ 
    hTrueMuEndpt->Fill( fTrue_MuTrackEnd_X, fTrue_MuTrackEnd_Y, fTrue_MuTrackEnd_Z );
    hTrueMuEndpt_ZX->Fill( fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_X);
  }
  
  // ******************************************************************
  // If the michel electron was found, then simulate the scintillation signals
  // ******************************************************************
  //if( fMuonID > 0 && fElectronID > 0 ) {
   
    // This stage loops through the IDEs from SimChannels and 
    // calculates the scintillation light produced
    ParticleTracker(e, fMuonID, fElectronID, hPMT_phelTimes, hPMT_phelTimes_muon, hPMT_phelTimes_electron); 
    
    fTrue_MuPE = 0;
    fTrue_PE = 0.; 
    fTrue_PEPrompt = 0.;

    // Loop over PMTs
    for( auto & i : fSelectChannels ) {

      float muX1 = 0.;
      float elX1 = 0.;
      // Approximate the "hit time" T0 as the time of arrival of the first photon
      if( fMuonID > 0 ) {
        fMuT0[i] = (float)hPMT_phelTimes[i]->FindFirstBinAbove(0,1);
        muX1 = fMuT0[i]-1;
      }
      if( fElectronID > 0 ) {
        fT0[i]   = (float)hPMT_phelTimes_electron[i]->FindFirstBinAbove(0,1);
        elX1 = fT0[i]-1;
      }

      // Save prompt, total PEs by integrating the photon arrival 
      // time histograms
      if( fElectronID > 0 ) {
        fTrue_PE_prompt[i]      = Integrate(hPMT_phelTimes_electron[i],elX1,fT0[i]+fPromptWindow);
        fTrue_PE_total[i]       = Integrate(hPMT_phelTimes_electron[i],elX1,fT0[i]+fFullWindow);
        fTrue_PE                += fTrue_PE_total[i];
        fTrue_PEPrompt          += fTrue_PE_prompt[i];
        hTrue_PE_preTrig[i]     ->Fill( fTrue_PE_total[i] );
      }
      if( fMuonID > 0 ) {
        fTrue_MuPE_prompt[i]    = Integrate(hPMT_phelTimes_muon[i],     muX1,fMuT0[i]+fPromptWindow);
        fTrue_MuPE_total[i]     = Integrate(hPMT_phelTimes_muon[i],     muX1,fMuT0[i]+fFullWindow);
        fTrue_MuPE              += fTrue_MuPE_total[i];
        if( fElectronID > 0 ) {
          fTrue_MuContam_prompt[i]= Integrate(hPMT_phelTimes_muon[i],     elX1,fT0[i]+fPromptWindow);
          fTrue_MuContam_total[i] = Integrate(hPMT_phelTimes_muon[i],     elX1,fT0[i]+fFullWindow);
        }
      }
      hTrue_LightYield[i]     ->Fill( hPMT_phelTimes[i]->GetEntries() / fTrue_TotalEnergyDep );
    }

    LOG_VERBATIM("MichelAna")
    <<"  Total charge deposited in event: "<<fTrue_TotalChargeDep<<" electrons, total energy = "<<fTrue_TotalEnergyDep<<" MeV\n"
    <<"  Total energy deposited by electron: "<<fTrue_ElShowerEnergyDep<<"\n"
    <<"    #pe on HMM (prompt/total) = "<<fTrue_PE_prompt[0]<<" / "<<fTrue_PE_total[0]<<"\n"
    <<"    #pe on ETL (prompt/total) = "<<fTrue_PE_prompt[1]<<" / "<<fTrue_PE_total[1];
  //}
  
  
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
bool MichelAna::IsParticleDescendedFrom( const art::Handle< std::vector<simb::MCParticle>>& particleHandle, int particleID, int motherID){
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
// Add up total displacement from mother ID
float MichelAna::CalcParticleHeritageDisplacement( const art::Handle< std::vector<simb::MCParticle>>& particleHandle, const simb::MCParticle& particle, int motherID){
  if( particle.TrackId() < 0 || motherID <= 0 ) return -9.;
  
  if( particle.TrackId() == motherID ) 
    return 0.;
  else if (particle.TrackId() < motherID ) 
    return -9.;
  
  int last = particle.NumberTrajectoryPoints()-1;
  float dL = (particle.Position(0).Vect()-particle.Position(last).Vect()).Mag();

  // start at particle and work backwards until we get to the motherID
  int nextMother = particle.Mother();
  bool flag = true;
  while( flag ) {
    
    // if the next particle in the lineage is the motherID, then
    // end the loop and return the cumulative displacement
    if( nextMother == motherID ) break;

    // if the next particle in the lineage is lower than the 
    // motherID, then something went wrong and the input
    // particle isn't actually related to the motherID at
    // all! So return junk dL value
    if( nextMother < motherID ) { dL = -9.; break; };

    // else, find the closest ancestor and update dL
    for( auto const& p : (*particleHandle) ){
      if( p.TrackId() == nextMother ) {
        last = p.NumberTrajectoryPoints()-1;
        dL += (p.Position(0).Vect()-p.Position(last).Vect()).Mag();
        nextMother = p.Mother();
        break;
      }
    }
  }

  return dL;
}



//########################################################################################
// Calculate total distance traversed by a particle
float MichelAna::CalcParticleDist( const simb::MCParticle& p){
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
void MichelAna::ParticleTracker(const art::Event& e, int muID, int elID, TH1D* hPhelTimes[], TH1D* hPhelTimes_muon[], TH1D* hPhelTimes_electron[]){

  LOG_VERBATIM("MichelAna")<<"Particle tracker...";
 
  // Get IDEs
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > simList;
  if(e.getByLabel(fSimProducerModule,SimListHandle)) { 
    art::fill_ptr_vector(simList, SimListHandle); 
  }
  
  // Get MCParticle information
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  e.getByLabel(fSimProducerModule,particleHandle);

  fTrue_NumBremmPhotons = 0;
  fTrue_TotalChargeDep  =0.;
  fTrue_TotalEnergyDep  =0.;
  fTrue_PE0             =0.;
  
  // Keep track of total charge (electrons escaping recomb.) 
  // for Michel electrons and their shower products, and
  // total number of electron-ion pairs. These are used to 
  // calculate recombination factors.
  int   Ni_trk  = 0;
  int   Ni_shwr = 0;
  float Q_trk   = 0.;
  float Q_shwr  = 0.;
  float Qc_trk   = 0.;
  float Qc_shwr  = 0.;

  // for calculating shower centroid
  float wsum_x = 0.;
  float wsum_y = 0.;
  float wsum_z = 0.;
 
  if( muID > 0 ) {
    fTrue_MuCharge      = 0.;
    fTrue_MuPhotons     = 0.;
    fTrue_MuEnergyDep   = 0.;
  }
  
  if( elID > 0 ) {
    fTrue_ElShowerVis     = 0.; // total frac. vis
    fTrue_ElShowerVisCh[0] = 0.; // total frac. vis: ETL
    fTrue_ElShowerVisCh[1] = 0.; // total frac. vis: HMM
    fTrue_ElTrackPhotons = 0;
    fTrue_ElShowerPhotons   = 0;
    fTrue_ElShowerPhotonsPrompt   = 0;
    fTrue_ElShowerPhotonsLate   = 0;
    //fTrue_PE[0] = 0;
    //fTrue_PE[1] = 0;
    fTrue_ElTrackChargeDep  = 0.;
    fTrue_ElTrackEnergyDep  = 0.;
    fTrue_ElShowerChargeDep = 0.;
    fTrue_ElShowerTotalTrajLength = 0.;
    fTrue_ElShowerIons = 0.;
    fTrue_ElShowerEnergyDep = 0.;
  }

  float sumMuEnergyDep_4mm = 0.;
  float sumMuEnergyDep_1cm = 0.;
  float sumMuEnergyDep_2cm = 0.;
  float sumMuEnergyDep_5cm = 0.;
  float sumMuEnergyDep_10cm = 0.;

  TH1D h_tmp("h_tmp","h_tmp",
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmax());
  TH1D h_tmp_n("h_tmp_n","h_tmp_n",
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx->GetXaxis()->GetXmax());
  TH1D h_tmp_HR("h_tmp_HR","h_tmp_HR",
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmax());
  TH1D h_tmp_HRn("h_tmp_HRn","h_tmp_HRn",
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetNbins(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmin(),
    hTrue_MuResRangeVsdEdx_HR->GetXaxis()->GetXmax());

  fTrue_IsElTrkContained = true;
  fTrue_IsElShwrContained = true;

  // Loop through particle list, and for each particle with the common 
  // ancestor, determine energy loss deposited INTO LAr at each step.
  for( auto const& particle : (*particleHandle) )
  {
    
    // Is this particle descended from electron, or is it the muon?
    //   - particle is considered part of the electron "track" if
    //     it is the Michel electron *or* if it's an ionization electron
    //     coming directly off of the Michel electron *or* if its total
    //     displacement from the Michel electron is < 0.4 cm
    //   - all particles descending from the Michel electron are
    //     considered a part of the "shower"
    //bool isElectron     = ( particle.TrackId() == elID || (abs(particle.PdgCode())==11 && particle.Mother() == elID) );
    float displacement  = CalcParticleHeritageDisplacement(particleHandle,particle,elID);
    bool isElectron     = ( particle.TrackId() == elID || (abs(particle.PdgCode())==11 && particle.Mother() == elID) || (displacement > 0 && displacement < 0.4 ));
    bool isFromElectron = IsParticleDescendedFrom(particleHandle,particle.TrackId(),elID);
    bool isElectronShwr = (isFromElectron || isElectron );
    bool isMuon         = (particle.TrackId() == muID);
    bool isFromMuon     = ( IsParticleDescendedFrom(particleHandle,particle.TrackId(),muID) && !isFromElectron );
    
    hTrue_PartDispFromElectron->Fill(displacement);
   
    // Bremm photon plots (not sure if this is right...)
    if( isFromElectron && particle.PdgCode() == 22 ) {
      fTrue_NumBremmPhotons++;
      hTrue_BremmPhotonLength->Fill( (particle.EndPosition().Vect() - particle.Position(0).Vect()).Mag() );
    }
    
    // Skip photons and unrelated particles
    //if( particle.PdgCode() == 22 ) continue;

    // Particle start time
    float particleTime = particle.T(0);

    // Particle's deposited charge and energy
    float particleDepEnergy = 0.;
    float particleDepCharge = 0.;
   
    // Loop over track's IDEs and convert each to scintillation 
    for(size_t nChan = 0; nChan < simList.size(); nChan++) {
      // Restrict to induction plane channels
      if(simList.at(nChan)->Channel() >= 240) break;
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
          hTrue_IDE_X   ->Fill(loc.X());
          hTrue_IDE_Y   ->Fill(loc.Y());
          hTrue_IDE_Z   ->Fill(loc.Z());
          hTrue_IDE_NumElectrons->Fill( it->second.at(i).numElectrons );
          
          // Correct numElectrons for drift attenuation. Account for the
          // transit times between the shield and induction plane.
          // Assuming induction plane is at X = -0.2, shield at X = 0.2
          float T = T_from_X(x,0);
          int numElectrons = it->second.at(i).numElectrons * TMath::Exp( T / fElectronLifetimeFromDB );
          
          // Use deposited energy to back-calculate scintillation.
          int numIon      = dE/fWion;
          int numPhotons  = std::max( (1.+fExcRatio)*numIon - numElectrons, 0.);

          // add up the energy
          fTrue_TotalChargeDep        += numElectrons; 
          fTrue_TotalEnergyDep        += dE;
          particleDepCharge           += numElectrons;
          particleDepEnergy           += dE;
        
//          std::cout<<"ch"<<simList.at(nChan)->Channel()<<"    isMuon? "<<isMuon<<"   dE "<<dE<<"   Q "<<numElectrons<<"\n";
         
          // Keep track of total charge, light, and quanta for the different components
          // of the Michel shower (electron ionization track and displaced shower products)
          if( isMuon ) {
            if( fTrue_IsMuStopping ) {
                float rr = (loc-fTrue_MuTrackEnd).Mag();
                h_tmp.Fill( rr , dE );
                h_tmp_n.Fill(rr);
                h_tmp_HR.Fill( rr , dE );
                h_tmp_HRn.Fill( rr );
            }
          }
          if( isMuon || (!isFromElectron && isFromMuon ) ){
            fTrue_MuPhotons   += numPhotons;
            fTrue_MuCharge    += numElectrons;
            fTrue_MuEnergyDep += dE;
            if( !fTrue_IsMuStopping ) {
              float d = (loc-fTrue_MuTrackVertex).Mag();
              if( d < 0.4 ) sumMuEnergyDep_4mm += dE;
              if( d < 1.0 ) sumMuEnergyDep_1cm += dE;      
              if( d < 2.0 ) sumMuEnergyDep_2cm += dE;      
              if( d < 5.0 ) sumMuEnergyDep_5cm += dE;      
              if( d < 10.0 ) sumMuEnergyDep_10cm += dE;      
            }
          }
          if( isElectron ){
            Ni_trk  += numIon; 
            Q_trk   += numElectrons;
            Qc_trk  += it->second.at(i).numElectrons;
            fTrue_ElTrackPhotons      += numPhotons;
            fTrue_ElTrackEnergyDep    += dE;
          }
          if( isElectronShwr ){
            fTrue_ElShowerEnergyDep   += dE;
            fTrue_ElShowerPhotons     += numPhotons;
            hTrue_ElShowerDepVsDistFromMuTrackEnd->Fill( (loc-fTrue_MuTrackEnd).Mag() );
            wsum_x += loc.X()*numElectrons;
            wsum_y += loc.Y()*numElectrons;
            wsum_z += loc.Z()*numElectrons;
          }
          if( !isElectron && isFromElectron ) {
              Ni_shwr += numIon;
              Q_shwr  += numElectrons;
              Qc_shwr  += it->second.at(i).numElectrons;
          }

          //if( numElectrons == 0 ) {
          //  std::cout<<"numElectrons "<<numElectrons<<"   dE= "<<dE<<"  numIon "<<numIon<<"  XYZ = "<<loc.X()<<", "<<loc.Y()<<", "<<loc.Z()<<"\n";
          //}

          // Divide photons into singlet- and triplet-produced populations
          std::binomial_distribution<int> distPhotonsFast(numPhotons, fSingletToTripletRatio/(1.+fSingletToTripletRatio));
          int numPhotons_fast = distPhotonsFast(generator);
          int numPhotons_late = numPhotons - numPhotons_fast;
          if( isElectronShwr ) {
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
          float numPhotons0 = numPhotons; 
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
            if( isElectronShwr ){ 
              fTrue_PE0 += Ndet[ch];
              // weighted vis
              fTrue_ElShowerVisCh[ch] += total_vis_ch * float(numPhotons0);
              fTrue_ElShowerVis       +=  total_vis_ch * float(numPhotons0); 
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
              PropagatePhoton( particleTime,  tauLAr, tmin[isvis], tmean[isvis], trms[isvis], isFromMuon, isFromElectron, hPhelTimes[ch], hPhelTimes_muon[ch], hPhelTimes_electron[ch]);
              
            }// end loop over photons
          }// end loop over PMTs
        }//<-- endloop over IDEs
      }//<-- endloop over time ticks
    }//<-- loop over channels
    // ----------------------------------------------------------------------------------
      

      // traj. length (if fully contained in TPC)
      float particleTrajLength = particle.Trajectory().TotalLength();
      bool containedInTPC = IsParticleContained(particle);
      
      // Set booleans for electron track/shower containment
      if( isElectron && !containedInTPC )       fTrue_IsElTrkContained  = false;
      if( isElectronShwr && !containedInTPC )   fTrue_IsElShwrContained = false;
      
      if( containedInTPC && particleTrajLength > 0 ) {
       
        float dEdx = particleDepEnergy / particleTrajLength;
        float dQdx = particleDepCharge / particleTrajLength;
        
        // -----------------------------------------
        // Electrons
        if( (isElectron || isFromElectron) ) {
          if( particleTrajLength > 0 ) {
          
            if( particleDepEnergy > 0. ) {
              hTrue_dEdx          ->Fill( dEdx, particleDepEnergy );
              hTrue_dQdx          ->Fill( dQdx, particleDepEnergy );
              hTrue_dEdx_vs_dQdx  ->Fill( dEdx, dQdx, particleDepEnergy );
              if( isElectron ) {
                hTrue_dEdx_ElTrk  ->Fill( dEdx, particleDepEnergy ); 
                hTrue_dQdx_ElTrk  ->Fill( dQdx, particleDepEnergy );
                hTrue_dEdx_vs_dQdx_ElTrk  ->Fill( dEdx, dQdx, particleDepEnergy );
              }
              if( !isElectron && isFromElectron ) {
                hTrue_dEdx_ElShw  ->Fill( dEdx, particleDepEnergy );
                hTrue_dQdx_ElShw  ->Fill( dQdx, particleDepEnergy );
                hTrue_dEdx_vs_dQdx_ElShw  ->Fill( dEdx, dQdx, particleDepEnergy );
              }
            }// endif L > 0.1 && dE > 0
          }// endif > 0
        }// endif electron shower product
      }// traj length cut

  }//<-- end loop over particles

  // Record containment of electron shower
  if( fTrue_ElEnergy > 0 ) 
    fTrue_ContainmentFrac = fTrue_ElShowerEnergyDep / fTrue_ElEnergy;

  // ---------------------------------------------------------
  // Compute stopping muon dE/dx profile
  if( h_tmp.GetEntries() > 0 ) {
    
    float bin_width = 1.;
  
    // first for the less finely-binned case
    bin_width = h_tmp.GetBinWidth(1);
    for(int i=1; i<=h_tmp.GetXaxis()->GetNbins()-1; i++){
//      std::cout<<"  h_tmp bin "<<i<<" counts = "<<h_tmp_n.GetBinContent(i)<<"\n";
      if( h_tmp_n.GetBinContent(i) == 0 ) break;
      float dE   = h_tmp.GetBinContent(i);
      float dEdx = dE / bin_width;
      if( dEdx > 0. && h_tmp_n.GetBinContent(i+1) > 0 && h_tmp_n.GetBinContent(i) >= 5 ) {
        hTrue_MuResRangeVsdEdx->Fill( h_tmp.GetBinCenter(i), dEdx );
      }
    }
    
    // now fill the "high resolution" histogram
    bin_width = h_tmp_HR.GetBinWidth(1);
    for(int i=1; i<=h_tmp_HR.GetXaxis()->GetNbins()-1; i++){
      if( h_tmp_HRn.GetBinContent(i) == 0 ) break;
      //std::cout<<"  h_tmp_HR bin "<<i<<" counts = "<<h_tmp_HRn.GetBinContent(i)<<"\n";
      if( h_tmp_HRn.GetBinContent(i) < 3 ) continue;
      float dEdx = h_tmp_HR.GetBinContent(i) / bin_width;
      if( dEdx > 0. && h_tmp_HRn.GetBinContent(i+1) > 0) hTrue_MuResRangeVsdEdx_HR->Fill( h_tmp_HR.GetBinCenter(i), dEdx );
    }
  
  }
 
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
  fTrue_ElShowerChargeDep   = Q_trk + Q_shwr;
  fTrue_ElShowerIons        = float(Ni_trk + Ni_shwr);
  hTrue_RecombFactor_ElTrk  ->Fill( Q_trk / float(Ni_trk) );
  hTrue_RecombFactor_ElShw  ->Fill( Q_shwr / float(Ni_shwr) );
      
  // Calculate shower centroid
  fTrue_ElShowerCentroid_X = wsum_x / fTrue_ElShowerChargeDep;
  fTrue_ElShowerCentroid_Y = wsum_y / fTrue_ElShowerChargeDep;
  fTrue_ElShowerCentroid_Z = wsum_z / fTrue_ElShowerChargeDep;
        
  // -----------------------------------------
  // Crossing muons
  if(     !fTrue_IsMuStopping 
      &&  fTrue_MuCharge > 0. 
      &&  fTrue_MuTrackLength > 10. 
      &&  fTrue_MuEnergyDep > 0. ) {
    hTrue_dQdx_CrsMu->Fill( fTrue_MuCharge / fTrue_MuTrackLength );
    hTrue_dEdx_CrsMu->Fill( fTrue_MuEnergyDep / fTrue_MuTrackLength );
    if( fabs(fTrue_MuTrackVertex.Y()-20.)<0.01 && sumMuEnergyDep_4mm > 0. ) {
      hTrue_dEdx_CrsMu_4mm->Fill( sumMuEnergyDep_4mm / 0.4 );
      hTrue_dEdx_CrsMu_1cm->Fill( sumMuEnergyDep_1cm / 1. );
      hTrue_dEdx_CrsMu_2cm->Fill( sumMuEnergyDep_2cm / 2. );
      hTrue_dEdx_CrsMu_5cm->Fill( sumMuEnergyDep_5cm / 5. );
      hTrue_dEdx_CrsMu_10cm->Fill(sumMuEnergyDep_10cm/ 10. );
    }
  }

  LOG_VERBATIM("MichelAna")<<"Done with particle tracker.";

}



//########################################################################################
// Fast-propagate a photon to a PMT, given some visibility, propagation time info
void MichelAna::PropagatePhoton(float t0, float tlar, float tmin, float tmean, float trms, bool isMuon, bool isElectron, TH1D* h, TH1D* h_mu, TH1D* h_el){
  
  if( tmean < tmin ) tmin = 0;
  double time = t0;
    
  // add the LAr time component
  time  += fRand->Exp(tlar);
    
  // add physical propagation time
  double dt = 0.;
  while( dt <= tmin ) { 
    dt = fRand->Gaus(tmean,trms);
  }
  time += dt;

  //time += tmin + fRand->Exp(tmean-tmin);
      
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
  if(isMuon)     h_mu ->Fill(time);

}



//########################################################################################
// Given a vector of hits and some reference 'X' coordinate, find the hit that is best
// matched to this reference. 
void MichelAna::FindBestMatchedHit( int plane, float x, int &bestMatchedHit, float &minDist){
  std::vector<int> hits;
  FindBestMatchedHit(plane, x, bestMatchedHit, minDist, hits);
}
void MichelAna::FindBestMatchedHit( int plane, float x, int &bestMatchedHit, float &minDist, std::vector<int>& hits){
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
// The primary all-in-one clustering routine which takes in some empty MichelCluster data
// object and performs a proximity clustering procedure. A 2D shower is reconstructed if the
// mu-electron boundary is found. This function requires that the candidate stopping mu track
// was already identified and the class member variables fMuTrackVertex,fMuTrackHits, and 
// fMuTrackIndex are defined.

void MichelAna::Clustering( int plane, MichelCluster& cl) { 
  if( fMuTrackIndex < 0 ) return;
  Clustering( plane, cl, float(fMuTrackVertex.X()), fMuTrackHits, -99.);
}

void MichelAna::Clustering( int plane, MichelCluster& cl, float xbnd ) { 
  if( fMuTrackIndex < 0 ) return;
  Clustering( plane, cl, float(fMuTrackVertex.X()), fMuTrackHits, xbnd );
}

void MichelAna::Clustering( int plane, MichelCluster& cl, float muTrackX, std::vector<int> &muTrackHits, float xbnd ) { 

  // Must have already done 3D tracking and tagged a candidate muon 
  // track for determination of seed, and also calculated X/W positions for all hits.
  if( fHitX.size() == 0 ) return;

  // Assign plane
  cl.plane = plane;
  LOG_VERBATIM("MichelAna")
  <<"Begin clustering on plane "<<plane<<"\n"
  <<"Looking for hit best matched to mu track startX = "<<muTrackX<<" cm .... "<<muTrackHits.size()<<" total hits associated with mu trk";
 
  // -------------------------------------------------------------
  // Find hit that is best matched to the muon's vertex X coordinate
  int muStartHit = -9;
  float minDistStart = 9999.;
  FindBestMatchedHit( plane, muTrackX, muStartHit, minDistStart, muTrackHits );
  LOG_VERBATIM("MichelAna")
  <<"Mu start hit dist = "<<minDistStart<<", key = "<<muStartHit;
  if( muStartHit < 0 ) return; 
  
  // --------------------------------------------------------------- 
  // Perform clustering of the muon + electron hits, using the hit
  // best matched to the reco vertex of our stopping track as the seed.
  std::vector<int> cluster;
  int seedHit = muStartHit;
  if( minDistStart  < 5.0 ) ProximityCluster(fHitX, fHitW, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  
  // Starting again at the seedpoint, try clustering again (excluding hits 
  // previously clustered) to see if the initial clustering was incomplete.
  std::vector<int> cluster2 = cluster;
  ProximityCluster(fHitX, fHitW, fHitPlane, cluster2, plane, seedHit, fMaxHitSeparation);

  // If cluster2 is same size as cluster, then no new hits were added in the
  // reclustering and we should just use the original cluster.
  // If cluster2 is larger, it means the original cluster was incomplete and 
  // we need to redefine the cluster using new endpoint
  if( cluster2.size() > cluster.size() ) {
    seedHit = cluster2.at(cluster2.size()-1);
    cluster.clear(); 
    ProximityCluster(fHitX, fHitW, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  } 
 
  // One final check to ensure the endpoint closest to the mu start X is
  // the one used as the seed    
  float x_end   = fHitX.at(cluster[cluster.size()-1]);
  float x_start = fHitX.at(cluster[0]);
  if( fabs(x_end - fMuTrackVertex.X()) < fabs(x_start - fMuTrackVertex.X()) ) {
    seedHit = cluster.at(cluster.size()-1);
    cluster.clear(); 
    ProximityCluster(fHitX, fHitW, fHitPlane, cluster, plane, seedHit, fMaxHitSeparation);
  }

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
      //profile_dQ.at(i-1)    = fHitIntegral[cluster.at(i)];
      profile_X.at(i-1)     = fHitX[cluster.at(i)];
      profile_W.at(i-1)     = fHitW[cluster.at(i)];
      //profile_dQds.at(i-1)  = fHitCharge[cluster.at(i)] / ds;
      profile_dQds.at(i-1)  = fHitCharge[cluster.at(i)];
      //profile_dQds.at(i-1)  = fHitIntegral[cluster.at(i)] / ds;
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
    float dropmax = -9.;
    size_t k1 = 0;
    size_t k2 = (nprof-1);
    for(size_t i=k1; i<k2; i++){
      if( profile_dQds_t.at(i) > dQtmax ) {
        dQtmax = profile_dQds_t.at(i);
        bnd1 = i;
      }
    }

    // reject if the max is the first sample
    if( (bnd1 >= 0 && bnd1 < 2) ) bnd1 = -99; 

    LOG_VERBATIM("MichelAna") 
    <<"Determining boundary point in dQds, dQds_t_max = "<<dQtmax<<" at index "<<bnd1;

    // ---------------------------------------------------------------------
    // Calculate slope of truncated mean preceeding this max
    if( bnd1 > 0 ) {
      TF1 fLine("fLine","[0]+[1]*x");
      TGraphErrors g;
      k1 = std::max(bnd1-15,0);
      k2 = bnd1;
      for(size_t i=k1; i<=k2; i++){
        //g.SetPoint(g.GetN(),  i, profile_dQds_t.at(i));
        //g.SetPointError(g.GetN()-1,0.,0.05*profile_dQds_t.at(i));
        g.SetPoint(g.GetN(), profile_s.at(i), profile_dQds_t.at(i));
        g.SetPointError(g.GetN()-1,0.,0.05*profile_dQds_t.at(i));
      }
      fLine.SetParameter(0,1);
      fLine.SetParameter(1,10);
      g.Fit("fLine","QN0");
      cl.truncMeanSlope = fLine.GetParameter(1); 
      LOG_VERBATIM("MichelAna") 
      <<"Slope preceding max: "<<cl.truncMeanSlope<<"  chi2 = "<<fLine.GetChisquare()<<"   NDF: "<<fLine.GetNDF();
      hTruncMeanSlope->Fill(cl.truncMeanSlope);
      hTruncMeanSlopeChi2->Fill(fLine.GetChisquare()/(fLine.GetNDF()-1));
    
    
      // ----------------------------------------------------------------------
      // The truncated mean peak will precede the charge drop-off boundary, 
      // so only look forward
      int i_maxQ = -99;
      int i_maxDrop = -99;
      for(size_t i=0; i< nprof; i++){
        int dist = int(i) - bnd1;
        float drop = -9.;
        if( i < nprof-1 ) drop = (profile_dQds.at(i) - profile_dQds.at(i+1));
        //if( ( dist >=0 ) && ( dist <= fMaxMatchDistance ) ) {
        if( abs(dist) <= fMaxMatchDistance ) {
          if( profile_dQds.at(i) > dQmax ) {
            dQmax = profile_dQds.at(i);
            i_maxQ = i;
          }
          if( drop > dropmax ) {
            dropmax = drop;
            i_maxDrop = i;
          }
        }
      }
        
      hDistMaxTdQdsToMaxQ->Fill(i_maxQ - bnd1);
      hDistMaxTdQdsToMaxDrop->Fill(i_maxDrop - bnd1);
    
      if( i_maxQ >= 0 ) bnd2 = i_maxQ;
      if( fUseMaxDropAsBnd && i_maxDrop >= 0 ) bnd2 = i_maxDrop;
     
      // Make sure there are no other hits of comparable size along the electron trajectory 
      // (if so, it means we were tricked and this might not be the true muon end)
      if( bnd2 >= 0 ) {
        float dQbnd = profile_dQds.at(bnd2);
        for(size_t i=bnd2+1; i< nprof; i++){
          float ratio = profile_dQds.at(i) / dQbnd;
          hElHitChargeRatio->Fill( i-bnd2, ratio );
          //if( i-bnd2 <= 1 && ratio > 0.5 ) {
          //  // re-assign boundary point
          //  dQmax = profile_dQds.at(i);
          //  bnd2 = i;
          //}
        }
      }
    
    }
    LOG_VERBATIM("MichelAna") 
    <<"dQds_max = "<<dQmax<<" at index "<<bnd2;


    if( bnd2 >= 0 ) {
      cl.maxQ_i = bnd2;
      cl.maxQ_s = profile_s.at(bnd2);
    

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
      
      LOG_VERBATIM("MichelAna")
      <<"Using linearity threshold: "<<linearityThresh;

      // -------------------------------------
      // Find the covariance at the boundary
      float minCov = 9;
      //if( bnd2 >= 0 ) minCov = profile_linearity.at(bnd2);
      for(int i = std::max(int(0), bnd2-1); i<= std::min(int(nprof-1), bnd2+1); i++){
        if( profile_linearity.at(i) < minCov ) minCov = profile_linearity.at(i);
      }
      cl.minCovAtBnd = minCov;
    
      // ------------------------------------------------------
      // If an allowable 'x' coordinate for the boundary was provided
      // as an argument (xbnd), then merely require that this boundary
      // point be close to it. Otherwise, look for kink using linearity.
      if( xbnd <= 0 ) {
        LOG_VERBATIM("MichelAna")<<"Lin threshold = "<<linearityThresh<<", minCov = "<<cl.minCovAtBnd;
        if(  minCov > linearityThresh ) bnd2 = -9;
      } else {
        LOG_VERBATIM("MichelAna")<<"Allowable X-coordinate provided ("<<xbnd<<"), ours is "<<fHitX.at(profile_hitKey[bnd2]);
        double dx = fSamplingRate*fHitRMS.at( profile_hitKey[bnd2] ) * fDriftVelocity[0];
        if( fabs(fHitX.at(profile_hitKey[bnd2]) - xbnd ) > dx ) bnd2 = -9; 
      }
    }


    if( bnd2 >= 0 ) {
      LOG_VERBATIM("MichelAna")<<"Boundary pt identified!  bnd2 = "<<bnd2;
      if( fBndOffset != 0 ){
        LOG_VERBATIM("MichelAna")<<"Offsetting bnd by "<<fBndOffset<<" hits";
        bnd2 += fBndOffset;
      }
      cl.bnd_i = bnd2;
      cl.bnd_s = profile_s.at(bnd2);
      cl.muEnd2D.SetXYZ(profile_W.at(bnd2), profile_X.at(bnd2), 0.);
      cl.muEndHit = profile_hitKey.at(bnd2); 
      cl.muEnd2D_X = profile_X.at(bnd2);
      cl.muEnd2D_W = profile_W.at(bnd2);
    }
   
    // Assign profiles to MichelCluster object
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
      cl.muCharge = 0.;
      float corr = 1.;
      if( fDecayTime > 0 ) corr = exp((fDecayTime/1000.)/fElectronLifetime);
      for(size_t i=0; i< nprof ; i++){
        if( (int)i <= bnd2 ) {
          cl.cluster_mu.push_back(profile_hitKey.at(i));
          cl.muCharge += fHitCharge.at(profile_hitKey.at(i)) * corr;
          if( profile_linearity.at(i) > fLinThreshMuFit && bnd2-(int)i < 20 ) {
            if( muEndDir_Pt1.X() == -9. ) muEndDir_Pt1.SetXYZ(profile_W.at(i), profile_X.at(i), 0.);
            muEndDir_Pt2.SetXYZ(profile_W.at(i), profile_X.at(i), 0.);
            g_mu.SetPoint(g_mu.GetN(),profile_W.at(i), profile_X.at(i));
          }
        } else { 
          //if( (int)i-bnd2 <= 1 || profile_ds.at(i) < 1.0 ){
            cl.cluster_el.push_back(profile_hitKey.at(i));
            g_el.SetPoint(g_el.GetN(),profile_W.at(i), profile_X.at(i));
          //}
        }
      }
     
      int muHits = cl.cluster_mu.size();
      int elHits = cl.cluster_el.size();
      cl.nPtsMuFit = g_mu.GetN();
      hMuClusterHitsEndFit->Fill(cl.nPtsMuFit);
     
      // If we have > 1 electron track hits in the cluster, define the electron
      // start-point as the first hit after the boundary point 
      if( elHits > 1 ) {
        elStart2D.SetXYZ( fHitW.at(cl.cluster_el[0]), fHitX.at(cl.cluster_el[0]), 0.);
      } else {
        elStart2D = cl.muEnd2D;
      }

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
      
      LOG_VERBATIM("MichelAna")
      <<"Linearity check... good mu hits = "<<numMuHitsGood<<" (frac = "<<cl.fracMuHitsLinear<<")";
    
      // --------------------------------------------------------------------
      // Count extra unclustered hits around muon endpoint and electron startpoint
      /*
      cl.extraHits   = 0;
      for(size_t i=0; i<fHitX.size(); i++){
        if( fHitPlane.at(i) != plane ) continue; 
        if( std::find(cluster.begin(), cluster.end(), i ) != cluster.end() ) continue;
        TVector3 loc( fHitW.at(i), fHitX.at(i), 0. );
        float distMu = (loc - cl.muEnd2D).Mag();
        float distEl = (loc - elStart2D).Mag();
        if( distMu <= 2.0 || distEl <= 2.0 ) cl.extraHits++;
      }
      */
      
      // ----------------------------------------------------------------
      // Count up extra unclustered hits along entire cluster
      cl.extraHits   = 0;
      for(size_t i=0; i<fHitX.size(); i++){
        if( fHitPlane.at(i) != plane || std::find(cluster.begin(), cluster.end(), i ) != cluster.end() ) continue;
        //TVector3 loc( fHitW.at(i), fHitX.at(i), 0. );
        for(size_t j=0; j<cluster.size(); j++){
          //TVector3 clstrHitLoc( fHitW.at(cluster.at(j)), fHitX.at(cluster.at(j)), 0.);
          float d = sqrt( pow(fHitW.at(cluster.at(j))-fHitW.at(i),2) + pow(fHitX.at(cluster.at(j))-fHitX.at(i),2)  );
          if( d <= 2.0 ) {
            cl.extraHits++;
            break;
          }
        }
      }
      LOG_VERBATIM("MichelAna")
      <<"Number of points to be used in muon fit: "<<cl.nPtsMuFit<<" (min threshold: "<<fMinMuClusterHitsEndFit<<")";

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
        if( elHits > 0 ) {
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
          elDir2D.SetMag(1.);
          cl.decayAngle2D = cl.muEndDir2D.Angle(elDir2D) * RAD_TO_DEG;
          cl.elDir2D_W = elDir2D.X();
          cl.elDir2D_X = elDir2D.Y();
          LOG_VERBATIM("MichelAna")
          <<"Angle between muon and electron clusters: "<<cl.decayAngle2D;
        }

      }//<-- endif enough muon points for directional fit

      /* 
      // Before we start the showering, bifurcate the electron cluster 
      // to better separate track and shower-like components by cutting
      // the cluster where the separation exceeds 2 cm.
      std::cout<<"Trimming the electron cluster... bnd2= "<<bnd2<<"   nprof= "<<nprof<<"\n";
      std::vector<int> tmp_vec;
      for(size_t i=bnd2+1; i< nprof; i++){
        std::cout<<"  "<<i<<"   ds= "<<profile_ds.at(i)<<"\n";   
        if( int(i)-bnd2 < 2 || profile_ds.at(i) < 2.0 ) {
          tmp_vec.push_back( profile_hitKey.at(i) );
        } else {
          std::cout<<"*** breaking off here ***\n";
          break;
        }
      }
      cl.cluster_el = tmp_vec;
      */

    }//<-- endif a boundary point was determined
  
  }//<-- endif cluster size cut
   
   
  

  // =====================================================================
  // Electron shower clustering: if a decay angle was calculated above, then
  // proceed to cluster all hits along the reconstructed direction in 2D.
  
  if( cl.decayAngle2D > 0 ){
    float sum_t = 0.; 
    float sum_x = 0.;
    float sum_w = 0.;
    
    // ---------------------------------------------------------------
    // Extend the collection of hits associated with the muon cluster
    // by adding extra unclustered hits within proximity
    std::vector<int> muHits = cl.cluster_mu;
    /*
    std::vector<int> extraHits;
    for(int i=0; i<(int)fHitX.size(); i++){
      int hitPlane = fHitPlane.at(i);
      if( hitPlane != plane ) continue; 
      // if already one of the clustered hits, either muon or electron, skip
      if( std::find( cl.cluster.begin(), cl.cluster.end(), i ) != cl.cluster.end() ) continue;
      // loop through muon hits
      for(int j=0; j<(int)cl.cluster_mu.size(); j++){
        int hitID = cl.cluster_mu.at(j);
        float d = sqrt(pow(fHitX.at(hitID)-fHitX.at(i),2) + pow(fHitW.at(hitID)-fHitW.at(i),2));
          if( d <= 1. ) extraHits.push_back(i); 
      }
    }
    //std::cout<<"Expanding the muon hit cluster. Number hits before: "<<cl.cluster_mu.size()<<" --> extra hits added: "<<extraHits.size()<<"\n";
    for(size_t i=0; i<extraHits.size(); i++){
      muHits.push_back(extraHits.at(i) );
    }
    */
  
//    std::cout<<"Forming electron shower on plane "<<plane<<"\n";
    for(size_t i=0; i<fHitX.size(); i++){
       
        if( fHitPlane.at(i) != plane ) continue;
        //if( fHitAmplitude.at(i) < fMinHitAmp[plane] ) continue; 
        if( fHitCharge.at(i) < fMinHitCharge[plane] ) continue; 
        
        bool isInElCluster = false;
        bool isMuonHit     = false;
        bool isBgHit       = false;
        if( std::find(cl.cluster_el.begin(), cl.cluster_el.end(), i ) != cl.cluster_el.end() ) isInElCluster = true;
        if( std::find(muHits.begin(), muHits.end(), i ) != muHits.end() ) isMuonHit = true;
        //if( std::find(cl.cluster_mu.begin(), cl.cluster_mu.end(), i ) != cl.cluster_mu.end() ) isMuonHit = true;
//        if( std::find(fBgTrackHits.begin(), fBgTrackHits.end(), i ) != fBgTrackHits.end() ) isBgHit = true;
       
        if( isMuonHit || (!isInElCluster && isBgHit) ) continue;

        // Find angle of this hit relative to Michel direction (2D)
        TVector3 loc(fHitW.at(i), fHitX.at(i), 0.);
        TVector3 hv = loc - elStart2D; 
        float dist = (loc-elStart2D).Mag();
        float ang = hv.Angle(elDir2D) * RAD_TO_DEG;
        
        if( cl.plane == 1 ) {
          hElShowerDepDist    ->Fill(dist);
          if( dist > 1.5 ) hElShowerDepAngle2D ->Fill(ang);
        }
   
    
        // Add this hit to the shower if it:
        //  (a) is part of the original electron cluster
        //  (b) lies within proximity to the boundary
        //  (c) is within 2D acceptance angle of electron direction
        if( isInElCluster || ( dist < 1.5 ) || (!isMuonHit && ang <= fShowerAcceptanceAngle) ) {
        
          //std::cout<<"   "<<i<<"  X: "<<fHitX.at(i)<<"  charge: "<<fHitCharge.at(i)<<"   amp: "<<fHitAmplitude.at(i)<<"  RMS: "<<fHitRMS.at(i)<<"   isInTrk? "<<isInElCluster<<"\n";
          cl.shower     .push_back( i );
          cl.isInTrk    .push_back(isInElCluster);
          sum_t         += fHitT.at(i) * fHitCharge.at(i);
          sum_x         += fHitX.at(i) * fHitCharge.at(i);
          sum_w         += fHitCharge.at(i);
          if( !isInElCluster ) {
            hHitRMS_ElShower[plane]->Fill( fHitRMS.at(i) );
            hHitAmp_ElShower[plane]->Fill( fHitAmplitude.at(i) );
            hHitIntegral_ElShower[plane]->Fill( fHitIntegral.at(i) );
            //if( fHitPlane.at(i) == 1 ) {
            //  hHitIntegral_vs_DiffSummedADC_ElShower->Fill( fHitIntegral.at(i), (fHitIntegral.at(i)-fHitSummedADC.at(i))/fHitSummedADC.at(i) );
            //  hHitIntegral_vs_SummedADC_ElShower->Fill( fHitIntegral.at(i), fHitSummedADC.at(i) );
            //  hHitAmp_vs_DiffSummedADC_ElShower     ->Fill( fHitAmplitude.at(i), (fHitIntegral.at(i)-fHitSummedADC.at(i))/fHitSummedADC.at(i) );
            //}
          } //else {
            //if( fHitPlane.at(i) == 1 ) {
            //  hHitIntegral_vs_DiffSummedADC_ElTrk->Fill( fHitIntegral.at(i), (fHitIntegral.at(i)-fHitSummedADC.at(i))/fHitSummedADC.at(i) );
            //  hHitAmp_vs_DiffSummedADC_ElTrk     ->Fill( fHitAmplitude.at(i), (fHitIntegral.at(i)-fHitSummedADC.at(i))/fHitSummedADC.at(i) );
            //}
          //}
        }//<-- end if in 2D cone
        
    }//<-- end loop over hits

    if( sum_w > 0. ) {
      cl.aveDriftTime  = sum_t / sum_w;
      cl.aveX          = sum_x / sum_w;
    }
  
    cl.elShowerFrac = cl.shower.size() / float( fNumPlaneHits[plane] - (int)cl.cluster_mu.size());
  
  } // end Michel shower 2D clustering (if decay angle found)

}



//########################################################################################
void MichelAna::CalcMichelShowerEnergy( MichelCluster& cl ) { 

  cl.elShowerCharge     = 0.;
  cl.elShowerEnergy     = 0.;
  cl.elTrackEnergy      = 0.;
  cl.elTrackCharge      = 0.;

  for(size_t i=0; i<cl.shower.size(); i++){
    
    size_t  hitID       = cl.shower.at(i);

    if( fHitX.at(hitID) < fEdgeMarginX ) continue;

    bool    isInTrack   = cl.isInTrk.at(i);
    float   q           = fHitCharge.at(hitID);
    float   rfac        = fRecombFactorShower;
    if(isInTrack) rfac  = fRecombFactorTrack;

    // Add to Michel charge
    cl.elShowerCharge         += q;
    if( isInTrack )
      cl.elTrackCharge        += q;
    
    hHitCharge_ElShower[cl.plane]->Fill(q);
     
    // Scale charge to energy
    //   dE = N_i * W_ion
    //   N_i = dQ/R
    //   W_ion = W_ph * (1+alpha)
    float dE = ( q / rfac ) * fWion;
    cl.elShowerEnergy   += dE;
    if( isInTrack ) 
      cl.elTrackEnergy  += dE;

  }//<-- end loop over hits
 
  LOG_VERBATIM("MichelAna")<<"Plane "<<cl.plane
  <<": Michel electron charge= "<<cl.elShowerCharge
  <<", trk energy = "<<cl.elTrackEnergy
  <<", shwr energy= "<<cl.elShowerEnergy;

}


//########################################################################################
// A simple clustering algorithm that starts at some "seed" hit and clusters hits based on
// proximity in 2D 'WX' space.
void MichelAna::ProximityCluster(
std::vector<float>& hitX, std::vector<float>& hitW, std::vector<int>& hitPlane, 
std::vector<int>& cluster, int plane, int seedHit, float distThresh)
{
  LOG_VERBATIM("MichelAna")
  <<"Beginning proximity-based clustering; input cluster has "<<cluster.size()<<" hits to start with.\n"
  <<"  seed hit = "<<hitW[seedHit]<<","<<hitX[seedHit]<<"\n"
  <<"  looking at plane "<<plane;
  
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
 
  LOG_VERBATIM("MichelAna")
  <<"  total hits in plane "<<plane<<": "<<fNumPlaneHits[plane]<<", cluster size: "<<cluster.size();
}


/*
//########################################################################################
void MichelAna::SuperSonicCluster(
std::vector<float>& hitX, std::vector<float>& hitW, std::vector<int>& hitPlane, 
std::vector<int>& cluster, const std::vector<int>& clusterExl, int plane, int seedHit, float distThresh)
{
  LOG_VERBATIM("MichelAna")
  <<"Beginning supersonic clustering; input cluster has "<<cluster.size()<<" hits to start with.\n"
  <<"  seed hit = "<<hitW[seedHit]<<","<<hitX[seedHit]<<"\n"
  <<"  looking at plane "<<plane;
 
  // Check if the seed hit is already in cluster; if not, add it
  if( std::find(cluster.begin(), cluster.end(), seedHit ) == cluster.end() ) 
    cluster.push_back(seedHit);

  // Starting from the seed hit, look for nearest hit not already in the cluster
  float radius = 0.;
  for(size_t i=0; i<hitX.size(); i++){
    // Skip hits not in plane, or already in cluster
    if( hitPlane.at(i) != plane ) continue;
    if( std::find(cluster.begin(), cluster.end(), i ) != cluster.end() ) continue;
    if( std::find(clusterExl.begin(), clusterExl.end(), i ) != clusterExl.end() ) continue;
    TVector3 hit_loc(hitX[i], hitW[i], 0.);
    TVector3 hit_loc_prev(hitX[prevHit], hitW[prevHit], 0.);
    TVector3 dir = hit_loc - hit_loc_prev;
    float dist = dir.Mag();
    if( dist < minDist ) { minDist = dist; minKey = i; }
  }
  
  // UNDER CONSTRUCTION  
 
}
*/



//########################################################################################
// Calculate the truncated mean in the neighborhood of a sample within a clustered hit profile.
float MichelAna::CalcTruncatedMeanInProfile(std::vector<float>& v, size_t index, int nb, float p){
  
  //???
  // Resize neighborhood depending on how close to edge we are
  //nb = std::min( nb, (int)index );
  //nb = std::min( nb, (int)(v.size()-index-1));
  //???
  
  // Restrict neighborhood (nb) at edges
  //if( nb > 1 && (index < 2 || index >= v.size() - 2 )) nb = 1;

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
float MichelAna::CalcLocalLinearity(std::vector<float>& vx, std::vector<float>& vy, size_t index, int nb){

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

  // ................................
  aveX = vx.at(index);
  aveY = vy.at(index);
  // ................................

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
    return fabs(cov) / (stDevX * stDevY);
  else 
    return 1.;
}



//########################################################################################
// This function performs a fit to the (binned) muon late-light and extrapolates this fit
// forward to estimate how much of this light spills over into the Michel electron integration
// window. It then subtracts this off of the integrated PE for the Michel.
void MichelAna::MuContamCorr(int ch) {
  
  LOG_VERBATIM("MichelAna")
  <<"Correcting for muon late-light... ch "<<ch<<"  dT = "<<fdT[ch];
  // Check that there are two hits
  if( fvHitTimes[ch].size() != 2 || fSPE[ch] <= 0 || fPE_totalRaw[ch] < 0) return;

//  const int N = 7;
//  float binW[N] = { 100., 200., 200., 200., 300.,   300., 300.}; 
//  float x[N]    = { 250., 400., 600., 800., 1050., 1350., 1650.};
//  400-500
//  500-700
//  700-900
//  900-1200
//  1200-1500
//  1500-1800
  const int N = 6;
  float binW[N] = { 100., 200., 200., 300.,   300., 300.}; 
  float x[N]    = { 450., 600., 800., 1050., 1350., 1650.};
  float y[N];
  float dx[N], dy[N];
  float sum_100 = 0.;
  float sum_total = 0.;
  float sum_total_mu = -99.;
 
  // there must be at least 2 bins
  if( fdT[ch] >= fMuContamCorr_dTmin && fdT[ch] >= x[1]+binW[1]/2. ) {
    
    sum_total_mu = 0.;
  
    //float adc100 = fvHitADC_100ns[ch].at(0);
    //float adc200 = fvHitADC_200ns[ch].at(0);
    //float adc300 = fvHitADC_300ns[ch].at(0);
    float adc400 = fvHitADC_400ns[ch].at(0); 
    float adc500 = fvHitADC_500ns[ch].at(0);
    float adc700 = fvHitADC_700ns[ch].at(0);
    float adc900 = fvHitADC_900ns[ch].at(0);
    float adc1200 = fvHitADC_1200ns[ch].at(0);
    float adc1500 = fvHitADC_1500ns[ch].at(0);
    float adc1800 = fvHitADC_1800ns[ch].at(0);

    //y[0] = ((adc300-adc200) / fSPE[ch])/binW[0];
    y[0] = ((adc500-adc400) / fSPE[ch])/binW[0];
    y[1] = ((adc700-adc500) / fSPE[ch])/binW[1];
    y[2] = ((adc900-adc700) / fSPE[ch])/binW[2];
    y[3] = ((adc1200-adc900) / fSPE[ch])/binW[3];
    y[4] = ((adc1500-adc1200) / fSPE[ch])/binW[4];
    y[5] = ((adc1800-adc1500) / fSPE[ch])/binW[5];
  
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
    
    TGraphErrors gr;

    for(int i=0; i<N; i++){
      if( fdT[ch] < x[i] - binW[i]/2. ) continue;
      dx[i] = binW[i]/2.;
      if( y[i] >= 0 ) dy[i] = sqrt(y[i]*binW[i]) / binW[i];
      else            dy[i] = 0.;
      gr.SetPoint(gr.GetN(), x[i], y[i] );
      gr.SetPointError(gr.GetN()-1, dx[i], dy[i]);
      if( fdT[ch] > x[N-1] + binW[N-1]/2. ) {
        hAvePhelProfile_mu[ch]->Fill(x[i],y[i]);
        avePhelProfileCounter_mu[ch]++;
      }
    }
 
    /* 
    TF1 fit("fit","[0]*exp(-x/[1])",100.,1800.);
    fit.SetParameter(0, y[0] );  
    fit.SetParameter(1, 1300. );
    fit.SetParLimits(1, 0., 9999. );
    if( fMuContamCorr_EffTau[ch] > 0 ) {// SetTau 
      fit.SetParameter(1, TauConversion(fMuContamCorr_EffTau[ch]) );
//      fit.FixParameter(1, TauConversion(fMuContamCorr_EffTau[ch]) );
    }
    */
   
//    float tau = TauConversion(fMuContamCorr_EffTau[ch]);
    float tau = fMuContamCorr_EffTau[ch];

    // The TPB component is ~11.5% of late-light component + prompt light
    // according to MC studies performed on 5/28/2018.
    //   let PE_prompt = param [2]
    //   norm*[0]*3550. = 0.115 * ( [2] + [0]*[1] )
    //   --> norm = (0.115*([2]+[0]*[1])/([0]*3550.))
    TF1 fit("fit","[0]*exp(-x/[1])+0.115*(([2]+[0]*[1])/3550.)*exp(-x/3550)",100.,1800.);
    fit.SetParameter(0, y[0]);
    fit.FixParameter(1, tau);
    fit.FixParameter(2, 0.);


    gr.Fit("fit","RN0Q");

    tau = fit.GetParameter(1);
    float Chi2 = fit.GetChisquare() / fit.GetNDF(); 
    LOG_VERBATIM("MichelAna")<<"  tau = "<<tau<<", red chi2 = "<<Chi2;
    hMuTau[ch]->Fill(tau);
  
    // Save info to be plotted later with waveforms
    fvPrepulseX1[ch].at(1)        = fvHitTimes[ch].at(0) + 400;
    fvPrepulseSlowNorm[ch].at(1)  = fit.GetParameter(0)*fSPE[ch]*-1.;
    fvPrepulseZeroPoint[ch].at(1) = fvHitTimes[ch].at(0);
    fvPrepulseSlowTau[ch].at(1)   = fit.GetParameter(1);
  
    // Calculate integral of this component
    // subtract 8%
    /*
    float A = fit.GetParameter(0);
    fit.SetParameter(0, A*0.92);
    // Add 8% TPB component
    float tau2  = 3550.;
    float A2    = (0.08*A*tau)/tau2;
    TF1 fit2("fit2","[0]*exp(-x/[1])",100.,1800.);
    fit2.SetParameter(0, A2);
    fit2.SetParameter(1, tau2);
    */

    int k1 = fdT[ch]-10;
    int k2a = fdT[ch] + fPromptWindow;
    int k2 = fdT[ch] + fFullWindow;
    for(int i=400; i<k2; i++){
      float val = fit.Eval(i);// + fit2.Eval(i);
      if( i < fFullWindow ) sum_total_mu  += val;
      if( i >= k1 ) {
        sum_total += val;
        if( i < k2a ) sum_100   += val;
      }
    }
  
  } else {
    
    sum_100   = fPE_promptRaw[ch]  - fvHitADCpf_100ns[ch].at(1)/fSPE[ch];
    sum_total = fPE_totalRaw[ch]   - fvHitADCpf_total[ch].at(1)/fSPE[ch];
  
  }
  

  //TGraphErrors gr(N,x,y,dx,dy);

  fMuContam_prompt[ch]  = sum_100;
  fMuContam_total[ch]   = sum_total;

  // Assign total muon light estimated from fit
  if( sum_total_mu > 0. ) fMuPE_total[ch] = fvHitADC_400ns[ch].at(0)/fSPE[ch] + sum_total_mu; 

  // Correct off the contamination
  fPE_prompt[ch] = fPE_promptRaw[ch] - sum_100;
  fPE_total[ch]  = fPE_totalRaw[ch] - sum_total;
  
  LOG_VERBATIM("MichelAna")
  <<"  pe100 contam: "<<sum_100<<", total contam: "<<sum_total<<"\n"
  <<"  raw prompt: "<<fPE_promptRaw[ch]<<", total: "<<fPE_totalRaw[ch]<<"\n"
  <<"  corrected prompt: "<<fPE_prompt[ch]<<", total: "<<fPE_total[ch];
          
}



//########################################################################################
int MichelAna::GetGlobalBin(const TH3D* h, double x, double y, double z){
  int xbin = h->GetXaxis()->FindBin(x);
  int ybin = h->GetYaxis()->FindBin(y);
  int zbin = h->GetZaxis()->FindBin(z);
  return h->GetBin(xbin,ybin,zbin);
}
int MichelAna::GetGlobalBin(const TH3D* h, const TLorentzVector& v){
  return GetGlobalBin(h, v.X(), v.Y(), v.Z() );
}
int MichelAna::GetGlobalBin(const TH3D* h, const TVector3& v){
  return GetGlobalBin(h, v.X(), v.Y(), v.Z() );
}



//########################################################################################
float MichelAna::GetVisibility(TVector3 loc, int ch){
  return GetVisibility(loc, ch, 0) + GetVisibility(loc, ch, 1);
}



//########################################################################################
float MichelAna::GetVisibility(TVector3 loc, int ch, int isVis ){
  float tmin = -999; // -999 = don't fetch time info
  float tmean = 0.;
  float trms = 0.;
  float vis = 0.;
  ReadPhotonLibrary(loc, ch, isVis, vis, tmin, tmean, trms);
  return vis;
}



//########################################################################################
void MichelAna::ReadPhotonLibrary(TVector3 &loc, int ch, int isVis, float& vis, float& tmin, float& tmean, float& trms){
  
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
//float MichelAna::ElectronsFromADCArea( float area, int plane ) {
//  return area / fCalAreaConstants[plane];
//}

// ########################################
float MichelAna::dQdx_from_dEdx(float dEdx ){
  float rho    = 1.383; // LAr density in g/cm^3 (from ArgoNeuT)
  float Beta    = util::kModBoxB / (rho * fEfield);
  float Alpha   = util::kModBoxA;
  return log( Beta*dEdx + Alpha ) / (Beta*fWion); 
}

// ########################################################################################
float MichelAna::Integrate(const TH1D* h, float x1, float x2){
  int Nbins = h->GetNbinsX();
  double integral = 0.;
  for(int i=1; i<Nbins; i++ ){
    double x = h->GetBinLowEdge(i);
    if( x >= x1 and x <= x2 ) integral += h->GetBinContent(i);
  }
  return integral;
}



//########################################################################################
float MichelAna::Integrate(const std::vector<float> &vec, int x1, int x2){
  float integral = 0.;
  for(int i=x1; i<x2; i++ ) integral += vec.at(i);
  return integral;
}



//########################################################################################
void MichelAna::AddToAverageWfm( std::vector<float> &wfm, short t1, short t2, TH1D* h, int &counter, int polarity){
  int nbins = h->GetNbinsX();
  short T = t2 - t1;
  if( nbins ==  T ) {
    float xmin = h->GetXaxis()->GetXmin();
    for(short i=0; i<T; i++) h->Fill(xmin+i, polarity*wfm.at(t1+i) );
    counter++;
  } else {
    LOG_VERBATIM("MichelAna")<<"MichelAna::AddToAverageWfm ERROR: Number of bins ("<<nbins<<") does not match range given ("<<T<<")";
  }
}

//########################################################################################
void MichelAna::AddToAverageWfm( const TH1D* hin, short t1, short t2, short thit, TH1D* h, int &counter, int polarity){
  LOG_VERBATIM("MichelAna")<<"Adding to average MC waveform: "<<t1<<", "<<t2<<", thit = "<<thit;
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
void MichelAna::MakeWfmGraphs(art::Event& e){

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
        g_bs.SetMarkerColor(kBlue);
        g_bs.SetLineColor(kBlue);
        TGraph g_grad;
        g_grad.SetMarkerColor(kViolet+2);
        g_grad.SetLineColor(kViolet+2);
        std::vector<float> g = fOpHitBuilderAlg.MakeGradient(fPMT_wfm[ch]);
        for(size_t i=0; i<fPMT_wfm[ch].size(); i++) {
          g_wfm .SetPoint(  g_wfm.GetN(),  i, -1.*fPMT_wfm[ch].at(i));
          g_grad.SetPoint(  g_grad.GetN(), i, -1.*g.at(i) );
          if( fvbs[ch].size() > 0 ){
            g_bs.SetPoint( g_bs.GetN(), i, -1.*fvbs[ch].at(i) );
          }
        }
         
        sprintf(histName, "wfm_ch%lu_%i_r%i_sr%i_e%i", ch, fNumSavedWfmGraphs, e.run(), e.subRun(), e.id().event());
        fTCanvas     = graphDir_wfms.make<TCanvas>(histName,"c1",700,600);
        TPad c1("c1","c1",0.0,0.0,1.0,1.0);
        c1.Draw(); 
        c1.Divide(1,2);
          
        c1.cd(1);

        TMultiGraph mg;
        if(g_wfm.GetN() > 0 )    mg.Add(&g_wfm,"APL");
        //if(g_bs.GetN() > 0 ) mg.Add(&g_bs,"L");
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
void MichelAna::MakeClusteringGraphs( MichelCluster& cl ) {
  
  if( fNumSavedHitGraphs >= fMaxSavedHitGraphs ) return;
  
  LOG_VERBATIM("MichelAna")<<"Making clustering graph";
  
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
    sprintf(buffer, "clstr%i_r%i_sr%i_e%i_pl%d_Reco-%3.0fMeV_True-%3.0fMeV_%3.0fdeg", fNumSavedHitGraphs, fRunNumber, fSubRunNumber, fEventNumber, plane, elShowerEnergy, fTrue_ElShowerEnergyDep, fDecayAngle2D);
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
//  g2.GetYaxis()->SetTitle("Charge density [e-/cm]");
  g2.GetYaxis()->SetTitle("Charge [e-]");
  g2.GetYaxis()->SetNoExponent(false);
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
  g3.GetXaxis()->SetTitle("Projected 2D distance [cm]");
  g3.GetYaxis()->SetTitle("Local linearity #chi^{2}");
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

//########################################################################################
void MichelAna::Display3DEvent( std::vector<TVector3> &muHits, std::vector<TVector3> &elHits ) {
}

float MichelAna::T_from_X(float x, int pl){
  float t1 = fSamplingRate*(fXTicksOffset[0] - fTriggerOffset );
  float T  = fSamplingRate*(fXTicksOffset[pl] - fTriggerOffset );
  if( x > 0.2 ) {
    T += t1 + (x-0.2)/fDriftVelocity[0];
  } else {
    T = t1*(x/0.2);
  }
  return T;
}

DEFINE_ART_MODULE(MichelAna)
