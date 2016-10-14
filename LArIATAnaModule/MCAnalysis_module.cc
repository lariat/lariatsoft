////////////////////////////////////////////////////////////////////////
// Class:       MCAnalysis
// Module Type: analyzer
// File:        MCAnalysis_module.cc
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
#include <TStyle.h>
#include <THStack.h>
const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxTrackHits  = 1000;  //maximum number of space points
const int kMaxTrajHits   = 1000;  //maximum number of trajectory points
const int kMaxPrimaries  = 20000;  //maximum number of primary particles
const int kMaxTruePrimaryPts = 100000; //maximum number of truth trajectory points
// ### Number of centimeters in Z we require a track ###
// ### to have a space point within (default = 2 cm) ###
const double FirstSpacePointZPos = 2.0;


// ##############################################################
// ### Delta X Between Extrapolated MC Particle and TPC Track ###
// ##############################################################
const double DeltaXLowerBound = -2.0;
const double DeltaXUpperBound = 6.0;

// ##############################################################
// ### Delta Y Between Extrapolated MC Particle and TPC Track ###
// ##############################################################
const double DeltaYLowerBound = -3.0;
const double DeltaYUpperBound = 6.0;


// ########################################################################
// ### Fiducial Boundry Cuts (used to determine if a track is stopping) ###
// ########################################################################
const double XLowerFid = 0;
const double XUpperFid = 47;

const double YLowerFid = -20;
const double YUpperFid = 20;

const double ZLowerFid = 0;
const double ZUpperFid = 90;


// ########################################################################
// ### Definition of the upstream part of the TPC where we restrict the ###
// ###             number of tracks which can be present                ###
// ########################################################################
const int UpperPartOfTPC = 14.0;

// #####################################################
// ### Number of tracks allowed in the upstream part ###
// #####################################################
const int nLowZTracksAllowed = 4;


// ############################
// ### Alpha Cut in degrees ###
// ############################
const double alphaCut = 10;

//Create the cross section from the incident and interaction plots
const float rho = 1400; //kg/m^3
//  float cm_per_m = 100;
const float molar_mass = 39.9; //g/mol
const float g_per_kg = 1000; 
const float avogadro = 6.02e+23; //number/mol
const float number_density = rho*g_per_kg/molar_mass*avogadro;
const float slab_width = 0.0045;//in m



namespace lariat 
{
   class MCAnalysis;
}

class lariat::MCAnalysis : public art::EDAnalyzer 
{
public:
   explicit MCAnalysis(fhicl::ParameterSet const & p);
   virtual ~MCAnalysis();

   // Required functions.
   void analyze(art::Event const & e) override;

   // Selected optional functions.
   void beginJob();
   void reconfigure(fhicl::ParameterSet const & p);
   void endJob() override;
private:

   // === Function used to reset all the variables  ===
   void ResetVars();
  
   //=== Storing Run Information ===
   int run;			//<---Run Number
   int subrun;			//<---SubRun Number
   int event;			//<---Event Number

   double trkidpri;
   double trkidmot;
   
   //####### dummy variable ##############//
   double dummyTrkX[kMaxTrack];
   double dummyTrkY[kMaxTrack];
   double dummyTrkZ[kMaxTrack];
   double dummyTrk_pHat0X[kMaxTrack];
   double dummyTrk_pHat0Y[kMaxTrack];
   double dummyTrk_pHat0Z[kMaxTrack];
   double dummyTrk_Index[kMaxTrack];
   double dummyTrk_Theta[kMaxTrack];
   double dummyTrk_Phi[kMaxTrack];
   // === Storing the tracks Calorimetry Information
   int    trkhits[kMaxTrack][2];
   //  double trkpida[kMaxTrack][2];
   double trkke[kMaxTrack][2];
   double trkdedx[kMaxTrack][2][1000];
   double trkrr[kMaxTrack][2][1000];
   double trkpitchhit[kMaxTrack][2][1000];
   
   
   double trjPt_dmX[kMaxTrack][kMaxTrajHits];     //<---Storing the trajector point location in X
   double trjPt_dmY[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Y
   double trjPt_dmZ[kMaxTrack][kMaxTrajHits];	//<---Storing the trajector point location in Z

   // === dummy variable ===


   double dummyXpoint, dummyYpoint, dummyZpoint;
   double dummypoint_TempTrjX, dummypoint_TempTrjY,dummypoint_TempTrjZ,dummyPointTrkInd;
   int nEvtsMCTrackMatch=0;
   int nEventsPassingAlpha=0;
   int nEvtsTrackZPos =0;
   int nEvtsGoodMC = 0;
   int nLowZTrkEvents =0;
   float mcPhi = 0;
   float mcTheta = 0;
   double NTpts;
   double alpha, DeltaX, DeltaY, DeltaZ;
   
   // === Storing Geant4 MC Truth Information ===
   double g4Primary_X0[kMaxPrimaries];
   double g4Primary_Y0[kMaxPrimaries];
   double g4Primary_Z0[kMaxPrimaries];
   double g4Primary_Px[kMaxPrimaries];
   double g4Primary_Py[kMaxPrimaries];
   double g4Primary_Pz[kMaxPrimaries];
   double g4Primary_Xf[kMaxPrimaries];
   double g4Primary_Yf[kMaxPrimaries];
   double g4Primary_Zf[kMaxPrimaries];
   double g4Primary_ProjX0[kMaxPrimaries];
   double g4Primary_ProjY0[kMaxPrimaries];
   double g4Primary_ProjZ0[kMaxPrimaries];
   double g4PrimaryProcess[kMaxPrimaries];


	// === Storing additionnal Geant4 MC Truth Information for the primary track only ===	   
	//int NTrTrajPts[kMaxPrimaries];							 //<--Nb. of true points in the true primary trajectories
   double MidE[kMaxTruePrimaryPts];  //<--E Energy of a point in the true primary trajectory 

//################# Histogram for the Analysis ####################

// ###################################   
// ### MC Starting Point Histogram ###
// ###################################
TH1D *fStartXMC; 
TH1D *fStartYMC;
TH1D *fStartZMC;

// ######################################   
// ### MC Starting Momentum Histogram ###
// ######################################
TH1D *fStartPxMC;
TH1D *fStartPyMC;
TH1D *fStartPzMC;

// #################################   
// ### MC Ending Point Histogram ###
// #################################
TH1D *fEndXMC;
TH1D *fEndYMC;
TH1D *fEndZMC; 

// ################################################   
// ### Extrapolated MC Starting Point Histogram ###
// ################################################   
TH1D *fProjX0;
TH1D *fProjY0;
TH1D *fProjZ0;
TH1D *fProcess;
   

// ###################################
// ### Primary Particle Z_f vs X_f ###
// ###################################
TH2D *fPrimaryEndXvsZ;

// ###################################
// ### Primary Particle Y_f vs Z_f ###
// ###################################
TH2D *fPrimaryEndYvsZ;

// ##########################
// ### Histogram for cuts ###
// ##########################
TH1D *fAlpha;
TH1D *fDeltaX;
TH1D *fDeltaY;
TH1D *fDeltaZ;
TH1D *fUpstZpts;
TH1D *fMCInitalKE;
TH1D *fDeltaEndX;
TH1D *fDeltaEndY;
TH1D *fDeltaEndZ;
TH1F *fCutHistogram;
TH1D *fnUpstmTrk;

// ############################################
// ### Histogram for Delta End Z vs process ###
// ############################################
TH1D *fDeltaEndZInElastic;
TH1D *fDeltaEndZNeutronInElastic;
TH1D *fDeltaEndZHadElastic;
TH1D *fDeltaEndZnCap;
TH1D *fDeltaEndZnuclearCapatureAtRest;
TH1D *fDeltaEndZDecay;
TH1D *fDeltaEndZKaonZeroInElastic;
TH1D *fDeltaEndZCoulombScat;
TH1D *fDeltaEndZMuMinusCapture;
TH1D *fDeltaEndZProtonInelastic;
TH1D *fDeltaEndZPiMinusAbsorptionAtRest;


// ############################################
// ### Histogram for Delta End Y vs process ###
// ############################################
TH1D *fDeltaEndYInElastic;
TH1D *fDeltaEndYNeutronInElastic;
TH1D *fDeltaEndYHadElastic;
TH1D *fDeltaEndYnCap;
TH1D *fDeltaEndYnuclearCapatureAtRest;
TH1D *fDeltaEndYDecay;
TH1D *fDeltaEndYKaonZeroInElastic;
TH1D *fDeltaEndYCoulombScat;
TH1D *fDeltaEndYMuMinusCapture;
TH1D *fDeltaEndYProtonInelastic;
TH1D *fDeltaEndYPiMinusAbsorptionAtRest;

// ############################################
// ### Histogram for Delta End X vs process ###
// ############################################
TH1D *fDeltaEndXInElastic;
TH1D *fDeltaEndXNeutronInElastic;
TH1D *fDeltaEndXHadElastic;
TH1D *fDeltaEndXnCap;
TH1D *fDeltaEndXnuclearCapatureAtRest;
TH1D *fDeltaEndXDecay;
TH1D *fDeltaEndXKaonZeroInElastic;
TH1D *fDeltaEndXCoulombScat;
TH1D *fDeltaEndXMuMinusCapture;
TH1D *fDeltaEndXProtonInelastic;
TH1D *fDeltaEndXPiMinusAbsorptionAtRest;


// #########################################################
// ### Histograms tracks which go into the Cross-Section ###
// #########################################################
TH1D *fdataPiondEdX;
TH1D *fdataPionRR;
TH1D *fdataPionTrkPitch;
TH2D *fdataPiondEdXvsRR;

//################### Cross-section ####################################//
   TH1D *fdataPionIncidentKE;
   TH1D *fPionInteractions;
   TH1F *fCrossSection;
   TH1D *fTruthIncidentKE;
   TH1D *fTruthInteractingKE;
   
   
  //std::string fTrigModuleLabel;
  std::string fClusterModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fG4ModuleLabel;
  std::string fShowerModuleLabel;       // Producer that makes showers from clustering
  std::string fMCShowerModuleLabel;	// Producer name that makes MCShower Object

  calo::CalorimetryAlg fCalorimetryAlg;

};


lariat::MCAnalysis::MCAnalysis(fhicl::ParameterSet const & pset) 
  : EDAnalyzer(pset)
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
   this->reconfigure(pset);
}

lariat::MCAnalysis::~MCAnalysis()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::MCAnalysis::reconfigure(fhicl::ParameterSet const & pset)
{
  //fTrigModuleLabel	 	= pset.get< std::string >("TriggerUtility");
   fHitsModuleLabel      	= pset.get< std::string >("HitsModuleLabel");
   fTrackModuleLabel		= pset.get< std::string >("TrackModuleLabel");
   fCalorimetryModuleLabel 	= pset.get< std::string >("CalorimetryModuleLabel");
   fParticleIDModuleLabel  	= pset.get< std::string >("ParticleIDModuleLabel");
   fClusterModuleLabel          = pset.get< std::string >("ClusterModuleLabel");
   fG4ModuleLabel               = pset.get< std::string >("G4ModuleLabel");
   fShowerModuleLabel           = pset.get< std::string >("ShowerModuleLabel");
   fMCShowerModuleLabel		= pset.get< std::string >("MCShowerModuleLabel");
   return;
}

void lariat::MCAnalysis::analyze(art::Event const & evt)
{
// #############################################
// ### Reset variables before we get started ###
// #############################################
ResetVars();

// #######################################
// ### Get potentially useful services ###
// #######################################
// === Geometry Service ===
art::ServiceHandle<geo::Geometry> geom;
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
//						FILLING THE MCTruth Geant4 INFORMATION
// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------
int primary=0;

if(!isdata)
   {
   
   // ### Filling all MC Events in the Cut histogram ###
   fCutHistogram->Fill(0);
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
	
      
   // ### Setting a string for primary ###
   std::string pri("primary");
      
   // ### Setting a string for PionMinusInelastic ###
   std::string PionMinusInelastic("PionMinusInelastic");
      
   // ### Setting a string for NeutronInelastic ###
   std::string NeutronInelastic("NeutronInelastic");
      
    // ### Setting a string for hadElastic ###
   std::string hadElastic("hadElastic");
      
   // ### Setting a string for nCapture ###
   std::string nCapture("nCapture");
      
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
   std::string ProtonInelastic("ProtonInelastic");
      
   // ### Setting a string for PiMinusAbsorptionAtRest ###
   std::string PiMinusAbsorptionAtRest("PiMinusAbsorptionAtRest");


   int geant_particle=0;
   // float g4PrimaryProcess[100] = {0};
   int g4Primary_TrkID[100] = {999};
   // ############################################################
   // ### Determine the number of primary particles from geant ###
   // ############################################################
   for( unsigned int i = 0; i < geant_part.size(); ++i )
      {
      geant_particle++;
      
      // ##################################################
      // ### Grabbing the primary particles information ###
      // ##################################################
      if(geant_part[i]->Process()==pri)
         {
	 
	 // ##############################################################
	 // ### Getting the trajectory points for the primary particle ###
	 // ##############################################################
	 simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
	 
	 // ### Number of MC trajectory points ###
	 int iPrimPt = 0;	
	 for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	    {
	    MidE[iPrimPt]  = truetraj.E(iPrimPt)*1000;
	    iPrimPt++;
	    }//<--End loop on true trajectory points
         NTpts = iPrimPt;
         
	 // ### Getting the primary particles starting information ###
	 g4Primary_X0[primary] = geant_part[i]->Vx();
         g4Primary_Y0[primary] = geant_part[i]->Vy();
         g4Primary_Z0[primary] = geant_part[i]->Vz();

         g4Primary_Px[primary] = geant_part[i]->Px()*1000;
         g4Primary_Py[primary] = geant_part[i]->Py()*1000;
         g4Primary_Pz[primary] = geant_part[i]->Pz()*1000;
	 
	 // #########################################################################
	 // ## Extrapolating the MC Primary particle to the front face of the TPC ###
	 // #########################################################################
	 
	 double PX0 = -999, PY0 = -999; 
	 if(geant_part[i]->Vz() < 0)
	    {
	    double b1 = geant_part[i]->Vz() - geant_part[i]->Vx()*geant_part[i]->Pz()/geant_part[i]->Px();
	    double b2 = geant_part[i]->Vz() - geant_part[i]->Vy()*geant_part[i]->Pz()/geant_part[i]->Py();

            PX0 = -b1*geant_part[i]->Px()/geant_part[i]->Pz();
	    PY0 = -b2*geant_part[i]->Py()/geant_part[i]->Pz();

            g4Primary_ProjX0[primary] = PX0;
	    g4Primary_ProjY0[primary] = PY0;
	    g4Primary_ProjZ0[primary] = 0.0;
	    }//<---End projecting the particle if it started outside the TPC
 
       
         // ### Getting the primary particles ending information ###
	 g4Primary_TrkID[primary] = geant_part[i]->TrackId();
	 g4Primary_Xf[primary] =geant_part[i]->EndPosition()[0];
         g4Primary_Yf[primary] =geant_part[i]->EndPosition()[1];
         g4Primary_Zf[primary] =geant_part[i]->EndPosition()[2];
	 
	 // ##########################
	 // ### Filling Histograms ###
	 // ##########################
	 fProjX0->Fill(PX0);
	 fProjY0->Fill(PY0);
	 fProjZ0->Fill(0.0);
	 fStartXMC->Fill(geant_part[i]->Vx());
	 fStartYMC->Fill(geant_part[i]->Vy());
	 fStartZMC->Fill(geant_part[i]->Vz());
	 fStartPxMC->Fill(geant_part[i]->Px()*1000);
	 fStartPyMC->Fill(geant_part[i]->Py()*1000);
	 fStartPzMC->Fill(geant_part[i]->Pz()*1000);

	 fEndXMC->Fill(geant_part[i]->EndPosition()[0]);
	 fEndYMC->Fill(geant_part[i]->EndPosition()[1]);
	 fEndZMC->Fill(geant_part[i]->EndPosition()[2]);
         fPrimaryEndXvsZ->Fill(geant_part[i]->EndPosition()[2],geant_part[i]->EndPosition()[0]);
         fPrimaryEndYvsZ->Fill(geant_part[i]->EndPosition()[2],geant_part[i]->EndPosition()[1]);

	 // ### Counting the number of primaries
	 primary++;
	 
	 
	 } // End of if particle is primary
      trkidmot = geant_part[i]->Mother();
      }  // End of geant_part size


   // ############################################
   // ### Storing the process name for primary ###
   // ############################################
   for( unsigned int i = 0; i < geant_part.size(); ++i )
      {
      for( int jprime = 0; jprime < primary; jprime++)
         {
         if(geant_part[i]->Mother() == g4Primary_TrkID[jprime])
            {
	 
	    if(geant_part[i]->Process() == PionMinusInelastic)
	       {g4PrimaryProcess[primary -1] = 1;}
	     
	    if(geant_part[i]->Process() == NeutronInelastic)
	       {g4PrimaryProcess[primary -1] = 2;}
	     
	    if(geant_part[i]->Process() == hadElastic)
	       {g4PrimaryProcess[primary -1] = 3;}
	  
	    if(geant_part[i]->Process() == nCapture)
	       {g4PrimaryProcess[primary -1] = 4;}
	     
	    if(geant_part[i]->Process() == CHIPSNuclearCaptureAtRest)
	       {g4PrimaryProcess[primary -1] = 5;}
	  
	    if(geant_part[i]->Process() == Decay)
	       {g4PrimaryProcess[primary -1] = 6;}
	  
	    if(geant_part[i]->Process() == KaonZeroLInelastic)
	       {g4PrimaryProcess[primary -1] = 7;}
	     
	    if(geant_part[i]->Process() == CoulombScat)
	       {g4PrimaryProcess[primary -1] = 8;}
	     
	    if(geant_part[i]->Process() == muMinusCaptureAtRest)
	       {g4PrimaryProcess[primary -1] = 9;}
	     
	    if(geant_part[i]->Process() == ProtonInelastic)
	       {g4PrimaryProcess[primary -1] = 10;}

            if(geant_part[i]->Process()==PiMinusAbsorptionAtRest)
               {g4PrimaryProcess[primary -1] = 0;}

            }//<---End getting the particles associated with the primary
         fProcess->Fill(g4PrimaryProcess[jprime]);
         }//<---End jprime loop
      }//<--End i loop over all geant4 particles
       
}//<---End checking if this is data   
  
  
  
//=======================================================================================================================
//				Only looking at events where the primary particle enters the TPC
//=======================================================================================================================
   
bool GoodMCEventInTPC = true;
   
// ##############################################
// ### Looping over all the primary particles ###
// ##############################################
for(int npri = 0; npri < primary; npri++)
   {
   if(g4Primary_Zf[npri] < 0){GoodMCEventInTPC = false;}
   
   }//<---End npri loop
   
// ### Filling the cut histogram if the primary enters the TPC ###
if(GoodMCEventInTPC){fCutHistogram->Fill(1);}



//=======================================================================================================================
//						Low Z Spacepoint Track Cut
//=======================================================================================================================
bool TrackTrjPtsZCut = false;
int nUpStreamTrk = 0;

// ##################################################
// ### Loop over all the Reconstructed TPC Tracks ###
// ##################################################
for(size_t iTrk = 0; iTrk<tracklist.size(); iTrk++)
   {
   TVector3 p_hat_dm0;
   // ### Resetting the variables for each track ###
   dummyXpoint = 999, dummyYpoint = 999, dummyZpoint = 999;
   // ########################################################
   // ### Looping over the trajectory points for the track ###
   // ########################################################
   for(size_t iTrjPt = 0; iTrjPt<tracklist[iTrk]->NumberTrajectoryPoints(); iTrjPt++)
      {
      p_hat_dm0 = tracklist[iTrk]->DirectionAtPoint(iTrjPt);
      
      // ### Need to understand this .... ###
      if( p_hat_dm0.Z() < 0 )
         {
	 p_hat_dm0.SetX(p_hat_dm0.X()*-1);
	 p_hat_dm0.SetY(p_hat_dm0.Y()*-1);
	 p_hat_dm0.SetZ(p_hat_dm0.Z()*-1);
	 }
      trjPt_dmX[iTrk][iTrjPt] = tracklist[iTrk]->LocationAtPoint(iTrjPt).X();
      trjPt_dmY[iTrk][iTrjPt] = tracklist[iTrk]->LocationAtPoint(iTrjPt).Y();
      trjPt_dmZ[iTrk][iTrjPt] = tracklist[iTrk]->LocationAtPoint(iTrjPt).Z();
      
      // ###########################################################################
      // ### Setting our dummypoints if this is the lowest Z point on this track ###
      // ###           and still within the active volume of the TPC             ###
      // ###########################################################################
      if(trjPt_dmZ[iTrk][iTrjPt] < dummyZpoint && trjPt_dmZ[iTrk][iTrjPt] > 0.0 && 
	 trjPt_dmY[iTrk][iTrjPt] > -20.0 && trjPt_dmY[iTrk][iTrjPt] < 20.0 &&       
	 trjPt_dmX[iTrk][iTrjPt] > 0.0 && trjPt_dmX[iTrk][iTrjPt] < 47 )
         {
	 dummyXpoint = trjPt_dmX[iTrk][iTrjPt];
	 dummyYpoint = trjPt_dmY[iTrk][iTrjPt];
	 dummyZpoint = trjPt_dmZ[iTrk][iTrjPt];
	 
	 dummypoint_TempTrjX = p_hat_dm0.X();
	 dummypoint_TempTrjY = p_hat_dm0.Y();
	 dummypoint_TempTrjZ = p_hat_dm0.Z();
	 
	 dummyPointTrkInd = iTrk;
	  
	 } //--End sorting for the lowest point in Z
 


      }  //---End iTrjPt loop

      // ################################################################
      // ### Record this track if the upstream point is less than 2cm ###
      // ################################################################
     
      if(dummyZpoint < 2)
       {
       	 dummyTrkX[nUpStreamTrk] = dummyXpoint;
	 dummyTrkY[nUpStreamTrk] = dummyYpoint;
	 dummyTrkZ[nUpStreamTrk] = dummyZpoint;


	  dummyTrk_pHat0X[nUpStreamTrk] = dummypoint_TempTrjX;
	  dummyTrk_pHat0Y[nUpStreamTrk] = dummypoint_TempTrjY;
	  dummyTrk_pHat0Z[nUpStreamTrk] = dummypoint_TempTrjZ;  
	  dummyTrk_Index[nUpStreamTrk] = dummyPointTrkInd;

	  nUpStreamTrk++;
          TrackTrjPtsZCut = true;
          
	} // End of If loop 
        
      
        //} //End of Spacept Loop
       fUpstZpts->Fill(dummyZpoint);


   }  // End iTrk


   fnUpstmTrk->Fill(nUpStreamTrk);

   // ###############################################
   // ### Skipping events that don't have a track ###
   // ###   in the front of the TPC (Z) Position  ###
   // ###############################################
   if(TrackTrjPtsZCut)
   {
   // ### Counting Events w/ front face TPC Track ###
   nEvtsTrackZPos++;
   fCutHistogram->Fill(2);
   
   
      
// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------
//						FILLING THE 3-D TRACK INFORMATION
// ----------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------
   
//####################Filling Trk Pitchhit, dE/dX info for all tracks###########################################


   // ### Looping over tracks ###
    for(size_t i=0; i<tracklist.size();++i)
      {

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
    

      }



   //=======================================================================================================================
   //					Cutting on the number of tracks in the upstream TPC
   //=======================================================================================================================
   
   int nLowZTracksInTPC = 0;
   // ################################################################
   // ### Initializing variables for study of low Z track location ###
   // ################################################################
   bool LowZTrackInTPC = false;
   
  // float templowz1 = 0;
   
   
   
    // #################################################################
    // ### Only keeping events if there is less than N tracks in the ###
    // ###    first ## cm of the TPC (to help cut out EM Showers     ###
    // #################################################################
    for(size_t iTrk = 0; iTrk<tracklist.size(); iTrk++)
      {
      // ### Start by assuming this track is not in the ###
      // ###          low Z part of the TPC             ###
      LowZTrackInTPC = false;
           
      // ##################################################
      // ### Looping over the spacepoints for the track ###
      // ##################################################
      for(size_t iTrjPt = 0; iTrjPt<tracklist[iTrk]->NumberTrajectoryPoints(); iTrjPt++)
         {
	 
	 // ##################################################
	 // ### Count this track if it has a spacepoint in ###
	 // ###       the low Z region of the TPC          ###
	 // ##################################################
	 if(trjPt_dmZ[iTrk][iTrjPt] < UpperPartOfTPC)
	    {  
	    if(trjPt_dmY[iTrk][iTrjPt] > YLowerFid && trjPt_dmY[iTrk][iTrjPt] < YUpperFid && 
	       trjPt_dmX[iTrk][iTrjPt] > XLowerFid && trjPt_dmX[iTrk][iTrjPt] < XUpperFid)
	        {LowZTrackInTPC = true; }
		
            }//<---End counting if 
	
         }//<---End nspts loop
      
      // ##################################################################
      // ### If the track was in the "UpperPartOfTPC", bump the counter ###
      // ##################################################################
      if(LowZTrackInTPC)
         {
	 
	 nLowZTracksInTPC++;
	 }//<---End counting track in the Upstream part

      }//<---End nTPCtrk
    
     
    // ### Skipping the event if there are too many ###
    // ###       low Z tracks in the event          ###
    if((nLowZTracksInTPC <= nLowZTracksAllowed) && (nLowZTracksInTPC > 0))
      {  
    fCutHistogram->Fill(3);
    // ### Counting the event if it passes ###
    nLowZTrkEvents++;
    
   //=======================================================================================================================
   //				Matching the MC Particle to the Reco Track (similar to the WC Track match)
   //=======================================================================================================================


   // ################################################
   // ### Calculating the angles for the Geant4 MC ###
   // ################################################
   TVector3 z_hat_MC(0,0,1);
   TVector3 p_hat_0_MC;

   // ### Setting the vector for the MC using the ###
   // ###  extrapolated Momentum vector   ###
   p_hat_0_MC.SetX(g4Primary_Px[0]);
   p_hat_0_MC.SetY(g4Primary_Py[0]);
   p_hat_0_MC.SetZ(g4Primary_Pz[0]); 


   // ### Getting everything in the same convention ###
   float mcPhi = 0;
   float mcTheta = 0;

   // === Calculating Theta for MC ===
   mcTheta = acos(z_hat_MC.Dot(p_hat_0_MC)/p_hat_0_MC.Mag());
 // std::cout<<"Truth Theta"<<mcTheta<<std::endl;
   // === Calculating Phi for MC ===
   if( p_hat_0_MC.Y() > 0 && p_hat_0_MC.X() > 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X()); }
   else if( p_hat_0_MC.Y() > 0 && p_hat_0_MC.X() < 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+3.141592654; }
   else if( p_hat_0_MC.Y() < 0 && p_hat_0_MC.X() < 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+3.141592654; }
   else if( p_hat_0_MC.Y() < 0 && p_hat_0_MC.X() > 0 ){ mcPhi = atan(p_hat_0_MC.Y()/p_hat_0_MC.X())+6.28318; }
   else if( p_hat_0_MC.Y() == 0 && p_hat_0_MC.X() == 0 ){ mcPhi = 0; }//defined by convention
   else if( p_hat_0_MC.Y() == 0 )
      {
      if( p_hat_0_MC.X() > 0 ){ mcPhi = 0; }

      else{ mcPhi = 3.141592654; }

      }
   else if( p_hat_0_MC.X() == 0 )
      {
      if( p_hat_0_MC.Y() > 0 ){ mcPhi = 3.141592654/2; }
      else{ mcPhi = 3.141592654*3/2; }

      }

    std::cout<<"Truth Phi"<<mcPhi<<std::endl;
   // ######################################################
   // ### Calculating the angles for the Upstream Tracks ###
   // ######################################################

   TVector3 z_hat_TPC(0,0,1);
   TVector3 p_hat_0_TPC;
   std::cout<<"nUpStreamTrk"<<nUpStreamTrk<<std::endl;
   for(int aa = 0; aa < nUpStreamTrk; aa++)
      {
      // ### Setting the TVector ###
      p_hat_0_TPC.SetX(dummyTrk_pHat0X[aa]);
      p_hat_0_TPC.SetY(dummyTrk_pHat0Y[aa]);
      p_hat_0_TPC.SetZ(dummyTrk_pHat0Z[aa]);
      
      
      // ### Calculating TPC track theta ###
      dummyTrk_Theta[aa] = acos(z_hat_TPC.Dot(p_hat_0_TPC)/p_hat_0_TPC.Mag());
      
      // ### Calculating TPC track phi ###
      if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X()); }
      else if( p_hat_0_TPC.Y() > 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() < 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+3.141592654; }
      else if( p_hat_0_TPC.Y() < 0 && p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = atan(p_hat_0_TPC.Y()/p_hat_0_TPC.X())+6.28318; }
      else if( p_hat_0_TPC.Y() == 0 && p_hat_0_TPC.X() == 0 ){ dummyTrk_Phi[aa] = 0; }//defined by convention
      else if( p_hat_0_TPC.Y() == 0 )
         {
         if( p_hat_0_TPC.X() > 0 ){ dummyTrk_Phi[aa] = 0; }

         else{ dummyTrk_Phi[aa] = 3.141592654; }

         }
      else if( p_hat_0_TPC.X() == 0 )
         {
         if( p_hat_0_TPC.Y() > 0 ){ dummyTrk_Phi[aa] = 3.141592654/2; }
         else{ dummyTrk_Phi[aa] = 3.141592654*3/2; }

         }
      
      std::cout<<"Track theta"<<dummyTrk_Theta[aa]<<"Track Phi"<<dummyTrk_Phi[aa]<<std::endl;
      }//<---End aa loop


   // =============================================================================================================
   // 				Cutting on the DeltaX and Delta Y between MC and Reco Track
   // =============================================================================================================

   bool MC_TPCMatch = false;
   int nDeltaMatch = 0;
   for(int bb = 0; bb < nUpStreamTrk; bb++)
       {

        DeltaX = dummyTrkX[bb] - g4Primary_ProjX0[0];
        DeltaY = dummyTrkY[bb] - g4Primary_ProjY0[0];
        DeltaZ = dummyTrkZ[bb] - g4Primary_ProjZ0[0];

        fDeltaX->Fill(DeltaX);
        fDeltaY->Fill(DeltaY);
        fDeltaZ->Fill(DeltaZ);

      // ### Matching in Delta X and Delta Y  ###
      if(DeltaX > DeltaXLowerBound && DeltaX < DeltaXUpperBound && DeltaY > DeltaYLowerBound && DeltaY < DeltaYUpperBound)
         {
	 MC_TPCMatch = true;
	 nDeltaMatch++;
	 }//<---End matching



	}  // End of bb Loop

       if((MC_TPCMatch) || (nDeltaMatch==1))

       {
       nEvtsMCTrackMatch++;
       fCutHistogram->Fill(4);

   // =============================================================================================================
   // 				Cutting on the Alpha Angle between MC and Reco Track
   // =============================================================================================================
   
   bool AlphaMatch = false;
   size_t RecoTrackIndex = -1;
   // #########################################################
   // ### Define the unit vectors for the MC and TPC Tracks ###
   // #########################################################
   TVector3 theUnitVector_MC;
   TVector3 theUnitVector_TPCTrack;
   
   // ###################################
   // ### Loop over all the US Tracks ###
   // ###################################
   for(int bb = 0; bb < nUpStreamTrk; bb++)
      {
      
      // === MC Unit Vector ===
      theUnitVector_MC.SetX(sin(mcTheta)*cos(mcPhi));
      theUnitVector_MC.SetY(sin(mcTheta)*sin(mcPhi));
      theUnitVector_MC.SetZ(cos(mcTheta));
   
      // === TPC Track Unit Vector ===
      theUnitVector_TPCTrack.SetX(sin(dummyTrk_Theta[bb])*cos(dummyTrk_Phi[bb]));
      theUnitVector_TPCTrack.SetY(sin(dummyTrk_Theta[bb])*sin(dummyTrk_Phi[bb]));
      theUnitVector_TPCTrack.SetZ(cos(dummyTrk_Theta[bb]));
      
      // ###########################################################
      // ### Calculating the angle between WCTrack and TPC Track ###
      // ###########################################################
      float alpha = ( acos(theUnitVector_MC.Dot(theUnitVector_TPCTrack)) )* (180.0/3.141592654);
      
      fAlpha->Fill(alpha);
      
      // ### Setting the boolian for the true match ###
      if(alpha < alphaCut)
         {
	 AlphaMatch = true;
	 RecoTrackIndex = dummyTrk_Index[bb];
	 }
      
      }//<---End bb loop


     if(AlphaMatch)
      {

      nEventsPassingAlpha++;
      fCutHistogram->Fill(5);


//===========================================================================================================================================
   // 						Calculate energy loss and fill histogram for cross-section analysis 
   // ===========================================================================================================================================   

   // ### The assumed energy loss between the cryostat and the TPC ###
   float entryTPCEnergyLoss = 8.6; //MeV

   // ### The assumed mass of the incident particle (here we assume a pion) ###
    float mass = 139.57;
   
   // #############################################################
   // ### Calculating the momentum from the MC Primary Particle ###
   // #############################################################
   float momentum = sqrt( (g4Primary_Px[0]*g4Primary_Px[0]) + (g4Primary_Py[0]*g4Primary_Py[0]) + (g4Primary_Pz[0]*g4Primary_Pz[0]) );
   
   // ###   Calculating the initial Kinetic Energy    ###
   // ### KE = Energy - mass = (p^2 + m^2)^1/2 - mass ###   
   float kineticEnergy = pow( (momentum*momentum) + (mass*mass) ,0.5) - mass;
   
   // ### The kinetic energy is that which we calculated ###
   // ###       minus the calculated energy loss         ###
   kineticEnergy -= entryTPCEnergyLoss;
   fMCInitalKE->Fill(kineticEnergy);

       // ########################################################################
   // ### Variables for the track we are calculating the cross-section for ###
   // ########################################################################
   double Piondedx[1000]={0.};
   double Pionresrange[1000]={0.};
   double Pionpitchhit[1000]={0.};
   size_t nPionSpts = 0;
   double PionSumEnergy = 0;
   
   // ### Variables for determining the matching of the end point ###
   float TrackEndX = 999, TrackEndY = 999, TrackEndZ = 999;
   std::vector<double> trackStart0;
   std::vector<double> trackEnd0;

   for(size_t iTrk = 0; iTrk<tracklist.size(); iTrk++)
      {

      tracklist[iTrk]->Extent(trackStart0,trackEnd0); 
      if(iTrk == RecoTrackIndex)
	{
      	TrackEndX = trackEnd0[0];
      	TrackEndY = trackEnd0[1];
      	TrackEndZ = trackEnd0[2];
      

      	nPionSpts = 0;

      	std::vector<art::Ptr<recob::SpacePoint> > spts0 = fmsp.at(iTrk);

      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
      	for(size_t nspts = 0; nspts<spts0.size(); ++nspts)
       	{
        
        	Piondedx[nPionSpts]     = trkdedx[iTrk][1][nspts];

          // ### Putting in a fix in the case that the dE/dX is negative in this step ### 
	 // ###  then take the point before and the point after and average them
		 if(Piondedx[nPionSpts] < 0 && nspts < spts0.size() && nspts > 0)
		    {Piondedx[nPionSpts] = ( (trkdedx[iTrk][1][nspts - 1] + trkdedx[iTrk][1][nspts + 1]) / 2);}
	 
	 // ### If this didn't fix it, then just put in a flat 2.4 MeV / cm fix ###
		 if(Piondedx[nPionSpts] < 0)
		    {
		    Piondedx[nPionSpts] = 2.4;
		    continue;
		    }

        	 Pionresrange[nPionSpts] = trkrr[iTrk][1][nspts];
        	 Pionpitchhit[nPionSpts] = trkpitchhit[iTrk][1][nspts];
         
		 PionSumEnergy = (Piondedx[nPionSpts] * Pionpitchhit[nPionSpts]) + PionSumEnergy;
	 
		 // ### Recording the dE/dX ###
		 fdataPiondEdX->Fill(Piondedx[nPionSpts]);
		 // ### Recording the residual range ###
		 fdataPionRR->Fill(Pionresrange[nPionSpts]);
		 // ### Recording the Pitch ###
		 fdataPionTrkPitch->Fill(Pionpitchhit[nPionSpts]);
	 
		 // ### Filling 2d dE/dX vs RR ###
		 fdataPiondEdXvsRR->Fill(Pionresrange[nPionSpts], Piondedx[nPionSpts]);
	 
		 nPionSpts++;

		} //End of Spacept loop

        	} // End of If Loop


      	} // End of Trk Loop

   // ###################################
   // ### Filling the Delta End Point ###
   // ###################################
   float DeltaEndX = g4Primary_Xf[0] - TrackEndX;
   float DeltaEndY = g4Primary_Yf[0] - TrackEndY;
   float DeltaEndZ = g4Primary_Zf[0] - TrackEndZ;


   fDeltaEndX->Fill(DeltaEndX);
   fDeltaEndY->Fill(DeltaEndY);
   fDeltaEndZ->Fill(DeltaEndZ);


   if(g4PrimaryProcess[0] == 1){fDeltaEndZInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 2){fDeltaEndZNeutronInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 3){fDeltaEndZHadElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 4){fDeltaEndZnCap->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 5){fDeltaEndZnuclearCapatureAtRest->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 6){fDeltaEndZDecay->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 7){fDeltaEndZKaonZeroInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 8){fDeltaEndZCoulombScat->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 9){fDeltaEndZMuMinusCapture->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 10){fDeltaEndZProtonInelastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 0){fDeltaEndZPiMinusAbsorptionAtRest->Fill(DeltaEndZ);}


    // ###########################################################
   // ### Looping over the Truth spacepoints to fill the histograms ###
   // ###########################################################

for(int i =0; i<NTpts; i++)
{
  
  fTruthIncidentKE->Fill(MidE[i]-139.57);
  if(i==NTpts-1)
  {
  fTruthInteractingKE->Fill(MidE[i]-139.57);
   }
}



   // ###########################################################
   // ### Looping over the spacepoints to fill the histograms ###
   // ###########################################################
   for(size_t npoints = 0; npoints < nPionSpts; npoints++)
      {
      // ### Filling the incidient histogram ###
      fdataPionIncidentKE->Fill(kineticEnergy);
      
      // ### Filling the interaction histogram for the last spt ###
      if(npoints == nPionSpts -1)
         {fPionInteractions->Fill(kineticEnergy);}
      
      //float energyLossInStep = Piondedx[npoints] * Pionresrange[npoints] * RecombinationFactor;
      float energyLossInStep = Piondedx[npoints] * Pionpitchhit[npoints];
      
      kineticEnergy -= energyLossInStep;
      
      
      }//<---End npoints loop

      }  // End of Alpha cut if loop

       } //End of Delta cut if loop

      } // End of Number of Upstream Track if loop 


   } //End of TrackTrjPtsZCut if loop





} // End of analyzer


void lariat::MCAnalysis::beginJob()
{
   
   	art::ServiceHandle<art::TFileService> tfs;
	fStartXMC	= tfs->make<TH1D>("fStartXMC","Primary X0",100,-50,50);
	fStartYMC	= tfs->make<TH1D>("fStartYMC","Primary Y0",100,-50,50);
	fStartZMC	= tfs->make<TH1D>("fStartZMC","Primary Z0",100,-20,100);

	fStartPxMC     = tfs->make<TH1D>("fStartPxMC","Primary Px",100,-150,150);
	fStartPyMC     = tfs->make<TH1D>("fStartPyMC","Primary Py",100,-150,150);
	fStartPzMC     = tfs->make<TH1D>("fStartPzMC","Primary Pz",100,-500,2500);

	fEndXMC        = tfs->make<TH1D>("fEndXMC","Primary Xf",100,-200,200);
	fEndYMC        = tfs->make<TH1D>("fEndYMC","Primary Yf",100,-200,200);
	fEndZMC        = tfs->make<TH1D>("fEndZMC","Primary Zf",100,-50,500);

        fProjX0        = tfs->make<TH1D>("fProjX0","Primary Extrapolated X0",100,-50,50);
	fProjY0        = tfs->make<TH1D>("fProjY0","Primary Extrapolated Y0",100,-50,50);
	fProjZ0        = tfs->make<TH1D>("fProjZ0","Primary Extrapolated Z0",100,-20,100);
        fProcess       = tfs->make<TH1D>("fProcess","Process name", 100, 0 , 10);
	
	fPrimaryEndXvsZ= tfs->make<TH2D>("fPrimaryEndXvsZ","Xf vs Zf", 500, -50, 450, 400, -200, 200);
        fPrimaryEndYvsZ= tfs->make<TH2D>("fPrimaryEndYvsZ","Yf vs Zf", 500, -50, 450, 200, -200, 200);

	fAlpha            = tfs->make<TH1D>("fAlpha","Angle between Truth and Reco", 90, 0,90);
        fDeltaX           = tfs->make<TH1D>("fDeltaX","Delta_X0 of Reco upstream pt and Truth Projected X0", 200, -50, 50);
        fDeltaY           = tfs->make<TH1D>("fDeltaY","Delta_Y0 of Reco upstream pt and Truth Projected Y0", 200, -50, 50);
        fDeltaZ           = tfs->make<TH1D>("fDeltaZ","Delta_Z0 of Reco upstream pt and Truth Projected Z0", 200, -50, 50);
        fUpstZpts         = tfs->make<TH1D>("fUpstZpts","Most upstream spacepoint of all TPC Tracks", 20, 0, 10);
        fMCInitalKE       = tfs->make<TH1D>("fMCInitalKE","Pion Initial KE ", 500, 0,2500);
        fDeltaEndX        = tfs->make<TH1D>("fDeltaEndX"," #Delta X_{f} of the most upstream Reco Track and the Projected Primary Particle 							X_{f}",200,-50,50);
        fDeltaEndY        = tfs->make<TH1D>("fDeltaEndY"," #Delta Y_{f} of the most upstream Reco Track and the Projected Primary Particle 							Y_{f}",200,-50,50);
        fDeltaEndZ        = tfs->make<TH1D>("fDeltaEndZ"," #Delta Z_{f} of the most upstream Reco Track and the Projected Primary  Particle 							Z_{f}",200,-50,50);

        fDeltaEndZInElastic        = tfs->make<TH1D>("fDeltaEndZInElastic","#Delta Z_{f} InElastic", 200, -50 , 50);
	fDeltaEndZNeutronInElastic = tfs->make<TH1D>("fDeltaEndZNeutronInElastic", "#Delta Z_{f} Neutron InElastic", 200, -50 , 50);
	fDeltaEndZHadElastic       = tfs->make<TH1D>("fDeltaEndZHadElastic","#Delta Z_{f} Hadronic Elastic", 200, -50 , 50);

        fDeltaEndZnCap             = tfs->make<TH1D>("fDeltaEndZnCap","#Delta Z_{f} Neutron Capture", 200, -50 , 50);
        fDeltaEndZnuclearCapatureAtRest = tfs->make<TH1D>("fDeltaEndZnuclearCapatureAtRest","#Delta Z_{f} Nuclear capture at rest ", 200, -50 , 50);

	fDeltaEndZDecay            = tfs->make<TH1D>("fDeltaEndZDecay","#Delta Z_{f} Decay ", 200, -50 , 50);
        fDeltaEndZKaonZeroInElastic = tfs->make<TH1D>("fDeltaEndZKaonZeroInElastic","#Delta Z_{f} Kaon Zero InElastic ", 200, -50 , 50);
        fDeltaEndZCoulombScat       = tfs->make<TH1D>("fDeltaEndZCoulombScat","#Delta Z_{f} Coulomb Scattering ", 200, -50 , 50);
        fDeltaEndZMuMinusCapture    = tfs->make<TH1D>("fDeltaEndZMuMinusCapture","#Delta Z_{f} MuMinus Capture ", 200, -50 , 50);
        fDeltaEndZProtonInelastic   = tfs->make<TH1D>("fDeltaEndZProtonInelastic","#Delta Z_{f} Proton Inelastic ", 200, -50 , 50);
        
	fDeltaEndZPiMinusAbsorptionAtRest = tfs->make<TH1D>("fDeltaEndZPiMinusAbsorptionAtRest","#Delta Z_{f} fDeltaEndZPiMinusAbsorptionAtRest ", 200, -50, 50);

        fnUpstmTrk    = tfs->make<TH1D>("fnUpstmTrk","# of Upstream Trk", 10,0,10);
        fCutHistogram = tfs->make<TH1F>("CutHistogram","Number of events left after each cut",6,0,6);

        fCutHistogram->GetXaxis()->SetTitle("0:NoCut, 1:Good Events, 2:2cm , 3:14cm,4:Delta , 5:Alpha ");
        fCutHistogram->GetYaxis()->SetTitle("Remaining events");
        fdataPiondEdX     = tfs->make<TH1D>("fdataPiondEdX","Pion dE/dX", 200, 0,50);
	fdataPionRR       = tfs->make<TH1D>("fdataPionRR","Pion Residual range", 240, -10,110);
        fdataPionTrkPitch = tfs->make<TH1D>("fdataPionTrkPitch","Track Pitch", 1000, 0,5);
        fdataPiondEdXvsRR = tfs->make<TH2D>("fdataPiondEdXvsRR","dE/dX vs Residual Range",200, 0, 100, 200, 0, 50);
        fdataPionIncidentKE = tfs->make<TH1D>("fdataPionIncidentKE","Pion Incident Kinetic Energy (MeV)", 40, 0, 2000);
        fPionInteractions   = tfs->make<TH1D>("fPionInteractions", "Pion Out Kinetic Energy (MeV)", 40, 0, 2000);
        fCrossSection       = tfs->make<TH1F>("CrossSection","Cross Section",20,0.,1000.);
        fTruthIncidentKE    = tfs->make<TH1D>("fTruthIncidentKE"," Truth Pion Incident Kinetic Energy (MeV)", 40, 0, 2000);
        fTruthInteractingKE = tfs->make<TH1D>("fTruthInteractingKE"," Truth Pion Interacting Kinetic Energy (MeV)", 40, 0, 2000);

        
	


	fStartXMC->GetXaxis()->SetTitle("cm");
	fStartYMC->GetXaxis()->SetTitle("cm");
        fStartZMC->GetXaxis()->SetTitle("cm");
        fStartPxMC->GetXaxis()->SetTitle("MeV");
        fStartPyMC->GetXaxis()->SetTitle("MeV");
        fStartPzMC->GetXaxis()->SetTitle("MeV");
        fEndXMC->GetXaxis()->SetTitle("cm");
	fEndYMC->GetXaxis()->SetTitle("cm");
        fEndZMC->GetXaxis()->SetTitle("cm");
        fProjX0->GetXaxis()->SetTitle("cm");
        fProjY0->GetXaxis()->SetTitle("cm");
	fProjZ0->GetXaxis()->SetTitle("cm");
	fProcess->GetXaxis()->SetTitle("process number");
	fPrimaryEndXvsZ->GetXaxis()->SetTitle("Z_f(cm)");
	fPrimaryEndXvsZ->GetYaxis()->SetTitle("X_f(cm)");
	fPrimaryEndYvsZ->GetXaxis()->SetTitle("Z_f(cm)");
	fPrimaryEndYvsZ->GetYaxis()->SetTitle("Y_f(cm)");
	fAlpha->GetXaxis()->SetTitle("Deg");
	fDeltaX->GetXaxis()->SetTitle("cm");
	fDeltaY->GetXaxis()->SetTitle("cm");
	fDeltaZ->GetXaxis()->SetTitle("cm");
	fMCInitalKE->GetXaxis()->SetTitle("MeV");
	fDeltaEndX->GetXaxis()->SetTitle("cm");
	fDeltaEndY->GetXaxis()->SetTitle("cm");
	fDeltaEndZ->GetXaxis()->SetTitle("cm");
	fdataPiondEdX->GetXaxis()->SetTitle("MeV/cm");
	fdataPionRR->GetXaxis()->SetTitle("cm");
	fdataPiondEdXvsRR->GetXaxis()->SetTitle("cm");
	fdataPiondEdXvsRR->GetXaxis()->SetTitle("MeV/cm");
	fdataPionIncidentKE->GetXaxis()->SetTitle("MeV"); 
	fPionInteractions->GetXaxis()->SetTitle("MeV");
        fCrossSection->GetXaxis()->SetTitle("Kinetic Energy, MeV");
        fCrossSection->GetYaxis()->SetTitle("Cross Section, barns");
 
   
}  //End of beginJob loop


void lariat::MCAnalysis::endJob()
{

  // Implementation of optional member function here.
  //Create the cross section from the incident and interaction plots
  float rho = 1400; //kg/m^3
  //  float cm_per_m = 100;
  float molar_mass = 39.9; //g/mol
  float g_per_kg = 1000; 
  float avogadro = 6.02e+23; //number/mol
  float number_density = rho*g_per_kg/molar_mass*avogadro;
  float slab_width = 0.0045;//in m


// ###################################################################
// #### Looping over the exiting bins to extract the cross-section ###
// ###################################################################


   for( int iBin = 1; iBin <= fPionInteractions->GetNbinsX(); ++iBin )
   {

    if( fdataPionIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad

       // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = (fPionInteractions->GetBinContent(iBin)/fdataPionIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);
   
   //std::cout<<"Cross-Section before conversion to barns = "<<TempCrossSection<<std::endl;
   
   float crossSection = TempCrossSection * (1/1e-28); //To put this into barns
   //std::cout<<"Cross-Section = "<<crossSection<<std::endl;
   
   
   fCrossSection->SetBinContent(iBin,crossSection);
   
   float numError = pow(fPionInteractions->GetBinContent(iBin),0.5);
   float num = fPionInteractions->GetBinContent(iBin);

   
   if(num == 0){num = 1;}
   float term1 = numError/num;
   //std::cout<<"term1 = "<<term1<<std::endl;
   
   float denomError = pow(fdataPionIncidentKE->GetBinContent(iBin),0.5);
   float denom = fdataPionIncidentKE->GetBinContent(iBin);
   if(denom == 0){denom = 1;}
   float term2 = denomError/denom;
   //std::cout<<"term2 = "<<term2<<std::endl;
   float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
   //std::cout<<"totalError = "<<totalError<<std::endl;
   fCrossSection->SetBinError(iBin,totalError);



   }

} // End of endJob loop




void lariat::MCAnalysis::ResetVars()
{

  run = -99999;
  subrun = -99999;
  event = -99999;
  dummyXpoint =-99999;
  dummyYpoint = -99999;
  dummyZpoint = -99999;

  for (int i = 0; i < kMaxTrack; ++i)

{
    dummyTrkX[i] = -99999;
    dummyTrkY[i] = -99999;
    dummyTrkZ[i] = -99999;
    dummyTrk_pHat0X[i] = -99999;
    dummyTrk_pHat0Y[i] = -999999;
    dummyTrk_pHat0Z[i] =-999999;
    dummyTrk_Index[i] =-999999;
    dummyTrk_Theta[i] = -99999;
    dummyTrk_Phi[i]  = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trjPt_dmX[i][j] = -99999;
      trjPt_dmY[i][j] = -99999;
      trjPt_dmZ[i][j] = -99999;

    }

    for (int j = 0; j<2; ++j){
 
      trkhits[i][j] = -99999; 
      for (int k = 0; k<1000; ++k){
	trkdedx[i][j][k] = -99999;
	trkrr[i][j][k] = -99999;
	trkpitchhit[i][j][k] = -99999;
      }
    }
    
  }
  for (int i = 0; i<kMaxPrimaries; ++i){
 
    g4Primary_X0[i] =-9999;
    g4Primary_Y0[i] =-9999;
    g4Primary_Z0[i] =-9999;
    g4Primary_Px[i] =-99999;
    g4Primary_Py[i] =-99999;
    g4Primary_Pz[i] =-99999;
    g4Primary_Xf[i] =-99999; 
    g4Primary_Yf[i] =-99999; 
    g4Primary_Zf[i] =-99999;
    g4Primary_ProjX0[i] =-99999;
    g4Primary_ProjY0[i] =-99999; 
    g4Primary_ProjZ0[i] =-99999;
    g4PrimaryProcess[i] =-99999;
}
 

		for(int j = 0; j<kMaxTruePrimaryPts; j++){	
   		MidE[j] = -99999; 
		}
   

}

DEFINE_ART_MODULE(lariat::MCAnalysis)
