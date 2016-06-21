#define Proton_cxx
#include "Proton.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


// ===================================================================================================================
// ====================================       PUT HISTOGRAMS HERE           ==========================================
// ===================================================================================================================

/////////////////////////////////// Primary Particle Start X Position //////////////////////////////////////////
TH1D *hMCPrimaryStartX = new TH1D("hMCPrimaryStartX", "Primary Particle X_{0}", 200, -50 , 50);
/////////////////////////////////// Primary Particle Start Y Position //////////////////////////////////////////
TH1D *hMCPrimaryStartY = new TH1D("hMCPrimaryStartY", "Primary Particle Y_{0}", 200, -50 , 50);
/////////////////////////////////// Primary Particle Start Z Position //////////////////////////////////////////
TH1D *hMCPrimaryStartZ = new TH1D("hMCPrimaryStartZ", "Primary Particle Z_{0}", 600, -120 , 180);

/////////////////////////////////// Primary Projected Particle Start X Position //////////////////////////////////////////
TH1D *hMCPrimaryProjectedStartX = new TH1D("hMCPrimaryProjectedStartX", "Primary Particle X_{0}", 200, -50 , 50);
/////////////////////////////////// Primary Projected Particle Start Y Position //////////////////////////////////////////
TH1D *hMCPrimaryProjectedStartY = new TH1D("hMCPrimaryProjectedStartY", "Primary Particle Y_{0}", 200, -50 , 50);
/////////////////////////////////// Primary Projected Particle Start Z Position //////////////////////////////////////////
TH1D *hMCPrimaryProjectedStartZ = new TH1D("hMCPrimaryProjectedStartZ", "Primary Particle Z_{0}", 400, -50 , 150);

/////////////////////////////////// Primary Particle End X Position //////////////////////////////////////////
TH1D *hMCPrimaryEndX = new TH1D("hMCPrimaryEndX", "Primary Particle X_{f}", 400, -200 , 200);
/////////////////////////////////// Primary Particle End Y Position //////////////////////////////////////////
TH1D *hMCPrimaryEndY = new TH1D("hMCPrimaryEndY", "Primary Particle Y_{f}", 400, -200 , 200);
/////////////////////////////////// Primary Particle End Z Position //////////////////////////////////////////
TH1D *hMCPrimaryEndZ = new TH1D("hMCPrimaryEndZ", "Primary Particle Z_{f}", 600, -120 , 480);

/////////////////////////////////// Primary Particle Px  //////////////////////////////////////////
TH1D *hMCPrimaryPx = new TH1D("hMCPrimaryPx", "Primary Particle P_{x}", 300, -150 , 150);
/////////////////////////////////// Primary Particle Py  //////////////////////////////////////////
TH1D *hMCPrimaryPy = new TH1D("hMCPrimaryPy", "Primary Particle P_{y}", 300, -159 , 150);
/////////////////////////////////// Primary Particle Pz //////////////////////////////////////////
TH1D *hMCPrimaryPz = new TH1D("hMCPrimaryPz", "Primary Particle P_{z}", 3000, -500 , 2500);

/////////////////////////////////// Primary Particle Process //////////////////////////////////////////
TH1D *hMCPrimaryProcess = new TH1D("hMCPrimaryProcess", "Primary Particle Process", 100, 0 , 10);

/////////////////////////////////// Primary End X vs Z Position //////////////////////////////////////////////
TH2D *hMCPrimaryEndXvsZ = new TH2D("hMCPrimaryEndXvsZ", "X_{f} vs Z_{f}", 600, -150, 450, 400, -200, 200);
/////////////////////////////////// Primary End Y vs Z Position //////////////////////////////////////////////
TH2D *hMCPrimaryEndYvsZ = new TH2D("hMCPrimaryEndYvsZ", "Y_{f} vs Z_{f}", 600, -150, 450, 200, -200, 200);

/////////////////////////////////// Energy Loss in the upstream region of the beamline ///////////////////////
TH1D *hMCELossUpstream = new TH1D("hMCELossUpstream", "Energy loss prior to entering the TPC", 1000, 0, 1000);

/////////////////////////////////// True Length //////////////////////////////////////////
TH1D *hTrueLength = new TH1D("hTrueLength", "#True Length of the Primary Particle inside the TPC", 200, 0 , 100);



// ================================================================================================================
// 						Reconstructed Information 
// ================================================================================================================

/////////////////////////////////// Most Upstream Z point of tracks //////////////////////////////////////////
TH1D *hdataUpstreamZPos = new TH1D("hdataUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);

/////////////////////////////////// Number of Tracks in the TPC versus distance //////////////////////////////////////////
TH2D *hdataNTracksvsZpos = new TH2D("hdataNTracksvsZpos", "Number of TPC tracks vs Z ", 30, 0, 30, 20, 0, 10);

/////////////////////////////////// Alpha Between WC and TPC Tracks //////////////////////////////////////////
TH1D *hAlpha = new TH1D("hAlpha", "#alpha between MC Particle and TPC Track", 90, 0, 90);

/////////////////////////////////// Delta Start X Position //////////////////////////////////////////
TH1D *hDeltaX = new TH1D("hDeltaX", "#Delta X_{0} of the most upstream Reco Track and the Projected Primary Particle X_{0}", 200, -50 , 50);
/////////////////////////////////// Delta Start Y Position //////////////////////////////////////////
TH1D *hDeltaY = new TH1D("hDeltaY", "#Delta Y_{0} of the most upstream Reco Track and the Projected Primary Particle Y_{0}", 200, -50 , 50);
/////////////////////////////////// Delta Start Z Position //////////////////////////////////////////
TH1D *hDeltaZ = new TH1D("hDeltaZ", "#Delta Z_{0} of the most upstream Reco Track and the Projected Primary Particle Z_{0}", 200, -50 , 50);


/////////////////////////////////// Reconstructed Length //////////////////////////////////////////
TH1D *hRecoLength = new TH1D("hRecoLength", "#Reconstructed Length of the Primary Particle inside the TPC", 200, 0 , 100);

/////////////////////////////////// Delta Length Position //////////////////////////////////////////
TH1D *hDeltaLength = new TH1D("hDeltaLength", "#Delta Length of the most upstream Reco Track and the Primary Particle Z_{0}", 200, -100 , 100);





/////////////////////////////////// "Proton" initial Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hMCInitalKE = new TH1D("hMCInitalKE", "Proton Initial Kinetic Energy (MeV)", 500, 0, 2500);
/////////////////////////////////// "Proton" dE/dX //////////////////////////////////////////
TH1D *hdataProtondEdX = new TH1D("hdataProtondEdX", "Proton dE/dX", 200, 0, 50);
/////////////////////////////////// "Proton" Residual Range //////////////////////////////////////////
TH1D *hdataProtonRR = new TH1D("hdataProtonRR", "Proton Residual Range", 240, -10, 110);
/////////////////////////////////// "Proton" Track Pitch //////////////////////////////////////////
TH1D *hdataProtonTrkPitch = new TH1D("hdataProtonTrkPitch", "Track Pitch", 1000, 0, 5);
///////////////////////////////// "Proton dE/dX vs RR ///////////////////////////////////////////
TH2D *hdataProtondEdXvsRR = new TH2D("", "dE/dX vs Residual Range",200, 0, 100, 200, 0, 50);
/////////////////////////////////// "Proton" Incident to the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hdataProtonIncidentKE = new TH1D("hdataProtonIncidentKE", "Proton Incident Kinetic Energy (MeV)", 40, 0, 2000);
/////////////////////////////////// "Proton" Exiting the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *fProtonInteractions = new TH1D("fProtonInteractions", "Proton Out Kinetic Energy (MeV)", 40, 0, 2000);
/////////////////////////////////// Cross-Section //////////////////////////////////////////////////////////////////////
TH1F *fCrossSection = new TH1F("fCrossSection", "Cross-Section", 40, 0, 2000); 

////////////////////////////////// Proton Track Start and End Positions //////////////////////////////////////////////////
TH1D *hdataProtonTrackEndX = new TH1D("hdataProtonTrackEndX", "Proton Track End X Position", 50, 0, 50);
TH1D *hdataProtonTrackEndY = new TH1D("hdataProtonTrackEndY", "Proton Track End Y Position", 40, -20, 20);
TH1D *hdataProtonTrackEndZ = new TH1D("hdataProtonTrackEndZ", "Proton Track End Z Position", 100, 0, 100);

TH1D *hdataProtonTrackStartX = new TH1D("hdataProtonTrackStartX", "Proton Track Start X Position", 50, 0, 50);
TH1D *hdataProtonTrackStartY = new TH1D("hdataProtonTrackStartY", "Proton Track Start Y Position", 40, -20, 20);
TH1D *hdataProtonTrackStartZ = new TH1D("hdataProtonTrackStartZ", "Proton Track Start Z Position", 100, 0, 100);


/////////////////////////////////// Delta End X Position //////////////////////////////////////////
TH1D *hDeltaEndX = new TH1D("hDeltaEndX", "#Delta X_{f} of the most upstream Reco Track and the Projected Primary Particle X_{f}", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position //////////////////////////////////////////
TH1D *hDeltaEndY = new TH1D("hDeltaEndY", "#Delta Y_{f} of the most upstream Reco Track and the Projected Primary Particle Y_{f}", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position //////////////////////////////////////////
TH1D *hDeltaEndZ = new TH1D("hDeltaEndZ", "#Delta Z_{f} of the most upstream Reco Track and the Projected Primary Particle Z_{f}", 200, -50 , 50);




//-----------------------------------------------------------------------------------------------------------------------------------------------
//								Stacked Histograms

/////////////////////////////////// Delta End Z Position InElastic //////////////////////////////////////////
TH1D *hDeltaEndZInElastic = new TH1D("hDeltaEndZInElastic", "#Delta Z_{f} InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Neutron InElastic //////////////////////////////////////////
TH1D *hDeltaEndZNeutronInElastic = new TH1D("hDeltaEndZNeutronInElastic", "#Delta Z_{f} Neutron InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Hadrondic Elastic //////////////////////////////////////////
TH1D *hDeltaEndZHadElastic = new TH1D("hDeltaEndZHadElastic", "#Delta Z_{f} Hadronic Elastic", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Neutron Capture //////////////////////////////////////////
TH1D *hDeltaEndZnCap = new TH1D("hDeltaEndZnCap", "#Delta Z_{f} Neutron Capture", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Nuclear Capture at Rest //////////////////////////////////////////
TH1D *hDeltaEndZnuclearCapatureAtRest = new TH1D("hDeltaEndZnuclearCapatureAtRest", "#Delta Z_{f} Nuclear capture at rest ", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Decay //////////////////////////////////////////
TH1D *hDeltaEndZDecay = new TH1D("hDeltaEndZDecay", "#Delta Z_{f} Decay ", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position KaonZero InElastic //////////////////////////////////////////
TH1D *hDeltaEndZKaonZeroInElastic = new TH1D("hDeltaEndZKaonZeroInElastic", "#Delta Z_{f} Kaon Zero InElastic ", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Coulomb Scattering //////////////////////////////////////////
TH1D *hDeltaEndZCoulombScat = new TH1D("hDeltaEndZCoulombScat", "#Delta Z_{f} Coulomb Scattering ", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position MuMinus Capture //////////////////////////////////////////
TH1D *hDeltaEndZMuMinusCapture = new TH1D("hDeltaEndZMuMinusCapture", "#Delta Z_{f} MuMinus Capture ", 200, -50 , 50);
/////////////////////////////////// Delta End Z Position Proton Inelastic //////////////////////////////////////////
TH1D *hDeltaEndZProtonInelastic = new TH1D("hDeltaEndZProtonInelastic", "#Delta Z_{f} Proton Inelastic ", 200, -50 , 50);





/////////////////////////////////// Delta End Y Position InElastic //////////////////////////////////////////
TH1D *hDeltaEndYInElastic = new TH1D("hDeltaEndYInElastic", "#Delta Y_{f} InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Neutron InElastic //////////////////////////////////////////
TH1D *hDeltaEndYNeutronInElastic = new TH1D("hDeltaEndYNeutronInElastic", "#Delta Y_{f} Neutron InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Hadrondic Elastic //////////////////////////////////////////
TH1D *hDeltaEndYHadElastic = new TH1D("hDeltaEndYHadElastic", "#Delta Y_{f} Hadronic Elastic", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Neutron Capture //////////////////////////////////////////
TH1D *hDeltaEndYnCap = new TH1D("hDeltaEndYnCap", "#Delta Y_{f} Neutron Capture", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Nuclear Capture at Rest //////////////////////////////////////////
TH1D *hDeltaEndYnuclearCapatureAtRest = new TH1D("hDeltaEndYnuclearCapatureAtRest", "#Delta Y_{f} Nuclear capture at rest ", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Decay //////////////////////////////////////////
TH1D *hDeltaEndYDecay = new TH1D("hDeltaEndYDecay", "#Delta Y_{f} Decay ", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position KaonZero InElastic //////////////////////////////////////////
TH1D *hDeltaEndYKaonZeroInElastic = new TH1D("hDeltaEndYKaonZeroInElastic", "#Delta Y_{f} Kaon Zero InElastic ", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Coulomb Scattering //////////////////////////////////////////
TH1D *hDeltaEndYCoulombScat = new TH1D("hDeltaEndYCoulombScat", "#Delta Y_{f} Coulomb Scattering ", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position MuMinus Capture //////////////////////////////////////////
TH1D *hDeltaEndYMuMinusCapture = new TH1D("hDeltaEndYMuMinusCapture", "#Delta Y_{f} MuMinus Capture ", 200, -50 , 50);
/////////////////////////////////// Delta End Y Position Proton Inelastic //////////////////////////////////////////
TH1D *hDeltaEndYProtonInelastic = new TH1D("hDeltaEndYProtonInelastic", "#Delta Y_{f} Proton Inelastic ", 200, -50 , 50);




/////////////////////////////////// Delta End X Position InElastic //////////////////////////////////////////
TH1D *hDeltaEndXInElastic = new TH1D("hDeltaEndXInElastic", "#Delta X_{f} InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Neutron InElastic //////////////////////////////////////////
TH1D *hDeltaEndXNeutronInElastic = new TH1D("hDeltaEndXNeutronInElastic", "#Delta X_{f} Neutron InElastic", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Hadrondic Elastic //////////////////////////////////////////
TH1D *hDeltaEndXHadElastic = new TH1D("hDeltaEndXHadElastic", "#Delta X_{f} Hadronic Elastic", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Neutron Capture //////////////////////////////////////////
TH1D *hDeltaEndXnCap = new TH1D("hDeltaEndXnCap", "#Delta X_{f} Neutron Capture", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Nuclear Capture at Rest //////////////////////////////////////////
TH1D *hDeltaEndXnuclearCapatureAtRest = new TH1D("hDeltaEndXnuclearCapatureAtRest", "#Delta X_{f} Nuclear capture at rest ", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Decay //////////////////////////////////////////
TH1D *hDeltaEndXDecay = new TH1D("hDeltaEndXDecay", "#Delta X_{f} Decay ", 200, -50 , 50);
/////////////////////////////////// Delta End X Position KaonZero InElastic //////////////////////////////////////////
TH1D *hDeltaEndXKaonZeroInElastic = new TH1D("hDeltaEndXKaonZeroInElastic", "#Delta X_{f} Kaon Zero InElastic ", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Coulomb Scattering //////////////////////////////////////////
TH1D *hDeltaEndXCoulombScat = new TH1D("hDeltaEndXCoulombScat", "#Delta X_{f} Coulomb Scattering ", 200, -50 , 50);
/////////////////////////////////// Delta End X Position MuMinus Capture //////////////////////////////////////////
TH1D *hDeltaEndXMuMinusCapture = new TH1D("hDeltaEndXMuMinusCapture", "#Delta X_{f} MuMinus Capture ", 200, -50 , 50);
/////////////////////////////////// Delta End X Position Proton Inelastic //////////////////////////////////////////
TH1D *hDeltaEndXProtonInelastic = new TH1D("hDeltaEndXProtonInelastic", "#Delta X_{f} Proton Inelastic ", 200, -50 , 50);

// ### Delta Z Stacked  ###
THStack *hStackDeltaEndZ = new THStack("hStackDeltaEndZ", "#Delta Z_{f}");

// ### Delta Y Stacked  ###
THStack *hStackDeltaEndY = new THStack("hStackDeltaEndY", "#Delta Y_{f}");

// ### Delta X Stacked  ###
THStack *hStackDeltaEndX = new THStack("hStackDeltaEndX", "#Delta X_{f}");


void Proton::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Proton.C
//      root> Proton t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
// 					  Putting Flexible Cuts here
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------


// ### Putting in the particle mass that is being simulated here ###
// ###    which is used when calculating the energy loss before  ###
// ###                       entering the TPC                    ###

float particle_mass = 139.57; //<---Mass of Proton in MeV


// ### Number of centimeters in Z we require a track ###
// ### to have a space point within (default = 2 cm) ###
double FirstSpacePointZPos = 2.0;


// ########################################################
// ### Delta X Between Wire Chamber Track and TPC Track ###
// ########################################################
double DeltaXLowerBound = -2.0;
double DeltaXUpperBound = 6.0;

// ########################################################
// ### Delta Y Between Wire Chamber Track and TPC Track ###
// ########################################################
double DeltaYLowerBound = -3.0;
double DeltaYUpperBound = 6.0;


// ########################################################################
// ### Fiducial Boundry Cuts (used to determine if a track is stopping) ###
// ########################################################################
double XLowerFid = 0;
double XUpperFid = 47;

double YLowerFid = -20;
double YUpperFid = 20;

double ZLowerFid = 0;
double ZUpperFid = 90;


// ########################################################################
// ### Definition of the upstream part of the TPC where we restrict the ###
// ###             number of tracks which can be present                ###
// ########################################################################
int UpperPartOfTPC = 14.0;

// #####################################################
// ### Number of tracks allowed in the upstream part ###
// #####################################################
int nLowZTracksAllowed = 4;


// #################################################################################
// ### Making shower Cut (ShortTkLength) and the number of short tracks we allow ###
// #################################################################################
double ShortTkLength = 5.0;
int nShortTracksAllowed = 2;


// ############################
// ### Alpha Cut in degrees ###
// ############################
double alphaCut = 10;

// ### Setting the global event weight based on ###
// ###   open box WCTrack momentum spectrum     ###
   
// 100 MeV < P < 200 MeV = 0.02
// 200 MeV < P < 300 MeV = 0.10
// 300 MeV < P < 400 MeV = 0.535
// 400 MeV < P < 500 MeV = 0.84
// 500 MeV < P < 600 MeV = 0.965
// 600 MeV < P < 700 MeV = 1.0
// 700 MeV < P < 800 MeV = 0.62
// 800 MeV < P < 900 MeV = 0.225
// 900 MeV < P < 1000MeV = 0.094
// 1000MeV < P < 1100MeV = 0.0275
// 1100MeV < P < 1500MeV = 0.01
   
float EventWeight = 1.0;

// #################################################   
// ### True  = Use the momentum based weighting  ###
// ### False = Don't weight events               ###
// #################################################
bool UseEventWeight = true;




// ----------------------------------------------------------------
// Create the cross section from the incident and interaction plots
// ----------------------------------------------------------------
float rho = 1400; //kg/m^3
//  float cm_per_m = 100;
float molar_mass = 39.9; //g/mol
float g_per_kg = 1000; 
float avogadro = 6.02e+23; //number/mol
float number_density = rho*g_per_kg/molar_mass*avogadro;
float slab_width = 0.0045;//in m

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;


// ##########################################################
// ### Putting in some counters for event reduction table ###
// ##########################################################
int nTotalEvents = 0, nEvtsGoodMC = 0, nEvtsTrackZPos = 0, nEvtsMCTrackMatch = 0;
int nEvtsTOF = 0, nEvtsPID = 0, nEventsPassingAlpha = 0, nLowZTrkEvents = 0;
int nNonShowerEvents = 0;

int counter = 0;

// ###############################
// ### Looping over all Events ###
// ###############################
//for (Long64_t jentry=0; jentry<nentries;jentry++)
for (Long64_t jentry=0; jentry<5000;jentry++)
   {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   // if (Cut(ientry) < 0) continue;
   
   // #############################
   // ### Counting Total Events ###
   // #############################
   nTotalEvents++;
   
   
   // === Outputting every nEvents to the screen ===
   if(nTotalEvents % 1000 == 0){std::cout<<"Event = "<<nTotalEvents<<std::endl;} 
   

   //=======================================================================================================================
   //					      GEANT 4 Information
   //=======================================================================================================================
   // ##################################################
   // ### Defining some useful variables we will use ###
   // ##################################################
   int nG4Primary = 0;
   int nG4TrajPoints = 0;
   
   float g4Primary_X0[100] = {0.}, g4Primary_Y0[100] = {0.}, g4Primary_Z0[100] = {0.};
   float g4Primary_ProjX0[100] = {0.}, g4Primary_ProjY0[100] = {0.}, g4Primary_ProjZ0[100] = {0.};
   
   float g4Primary_Xf[100] = {0.}, g4Primary_Yf[100] = {0.}, g4Primary_Zf[100] = {0.};
   float g4Primary_Px[100] = {0.}, g4Primary_Py[100] = {0.}, g4Primary_Pz[100] = {0.};
   
   int g4Primary_TrkID[100] = {999}, g4PrimaryProcess[100] = {0};
   
   int nG4PriTrj = 0;
   float g4Primary_TrueTrjX[10][10000] = {0.}, g4Primary_TrueTrjY[10][10000] = {0.}, g4Primary_TrueTrjZ[10][10000] = {0.};
   float g4Primary_TrueTrjPx[10][10000] = {0.}, g4Primary_TrueTrjPy[10][10000] = {0.}, g4Primary_TrueTrjPz[10][10000] = {0.};
   
   
   float TrueLength = 0;
   float RecoLength = 0;
   // ##########################################
   // ### Loop over all the GEANT4 particles ###
   // ##########################################
   for (int iG4 = 0; iG4 < geant_list_size; iG4++)
      {
      // #####################################################
      // ### If this is a primary particle then look at it ###
      // #####################################################
      if(process_primary[iG4] == 1)
         {

	 
	 // ### Recording information for use later ###
	 g4Primary_X0[nG4Primary] = StartPointx[iG4];
	 g4Primary_Y0[nG4Primary] = StartPointy[iG4];
	 g4Primary_Z0[nG4Primary] = StartPointz[iG4];
	 
	 g4Primary_Xf[nG4Primary] = EndPointx[iG4];
	 g4Primary_Yf[nG4Primary] = EndPointy[iG4];
	 g4Primary_Zf[nG4Primary] = EndPointz[iG4];
	 
	 g4Primary_Px[nG4Primary] = Px[iG4] * 1000; //<---Converting to MeV
	 g4Primary_Py[nG4Primary] = Py[iG4] * 1000; //<---Converting to MeV
	 g4Primary_Pz[nG4Primary] = Pz[iG4] * 1000; //<---Converting to MeV
	 
	 // ### Setting a global event weight ###
	 // 100 MeV < P < 200 MeV = 0.02
   	 // 200 MeV < P < 300 MeV = 0.10
   	 // 300 MeV < P < 400 MeV = 0.535
   	 // 400 MeV < P < 500 MeV = 0.84
   	 // 500 MeV < P < 600 MeV = 0.965
   	 // 600 MeV < P < 700 MeV = 1.0
   	 // 700 MeV < P < 800 MeV = 0.62
   	 // 800 MeV < P < 900 MeV = 0.225
   	 // 900 MeV < P < 1000MeV = 0.094
   	 // 1000MeV < P < 1100MeV = 0.0275
   	 // 1100MeV < P < 1500MeV = 0.01
	 
	 // ### Setting Event weight ### 
	 if(UseEventWeight)
	    {
	    if(g4Primary_Pz[nG4Primary] > 0   && g4Primary_Pz[nG4Primary] < 100){EventWeight = 0.010;}
	    if(g4Primary_Pz[nG4Primary] > 100 && g4Primary_Pz[nG4Primary] < 200){EventWeight = 0.020;}
	    if(g4Primary_Pz[nG4Primary] > 200 && g4Primary_Pz[nG4Primary] < 300){EventWeight = 0.100;}
	    if(g4Primary_Pz[nG4Primary] > 300 && g4Primary_Pz[nG4Primary] < 400){EventWeight = 0.535;}
	    if(g4Primary_Pz[nG4Primary] > 400 && g4Primary_Pz[nG4Primary] < 500){EventWeight = 0.840;}
	    if(g4Primary_Pz[nG4Primary] > 500 && g4Primary_Pz[nG4Primary] < 600){EventWeight = 0.965;}
	    if(g4Primary_Pz[nG4Primary] > 600 && g4Primary_Pz[nG4Primary] < 700){EventWeight = 1.000;}
	    if(g4Primary_Pz[nG4Primary] > 700 && g4Primary_Pz[nG4Primary] < 800){EventWeight = 0.620;}
	    if(g4Primary_Pz[nG4Primary] > 800 && g4Primary_Pz[nG4Primary] < 900){EventWeight = 0.225;}
	    if(g4Primary_Pz[nG4Primary] > 900 && g4Primary_Pz[nG4Primary] <1000){EventWeight = 0.094;}
	    if(g4Primary_Pz[nG4Primary] >1000 && g4Primary_Pz[nG4Primary] <1100){EventWeight = 0.0275;}
	    if(g4Primary_Pz[nG4Primary] >1100){EventWeight = 0.010;}
	    
	    }//<---End Assigning Event weight
	 
	 TrueLength = sqrt( ((EndPointz[iG4]-StartPointz[iG4])*(EndPointz[iG4]-StartPointz[iG4])) + 
	                    ((EndPointy[iG4]-StartPointy[iG4])*(EndPointy[iG4]-StartPointy[iG4])) + 
	                    ((EndPointx[iG4]-StartPointx[iG4])*(EndPointx[iG4]-StartPointx[iG4])) );
			    
	 hTrueLength->Fill(TrueLength);
	 
	 // ------------------------------------------------------------------------------------
	 // ---------------        Extrapolate the X, Y, Z position of the primary         -----
	 // ---------------     if it started upstream of the front face of the TPC        -----
	 // ------------------------------------------------------------------------------------
	 
	 double b1 = StartPointz[iG4] - StartPointx[iG4]*Pz[iG4]/Px[iG4];
	 double b2 = StartPointz[iG4] - StartPointy[iG4]*Pz[iG4]/Py[iG4];
	 
	 g4Primary_ProjX0[nG4Primary] = -b1*Px[iG4]/Pz[iG4];
	 g4Primary_ProjY0[nG4Primary] = -b2*Py[iG4]/Pz[iG4];
	 g4Primary_ProjZ0[nG4Primary] = 0.0;
	 
	 // ### Setting the primary particles Track ID ###
	 g4Primary_TrkID[nG4Primary] = TrackId[iG4];
	 
	 
	 hMCPrimaryEndXvsZ->Fill(EndPointz[iG4], EndPointx[iG4]);
	 hMCPrimaryEndYvsZ->Fill(EndPointz[iG4], EndPointy[iG4]);
	 
	 // ##############################
	 // ### Filling the histograms ###
	 // ##############################
	 hMCPrimaryPx->Fill(g4Primary_Px[nG4Primary], EventWeight);
	 hMCPrimaryPy->Fill(g4Primary_Py[nG4Primary], EventWeight);
	 hMCPrimaryPz->Fill(g4Primary_Pz[nG4Primary], EventWeight);

	 hMCPrimaryStartX->Fill(StartPointx[iG4]);
	 hMCPrimaryStartY->Fill(StartPointy[iG4]);
	 hMCPrimaryStartZ->Fill(StartPointz[iG4]);
	 
	 hMCPrimaryEndX->Fill(EndPointx[iG4]);
	 hMCPrimaryEndY->Fill(EndPointy[iG4]);
	 hMCPrimaryEndZ->Fill(EndPointz[iG4]);

	 hMCPrimaryProjectedStartX->Fill(g4Primary_ProjX0[nG4Primary]);
	 hMCPrimaryProjectedStartY->Fill(g4Primary_ProjY0[nG4Primary]);
	 hMCPrimaryProjectedStartZ->Fill(g4Primary_ProjZ0[nG4Primary]);
	 
	 nG4TrajPoints = 0;
	 
	 // ### Recording the primary particles trajectory points ###
	 for(int iG4Tr = 0; iG4Tr < NTrTrajPts[iG4]; iG4Tr++)
	    {
	    g4Primary_TrueTrjX[nG4Primary][nG4PriTrj] = MidPosX[iG4][iG4Tr];
	    g4Primary_TrueTrjY[nG4Primary][nG4PriTrj] = MidPosY[iG4][iG4Tr];
	    g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj] = MidPosZ[iG4][iG4Tr];
	    
	    //std::cout<<"g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj] = "<<g4Primary_TrueTrjZ[nG4Primary][nG4PriTrj]<<std::endl;
	    
	    g4Primary_TrueTrjPx[nG4Primary][nG4PriTrj] = MidPx[iG4][iG4Tr]*1000;//<---Converting to MeV
	    g4Primary_TrueTrjPy[nG4Primary][nG4PriTrj] = MidPy[iG4][iG4Tr]*1000;//<---Converting to MeV
	    g4Primary_TrueTrjPz[nG4Primary][nG4PriTrj] = MidPz[iG4][iG4Tr]*1000;//<---Converting to MeV
	    
	    nG4PriTrj++;
	    }//<---end looping over this primary particles true trajectory points
	 
	 
	 // ### Bumping the counters ###
	 nG4Primary++;
	 
	 }//<---End Looking only at primaries

      }// <---End iG4 Loop
   
   // ################################################
   // ### Loop over all the GEANT4 particles again ###
   // ###  to get the process from the daughters   ###
   // ################################################
   for (int iG4 = 0; iG4 < geant_list_size; iG4++)
      {
      
      // ### Looking for the Daughters of the primary ###
      if(Mother[iG4] == g4Primary_TrkID[nG4Primary - 1])
	 {
	 g4PrimaryProcess[nG4Primary - 1] = Process[iG4];
	 
	 
	 
	 
	 }//<---End matching daughters
      
      }//<---end iG4 loop
   hMCPrimaryProcess->Fill(g4PrimaryProcess[nG4Primary - 1]);
   
   
   
   //=======================================================================================================================
   //				Only looking at events where the primary particle enters the TPC
   //=======================================================================================================================
   
   bool GoodMCEventInTPC = true;
   
   // ##############################################
   // ### Looping over all the primary particles ###
   // ##############################################
   for(int npri = 0; npri < nG4Primary; npri++)
      {
      if(g4Primary_Zf[npri] < 0){GoodMCEventInTPC = false;}
      
      // ####################################################################
      // ### Calculating the energy loss for particles that enter the TPC ###
      // ####################################################################
      if(GoodMCEventInTPC)
         {
	 float DifferenceInEnergy = 0;
	 // ### Loop over the true trajectory points ###
	 for(int ntrj = 0; ntrj < nG4PriTrj; ntrj++)
	    {
\
	    // ### Only looking at point which are upstream of the TPC ###
	    if(g4Primary_TrueTrjZ[npri][ntrj] < 0)
	       {
	       
	       float Momentum_Point1 = sqrt((g4Primary_TrueTrjPx[npri][ntrj]*g4Primary_TrueTrjPx[npri][ntrj]) + 
	                               (g4Primary_TrueTrjPy[npri][ntrj]*g4Primary_TrueTrjPy[npri][ntrj]) +
				       (g4Primary_TrueTrjPz[npri][ntrj]*g4Primary_TrueTrjPz[npri][ntrj]));
				       
	       float Momentum_Point2 = sqrt((g4Primary_TrueTrjPx[npri][ntrj+1]*g4Primary_TrueTrjPx[npri][ntrj+1]) + 
	                               (g4Primary_TrueTrjPy[npri][ntrj+1]*g4Primary_TrueTrjPy[npri][ntrj+1]) +
				       (g4Primary_TrueTrjPz[npri][ntrj+1]*g4Primary_TrueTrjPz[npri][ntrj+1]));
				       
	       float Energy_Point1 = sqrt( (Momentum_Point1*Momentum_Point1) + (particle_mass*particle_mass)  );
	       
	       float Energy_Point2 = sqrt( (Momentum_Point2*Momentum_Point2) + (particle_mass*particle_mass)  );
	       
	       DifferenceInEnergy +=  Energy_Point1 - Energy_Point2;
	       
	       std::cout<<"z = "<<g4Primary_TrueTrjZ[npri][ntrj]<<", DifferenceInEnergy = "<<DifferenceInEnergy<<std::endl;
	       
	       
	       }//<---End only look at points which are upstream of the TPC
	    
	    
	    }//<---End ntrj for loop
	 
	 
	 hMCELossUpstream->Fill(DifferenceInEnergy);
	 }//<---Only looking at events that actually make it into the TPC
      
      
      }//<---End npri loop
   
   if(!GoodMCEventInTPC){continue;}
   nEvtsGoodMC++;


   
   
   //=======================================================================================================================
   //						Low Z Spacepoint Track Cut
   //=======================================================================================================================
   
   // ### Boolian for events w/ track which ###
   // ###     starts at the front face      ###
   bool TrackTrjPtsZCut = false;
   
   // ### Recording the index of the track which ###
   // ###   starts at the front face of the TPC  ###
   bool PreLimTrackIndex[500] = {false};
   
   // ##################################################
   // ### Defining a dummy variable used for sorting ###
   // ##################################################
   double dummyXpoint = 999, dummyYpoint = 999, dummyZpoint = 100, dummyPointTrkInd = -1;
   double dummypoint_TempTrjX = 999, dummypoint_TempTrjY = 999, dummypoint_TempTrjZ = 999;
   double dummyTrkX[200] = {0.}, dummyTrkY[200] = {0.}, dummyTrkZ[200] = {100.};
   double dummyTrk_pHat0X[200] = {0.}, dummyTrk_pHat0Y[200] = {0.}, dummyTrk_pHat0Z[200] = {0.};
   double dummyTrk_Theta[200] = {0.}, dummyTrk_Phi[200] = {0.};
   double dummyTrk_Index[200] = {0.};
   
   int nUpStreamTrk = 0;
   
   // ###########################
   // ### Looping over tracks ###
   // ###########################
   for(int iTrk = 0; iTrk < ntracks_reco; iTrk++)
      {
      
      float tempZpoint = 100;
      // ########################################################
      // ### Looping over the trajectory points for the track ###
      // ########################################################
      for(int iTrjPt = 0; iTrjPt < nTrajPoint[iTrk]; iTrjPt++)
         {
	 
	 // ################################################################################ 
         // ### Tracking the lowest Z point that is inside fiducial boundries of the TPC ###
	 // ################################################################################
	 // ### Resetting the variables for each track ###
         dummyXpoint = 999, dummyYpoint = 999, dummyZpoint = 999;
	 
	 // ###########################################################################
	 // ### Setting our dummypoints if this is the lowest Z point on this track ###
	 // ###           and still within the active volume of the TPC             ###
	 // ###########################################################################
	 if(trjPt_Z[iTrk][iTrjPt] < tempZpoint && trjPt_Z[iTrk][iTrjPt] > ZLowerFid && //<---Note: the variable "trjPt_Z[iTrk][iTrjPt]"
	    trjPt_Y[iTrk][iTrjPt] > YLowerFid && trjPt_Y[iTrk][iTrjPt] < YUpperFid &&       // is the z position of the 3d point for the iTrk
	    trjPt_X[iTrk][iTrjPt] > XLowerFid && trjPt_X[iTrk][iTrjPt] < XUpperFid )
	     
	    {
	    dummyXpoint = trjPt_X[iTrk][iTrjPt];
	    dummyYpoint = trjPt_Y[iTrk][iTrjPt];
	    dummyZpoint = trjPt_Z[iTrk][iTrjPt];
	    
	    dummypoint_TempTrjX = pHat0_X[iTrk][iTrjPt];
	    dummypoint_TempTrjY = pHat0_Y[iTrk][iTrjPt];
	    dummypoint_TempTrjZ = pHat0_Z[iTrk][iTrjPt];
	    
	    dummyPointTrkInd = iTrk;
	    
	    tempZpoint = trjPt_Z[iTrk][iTrjPt];
	    }//<---End looking for the most upstream point
        
	 // ### Only passing events with a track that has ###
	 // ###  a spacepoint within the first N cm in Z  ### 
	 // ###    And requiring it to be inside the TPC  ###
	 if(dummyZpoint < FirstSpacePointZPos)
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
	    }
	 }//<---End looping over nspts	
	 
      // ### Filling the most upstream spacepoint for this track ###
      hdataUpstreamZPos->Fill(tempZpoint);
      
      // ### Recording that this track is a "good Track if ###
      // ###  it has a space point in the first N cm in Z  ###
      if(TrackTrjPtsZCut){ PreLimTrackIndex[iTrk] = true;}
      	 
      }//<---End nTPCtrk loop
      
   // ###############################################
   // ### Skipping events that don't have a track ###
   // ###   in the front of the TPC (Z) Position  ###
   // ###############################################
   if(!TrackTrjPtsZCut){continue;}
   // ### Counting Events w/ front face TPC Track ###
   nEvtsTrackZPos++;
   
   // ###############################################################
   // ### Making plot of number of tracks as a function of length ###
   // ###############################################################
   LowZCut();
   
   
   //=======================================================================================================================
   //					Cutting on the number of tracks in the upstream TPC
   //=======================================================================================================================
   
   int nLowZTracksInTPC = 0;
   // ################################################################
   // ### Initializing variables for study of low Z track location ###
   // ################################################################
   bool LowZTrackInTPC = false;
   
   float templowz1 = 0;
   float templowz2 = 0;
   
   
    // #################################################################
    // ### Only keeping events if there is less than N tracks in the ###
    // ###    first ## cm of the TPC (to help cut out EM Showers     ###
    // #################################################################
    for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      // ### Start by assuming this track is not in the ###
      // ###          low Z part of the TPC             ###
      LowZTrackInTPC = false;
           
      // ##################################################
      // ### Looping over the spacepoints for the track ###
      // ##################################################
      for(size_t nspts = 0; nspts < nTrajPoint[nTPCtrk]; nspts++)
         {
	 
	 // ##################################################
	 // ### Count this track if it has a spacepoint in ###
	 // ###       the low Z region of the TPC          ###
	 // ##################################################
	 if(trkz[nTPCtrk][nspts] < UpperPartOfTPC)
	    {  
	    if(trky[nTPCtrk][nspts] > YLowerFid && trky[nTPCtrk][nspts] < YUpperFid && 
	       trkx[nTPCtrk][nspts] > XLowerFid && trkx[nTPCtrk][nspts] < XUpperFid)
	        {LowZTrackInTPC = true; templowz1 = trkz[nTPCtrk][nspts];}
		
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
    if(nLowZTracksInTPC > nLowZTracksAllowed || nLowZTracksInTPC == 0){continue;}  
    
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
   
   // === Calculating Phi for MC ===
   //---------------------------------------------------------------------------------------------------------------------
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
   //---------------------------------------------------------------------------------------------------------------------
   
   
   // ######################################################
   // ### Calculating the angles for the Upstream Tracks ###
   // ######################################################
   
   TVector3 z_hat_TPC(0,0,1);
   TVector3 p_hat_0_TPC;
   for(int aa = 0; aa < nUpStreamTrk; aa++)
      {
      // ### Setting the TVector ###
      p_hat_0_TPC.SetX(dummyTrk_pHat0X[aa]);
      p_hat_0_TPC.SetY(dummyTrk_pHat0Y[aa]);
      p_hat_0_TPC.SetZ(dummyTrk_pHat0Z[aa]);
      
      
      // ### Calculating TPC track theta ###
      dummyTrk_Theta[aa] = acos(z_hat_TPC.Dot(p_hat_0_TPC)/p_hat_0_TPC.Mag());
      
      // ### Calculating TPC track phi ###
      //---------------------------------------------------------------------------------------------------------------------
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
      //---------------------------------------------------------------------------------------------------------------------
      
      }//<---End aa loop
   
   // =============================================================================================================
   // 				Cutting on the DeltaX and Delta Y between MC and Reco Track
   // =============================================================================================================
   
   bool MC_TPCMatch = false;
   int nDeltaMatch = 0;
   // ###################################
   // ### Loop over all the US Tracks ###
   // ###################################
   for(int bb = 0; bb < nUpStreamTrk; bb++)
      {
      float DeltaX = dummyTrkX[bb] - g4Primary_ProjX0[0];
      float DeltaY = dummyTrkY[bb] - g4Primary_ProjY0[0];
      float DeltaZ = dummyTrkZ[bb] - g4Primary_ProjZ0[0];
      
      hDeltaX->Fill(DeltaX);
      hDeltaY->Fill(DeltaY);
      hDeltaZ->Fill(DeltaZ);
      
      // ### Matching in Delta X and Delta Y and alpha ###
      if(DeltaX > DeltaXLowerBound && DeltaX < DeltaXUpperBound && DeltaY > DeltaYLowerBound && DeltaY < DeltaYUpperBound)
         {
	 MC_TPCMatch = true;
	 nDeltaMatch++;
	 }//<---End matching
      
      }//<---End bb loop
   
   if(!MC_TPCMatch || nDeltaMatch < 1 || nDeltaMatch > 1){continue;}
   nEvtsMCTrackMatch++;
   
   // =============================================================================================================
   // 				Cutting on the Alpha Angle between MC and Reco Track
   // =============================================================================================================
   
   bool AlphaMatch = false;
   int RecoTrackIndex = -1;
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
      
      hAlpha->Fill(alpha);
      
      // ### Setting the boolian for the true match ###
      if(alpha < alphaCut)
         {
	 AlphaMatch = true;
	 RecoTrackIndex = dummyTrk_Index[bb];
	 }
      
      }//<---End bb loop
   
   // #############################################################
   // ### Skipping the event if there is no good match in alpha ###
   // #############################################################
   if(!AlphaMatch){continue;}
   nEventsPassingAlpha++;
   
   
   
   // ===========================================================================================================================================
   // 						Vetoing Shower Like Events 
   // ===========================================================================================================================================   
   
   int nShortTrks = 0;
   // ###############################
   // ### Looping over TPC tracks ###
   // ###############################
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      
      // ### If the track is shorter than our cut ###
      if(trklength[nTPCtrk] < ShortTkLength)
         {
	 nShortTrks++;
	 
	 }
	 
      }//<---End nTPCtrk
   
   // ### Skipping the event if there are too many short tracks ###
   if(nShortTrks > nShortTracksAllowed){continue;}
   
   // ### Bumping the counter ###
   nNonShowerEvents++;
   
   
   
   // ===========================================================================================================================================
   // 						Calculating the cross-section 
   // ===========================================================================================================================================   
   
   // ### The assumed energy loss between the cryostat and the TPC ###
   float entryTPCEnergyLoss = 8.6; //MeV

   // ### The assumed mass of the incident particle (here we assume a Proton) ###
   float mass = 139.57;
   
   // #################################################
   // ### Setting a flag to exclude exiting tracks  ###
   // ###  from the numerator of the cross-section  ###
   // #################################################
   
   bool ExitingTrack = false;
   
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
   
   // ### Filling the initial kinetic energy plot ###
   hMCInitalKE->Fill(kineticEnergy);
   
   
   // ########################################################################
   // ### Variables for the track we are calculating the cross-section for ###
   // ########################################################################
   double Protondedx[1000]={0.};
   double Protonresrange[1000]={0.};
   double Protonpitchhit[1000]={0.};
   int nProtonSpts = 0;
   double ProtonSumEnergy = 0;
   
   // ### Variables for determining the matching of the end point ###
   float TrackEndX = 999, TrackEndY = 999, TrackEndZ = 999;
   
   // ################################
   // ### Loop over all TPC Tracks ###
   // ################################
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      if(nTPCtrk != RecoTrackIndex){continue;}
      
      // ### Recording the end-point of this track ###
      TrackEndX = trkendx[nTPCtrk];
      TrackEndY = trkendy[nTPCtrk];
      TrackEndZ = trkendz[nTPCtrk];
      
      hdataProtonTrackEndX->Fill(TrackEndX);
      hdataProtonTrackEndY->Fill(TrackEndY);
      hdataProtonTrackEndZ->Fill(TrackEndZ);
      
      // ### Recording the start-point of this track ###
      
      hdataProtonTrackStartX->Fill(trkvtxx[nTPCtrk]);
      hdataProtonTrackStartY->Fill(trkvtxy[nTPCtrk]);
      hdataProtonTrackStartZ->Fill(trkvtxz[nTPCtrk]);
      
      RecoLength = sqrt( ((trkendz[nTPCtrk]-trkvtxz[nTPCtrk])*(trkendz[nTPCtrk]-trkvtxz[nTPCtrk])) + 
	                    ((trkendy[nTPCtrk]-trkvtxy[nTPCtrk])*(trkendy[nTPCtrk]-trkvtxy[nTPCtrk])) + 
	                    ((trkendx[nTPCtrk]-trkvtxx[nTPCtrk])*(trkendx[nTPCtrk]-trkvtxx[nTPCtrk])) );
      
      hRecoLength->Fill(RecoLength);
      
      // #################################################
      // ### If this tracks end point is at a boundary ###
      // ###   then tag this track as "through-going"  ###
      // #################################################
      if(TrackEndX < 1   || TrackEndX > 42.0 || TrackEndY > 19 ||
         TrackEndY < -19 || TrackEndZ > 89.0)
	 {ExitingTrack = true;}
      
      nProtonSpts = 0;
      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
      for(size_t nspts = 0; nspts < ntrkhits[nTPCtrk]; nspts++)
         {
	 // ###                 Note: Format for this variable is:             ###
	 // ### [trk number][plane 0 = induction, 1 = collection][spts number] ###
         Protondedx[nProtonSpts]     = trkdedx[nTPCtrk][1][nspts];
	 
	 // ### Putting in a fix in the case that the dE/dX is negative in this step ### 
	 // ###  then take the point before and the point after and average them
	 if(Protondedx[nProtonSpts] < 0 && nspts < ntrkhits[nTPCtrk] && nspts > 0)
	    {Protondedx[nProtonSpts] = ( (trkdedx[nTPCtrk][1][nspts - 1] + trkdedx[nTPCtrk][1][nspts + 1]) / 2);}
	 
	 // ### If this didn't fix it, then just put in a flat 2.4 MeV / cm fix ###
	 if(Protondedx[nProtonSpts] < 0)
	    {
	    Protondedx[nProtonSpts] = 2.4;
	    continue;
	    }
	 
         Protonresrange[nProtonSpts] = trkrr[nTPCtrk][1][nspts];
         Protonpitchhit[nProtonSpts] = trkpitchhit[nTPCtrk][1][nspts];
         
	 ProtonSumEnergy = (Protondedx[nProtonSpts] * Protonpitchhit[nProtonSpts]) + ProtonSumEnergy;
	 
	 // ### Recording the dE/dX ###
	 hdataProtondEdX->Fill(Protondedx[nProtonSpts]);
	 // ### Recording the residual range ###
	 hdataProtonRR->Fill(Protonresrange[nProtonSpts]);
	 // ### Recording the Pitch ###
	 hdataProtonTrkPitch->Fill(Protonpitchhit[nProtonSpts]);
	 
	 // ### Filling 2d dE/dX vs RR ###
	 hdataProtondEdXvsRR->Fill(Protonresrange[nProtonSpts], Protondedx[nProtonSpts]);
	 
	 nProtonSpts++;
         }//<---End nspts loop
      
      
      }//<---End nTPCtrk loop
   
   
   // ###################################
   // ### Filling the Delta End Point ###
   // ###################################
   float DeltaEndX = g4Primary_Xf[0] - TrackEndX;
   float DeltaEndY = g4Primary_Yf[0] - TrackEndY;
   float DeltaEndZ = g4Primary_Zf[0] - TrackEndZ;
      
   hDeltaEndX->Fill(DeltaEndX);
   hDeltaEndY->Fill(DeltaEndY);
   hDeltaEndZ->Fill(DeltaEndZ);
   
   if(DeltaEndZ > 3 || DeltaEndZ < 3){counter++;}
   
   if(g4PrimaryProcess[0] == 1){hDeltaEndZInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 2){hDeltaEndZNeutronInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 3){hDeltaEndZHadElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 4){hDeltaEndZnCap->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 5){hDeltaEndZnuclearCapatureAtRest->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 6){hDeltaEndZDecay->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 7){hDeltaEndZKaonZeroInElastic->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 8){hDeltaEndZCoulombScat->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 9){hDeltaEndZMuMinusCapture->Fill(DeltaEndZ);}
   if(g4PrimaryProcess[0] == 10){hDeltaEndZProtonInelastic->Fill(DeltaEndZ);}
   
   if(g4PrimaryProcess[0] == 1){hDeltaEndYInElastic->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 2){hDeltaEndYNeutronInElastic->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 3){hDeltaEndYHadElastic->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 4){hDeltaEndYnCap->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 5){hDeltaEndYnuclearCapatureAtRest->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 6){hDeltaEndYDecay->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 7){hDeltaEndYKaonZeroInElastic->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 8){hDeltaEndYCoulombScat->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 9){hDeltaEndYMuMinusCapture->Fill(DeltaEndY);}
   if(g4PrimaryProcess[0] == 10){hDeltaEndYProtonInelastic->Fill(DeltaEndY);}
   
   if(g4PrimaryProcess[0] == 1){hDeltaEndXInElastic->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 2){hDeltaEndXNeutronInElastic->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 3){hDeltaEndXHadElastic->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 4){hDeltaEndXnCap->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 5){hDeltaEndXnuclearCapatureAtRest->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 6){hDeltaEndXDecay->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 7){hDeltaEndXKaonZeroInElastic->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 8){hDeltaEndXCoulombScat->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 9){hDeltaEndXMuMinusCapture->Fill(DeltaEndX);}
   if(g4PrimaryProcess[0] == 10){hDeltaEndXProtonInelastic->Fill(DeltaEndX);}
   
   // ###########################################################
   // ### Looping over the spacepoints to fill the histograms ###
   // ###########################################################
   for(size_t npoints = 0; npoints < nProtonSpts; npoints++)
      {
      // ### Filling the incidient histogram ###
      hdataProtonIncidentKE->Fill(kineticEnergy, EventWeight);
      
      // ### Filling the interaction histogram for the last spt ###
      if(npoints == nProtonSpts -1 && !ExitingTrack)
         {fProtonInteractions->Fill(kineticEnergy, EventWeight);}
      
      //float energyLossInStep = Protondedx[npoints] * Protonresrange[npoints] * RecombinationFactor;
      float energyLossInStep = Protondedx[npoints] * Protonpitchhit[npoints];
      
      kineticEnergy -= energyLossInStep;
      
      
      }//<---End npoints loop
    
   
   }//<---End jentry loop
//-------------------------------------------------------------------------------------------------------------------------------------------------





//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//						Coloring the histograms
hDeltaEndZInElastic->SetFillColor(kYellow+2);
hDeltaEndZInElastic->SetLineColor(kYellow+2);
hDeltaEndZInElastic->SetFillStyle(3006);
hStackDeltaEndZ->Add(hDeltaEndZInElastic);


hDeltaEndZNeutronInElastic->SetFillColor(kGreen+2);
hDeltaEndZNeutronInElastic->SetLineColor(kGreen+2);
hDeltaEndZNeutronInElastic->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZNeutronInElastic);

hDeltaEndZHadElastic->SetFillColor(kCyan+3);
hDeltaEndZHadElastic->SetLineColor(kCyan+3);
hDeltaEndZHadElastic->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZHadElastic);

hDeltaEndZnCap->SetFillColor(kBlue);
hDeltaEndZnCap->SetLineColor(kBlue);
hDeltaEndZnCap->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZnCap);

hDeltaEndZnuclearCapatureAtRest->SetFillColor(kMagenta);
hDeltaEndZnuclearCapatureAtRest->SetLineColor(kMagenta);
hDeltaEndZnuclearCapatureAtRest->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZnuclearCapatureAtRest);

hDeltaEndZDecay->SetFillColor(kRed+2);
hDeltaEndZDecay->SetLineColor(kRed+2);
hDeltaEndZDecay->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZDecay);

hDeltaEndZKaonZeroInElastic->SetFillColor(kOrange+7);
hDeltaEndZKaonZeroInElastic->SetLineColor(kOrange+7);
hDeltaEndZKaonZeroInElastic->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZKaonZeroInElastic);

hDeltaEndZCoulombScat->SetFillColor(kAzure+1);
hDeltaEndZCoulombScat->SetLineColor(kAzure+1);
hDeltaEndZCoulombScat->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZCoulombScat);

hDeltaEndZMuMinusCapture->SetFillColor(kViolet-1);
hDeltaEndZMuMinusCapture->SetLineColor(kViolet-1);
hDeltaEndZMuMinusCapture->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZMuMinusCapture);

hDeltaEndZProtonInelastic->SetFillColor(kBlack);
hDeltaEndZProtonInelastic->SetLineColor(kBlack);
hDeltaEndZProtonInelastic->SetFillStyle(3002);
hStackDeltaEndZ->Add(hDeltaEndZProtonInelastic);

////////////////////////////////////////////////////////////////////////

hDeltaEndYInElastic->SetFillColor(kYellow+2);
hDeltaEndYInElastic->SetLineColor(kYellow+2);
hDeltaEndYInElastic->SetFillStyle(3006);
hStackDeltaEndY->Add(hDeltaEndYInElastic);


hDeltaEndYNeutronInElastic->SetFillColor(kGreen+2);
hDeltaEndYNeutronInElastic->SetLineColor(kGreen+2);
hDeltaEndYNeutronInElastic->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYNeutronInElastic);

hDeltaEndYHadElastic->SetFillColor(kCyan+3);
hDeltaEndYHadElastic->SetLineColor(kCyan+3);
hDeltaEndYHadElastic->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYHadElastic);

hDeltaEndYnCap->SetFillColor(kBlue);
hDeltaEndYnCap->SetLineColor(kBlue);
hDeltaEndYnCap->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYnCap);

hDeltaEndYnuclearCapatureAtRest->SetFillColor(kMagenta);
hDeltaEndYnuclearCapatureAtRest->SetLineColor(kMagenta);
hDeltaEndYnuclearCapatureAtRest->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYnuclearCapatureAtRest);

hDeltaEndYDecay->SetFillColor(kRed+2);
hDeltaEndYDecay->SetLineColor(kRed+2);
hDeltaEndYDecay->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYDecay);

hDeltaEndYKaonZeroInElastic->SetFillColor(kOrange+7);
hDeltaEndYKaonZeroInElastic->SetLineColor(kOrange+7);
hDeltaEndYKaonZeroInElastic->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYKaonZeroInElastic);

hDeltaEndYCoulombScat->SetFillColor(kAzure+1);
hDeltaEndYCoulombScat->SetLineColor(kAzure+1);
hDeltaEndYCoulombScat->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYCoulombScat);

hDeltaEndYMuMinusCapture->SetFillColor(kViolet-1);
hDeltaEndYMuMinusCapture->SetLineColor(kViolet-1);
hDeltaEndYMuMinusCapture->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYMuMinusCapture);

hDeltaEndYProtonInelastic->SetFillColor(kBlack);
hDeltaEndYProtonInelastic->SetLineColor(kBlack);
hDeltaEndYProtonInelastic->SetFillStyle(3002);
hStackDeltaEndY->Add(hDeltaEndYProtonInelastic);


////////////////////////////////////////////////////////////////////////

hDeltaEndXInElastic->SetFillColor(kYellow+2);
hDeltaEndXInElastic->SetLineColor(kYellow+2);
hDeltaEndXInElastic->SetFillStyle(3006);
hStackDeltaEndX->Add(hDeltaEndXInElastic);


hDeltaEndXNeutronInElastic->SetFillColor(kGreen+2);
hDeltaEndXNeutronInElastic->SetLineColor(kGreen+2);
hDeltaEndXNeutronInElastic->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXNeutronInElastic);

hDeltaEndXHadElastic->SetFillColor(kCyan+3);
hDeltaEndXHadElastic->SetLineColor(kCyan+3);
hDeltaEndXHadElastic->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXHadElastic);

hDeltaEndXnCap->SetFillColor(kBlue);
hDeltaEndXnCap->SetLineColor(kBlue);
hDeltaEndXnCap->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXnCap);

hDeltaEndXnuclearCapatureAtRest->SetFillColor(kMagenta);
hDeltaEndXnuclearCapatureAtRest->SetLineColor(kMagenta);
hDeltaEndXnuclearCapatureAtRest->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXnuclearCapatureAtRest);

hDeltaEndXDecay->SetFillColor(kRed+2);
hDeltaEndXDecay->SetLineColor(kRed+2);
hDeltaEndXDecay->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXDecay);

hDeltaEndXKaonZeroInElastic->SetFillColor(kOrange+7);
hDeltaEndXKaonZeroInElastic->SetLineColor(kOrange+7);
hDeltaEndXKaonZeroInElastic->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXKaonZeroInElastic);

hDeltaEndXCoulombScat->SetFillColor(kAzure+1);
hDeltaEndXCoulombScat->SetLineColor(kAzure+1);
hDeltaEndXCoulombScat->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXCoulombScat);

hDeltaEndXMuMinusCapture->SetFillColor(kViolet-1);
hDeltaEndXMuMinusCapture->SetLineColor(kViolet-1);
hDeltaEndXMuMinusCapture->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXMuMinusCapture);

hDeltaEndXProtonInelastic->SetFillColor(kBlack);
hDeltaEndXProtonInelastic->SetLineColor(kBlack);
hDeltaEndXProtonInelastic->SetFillStyle(3002);
hStackDeltaEndX->Add(hDeltaEndXProtonInelastic);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



// ===============================================================================================================
//					MAKING THE CROSS-SECTION PLOT
// ===============================================================================================================


// ###################################################################
// #### Looping over the exiting bins to extract the cross-section ###
// ###################################################################
for( int iBin = 1; iBin <= fProtonInteractions->GetNbinsX(); ++iBin )
   {
   
   std::cout<<std::endl;
   // ### If an incident bin is equal to zero then skip that bin ###
   if( hdataProtonIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
   
   // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = (fProtonInteractions->GetBinContent(iBin)/hdataProtonIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);
   
   //std::cout<<"Cross-Section before conversion to barns = "<<TempCrossSection<<std::endl;
   
   float crossSection = TempCrossSection * (1/1e-28); //To put this into barns
   //std::cout<<"Cross-Section = "<<crossSection<<std::endl;
   
   
   fCrossSection->SetBinContent(iBin,crossSection);
   
   float numError = pow(fProtonInteractions->GetBinContent(iBin),0.5);
   float num = fProtonInteractions->GetBinContent(iBin);

   
   if(num == 0){num = 1;}
   float term1 = numError/num;
   //std::cout<<"term1 = "<<term1<<std::endl;
   
   float denomError = pow(hdataProtonIncidentKE->GetBinContent(iBin),0.5);
   float denom = hdataProtonIncidentKE->GetBinContent(iBin);
   if(denom == 0){denom = 1;}
   float term2 = denomError/denom;
   //std::cout<<"term2 = "<<term2<<std::endl;
   float totalError = (TempCrossSection) * (pow( ( (term1*term1) + (term2*term2) ),0.5)) * (1/number_density) * (1/slab_width) * (1/1e-28) *(1e26);
   //std::cout<<"totalError = "<<totalError<<std::endl;
   fCrossSection->SetBinError(iBin,totalError);
       
   }//<---End iBin loop



// ========================================================================================================
// ===					EVENT REDUCTION TABLE						===
// ========================================================================================================
std::cout<<std::endl;
std::cout<<"###################################################################################"<<std::endl;
std::cout<<"### Number of Events in AnaModule                                         = "<<nTotalEvents<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ Primary MC which enters the TPC                   = "<<nEvtsGoodMC<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ Trk Z < "<<FirstSpacePointZPos<<"                                = "<<nEvtsTrackZPos<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ < "<<nLowZTracksAllowed<<" tracks in the first "<<UpperPartOfTPC<<" cm of the TPC     = "<<nLowZTrkEvents<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ ONE MC Particle Matched                               = "<<nEvtsMCTrackMatch<<" ###"<<std::endl;
std::cout<<"###              ( "<<DeltaXLowerBound<<" < Delta X < "<<DeltaXUpperBound<<" , "<<DeltaYLowerBound<<" < Delta Y < "<<DeltaYUpperBound<<" )             ###"<<std::endl;
std::cout<<"### Number of Events w/ angle alpha less the "<<alphaCut<<" degrees  = "<<nEventsPassingAlpha<<"         ###"<<std::endl;
std::cout<<"### Number of Events that are not Shower Like                        = "<<nNonShowerEvents<<std::endl;
std::cout<<"###################################################################################"<<std::endl;
std::cout<<std::endl;

std::cout<<"counter = "<<counter<<std::endl;



// ===================================================================================================================
// =====================================   Drawing Histograms   ======================================================
// ===================================================================================================================
// ### Making a TCanvas ###
TCanvas *c1= new TCanvas("c1","Initial X Pos");
c1->SetTicks();
c1->SetLogy();
c1->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryStartX->SetLineColor(kRed);
hMCPrimaryStartX->SetLineStyle(0);
hMCPrimaryStartX->SetLineWidth(3);
hMCPrimaryStartX->SetMarkerStyle(8);
hMCPrimaryStartX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryStartX->GetXaxis()->SetTitle("Primary Particle X_{0} (cm)");
hMCPrimaryStartX->GetXaxis()->CenterTitle();

hMCPrimaryStartX->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryStartX->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryStartX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c2= new TCanvas("c2","Initial Y Pos");
c2->SetTicks();
c2->SetLogy();
c2->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryStartY->SetLineColor(kRed);
hMCPrimaryStartY->SetLineStyle(0);
hMCPrimaryStartY->SetLineWidth(3);
hMCPrimaryStartY->SetMarkerStyle(8);
hMCPrimaryStartY->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryStartY->GetXaxis()->SetTitle("Primary Particle Y_{0} (cm)");
hMCPrimaryStartY->GetXaxis()->CenterTitle();

hMCPrimaryStartY->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryStartY->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryStartY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c3= new TCanvas("c3","Initial Z Pos");
c3->SetTicks();
c3->SetLogy();
c3->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryStartZ->SetLineColor(kRed);
hMCPrimaryStartZ->SetLineStyle(0);
hMCPrimaryStartZ->SetLineWidth(3);
hMCPrimaryStartZ->SetMarkerStyle(8);
hMCPrimaryStartZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryStartZ->GetXaxis()->SetTitle("Primary Particle Z_{0} (cm)");
hMCPrimaryStartZ->GetXaxis()->CenterTitle();

hMCPrimaryStartZ->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryStartZ->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryStartZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TCanvas *c4= new TCanvas("c4","Final X Pos");
c4->SetTicks();
c4->SetLogy();
c4->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryEndX->SetLineColor(kRed);
hMCPrimaryEndX->SetLineStyle(0);
hMCPrimaryEndX->SetLineWidth(3);
hMCPrimaryEndX->SetMarkerStyle(8);
hMCPrimaryEndX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryEndX->GetXaxis()->SetTitle("Primary Particle X_{F} (cm)");
hMCPrimaryEndX->GetXaxis()->CenterTitle();

hMCPrimaryEndX->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryEndX->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryEndX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TCanvas *c5= new TCanvas("c5","Final Y Pos");
c5->SetTicks();
c5->SetLogy();
c5->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryEndY->SetLineColor(kRed);
hMCPrimaryEndY->SetLineStyle(0);
hMCPrimaryEndY->SetLineWidth(3);
hMCPrimaryEndY->SetMarkerStyle(8);
hMCPrimaryEndY->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryEndY->GetXaxis()->SetTitle("Primary Particle Y_{F} (cm)");
hMCPrimaryEndY->GetXaxis()->CenterTitle();

hMCPrimaryEndY->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryEndY->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryEndY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




TCanvas *c6= new TCanvas("c6","Final Z Pos");
c6->SetTicks();
c6->SetLogy();
c6->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryEndZ->SetLineColor(kRed);
hMCPrimaryEndZ->SetLineStyle(0);
hMCPrimaryEndZ->SetLineWidth(3);
hMCPrimaryEndZ->SetMarkerStyle(8);
hMCPrimaryEndZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryEndZ->GetXaxis()->SetTitle("Primary Particle Z_{F} (cm)");
hMCPrimaryEndZ->GetXaxis()->CenterTitle();

hMCPrimaryEndZ->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryEndZ->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryEndZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c7= new TCanvas("c7","Primary Px");
c7->SetTicks();
c7->SetLogy();
c7->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryPx->SetLineColor(kRed);
hMCPrimaryPx->SetLineStyle(0);
hMCPrimaryPx->SetLineWidth(3);
hMCPrimaryPx->SetMarkerStyle(8);
hMCPrimaryPx->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryPx->GetXaxis()->SetTitle("Primary Particle P_{X} (MeV)");
hMCPrimaryPx->GetXaxis()->CenterTitle();

hMCPrimaryPx->GetYaxis()->SetTitle("Events / MeV");
hMCPrimaryPx->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryPx->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TCanvas *c8= new TCanvas("c8","Primary Py");
c8->SetTicks();
c8->SetLogy();
c8->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryPy->SetLineColor(kRed);
hMCPrimaryPy->SetLineStyle(0);
hMCPrimaryPy->SetLineWidth(3);
hMCPrimaryPy->SetMarkerStyle(8);
hMCPrimaryPy->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryPy->GetXaxis()->SetTitle("Primary Particle P_{Y} (MeV)");
hMCPrimaryPy->GetXaxis()->CenterTitle();

hMCPrimaryPy->GetYaxis()->SetTitle("Events / MeV");
hMCPrimaryPy->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryPy->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TCanvas *c9= new TCanvas("c9","Primary Pz");
c9->SetTicks();
c9->SetLogy();
c9->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryPz->SetLineColor(kRed);
hMCPrimaryPz->SetLineStyle(0);
hMCPrimaryPz->SetLineWidth(3);
hMCPrimaryPz->SetMarkerStyle(8);
hMCPrimaryPz->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryPz->GetXaxis()->SetTitle("Primary Particle P_{Z} (MeV)");
hMCPrimaryPz->GetXaxis()->CenterTitle();

hMCPrimaryPz->GetYaxis()->SetTitle("Events / MeV");
hMCPrimaryPz->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryPz->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// ### Making a TCanvas ###
TCanvas *c10= new TCanvas("c10","Projected Initial X Pos");
c10->SetTicks();
//c10->SetLogy();
c10->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryProjectedStartX->SetLineColor(kGreen);
hMCPrimaryProjectedStartX->SetLineStyle(0);
hMCPrimaryProjectedStartX->SetLineWidth(3);
hMCPrimaryProjectedStartX->SetMarkerStyle(8);
hMCPrimaryProjectedStartX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryProjectedStartX->GetXaxis()->SetTitle("Primary Particle Projected X_{0} (cm)");
hMCPrimaryProjectedStartX->GetXaxis()->CenterTitle();

hMCPrimaryProjectedStartX->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryProjectedStartX->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryProjectedStartX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c11= new TCanvas("c11","Projected Y Pos");
c11->SetTicks();
c11->SetLogy();
c11->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryProjectedStartY->SetLineColor(kGreen);
hMCPrimaryProjectedStartY->SetLineStyle(0);
hMCPrimaryProjectedStartY->SetLineWidth(3);
hMCPrimaryProjectedStartY->SetMarkerStyle(8);
hMCPrimaryProjectedStartY->SetMarkerSize(0.9);


// ### Labeling the axis ###
hMCPrimaryProjectedStartY->GetXaxis()->SetTitle("Primary Particle Projected Y_{0} (cm)");
hMCPrimaryProjectedStartY->GetXaxis()->CenterTitle();

hMCPrimaryProjectedStartY->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryProjectedStartY->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryProjectedStartY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c12= new TCanvas("c12","Projected Initial Z Pos");
c12->SetTicks();
c12->SetLogy();
c12->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryProjectedStartZ->SetLineColor(kGreen);
hMCPrimaryProjectedStartZ->SetLineStyle(0);
hMCPrimaryProjectedStartZ->SetLineWidth(3);
hMCPrimaryProjectedStartZ->SetMarkerStyle(8);
hMCPrimaryProjectedStartZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryProjectedStartZ->GetXaxis()->SetTitle("Primary Particle Projected Z_{0} (cm)");
hMCPrimaryProjectedStartZ->GetXaxis()->CenterTitle();

hMCPrimaryProjectedStartZ->GetYaxis()->SetTitle("Events / 0.5 cm");
hMCPrimaryProjectedStartZ->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryProjectedStartZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c13= new TCanvas("c13","Final X vs Z");
c13->SetTicks();
//c13->SetLogy();
c13->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryEndXvsZ->SetLineColor(kRed);
hMCPrimaryEndXvsZ->SetLineStyle(0);
hMCPrimaryEndXvsZ->SetLineWidth(3);
hMCPrimaryEndXvsZ->SetMarkerStyle(8);
hMCPrimaryEndXvsZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryEndXvsZ->GetXaxis()->SetTitle("Primary Particle Z_{f} (cm)");
hMCPrimaryEndXvsZ->GetXaxis()->CenterTitle();

hMCPrimaryEndXvsZ->GetYaxis()->SetTitle("Primary Particle X_{f} (cm)");
hMCPrimaryEndXvsZ->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryEndXvsZ->Draw("colz");

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TCanvas *c14= new TCanvas("c14","Final Y vs Z");
c14->SetTicks();
//c14->SetLogy();
c14->SetFillColor(kWhite);  

// === Histogram Drawing Settings ===
hMCPrimaryEndYvsZ->SetLineColor(kRed);
hMCPrimaryEndYvsZ->SetLineStyle(0);
hMCPrimaryEndYvsZ->SetLineWidth(3);
hMCPrimaryEndYvsZ->SetMarkerStyle(8);
hMCPrimaryEndYvsZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCPrimaryEndYvsZ->GetXaxis()->SetTitle("Primary Particle Z_{f} (cm)");
hMCPrimaryEndYvsZ->GetXaxis()->CenterTitle();

hMCPrimaryEndYvsZ->GetYaxis()->SetTitle("Primary Particle Y_{f} (cm)");
hMCPrimaryEndYvsZ->GetYaxis()->CenterTitle();

// ### Drawing ###
hMCPrimaryEndYvsZ->Draw("colz");

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,"");

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c15= new TCanvas("c15","Most upstream Spacepoint");
c15->SetTicks();
//c9->SetLogy();
c15->SetFillColor(kWhite); 

// ### Formatting the histograms ###
hdataUpstreamZPos->SetLineColor(kBlack);
hdataUpstreamZPos->SetLineStyle(0);
hdataUpstreamZPos->SetLineWidth(3);
hdataUpstreamZPos->SetMarkerStyle(8);
hdataUpstreamZPos->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataUpstreamZPos->GetXaxis()->SetTitle("Upstream Z position of the track (cm)");
hdataUpstreamZPos->GetXaxis()->CenterTitle();

hdataUpstreamZPos->GetYaxis()->SetTitle("Events / 0.5 cm");
hdataUpstreamZPos->GetYaxis()->CenterTitle();

// ### Drawing the histograms ###
hdataUpstreamZPos->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();



// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c16= new TCanvas("c16","Number of tracks vs Z");
c16->SetTicks();
c16->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataNTracksvsZpos->SetLineColor(kBlack);
hdataNTracksvsZpos->SetLineStyle(0);
hdataNTracksvsZpos->SetLineWidth(3);
hdataNTracksvsZpos->SetMarkerStyle(8);
hdataNTracksvsZpos->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataNTracksvsZpos->GetXaxis()->SetTitle("Distance from the upstream face of the TPC (cm)");
hdataNTracksvsZpos->GetXaxis()->CenterTitle();

hdataNTracksvsZpos->GetYaxis()->SetTitle("Number of Tracks");
hdataNTracksvsZpos->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataNTracksvsZpos->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c17= new TCanvas("c17","Delta X0");
c17->SetTicks();
c17->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaX->SetLineColor(kBlue);
hDeltaX->SetLineStyle(0);
hDeltaX->SetLineWidth(3);
hDeltaX->SetMarkerStyle(8);
hDeltaX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaX->GetXaxis()->SetTitle("#Delta X_{0} (Reco - True) (cm)");
hDeltaX->GetXaxis()->CenterTitle();

hDeltaX->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaX->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c18= new TCanvas("c18","Delta Y0");
c18->SetTicks();
c18->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaY->SetLineColor(kBlue);
hDeltaY->SetLineStyle(0);
hDeltaY->SetLineWidth(3);
hDeltaY->SetMarkerStyle(8);
hDeltaY->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaY->GetXaxis()->SetTitle("#Delta Y_{0} (Reco - True) (cm)");
hDeltaY->GetXaxis()->CenterTitle();

hDeltaY->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaY->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c19= new TCanvas("c19","Delta Z0");
c19->SetTicks();
c19->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaZ->SetLineColor(kBlue);
hDeltaZ->SetLineStyle(0);
hDeltaZ->SetLineWidth(3);
hDeltaZ->SetMarkerStyle(8);
hDeltaZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaZ->GetXaxis()->SetTitle("#Delta Z_{0} (Reco - True) (cm)");
hDeltaZ->GetXaxis()->CenterTitle();

hDeltaZ->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaZ->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c20= new TCanvas("c20","Alpha");
c20->SetTicks();
c20->SetFillColor(kWhite);

// ### Formatting the histograms ###
hAlpha->SetLineColor(kBlue);
hAlpha->SetLineStyle(0);
hAlpha->SetLineWidth(3);
hAlpha->SetMarkerStyle(8);
hAlpha->SetMarkerSize(0.9);

// ### Labeling the axis ###
hAlpha->GetXaxis()->SetTitle("#alpha (degrees)");
hAlpha->GetXaxis()->CenterTitle();

hAlpha->GetYaxis()->SetTitle("Events / degree");
hAlpha->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hAlpha->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------


// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c21= new TCanvas("c21","dEdX");
c21->SetTicks();
c21->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataProtondEdX->SetLineColor(kBlack);
hdataProtondEdX->SetLineStyle(0);
hdataProtondEdX->SetLineWidth(3);
hdataProtondEdX->SetMarkerStyle(8);
hdataProtondEdX->SetMarkerSize(0.9);


// ### Labeling the axis ###
hdataProtondEdX->GetXaxis()->SetTitle("dE/dX (MeV/cm)");
hdataProtondEdX->GetXaxis()->CenterTitle();

hdataProtondEdX->GetYaxis()->SetTitle("");
hdataProtondEdX->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataProtondEdX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------


// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c22= new TCanvas("c22","RR");
c22->SetTicks();
c22->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataProtonRR->SetLineColor(kBlack);
hdataProtonRR->SetLineStyle(0);
hdataProtonRR->SetLineWidth(3);
hdataProtonRR->SetMarkerStyle(8);
hdataProtonRR->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataProtonRR->GetXaxis()->SetTitle("Residual Range (cm)");
hdataProtonRR->GetXaxis()->CenterTitle();

hdataProtonRR->GetYaxis()->SetTitle("Events / 0.5 cm");
hdataProtonRR->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataProtonRR->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c23= new TCanvas("c23","Pitch");
c23->SetTicks();
c23->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataProtonTrkPitch->SetLineColor(kBlack);
hdataProtonTrkPitch->SetLineStyle(0);
hdataProtonTrkPitch->SetLineWidth(3);
hdataProtonTrkPitch->SetMarkerStyle(8);
hdataProtonTrkPitch->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataProtonTrkPitch->GetXaxis()->SetTitle("Track Pitch (cm)");
hdataProtonTrkPitch->GetXaxis()->CenterTitle();

hdataProtonTrkPitch->GetYaxis()->SetTitle("Events / 0.005 cm");
hdataProtonTrkPitch->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataProtonTrkPitch->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c24= new TCanvas("c24","Incident KE");
c24->SetTicks();
c24->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataProtonIncidentKE->SetLineColor(kBlack);
hdataProtonIncidentKE->SetLineStyle(0);
hdataProtonIncidentKE->SetLineWidth(3);
hdataProtonIncidentKE->SetMarkerStyle(8);
hdataProtonIncidentKE->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataProtonIncidentKE->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
hdataProtonIncidentKE->GetXaxis()->CenterTitle();

hdataProtonIncidentKE->GetYaxis()->SetTitle("Events / 50 MeV");
hdataProtonIncidentKE->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataProtonIncidentKE->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c25= new TCanvas("c25","Interaction KE");
c25->SetTicks();
c25->SetFillColor(kWhite);

// ### Formatting the histograms ###
fProtonInteractions->SetLineColor(kBlack);
fProtonInteractions->SetLineStyle(0);
fProtonInteractions->SetLineWidth(3);
fProtonInteractions->SetMarkerStyle(8);
fProtonInteractions->SetMarkerSize(0.9);

// ### Labeling the axis ###
fProtonInteractions->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
fProtonInteractions->GetXaxis()->CenterTitle();

fProtonInteractions->GetYaxis()->SetTitle("Events / 50 MeV");
fProtonInteractions->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
fProtonInteractions->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------


// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c26= new TCanvas("c26","Cross-Section");
c26->SetTicks();
c26->SetFillColor(kWhite);

// ### Formatting the histograms ###
fCrossSection->SetLineColor(kBlack);
fCrossSection->SetLineStyle(0);
fCrossSection->SetLineWidth(3);
fCrossSection->SetMarkerStyle(8);
fCrossSection->SetMarkerSize(0.9);

// ### Labeling the axis ###
fCrossSection->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
fCrossSection->GetXaxis()->CenterTitle();

fCrossSection->GetYaxis()->SetTitle("#sigma(Inclusive #pi) (barns)");
fCrossSection->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
fCrossSection->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c27= new TCanvas("c27","Cross-Section");
c27->SetTicks();
c27->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataProtondEdXvsRR->SetLineColor(kBlack);
hdataProtondEdXvsRR->SetLineStyle(0);
hdataProtondEdXvsRR->SetLineWidth(3);
hdataProtondEdXvsRR->SetMarkerStyle(8);
hdataProtondEdXvsRR->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataProtondEdXvsRR->GetXaxis()->SetTitle("Residual Range (cm)");
hdataProtondEdXvsRR->GetXaxis()->CenterTitle();

hdataProtondEdXvsRR->GetYaxis()->SetTitle("dE/dX (MeV/cm)");
hdataProtondEdXvsRR->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataProtondEdXvsRR->Draw("colz");

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c28= new TCanvas("c28","Delta XF");
c28->SetTicks();
c28->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaEndX->SetLineColor(kBlue);
hDeltaEndX->SetLineStyle(0);
hDeltaEndX->SetLineWidth(3);
hDeltaEndX->SetMarkerStyle(8);
hDeltaEndX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaEndX->GetXaxis()->SetTitle("#Delta X_{F} (Reco - True) (cm)");
hDeltaEndX->GetXaxis()->CenterTitle();

hDeltaEndX->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaEndX->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaEndX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c29= new TCanvas("c29","Delta YF");
c29->SetTicks();
c29->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaEndY->SetLineColor(kBlue);
hDeltaEndY->SetLineStyle(0);
hDeltaEndY->SetLineWidth(3);
hDeltaEndY->SetMarkerStyle(8);
hDeltaEndY->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaEndY->GetXaxis()->SetTitle("#Delta Y_{F} (Reco - True) (cm)");
hDeltaEndY->GetXaxis()->CenterTitle();

hDeltaEndY->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaEndY->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaEndY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c30= new TCanvas("c30","Delta ZF");
c30->SetTicks();
c30->SetFillColor(kWhite);

// ### Formatting the histograms ###
hDeltaEndZ->SetLineColor(kBlue);
hDeltaEndZ->SetLineStyle(0);
hDeltaEndZ->SetLineWidth(3);
hDeltaEndZ->SetMarkerStyle(8);
hDeltaEndZ->SetMarkerSize(0.9);

// ### Labeling the axis ###
hDeltaEndZ->GetXaxis()->SetTitle("#Delta Z_{F} (Reco - True) (cm)");
hDeltaEndZ->GetXaxis()->CenterTitle();

hDeltaEndZ->GetYaxis()->SetTitle("Events / 0.5 cm");
hDeltaEndZ->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hDeltaEndZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c31= new TCanvas("c31","Delta ZF Stacked");
c31->SetTicks();
c31->SetFillColor(kWhite);

// ### Formatting the histograms ###

// ### Labeling the axis ###
//hStackDeltaEndZ->GetXaxis()->SetTitle("#Delta Z_{F} (Reco - True) (cm)");
//hStackDeltaEndZ->GetXaxis()->CenterTitle();

//hStackDeltaEndZ->GetYaxis()->SetTitle("Events / 0.5 cm");
//hStackDeltaEndZ->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hStackDeltaEndZ->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->AddEntry(hDeltaEndZInElastic, "InElastic");
leg->AddEntry(hDeltaEndZNeutronInElastic, "Neutron InElastic");
leg->AddEntry(hDeltaEndZHadElastic, "Hadronic Elastic");
leg->AddEntry(hDeltaEndZnCap, "Neutron Capture");
leg->AddEntry(hDeltaEndZnuclearCapatureAtRest, "Nuclear Capture at Rest");
leg->AddEntry(hDeltaEndZDecay, "Decay");
leg->AddEntry(hDeltaEndZKaonZeroInElastic, "K0 InElastic");
leg->AddEntry(hDeltaEndZCoulombScat, "Coulomb Scatter");
leg->AddEntry(hDeltaEndZMuMinusCapture, "#mu^{-} Capture");
leg->AddEntry(hDeltaEndZProtonInelastic, "Proton InElastic");
leg->Draw();
// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c32= new TCanvas("c32","Delta YF Stacked");
c32->SetTicks();
c32->SetFillColor(kWhite);

// ### Formatting the histograms ###

// ### Labeling the axis ###
//hStackDeltaEndZ->GetXaxis()->SetTitle("#Delta Z_{F} (Reco - True) (cm)");
//hStackDeltaEndZ->GetXaxis()->CenterTitle();

//hStackDeltaEndZ->GetYaxis()->SetTitle("Events / 0.5 cm");
//hStackDeltaEndZ->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hStackDeltaEndY->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->AddEntry(hDeltaEndYInElastic, "InElastic");
leg->AddEntry(hDeltaEndYNeutronInElastic, "Neutron InElastic");
leg->AddEntry(hDeltaEndYHadElastic, "Hadronic Elastic");
leg->AddEntry(hDeltaEndYnCap, "Neutron Capture");
leg->AddEntry(hDeltaEndYnuclearCapatureAtRest, "Nuclear Capture at Rest");
leg->AddEntry(hDeltaEndYDecay, "Decay");
leg->AddEntry(hDeltaEndYKaonZeroInElastic, "K0 InElastic");
leg->AddEntry(hDeltaEndYCoulombScat, "Coulomb Scatter");
leg->AddEntry(hDeltaEndYMuMinusCapture, "#mu^{-} Capture");
leg->AddEntry(hDeltaEndYProtonInelastic, "Proton InElastic");
leg->Draw();
// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c33= new TCanvas("c33","Delta XF Stacked");
c33->SetTicks();
c33->SetFillColor(kWhite);

// ### Formatting the histograms ###

// ### Labeling the axis ###
//hStackDeltaEndZ->GetXaxis()->SetTitle("#Delta Z_{F} (Reco - True) (cm)");
//hStackDeltaEndZ->GetXaxis()->CenterTitle();

//hStackDeltaEndZ->GetYaxis()->SetTitle("Events / 0.5 cm");
//hStackDeltaEndZ->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hStackDeltaEndX->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->AddEntry(hDeltaEndXInElastic, "InElastic");
leg->AddEntry(hDeltaEndXNeutronInElastic, "Neutron InElastic");
leg->AddEntry(hDeltaEndXHadElastic, "Hadronic Elastic");
leg->AddEntry(hDeltaEndXnCap, "Neutron Capture");
leg->AddEntry(hDeltaEndXnuclearCapatureAtRest, "Nuclear Capture at Rest");
leg->AddEntry(hDeltaEndXDecay, "Decay");
leg->AddEntry(hDeltaEndXKaonZeroInElastic, "K0 InElastic");
leg->AddEntry(hDeltaEndXCoulombScat, "Coulomb Scatter");
leg->AddEntry(hDeltaEndXMuMinusCapture, "#mu^{-} Capture");
leg->AddEntry(hDeltaEndXProtonInelastic, "Proton InElastic");
leg->Draw();


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c34= new TCanvas("c34","Energy Loss");
c34->SetTicks();
c34->SetFillColor(kWhite);

// ### Formatting the histograms ###
// ### Formatting the histograms ###
hMCELossUpstream->SetLineColor(kBlue);
hMCELossUpstream->SetLineStyle(0);
hMCELossUpstream->SetLineWidth(3);
hMCELossUpstream->SetMarkerStyle(8);
hMCELossUpstream->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCELossUpstream->GetXaxis()->SetTitle("Energy Loss Prior to Entering the TPC (MeV)");
hMCELossUpstream->GetXaxis()->CenterTitle();

hMCELossUpstream->GetYaxis()->SetTitle("Events / MeV");
hMCELossUpstream->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hMCELossUpstream->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->Draw();



// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c35= new TCanvas("c35","True Length");
c35->SetTicks();
c35->SetFillColor(kWhite);

// ### Formatting the histograms ###
// ### Formatting the histograms ###
hTrueLength->SetLineColor(kBlue);
hTrueLength->SetLineStyle(0);
hTrueLength->SetLineWidth(3);
hTrueLength->SetMarkerStyle(8);
hTrueLength->SetMarkerSize(0.9);

// ### Labeling the axis ###
hTrueLength->GetXaxis()->SetTitle("Length of the particle (cm)");
hTrueLength->GetXaxis()->CenterTitle();

hTrueLength->GetYaxis()->SetTitle("Events / 0.5 cm");
hTrueLength->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hTrueLength->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->Draw();


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c36= new TCanvas("c36","Reco Length");
c36->SetTicks();
c36->SetFillColor(kWhite);

// ### Formatting the histograms ###
// ### Formatting the histograms ###
hRecoLength->SetLineColor(kBlue);
hRecoLength->SetLineStyle(0);
hRecoLength->SetLineWidth(3);
hRecoLength->SetMarkerStyle(8);
hRecoLength->SetMarkerSize(0.9);

// ### Labeling the axis ###
hRecoLength->GetXaxis()->SetTitle("Length of the particle (cm)");
hRecoLength->GetXaxis()->CenterTitle();

hRecoLength->GetYaxis()->SetTitle("Events / 0.5 cm");
hRecoLength->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hRecoLength->Draw();

// ############################
// # Setting the Latex Header #
// ############################
//TLatex *t = new TLatex();
t->SetNDC();
t->SetTextFont(62);
t->SetTextSize(0.04);
t->SetTextAlign(40);
t->DrawLatex(0.1,0.90,"LArIAT Preliminary");
t->DrawLatex(0.13,0.84,""); 

// ######################
// # Setting the Legend #
// ######################
//TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Proton MC");
leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ===========================================================================================
// ============================  Writing out histograms to ROOT File =========================
// ===========================================================================================
// ###############################################
// ### Creating a file to output my histograms ###
// ###############################################
TFile myfile("ProtonMCXSection_weighted_histos.root","RECREATE");

// ### Reco Info ###
hdataUpstreamZPos->Write();
hdataNTracksvsZpos->Write();

// ### MC Info ###
hMCPrimaryStartX->Write();
hMCPrimaryStartY->Write();
hMCPrimaryStartZ->Write();
hMCPrimaryProjectedStartX->Write();
hMCPrimaryProjectedStartY->Write();
hMCPrimaryProjectedStartZ->Write();
hMCPrimaryEndX->Write();
hMCPrimaryEndY->Write();
hMCPrimaryEndZ->Write();
hMCPrimaryEndXvsZ->Write();
hMCPrimaryEndYvsZ->Write();
hMCPrimaryProcess->Write();
hMCPrimaryPx->Write();
hMCPrimaryPy->Write();
hMCPrimaryPz->Write();
hAlpha->Write();
hDeltaX->Write();
hDeltaY->Write();
hDeltaZ->Write();  
hMCInitalKE->Write(); 
hdataProtondEdX->Write();
hdataProtonRR->Write();
hdataProtonTrkPitch->Write();
hdataProtonIncidentKE->Write();
fProtonInteractions->Write();
fCrossSection->Write();
hdataProtondEdXvsRR->Write();
hDeltaEndX->Write();
hDeltaEndY->Write();
hDeltaEndZ->Write();

hTrueLength->Write();
hRecoLength->Write();

hDeltaEndZInElastic->Write();
hDeltaEndZNeutronInElastic->Write();
hDeltaEndZHadElastic->Write();
hDeltaEndZnCap->Write();
hDeltaEndZnuclearCapatureAtRest->Write();
hDeltaEndZDecay->Write();
hDeltaEndZKaonZeroInElastic->Write();
hDeltaEndZCoulombScat->Write();
hDeltaEndZMuMinusCapture->Write();
hDeltaEndZProtonInelastic->Write();

hDeltaEndYInElastic->Write();
hDeltaEndYNeutronInElastic->Write();
hDeltaEndYHadElastic->Write();
hDeltaEndYnCap->Write();
hDeltaEndYnuclearCapatureAtRest->Write();
hDeltaEndYDecay->Write();
hDeltaEndYKaonZeroInElastic->Write();
hDeltaEndYCoulombScat->Write();
hDeltaEndYMuMinusCapture->Write();
hDeltaEndYProtonInelastic->Write();

hDeltaEndXInElastic->Write();
hDeltaEndXNeutronInElastic->Write();
hDeltaEndXHadElastic->Write();
hDeltaEndXnCap->Write();
hDeltaEndXnuclearCapatureAtRest->Write();
hDeltaEndXDecay->Write();
hDeltaEndXKaonZeroInElastic->Write();
hDeltaEndXCoulombScat->Write();
hDeltaEndXMuMinusCapture->Write();
hDeltaEndXProtonInelastic->Write();

hStackDeltaEndZ->Write();
hStackDeltaEndY->Write();
hStackDeltaEndX->Write();

hMCELossUpstream->Write();

hdataProtonTrackEndX->Write();
hdataProtonTrackEndY->Write();
hdataProtonTrackEndZ->Write();
hdataProtonTrackStartX->Write();
hdataProtonTrackStartY->Write();
hdataProtonTrackStartZ->Write();

}//<---End Loop Function






//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
// 				### Function for plotting the Low Z track location ###
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
void Proton::LowZCut()
{

// ################################################################
// ### Initializing variables for study of low Z track location ###
// ################################################################
bool LowZTrackInTPC = false;
bool trackwZLT1 = false, trackwZLT2 = false, trackwZLT3 = false;
bool trackwZLT4 = false, trackwZLT5 = false, trackwZLT6 = false;
bool trackwZLT7 = false, trackwZLT8 = false, trackwZLT9 = false;
bool trackwZLT10 = false, trackwZLT11 = false, trackwZLT12 = false;
bool trackwZLT13 = false, trackwZLT14 = false, trackwZLT15 = false;
bool trackwZLT16 = false, trackwZLT17 = false, trackwZLT18 = false;
bool trackwZLT19 = false, trackwZLT20 = false, trackwZLT21 = false;
bool trackwZLT22 = false, trackwZLT23 = false, trackwZLT24 = false;
bool trackwZLT25 = false, trackwZLT26 = false, trackwZLT27 = false;
bool trackwZLT28 = false, trackwZLT29 = false, trackwZLT30 = false;
    
// ########################################################## 
// ### Initializing the counters for the Z track location ###
// ##########################################################
int trkwZLT1 = 0, trkwZLT2 = 0, trkwZLT3 = 0;
int trkwZLT4 = 0, trkwZLT5 = 0, trkwZLT6 = 0;
int trkwZLT7 = 0, trkwZLT8 = 0, trkwZLT9 = 0;
int trkwZLT10 = 0, trkwZLT11 = 0, trkwZLT12 = 0;
int trkwZLT13 = 0, trkwZLT14 = 0, trkwZLT15 = 0;
int trkwZLT16 = 0, trkwZLT17 = 0, trkwZLT18 = 0;
int trkwZLT19 = 0, trkwZLT20 = 0, trkwZLT21 = 0;
int trkwZLT22 = 0, trkwZLT23 = 0, trkwZLT24 = 0;
int trkwZLT25 = 0, trkwZLT26 = 0, trkwZLT27 = 0;
int trkwZLT28 = 0, trkwZLT29 = 0, trkwZLT30 = 0;

// ###########################
// ### Looping over tracks ###
// ###########################
for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
   {
   // === Start each track with the boolians set to false ===
   trackwZLT1 = false; trackwZLT2 = false; trackwZLT3 = false;
   trackwZLT4 = false; trackwZLT5 = false; trackwZLT6 = false;
   trackwZLT7 = false; trackwZLT8 = false; trackwZLT9 = false;
   trackwZLT10 = false; trackwZLT11 = false; trackwZLT12 = false;
   trackwZLT13 = false; trackwZLT14 = false; trackwZLT15 = false;
   trackwZLT16 = false; trackwZLT17 = false; trackwZLT18 = false;
   trackwZLT19 = false; trackwZLT20 = false; trackwZLT21 = false;
   trackwZLT22 = false; trackwZLT23 = false; trackwZLT24 = false;
   trackwZLT25 = false; trackwZLT26 = false; trackwZLT27 = false;
   trackwZLT28 = false; trackwZLT29 = false; trackwZLT30 = false;
   
   // ##################################################
   // ### Looping over the spacepoints for the track ###
   // ##################################################
   for(size_t nspts = 0; nspts < nTrajPoint[nTPCtrk]; nspts++)
      {
      
      // ### Putting the number of tracks w/ spts in this area of Z ###
	 if(trkz[nTPCtrk][nspts] < 1 && !trackwZLT1)
	 	{trkwZLT1++; trackwZLT1 = true;} 
	 if(trkz[nTPCtrk][nspts] < 2 && !trackwZLT2)
	 	{trkwZLT2++; trackwZLT2 = true;} 
	 if(trkz[nTPCtrk][nspts] < 3 && !trackwZLT3)
	 	{trkwZLT3++; trackwZLT3 = true;} 
	 if(trkz[nTPCtrk][nspts] < 4 && !trackwZLT4)
	 	{trkwZLT4++; trackwZLT4 = true;} 
	 if(trkz[nTPCtrk][nspts] < 5 && !trackwZLT5)
	 	{trkwZLT5++; trackwZLT5 = true;} 
	 if(trkz[nTPCtrk][nspts] < 6 && !trackwZLT6)
	 	{trkwZLT6++; trackwZLT6 = true;} 
	 if(trkz[nTPCtrk][nspts] < 7 && !trackwZLT7)
	 	{trkwZLT7++; trackwZLT7 = true;}
	 if(trkz[nTPCtrk][nspts] < 8 && !trackwZLT8)
	 	{trkwZLT8++; trackwZLT8 = true;}
	 if(trkz[nTPCtrk][nspts] < 9 && !trackwZLT9)
	 	{trkwZLT9++; trackwZLT9 = true;}
	 if(trkz[nTPCtrk][nspts] < 10 && !trackwZLT10)
	 	{trkwZLT10++; trackwZLT10 = true;}
	 if(trkz[nTPCtrk][nspts] < 11 && !trackwZLT11)
	 	{trkwZLT11++; trackwZLT11 = true;}
	 if(trkz[nTPCtrk][nspts] < 12 && !trackwZLT12)
	 	{trkwZLT12++; trackwZLT12 = true;} 
	 if(trkz[nTPCtrk][nspts] < 13 && !trackwZLT13)
	 	{trkwZLT13++; trackwZLT13 = true;}
	 if(trkz[nTPCtrk][nspts] < 14 && !trackwZLT14)
	 	{trkwZLT14++; trackwZLT14 = true;}
	 if(trkz[nTPCtrk][nspts] < 15 && !trackwZLT15)
	 	{trkwZLT15++; trackwZLT15 = true;}
	 if(trkz[nTPCtrk][nspts] < 16 && !trackwZLT16)
	 	{trkwZLT16++; trackwZLT16 = true;}
	 if(trkz[nTPCtrk][nspts] < 17 && !trackwZLT17)
	 	{trkwZLT17++; trackwZLT17 = true;}
	 if(trkz[nTPCtrk][nspts] < 18 && !trackwZLT18)
	 	{trkwZLT18++; trackwZLT18 = true;}		
	 if(trkz[nTPCtrk][nspts] < 19 && !trackwZLT19)
	 	{trkwZLT19++; trackwZLT19 = true;}
	 if(trkz[nTPCtrk][nspts] < 20 && !trackwZLT20)
	 	{trkwZLT20++; trackwZLT20 = true;}
	 if(trkz[nTPCtrk][nspts] < 21 && !trackwZLT21)
	 	{trkwZLT21++; trackwZLT21 = true;} 
	 if(trkz[nTPCtrk][nspts] < 22 && !trackwZLT22)
	 	{trkwZLT22++; trackwZLT22 = true;}
	 if(trkz[nTPCtrk][nspts] < 23 && !trackwZLT23)
	 	{trkwZLT23++; trackwZLT23 = true;}
	 if(trkz[nTPCtrk][nspts] < 24 && !trackwZLT24)
	 	{trkwZLT24++; trackwZLT24 = true;}
	 if(trkz[nTPCtrk][nspts] < 25 && !trackwZLT25)
	 	{trkwZLT25++; trackwZLT25 = true;}
	 if(trkz[nTPCtrk][nspts] < 26 && !trackwZLT26)
	 	{trkwZLT26++; trackwZLT26 = true;}
	 if(trkz[nTPCtrk][nspts] < 27 && !trackwZLT27)
	 	{trkwZLT27++; trackwZLT27 = true;}		
	 if(trkz[nTPCtrk][nspts] < 28 && !trackwZLT28)
	 	{trkwZLT28++; trackwZLT28 = true;}
	 if(trkz[nTPCtrk][nspts] < 29 && !trackwZLT29)
	 	{trkwZLT29++; trackwZLT29 = true;}
	 if(trkz[nTPCtrk][nspts] < 30 && !trackwZLT30)
	 	{trkwZLT30++; trackwZLT30 = true;}
   
   
      }//<---End nspts
   
   }///<---End nTPCtrk

hdataNTracksvsZpos->Fill(1, trkwZLT1);
hdataNTracksvsZpos->Fill(2, trkwZLT2);
hdataNTracksvsZpos->Fill(3, trkwZLT3);
hdataNTracksvsZpos->Fill(4, trkwZLT4);
hdataNTracksvsZpos->Fill(5, trkwZLT5);
hdataNTracksvsZpos->Fill(6, trkwZLT6);
hdataNTracksvsZpos->Fill(7, trkwZLT7);
hdataNTracksvsZpos->Fill(8, trkwZLT8);
hdataNTracksvsZpos->Fill(9, trkwZLT9);
hdataNTracksvsZpos->Fill(10, trkwZLT10);
hdataNTracksvsZpos->Fill(11, trkwZLT11);
hdataNTracksvsZpos->Fill(12, trkwZLT12);
hdataNTracksvsZpos->Fill(13, trkwZLT13);
hdataNTracksvsZpos->Fill(14, trkwZLT14);
hdataNTracksvsZpos->Fill(15, trkwZLT15);
hdataNTracksvsZpos->Fill(16, trkwZLT16);
hdataNTracksvsZpos->Fill(17, trkwZLT17);
hdataNTracksvsZpos->Fill(18, trkwZLT18); 
hdataNTracksvsZpos->Fill(19, trkwZLT19);
hdataNTracksvsZpos->Fill(20, trkwZLT20);
hdataNTracksvsZpos->Fill(21, trkwZLT21);
hdataNTracksvsZpos->Fill(22, trkwZLT22);
hdataNTracksvsZpos->Fill(23, trkwZLT23);
hdataNTracksvsZpos->Fill(24, trkwZLT24);
hdataNTracksvsZpos->Fill(25, trkwZLT25);
hdataNTracksvsZpos->Fill(26, trkwZLT26);
hdataNTracksvsZpos->Fill(27, trkwZLT27); 
hdataNTracksvsZpos->Fill(28, trkwZLT28);
hdataNTracksvsZpos->Fill(29, trkwZLT29);
hdataNTracksvsZpos->Fill(30, trkwZLT30);



}//<----End LowZCut() Function*/
