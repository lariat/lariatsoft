#define DataAnalysis_cxx
#include "DataAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ===================================================================================================================
// ====================================       PUT HISTOGRAMS HERE           ==========================================
// ===================================================================================================================

/////////////////////////////////// Events with and without WC Tracks //////////////////////////////////////////
TH1D *hdataWCTrackExist = new TH1D("hdataWCTrackExist", "Existance of WCTrack", 2, 0, 2);


/////////////////////////////////// Most Upstream Z point of tracks //////////////////////////////////////////
TH1D *hdataUpstreamZPos = new TH1D("hdataUpstreamZPos", "Most upstream spacepoint of all TPC Tracks", 20, 0, 10);

/////////////////////////////////// Number of Tracks in the TPC versus distance //////////////////////////////////////////
TH2D *hdataNTracksvsZpos = new TH2D("hdataNTracksvsZpos", "Number of TPC tracks vs Z ", 30, 0, 30, 20, 0, 10);

/////////////////////////////////// Number of Tracks in the TPC versus distance //////////////////////////////////////////
TH2D *hdataWCTrackMomentumVSTOF = new TH2D("hdataWCTrackMomentumVSTOF", "TOF vs WCTrack Momentum", 250, 0, 2500, 200, 0, 100);

/////////////////////////////////// Delta WCTrack X //////////////////////////////////////////
TH1D *hdataDeltaWCTrkX = new TH1D("hdataDeltaWCTrkX", "#Delta X TPC/WC Track", 160, -40, 40);

/////////////////////////////////// Delta WCTrack Y //////////////////////////////////////////
TH1D *hdataDeltaWCTrkY = new TH1D("hdataDeltaWCTrkY", "#Delta Y TPC/WC Track", 160, -40, 40);

/////////////////////////////////// Number of Matched Tracks //////////////////////////////////////////
TH1D *hdataNMatchTPCWCTrk = new TH1D("hdataNMatchTPCWCTrk", "Number of matched TPC/WC Tracks", 20, 0, 10);

/////////////////////////////////// Alpha Between WC and TPC Tracks //////////////////////////////////////////
TH1D *hdataAlpha = new TH1D("hdataAlpha", "#alpha between WC and TPC Track", 90, 0, 90);

/////////////////////////////////// TPC Track Theta //////////////////////////////////////////
TH1D *hdataTPCTheta = new TH1D("hdataTPCTheta", "TPC Track Theta", 180, 0, 90);

/////////////////////////////////// TPC Track Phi //////////////////////////////////////////
TH1D *hdataTPCPhi = new TH1D("hdataTPCPhi", "TPC Track Phi", 360, 0, 360);

/////////////////////////////////// Wire Chamber Theta //////////////////////////////////////////
TH1D *hdataWCTheta = new TH1D("hdataWCTheta", "WC Track Theta", 180, 0, 90);

/////////////////////////////////// Wire Chamber Phi //////////////////////////////////////////
TH1D *hdataWCPhi = new TH1D("hdataWCPhi", "WC Track Phi", 360, 0, 360);

/////////////////////////////////// WCTRK Momentum Histogram (MeV) //////////////////////////////////////////
TH1D *hdataWCTRKMomentum = new TH1D("hdataWCTRKMomentum", "WCtrk Momentum (MeV)", 250, 0, 2500);

/////////////////////////////////// "Pion" initial Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hdataPionInitialMomentum = new TH1D("hdataPionInitialMomentum", "Pion Initial Momentum (MeV)", 500, 0, 2500);

/////////////////////////////////// "Pion" dE/dX //////////////////////////////////////////
TH1D *hdataPiondEdX = new TH1D("hdataPiondEdX", "Pion dE/dX", 200, 0, 50);

/////////////////////////////////// "Pion" Residual Range //////////////////////////////////////////
TH1D *hdataPionRR = new TH1D("hdataPionRR", "Pion Residual Range", 400, -100, 100);

/////////////////////////////////// "Pion" Track Pitch //////////////////////////////////////////
TH1D *hdataPionTrkPitch = new TH1D("hdataPionTrkPitch", "Track Pitch", 1000, 0, 5);

///////////////////////////////// "Pion dE/dX vs RR ///////////////////////////////////////////
TH2D *hdataPiondEdXvsRR = new TH2D("", "dE/dX vs Residual Range",200, 0, 100, 200, 0, 50);

/////////////////////////////////// "Pion" Incident to the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *hdataPionIncidentKE = new TH1D("hdataPionIncidentKE", "Pion Incident Kinetic Energy (MeV)", 40, 0, 2000);

/////////////////////////////////// "Pion" Exiting the slab Kinetic Energy (MeV) //////////////////////////////////////////
TH1D *fPionInteractions = new TH1D("fPionInteractions", "Pion Out Kinetic Energy (MeV)", 40, 0, 2000);

/////////////////////////////////// Cross-Section //////////////////////////////////////////////////////////////////////
TH1F *fCrossSection = new TH1F("fCrossSection", "Cross-Section", 40, 0, 2000);

////////////////////////////////// Pion Track Start and End Positions //////////////////////////////////////////////////
TH1D *hdataPionTrackEndX = new TH1D("hdataPionTrackEndX", "Pion Track End X Position", 50, 0, 50);
TH1D *hdataPionTrackEndY = new TH1D("hdataPionTrackEndY", "Pion Track End Y Position", 40, -20, 20);
TH1D *hdataPionTrackEndZ = new TH1D("hdataPionTrackEndZ", "Pion Track End Z Position", 100, 0, 100);

TH1D *hdataPionTrackStartX = new TH1D("hdataPionTrackStartX", "Pion Track Start X Position", 50, 0, 50);
TH1D *hdataPionTrackStartY = new TH1D("hdataPionTrackStartY", "Pion Track Start Y Position", 40, -20, 20);
TH1D *hdataPionTrackStartZ = new TH1D("hdataPionTrackStartZ", "Pion Track Start Z Position", 100, 0, 100);


/////////////////////////////////// Reconstructed Length //////////////////////////////////////////
TH1D *hRecoLength = new TH1D("hRecoLength", "#Reconstructed Length of the Primary Particle inside the TPC", 200, 0 , 100);


void DataAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DataAnalysis.C
//      Root > DataAnalysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

float particle_mass = 139.57; //<---Mass of Pion in MeV


// #################################################
// ### Setting a flag to exclude exiting tracks  ###
// ###  from the numerator of the cross-section  ###
// #################################################
   
bool ExitingTrack = false;

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


// ########################################################
// ###      Scaling factor for through going muon       ###
// ### Note: If you have a 10% contamination than your  ###
// ###            scale factor should be 0.90           ###
// ########################################################
double MuonContaminationScaleFactor = 0.90;


// ----------------------------------------------------------------
// Create the cross section from the incident and interaction plots
// ----------------------------------------------------------------
int MatchWCTrackIndex[10] = {0};

// ### The assumed energy loss between the cryostat and the TPC ###
float entryTPCEnergyLoss = 36; //MeV

//float entryTPCEnergyLoss = 45; //MeV


// ### The assumed mass of the incident particle (here we assume a pion) ###
float mass = 139.57;

float rho = 1400; //kg/m^3
//  float cm_per_m = 100;
float molar_mass = 39.9; //g/mol
float g_per_kg = 1000; 
float avogadro = 6.02e+23; //number/mol
float number_density = rho*g_per_kg/molar_mass*avogadro;
float slab_width = 0.0045;//in m

float RecoLength = 0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;

// ##########################################################
// ### Putting in some counters for event reduction table ###
// ##########################################################
int nTotalEvents = 0, nEvtsWCTrack = 0, nEvtsTrackZPos = 0, nEvtsWCTrackMatch = 0;
int nEvtsTOF = 0, nEvtsPID = 0, nEventsPassingAlpha = 0, nLowZTrkEvents = 0;

int nEvntsTPC = 0;

int nNonShowerEvents = 0;

// ################################
// ### Loop over events (jentry)###
// ################################   
for (Long64_t jentry=0; jentry<nentries;jentry++)
//for (Long64_t jentry=0; jentry<20000;jentry++)  
   {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
      
   // #############################
   // ### Counting Total Events ###
   // #############################
   nTotalEvents++;
   
   // === Outputting every nEvents to the screen ===
   if(nTotalEvents % 1000 == 0){std::cout<<"Event = "<<nTotalEvents<<std::endl;}   
   
   ExitingTrack = false;
   
   //=======================================================================================================================
   //						Wire Chamber Track Cuts
   //=======================================================================================================================
   // ########################################
   // ### Skipping events with no WC Track ###
   // ########################################
   if(nwctrks < 1){hdataWCTrackExist->Fill(0); continue;}
   // ### Counting Events w/ WC Track ###
   hdataWCTrackExist->Fill(1);
   nEvtsWCTrack++;
   
   //=======================================================================================================================
   //						 TOF Event Selection (ns)
   //=======================================================================================================================
   
   bool tofGood = true;
   if(ntof < 1){continue;}
   for(int mmtof = 0; mmtof < ntof; mmtof++)
      {
      if(tofObject[mmtof] < 0 && tofObject[mmtof] > 30)
         {tofGood = false;}
      
      }//<---End mmtof
   
   if(!tofGood){continue;}
   nEvtsTOF++;
   
   
   //=======================================================================================================================
   //						 Putting in a GOOD TPC Cut (looking for nhits > 0
   //=======================================================================================================================
   if(nhits < 1){continue;}
   nEvntsTPC++;
   
   
   
   // ======================================================================================================================
   //						  Particle ID Filter
   // ======================================================================================================================
   
   bool GoodPID = false;

   // ### Loop over the WCTracks and TOF Objects ###
   for (int numWCTrk = 0; numWCTrk < nwctrks; numWCTrk++)
      {
      
      // ### Checking the number of TOF objects ###
      int TOFObject = numWCTrk;
      hdataWCTrackMomentumVSTOF->Fill(wctrk_momentum[numWCTrk] , tofObject[TOFObject]);
      // ### If we have more WCObjects 
      if(TOFObject > ntof){TOFObject = ntof -1;}
      
      if(wctrk_momentum[numWCTrk] > 100 && wctrk_momentum[numWCTrk] < 1500 && 
         tofObject[TOFObject] > 10 && tofObject[TOFObject] < 25)
         {GoodPID = true;}
      
      
      
      }//<---end numWCTrk
   if(!GoodPID){continue;}
   
   nEvtsPID++;
   
   
   
   //=======================================================================================================================
   //						Low Z Spacepoint Track Cut
   //=======================================================================================================================
   
   // ### Boolian for events w/ track which ###
   // ###     starts at the front face      ###
   bool TrackSptsZCut = false;
   
   // ### Recording the index of the track which ###
   // ###   starts at the front face of the TPC  ###
   bool PreLimTrackIndex[500] = {false};
   
   // ###########################
   // ### Looping over tracks ###
   // ###########################
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      
      float tempZpoint = 100;
      
      // ##################################################
      // ### Looping over the spacepoints for the track ###
      // ##################################################
      for(size_t ntrjpts = 0; ntrjpts < nTrajPoint[nTPCtrk]; ntrjpts++)
         {
	 // ################################################################################ 
         // ### Tracking the lowest Z point that is inside fiducial boundries of the TPC ###
	 // ################################################################################
	 if( trjPt_Z[nTPCtrk][ntrjpts] < tempZpoint  && trjPt_Z[nTPCtrk][ntrjpts] > ZLowerFid && 
	     trjPt_Z[nTPCtrk][ntrjpts] < ZUpperFid && trjPt_X[nTPCtrk][ntrjpts] > XLowerFid && trjPt_X[nTPCtrk][ntrjpts] < XUpperFid &&
	     trjPt_Y[nTPCtrk][ntrjpts] < YUpperFid && trjPt_Y[nTPCtrk][ntrjpts] > YLowerFid )
	     
	    {tempZpoint = trjPt_Z[nTPCtrk][ntrjpts];}//<---End looking for the most upstream point
        
	 // ### Only passing events with a track that has ###
	 // ###  a spacepoint within the first N cm in Z  ### 
	 // ###    And requiring it to be inside the TPC  ###
	 if(tempZpoint < FirstSpacePointZPos)
	    {TrackSptsZCut = true;}
	 }//<---End looping over ntrjpts	
	 
      // ### Filling the most upstream spacepoint for this track ###
      hdataUpstreamZPos->Fill(tempZpoint);
      
      // ### Recording that this track is a "good Track if ###
      // ###  it has a space point in the first N cm in Z  ###
      if(TrackSptsZCut){ PreLimTrackIndex[nTPCtrk] = true;}
      	 
      }//<---End nTPCtrk loop
      
   // ###############################################
   // ### Skipping events that don't have a track ###
   // ###   in the front of the TPC (Z) Position  ###
   // ###############################################
   if(!TrackSptsZCut){continue;}
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
      for(size_t ntrjpts = 0; ntrjpts < nTrajPoint[nTPCtrk]; ntrjpts++)
         {
	 
	 // ##################################################
	 // ### Count this track if it has a spacepoint in ###
	 // ###       the low Z region of the TPC          ###
	 // ##################################################
	 if(trjPt_Z[nTPCtrk][ntrjpts] < UpperPartOfTPC)
	    {  
	    if(trjPt_Y[nTPCtrk][ntrjpts] > YLowerFid && trjPt_Y[nTPCtrk][ntrjpts] < YUpperFid && 
	       trjPt_X[nTPCtrk][ntrjpts] > XLowerFid && trjPt_X[nTPCtrk][ntrjpts] < XUpperFid)
	        {LowZTrackInTPC = true; templowz1 = trjPt_Z[nTPCtrk][ntrjpts];}
		
            }//<---End counting if 
	
         }//<---End ntrjpts loop
      
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
   //						Uniquely matching one WC Track to TPC Track
   //=======================================================================================================================
   
   // ### Keeping track of the number of matched tracks ###
   int nMatchedTracks = 0;
   
   // ### Setting the index for the track which is ###
   // ### uniquely matched to a wire chamber track ###
   bool MatchTPC_WVTrack[500] = {false};
   
   MatchWCTrackIndex[0] = 0;
   MatchWCTrackIndex[1] = 0;
   MatchWCTrackIndex[2] = 0;
   MatchWCTrackIndex[3] = 0;
   MatchWCTrackIndex[4] = 0;
   MatchWCTrackIndex[5] = 0;
   MatchWCTrackIndex[6] = 0;
   MatchWCTrackIndex[7] = 0;
   MatchWCTrackIndex[8] = 0;
   
   // ### Loop over tracks again for WCTrack / TPC Track Matching ###
   
   float FirstSpacePointIndex[500] = {0.};
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      
      // ### Only looking at tracks that pass low Z cut ###
      if(!PreLimTrackIndex[nTPCtrk]){continue;}
      
      float FirstSpacePointZ = 999;
      float FirstSpacePointY = 999;
      float FirstSpacePointX = 999;
      // #########################################################
      // ### Looping over the spacepoints for the prelim-track ###
      // #########################################################
      for(size_t ntrjpts = 0; ntrjpts < nTrajPoint[nTPCtrk]; ntrjpts++)
         {
	 
	 // ### Recording this tracks upstream most X, Y, Z location, ###
	 // ###       which is inside the fiducial boundary           ###
	 if(trjPt_Z[nTPCtrk][ntrjpts] < FirstSpacePointZ && trjPt_Y[nTPCtrk][ntrjpts] > YLowerFid && 
	    trjPt_Y[nTPCtrk][ntrjpts] < YUpperFid && trjPt_X[nTPCtrk][ntrjpts] > XLowerFid && 
	    trjPt_X[nTPCtrk][ntrjpts] < XUpperFid && trjPt_Z[nTPCtrk][ntrjpts] < UpperPartOfTPC)
	    {
	    FirstSpacePointIndex[nTPCtrk] = ntrjpts;
	    FirstSpacePointZ = trjPt_Z[nTPCtrk][ntrjpts];
	    FirstSpacePointY = trjPt_Y[nTPCtrk][ntrjpts];
	    FirstSpacePointX = trjPt_X[nTPCtrk][ntrjpts];
	    }//<---End Recording the tracks upstream X,Y,Z coordinate
         }//<---end ntrjpts loop
      
      // ### Skipping if this isn't the going to be matched ###
      if(FirstSpacePointZ == 999){continue;}
      // ### Now that I've gotten the earliest X, Y, Z spacepoint for this track ###
      // ###    I record the Delta X, Y, Z between this track and the WCTrack    ###
      
      
      // ### Variables for Delta WC and TPC tracks ###
      float DeltaX_WC_TPC_Track = 999;
      float DeltaY_WC_TPC_Track = 999;
            
      // ### Looping over the WCTrack ####
      for(size_t mwctrk = 0; mwctrk < nwctrks; mwctrk++)
         {
	 // ### Calculating the Delta X and Delta Y between WC track and TPC track ###
	 DeltaX_WC_TPC_Track = FirstSpacePointX - (wctrk_XFaceCoor[mwctrk]* 0.1);//<---Note: *0.1 to convert to cm
	 DeltaY_WC_TPC_Track = FirstSpacePointY - (wctrk_YFaceCoor[mwctrk]* 0.1);
	 
	 // ### Filling the Delta X and Delta Y between WC tracks and TPC Tracks ###
	 hdataDeltaWCTrkY->Fill(DeltaY_WC_TPC_Track);
	 hdataDeltaWCTrkX->Fill(DeltaX_WC_TPC_Track);
	 
	 // ### Counting the number of matched tracks ###
	 if( (DeltaX_WC_TPC_Track > DeltaXLowerBound && DeltaX_WC_TPC_Track < DeltaXUpperBound) &&
	     (DeltaY_WC_TPC_Track > DeltaYLowerBound && DeltaY_WC_TPC_Track < DeltaYUpperBound) )
	     {
	     // ### Setting the index of this track to true ###
	     MatchTPC_WVTrack[nTPCtrk] = true;
	     // ### Counting the matched tracks ###
	     nMatchedTracks++;
	     
	     // ### Setting the WCTrack Index = 1 if this WCTrack was matched ###
	     MatchWCTrackIndex[mwctrk] = 1;
	     }
	 
	 }//<---End mwctrk loop
      
      }//<---End nTPCtrk loop
    
    
       //=======================================================================================================================
   //				Calculating and cutting on the angle between the WC and TPC Track (alpha)
   //=======================================================================================================================
   
   
   // ###             Loop over tracks again only looking for              ###
   // ###                the one which matches the WC Track                ###
   // ###   and calculate the angle between the WCTrack and the TPC Track  ###
   
   // ###################################################
   // ### Vectors for angles between TPC and WC Track ###
   // ###################################################
   TVector3 z_hat(0,0,1);
   TVector3 p_hat_0;
    
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      // ### Skipping all the tracks which aren't well matched to a WC Track ###
      if(!MatchTPC_WVTrack[nTPCtrk]){continue;}
      if(!PreLimTrackIndex[nTPCtrk]){continue;}
      // ### Storing temp trajectory points to find the most upstream point ###
      float TempTrj_X = 900, TempTrj_Y = 900, TempTrj_Z = 900;
      
     
      // ### looping over the trajectory points for the track ###
      for(size_t nTrj = 0; nTrj < nTrajPoint[nTPCtrk]; nTrj++)
         {
         // ### Recording this tracks upstream most X, Y, Z location
	 if(FirstSpacePointIndex[nTPCtrk] == nTrj)
	    {
	    
	    // #################################################
	    // ### Saving the most upstream trajectory point ###
	    // #################################################
	    TempTrj_X = pHat0_X[nTPCtrk][nTrj];
	    TempTrj_Y = pHat0_Y[nTPCtrk][nTrj];
	    TempTrj_Z = pHat0_Z[nTPCtrk][nTrj];
	    	    
	    }//<---End Recording the tracks upstream X,Y,Z coordinate

         }//<---End nTrj loop
      
      // ### Setting the vector for the matched track ###
      // ###      most upstream trajectory point      ###
      p_hat_0.SetX(TempTrj_X);
      p_hat_0.SetY(TempTrj_Y);
      p_hat_0.SetZ(TempTrj_Z); //<--Note: since at this point we only have one unique match
      			       //         only having one entry should be fine
      
      }//<---End nTPCtrk loop


   // ### Grabbing the WCTrack Theta ###
   float wcTheta = wctrk_theta[0];
   hdataWCTheta->Fill(wcTheta* (180.0/3.141592654));
   
   // ### Grabbing the WCTRack Phi ###
   float wcPhi = wctrk_phi[0];
   hdataWCPhi->Fill(wcPhi* (180.0/3.141592654));
   
   // ### Calculating the Theta for the TPC Track ###
   float tpcTheta = acos(z_hat.Dot(p_hat_0)/p_hat_0.Mag());  
   hdataTPCTheta->Fill(tpcTheta* (180.0/3.141592654));
   
   // ### Using same convention as WCTrack to calculate phi ###
   float phi = 0;
   //Calculating phi (degeneracy elimination for the atan function)
   //----------------------------------------------------------------------------------------------
   if( p_hat_0.Y() > 0 && p_hat_0.X() > 0 ){ phi = atan(p_hat_0.Y()/p_hat_0.X()); }
   else if( p_hat_0.Y() > 0 && p_hat_0.X() < 0 ){ phi = atan(p_hat_0.Y()/p_hat_0.X())+3.141592654; }
   else if( p_hat_0.Y() < 0 && p_hat_0.X() < 0 ){ phi = atan(p_hat_0.Y()/p_hat_0.X())+3.141592654; }
   else if( p_hat_0.Y() < 0 && p_hat_0.X() > 0 ){ phi = atan(p_hat_0.Y()/p_hat_0.X())+6.28318; }
   else if( p_hat_0.Y() == 0 && p_hat_0.X() == 0 ){ phi = 0; }//defined by convention
   else if( p_hat_0.Y() == 0 )
      {
      if( p_hat_0.X() > 0 ){ phi = 0; }

      else{ phi = 3.141592654; }

      }
   else if( p_hat_0.X() == 0 )
      {
      if( p_hat_0.Y() > 0 ){ phi = 3.141592654/2; }
      else{ phi = 3.141592654*3/2; }

      }
   //----------------------------------------------------------------------------------------------
   
   // ### Using TPC Phi ###
   float tpcPhi = phi; 
   hdataTPCPhi->Fill(tpcPhi* (180.0/3.141592654));
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
   theUnitVector_TPCTrack.SetX(sin(tpcTheta)*cos(tpcPhi));
   theUnitVector_TPCTrack.SetY(sin(tpcTheta)*sin(tpcPhi));
   theUnitVector_TPCTrack.SetZ(cos(tpcTheta));
   
   // ###########################################################
   // ### Calculating the angle between WCTrack and TPC Track ###
   // ###########################################################
   float alpha = ( acos(theUnitVector_WCTrack.Dot(theUnitVector_TPCTrack)) )* (180.0/3.141592654);
   
   //std::cout<<"alpha = "<<alpha<<std::endl;
   hdataAlpha->Fill(alpha);
   
   
    
    // ### Filling the number of matched WC and TPC Tracks ###
    hdataNMatchTPCWCTrk->Fill(nMatchedTracks); 
    
    // ### Skipping this event if no WC track is matched ###
    // ###    OR if more than one WC track is matched    ###
    if( (nMatchedTracks < 1 || nMatchedTracks > 1) && alpha > alphaCut){continue;}
    
    //if(alpha > alphaCut){continue;}
   
    
    // ### Counting the number of events with ONE WC track matched ###
    nEvtsWCTrackMatch++;
    
    // ### Counting if we have a well matched stopping track ###
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






   // ---   First grab the tracks "initial" momentum which we take from ---
   // --- the momentum of the wire chamber track which has been matched ---
   // ---  and correct for the "typical" energy loss for a track in the ---
   // ---   argon between the cryostat and the front face of the TPC    ---
   
   
   float momentum = 999;
   
   //Vectors with calo info of the matched tpc track
   double Piondedx[1000]={0.};
   double Pionresrange[1000]={0.};
   double Pionpitchhit[1000]={0.};
   int nPionSpts = 0;
   double PionSumEnergy = 0;
   
   // ### Variables for determining the matching of the end point ###
   float TrackEndX = 999, TrackEndY = 999, TrackEndZ = 999;
   
   
   
   
   // ###################################################
   // ### Grabbing the Wire Chamber Track Information ###
   // ###################################################
   for(size_t mwctrk = 0; mwctrk < nwctrks; mwctrk++)
      {
      // ### Skip this WCTrack if it isn't the matched one ###
      if(MatchWCTrackIndex[mwctrk] < 1 || MatchWCTrackIndex[mwctrk] > 1){continue;}
      
      
      hdataWCTRKMomentum->Fill(wctrk_momentum[mwctrk]);
      momentum =wctrk_momentum[mwctrk];
      
      }//<---End mwctrk
   
   // ###   Calculating the initial Kinetic Energy    ###
   // ### KE = Energy - mass = (p^2 + m^2)^1/2 - mass ###
   float kineticEnergy = pow( (momentum*momentum) + (mass*mass) ,0.5) - mass;
   
   // ### The kinetic energy is that which we calculated ###
   // ###       minus the calculated energy loss         ###
   kineticEnergy -= entryTPCEnergyLoss;
   
   // ### Filling the initial kinetic energy plot ###
   hdataPionInitialMomentum->Fill(kineticEnergy);
   
   int howManyTracks = 0;
   
   // ### Loop over the tracks ###
   for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
      {
      // ### Skipping all the tracks which aren't well matched ###
      if(!MatchTPC_WVTrack[nTPCtrk]){continue;}
      
      // ### Recording the end-point of this track ###
      TrackEndX = trkendx[nTPCtrk];
      TrackEndY = trkendy[nTPCtrk];
      TrackEndZ = trkendz[nTPCtrk];
      
      hdataPionTrackEndX->Fill(TrackEndX);
      hdataPionTrackEndY->Fill(TrackEndY);
      hdataPionTrackEndZ->Fill(TrackEndZ);
      
      RecoLength = trklength[nTPCtrk];
      
      hRecoLength->Fill(RecoLength);
      
      // ### Recording the start-point of this track ###
      
      hdataPionTrackStartX->Fill(trkvtxx[nTPCtrk]);
      hdataPionTrackStartY->Fill(trkvtxy[nTPCtrk]);
      hdataPionTrackStartZ->Fill(trkvtxz[nTPCtrk]);
      
      // #################################################
      // ### If this tracks end point is at a boundary ###
      // ###   then tag this track as "through-going"  ###
      // #################################################
      if(TrackEndX < 1   || TrackEndX > 42.0 || TrackEndY > 19 ||
         TrackEndY < -19 || TrackEndZ > 89.0)
	 {ExitingTrack = true;}
      
      nPionSpts = 0;
      // ###################################################
      // ### Looping over the spacepoints for this track ###
      // ###################################################
      for(size_t nspts = 0; nspts < ntrkhits[nTPCtrk]; nspts++)
         {
	 // ###                 Note: Format for this variable is:             ###
	 // ### [trk number][plane 0 = induction, 1 = collection][spts number] ###
         Piondedx[nPionSpts]     = trkdedx[nTPCtrk][1][nspts] * 1.105; //<---dEdX correction factor
	 
	 // ### Putting in a fix in the case that the dE/dX is negative in this step ### 
	 // ###  then take the point before and the point after and average them
	 if(Piondedx[nPionSpts] < 0 && nspts < ntrkhits[nTPCtrk] && nspts > 0)
	    {Piondedx[nPionSpts] = ( (trkdedx[nTPCtrk][1][nspts - 1] + trkdedx[nTPCtrk][1][nspts + 1]) / 2);}
	 
	 // ### If this didn't fix it, then just put in a flat 2.4 MeV / cm fix ###
	 if(Piondedx[nPionSpts] < 0)
	    {
	    Piondedx[nPionSpts] = 2.4;
	    continue;}
	 
         Pionresrange[nPionSpts] = trkrr[nTPCtrk][1][nspts];
         Pionpitchhit[nPionSpts] = trkpitchhit[nTPCtrk][1][nspts];
         
	 PionSumEnergy = (Piondedx[nPionSpts] * Pionpitchhit[nPionSpts]) + PionSumEnergy;
	 
	 // ### Recording the dE/dX ###
	 hdataPiondEdX->Fill(Piondedx[nPionSpts]);
	 // ### Recording the residual range ###
	 hdataPionRR->Fill(Pionresrange[nPionSpts]);
	 // ### Recording the Pitch ###
	 hdataPionTrkPitch->Fill(Pionpitchhit[nPionSpts]);
	 
	 // ### Filling 2d dE/dX vs RR ###
	 hdataPiondEdXvsRR->Fill(Pionresrange[nPionSpts], Piondedx[nPionSpts]);
	 
	 nPionSpts++;
         }//<---End nspts loop
      howManyTracks++;
      }//<---End nTPCtrk loop    
   
   
   // ### Looping over the spacepoints to fill the histograms ###
   
for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
   {
   // ### Skipping all the tracks which aren't well matched ###
   if(!MatchTPC_WVTrack[nTPCtrk]){continue;}
   
   for(size_t npoints = 0; npoints < nPionSpts; npoints++)
      {
      // ### Filling the incidient histogram ###
      hdataPionIncidentKE->Fill(kineticEnergy);
      
      // ### Filling the interaction histogram for the last spt ###
      if(npoints == nPionSpts -1 && !ExitingTrack)
         {fPionInteractions->Fill(kineticEnergy);}
      
      //float energyLossInStep = Piondedx[npoints] * Pionresrange[npoints] * RecombinationFactor;
      float energyLossInStep = Piondedx[npoints] * Pionpitchhit[npoints];
      
      kineticEnergy -= energyLossInStep;
      
      
      }//<---End npoints loop
      
   }//<---End nTPC trk

      
   }//<---End jentry loop


// ===============================================================================================================
// 					SCALING FOR THE MUON CONTAMINATION
// ===============================================================================================================

hdataPionIncidentKE->Scale(MuonContaminationScaleFactor);



// ===============================================================================================================
//					MAKING THE CROSS-SECTION PLOT
// ===============================================================================================================

// ###################################################################
// #### Looping over the exiting bins to extract the cross-section ###
// ###################################################################
for( int iBin = 1; iBin <= fPionInteractions->GetNbinsX(); ++iBin )
   {
   
   std::cout<<std::endl;
   // ### If an incident bin is equal to zero then skip that bin ###
   if( hdataPionIncidentKE->GetBinContent(iBin) == 0 )continue; //Temporary fix to ensure that no Infinities are propagated to pad
   
   // ### Cross-section = (Exit Bins / Incident Bins) * (1/Density) * (1/slab width) ###
   float TempCrossSection = (fPionInteractions->GetBinContent(iBin)/hdataPionIncidentKE->GetBinContent(iBin)) * (1/number_density) * (1/slab_width);
   
   //std::cout<<"Cross-Section before conversion to barns = "<<TempCrossSection<<std::endl;
   
   float crossSection = TempCrossSection * (1/1e-28); //To put this into barns
   //std::cout<<"Cross-Section = "<<crossSection<<std::endl;
   
   
   fCrossSection->SetBinContent(iBin,crossSection);
   
   float numError = pow(fPionInteractions->GetBinContent(iBin),0.5);
   float num = fPionInteractions->GetBinContent(iBin);

   
   if(num == 0){num = 1;}
   float term1 = numError/num;
   //std::cout<<"term1 = "<<term1<<std::endl;
   
   float denomError = pow(hdataPionIncidentKE->GetBinContent(iBin),0.5);
   float denom = hdataPionIncidentKE->GetBinContent(iBin);
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
std::cout<<"########################################################################"<<std::endl;
std::cout<<"### Number of Events in AnaModule                                = "<<nTotalEvents<<" ###"<<std::endl;
std::cout<<"-------------------------------   Stage 0   ----------------------------"<<std::endl;
std::cout<<"### Number of Events w/ WC Track                                 = "<<nEvtsWCTrack<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ TOF > 0 ns and < 30 ns                   = "<<nEvtsTOF<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ Good TPC info (nHits > 0)		     = "<<nEvntsTPC<<" ###"<<std::endl;
std::cout<<"-------------------------------   Stage 1   ----------------------------"<<std::endl;
std::cout<<"### Number of Events w/ PID consistent with Pi/Mu                = "<<nEvtsPID<<" ###"<<std::endl;
std::cout<<"-------------------------------   Stage 2   ----------------------------"<<std::endl;
std::cout<<"### Number of Events w/ Trk Z < "<<FirstSpacePointZPos<<"                                = "<<nEvtsTrackZPos<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ < "<<nLowZTracksAllowed<<" tracks in the first "<<UpperPartOfTPC<<" cm of the TPC = "<<nLowZTrkEvents<<" ###"<<std::endl;
std::cout<<"### Number of Events w/ ONE WC Track Matched                     = "<<nEvtsWCTrackMatch<<" ###"<<std::endl;
std::cout<<"###              ( "<<DeltaXLowerBound<<" < Delta X < "<<DeltaXUpperBound<<" , "<<DeltaYLowerBound<<" < Delta Y < "<<DeltaYUpperBound<<" )             ###"<<std::endl;
std::cout<<"### Number of Events w/ angle alpha less the "<<alphaCut<<" degrees  = "<<nEventsPassingAlpha<<"         ###"<<std::endl;
std::cout<<"### Number of Events that are not Shower Like                        = "<<nNonShowerEvents<<std::endl;
std::cout<<"########################################################################"<<std::endl;
std::cout<<std::endl;

//std::cout<<"nTrajPoints = "<<nTrajPoints<<std::endl;
//std::cout<<"nSpacePoints = "<<nSpacePoints<<std::endl; 

// ===========================================================================================
// ============================  Writing out histograms to ROOT File =========================
// ===========================================================================================
// ###############################################
// ### Creating a file to output my histograms ###
// ###############################################
TFile myfile("PionXSection_histos.root","RECREATE");
hdataWCTrackExist->Write();
hdataWCTRKMomentum->Write();
hdataUpstreamZPos->Write();
hdataDeltaWCTrkX->Write();
hdataDeltaWCTrkY->Write();
hdataNMatchTPCWCTrk->Write();
hdataNTracksvsZpos->Write();
hdataAlpha->Write();
hdataPionInitialMomentum->Write();
hdataPiondEdX->Write();
hdataPionRR->Write();
hdataPionTrkPitch->Write();
hdataPiondEdXvsRR->Write();
hdataPionIncidentKE->Write();
fPionInteractions->Write();
fCrossSection->Write();
hdataTPCTheta->Write();
hdataTPCPhi->Write();
hdataWCTheta->Write();
hdataWCPhi->Write();
hdataWCTrackMomentumVSTOF->Write();
hdataPionTrackEndX->Write();
hdataPionTrackEndY->Write();
hdataPionTrackEndZ->Write();
hdataPionTrackStartX->Write();
hdataPionTrackStartY->Write();
hdataPionTrackStartZ->Write();
hRecoLength->Write();





// ===================================================================================================================
// =====================================   Drawing Histograms   ======================================================
// ===================================================================================================================


// ### Making a TCanvas ###
TCanvas *c2= new TCanvas("c2","WC Momentum");
c2->SetTicks();
c2->SetLogy();
c2->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataWCTRKMomentum->SetLineColor(kBlack);
hdataWCTRKMomentum->SetLineStyle(0);
hdataWCTRKMomentum->SetLineWidth(3);
hdataWCTRKMomentum->SetMarkerStyle(8);
hdataWCTRKMomentum->SetMarkerSize(0.9);

hdataWCTRKMomentum->Draw();  

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
TLegend *leg2 = new TLegend();
leg2 = new TLegend(0.58,0.65,0.88,0.88);
leg2->SetTextSize(0.04);
leg2->SetTextAlign(12);
leg2->SetFillColor(kWhite);
leg2->SetLineColor(kWhite);
leg2->SetShadowColor(kWhite);
leg2->SetHeader("Good Run List");
leg2->AddEntry(hdataWCTRKMomentum, "WC Track Momentum");
leg2->Draw();






// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c3= new TCanvas("c3","Delta Y WC/TPC Track");
c3->SetTicks();
//c3->SetLogy();
c3->SetFillColor(kWhite);  

// ### Drawing Histo ###

hdataDeltaWCTrkY->SetLineColor(kBlack);
hdataDeltaWCTrkY->SetLineStyle(0);
hdataDeltaWCTrkY->SetLineWidth(3);
hdataDeltaWCTrkY->SetMarkerStyle(8);
hdataDeltaWCTrkY->SetMarkerSize(0.9); 
hdataDeltaWCTrkY->Draw();  

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
TLegend *leg3 = new TLegend();
leg3 = new TLegend(0.58,0.65,0.88,0.88);
leg3->SetTextSize(0.04);
leg3->SetTextAlign(12);
leg3->SetFillColor(kWhite);
leg3->SetLineColor(kWhite);
leg3->SetShadowColor(kWhite);
leg3->SetHeader("Good Run List");
leg3->AddEntry(hdataDeltaWCTrkY, "#Delta Y WC/TPC Track");
leg3->Draw();








// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c4= new TCanvas("c4","Delta X WC/TPC Track");
c4->SetTicks();
//c4->SetLogy();
c4->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataDeltaWCTrkX->SetLineColor(kBlack);
hdataDeltaWCTrkX->SetLineStyle(0);
hdataDeltaWCTrkX->SetLineWidth(3);
hdataDeltaWCTrkX->SetMarkerStyle(8);
hdataDeltaWCTrkX->SetMarkerSize(0.9);

hdataDeltaWCTrkX->Draw();  

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
TLegend *leg4 = new TLegend();
leg4 = new TLegend(0.58,0.65,0.88,0.88);
leg4->SetTextSize(0.04);
leg4->SetTextAlign(12);
leg4->SetFillColor(kWhite);
leg4->SetLineColor(kWhite);
leg4->SetShadowColor(kWhite);
leg4->SetHeader("Open Box Data");
leg4->AddEntry(hdataDeltaWCTrkX, "#Delta X WC/TPC Track");
leg4->Draw();








// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c5= new TCanvas("c5","Number of matched WC/TPC Track");
c5->SetTicks();
//c5->SetLogy();
c5->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataNMatchTPCWCTrk->SetLineColor(kBlack);
hdataNMatchTPCWCTrk->SetLineStyle(0);
hdataNMatchTPCWCTrk->SetLineWidth(3);
hdataNMatchTPCWCTrk->SetMarkerStyle(8);
hdataNMatchTPCWCTrk->SetMarkerSize(0.9);

hdataNMatchTPCWCTrk->Draw();  

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
TLegend *leg5 = new TLegend();
leg5 = new TLegend(0.58,0.65,0.88,0.88);
leg5->SetTextSize(0.04);
leg5->SetTextAlign(12);
leg5->SetFillColor(kWhite);
leg5->SetLineColor(kWhite);
leg5->SetShadowColor(kWhite);
leg5->SetHeader("Open Box Data");
leg5->AddEntry(hdataNMatchTPCWCTrk, "Number of Matched WC/TPC Track");
leg5->Draw();




// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c6= new TCanvas("c6","Pion dE/dX");
c6->SetTicks();
//c5->SetLogy();
c6->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataPiondEdX->SetLineColor(kBlack);
hdataPiondEdX->SetLineStyle(0);
hdataPiondEdX->SetLineWidth(3);
hdataPiondEdX->SetMarkerStyle(8);
hdataPiondEdX->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataPiondEdX->GetXaxis()->SetTitle("dE/dX (MeV/cm)");
hdataPiondEdX->GetXaxis()->CenterTitle();

hdataPiondEdX->GetYaxis()->SetTitle("");
hdataPiondEdX->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataPiondEdX->Draw(); 

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
TLegend *leg6 = new TLegend();
leg6 = new TLegend(0.58,0.65,0.88,0.88);
leg6->SetTextSize(0.04);
leg6->SetTextAlign(12);
leg6->SetFillColor(kWhite);
leg6->SetLineColor(kWhite);
leg6->SetShadowColor(kWhite);
leg6->SetHeader("Open Box Data");
leg6->AddEntry(hdataPiondEdX, "dE/dX");
leg6->Draw();



// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c7= new TCanvas("c7","Pion Residual Range");
c7->SetTicks();
//c7->SetLogy();
c7->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataPionRR->SetLineColor(kBlack);
hdataPionRR->SetLineStyle(0);
hdataPionRR->SetLineWidth(3);
hdataPionRR->SetMarkerStyle(8);
hdataPionRR->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataPionRR->GetXaxis()->SetTitle("Residual Range (cm)");
hdataPionRR->GetXaxis()->CenterTitle();

hdataPionRR->GetYaxis()->SetTitle("Events / 0.5 cm");
hdataPionRR->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataPionRR->Draw();  

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
TLegend *leg7 = new TLegend();
leg7 = new TLegend(0.58,0.65,0.88,0.88);
leg7->SetTextSize(0.04);
leg7->SetTextAlign(12);
leg7->SetFillColor(kWhite);
leg7->SetLineColor(kWhite);
leg7->SetShadowColor(kWhite);
leg7->SetHeader("Open Box Data");
leg7->AddEntry(hdataDeltaWCTrkY, "Residual Range");
leg7->Draw();




// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c8= new TCanvas("c8","dEdX vs RR");
c8->SetTicks();
//c8->SetLogy();
c8->SetFillColor(kWhite);  

// ### Labeling the axis ###
hdataPiondEdXvsRR->GetXaxis()->SetTitle("Residual Range (cm)");
hdataPiondEdXvsRR->GetXaxis()->CenterTitle();

hdataPiondEdXvsRR->GetYaxis()->SetTitle("dE/dX (MeV/cm)");
hdataPiondEdXvsRR->GetYaxis()->CenterTitle();

// ### Drawing Histo ### 
hdataPiondEdXvsRR->Draw("colz");  

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
TLegend *leg8 = new TLegend();
leg8 = new TLegend(0.58,0.65,0.88,0.88);
leg8->SetTextSize(0.04);
leg8->SetTextAlign(12);
leg8->SetFillColor(kWhite);
leg8->SetLineColor(kWhite);
leg8->SetShadowColor(kWhite);
leg8->SetHeader("Open Box Data");
leg8->AddEntry(hdataPiondEdXvsRR, "");
leg8->Draw();





// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c9= new TCanvas("c9","Incident Kinetic Energy");
c9->SetTicks();
//c9->SetLogy();
c9->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataPionIncidentKE->SetLineColor(kBlack);
hdataPionIncidentKE->SetLineStyle(0);
hdataPionIncidentKE->SetLineWidth(3);
hdataPionIncidentKE->SetMarkerStyle(8);
hdataPionIncidentKE->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataPionIncidentKE->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
hdataPionIncidentKE->GetXaxis()->CenterTitle();

hdataPionIncidentKE->GetYaxis()->SetTitle("Events / 50 MeV");
hdataPionIncidentKE->GetYaxis()->CenterTitle();

hdataPionIncidentKE->Draw("");  


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
TLegend *leg9 = new TLegend();
leg9 = new TLegend(0.58,0.65,0.88,0.88);
leg9->SetTextSize(0.04);
leg9->SetTextAlign(12);
leg9->SetFillColor(kWhite);
leg9->SetLineColor(kWhite);
leg9->SetShadowColor(kWhite);
leg9->SetHeader("Open Box Data");
leg9->AddEntry(hdataPionIncidentKE, "");
leg9->Draw();


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c25= new TCanvas("c25","Interaction KE");
c25->SetTicks();
c25->SetFillColor(kWhite);

// ### Formatting the histograms ###
fPionInteractions->SetLineColor(kBlack);
fPionInteractions->SetLineStyle(0);
fPionInteractions->SetLineWidth(3);
fPionInteractions->SetMarkerStyle(8);
fPionInteractions->SetMarkerSize(0.9);

// ### Labeling the axis ###
fPionInteractions->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
fPionInteractions->GetXaxis()->CenterTitle();

fPionInteractions->GetYaxis()->SetTitle("Events / 50 MeV");
fPionInteractions->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
fPionInteractions->Draw();

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
leg->SetHeader("Open Box Data");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();





// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c10= new TCanvas("c10","Track Pitch");
c10->SetTicks();
//c10->SetLogy();
c10->SetFillColor(kWhite);  

// ### Drawing Histo ### 
hdataPionTrkPitch->SetLineColor(kBlack);
hdataPionTrkPitch->SetLineStyle(0);
hdataPionTrkPitch->SetLineWidth(3);
hdataPionTrkPitch->SetMarkerStyle(8);
hdataPionTrkPitch->SetMarkerSize(0.9);

hdataPionTrkPitch->Draw("");  

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
//leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Open Box Data");
//leg->AddEntry(hdataPionTrkPitch, "Track Pitch");

leg->Draw();






// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c14= new TCanvas("c14","WC Track Events");
c14->SetTicks();
//c9->SetLogy();
c14->SetFillColor(kWhite); 

// ### Drawing Histo ### 
hdataWCTrackExist->SetLineColor(kBlack);
hdataWCTrackExist->SetLineStyle(0);
hdataWCTrackExist->SetLineWidth(3);
hdataWCTrackExist->SetMarkerStyle(8);
hdataWCTrackExist->SetMarkerSize(0.9);

hdataWCTrackExist->Draw();

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
//leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Good Run List");
leg->AddEntry(hdataWCTrackExist, "Events w/ WC Track");

leg->Draw();

// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c15= new TCanvas("c15","Most upstream Spacepoint");
c15->SetTicks();
//c9->SetLogy();
c15->SetFillColor(kWhite); 

hdataUpstreamZPos->SetLineColor(kBlack);
hdataUpstreamZPos->SetLineStyle(0);
hdataUpstreamZPos->SetLineWidth(3);
hdataUpstreamZPos->SetMarkerStyle(8);
hdataUpstreamZPos->SetMarkerSize(0.9);

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
//leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Good Run List");
leg->AddEntry(hdataWCTrackExist, "Most Upstream Spacepoint");

leg->Draw();



// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c16= new TCanvas("c16","Number of tracks vs Z");
c16->SetTicks();
//c9->SetLogy();
c16->SetFillColor(kWhite); 

hdataNTracksvsZpos->SetLineColor(kBlack);
hdataNTracksvsZpos->SetLineStyle(0);
hdataNTracksvsZpos->SetLineWidth(3);
hdataNTracksvsZpos->SetMarkerStyle(8);
hdataNTracksvsZpos->SetMarkerSize(0.9);
hdataNTracksvsZpos->Draw("colz");

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
//leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
leg->SetHeader("Good Run List");
leg->AddEntry(hdataWCTrackExist, "Most Upstream Spacepoint");


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------
// ### Making a TCanvas ###
TCanvas *c17= new TCanvas("c17","#alpha between WC and TPC Track");
c17->SetTicks();
//c9->SetLogy();
c17->SetFillColor(kWhite); 

hdataAlpha->SetLineColor(kBlack);
hdataAlpha->SetLineStyle(0);
hdataAlpha->SetLineWidth(3);
hdataAlpha->SetMarkerStyle(8);
hdataAlpha->SetMarkerSize(0.9);

hdataAlpha->Draw();   
 
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
leg->SetHeader("Open Box Data");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();


// ---------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------


// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c27= new TCanvas("c27","TOF vs Pz");
c27->SetTicks();
c27->SetFillColor(kWhite);

// ### Formatting the histograms ###
hdataWCTrackMomentumVSTOF->SetLineColor(kBlack);
hdataWCTrackMomentumVSTOF->SetLineStyle(0);
hdataWCTrackMomentumVSTOF->SetLineWidth(3);
hdataWCTrackMomentumVSTOF->SetMarkerStyle(8);
hdataWCTrackMomentumVSTOF->SetMarkerSize(0.9);

// ### Labeling the axis ###
hdataWCTrackMomentumVSTOF->GetXaxis()->SetTitle("WC Track P_{z} (MeV)");
hdataWCTrackMomentumVSTOF->GetXaxis()->CenterTitle();

hdataWCTrackMomentumVSTOF->GetYaxis()->SetTitle("TOF (ns)");
hdataWCTrackMomentumVSTOF->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hdataWCTrackMomentumVSTOF->Draw("colz");

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
leg->SetHeader("Open Box Data");
//leg->AddEntry(hMCPrimaryStartX, "X_{0}");
leg->Draw();



   
   
}//<---End Loop() Function




//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
// 				### Function for plotting the Low Z track location ###
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
void DataAnalysis::LowZCut()
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
   for(size_t ntrjpts = 0; ntrjpts < nTrajPoint[nTPCtrk]; ntrjpts++)
      {
      
      // ### Putting the number of tracks w/ spts in this area of Z ###
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 1 && !trackwZLT1)
	 	{trkwZLT1++; trackwZLT1 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 2 && !trackwZLT2)
	 	{trkwZLT2++; trackwZLT2 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 3 && !trackwZLT3)
	 	{trkwZLT3++; trackwZLT3 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 4 && !trackwZLT4)
	 	{trkwZLT4++; trackwZLT4 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 5 && !trackwZLT5)
	 	{trkwZLT5++; trackwZLT5 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 6 && !trackwZLT6)
	 	{trkwZLT6++; trackwZLT6 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 7 && !trackwZLT7)
	 	{trkwZLT7++; trackwZLT7 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 8 && !trackwZLT8)
	 	{trkwZLT8++; trackwZLT8 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 9 && !trackwZLT9)
	 	{trkwZLT9++; trackwZLT9 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 10 && !trackwZLT10)
	 	{trkwZLT10++; trackwZLT10 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 11 && !trackwZLT11)
	 	{trkwZLT11++; trackwZLT11 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 12 && !trackwZLT12)
	 	{trkwZLT12++; trackwZLT12 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 13 && !trackwZLT13)
	 	{trkwZLT13++; trackwZLT13 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 14 && !trackwZLT14)
	 	{trkwZLT14++; trackwZLT14 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 15 && !trackwZLT15)
	 	{trkwZLT15++; trackwZLT15 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 16 && !trackwZLT16)
	 	{trkwZLT16++; trackwZLT16 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 17 && !trackwZLT17)
	 	{trkwZLT17++; trackwZLT17 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 18 && !trackwZLT18)
	 	{trkwZLT18++; trackwZLT18 = true;}		
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 19 && !trackwZLT19)
	 	{trkwZLT19++; trackwZLT19 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 20 && !trackwZLT20)
	 	{trkwZLT20++; trackwZLT20 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 21 && !trackwZLT21)
	 	{trkwZLT21++; trackwZLT21 = true;} 
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 22 && !trackwZLT22)
	 	{trkwZLT22++; trackwZLT22 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 23 && !trackwZLT23)
	 	{trkwZLT23++; trackwZLT23 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 24 && !trackwZLT24)
	 	{trkwZLT24++; trackwZLT24 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 25 && !trackwZLT25)
	 	{trkwZLT25++; trackwZLT25 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 26 && !trackwZLT26)
	 	{trkwZLT26++; trackwZLT26 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 27 && !trackwZLT27)
	 	{trkwZLT27++; trackwZLT27 = true;}		
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 28 && !trackwZLT28)
	 	{trkwZLT28++; trackwZLT28 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 29 && !trackwZLT29)
	 	{trkwZLT29++; trackwZLT29 = true;}
	 if(trjPt_Z[nTPCtrk][ntrjpts] < 30 && !trackwZLT30)
	 	{trkwZLT30++; trackwZLT30 = true;}
   
   
      }//<---End ntrjpts
   
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

 

}//<----End LowZCut() Function
