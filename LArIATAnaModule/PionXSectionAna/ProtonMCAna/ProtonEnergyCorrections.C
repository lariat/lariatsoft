#define ProtonEnergyCorrections_cxx
#include "ProtonEnergyCorrections.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ===================================================================================================================
// ====================================       PUT HISTOGRAMS HERE           ==========================================
// ===================================================================================================================

/////////////////////////////////// Energy Loss in the upstream region of the beamline ///////////////////////
TH1D *hMCELossUpstream = new TH1D("hMCELossUpstream", "Energy loss prior to entering the TPC", 1000, 0, 1000);

/////////////////////////////////// Energy Loss in the TPC  ///////////////////////
TH1D *hMCELossInTPC = new TH1D("hMCELossInTPC", "Energy loss inside the TPC", 1000, 0, 1000);


void ProtonEnergyCorrections::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ProtonEnergyCorrections.C
//      Root > ProtonEnergyCorrections t
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
// ##########################################################
// ### Putting in some counters for event reduction table ###
// ##########################################################
int nTotalEvents = 0;


float particle_mass = 938.28 ;

// #################################
// ### Variables for Energy Loss ###
// #################################

float EnergyLossOutsideTPC = 0;
float EnergyLossInsideTPC = 0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;


// ###########################
// ### Looping over events ###
// ###########################
for (Long64_t jentry=0; jentry<nentries;jentry++) 
//for (Long64_t jentry=0; jentry<2000;jentry++)
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
   // #######################################
   // ### Loop over all the G4  particles ###
   // #######################################
   for (int iG4 = 0; iG4 < geant_list_size; iG4++)
      {
      
      EnergyLossOutsideTPC = 0;
      EnergyLossInsideTPC = 0;
      
      // #########################################################
      // ### Skipping any particle which doesn't enter the TPC ###
      // #########################################################
      if(EndPointz[iG4] < 0){continue;}
      
      // #############################################
      // ### Loop over all the primary true points ###
      // #############################################
      if(process_primary[iG4] == 1)
         {
	 //std::cout<<"NTrTrajPts[iG4] = "<<NTrTrajPts[iG4]<<std::endl;
	 for(int iPriTrj = 1; iPriTrj < NTrTrajPts[iG4]; iPriTrj++)
	    {
	    
	    
	    //if(iPriTrj > 999){continue;}
	    
	    // ### Only looking at points which are upstream of the TPC ###
	    if(MidPosZ[iG4][iPriTrj] < 0)
	       {
	    
	       float Momentum_Point1 = (sqrt((MidPx[iG4][iPriTrj-1]*MidPx[iG4][iPriTrj-1]) + 
	                               (MidPy[iG4][iPriTrj-1]*MidPy[iG4][iPriTrj-1]) +
				       (MidPz[iG4][iPriTrj-1]*MidPz[iG4][iPriTrj-1])))*1000;
				       
	       float Momentum_Point2 = (sqrt((MidPx[iG4][iPriTrj]*MidPx[iG4][iPriTrj]) + 
	                               (MidPy[iG4][iPriTrj]*MidPy[iG4][iPriTrj]) +
				       (MidPz[iG4][iPriTrj]*MidPz[iG4][iPriTrj])))*1000;
				       
	       float Energy_Point1 = sqrt( (Momentum_Point1*Momentum_Point1) + (particle_mass*particle_mass)  );
	       
	       float Energy_Point2 = sqrt( (Momentum_Point2*Momentum_Point2) + (particle_mass*particle_mass)  );
	       
	       EnergyLossOutsideTPC +=  Energy_Point1 - Energy_Point2;
	       }//<---End Looking at energy loss upstream of TPC
	       
	    // ### Only looking at points which are inside of the TPC ###
	    if(MidPosZ[iG4][iPriTrj] > 0 && MidPosZ[iG4][iPriTrj] < 90 &&
	       MidPosX[iG4][iPriTrj] > 0 && MidPosX[iG4][iPriTrj] < 43 &&
	       MidPosY[iG4][iPriTrj] > -20 && MidPosY[iG4][iPriTrj] < 20 )
	       {
	       
	       float Momentum_Point1 = (sqrt((MidPx[iG4][iPriTrj-1]*MidPx[iG4][iPriTrj-1]) + 
	                               (MidPy[iG4][iPriTrj-1]*MidPy[iG4][iPriTrj-1]) +
				       (MidPz[iG4][iPriTrj-1]*MidPz[iG4][iPriTrj-1])))*1000;
				       
	       float Momentum_Point2 = (sqrt((MidPx[iG4][iPriTrj]*MidPx[iG4][iPriTrj]) + 
	                               (MidPy[iG4][iPriTrj]*MidPy[iG4][iPriTrj]) +
				       (MidPz[iG4][iPriTrj]*MidPz[iG4][iPriTrj])))*1000;
				       
	       float Energy_Point1 = sqrt( (Momentum_Point1*Momentum_Point1) + (particle_mass*particle_mass)  );
	       
	       //std::cout<<"Energy_Point1 = "<<Energy_Point1<<std::endl;
	       float Energy_Point2 = sqrt( (Momentum_Point2*Momentum_Point2) + (particle_mass*particle_mass)  );
	       //std::cout<<"Energy_Point2 = "<<Energy_Point2<<std::endl;
	       
	       //std::cout<<std::endl;
	       EnergyLossInsideTPC +=  Energy_Point1 - Energy_Point2;
	   
	       }//<---End Looking at energy loss upstream of TPC
	    
	    
	    
	    
	    }//<---iPriTrj loop
	    
	 hMCELossUpstream->Fill(EnergyLossOutsideTPC);
	 hMCELossInTPC->Fill(EnergyLossInsideTPC);
	 }//<---Only looking at primary particles
      
      
      
      }//<---End iG4Pri loop
   
   
   }//<---End jentry loop

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c01= new TCanvas("c01","Energy Loss Upstream");
c01->SetTicks();
c01->SetFillColor(kWhite);

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
leg->Draw();

//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c02= new TCanvas("c02","Energy Loss Inside TPC");
c02->SetTicks();
c02->SetFillColor(kWhite);

// ### Formatting the histograms ###
// ### Formatting the histograms ###
hMCELossInTPC->SetLineColor(kBlue);
hMCELossInTPC->SetLineStyle(0);
hMCELossInTPC->SetLineWidth(3);
hMCELossInTPC->SetMarkerStyle(8);
hMCELossInTPC->SetMarkerSize(0.9);

// ### Labeling the axis ###
hMCELossInTPC->GetXaxis()->SetTitle("Energy Loss inside the TPC (MeV)");
hMCELossInTPC->GetXaxis()->CenterTitle();

hMCELossInTPC->GetYaxis()->SetTitle("Events / MeV");
hMCELossInTPC->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hMCELossInTPC->Draw();

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


}//<---End Loop()
