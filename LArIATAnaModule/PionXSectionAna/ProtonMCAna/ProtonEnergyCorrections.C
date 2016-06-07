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
TH1D *hMCELossInTPC = new TH1D("hMCELossInTPC", "Energy loss inside the TPC", 110, -100, 1000);

/////////////////////////////////// Initial Kinetic Energy /////////////////////////
TH1D *hInitialKE = new TH1D("hInitialKE", "Initial Kinetic Energy", 110, -100, 1000);

/////////////////////////////////// Initial Kinetic Energy @ the TPC /////////////////////////
TH1D *hInitialKEAtTPC = new TH1D("hInitialKEAtTPC", "Initial Kinetic Energy at the TPC", 110, -100, 1000);

/////////////////////////////////// Final Kinetic Energy in the TPC /////////////////////////
TH1D *hFinalKEInTPC = new TH1D("hFinalKEInTPC", "Final Kinetic Energy in the TPC", 110, -100, 1000);

/////////////////////////////////// E Loss upstream of the TPC /////////////////////////
TH2D *hELossXvsY = new TH2D("hELossXvsY", "Energy Loss X vs Y", 86, 0, 43, 80, -20, 20);

/////////////////////////////////// E Loss upstream of the Flux TPC /////////////////////////
TH2D *hELossXvsYFlux = new TH2D("hELossXvsYFlux", "Energy Loss X vs Y", 86, 0, 43, 80, -20, 20);

/////////////////////////////////// Divided E Loss /////////////////////////
TH2D *hELossXvsYDivide = new TH2D("hELossXvsYDivide", "Energy Loss X vs Y", 86, 0, 43, 80, -20, 20);


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
int nTotalEvents = 0, nEventsReachTPC = 0, nParticlesStop = 0;



bool EventsReachTPC = false;
bool EventsWhereParticleStops = false;

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
   
   // ### Assume the event doesn't reach the TPC ###
   EventsReachTPC = false;
   
   // ### Assume the particle does stop ###
   EventsWhereParticleStops = true;
      
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
      
      
      float FinalKEInTPC = 0;
      
      // #########################################################
      // ### Skipping any particle which doesn't enter the TPC ###
      // #########################################################
      if(EndPointz[iG4] < 0){continue;}
      
      // ### If we are here, the event must pass through the TPC ###
      EventsReachTPC = true;
      
      // #######################################################
      // ### Skipping the event if the particle doesn't stop ###
      // #######################################################
      if(EndPointz[iG4] < 0 || EndPointz[iG4] > 90 || EndPointx[iG4] < 0 || 
         EndPointx[iG4] > 43 || EndPointy[iG4] < -20 || 
	 EndPointy[iG4] > 20 || NumberDaughters[iG4] > 0) {EventsWhereParticleStops = false; continue; }
      
      
      float g4Primary_ProjX0 = 0.0;
      float g4Primary_ProjY0 = 0.0;
      float g4Primary_ProjZ0 = 0.0;      
      // #############################################
      // ### Loop over all the primary true points ###
      // #############################################
      if(process_primary[iG4] == 1 && NumberDaughters[iG4] == 0)
         {
	 
	 // ############################################
	 // ### Calculate the initial Kinetic Energy ###
	 // ############################################
	 
	 float momentum = sqrt( (Px[iG4] * Px[iG4]) + (Py[iG4] * Py[iG4]) + (Pz[iG4] * Pz[iG4]) ) * 1000;
	 
	 float InitialKE = sqrt( (momentum*momentum) + (particle_mass*particle_mass)) - particle_mass;
	 
	 hInitialKE->Fill(InitialKE);
	 
	 for(int iPriTrj = 1; iPriTrj < NTrTrajPts[iG4]; iPriTrj++)
	    {
	    	    
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
	    if(MidPosZ[iG4][iPriTrj] >= 0 && MidPosZ[iG4][iPriTrj] < 90 &&
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
	 
	 // ###############################################################
	 // ### Calculating the initial KE at the front face of the TPC ###
	 // ###############################################################
	 float InitialKEAtTPC = InitialKE - EnergyLossOutsideTPC;
	 
	 hInitialKEAtTPC->Fill(InitialKEAtTPC);
	 
	 // ###############################################
	 // ### Calculating the Final KE inside the TPC ###
	 // ###############################################
	 FinalKEInTPC = InitialKE - EnergyLossOutsideTPC - EnergyLossInsideTPC;
	 hFinalKEInTPC->Fill(FinalKEInTPC);
	 
	 
	 // ------------------------------------------------------------------------------------
	 // ---------------        Extrapolate the X, Y, Z position of the primary         -----
	 // ---------------     if it started upstream of the front face of the TPC        -----
	 // ------------------------------------------------------------------------------------
	 
	 double b1 = StartPointz[iG4] - StartPointx[iG4]*Pz[iG4]/Px[iG4];
	 double b2 = StartPointz[iG4] - StartPointy[iG4]*Pz[iG4]/Py[iG4];
	 
	 g4Primary_ProjX0 = -b1*Px[iG4]/Pz[iG4];
	 g4Primary_ProjY0 = -b2*Py[iG4]/Pz[iG4];
	 g4Primary_ProjZ0 = 0.0;
	 
	 hELossXvsY->Fill(g4Primary_ProjX0, g4Primary_ProjY0, EnergyLossOutsideTPC);
	 
	 hELossXvsYFlux->Fill(g4Primary_ProjX0, g4Primary_ProjY0);
	 
	 }//<---Only looking at primary particles
      
      float P = sqrt( (Px[iG4]*Px[iG4]) + (Py[iG4]*Py[iG4]) + (Pz[iG4]*Pz[iG4])) * 1000;
      //if(FinalKEInTPC == 0){std::cout<<"Run "<<run<<", SubRun "<<subrun<<", Event "<<event<<" P "<<P<<std::endl;}
      }//<---End iG4Pri loop
   
   
   // ### Counting the events which pass through the TPC ###
   if(EventsReachTPC){nEventsReachTPC++;}
   if(EventsWhereParticleStops){nParticlesStop++;}
   
   }//<---End jentry loop


std::cout<<"========================================================================================"<<std::endl;
std::cout<<"Total Number of Events 				"<<nTotalEvents<<std::endl;
std::cout<<"Particles reach the TPC				"<<nEventsReachTPC<<std::endl;
std::cout<<"Particle stops inside the TPC			"<<nParticlesStop<<std::endl;
std::cout<<"========================================================================================"<<std::endl;


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


//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c03= new TCanvas("c03","Initial KE");
c03->SetTicks();
c03->SetFillColor(kWhite);

// ### Formatting the histograms ###
// ### Formatting the histograms ###
hInitialKE->SetLineColor(kBlue);
hInitialKE->SetLineStyle(0);
hInitialKE->SetLineWidth(3);
hInitialKE->SetMarkerStyle(8);
hInitialKE->SetMarkerSize(0.9);

// ### Labeling the axis ###
hInitialKE->GetXaxis()->SetTitle("K.E._{inital} (MeV)");
hInitialKE->GetXaxis()->CenterTitle();

hInitialKE->GetYaxis()->SetTitle("Events / 10 MeV");
hInitialKE->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hInitialKE->Draw();

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



//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c04= new TCanvas("c04","Initial KE at TPC");
c04->SetTicks();
c04->SetFillColor(kWhite);

// ### Formatting the histograms ###
hInitialKEAtTPC->SetLineColor(kBlue);
hInitialKEAtTPC->SetLineStyle(0);
hInitialKEAtTPC->SetLineWidth(3);
hInitialKEAtTPC->SetMarkerStyle(8);
hInitialKEAtTPC->SetMarkerSize(0.9);

// ### Labeling the axis ###
hInitialKEAtTPC->GetXaxis()->SetTitle("K.E._{inital} @ TPC (MeV)");
hInitialKEAtTPC->GetXaxis()->CenterTitle();

hInitialKEAtTPC->GetYaxis()->SetTitle("Events / 10 MeV");
hInitialKEAtTPC->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hInitialKEAtTPC->Draw();

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


//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c05= new TCanvas("c05","Final KE in TPC");
c05->SetTicks();
c05->SetFillColor(kWhite);

// ### Formatting the histograms ###
hFinalKEInTPC->SetLineColor(kBlue);
hFinalKEInTPC->SetLineStyle(0);
hFinalKEInTPC->SetLineWidth(3);
hFinalKEInTPC->SetMarkerStyle(8);
hFinalKEInTPC->SetMarkerSize(0.9);

// ### Labeling the axis ###
hFinalKEInTPC->GetXaxis()->SetTitle("K.E._{Final} (MeV)");
hFinalKEInTPC->GetXaxis()->CenterTitle();

hFinalKEInTPC->GetYaxis()->SetTitle("Events / 10 MeV");
hFinalKEInTPC->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hFinalKEInTPC->Draw();

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



//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c06= new TCanvas("c06","Energy Loss vs X and Y");
c06->SetTicks();
c06->SetFillColor(kWhite);

// ### Formatting the histograms ###

// ### Labeling the axis ###
hELossXvsY->GetXaxis()->SetTitle("X Postion at TPC Face (cm)");
hELossXvsY->GetXaxis()->CenterTitle();

hELossXvsY->GetYaxis()->SetTitle("Y Position at TPC Face (cm)");
hELossXvsY->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hELossXvsY->Draw("colz");

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
//leg->SetTextSize(0.04);
//leg->SetTextAlign(12);
//leg->SetFillColor(kWhite);
//leg->SetLineColor(kWhite);
//leg->SetShadowColor(kWhite);
//leg->SetHeader("Proton MC");
//leg->Draw();

//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c07= new TCanvas("c07","Energy Loss vs X and Y Flux");
c07->SetTicks();
c07->SetFillColor(kWhite);

// ### Formatting the histograms ###

// ### Labeling the axis ###
hELossXvsYFlux->GetXaxis()->SetTitle("X Postion at TPC Face (cm)");
hELossXvsYFlux->GetXaxis()->CenterTitle();

hELossXvsYFlux->GetYaxis()->SetTitle("Y Position at TPC Face (cm)");
hELossXvsYFlux->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hELossXvsYFlux->Draw("colz");

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
//leg->SetTextSize(0.04);
//leg->SetTextAlign(12);
//leg->SetFillColor(kWhite);
//leg->SetLineColor(kWhite);
//leg->SetShadowColor(kWhite);
//leg->SetHeader("Proton MC");
//leg->Draw();

//----------------------------------------------------------------------------------------------

// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c08= new TCanvas("c08","Energy Loss vs X and Y Flux");
c08->SetTicks();
c08->SetFillColor(kWhite);

// ### Formatting the histograms ###

hELossXvsYDivide->Divide(hELossXvsY, hELossXvsYFlux);

// ### Labeling the axis ###
hELossXvsYDivide->GetXaxis()->SetTitle("X Postion at TPC Face (cm)");
hELossXvsYDivide->GetXaxis()->CenterTitle();

hELossXvsYDivide->GetYaxis()->SetTitle("Y Position at TPC Face (cm)");
hELossXvsYDivide->GetYaxis()->CenterTitle();

// ### Drawing the histogram ### 
hELossXvsYDivide->Draw("colz");

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
//leg->SetTextSize(0.04);
//leg->SetTextAlign(12);
//leg->SetFillColor(kWhite);
//leg->SetLineColor(kWhite);
//leg->SetShadowColor(kWhite);
//leg->SetHeader("Proton MC");
//leg->Draw();

}//<---End Loop()
