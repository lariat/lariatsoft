{
// #######################
// ### Load Data Plots ###
// #######################
TFile *f1 = new TFile("./PionXSection_histos_05_12_00.root");


// ###################################
// ### Load Pion Monte Carlo Plots ###
// ###################################
TFile *f2 = new TFile("./PionMCXSection_weighted_histos.root");



//--------------------------------------------------------------------------------------------------------------
//						dE/dX Plots
//--------------------------------------------------------------------------------------------------------------

// ### Getting the data dE/dX plot ###
TH1F *hDatadEdX = (TH1F*)f1->Get("hdataPiondEdX");

// ### Getting the MC dE/dX plot ###
TH1F *hMCdEdX = (TH1F*)f2->Get("hdataPiondEdX");

// ### Normalizing MC to Data ###
double MCIntegral = hMCdEdX->Integral();
double DataIntegral = hDatadEdX->Integral();

double scaleMC = DataIntegral/MCIntegral;

// ### Scaling MC ###
hMCdEdX->Sumw2();
hMCdEdX->Scale(scaleMC);



// ########################
// ### Making a TCanvas ###
// ########################
TCanvas *c01 = new TCanvas("c01", "dE/dX");
c01->SetTicks();
c01->SetFillColor(kWhite);

hMCdEdX->SetLineColor(kBlue);
hMCdEdX->SetLineStyle(0);
hMCdEdX->SetLineWidth(3);

hDatadEdX->SetLineColor(kBlack);
hDatadEdX->SetLineStyle(0);
hDatadEdX->SetLineWidth(3);
hDatadEdX->SetMarkerStyle(8);
hDatadEdX->SetMarkerSize(0.9);


// ### Drawing the histograms ###

hDatadEdX->Draw("E1");
hMCdEdX->Draw("histosame");

// ### Fitting the MC ###
//hMCdEdX->Fit("landau","","",0, 50);

// ### Fitting the Data ###
hDatadEdX->Fit("landau","","",1.0, 10);

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

}
