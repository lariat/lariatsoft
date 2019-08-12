#define BeamComp64GeVNeg100ARound2_cxx
#include "BeamComp64GeVNeg100ARound2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>
double Polarity=-1; //Sets the pdgs to use for the plots -1 for negative, +1 for positive.

//Setting WCQuality and WCTrack Residual cut values
double MidCut=21; 
double MidCutX=15; 
double MidCutY=15;
double MidCutXMean=3;
double MidCutYMean=0;
double MinResCut=12;

//If you want to smear the momentum or TOF by a percentage, define it here.
double smearfactorp=0.00; //Sigma for the smearing gaussian for the momentum 
double smearfactorlowTOF=0.00; //Sigma for the smearing gaussian for the momentum 
double smearfactorhighTOF=.00; //Sigma for the smearing gaussian for the momentum

//Allows separation of particles by TOF. Was used for TOF smearing in certain ranges, but doesnt do anything right now.
double TOFlowcut=0E-9;
double TOFlowmidcut=24E-9;
double TOFhighmidcut=29E-9;
double TOFhighcut=100E-9;


//Used for "LArIAT" Stats Box.

//float mar_l=0.15;
//float mar_r=0.05;
//float mar_b=0.15;
//float mar_t=0.05;
float mar_l=0.15;
float mar_r=0.05;
float mar_b=0.15;
float mar_t=0.05;
//float hd_x1=mar_l+0.02;
float hd_x1=.58;

//float hd_y2=1.-mar_t-0.02;
float hd_y2=.92;
float axisTitleSize=0.045;
float axisLabelSize=0.045;
float textSize=0.025;

//Picky hist Declare
TH1F *pip=new TH1F("Picky Pion","Picky Pion",100,0,2000);
TH1F *mup=new TH1F("Picky Muon","Picky Muon",100,0,2000);
TH1F *ep=new TH1F("Picky Electron","Picky Electron",100,0,2000);
TH1F *kp=new TH1F("Picky Kaon","Picky Kaon",100,0,2000);
TH1F *pp=new TH1F("Picky Proton","Picky Proton",100,0,2000);
THStack *hp=new THStack("Picky Beam Comp 64GeV Neg 100A","Picky Beam Comp 64GeV Neg 100A");
TH1F *momp=new TH1F("Picky Sim Momentum","Picky Sim Momentum",100,0,2000);
TH1F *dp=new TH1F("Picky Data Momentum","Picky Data Momentum",100,0,2000);

//HY Four Point hist Declare
TH1F *pi4=new TH1F("Four Point HY Pion","Four Point HY Pion",100,0,2000);
TH1F *mu4=new TH1F("Four Point HY Muon","Four Point HY Muon",100,0,2000);
TH1F *e4=new TH1F("Four Point HY Electron","Four Point HY Electron",100,0,2000);
TH1F *k4=new TH1F("Four Point HY Kaon","Four Point HY Kaon",100,0,2000);
TH1F *p4=new TH1F("Four Point HY Proton","Four Point HY Proton",100,0,2000);
THStack *h4=new THStack("","");
TH1F *mom4=new TH1F("Four Point HY Sim Momentum","Four Point HY Sim Momentum",100,0,2000);
TH1F *d4=new TH1F("Four Point HY Data Momentum","Four Point HY Data Momentum",100,0,2000);
TH1F *mom4Low=new TH1F("Low Mass Four Point HY Sim Momentum","Low Mass Four Point HY Sim Momentum",100,0,2000);
//HY S2 hist Declare
TH1F *pi2=new TH1F("S2 HY Pion","S2 HY Pion",100,0,2000);
TH1F *mu2=new TH1F("S2 HY Muon","S2 HY Muon",100,0,2000);
TH1F *e2=new TH1F("S2 HY Electron","S2 HY Electron",100,0,2000);
TH1F *k2=new TH1F("S2 HY Kaon","S2 HY Kaon",100,0,2000);
TH1F *p2=new TH1F("S2 HY Proton","S2 HY Proton",100,0,2000);
THStack *h2=new THStack("S2 HY Beam Comp 64GeV Neg 100A","S2 HY Beam Comp 64GeV Neg 100A");
TH1F *mom2=new TH1F("S2 HY Sim Momentum","S2 HY Sim Momentum",100,0,2000);
TH1F *d2=new TH1F("S2 HY Data Momentum","S2 HY Data Momentum",100,0,2000);

//HY S3 hist Declare
TH1F *pi3=new TH1F("S3 HY Pion","S3 HY Pion",100,0,2000);
TH1F *mu3=new TH1F("S3 HY Muon","S3 HY Muon",100,0,2000);
TH1F *e3=new TH1F("S3 HY Electron","S3 HY Electron",100,0,2000);
TH1F *k3=new TH1F("S3 HY Kaon","S3 HY Kaon",100,0,2000);
TH1F *p3=new TH1F("S3 HY Proton","S3 HY Proton",100,0,2000);
THStack *h3=new THStack("S3 HY Beam Comp 64GeV Neg 100A","S3 HY Beam Comp 64GeV Neg 100A");
TH1F *mom3=new TH1F("S3 HY Sim Momentum","S3 HY Sim Momentum",100,0,2000);
TH1F *d3=new TH1F("S3 HY Data Momentum","S3 HY Data Momentum",100,0,2000);

//Declare the error plots
TH2F *ratp= new TH2F("Fractional Error of Truth and Sim: Picky Tracks","Fractional Error of Truth and Sim: Picky Tracks",100,0,2000,400,-1,1);
TH1F *diffp=new TH1F("Raw Error of Truth and Sim: Picky Tracks","Raw Error of Truth and Sim: Picky Tracks",400,-.02,.02);

TH2F *rat4= new TH2F("Fractional Error of Truth and Sim: Four Point HY Tracks","Fractional Error of Truth and Sim: Four Point HY Tracks",100,0,2000,400,-1,1);
TH1F *dif4=new TH1F("Raw Error of Truth and Sim: Four Point HY Tracks","Raw Error of Truth and Sim: Four Point HY Tracks",400,-.02,.02);

TH2F *rat2= new TH2F("Fractional Error of Truth and Sim: S2 Tracks","Fractional Error of Truth and Sim: S2 Tracks",100,0,2000,400,-1,1);
TH1F *diff2=new TH1F("Raw Error of Truth and Sim: S2 Tracks","Raw Error of Truth and Sim: S2 Tracks",400,-.02,.02);

TH2F *rat3= new TH2F("Fractional Error of Truth and Sim: S3 Tracks","Fractional Error of Truth and Sim: S3 Tracks",100,0,2000,400,-1,1);
TH1F *diff4=new TH1F("Raw Error of Truth and Sim: S3 Tracks","Raw Error of Truth and Sim: S3 Tracks",400,-.02,.02);


//Plots Only in sim: "pure" and "frankentrack"
TH1F *mompurep=new TH1F("Pure Picky Sim Momentum","Pure Picky Sim Momentum",100,0,2000);
TH1F *momfrankp=new TH1F("Frankentrack Picky Sim Momentum","Frankentrack Picky Sim Momentum",100,0,2000);

TH1F *mompure4=new TH1F("Pure Four Point HY Sim Momentum","Pure Four Point HY Sim Momentum",100,0,2000);
TH1F *momfrank4=new TH1F("Frankentrack Four Point HY Sim Momentum","Frankentrack Four Point HY Sim Momentum",100,0,2000);

TH1F *mompure2=new TH1F("Pure S2 Sim Momentum","Pure S2 Sim Momentum",100,0,2000);
TH1F *momfrank2=new TH1F("Frankentrack S2 Sim Momentum","Frankentrack S2 Sim Momentum",100,0,2000);

TH1F *mompure3=new TH1F("Pure S3 Sim Momentum","Pure S3 Sim Momentum",100,0,2000);
TH1F *momfrank3=new TH1F("Frankentrack S3 Sim Momentum","Frankentrack S3 Sim Momentum",100,0,2000);


//HY Mass Plots
TH1F *pi4m=new TH1F("Four Point HY Pion Mass","Four Point HY Pion Mass",400,-2000,2000);
TH1F *mu4m=new TH1F("Four Point HY Muon Mass","Four Point HY Muon Mass",400,-2000,2000);
TH1F *e4m=new TH1F("Four Point HY Electron Mass","Four Point HY Electron Mass",400,-2000,2000);
TH1F *k4m=new TH1F("Four Point HY Kaon Mass","Four Point HY Kaon Mass",400,-2000,2000);
TH1F *p4m=new TH1F("Four Point HY Proton Mass","Four Point HY Proton Mass",400,-2000,2000);
TH1F *d4m=new TH1F("Data Mass","Data Mass",400,-2000,2000);
//THStack *h4m=new THStack("Four Point HY Beam Comp 64GeV Pos 100A Mass","Four Point HY Beam Comp 64GeV Pos 100A Mass");
THStack *h4m=new THStack("","");
//HY Mass Plots
TH1F *pi4t=new TH1F("Four Point HY Pion TOF","Four Point HY Pion TOF",400,0,100);
TH1F *mu4t=new TH1F("Four Point HY Muon TOF","Four Point HY Muon TOF",400,0,100);
TH1F *e4t=new TH1F("Four Point HY Electron TOF","Four Point HY Electron TOF",400,0,100);
TH1F *k4t=new TH1F("Four Point HY Kaon TOF","Four Point HY Kaon TOF",400,0,100);
TH1F *p4t=new TH1F("Four Point HY Proton TOF","Four Point HY Proton TOF",400,0,100);
TH1F *d4t=new TH1F("Data TOF","Data TOF",400,0,100);
//THStack *h4t=new THStack("Four Point HY Beam Comp 64GeV Pos 100A TOF","Four Point HY Beam Comp 64GeV Pos 100A TOF");
THStack *h4t=new THStack("","");
//Plots people want

TH1F *binfracpi= new TH1F("Pion Bin by Bin Content","Pion Bin By Bin Content",100,0,2000);
TH1F *binfracmu= new TH1F("Muon Bin by Bin Content","Muon Bin By Bin Content",100,0,2000);
TH1F *binfrace= new TH1F("Electron Bin by Bin Content","Electron Bin By Bin Content",100,0,2000);
TH1F *binfrack= new TH1F("Kaon Bin by Bin Content","Kaon Bin By Bin Content",100,0,2000);
TH1F *binfracp= new TH1F("Proton Bin by Bin Content","Proton Bin By Bin Content",100,0,2000);
void BeamComp64GeVNeg100ARound2::Loop()
{
//Get the data plots
//This could be cleaned up once we have a "production" of beamline reco and a standardized file system.
TFile *fp=new TFile("../FinalBeamComp/WCTrackHistosNeg100APicky.root"); //Change this your data picky tracks after WCFilter
TH1F *dp=(TH1F*)fp->Get("FourPointMom");
TFile *fh=new TFile("../FinalBeamComp/JHoWCTrackHistosNeg100.root"); //Change this to your data 3 pt tracks after WCFilt 
TFile *fhy=new TFile("../FinalBeamComp/Data100Plots/Neg100A.root");  //File With the TOF and HY tracks
TH1F *d4=(TH1F*)fhy->Get("wcquality/MomAfterQuality");
TH1F *d2=(TH1F*)fh->Get("Skip2Mom");
TH1F *d3=(TH1F*)fh->Get("Skip3Mom");
TH1F *d4m=(TH1F*)fh->Get("Mass");
TH1F *d4t=(TH1F*)fhy->Get("tof/TOF");


TRandom *r3=new TRandom3(0);
//Now lets get the Sim plots
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   //use RNG to get a smear factor for each track....assuming you're smearing. 
   double arrayofsmearsp[nentries];
   double arrayofsmearslowTOF[nentries];
   double arrayofsmearshighTOF[nentries];
   for (int i=0; i<nentries; i++)
   {
     Double_t factorp=r3->Gaus(0,smearfactorp);
     arrayofsmearsp[i]=1+factorp;
     
     Double_t factorlowTOF=r3->Gaus(0,smearfactorlowTOF);
     arrayofsmearslowTOF[i]=1+factorlowTOF;
     
     Double_t factorhighTOF=r3->Gaus(0,smearfactorhighTOF);
     arrayofsmearshighTOF[i]=1+factorhighTOF;
   }    
   
   //Time to loop over the tracks
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //Skip events that fail WCQuality Cuts, or doesnt have a TOF reco'd
      if(!WCFilt || TOF<-1000 ||fabs(MidDiffX-MidCutXMean)>MidCutX ||fabs(MidDiffY-MidCutYMean)>MidCutY ||MinRes>MinResCut){continue;}
      
      // Fill Sim-only plots, which don't get separated by particle species
      if(PickyTrack && PureTrack && FourptTrack){mompurep->Fill(Momentum);}
      if(PickyTrack && Frankentrack && FourptTrack){momfrankp->Fill(Momentum);}

      if((PickyTrack||HYTrack) && FourptTrack && PureTrack){mompure4->Fill(Momentum);}
      if((PickyTrack||HYTrack) && FourptTrack && Frankentrack){momfrank4->Fill(Momentum);}

      if((PickyTrack||HYTrack) && S2Track && PureTrack){mompure2->Fill(Momentum);}
      if((PickyTrack||HYTrack) && S2Track && Frankentrack){momfrank2->Fill(Momentum);}
      
      if((PickyTrack||HYTrack) && S3Track && PureTrack){mompure3->Fill(Momentum);}
      if((PickyTrack||HYTrack) && S3Track && Frankentrack){momfrank3->Fill(Momentum);}
      
      //Easier to have these plots regardless of species. This will also exist in the form of a THStack
      if(PickyTrack && FourptTrack){momp->Fill(Momentum);}
      if((PickyTrack||HYTrack) && FourptTrack){mom4->Fill(Momentum*arrayofsmearsp[jentry]);}
      if((PickyTrack||HYTrack) && S2Track){mom2->Fill(Momentum);}
      if((PickyTrack||HYTrack) && S3Track){mom3->Fill(Momentum);}
      
      //Fill Sim plots that are separated by particle species. 
      //This is a lot of copy paste, so here's what happens to about line 436.
      //For each of the 5 particle species, fill a picky track mom hist. For HY 4pt tracks, fill a mom hist, using a smeared momentum, if you've smeared. 
      //For the various TOF ranges, fill the TOF plot for that particle, using the relevant smearing. With the smeared momentum and TOF, calculate the mass.
      //Fill histogram for the 3pt tracks by particle, if it's a 3pt event.   
      double ThisPDG=FindPDG(trackpdgs,MotherIDs,trackIDs);             
      if(ThisPDG==Polarity*211)
      {
        if(PickyTrack && FourptTrack){pip->Fill(Momentum);}
	if((PickyTrack || HYTrack) && FourptTrack)
	{
	  //pi4t->Fill(TOF*1E9);
	  pi4->Fill(Momentum*arrayofsmearsp[jentry]);
	  mom4Low->Fill(Momentum*arrayofsmearsp[jentry]);
	  if(TOF>TOFlowcut && TOF<TOFlowmidcut)
	  {
	  pi4t->Fill(TOF*arrayofsmearslowTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  pi4m->Fill(smass);
	  }
	  
	  else if(TOF>TOFhighmidcut && TOF<TOFhighcut)
	  {
	  pi4t->Fill(TOF*arrayofsmearshighTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	 // std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  pi4m->Fill(smass);
	  }
	  
	  else
	  {
	  pi4t->Fill(TOF*1E9);
	  double sunderroot=TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1)<0){smass=-smass;}
	 // std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  pi4m->Fill(smass);
	  }
	}
	if((PickyTrack || HYTrack) && S2Track){pi2->Fill(Momentum);}
	if((PickyTrack || HYTrack) && S3Track){pi3->Fill(Momentum);}
	
      }
      if(ThisPDG==Polarity*-13)
      {
        if(PickyTrack && FourptTrack){mup->Fill(Momentum);}
	if((PickyTrack || HYTrack) && FourptTrack)
	{
	  //mu4t->Fill(TOF*1E9); 
	  mu4->Fill(Momentum*arrayofsmearsp[jentry]);
	  mom4Low->Fill(Momentum*arrayofsmearsp[jentry]);
	  if(TOF>TOFlowcut && TOF<TOFlowmidcut)
	  {
	  mu4t->Fill(TOF*arrayofsmearslowTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  mu4m->Fill(smass);
	  }
	  
	  else if(TOF>TOFhighmidcut && TOF<TOFhighcut)
	  {
	  mu4t->Fill(TOF*arrayofsmearshighTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  mu4m->Fill(smass);
	  }
	  
	  else
	  {
	  mu4t->Fill(TOF*1E9);
	  double sunderroot=TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  mu4m->Fill(smass);
	  }
	}
	if((PickyTrack || HYTrack) && S2Track){mu2->Fill(Momentum);}
	if((PickyTrack || HYTrack) && S3Track){mu3->Fill(Momentum);}        
      }
      if(ThisPDG==Polarity*-11)
      {
        if(PickyTrack && FourptTrack){ep->Fill(Momentum);}
	if((PickyTrack || HYTrack) && FourptTrack)
	{
	  //e4t->Fill(TOF*1E9);
	  e4->Fill(Momentum*arrayofsmearsp[jentry]);
	   mom4Low->Fill(Momentum*arrayofsmearsp[jentry]);
	  if(TOF>TOFlowcut && TOF<TOFlowmidcut)
	  {
	  e4t->Fill(TOF*arrayofsmearslowTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  e4m->Fill(smass);
	  }
	  
	  else if(TOF>TOFhighmidcut && TOF<TOFhighcut)
	  {
	  e4t->Fill(TOF*arrayofsmearshighTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  e4m->Fill(smass);
	  }
	  
	  else
	  {
	  e4t->Fill(TOF*1E9);
	  double sunderroot=TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  e4m->Fill(smass);
	  }	
	}
	if((PickyTrack || HYTrack) && S2Track){e2->Fill(Momentum);}
	if((PickyTrack || HYTrack) && S3Track){e3->Fill(Momentum);}        
      }
      if(ThisPDG==Polarity*321)
      {
        if(PickyTrack && FourptTrack){kp->Fill(Momentum);}
	if((PickyTrack || HYTrack) && FourptTrack)
	{
	  //k4t->Fill(TOF*1E9);
	  k4->Fill(Momentum*arrayofsmearsp[jentry]);
	  if(TOF>TOFlowcut && TOF<TOFlowmidcut)
	  {
	  k4t->Fill(TOF*arrayofsmearslowTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  k4m->Fill(smass);
	  }
	  
	  else if(TOF>TOFhighmidcut && TOF<TOFhighcut)
	  {
	  k4t->Fill(TOF*arrayofsmearshighTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  k4m->Fill(smass);
	  }
	  
	  else
	  {
	  k4t->Fill(TOF*1E9);
	  double sunderroot=TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  k4m->Fill(smass);
	  }  
	}
	if((PickyTrack || HYTrack) && S2Track){k2->Fill(Momentum);}
	if((PickyTrack || HYTrack) && S3Track){k3->Fill(Momentum);}        
      }
      if(ThisPDG==Polarity*2212)
      {
        if(PickyTrack && FourptTrack){pp->Fill(Momentum);}
	if((PickyTrack || HYTrack) && FourptTrack)
	{
	  //p4t->Fill(TOF*1E9);
	  p4->Fill(Momentum*arrayofsmearsp[jentry]);
	  if(TOF>TOFlowcut && TOF<TOFlowmidcut)
	  {
	  p4t->Fill(TOF*arrayofsmearslowTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearslowTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearslowTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  p4m->Fill(smass);
	  }
	  
	  else if(TOF>TOFhighmidcut && TOF<TOFhighcut)
	  {
	  p4t->Fill(TOF*arrayofsmearshighTOF[jentry]*1E9);
	  double sunderroot=TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*arrayofsmearshighTOF[jentry]*0.299792458E9*0.299792458E9*TOF*arrayofsmearshighTOF[jentry]/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  p4m->Fill(smass);
	  }
	  
	  else
	  {
	  p4t->Fill(TOF*1E9);
	  double sunderroot=TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1;
	  if (sunderroot<0){sunderroot=-sunderroot;}
	  double smass = Momentum*arrayofsmearsp[jentry]*pow(sunderroot,0.5);
	  if((TOF*0.299792458E9*0.299792458E9*TOF/(6.652*6.652) - 1)<0){smass=-smass;}
	  //std::cout<<Momentum<<" "<<TOF<<" "<<smass<<std::endl;
	  p4m->Fill(smass);
	  }	  
	}//if track is 4pt
	if((PickyTrack || HYTrack) && S2Track){p2->Fill(Momentum);}
	if((PickyTrack || HYTrack) && S3Track){p3->Fill(Momentum);}        
      }//if track is proton                          
   } //jentry
   
   
//Now format some fills colors and lines.
   //Picky plots with bad color scheme.    
   pip->SetLineColor(9);
   mup->SetLineColor(41);
   ep->SetLineColor(15);
   kp->SetLineColor(7);
   pp->SetLineColor(2);
   pip->SetFillColor(9);
   mup->SetFillColor(41);
   ep->SetFillColor(15);
   kp->SetFillColor(7);
   pp->SetFillColor(2);

   //Four Point HY with "cool" color scheme

   char picolhex[10]="#b7dfcb";
   char pcolhex[10]="#5abad1";
   char mucolhex[10]="#3984b6";
   char kcolhex[10]="#161f63";
   char ecolhex[10]="#726499";

   
   TColor *TColGreg=new TColor();
   Int_t col1=TColGreg->GetColor(picolhex);
   Int_t col2=TColGreg->GetColor(pcolhex);
   Int_t col3=TColGreg->GetColor(mucolhex);
   Int_t col4=TColGreg->GetColor(ecolhex);
   Int_t col5=TColGreg->GetColor(kcolhex);
   pi4->SetLineColor(col1);
   mu4->SetLineColor(col3);
   e4->SetLineColor(col4);
   k4->SetLineColor(col5);
   p4->SetLineColor(col2);
   pi4->SetFillColor(col1);
   mu4->SetFillColor(col3);
   e4->SetFillColor(col4);
   k4->SetFillColor(col5);
   p4->SetFillColor(col2);
   pi4->SetMarkerColor(col1);
   mu4->SetMarkerColor(col3);
   e4->SetMarkerColor(col4);
   k4->SetMarkerColor(col5);
   p4->SetMarkerColor(col2); 
     
   //Skip 2
   pi2->SetLineColor(9);
   mu2->SetLineColor(41);
   e2->SetLineColor(15);
   k2->SetLineColor(7);
   p2->SetLineColor(2);
   pi2->SetFillColor(9);
   mu2->SetFillColor(41);
   e2->SetFillColor(15);
   k2->SetFillColor(7);
   p2->SetFillColor(2);      

   //Skip 3
   pi3->SetLineColor(9);
   mu3->SetLineColor(41);
   e3->SetLineColor(15);
   k3->SetLineColor(7);
   p3->SetLineColor(2);
   pi3->SetFillColor(9);
   mu3->SetFillColor(41);
   e3->SetFillColor(15);
   k3->SetFillColor(7);
   p3->SetFillColor(2);

   //Set Data markers to be more visible
   
   dp->SetMarkerStyle(20);
   d4->SetMarkerStyle(20); 
   d2->SetMarkerStyle(20); 
   d3->SetMarkerStyle(20);
   d4t->SetMarkerStyle(20);  
   d4m->SetMarkerStyle(20);
   //Get normalization constants   
   double dpint=dp->Integral();
   double d4int=d4->Integral();
   double d2int=d2->Integral();
   double d3int=d3->Integral();
   
   double mompint=momp->Integral();
   double mom4int=mom4->Integral();
   double mom4intlow=mom4Low->Integral();
   double mom2int=mom2->Integral();
   double mom3int=mom3->Integral();
   //Need these numbers for appropriate error propogation later on
   double alphap=mompint/dpint;
   double alpha4=mom4int/d4int;
   double alpha2=mom2int/d2int;
   double alpha3=mom3int/d3int;
   //Get total number of tracks in each sample
   double totalp=pip->Integral()+mup->Integral()+ep->Integral()+kp->Integral()+pp->Integral();
   double total4=pi4->Integral()+mu4->Integral()+e4->Integral()+k4->Integral()+p4->Integral();
   double total2=pi2->Integral()+mu2->Integral()+e2->Integral()+k2->Integral()+p2->Integral();
   double total3=pi3->Integral()+mu3->Integral()+e3->Integral()+k3->Integral()+p3->Integral();  
   double total4low=pi4->Integral()+mu4->Integral()+e4->Integral();
      
   //Normalize particle species sim plots
   pip->Scale(1/totalp);
   mup->Scale(1/totalp);
   ep->Scale(1/totalp);
   kp->Scale(1/totalp);
   pp->Scale(1/totalp);
   
   pi4->Scale(1/total4);
   mu4->Scale(1/total4);
   e4->Scale(1/total4);
   k4->Scale(1/total4);
   p4->Scale(1/total4);
   
   pi2->Scale(1/total2);
   mu2->Scale(1/total2);
   e2->Scale(1/total2);
   k2->Scale(1/total2);
   p2->Scale(1/total2);
   
   pi3->Scale(1/total3);
   mu3->Scale(1/total3);
   e3->Scale(1/total3);
   k3->Scale(1/total3);
   p3->Scale(1/total3);
   
   //Fill the THStacks
   hp->Add(pp,"hist");
   hp->Add(kp,"hist");
   hp->Add(ep,"hist");
   hp->Add(mup,"hist");
   hp->Add(pip,"hist");
   
   //h4->Add(p4,"hist");
   h4->Add(k4,"hist");
   h4->Add(e4,"hist");
   h4->Add(mu4,"hist");
   h4->Add(p4,"hist");
   h4->Add(pi4,"hist");
   
   //h2->Add(p2,"hist");
   //h2->Add(k2,"hist");
   //h2->Add(e2,"hist");
   //h2->Add(mu2,"hist");
   h2->Add(pi2,"hist");
   
   h3->Add(p3,"hist");
   h3->Add(k3,"hist");
   h3->Add(e3,"hist");
   h3->Add(mu3,"hist");
   h3->Add(pi3,"hist");
   
   //Optional. Print the composition for the 4 sim/data compared samples
   std::cout<<"====================== Picky Tracks ======================"<<std::endl;
   std::cout<<"Total Tracks: "<<totalp<<std::endl;
   std::cout<<"Pion: "<<100*pip->Integral()<<std::endl;
   std::cout<<"Muon: "<<100*mup->Integral()<<std::endl;
   std::cout<<"Electron: "<<100*ep->Integral()<<std::endl;
   std::cout<<"Kaon: "<<100*kp->Integral()<<std::endl;
   std::cout<<"Proton: "<<100*pp->Integral()<<std::endl;
   std::cout<<"=========================================================="<<std::endl;  
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<"====================== HY Four Point Tracks ======================"<<std::endl;
   std::cout<<"Total Tracks: "<<total4<<std::endl;
   std::cout<<"Pion: "<<100*pi4->Integral()<<std::endl;
   std::cout<<"Muon: "<<100*mu4->Integral()<<std::endl;
   std::cout<<"Electron: "<<100*e4->Integral()<<std::endl;
   std::cout<<"Kaon: "<<100*k4->Integral()<<std::endl;
   std::cout<<"Proton: "<<100*p4->Integral()<<std::endl;    
   std::cout<<"=================================================================="<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<"====================== HY Skip 2 Tracks ======================"<<std::endl;
   std::cout<<"Total Tracks: "<<total2<<std::endl;
   std::cout<<"Pion: "<<100*pi2->Integral()<<std::endl;
   std::cout<<"Muon: "<<100*mu2->Integral()<<std::endl;
   std::cout<<"Electron: "<<100*e2->Integral()<<std::endl;
   std::cout<<"Kaon: "<<100*k2->Integral()<<std::endl;
   std::cout<<"Proton: "<<100*p2->Integral()<<std::endl;    
   std::cout<<"=============================================================="<<std::endl; 
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<""<<std::endl;
   std::cout<<"====================== HY Skip 3 Tracks ======================"<<std::endl;
   std::cout<<"Total Tracks: "<<total3<<std::endl;
   std::cout<<"Pion: "<<100*pi3->Integral()<<std::endl;
   std::cout<<"Muon: "<<100*mu3->Integral()<<std::endl;
   std::cout<<"Electron: "<<100*e3->Integral()<<std::endl;
   std::cout<<"Kaon: "<<100*k3->Integral()<<std::endl;
   std::cout<<"Proton: "<<100*p3->Integral()<<std::endl;    
   std::cout<<"=============================================================="<<std::endl; 


   //Calculate the bin by bin errors
   
   //Picky Tracks
   int binsp=dp->GetNbinsX();
   double errorpy[binsp];
   double errorpx[binsp];
   double xp[binsp];
   double yp[binsp];
   for (int i=0; i<binsp; i++)
   {
     double dbin=dp->GetBinContent(i);
     double tbin=momp->GetBinContent(i);
     if(dbin>0 && tbin>0)
     {
     xp[i]=20*i;
     yp[i]=((tbin/mompint)-(dbin/dpint))/((tbin/mompint));
     errorpy[i]=alphap/tbin*pow(dbin*(1+dbin/tbin),.5);
     errorpx[i]=10;
     }
     else
     {
       xp[i]=0;
       yp[i]=0;
       errorpx[i]=0;
       errorpy[i]=0;
     }
   }
   
   //HY 4 Point
   int bins4=d4->GetNbinsX();
   double error4y[binsp];
   double error4x[binsp];
   double x4[binsp];
   double y4[binsp];
   for (int i=0; i<bins4; i++)
   {
     double dbin=d4->GetBinContent(i);
     double tbin=mom4->GetBinContent(i);
     if(dbin>0 && tbin>0)
     {
     x4[i]=20*i;
     y4[i]=((tbin/mom4int)-(dbin/d4int))/((tbin/mom4int));
     error4y[i]=alpha4/tbin*pow(dbin*(1+dbin/tbin),.5);
     error4x[i]=10;
     }
     else
     {
       x4[i]=0;
       y4[i]=0;
       error4x[i]=0;
       error4y[i]=0;
     }
   }
   
   //HY S2
   int bins2=d2->GetNbinsX();
   double error2y[binsp];
   double error2x[binsp];
   double x2[binsp];
   double y2[binsp];
   for (int i=0; i<bins2; i++)
   {
     double dbin=d2->GetBinContent(i);
     double tbin=mom2->GetBinContent(i);
     if(dbin>0 && tbin>0)
     {
     x2[i]=20*i;
     y2[i]=((tbin/mom2int)-(dbin/d2int))/((tbin/mom2int));
     error2y[i]=alpha2/tbin*pow(dbin*(1+dbin/tbin),.5);
     error2x[i]=10;
     }
     else
     {
       x2[i]=0;
       y2[i]=0;
       error2x[i]=0;
       error2y[i]=0;
     }
   }

   //HY S3 Point
   int bins3=d3->GetNbinsX();
   double error3y[binsp];
   double error3x[binsp];
   double x3[binsp];
   double y3[binsp];
   for (int i=0; i<bins3; i++)
   {
     double dbin=d3->GetBinContent(i);
     double tbin=mom3->GetBinContent(i);
     if(dbin>0 && tbin>0)
     {
     x3[i]=20*i;
     y3[i]=((tbin/mom3int)-(dbin/d3int))/((tbin/mom3int));
     error3y[i]=alpha3/tbin*pow(dbin*(1+dbin/tbin),.5);
     error3x[i]=10;
     }
     else
     {
       x3[i]=0;
       y3[i]=0;
       error3x[i]=0;
       error3y[i]=0;
     }
   }           
   
   //Scale data and species-independent plots. We have to do this after we calculate the error.
   dp->Scale(1/dpint);
   momp->Scale(1/mompint);

   d4->Scale(1/d4int);
   mom4->Scale(1/mom4int);
   
   d2->Scale(1/d2int);
   mom2->Scale(1/mom2int);
   
   d3->Scale(1/d3int);
   mom3->Scale(1/mom3int);
   
   mom4Low->Scale(1/mom4intlow);
   //Fill bin by bin content per particle
   
   TH1F* binfracpi=(TH1F*)pi4->Clone("binfracpi");
   TH1F* binfracmu=(TH1F*)mu4->Clone("binfracmu");
   TH1F* binfrace=(TH1F*) e4->Clone("binfrace");
   TH1F* binfrack=(TH1F*) k4->Clone("binfrack");
   TH1F* binfracp=(TH1F*) p4->Clone("binfracp");
   std::cout<<"mom4low: "<<mom4intlow<<std::endl;
   binfracpi->Scale(total4/total4low);
   binfracmu->Scale(total4/total4low);
   binfrace->Scale(total4/total4low);
   binfracpi->Divide(mom4Low);
   binfracmu->Divide(mom4Low);
   binfrace->Divide(mom4Low);
   //binfracpi->Divide(mom4);
   //binfracmu->Divide(mom4);
   //binfrace->Divide(mom4);
   binfrack->Divide(mom4);
   binfracp->Divide(mom4);
   
   
 //Make the Picky Track Momentum Comparion Plot       
   auto cp= new TCanvas("Picky Sim/Data Comparison","Picky Sim/Data Comparison");
   auto legp=new TLegend(0.5,0.4,0.9,0.9);
   legp->SetHeader("Picky Particle Species at WC4, Post WCFilter");
   legp->AddEntry(pip,TString::Format("#pi: %.1f %%",100*pip->Integral()));
   legp->AddEntry(mup,TString::Format("#mu: %.1f %%",100*mup->Integral()));
   legp->AddEntry(ep,TString::Format("e: %.1f %%",100*ep->Integral()));
   legp->AddEntry(kp,TString::Format("k: %.1f %%",100*kp->Integral()));
   legp->AddEntry(pp,TString::Format("p: %.1f %%",100*pp->Integral()));
   legp->AddEntry(dp,"+100A Picky Track, WCFiltered Data");
   hp->Draw();
   dp->Draw("same");
   legp->Draw();
   
   
   hp->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   hp->GetYaxis()->SetTitle("Entries per 20 MeV/c");
   cp->SetGridx();
   cp->SetGridy();
   
   //Make the Hy 4 Point Momentum Comparion Plot       
   auto c4= new TCanvas("HY 4 Point Sim to Data Comparison","HY 4 Point Sim to Data Comparison",800,800);
   auto leg4=new TLegend(0.6,0.5,0.94,0.8);
   leg4->SetHeader("64 GeV +100 Amps");
   leg4->AddEntry(pi4,TString::Format("#pi: %.1f%%",100*pi4->Integral()));
   leg4->AddEntry(p4,TString::Format("p: %.1f%%",100*p4->Integral()));
   leg4->AddEntry(mu4,TString::Format("#mu: %.1f%%",100*mu4->Integral()));
   leg4->AddEntry(e4,TString::Format("e: %.1f%%",100*e4->Integral()));
   leg4->AddEntry(k4,TString::Format("K: %.1f%%",100*k4->Integral()));
   
   leg4->AddEntry(d4,"Data");
   h4->Draw();
   d4->Draw("same");
   leg4->Draw();
   leg4->SetTextSize(textSize);
   h4->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   h4->GetXaxis()->SetRangeUser(200,1600);
   h4->GetYaxis()->SetTitle("Normalized Entries per 20 MeV/c");
   h4->GetXaxis()->SetTitleOffset(1.1);
   h4->GetYaxis()->SetTitleOffset(1.6);
   c4->SetGridx();
   c4->SetGridy();
   gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);
   TPaveText* header=MakeTextBox(hd_x1,hd_y2,textSize,3,0.3);
   header->AddText("#bf{LArIAT}");
   header->AddText(Form("%.0f Data Events",d4int));
   header->AddText(Form("%.0f Simulated Events",total4));
   header->Draw();

/*     //Make the S2 Point Momentum Comparion Plot       
   auto c2= new TCanvas("HY S2 Sim/Data Comparison","HY S2 Sim/Data Comparison");
   auto leg2=new TLegend(0.5,0.4,0.9,0.9);
   leg2->SetHeader("HY S2 Particle Species at WC4, Post WCFilter");
   leg2->AddEntry(pi2,"#pi: 78.82%");
   leg2->AddEntry(mu2,"#mu: 10.59%");
   leg2->AddEntry(e2,"e: 10.37%");
   leg2->AddEntry(k2,"k: .13%");
   leg2->AddEntry(p2,"p: .09%");
   leg2->AddEntry(d2,"-100A HY S2 Track, WCFiltered Data");
   h2->Draw();
   d2->Draw("same");
   leg2->Draw();
   h2->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   h2->GetYaxis()->SetTitle("Entries per 20 MeV/c");
   c2->SetGridx();
   c2->SetGridy();
    */
   
/*   //Make the S3 Point Momentum Comparion Plot       
   auto c3= new TCanvas("HY S3 Sim/Data Comparison","HY S3 Sim/Data Comparison");
   auto leg3=new TLegend(0.5,0.4,0.9,0.9);
   leg3->SetHeader("HY S3 Particle Species at WC4, Post WCFilter");
   leg3->AddEntry(pi3,"#pi: 82.93%");
   leg3->AddEntry(mu3,"#mu: 11.07%");
   leg3->AddEntry(e3,"e: 5.05%");
   leg3->AddEntry(k3,"k: .78%");
   leg3->AddEntry(p3,"p: .13%");
   leg3->AddEntry(d3,"-100A HY S3 Track, WCFiltered Data");
   h3->Draw();
   d3->Draw("same");
   leg3->Draw();
   h3->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   h3->GetYaxis()->SetTitle("Entries per 20 MeV/c");
   c3->SetGridx();
   c3->SetGridy();    */
   
   //Make a new TCanvas and plot the bin by bin ratio with errors
   
   //PickyTracks
   auto cerrorp= new TCanvas("Picky Tracks: Error Ratio","Picky Tracks: Error Ratio");
   TGraphErrors *grp = new TGraphErrors(binsp,xp,yp,errorpx,errorpy);
   grp->SetTitle("Picky Tracks: Bin by Bin Fractional Error Between Sim and Data");
   grp->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   grp->GetYaxis()->SetTitle("#frac{Truth-Data}{Truth}");
   cerrorp->SetGridy();
   grp->Draw("AP");
   
   //HY 4 Point
   auto cerror4= new TCanvas("HY 4 Point: Error Ratio","HY 4 Point: Error Ratio");
   TGraphErrors *gr4 = new TGraphErrors(bins4,x4,y4,error4x,error4y);
   gr4->SetTitle("Neg 100A HY 4 Point: Bin by Bin Fractional Error Between Sim and Data");
   gr4->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   gr4->GetYaxis()->SetTitle("#frac{Truth-Data}{Truth}");
   cerror4->SetGridy();
   gr4->Draw("AP");


/*    //S2
   auto cerror2= new TCanvas("HY S2: Error Ratio","HY S2: Error Ratio");
   TGraphErrors *gr2 = new TGraphErrors(bins2,x2,y2,error2x,error2y);
   gr2->SetTitle("HY S2: Bin by Bin Fractional Error Between Sim and Data");
   gr2->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   gr2->GetYaxis()->SetTitle("#frac{Truth-Data}{Truth}");
   cerror2->SetGridy();
   gr2->Draw("AP");

   //S3
   auto cerror3= new TCanvas("HY S3: Error Ratio","HY 4 S3: Error Ratio");
   TGraphErrors *gr3 = new TGraphErrors(bins3,x3,y3,error3x,error3y);
   gr3->SetTitle("HY S3: Bin by Bin Fractional Error Between Sim and Data");
   gr3->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   gr3->GetYaxis()->SetTitle("#frac{Truth-Data}{Truth}");
   cerror3->SetGridy();
   gr3->Draw("AP");     */   
   

//HY Four Point hist Declare
TH1F *pi4c=new TH1F("Corrected Four Point HY Pion","Corrected Four Point HY Pion",100,0,2000);
TH1F *mu4c=new TH1F("Corrected Four Point HY Muon","Corrected Four Point HY Muon",100,0,2000);
TH1F *e4c=new TH1F("Corrected Four Point HY Electron","Corrected Four Point HY Electron",100,0,2000);
TH1F *k4c=new TH1F("Corrected Four Point HY Kaon","Corrected Four Point HY Kaon",100,0,2000);
TH1F *p4c=new TH1F("Corrected Four Point HY Proton","Corrected Four Point HY Proton",100,0,2000);
THStack *h4c=new THStack("Corrected Four Point HY Beam Comp 64GeV Neg 100A","Corrected Four Point HY Beam Comp 64GeV Neg 100A");


//Once we have the original plots and errors, rerun the script, using the weight to fill the species histograms instead. These histograms will cause sim to match data perfectly.
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(!WCFilt || !(PickyTrack  || HYTrack) || !FourptTrack){continue;}
      Int_t flooredMom=TMath::Floor(Momentum*arrayofsmearsp[jentry]);
      Int_t remainder=flooredMom%20;
      Int_t weightbin=(flooredMom-remainder)/20;
      double ThisPDG=FindPDG(trackpdgs,MotherIDs,trackIDs);
      if(ThisPDG==Polarity*211)
      {
	pi4c->Fill(Momentum*arrayofsmearsp[jentry],1-y4[weightbin+1]);	
      }
      if(ThisPDG==Polarity*-13)
      {
        mu4c->Fill(Momentum*arrayofsmearsp[jentry],1-y4[weightbin+1]);      
      }
      if(ThisPDG==Polarity*-11)
      {
	e4c->Fill(Momentum*arrayofsmearsp[jentry],1-y4[weightbin+1]);       
      }
      if(ThisPDG==Polarity*321)
      {
	k4c->Fill(Momentum*arrayofsmearsp[jentry],1-y4[weightbin+1]);       
      }
      if(ThisPDG==Polarity*2212)
      {
	p4c->Fill(Momentum*arrayofsmearsp[jentry],1-y4[weightbin+1]);       
      }       
  }

   //Four Point HY
   pi4c->SetLineColor(9);
   mu4c->SetLineColor(41);
   e4c->SetLineColor(15);
   k4c->SetLineColor(7);
   p4c->SetLineColor(2);
   pi4c->SetFillColor(9);
   mu4c->SetFillColor(41);
   e4c->SetFillColor(15);
   k4c->SetFillColor(7);
   p4c->SetFillColor(2);
   double total4c=pi4c->Integral()+mu4c->Integral()+e4c->Integral()+k4c->Integral()+p4c->Integral();
   pi4c->Scale(1/total4c);
   mu4c->Scale(1/total4c);
   e4c->Scale(1/total4c);
   k4c->Scale(1/total4c);
   p4c->Scale(1/total4c);
   h4c->Add(p4c,"hist");
   h4c->Add(k4c,"hist");
   h4c->Add(e4c,"hist");
   h4c->Add(mu4c,"hist");
   h4c->Add(pi4c,"hist");
   
/*       //Make the Hy 4 Point Momentum Comparion Plot       
   auto c4c= new TCanvas("Corrected HY 4 Point Sim to Data Comparison","HY 4 Point Sim to Data Comparison");
   auto leg4c=new TLegend(0.5,0.4,0.9,0.9);
   leg4c->SetHeader("Corrected HY 4 Point Particle Species at WC4, Post WCFilter");
   leg4c->AddEntry(pi4c,TString::Format("#pi: %.1f %%",100*pi4c->Integral()));
   leg4c->AddEntry(mu4c,TString::Format("#mu: %.1f %%",100*mu4c->Integral()));
   leg4c->AddEntry(e4c,TString::Format("e: %.1f %%",100*e4c->Integral()));
   leg4c->AddEntry(k4c,TString::Format("k: %.1f %%",100*k4c->Integral()));
   leg4c->AddEntry(p4c,TString::Format("p: %.1f %%",100*p4c->Integral()));
   leg4c->AddEntry(d4,"+100A HY 4 Point Track, WCFiltered Data");
   h4c->Draw();
   d4->Draw("same");
   leg4c->Draw();
   h4c->GetXaxis()->SetTitle("Reconstructed Momentum (MeV/c)");
   h4c->GetYaxis()->SetTitle("Entries per 20 MeV/c");
   c4c->SetGridx();
   c4c->SetGridy();
    */
    
//Prepare the mass plot color scheme
   pi4m->SetLineColor(col1);
   mu4m->SetLineColor(col3);
   e4m->SetLineColor(col4);
   k4m->SetLineColor(col5);
   p4m->SetLineColor(col2);
   pi4m->SetFillColor(col1);
   mu4m->SetFillColor(col3);
   e4m->SetFillColor(col4);
   k4m->SetFillColor(col5);
   p4m->SetFillColor(col2);
   double total4m=pi4m->Integral()+mu4m->Integral()+e4m->Integral()+k4m->Integral()+p4m->Integral();
   
   pi4m->Scale(1/total4m);
   mu4m->Scale(1/total4m);
   e4m->Scale(1/total4m);
   k4m->Scale(1/total4m);
   p4m->Scale(1/total4m);
   d4m->Scale(1/d4m->Integral());
   h4m->Add(k4m,"hist");
   h4m->Add(e4m,"hist");
   h4m->Add(mu4m,"hist");
   h4m->Add(p4m,"hist");
   h4m->Add(pi4m,"hist");
      //Make the Hy 4 Point Mass comparison plot, using smeared TOF/mom      
   auto c4m= new TCanvas("HY 4 Point Sim Mass: Smeared Mom","HY 4 Point Sim Mass: Smeared Mom and TOF",800,800);
   auto leg4m=new TLegend(0.6,0.5,0.94,0.8);
   //leg4m->SetHeader("Mom and TOF Smeared, Post-WC Filter");
   leg4m->SetHeader("64 GeV +100 Amps");
   leg4m->AddEntry(pi4m,TString::Format("#pi: %.1f%%",100*pi4m->Integral()));
   leg4m->AddEntry(p4m,TString::Format("p: %.1f%%",100*p4m->Integral()));
   leg4m->AddEntry(mu4m,TString::Format("#mu: %.1f%%",100*mu4m->Integral()));
   leg4m->AddEntry(e4m,TString::Format("e: %.1f%%",100*e4m->Integral()));
   leg4m->AddEntry(k4m,TString::Format("K: %.1f%%",100*k4m->Integral()));
   leg4m->AddEntry(d4m,"Data");
   h4m->Draw();
   d4m->Draw("same");
   leg4m->Draw();
   leg4m->SetTextSize(textSize);
   h4m->GetYaxis()->SetRangeUser(0,1);
   h4m->GetXaxis()->SetTitle("Reconstructed Mass (MeV/c^{2})");
   h4m->GetYaxis()->SetTitle("Normalized Entries per 10 MeV/c^{2}");
   c4m->SetGridx();
   c4m->SetGridy();
   h4m->GetXaxis()->SetTitleOffset(1.1);
   h4m->GetYaxis()->SetTitleOffset(1.6);
   gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);   
   TPaveText* header4m=MakeTextBox(hd_x1,hd_y2,textSize,3,0.3);
   header4m->AddText("#bf{LArIAT}");
   header4m->AddText(Form("%.0f Data Events",d4int));
   header4m->AddText(Form("%.0f Simulated Events",total4));
   header4m->Draw();

   pi4t->SetLineColor(col1);
   mu4t->SetLineColor(col3);
   e4t->SetLineColor(col4);
   k4t->SetLineColor(col5);
   p4t->SetLineColor(col2);
   pi4t->SetFillColor(col1);
   mu4t->SetFillColor(col3);
   e4t->SetFillColor(col4);
   k4t->SetFillColor(col5);
   p4t->SetFillColor(col2);
   double total4t=pi4t->Integral()+mu4t->Integral()+e4t->Integral()+k4t->Integral()+p4t->Integral();
   std::cout<<total4t<<" "<<total4<<std::endl;
   pi4t->Scale(1/total4t);
   mu4t->Scale(1/total4t);
   e4t->Scale(1/total4t);
   k4t->Scale(1/total4t);
   p4t->Scale(1/total4t);
   d4t->Scale(1/d4t->Integral());
   h4t->Add(k4t,"hist");
   h4t->Add(e4t,"hist");
   h4t->Add(mu4t,"hist");
   h4t->Add(p4t,"hist");
   h4t->Add(pi4t,"hist");
      //Make the Hy 4 Point smeared tof plot       
   auto c4t= new TCanvas("HY 4 Point Sim TOF","HY 4 Point Sim TOF:",800,800);
   auto leg4t=new TLegend(0.6,0.5,0.94,0.8);
   leg4t->SetHeader("64 GeV +100 Amps");
   leg4t->AddEntry(pi4t,TString::Format("#pi: %.1f%%",100*pi4t->Integral()));
   leg4t->AddEntry(p4t,TString::Format("p: %.1f%%",100*p4t->Integral()));
   leg4t->AddEntry(mu4t,TString::Format("#mu: %.1f%%",100*mu4t->Integral()));
   leg4t->AddEntry(e4t,TString::Format("e: %.1f%%",100*e4t->Integral()));
   leg4t->AddEntry(k4t,TString::Format("K: %.1f%%",100*k4t->Integral()));
   leg4t->AddEntry(d4t,"Data");
   h4t->Draw();
   d4t->Draw("same");
   leg4t->Draw();
   leg4t->SetTextSize(textSize);
   h4t->GetXaxis()->SetTitle("Reconstructed TOF (ns)");
   h4t->GetYaxis()->SetTitle("Normalized entries per 0.25 ns");
   h4t->GetXaxis()->SetTitleOffset(1.1);
   h4t->GetYaxis()->SetTitleOffset(1.6);
   c4t->SetGridx();
   c4t->SetGridy();
   gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);   
   TPaveText* header4t=MakeTextBox(hd_x1,hd_y2,textSize,3,0.3);
   header4t->AddText("#bf{LArIAT}");
   header4t->AddText(Form("%.0f Data Events",d4int));
   header4t->AddText(Form("%.0f Simulated Events",total4));
   header4t->Draw();
   std::cout<<pi4->Integral()<<" "<<pi4t->Integral()<<std::endl;
   

  //Make the Bin by Bin Content Plot    
   binfracpi->SetName(""); 
   binfracpi->SetTitle("");
   gStyle->SetOptStat(0);
   auto c4cont= new TCanvas("pi/mu/e Content","pi/mu/e Content",800,800);
   auto leg4cont=new TLegend(0.6,0.5,0.94,0.8);
   leg4cont->SetHeader("64 GeV +100 Amps");
   leg4cont->AddEntry(binfracpi,TString::Format("#pi^{+}"));
   leg4cont->AddEntry(binfracmu,TString::Format("#mu^{+}"));
   leg4cont->AddEntry(binfrace,TString::Format("e^{+}"));
   leg4cont->SetTextSize(textSize);   
   binfracpi->GetXaxis()->SetTitle("Reco Momentum (MeV/c)");
   binfracpi->GetYaxis()->SetTitle("Fractional Content per 20 MeV/c"); 
   c4cont->SetGridx();
   c4cont->SetGridy();
   gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);   
   TPaveText* header4cont=MakeTextBox(hd_x1,hd_y2,textSize,3,0.3);
   header4cont->AddText("#bf{LArIAT}");
   header4cont->AddText(Form("%.0f Simulated Events",total4low));
   binfracpi->Draw();
   binfracmu->Draw("same");
   binfrace->Draw("same");
   header4cont->Draw("same");
   leg4cont->Draw("same");



   for (int i=0; i<bins4; i++)
   {
     if(binfracmu->GetBinContent(i)>1)
     {
       std::cout<<mu4->GetBinContent(i)<<" "<<mom4->GetBinContent(i)<<std::endl;
     }

   }


}
double BeamComp64GeVNeg100ARound2::FindPDG(double(&pdg)[4][2], double (&mother)[4][2], double (&id)[4][2])
{
  double pdg1=pdg[3][0];
  double pdg2=pdg[3][1];
  double id1=id[3][0];
  double id2=id[3][1];
  double mother1=mother[3][0];
  double mother2=mother[3][1];
  //If the particles are identical, just return the pdg.

  if (id1==id2){return pdg1;}
  //Consider whether the first pdg is a delta ray of the second particle
  else if(pdg1==11 && mother1==id2)
  {
    return pdg2;
  }
  else if(pdg2==11 && mother2==id1)
  {
    return pdg1;
  }  
  
  else
  {
    std::cout<<"Cannot find a causal relationship between these particles [trkid, pdg, mother] = ["<<id1<<", "<<pdg1<<", "<<mother1<<"], ["<<id2<<", "<<pdg2<<", "<<mother2<<"]. Returning -9999 and this event won't be included in the analysis."<<std::endl;
    return -9999;  
  }

}
