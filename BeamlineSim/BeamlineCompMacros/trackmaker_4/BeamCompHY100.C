#define BeamCompHY100SpillSPILLNUM_cxx
#include "BeamCompHY100SpillSPILLNUM.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom3.h>
int MaxHits=100;
//You must know the current setting used to make the BeamSim. Current should always be positive to give a positive B field.
double current=100;
double B=.3361;
double L=1143.7;

//Redefining the midplane not to have an angle, and instead be at a constant Z.
double magnet_midZ = (3801+4491)/2.;


/* double middiffXsigmascale=6.352/5.507; //From the data and sim distribution, the ratio of datasigma/simsigma. Will be used to smear the X distribution
double middiffXshift=2.871+.8505;  //From the data and sim distribution, meandata-meansim. Will be used to shift the X distribution before smearing.

double middiffYsigmascale=5.131/5.058; //From the data and sim distribution, the ratio of datasigma/simsigma. Will be used to smear the Y distribution
double middiffYshift=-.3327-.09;  //From the data and sim distribution, meandata-meansim. Will be used to shift the Y distribution before smearing. */

double middiffXsigmascale=1; //From the data and sim distribution, the ratio of datasigma/simsigma. Will be used to smear the X distribution
double middiffXshift=0;  //From the data and sim distribution, meandata-meansim. Will be used to shift the X distribution before smearing.

double middiffYsigmascale=1; //From the data and sim distribution, the ratio of datasigma/simsigma. Will be used to smear the Y distribution
double middiffYshift=0;  //From the data and sim distribution, meandata-meansim. Will be used to shift the Y distribution before smearing.

double TOFdeadtime=20E-9; //If a hit is registered in the TOF, how long to avoid another hit from being found. Standin for the resolution in pulse discrimination in the TOFhit finder.

double WCSmear=2; //How much smearing to apply to each hit in the WC to simulate noise. In mm.

//Declare/Initialize variables to be used/stored in the output tree.
std::vector<TVector3> wc1hits;
std::vector<TVector3> wc2hits;
std::vector<TVector3> wc3hits;
std::vector<TVector3> wc4hits;
std::vector<double>   USTOFhits;
std::vector<double>   DSTOFhits;
double XYZTrack[4][3]={{-9999,-9999,-9999},{-9999,-9999,-9999},{-9999,-9999,-9999},{-9999,-9999,-9999}};
std::vector<int> hitindex1;
std::vector<int> hitindex2;
std::vector<int> hitindex3;
std::vector<int> hitindex4;
std::vector<int> hitindexUS;
std::vector<int> hitindexDS;
int wcindex[4][2]={{-9999,-9999},{-9999,-9999},{-9999,-9999},{-9999,-9999}};
int hitindex[4][2]={{-9999,-9999},{-9999,-9999},{-9999,-9999},{-9999,-9999}};
int WCMissed=0;
double MinRes=-9999;
int perfectmatch=0;
int mismatched=0;
bool is4Pt=false;
bool isS2=false;
bool isS3=false;
int evtnum;
bool isFrankentrack=false;
bool isPure=false;
bool isPicky=false;
bool isHY=false;
bool wcfilt=false;
double momentum,tof,mass;
double error[2]={-9999,-9999};
int trackIDs[4][2]={{0,0},{0,0},{0,0},{0,0}};
int trackpdgs[4][2]={{0,0},{0,0},{0,0},{0,0}};
int MotherIDs[4][2]={{0,0},{0,0},{0,0},{0,0}};
double truemom[4][2]={{0,0},{0,0},{0,0},{0,0}};
double theta_x_us;
double theta_x_ds;
double yslope;
double MidDiff=-9999;
double MidDiffX=-9999;
double MidDiffY=-9999;
double WC4Diff=-9999;
double WC4DiffX=-9999;
double dshift=-9999;
double WC4DiffY=-9999;
double fDegtoRad=3.14159265359/180;
double WC1Mult, WC2Mult, WC3Mult,WC4Mult, WC1Miss;


//The bounds of the faces of the magnetic volumes. These will be used to remove tracks that project into the steel. [xlow_front, xhigh_front, xlow_back, xhigh_back]. Taken from input sim file.
 
double xboundMagnet1[4]={-1026.30, -725.524, -1111.326, -810.578};
double yboundMagnet1[4]={-70.2, 70.2, -70.2, 70.2};
  
double xboundMagnet2[4]={-1136.34,-831.847,-1181.073,-876.581};
double yboundMagnet2[4]={-70.0, 70.2, -70.2, 70.2};
 
double xboundDSCol[4]={-1170.41, -1019.59,-1170.41, -1019.59};
double yboundDSCol[4]={-76.53, 76.53, -76.53, 76.53}; 

//The z coordinate of the center of the face of the aperatures. [zcentus,zcentds]


double zcentMagnet1[2]={(3599.418+3543.672)/2, (4058.328+4002.582)/2};
double zcentMagnet2[2]={(4273.371+4244.052)/2, (4737.948+4708.629)/2};
double zcentDSCol[2]={(5595.379+5587.474)/2, (6508.526+6500.621)/2}; 


//The center of the WC, in G4BL coordinates

double xcent[4]={-403.0472, -738.0351, -1064.133, -1195.0827};
double ycent[4]={0.0508000, 0.0762000, -2.921000, -20.472400};
double zcent[4]={1730.3369, 3181.9215, 5167.5665, 7606.48720}; 

void BeamCompHY100SpillSPILLNUM::Loop()
{
TRandom3 *gRandom = new TRandom3();
std::string SpillNumber=("SpillSPILLNUM");
SpillNumber.erase(SpillNumber.begin(),SpillNumber.begin()+5);
std::cout<<"SpillNumber: "<<SpillNumber<<std::endl;
double SpillNumber1=std::stod(SpillNumber);
  TFile *outfile=new TFile("outdirbaseTrackTreeSpillSPILLNUM.root","recreate");

  TTree *tree= new TTree("TrackTree","TrackTree");
    tree->Branch("SpillNumber",&SpillNumber,"SpillNumber/D");
    tree->Branch("Momentum",&momentum,"momentum/D");
    tree->Branch("MomentumError",&error,"error[2]/D");
    tree->Branch("TrueHitMom",&truemom,"truemom[4][2]/D");
    tree->Branch("trackIDs",&trackIDs,"trackIDs[4][2]/I");
    tree->Branch("MotherIDs",&MotherIDs,"MotherIDs[4][2]/I");
    tree->Branch("trackpdgs",&trackpdgs,"trackpdgs[4][2]/I");
    tree->Branch("FourptTrack",&is4Pt,"is4Pt/O");
    tree->Branch("WCFilt",&wcfilt,"wcfilt/O");
    tree->Branch("S2Track",&isS2,"isS2/O");
    tree->Branch("S3Track",&isS3,"isS3/O");
    tree->Branch("PickyTrack",&isPicky,"isPicky/O");
    tree->Branch("HYTrack",&isHY,"isHY/O");
    tree->Branch("PureTrack",&isPure,"isPure/O");
    tree->Branch("Frankentrack",&isFrankentrack,"isFrankentrack/O");
    tree->Branch("thetaus",&theta_x_us,"theta_x_us/D");
    tree->Branch("thetads",&theta_x_ds,"theta_x_ds/D");
    tree->Branch("yslope",&yslope,"yslope/D");
    tree->Branch("MinRes",&MinRes,"MinRes/D");
    tree->Branch("MidDiff",&MidDiff,"MidDiff/D");
    tree->Branch("MidDiffX",&MidDiffX,"MidDiffX/D");
    tree->Branch("MidDiffY",&MidDiffY,"MidDiffY/D");
    tree->Branch("WC4Diff",&WC4Diff,"WC4Diff/D");
    tree->Branch("WC4DiffX",&WC4DiffX,"WC4DiffX/D");
    tree->Branch("WC4DiffY",&WC4DiffY,"WC4DiffY/D");    
    tree->Branch("WC1Mult",&WC1Mult,"WC1Mult/D");
    tree->Branch("WC2Mult",&WC2Mult,"WC2Mult/D");
    tree->Branch("WC3Mult",&WC3Mult,"WC3Mult/D");
    tree->Branch("WC4Mult",&WC4Mult,"WC4Mult/D");
    tree->Branch("WC1Miss",&WC1Miss,"WC1Miss/D");
    tree->Branch("XYZTrack",&XYZTrack,"XYZTrack[4][3]/D");
    tree->Branch("TOF",&tof, "tof/D");
    tree->Branch("Mass", &mass, "mass/D");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //For each event, get the WC and TOF hits and put them into a vector
      if(Event==jentry)
      {
        ResetVars();
        for(int i=0; i<MaxHits;i++)
	{
          if(tWC1[i]>0 && PDGWC1[i]!=22)
	  {
	  
	    double smearX1= gRandom->Gaus(0, WCSmear);
            double smearY1= gRandom->Gaus(0, WCSmear);
	    double smearZ1= gRandom->Gaus(0, WCSmear);
	    TVector3 v(xWC1[i]+smearX1,yWC1[i]+smearY1,zWC1[i]+smearZ1);
	    wc1hits.push_back(v);
	    hitindex1.push_back(i);
	  }
	  if(tWC2[i]>0 && PDGWC2[i]!=22)
	  {
	    double smearX2= gRandom->Gaus(0, WCSmear);
            double smearY2= gRandom->Gaus(0, WCSmear);
	    double smearZ2= gRandom->Gaus(0, WCSmear);
	  
	    TVector3 v(xWC2[i]+smearX2,yWC2[i]+smearY2,zWC2[i]+smearZ2);
	    wc2hits.push_back(v);
	    hitindex2.push_back(i);
	  }
	  if(tWC3[i]>0 && PDGWC3[i]!=22)
	  {
	    double smearX3= gRandom->Gaus(0, WCSmear);
            double smearY3= gRandom->Gaus(0, WCSmear);
	    double smearZ3= gRandom->Gaus(0, WCSmear);
	  
	    TVector3 v(xWC3[i]+smearX3,yWC3[i]+smearY3,zWC3[i]+smearZ3);
	    wc3hits.push_back(v);
	    hitindex3.push_back(i);
	  }
	  if(tWC4[i]>0 && PDGWC4[i]!=22)
	  {
	    double smearX4= gRandom->Gaus(0, WCSmear);
            double smearY4= gRandom->Gaus(0, WCSmear);
	    double smearZ4= gRandom->Gaus(0, WCSmear);
	  
	    TVector3 v(xWC4[i]+smearX4,yWC4[i]+smearY4,zWC4[i]+smearZ4);
	    wc4hits.push_back(v);
	    hitindex4.push_back(i);
	  }
	  if(tUSTOF[i]>0 && PDGUSTOF[i]!=22)
	  {
	    USTOFhits.push_back(tUSTOF[i]);
	    hitindexUS.push_back(i);
	  }
	  if(tDSTOF[i]>0 && PDGDSTOF[i]!=22)
	  {
	    DSTOFhits.push_back(tDSTOF[i]);
	    hitindexDS.push_back(i);
	  }		  	  	
 	}
	
	//See how many hits each WC has, and construct either a 3pt or 4pt track. If a 4pt track, reconstruct the TOF as well.
	WC1Mult=wc1hits.size();
	WC2Mult=wc2hits.size();
	WC3Mult=wc3hits.size();
	WC4Mult=wc4hits.size();
	if(wc2hits.size()==0 || wc3hits.size()==0)
	{
	  BuildThreePointTrack(wc1hits,wc2hits,wc3hits,wc4hits);
	}
	else
        {
	  momentum=BuildFourPointTrack(wc1hits,wc2hits,wc3hits,wc4hits);
	  tof=FindTOF(USTOFhits,DSTOFhits);
	  //If a momentum and tof are reconstructed, calculate the mass. If the TOF suggests a super luminal particle, the mass will be imaginary. Fill with negative mass instead of imaginary.
	  if (momentum>-9999 && tof>-9999)
	  {
	    double underroot=tof*0.299792458E9*0.299792458E9*tof/(6.652*6.652) - 1; 
	    if(underroot<0){underroot=-underroot;}
	    mass = momentum*pow(underroot ,0.5);
	    if((tof*0.299792458E9*0.299792458E9*tof/(6.652*6.652) - 1)<0)
	    {
	      mass=-mass;
	    }
	  }
	  else
	  {
	    mass=-9999;
	  }  
	}
      }
      if(isS2 || isS3 || is4Pt)
      {
        if(wc1hits.size()==1 && wc2hits.size()==1 && wc3hits.size()==1 && wc4hits.size()==1)
	{
	  isPicky=true;
	}
	else
	{
	  isHY=true;
	}
      }
      tree->Fill();
   }
  tree->Write();
  outfile->Close();
}
double BeamCompHY100SpillSPILLNUM::FindTOF(std::vector<double> ushits, std::vector<double> dshits)
{
  int nTOFs=0;
  double TOF=-9999;
  //sort the hits in each TOF by time
  std::sort(ushits.begin(),ushits.end());
  std::sort(dshits.begin(),dshits.end());
  //Loop over the TOF hits and remove any hits that are less than t=TOFdeadtime after another hit.
  for(int ius=ushits.size()-1; ius>0; ius--)
  {
    if((ushits.at(ius)-ushits.at(ius-1))<TOFdeadtime)
    {
      ushits.erase(ushits.begin()+ius-1);
    }    
  }
  for(int ids=dshits.size()-1; ids>0; ids--)
  {
    if((dshits.at(ids)-dshits.at(ids-1))<TOFdeadtime)
    {
      dshits.erase(dshits.begin()+ids-1);
    }    
  }
  for(int ius=0; ius<ushits.size(); ius++)
  {
    for(int ids=0; ids<dshits.size(); ids++)
    {
      if((dshits[ids]-ushits[ius]<100E-9) && (dshits[ids]-ushits[ius]>10E-9))
      {
        TOF=dshits[ids]-ushits[ius];
	nTOFs++;
      }
    }
  }
  if(nTOFs!=1)
  {
     TOF=-9999;
  }
  return TOF;
}
double BeamCompHY100SpillSPILLNUM::BuildFourPointTrack(std::vector<TVector3> wc1hits, std::vector<TVector3> wc2hits, std::vector<TVector3> wc3hits, std::vector<TVector3> wc4hits)
{
//Residual is the value we minimize to find the best combo of WC hits.
//After finding the line of best fit through a given combo of WC hits, find the average distance each point is from that best fit line. That's the Residual. If the residual is minimized, we've found most colinear set of points.
  std::vector<float> trackstats;
  double momentum=0;
  //Loop over every combination of points
  for(int iwc1z=0; iwc1z<wc1hits.size(); iwc1z++)
  {
    for(int iwc2z=0; iwc2z<wc2hits.size(); iwc2z++)
    {
      for(int iwc3z=0; iwc3z<wc3hits.size(); iwc3z++)
      {
        for(int iwc4z=0; iwc4z<wc4hits.size(); iwc4z++)
	{
	    for(int iwc1y=0; iwc1y<wc1hits.size(); iwc1y++)
            {
              for(int iwc2y=0; iwc2y<wc2hits.size(); iwc2y++)
              {
                for(int iwc3y=0; iwc3y<wc3hits.size(); iwc3y++)
                {
                  for(int iwc4y=0; iwc4y<wc4hits.size(); iwc4y++)
	          {
		    
	            //store the y,z hits in an array.
	            double y[4]={wc1hits[iwc1y].Y(),wc2hits[iwc2y].Y(),wc3hits[iwc3y].Y(),wc4hits[iwc4y].Y()};
	            double z[4]={wc1hits[iwc1z].Z(),wc2hits[iwc2z].Z(),wc3hits[iwc3z].Z(),wc4hits[iwc4z].Z()};
	            //Perform the linear regression, using "0" to indicate no missed WC. Using syntax from WCTrackBuilder code so I can copypaste
	            trackstats=Regression(y,z,WCMissed);
	            if(trackstats[2]<fabs(MinRes))
	            {
	              MinRes=trackstats[2];
	              yslope=trackstats[0];
	              //Store which hits are the ones to make a track
	              hitindex[0][0]=hitindex1[iwc1z];
	              hitindex[1][0]=hitindex2[iwc2z];
	              hitindex[2][0]=hitindex3[iwc3z];
	              hitindex[3][0]=hitindex4[iwc4z];
		      hitindex[0][1]=hitindex1[iwc1y];
	              hitindex[1][1]=hitindex2[iwc2y];
	              hitindex[2][1]=hitindex3[iwc3y];
	              hitindex[3][1]=hitindex4[iwc4y];
	              wcindex[0][0]=iwc1z;
	              wcindex[1][0]=iwc2z;
	              wcindex[2][0]=iwc3z;
	              wcindex[3][0]=iwc4z;
	              wcindex[0][1]=iwc1y;
	              wcindex[1][1]=iwc2y;
	              wcindex[2][1]=iwc3y;
	              wcindex[3][1]=iwc4y;		      
	            }
	          }
                }
              }
            }
	  }
	}
      }
    }  
    wcfilt=PassWCFilter(wcindex); //Check whether track would hit steel of magnets
    is4Pt=true;
    momentum=CalculateFourPointMomentum(wcindex,trackstats[0]);
  
    //We have the momentum, now lets see what particles actually went into the track.
    trackpdgs[0][0]=PDGWC1[hitindex[0][0]];
    trackpdgs[1][0]=PDGWC2[hitindex[1][0]];
    trackpdgs[2][0]=PDGWC3[hitindex[2][0]];
    trackpdgs[3][0]=PDGWC4[hitindex[3][0]];
    
    trackpdgs[0][1]=PDGWC1[hitindex[0][1]];
    trackpdgs[1][1]=PDGWC2[hitindex[1][1]];
    trackpdgs[2][1]=PDGWC3[hitindex[2][1]];
    trackpdgs[3][1]=PDGWC4[hitindex[3][1]];
        
    trackIDs[0][0]=TrackIDWC1[hitindex[0][0]];
    trackIDs[1][0]=TrackIDWC2[hitindex[1][0]];
    trackIDs[2][0]=TrackIDWC3[hitindex[2][0]];
    trackIDs[3][0]=TrackIDWC4[hitindex[3][0]];

    trackIDs[0][1]=TrackIDWC1[hitindex[0][1]];
    trackIDs[1][1]=TrackIDWC2[hitindex[1][1]];
    trackIDs[2][1]=TrackIDWC3[hitindex[2][1]];
    trackIDs[3][1]=TrackIDWC4[hitindex[3][1]];
        
    MotherIDs[0][0]=ParentIDWC1[hitindex[0][0]];
    MotherIDs[1][0]=ParentIDWC2[hitindex[1][0]]; 
    MotherIDs[2][0]=ParentIDWC3[hitindex[2][0]]; 
    MotherIDs[3][0]=ParentIDWC4[hitindex[3][0]]; 

    MotherIDs[0][1]=ParentIDWC1[hitindex[0][1]];
    MotherIDs[1][1]=ParentIDWC2[hitindex[1][1]]; 
    MotherIDs[2][1]=ParentIDWC3[hitindex[2][1]]; 
    MotherIDs[3][1]=ParentIDWC4[hitindex[3][1]]; 
        
    XYZTrack[0][0]=xWC1[hitindex[0][0]];
    XYZTrack[0][1]=yWC1[hitindex[0][1]];
    XYZTrack[0][2]=zWC1[hitindex[0][0]];
     
    XYZTrack[1][0]=xWC2[hitindex[1][0]];
    XYZTrack[1][1]=yWC2[hitindex[1][1]];
    XYZTrack[1][2]=zWC2[hitindex[1][0]];
     
    XYZTrack[2][0]=xWC3[hitindex[2][0]];
    XYZTrack[2][1]=yWC3[hitindex[2][1]];
    XYZTrack[2][2]=zWC3[hitindex[2][0]];
     
    XYZTrack[3][0]=xWC4[hitindex[3][0]];
    XYZTrack[3][1]=yWC4[hitindex[3][1]];
    XYZTrack[3][2]=zWC4[hitindex[3][0]];
    
   truemom[0][0]=pow(pow(PxWC1[hitindex[0][0]],2)+pow(PyWC1[hitindex[0][0]],2)+pow(PzWC1[hitindex[0][0]],2),0.5);
   truemom[1][0]=pow(pow(PxWC2[hitindex[1][0]],2)+pow(PyWC2[hitindex[1][0]],2)+pow(PzWC2[hitindex[1][0]],2),0.5);
   truemom[2][0]=pow(pow(PxWC3[hitindex[2][0]],2)+pow(PyWC3[hitindex[2][0]],2)+pow(PzWC3[hitindex[2][0]],2),0.5);
   truemom[3][0]=pow(pow(PxWC4[hitindex[3][0]],2)+pow(PyWC4[hitindex[3][0]],2)+pow(PzWC4[hitindex[3][0]],2),0.5);
   truemom[0][1]=pow(pow(PxWC1[hitindex[0][1]],2)+pow(PyWC1[hitindex[0][1]],2)+pow(PzWC1[hitindex[0][1]],2),0.5);
   truemom[1][1]=pow(pow(PxWC2[hitindex[1][1]],2)+pow(PyWC2[hitindex[1][1]],2)+pow(PzWC2[hitindex[1][1]],2),0.5);
   truemom[2][1]=pow(pow(PxWC3[hitindex[2][1]],2)+pow(PyWC3[hitindex[2][1]],2)+pow(PzWC3[hitindex[2][1]],2),0.5);
   truemom[3][1]=pow(pow(PxWC4[hitindex[3][1]],2)+pow(PyWC4[hitindex[3][1]],2)+pow(PzWC4[hitindex[3][1]],2),0.5);
   bool TrackStatus=PureFrankTrackStatus(trackIDs);
    if(TrackStatus)
    {
      isPure=true;
      perfectmatch++;
    }
    else
    { 
      isFrankentrack=true;
      mismatched++;
    }
    // Do I assess two errors, one per axis, or the average? One per axis for now.
     error[0]=100*(truemom[3][0]-momentum)/(truemom[3][0]);
     error[1]=100*(truemom[3][1]-momentum)/(truemom[3][1]);
    //Finding that tracks that used the same particle in WC2,3,4 but not WC1
/*     if(trackIDs[1]==trackIDs[2] && trackIDs[1]==trackIDs[3] && trackIDs[0]!=trackIDs[1])
    {
      std::cout<<"WC1 Breaking Frankentrack!"<<std::endl;
      double correctx=-9999;
      double incorrectx=-9999;
      //Find the WC1Xhit of the particle used for the WC1 hit in the track and the WC1Xhit of the particle used for the rest of the track, if it exists.
      for(int i=0; i<MaxHits;i++)
      {
        if(TrackIDWC1[i] == trackIDs[0])
	{
	  std::cout<<"Found the correct hit!"<<std::endl;
          correctx=xWC1[i];
	}
	if(TrackIDWC1[i]==trackIDs[1])
	{
	  std::cout<<"Found The Incorrect Hit!"<<std::endl;
	  incorrectx=xWC1[i];
	}
      }
      //If we found both particles present in WC1, find the difference in the xposition of these hits
      if(correctx!=-9999 && incorrectx!=-9999)
      {
        WC1Miss=incorrectx-correctx;
      }
      //If for some reason, the particle used in WC2,3,4 wasn't present in WC1 
      //(Like if a pion decayed after WC1 and the muon completed the track), store a nonsense number as a flag 
     
      else if(correctx!=-9999 && incorrectx==-9999)
      {
        WC1Miss=-123456;
      }	
      else{WC1Miss=-1234567;}//Please let us never get here.
    } */
  
  return momentum;
}


double BeamCompHY100SpillSPILLNUM::CalculateFourPointMomentum(int (&wcindex)[4][2], float yslope)
{
  double x[4]={wc1hits[wcindex[0][0]].X(),wc2hits[wcindex[1][0]].X(),wc3hits[wcindex[2][0]].X(),wc4hits[wcindex[3][0]].X()};
  double z[4]={wc1hits[wcindex[0][0]].Z(),wc2hits[wcindex[1][0]].Z(),wc3hits[wcindex[2][0]].Z(),wc4hits[wcindex[3][0]].Z()};
  float dx_us=x[1]-x[0];
  float dx_ds=x[3]-x[2];
  float dz_us=z[1]-z[0];
  float dz_ds=z[3]-z[2];
  theta_x_us= atan(dx_us/dz_us);
  theta_x_ds= atan(dx_ds/dz_ds);
  
  double mom=(fabs(B)*L)/((3.3356*(theta_x_ds-theta_x_us))*cos(atan(yslope)));
  return mom;
}

double BeamCompHY100SpillSPILLNUM::CalculateThreePointMomentum(int (&wcindex)[4][2], float yslope, int WCMissed)
{
  double momentum=0;
  if(WCMissed==2)
  {  
  
    std::vector<double> WC1Hit={wc1hits[wcindex[0][0]].X(),wc1hits[wcindex[0][1]].Y(),wc1hits[wcindex[0][0]].Z()};
    std::vector<double> WC3Hit={wc3hits[wcindex[2][0]].X(),wc3hits[wcindex[2][1]].Y(),wc3hits[wcindex[2][0]].Z()};  
    std::vector<double> WC4Hit={wc4hits[wcindex[3][0]].X(),wc4hits[wcindex[3][1]].Y(),wc4hits[wcindex[3][0]].Z()};
    
    float ds_dz=WC4Hit[2]-WC3Hit[2];
    float ds_dx=WC4Hit[0]-WC3Hit[0];

    std::vector<double> MidplaneHit=ProjToZ(WC3Hit,WC4Hit,magnet_midZ);
    float us_dz=(MidplaneHit[2]-WC1Hit[2]);
    float us_dx=(MidplaneHit[0]-WC1Hit[0]);       
    theta_x_us=atan(us_dx/us_dz);
    theta_x_ds=atan(ds_dx/ds_dz);
    double momentum=(fabs(B)*L)/((3.3356*(theta_x_ds-theta_x_us))*cos(atan(yslope)));
    double scalingfactor=0;
/*     if(B<0.3) //Only have correction factors for 60A and 100A. .3T will separate the samples, as that's like ~90A.
    {
       scalingfactor=1.21E-4*momentum+.0458;
       scalingfactor=0;
    }
    if(B>0.3)
    {
       scalingfactor=7.39E-5*momentum+.0479;
       scalingfactor=0;
    }
    momentum=momentum/(1+scalingfactor); */
  } 
    if(WCMissed==3)
  {
   
    
    std::vector<double> WC1Hit={wc1hits[wcindex[0][0]].X(),wc1hits[wcindex[0][1]].Y(),wc1hits[wcindex[0][0]].Z()};
    std::vector<double> WC2Hit={wc2hits[wcindex[1][0]].X(),wc2hits[wcindex[1][1]].Y(),wc2hits[wcindex[1][0]].Z()};  
    std::vector<double> WC4Hit={wc4hits[wcindex[3][0]].X(),wc4hits[wcindex[3][1]].Y(),wc4hits[wcindex[3][0]].Z()};
    float us_dz=WC2Hit[2]-WC1Hit[2];
    float us_dx=WC2Hit[0]-WC1Hit[0];

    std::vector<double> MidplaneHit=ProjToZ(WC1Hit,WC2Hit,magnet_midZ);	
    float ds_dz=-(MidplaneHit[2]-WC4Hit[2]);
    float ds_dx=-(MidplaneHit[0]-WC4Hit[0]);          
    theta_x_us=atan(us_dx/us_dz);
    theta_x_ds=atan(ds_dx/ds_dz);
    double momentum=(fabs(B)*L)/((3.3356*(theta_x_ds-theta_x_us))*cos(atan(yslope)));
    double scalingfactor=0;
/*     if(B<0.3) //Only have correction factors for 60A and 100A. .3T will separate the samples, as that's like ~90A.
    {
       scalingfactor=-5.71E-5*momentum-0.0483;
       scalingfactor=0;
    }
    if(B>0.3)
    {
       scalingfactor=-4.20E-5*momentum-0.0444;
       scalingfactor=0;
    }
    momentum=momentum/(1+scalingfactor); */
  }
  return momentum;
}
void BeamCompHY100SpillSPILLNUM::BuildThreePointTrack(std::vector<TVector3> wc1hits, std::vector<TVector3> wc2hits, std::vector<TVector3> wc3hits, std::vector<TVector3> wc4hits)
{
  std::vector<float> trackstats;
  int WCMissed;
  if(wc2hits.size()==0)
  {
    WCMissed=2;
    
  }
  if(wc3hits.size()==0)
  {
    WCMissed=3;
  }
  
  if(WCMissed==2)
  {
      //Loop over every combination of points
    for(int iwc1z=0; iwc1z<wc1hits.size(); iwc1z++)
    {
      for(int iwc3z=0; iwc3z<wc3hits.size(); iwc3z++)
      {
        for(int iwc4z=0; iwc4z<wc4hits.size(); iwc4z++)
  	{
	  for(int iwc1y=0; iwc1y<wc1hits.size(); iwc1y++)
          {
            for(int iwc3y=0; iwc3y<wc3hits.size(); iwc3y++)
            {
              for(int iwc4y=0; iwc4y<wc4hits.size(); iwc4y++)
  	      {
	        //store the y,z hits in an array.
	        double y[4]={wc1hits[iwc1y].Y(),-9999,wc3hits[iwc3y].Y(),wc4hits[iwc4y].Y()};
	        double z[4]={wc1hits[iwc1z].Z(),-9999,wc3hits[iwc3z].Z(),wc4hits[iwc4z].Z()};
	        //Perform the linear regression, using "0" to indicate no missed WC. Using syntax from WCTrackBuilder code so I can copypaste
	        trackstats=Regression(y,z,WCMissed);
	        if(trackstats[2]<fabs(MinRes))
	        {
	          MinRes=trackstats[2];
	          yslope=trackstats[0];
	      
	          //Store which hits to make a track, leaving WC2 as a nonsense number.
	          hitindex[0][0]=hitindex1[iwc1z];
	          hitindex[2][0]=hitindex3[iwc3z];
	          hitindex[3][0]=hitindex4[iwc4z];
		  hitindex[0][1]=hitindex1[iwc1y];
	          hitindex[2][1]=hitindex3[iwc3y];
	          hitindex[3][1]=hitindex4[iwc4y];
	          wcindex[0][0]=iwc1z;
	          wcindex[2][0]=iwc3z;
	          wcindex[3][0]=iwc4z;
	          wcindex[0][1]=iwc1y;
	          wcindex[2][1]=iwc3y;
	          wcindex[3][1]=iwc4y;		      	      
	        }
	      }
	    }
	  }
	}
      }
    }
  }
  if(WCMissed==3)
  {
      //Loop over every combination of points
    for(int iwc1z=0; iwc1z<wc1hits.size(); iwc1z++)
    {
      for(int iwc2z=0; iwc2z<wc2hits.size(); iwc2z++)
      {
        for(int iwc4z=0; iwc4z<wc4hits.size(); iwc4z++)
  	{
	  for(int iwc1y=0; iwc1y<wc1hits.size(); iwc1y++)
          {
            for(int iwc2y=0; iwc2y<wc2hits.size(); iwc2y++)
            {
              for(int iwc4y=0; iwc4y<wc4hits.size(); iwc4y++)
  	      {
	        //store the y,z hits in an array.
	        double y[4]={wc1hits[iwc1y].Y(),wc2hits[iwc2y].Y(),-9999,wc4hits[iwc4y].Y()};
	        double z[4]={wc1hits[iwc1z].Z(),wc2hits[iwc2z].Z(),-9999,wc4hits[iwc4z].Z()};
	        //Perform the linear regression, using "0" to indicate no missed WC. Using syntax from WCTrackBuilder code so I can copypaste
	        trackstats=Regression(y,z,WCMissed);
	        if(trackstats[2]<fabs(MinRes))
	        {
	          MinRes=trackstats[2];
	          yslope=trackstats[0];
	      
	          //Store which hits to make a track, leaving WC2 as a nonsense number.
	          hitindex[0][0]=hitindex1[iwc1z];
	          hitindex[1][0]=hitindex2[iwc2z];
	          hitindex[3][0]=hitindex4[iwc4z];
		  hitindex[0][1]=hitindex1[iwc1y];
	          hitindex[1][1]=hitindex2[iwc2y];
	          hitindex[3][1]=hitindex4[iwc4y];
	          wcindex[0][0]=iwc1z;
	          wcindex[1][0]=iwc2z;
	          wcindex[3][0]=iwc4z;
	          wcindex[0][1]=iwc1y;
	          wcindex[1][1]=iwc2y;
	          wcindex[3][1]=iwc4y;		      	      
	        }
	      }
	    }
	  }
	}
      }
    }
  }

     momentum=CalculateThreePointMomentum(wcindex, trackstats[0], WCMissed);
     wcfilt=PassWCFilter(wcindex, WCMissed, trackstats, momentum);//Check whether track would hit steel of magnets
      trackIDs[0][0]=TrackIDWC1[hitindex[0][0]];
      trackIDs[3][0]=TrackIDWC4[hitindex[3][0]];
      trackIDs[0][1]=TrackIDWC1[hitindex[0][1]];
      trackIDs[3][1]=TrackIDWC4[hitindex[3][1]];
      
      
      trackpdgs[0][0]=PDGWC1[hitindex[0][0]];
      trackpdgs[3][0]=PDGWC4[hitindex[3][0]];
      trackpdgs[0][1]=PDGWC1[hitindex[0][1]];
      trackpdgs[3][1]=PDGWC4[hitindex[3][1]];      
    if(WCMissed==2)
    {
      trackIDs[2][0]=TrackIDWC3[hitindex[2][0]];
      trackIDs[2][1]=TrackIDWC3[hitindex[2][1]];
      trackpdgs[2][0]=PDGWC3[hitindex[2][0]];
      trackpdgs[2][1]=PDGWC3[hitindex[2][1]];
      bool TrackStatus=PureFrankTrackStatus(trackIDs, 2);
      if(TrackStatus)
      {
       isPure=true;
       perfectmatch++;
      }
      else
      { 
        isFrankentrack=true; 
        mismatched++;
      }
      isS2=true;
      error[0]=100*(truemom[3][0]-momentum)/(truemom[3][0]);
      error[1]=100*(truemom[3][1]-momentum)/(truemom[3][1]);
 
    }  
    if(WCMissed==3)
    {
      trackIDs[1][0]=TrackIDWC2[hitindex[1][0]];
      trackIDs[1][1]=TrackIDWC2[hitindex[1][1]];
      trackpdgs[1][0]=PDGWC2[hitindex[1][0]];
      trackpdgs[1][1]=PDGWC2[hitindex[1][1]];
      bool TrackStatus=PureFrankTrackStatus(trackIDs, 3);
      if(TrackStatus)
      {
       isPure=true;
       perfectmatch++;
      }
      else
      { 
        isFrankentrack=true; 
        mismatched++;
      }
      isS3=true;
      error[0]=100*(truemom[3][0]-momentum)/(truemom[3][0]);
      error[1]=100*(truemom[3][1]-momentum)/(truemom[3][1]);
 
    }  
  
}

//Start here tomorrow.
std::vector<float> BeamCompHY100SpillSPILLNUM::Regression(double (&y)[4], double (&z)[4], int & WCMissed)
{
  std::vector<float> RegressionValues;
  int Npoints=0;
  float sum_z=0;
  float sum_zz=0;
  float sum_y=0;
  float sum_yz=0;
  float intercept=0;
  float slope=0;
  float residualsquare=0;
  float residual=0;
    for(int i=0; i<4; ++i){
      if(i != WCMissed-1){ //turning WC# to array index
        sum_y  += y[i];
        sum_zz += z[i]*z[i];
        sum_z  += z[i];
        sum_yz += y[i]*z[i];  
	++Npoints; //Residual is a function of number of points used.
      }
    }
  slope=(Npoints*sum_yz-sum_y*sum_z)/(Npoints*sum_zz-sum_z*sum_z);
  intercept=(sum_y-slope*sum_z)/Npoints;
  for(int i=0; i<4;++i){
    if(i != WCMissed-1){
    residual= (y[i]-slope*z[i]-intercept)/std::sqrt(1+slope*slope);
    residualsquare += (residual)*(residual); 
    }
  }
  float avgresidual= std::sqrt(residualsquare)/(Npoints-2);
  RegressionValues.push_back(slope);
  RegressionValues.push_back(intercept);
  RegressionValues.push_back(avgresidual);
  return RegressionValues;
} 
void BeamCompHY100SpillSPILLNUM::ResetVars()
{
  wc1hits.clear();
  wc2hits.clear();
  wc3hits.clear();
  wc4hits.clear();
  USTOFhits.clear();
  DSTOFhits.clear();
  WCMissed=0;
  MinRes=-9999;
  is4Pt=false;
  isS2=false;
  isS3=false;
  for(int i=0; i<2;i++)
  {
    error[i]=-9999;
  }
  momentum=-9999;
  tof=-9999;
  mass=-9999;
  isPicky=false;
  isHY=false;
  isPure=false;
  isFrankentrack=false;
  wcfilt=false;
  theta_x_us=-9999;
  theta_x_ds=-9999;
  yslope=-9999;
  MidDiff=-9999;
  MidDiffY=-9999;
  MidDiffX=-9999;  
  WC4Diff=-9999;
  WC4DiffY=-9999;
  WC4DiffX=-9999;  
  WC1Mult=-9999;
  WC2Mult=-9999;
  WC3Mult=-9999;
  WC4Mult=-9999;
  WC1Miss=-9999;
  for(int i=0; i<4; i++) 
  {
    for(int j=0; j<2; j++)
    {
      hitindex[i][j]=-9999;
      wcindex[i][j]=-9999;
      trackIDs[i][j]=-9999;
      trackpdgs[i][j]=-9999;
      MotherIDs[i][j]=-9999;
      truemom[i][j]=-9999;
    }  
    for(int j=0;j<3;j++)
    {
      XYZTrack[i][j]=-9999;   
    }
  }
}


bool BeamCompHY100SpillSPILLNUM::PassWCFilter(int (&wcindex)[4][2])
{

  std::vector<double> WC1Hit={wc1hits[wcindex[0][0]].X(),wc1hits[wcindex[0][1]].Y(),wc1hits[wcindex[0][0]].Z()};
  std::vector<double> WC2Hit={wc2hits[wcindex[1][0]].X(),wc2hits[wcindex[1][1]].Y(),wc2hits[wcindex[1][0]].Z()};
  std::vector<double> WC3Hit={wc3hits[wcindex[2][0]].X(),wc3hits[wcindex[2][1]].Y(),wc3hits[wcindex[2][0]].Z()};
  std::vector<double> WC4Hit={wc4hits[wcindex[3][0]].X(),wc4hits[wcindex[3][1]].Y(),wc4hits[wcindex[3][0]].Z()};

  bool PassUSMagnet=CheckUpstreamMagnetAperture(WC1Hit,WC2Hit);
  bool PassDSMagnet=CheckDownstreamMagnetAperture(WC3Hit,WC4Hit);
  bool PassDSCol=CheckDownstreamCollimatorAperture(WC3Hit,WC4Hit);
  bool PassMPtoWC4=MPtoWC4(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  bool PassMidPlane=MidplaneMatch(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  if(PassUSMagnet && PassDSMagnet && PassDSCol){return true;}
  else{return false;}
  
}
bool BeamCompHY100SpillSPILLNUM::PassWCFilter(int (&wcindex)[4][2], int WCMissed, std::vector<float> trackstats, double momentum)
{
  if(WCMissed==2)
  {
  std::vector<double> WC1Hit={wc1hits[wcindex[0][0]].X(),wc1hits[wcindex[0][1]].Y(),wc1hits[wcindex[0][0]].Z()};
  std::vector<double> WC3Hit={wc3hits[wcindex[2][0]].X(),wc3hits[wcindex[2][1]].Y(),wc3hits[wcindex[2][0]].Z()};  
  std::vector<double> WC4Hit={wc4hits[wcindex[3][0]].X(),wc4hits[wcindex[3][1]].Y(),wc4hits[wcindex[3][0]].Z()};
  
  float dx_ds=WC4Hit[0]-WC3Hit[0];
  float dz_ds=WC4Hit[2]-WC3Hit[2];
  float theta_x_ds= atan(dx_ds/dz_ds);
  float theta_x_us_reco=theta_x_ds-(fabs(B)*L/(3.3356*momentum*cos(atan(trackstats[0]))));
  float missedXwire=-(xcent[1]-WC1Hit[0]-(zcent[1]-WC1Hit[2])*tan(theta_x_us_reco))/(cos(fDegtoRad*13.0)-sin(fDegtoRad*13.0)*tan(theta_x_us_reco));
  float x2=xcent[1]+cos(fDegtoRad*13.0)*missedXwire;
  float z2=zcent[1]+sin(fDegtoRad*13.0)*missedXwire;
  float y2=trackstats[0]*z2+trackstats[1];
  //std::cout<<"extrapolated WC2 hit: "<<x2<<", "<<y2<<", "<<z2<<"."<<std::endl;
  std::vector<double> WC2Hit={x2,y2,z2};
  bool PassUSMagnet=CheckUpstreamMagnetAperture(WC1Hit,WC2Hit);
  bool PassDSMagnet=CheckDownstreamMagnetAperture(WC3Hit,WC4Hit);
  bool PassDSCol=CheckDownstreamCollimatorAperture(WC3Hit,WC4Hit);
  bool PassMPtoWC4=MPtoWC4(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  bool PassMidPlane=MidplaneMatch(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  //std::cout<<"This WC2 missed event....US: "<<PassUSMagnet<<" DS: "<<PassDSMagnet<<" Col: "<<PassDSCol<<std::endl;
  if(PassUSMagnet && PassDSMagnet && PassDSCol){return true;}
  else{return false;}
  } 
  
  
  else if(WCMissed==3)
  {
  std::vector<double> WC1Hit={wc1hits[wcindex[0][0]].X(),wc1hits[wcindex[0][1]].Y(),wc1hits[wcindex[0][0]].Z()};
  std::vector<double> WC2Hit={wc2hits[wcindex[1][0]].X(),wc2hits[wcindex[1][1]].Y(),wc2hits[wcindex[1][0]].Z()};
  std::vector<double> WC4Hit={wc4hits[wcindex[3][0]].X(),wc4hits[wcindex[3][1]].Y(),wc4hits[wcindex[3][0]].Z()};
  
  float dx_us=WC2Hit[0]-WC1Hit[0];
  float dz_us=WC2Hit[2]-WC1Hit[2];
  float theta_x_us= atan(dx_us/dz_us);
  float theta_x_ds_reco=theta_x_us+(fabs(B)*L/(3.3356*momentum*cos(atan(trackstats[0]))));
  float missedXwire=(WC4Hit[0]-xcent[2]-(WC4Hit[2]-zcent[2])*tan(theta_x_ds_reco))/(cos(fDegtoRad*3.0)-(sin(fDegtoRad*3.0))*tan(theta_x_ds_reco));
  float x3=xcent[2]+cos(fDegtoRad*3.0)*missedXwire;
  float z3=zcent[2]+sin(fDegtoRad*3.0)*missedXwire;
  float y3=trackstats[0]*z3+trackstats[1];

  //std::cout<<"extrapolated WC3 hit: "<<x3<<", "<<y3<<", "<<z3<<"."<<std::endl;
  std::vector<double> WC3Hit={x3,y3,z3};
  bool PassUSMagnet=CheckUpstreamMagnetAperture(WC1Hit,WC2Hit);
  bool PassDSMagnet=CheckDownstreamMagnetAperture(WC3Hit,WC4Hit);
  bool PassDSCol=CheckDownstreamCollimatorAperture(WC3Hit,WC4Hit);
  bool PassMPtoWC4=MPtoWC4(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  bool PassMidPlane=MidplaneMatch(WC1Hit,WC2Hit,WC3Hit,WC4Hit);
  //std::cout<<"This WC3 missed event....US: "<<PassUSMagnet<<" DS: "<<PassDSMagnet<<" Col: "<<PassDSCol<<std::endl;
  if(PassUSMagnet && PassDSMagnet && PassDSCol){return true;}
  else{return false;}
  }
 else {std::cout<<"This isn't right!"<<std::endl; return false;}
}
std::vector<double> BeamCompHY100SpillSPILLNUM::ProjToZ(std::vector<double> hit0, std::vector<double> hit1, double zpos) {

  //
  // This powerful little piece of code is what is used to find the point at z = zpos along 
  //   the line created by hit0 and hit1. Uses the parameterized vector form of a line
  //
  // <x, y, z> = <sx, sy, sz> * t + <startx, starty, startz>
  // sx, sy, sz are all slopes
  //
  // (z - startz) / sz = t
  // x = sx * t + startx
  // y = sx * t + srarty
  //
      
  double sx = hit1[0] - hit0[0];
  double sy = hit1[1] - hit0[1];
  double sz = hit1[2] - hit0[2];

  double t = (zpos - hit0[2]) / sz;

  std::vector<double> resulto {sx * t + hit0[0], sy * t + hit0[1], zpos};

  return resulto;

}
//COLLIMATOR CHECKS! MAGNET 1
//===================================================================================================
bool BeamCompHY100SpillSPILLNUM::CheckUpstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2)
{
   std::vector<double> USHit=ProjToZ(hit1,hit2,zcentMagnet1[0]);
   std::vector<double> DSHit=ProjToZ(hit1,hit2,zcentMagnet1[1]);
//Upstream Aperture X Check
   if ( USHit[0] < xboundMagnet1[0] || USHit[0] > xboundMagnet1[1] )
   {
     return false;
   }
//Upstream Aperture Y Check
   else if ( USHit[1] < yboundMagnet1[0] || USHit[1] > yboundMagnet1[1] )
   {
     return false;
   }
//Downstream Aperture Y Check   
   else if ( DSHit[1] < yboundMagnet1[2] || DSHit[1] > yboundMagnet1[3] )
   {
     return false;
   }
   else {return true;} //If all checks pass, then we're good for this collimator.
}

//MAGNET 2
//===================================================================================================
bool BeamCompHY100SpillSPILLNUM::CheckDownstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2)
{
   std::vector<double> USHit=ProjToZ(hit1,hit2,zcentMagnet2[0]);
   std::vector<double> DSHit=ProjToZ(hit1,hit2,zcentMagnet2[1]);
//Upstream Aperture Y Check
   if ( USHit[1] < yboundMagnet2[0] || USHit[1] > yboundMagnet2[1] )
   {
     return false;
   } 
//Downstream Aperture X Check    
   else if ( DSHit[0] < xboundMagnet2[2] || DSHit[0] > xboundMagnet2[3] )
   {
   
     return false;
   }   
//Downstream Aperture Y Check   
   else if ( DSHit[1] < yboundMagnet2[2] || DSHit[1] > yboundMagnet2[3] )
   {
   
     return false;
   }  

   else {return true;} //If all checks pass, then we're good for this collimator.
}
//DS Collimator
//===================================================================================================
bool BeamCompHY100SpillSPILLNUM::CheckDownstreamCollimatorAperture(std::vector<double> hit1, std::vector<double> hit2)
{
   std::vector<double> USHit=ProjToZ(hit1,hit2,zcentDSCol[0]);
   std::vector<double> DSHit=ProjToZ(hit1,hit2,zcentDSCol[1]);
 //Upstream Aperture X Check   
   if ( USHit[0] < xboundDSCol[0] || USHit[0] > xboundDSCol[1] )
   {

     return false;
   } 
 //Upstream Aperture Y Check   
   else if ( USHit[1] < yboundDSCol[0] || USHit[1] > yboundDSCol[1] )
   {

     return false;
   } 
 //Downstream Aperture X Check    
   else if ( DSHit[0] < xboundDSCol[2] || DSHit[0] > xboundDSCol[3] )
   {

     return false;
   }  
//Downstream Aperture Y Check   
   else if ( DSHit[1] < yboundDSCol[2] || DSHit[1] > yboundDSCol[3] )
   {

     return false;
   } 

   else {return true;} //If all checks pass, then we're good for this collimator.
}
//==================================================================================================================================
bool BeamCompHY100SpillSPILLNUM::MidplaneMatch(std::vector<double> hit1, std::vector<double> hit2, std::vector<double> hit3, std::vector<double> hit4)
{
  // This is horribly hardcoded stuff, but magnet_midZ is the average of the z positions of magnet 1 and 2
  std::vector<double> USLegProj=ProjToZ(hit1,hit2, magnet_midZ);
  std::vector<double> DSLegProj=ProjToZ(hit3,hit4, magnet_midZ);

  MidDiffX=middiffXsigmascale*(USLegProj[0]-DSLegProj[0])+middiffXshift;
  MidDiffY=middiffYsigmascale*(USLegProj[1]-DSLegProj[1])+middiffYshift; 
  
  MidDiff=pow(pow(MidDiffX,2)+pow(MidDiffY,2),.5);

  return true; //Until I figure out a cut value, always pass.

}
//=============================================================================================================================
bool BeamCompHY100SpillSPILLNUM::MPtoWC4(std::vector<double> hit1, std::vector<double> hit2, std::vector<double> hit3, std::vector<double> hit4)
{

  std::vector<double> USLegProj=ProjToZ(hit1,hit2, magnet_midZ);
  //Take MP hit from the US leg projection, use it with WC3 to project to the z coordinate where the wc4 center. 
  std::vector<double> MPtoWC4Proj=ProjToZ(USLegProj,hit3, zcent[3]);
  
  WC4Diff  =pow(pow(MPtoWC4Proj[0]-hit4[0],2)+pow(MPtoWC4Proj[1]-hit4[1],2),.5);
  WC4DiffX = MPtoWC4Proj[0]-hit4[0];
  WC4DiffY = MPtoWC4Proj[1]-hit4[1];
  return true; //Until I figure out a cut value, always pass.
}
//=======================================================================================================================

bool BeamCompHY100SpillSPILLNUM::PureFrankTrackStatus(int (&ID)[4][2])
{
  int FirstID=ID[0][0];
  for (int i=0; i<4; i++)
  {
    for (int j=0;j<2; j++)
    {
      if(ID[i][j]!=FirstID)
      {
        return false;
      } 
    }
  }
  return true;
}
//======================================================================================================================
bool BeamCompHY100SpillSPILLNUM::PureFrankTrackStatus(int (&ID)[4][2], int WCMissed)
{
 int FirstID=ID[0][0];
  for (int i=0; i<4; i++)
  {
    if(i!=WCMissed-1)
    {
      for (int j=0;j<2; j++)
      {
        if(ID[i][j]!=FirstID)
        {
          return false;
        } 
      }
    }
  }
  return true;
}
