#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <algorithm>
#include <vector>
//Declare the arrays that will put into the outtree of particles found within a readout buffer of a trigger.
double xWC1[100], yWC1[100], zWC1[100], tWC1[100], PDGWC1[100], PxWC1[100], PyWC1[100], PzWC1[100], ParentIDWC1[100], TrackIDWC1[100],PDGIDWC1[100], rtWC1[100], EventIDWC1[100];
double xWC2[100], yWC2[100], zWC2[100], tWC2[100], PDGWC2[100], PxWC2[100], PyWC2[100], PzWC2[100], ParentIDWC2[100], TrackIDWC2[100],PDGIDWC2[100], rtWC2[100], EventIDWC2[100];
double xWC3[100], yWC3[100], zWC3[100], tWC3[100], PDGWC3[100], PxWC3[100], PyWC3[100], PzWC3[100], ParentIDWC3[100], TrackIDWC3[100],PDGIDWC3[100], rtWC3[100], EventIDWC3[100];
double xWC4[100], yWC4[100], zWC4[100], tWC4[100], PDGWC4[100], PxWC4[100], PyWC4[100], PzWC4[100], ParentIDWC4[100], TrackIDWC4[100],PDGIDWC4[100], rtWC4[100], EventIDWC4[100];
double xUSTOF[100], yUSTOF[100], zUSTOF[100], tUSTOF[100], PDGUSTOF[100], PxUSTOF[100], PyUSTOF[100], PzUSTOF[100], ParentIDUSTOF[100], TrackIDUSTOF[100],PDGIDUSTOF[100], rtUSTOF[100], EventIDUSTOF[100];
double xDSTOF[100], yDSTOF[100], zDSTOF[100], tDSTOF[100], PDGDSTOF[100], PxDSTOF[100], PyDSTOF[100], PzDSTOF[100], ParentIDDSTOF[100], TrackIDDSTOF[100],PDGIDDSTOF[100], rtDSTOF[100], EventIDDSTOF[100];
double coincidencewindow[10][200];
int evt; 
Double_t timewindow=100E-9; //100 ns to allow a trigger to occur
Double_t buffertime=300E-6; //300 us readout to grab particles
Double_t offset=0;
Double_t WCbufferstart=-300*1.117E-9;
Double_t WCbufferend=724*1.117E-9;
Int_t wc1iter,wc2iter,wc3iter,wc4iter,USTOFiter,DSTOFiter;
double n45, n345, n2345, n12345, n012345, n1345, n01345, n245, n1245, n01245, spill;
//Declare the output root file and all the branches that will go into it.


using namespace std;

// Tree variables and branches definition

TTree          *fChain;
Int_t           fCurrent;

Float_t         x;
Float_t         y;
Float_t         z;
Float_t         Px;
Float_t         Py;
Float_t         Pz;
Float_t         t;
Float_t         PDGid;
Float_t         EventID;
Float_t         TrackID;
Float_t         ParentID;
Float_t         Weight;

TBranch        *b_x;   //!
TBranch        *b_y;   //!
TBranch        *b_z;   //!
TBranch        *b_Px;   //!
TBranch        *b_Py;   //!
TBranch        *b_Pz;   //!
TBranch        *b_t;   //!
TBranch        *b_PDGid;   //!
TBranch        *b_EventID;   //!
TBranch        *b_TrackID;   //!
TBranch        *b_ParentID;   //!
TBranch        *b_Weight;   //!




// struct will be used to keep all hits in the memory
struct beam_dethit {
  Double_t         x;
  Double_t         y;
  Double_t         z;
  Double_t         Px;
  Double_t         Py;
  Double_t         Pz;
  Double_t         t;
  Double_t         PDGid;
  Double_t         EventID;
  Double_t         TrackID;
  Double_t         ParentID;
  Double_t         spill_time;
  Double_t         Weight;
};

Double_t       beam_time[350000]={0}; // Keep random times for all events
Double_t       beam_timesave[350000]; // Keep random times for all events 
// vector containing all hits for all detectors
beam_dethit    detHits[6][1000000];  // vector of vectors nDet*nHit
 //blank hit list for one det


  
Int_t          detNumHits[6]; // Number of hits per detector

// Properties of the detector (4 WC and 2 TOF)  signals: width and delay

Double_t       detSigDelay[6]={18e-9, 15e-9, 10e-9, 5e-9, 0e-9, 0e-9};
Double_t       coincWidth = 10e-9;
Double_t       detSigWidth[6]={100e-9, 100e-9, 100e-9, 100e-9, 100e-9, 100e-9};
// Diagnostic histograms
TH1D *BucketRan = new TH1D("BucketRan", "", 84, 0, 84);
TH1D *BatchRan = new TH1D("BatchRan", "", 7, 0, 7);
TH1D *OrbitRan = new TH1D("OrbitRan", "", 379939, 0, 379939);
TH1D *TimeRan = new TH1D("TimeRan", "", 400, -2.2e-9, 2.2e-9);
TH1D *evid = new TH1D("evid", "", 400000, -100, 399900);
TH1D *TOFtime[6];// = new TH1D("TOFtime", "", 420000, 0, 4.2);
TH1D *SortTest = new TH1D("SortTest", "", 420000, 0, 4.2);

void Init(TTree *tree){
   // initializes the tree
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("Px", &Px, &b_Px);
   fChain->SetBranchAddress("Py", &Py, &b_Py);
   fChain->SetBranchAddress("Pz", &Pz, &b_Pz);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("PDGid", &PDGid, &b_PDGid);
   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("TrackID", &TrackID, &b_TrackID);
   fChain->SetBranchAddress("ParentID", &ParentID, &b_ParentID);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);

}
//A function that resets all the arrays for the tree before we process the next trigger
void ResetTriggerTree()
{  
  for(int i=0; i<100; i++)
  { 
    xWC1[i]=-9999;
    yWC1[i]=-9999;
    zWC1[i]=-9999;
    tWC1[i]=-9999;
    PDGWC1[i]=-9999;
    PxWC1[i]=-9999; 
    PyWC1[i]=-9999; 
    PzWC1[i]=-9999;
    ParentIDWC1[i]=-9999;
    TrackIDWC1[i]=-9999;
    EventIDWC1[i]=-9999;
    
    xWC2[i]=-9999;
    yWC2[i]=-9999;
    zWC2[i]=-9999;
    tWC2[i]=-9999;
    PDGWC2[i]=-9999;
    PxWC2[i]=-9999; 
    PyWC2[i]=-9999; 
    PzWC2[i]=-9999;
    ParentIDWC2[i]=-9999;
    TrackIDWC2[i]=-9999;
    EventIDWC2[i]=-9999;
    
    xWC3[i]=-9999;
    yWC3[i]=-9999;
    zWC3[i]=-9999;
    tWC3[i]=-9999;
    PDGWC3[i]=-9999;
    PxWC3[i]=-9999; 
    PyWC3[i]=-9999; 
    PzWC3[i]=-9999;
    ParentIDWC3[i]=-9999;
    TrackIDWC3[i]=-9999;
    EventIDWC3[i]=-9999;
    
    xWC4[i]=-9999;
    yWC4[i]=-9999;
    zWC4[i]=-9999;
    tWC4[i]=-9999;
    PDGWC4[i]=-9999;
    PxWC4[i]=-9999; 
    PyWC4[i]=-9999; 
    PzWC4[i]=-9999;
    ParentIDWC4[i]=-9999;
    TrackIDWC4[i]=-9999;
    EventIDWC4[i]=-9999;

    xUSTOF[i]=-9999;
    yUSTOF[i]=-9999;
    zUSTOF[i]=-9999;
    tUSTOF[i]=-9999;
    PDGUSTOF[i]=-9999;
    PxUSTOF[i]=-9999; 
    PyUSTOF[i]=-9999; 
    PzUSTOF[i]=-9999;
    ParentIDUSTOF[i]=-9999;
    TrackIDUSTOF[i]=-9999;
    EventIDUSTOF[i]=-9999;
    
    xDSTOF[i]=-9999;
    yDSTOF[i]=-9999;
    zDSTOF[i]=-9999;
    tDSTOF[i]=-9999;
    PDGDSTOF[i]=-9999;
    PxDSTOF[i]=-9999; 
    PyDSTOF[i]=-9999; 
    PzDSTOF[i]=-9999;
    ParentIDDSTOF[i]=-9999;
    TrackIDDSTOF[i]=-9999;
    EventIDDSTOF[i]=-9999;
            
    wc1iter=0;
    wc2iter=0;
    wc3iter=0;
    wc4iter=0;
    USTOFiter=0;
    DSTOFiter=0;    
    rtWC1[i]=-9999;
    rtWC2[i]=-9999;
    rtWC3[i]=-9999;
    rtWC4[i]=-9999;
    rtUSTOF[i]=-9999;
    rtDSTOF[i]=-9999;
    n45=-1;
    n345=-1;
    n2345=-1;
    n12345=-1;
    n012345=-1;
    n1345=-1;
    n01345=-1;
    n245=-1;
    n1245=-1;
    n01245=-1;
    spill=-1;
       
  }

}  
void ReadTree(Int_t idet, TTree *tree){
  // read a detector tree and save the info in the struct variable

  TOFtime[idet] = new TH1D(Form("TOFtime_%i",idet), "", 420000, 0, 4.2);
  Init(tree);

  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  std::cout << nentries << std::endl;
  detNumHits[idet] = nentries;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = tree->LoadTree(jentry);
    if(ientry < 0) break;
    nb = tree->GetEntry(jentry);
    Int_t evnum = (Int_t) EventID - (ShiftSPILLNUM-1)+1; //evnum must be between 0 and 350k, but EventID is between 350k(SPILLNUM-1) and 350k(SPILLNUM). Use the spillnum to scale down to the 0:350k range 
    if(evnum>350000+(ShiftSPILLNUM-1)+1){std::cout<<"Weird EventNum: "<<evnum<<std::endl;}
    if(evnum<0){std::cout<<"Weird EventNum: "<<evnum<<std::endl;}
    detHits[idet][jentry].x           = (Double_t) x;
    detHits[idet][jentry].y           = (Double_t) y;
    detHits[idet][jentry].z           = (Double_t) z;
    detHits[idet][jentry].t           = (Double_t) t;
    detHits[idet][jentry].Px          = (Double_t) Px;
    detHits[idet][jentry].Py          = (Double_t) Py;
    detHits[idet][jentry].Pz          = (Double_t) Pz;
    detHits[idet][jentry].PDGid       = (Double_t) PDGid;
    detHits[idet][jentry].EventID     = (Double_t) evnum;
    detHits[idet][jentry].TrackID     = (Double_t) TrackID;
    detHits[idet][jentry].ParentID    = (Double_t) ParentID;
    detHits[idet][jentry].Weight      = (Double_t) Weight; 
    
    
  
    // calculate the hit time based on event ID
    detHits[idet][jentry].spill_time  = (Double_t)(t*1e-9 + beam_time[evnum]);
    TOFtime[idet]->Fill(detHits[idet][jentry].spill_time);
    
  }
  std::cout << "Read detector " << idet <<std::endl;

}

bool timecomp(beam_dethit lhs, beam_dethit rhs){return lhs.spill_time < rhs.spill_time;} // needed to sort the hits

Int_t FindDetCoinc(Int_t detID, Double_t coincLow, Double_t coincHigh, Int_t detCoincCandID[100]){
   Double_t cut=0; //A lot of hits are at the beam left side of the detector, TOFUS excluded. Ignore hits in this weird peak by setting a detector dependent cut. 
   
   // Finding hits within a time window defined by coincLow and coincHigh 
   Int_t detCoincCandNum = 0;
   Bool_t detTrig = false;
   Int_t detTrigID = 0;

   for(Int_t idet = 0; idet < detNumHits[detID]; idet++){

     // skip looping over hits from photons, neutrinos, and neutrons
     Int_t PDGID = detHits[detID][idet].PDGid;
     Bool_t neutralCheck = (PDGID == 22) || (PDGID == 12) ||  (PDGID == 14);// || (PDGID == 2112);
   
     if(neutralCheck) continue;

     // skip all hits with times earlier than the signal end
     if(detTrig && (detHits[detID][idet].spill_time <= (detHits[detID][detTrigID].spill_time + detSigWidth[detID]))) continue;
     detTrig = false;

     // calculate the time of a hit including a delay
     Double_t detHitTime = detHits[detID][idet].spill_time + detSigDelay[detID];
  
     // if the hit has passed the time window boundary stop the loop
     if(detHitTime >= (coincHigh + (detSigWidth[detID] - coincWidth))) break;

     // chack if the hit signal overlaps with the time window for at least the coincidence width
     Bool_t detHitCand = (detHitTime >= (coincLow - (detSigWidth[detID] - coincWidth))) && (detHitTime <= (coincHigh + (detSigWidth[detID] - coincWidth) ));
     if(detHitCand){

        // we have a signal that overlaps so check the size of the overlap
        Double_t testLow  = std::max(coincLow, detHitTime);
        Double_t testHigh = std::min(coincHigh, detHitTime + detSigWidth[detID]);


        if((testHigh - testLow) >= coincWidth){

          // save this hit ID and set the flag to skip following hits within the signal width
          detCoincCandID[detCoincCandNum] = idet;
          detTrig = true;
          detTrigID = idet;
          detCoincCandNum++;

        }
     }

   }
   // return the number of hits overlapping with the time window
   return detCoincCandNum; 

}

void FindTriggers(){
 
  Int_t trig45Num = 0;
  Int_t trig345Num = 0;
  Int_t trig2345Num = 0;
  Int_t trig12345Num = 0;
  Int_t trig012345Num = 0;
  Int_t trig1345Num = 0;
  Int_t trig01345Num = 0;
  Int_t trigNum =  0;
  Int_t trig245Num = 0;
  Int_t trig1245Num = 0;
  Int_t trig01245Num = 0;
  std::vector<Double_t> triggertimes;
  std::vector<Double_t> TOFdshit; 
  //  loop over all DSTOF hits
  for(Int_t idet5 = 0; idet5 < detNumHits[5]; idet5++){

    if(detHits[5][idet5].PDGid == 22) continue;    

    // open 100ns gate starting with DSTOF hit
    // skip consequent hits within the 100ns gate
    if(idet5>0 && ((detHits[5][idet5].spill_time - detHits[5][idet5-1].spill_time) < detSigWidth[5])) continue;

    // to start coincidence window conincides with the width of the gate
    Double_t CoincLow5  = detHits[5][idet5].spill_time;
    Double_t CoincHigh5 = CoincLow5 + detSigWidth[5];    
    Int_t det4HitCoinc[100]={0};
    
    // find WC4 hits within this gate. returns number of hits and a vector with hit ID
    Double_t numDet4Coinc = FindDetCoinc(4, CoincLow5, CoincHigh5, det4HitCoinc);
    trig45Num += numDet4Coinc;
    for(Int_t idet4 = 0; idet4 < numDet4Coinc; idet4++){
      
      // reduce the coincidence window based 
      Double_t CoincLow4 = std::max(CoincLow5,detHits[4][det4HitCoinc[idet4]].spill_time+detSigDelay[4]);
      Double_t CoincHigh4 = std::min(CoincHigh5,detHits[4][det4HitCoinc[idet4]].spill_time + detSigWidth[4]+detSigDelay[4]);
      Int_t det3HitCoinc[100]={0};  
          
      // find WC3 hits within possible coincidence windows. returns number of hits and a vector with hit ID
      Double_t numDet3Coinc = FindDetCoinc(3, CoincLow4, CoincHigh4, det3HitCoinc);
      trig345Num += numDet3Coinc;
           
      for(Int_t idet3 = 0; idet3 < numDet3Coinc; idet3++){

        // reduce the coincidence window based
        Double_t CoincLow3 = std::max(CoincLow4, detHits[3][det3HitCoinc[idet3]].spill_time+detSigDelay[3]);
        Double_t CoincHigh3 = std::min(CoincHigh4, detHits[3][det3HitCoinc[idet3]].spill_time + detSigDelay[3]+detSigWidth[3]);
        Int_t det2HitCoinc[100]={0};
	
        // find WC2 hits within possible coincidence windows. returns number of hits and a vector with hit ID
        Double_t numDet2Coinc = FindDetCoinc(2, CoincLow3, CoincHigh3, det2HitCoinc);
        trig2345Num += numDet2Coinc;
        
        for(Int_t idet2 = 0; idet2 < numDet2Coinc; idet2++){
          // reduce the coincidence window based
          Double_t CoincLow2 = std::max(CoincLow3, detHits[2][det2HitCoinc[idet2]].spill_time+detSigDelay[2]);
          Double_t CoincHigh2 = std::min(CoincHigh3, detHits[2][det2HitCoinc[idet2]].spill_time + detSigWidth[2]+detSigDelay[2]);

          Int_t det1HitCoinc[100]={0};
          // find WC1 hits within possible coincidence windows. returns number of hits and a vector with hit ID
          Double_t numDet1Coinc = FindDetCoinc(1, CoincLow2, CoincHigh2, det1HitCoinc);
          trig12345Num += numDet1Coinc;
          
          for(Int_t idet1 = 0; idet1 < numDet1Coinc; idet1++){
            // reduce the coincidence window based
            Double_t CoincLow1 = std::max(CoincLow2, detHits[1][det1HitCoinc[idet1]].spill_time+detSigDelay[1]);
            Double_t CoincHigh1 = std::min(CoincHigh2, detHits[1][det1HitCoinc[idet1]].spill_time + detSigWidth[1]+detSigDelay[1]);

            Int_t det0HitCoinc[100]={0};
            // find USTOF hits within possible coincidence windows. returns number of hits and a vector with hit ID
            Double_t numDet0Coinc = FindDetCoinc(0, CoincLow1, CoincHigh1, det0HitCoinc);
            trig012345Num += numDet0Coinc;
           
            for(Int_t idet0 = 0; idet0 < numDet0Coinc; idet0++){
	      Double_t detHitTime0 = detHits[0][det0HitCoinc[idet0]].spill_time + detSigDelay[0]; 
              Double_t trigTimeLow = std::max(CoincLow1, detHitTime0);
              Double_t trigTimeHigh = std::min(CoincHigh1, detHitTime0 + detSigWidth[0]);

	      triggertimes.push_back(trigTimeLow);
	      TOFdshit.push_back(detHits[5][idet5].TrackID);
	    }
          }          
        }
      }
    }
  }
  
    //  S3 Trigger condition
  for(Int_t idet5 = 0; idet5 < detNumHits[5]; idet5++){

    if(detHits[5][idet5].PDGid == 22) continue;    

    // open 100ns gate starting with DSTOF hit
    // skip consequent hits within the 100ns gate
    if(idet5>0 && ((detHits[5][idet5].spill_time - detHits[5][idet5-1].spill_time) < detSigWidth[5])) continue;

    // to start coincidence window conincides with the width of the gate
    Double_t CoincLow5  = detHits[5][idet5].spill_time;
    Double_t CoincHigh5 = CoincLow5 + detSigWidth[5]; 
   
    Int_t det4HitCoinc[100]={0};
    // find WC4 hits within this gate. returns number of hits and a vector with hit ID
    Double_t numDet4Coinc = FindDetCoinc(4, CoincLow5, CoincHigh5, det4HitCoinc);
    for(Int_t idet4 = 0; idet4 < numDet4Coinc; idet4++){
      
      // reduce the coincidence window based 
      Double_t CoincLow4 = std::max(CoincLow5,detHits[4][det4HitCoinc[idet4]].spill_time+detSigDelay[4]);
      Double_t CoincHigh4 = std::min(CoincHigh5,detHits[4][det4HitCoinc[idet4]].spill_time + detSigWidth[4]+detSigDelay[4]);
      Int_t det2HitCoinc[100]={0};      
      // find WC3 hits within possible coincidence windows. returns number of hits and a vector with hit ID
      Double_t numDet2Coinc = FindDetCoinc(2, CoincLow4, CoincHigh4, det2HitCoinc);
      trig245Num += numDet2Coinc;
           
        
        for(Int_t idet2 = 0; idet2 < numDet2Coinc; idet2++){
          // reduce the coincidence window based
          Double_t CoincLow2 = std::max(CoincLow4, detHits[2][det2HitCoinc[idet2]].spill_time+detSigDelay[2]);
          Double_t CoincHigh2 = std::min(CoincHigh4, detHits[2][det2HitCoinc[idet2]].spill_time + detSigWidth[2]+detSigDelay[2]);

          Int_t det1HitCoinc[100]={0};
          // find WC1 hits within possible coincidence windows. returns number of hits and a vector with hit ID
          Double_t numDet1Coinc = FindDetCoinc(1, CoincLow2, CoincHigh2, det1HitCoinc);
          trig1245Num += numDet1Coinc;
          
          for(Int_t idet1 = 0; idet1 < numDet1Coinc; idet1++){
            // reduce the coincidence window based
            Double_t CoincLow1 = std::max(CoincLow2, detHits[1][det1HitCoinc[idet1]].spill_time+detSigDelay[1]);
            Double_t CoincHigh1 = std::min(CoincHigh2, detHits[1][det1HitCoinc[idet1]].spill_time + detSigWidth[1]+detSigDelay[1]);

            Int_t det0HitCoinc[100]={0};
            // find USTOF hits within possible coincidence windows. returns number of hits and a vector with hit ID
            Double_t numDet0Coinc = FindDetCoinc(0, CoincLow1, CoincHigh1, det0HitCoinc);
            trig01245Num += numDet0Coinc;
           
            for(Int_t idet0 = 0; idet0 < numDet0Coinc; idet0++){
	      Double_t detHitTime0 = detHits[0][det0HitCoinc[idet0]].spill_time + detSigDelay[0]; 
              Double_t trigTimeLow = std::max(CoincLow1, detHitTime0);
              Double_t trigTimeHigh = std::min(CoincHigh1, detHitTime0 + detSigWidth[0]);

	      triggertimes.push_back(trigTimeLow);
	      TOFdshit.push_back(detHits[5][idet5].TrackID);
	    }
          }          
        }
    }
  }    	
    //  S2 Trigger condition
  for(Int_t idet5 = 0; idet5 < detNumHits[5]; idet5++){

    if(detHits[5][idet5].PDGid == 22) continue;    

    // open 100ns gate starting with DSTOF hit
    // skip consequent hits within the 100ns gate
    if(idet5>0 && ((detHits[5][idet5].spill_time - detHits[5][idet5-1].spill_time) < detSigWidth[5])) continue;

    // to start coincidence window conincides with the width of the gate
    Double_t CoincLow5  = detHits[5][idet5].spill_time;
    Double_t CoincHigh5 = CoincLow5 + detSigWidth[5]; 
   
    Int_t det4HitCoinc[100]={0};
    // find WC4 hits within this gate. returns number of hits and a vector with hit ID
    Double_t numDet4Coinc = FindDetCoinc(4, CoincLow5, CoincHigh5, det4HitCoinc);
    for(Int_t idet4 = 0; idet4 < numDet4Coinc; idet4++){
      
      // reduce the coincidence window based 
      Double_t CoincLow4 = std::max(CoincLow5, detHits[4][det4HitCoinc[idet4]].spill_time+detSigDelay[4]);
      Double_t CoincHigh4 = std::min(CoincHigh5, detHits[4][det4HitCoinc[idet4]].spill_time + detSigWidth[4]+detSigDelay[4]);

      Int_t det3HitCoinc[100]={0};      
      // find WC3 hits within possible coincidence windows. returns number of hits and a vector with hit ID
      Double_t numDet3Coinc = FindDetCoinc(3, CoincLow4, CoincHigh4, det3HitCoinc);

           
        
        for(Int_t idet3 = 0; idet3 < numDet3Coinc; idet3++){

          // reduce the coincidence window based
          Double_t CoincLow3 = std::max(CoincLow4, detHits[3][det3HitCoinc[idet3]].spill_time+detSigDelay[3]);
          Double_t CoincHigh3 = std::min(CoincHigh4, detHits[3][det3HitCoinc[idet3]].spill_time + detSigWidth[3]+detSigDelay[3]);

          Int_t det1HitCoinc[100]={0};
          // find WC1 hits within possible coincidence windows. returns number of hits and a vector with hit ID
          Double_t numDet1Coinc = FindDetCoinc(1, CoincLow3, CoincHigh3, det1HitCoinc);
          trig1345Num += numDet1Coinc;
          
          for(Int_t idet1 = 0; idet1 < numDet1Coinc; idet1++){

            // reduce the coincidence window based
            Double_t CoincLow1 = std::max(CoincLow3, detHits[1][det1HitCoinc[idet1]].spill_time+detSigDelay[1]);
            Double_t CoincHigh1 = std::min(CoincHigh3, detHits[1][det1HitCoinc[idet1]].spill_time + detSigWidth[1]+detSigDelay[1]);

            Int_t det0HitCoinc[100]={0};
            // find USTOF hits within possible coincidence windows. returns number of hits and a vector with hit ID
            Double_t numDet0Coinc = FindDetCoinc(0, CoincLow1, CoincHigh1, det0HitCoinc);
            trig01345Num += numDet0Coinc;
           
            for(Int_t idet0 = 0; idet0 < numDet0Coinc; idet0++){

	      Double_t detHitTime0 = detHits[0][det0HitCoinc[idet0]].spill_time + detSigDelay[0]; 
              Double_t trigTimeLow = std::max(CoincLow1, detHitTime0);
              Double_t trigTimeHigh = std::min(CoincHigh1, detHitTime0 + detSigWidth[0]);

	      triggertimes.push_back(trigTimeLow);
	      TOFdshit.push_back(detHits[5][idet5].TrackID);
	    }
          }          
      }   
    }
  }        
   
  std::cout<<"trig45Num = "<<trig45Num<<std::endl;
  std::cout<<"trig345Num = "<<trig345Num<<std::endl;
  std::cout<<"trig2345Num = "<<trig2345Num<<std::endl;
  std::cout<<"trig12345Num = "<<trig12345Num<<std::endl;
  std::cout<<"trig012345Num = "<<trig012345Num<<std::endl;
  std::cout<<"trig1345Num = "<<trig1345Num<<std::endl;
  std::cout<<"trig01345Num = "<<trig01345Num<<std::endl;
  std::cout<<"trig245Num = "<<trig245Num<<std::endl;
  std::cout<<"trig1245Num = "<<trig1245Num<<std::endl;
  std::cout<<"trig01245Num = "<<trig01245Num<<std::endl;

  //Sort the trigger times, and remove times that are within a WC readout buffer away.
  std::cout<<"Found "<<triggertimes.size()<<" triggers across the three possible detector combinations. Clearing out triggers that are within a WC readout buffer of other triggers."<<std::endl;
  std::sort(triggertimes.begin(),triggertimes.end());
  for (int i=triggertimes.size()-1; i>0; i--)
  {
    if ((triggertimes.at(i)-triggertimes.at(i-1))<1024*1.117E-9)
    {
      triggertimes.erase(triggertimes.begin()+i-1);
    }
  }
  std::cout<<"Reduced to "<<triggertimes.size()<<" triggers. Collecting particles around the trigger."<<std::endl;
  
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!                                                                         !! 
//!!                                                                         !!
//!!  This is where your output tree will go. Make sure this line is right!  !!
//!!                                                                         !!
//!!                                                                         !!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TFile *outfile=new TFile("outdirbaseTriggerTreeSpillSPILLNUM.root","recreate");


TTree *outtree= new TTree("TriggerTree","TriggerTree");
TBranch *b_evt = outtree->Branch("Event",&evt,"evt/I");
TBranch *b_xUSTOF = outtree->Branch("xUSTOF",xUSTOF, "xUSTOF[100]/D");
TBranch *b_yUSTOF = outtree->Branch("yUSTOF",yUSTOF, "yUSTOF[100]/D");
TBranch *b_zUSTOF = outtree->Branch("zUSTOF",zUSTOF, "zUSTOF[100]/D");
TBranch *b_tUSTOF = outtree->Branch("tUSTOF",tUSTOF, "tUSTOF[100]/D");
TBranch *b_rtUSTOF = outtree->Branch("rtUSTOF",rtUSTOF, "rtUSTOF[100]/D");
TBranch *b_PDGUSTOF = outtree->Branch("PDGUSTOF",&PDGUSTOF, "PDGUSTOF[100]/D");
TBranch *b_PxUSTOF = outtree->Branch("PxUSTOF",&PxUSTOF, "PxUSTOF[100]/D");
TBranch *b_PyUSTOF = outtree->Branch("PyUSTOF",&PyUSTOF, "PyUSTOF[100]/D");
TBranch *b_PzUSTOF = outtree->Branch("PzUSTOF",&PzUSTOF, "PzUSTOF[100]/D");
TBranch *b_ParentIDUSTOF = outtree->Branch("ParentIDUSTOF",&ParentIDUSTOF, "ParentIDUSTOF[100]/D");
TBranch *b_TrackIDUSTOF = outtree->Branch("TrackIDUSTOF",&TrackIDUSTOF, "TrackIDUSTOF[100]/D");
TBranch *b_EventIDUSTOF = outtree->Branch("EventIDUSTOF",&EventIDUSTOF, "EventIDUSTOF[100]/D");

TBranch *b_xWC1 = outtree->Branch("xWC1",xWC1, "xWC1[100]/D");
TBranch *b_yWC1 = outtree->Branch("yWC1",yWC1, "yWC1[100]/D");
TBranch *b_zWC1 = outtree->Branch("zWC1",zWC1, "zWC1[100]/D");
TBranch *b_tWC1 = outtree->Branch("tWC1",tWC1, "tWC1[100]/D");
TBranch *b_rtWC1 = outtree->Branch("rtWC1",rtWC1, "rtWC1[100]/D");
TBranch *b_PDGWC1 = outtree->Branch("PDGWC1",&PDGWC1, "PDGWC1[100]/D");
TBranch *b_PxWC1 = outtree->Branch("PxWC1",&PxWC1, "PxWC1[100]/D");
TBranch *b_PyWC1 = outtree->Branch("PyWC1",&PyWC1, "PyWC1[100]/D");
TBranch *b_PzWC1 = outtree->Branch("PzWC1",&PzWC1, "PzWC1[100]/D");
TBranch *b_ParentIDWC1 = outtree->Branch("ParentIDWC1",&ParentIDWC1, "ParentIDWC1[100]/D");
TBranch *b_TrackIDWC1 = outtree->Branch("TrackIDWC1",&TrackIDWC1, "TrackIDWC1[100]/D");
TBranch *b_EventIDWC1 = outtree->Branch("EventIDWC1",&EventIDWC1, "EventIDWC1[100]/D");


TBranch *b_xWC2 = outtree->Branch("xWC2",xWC2, "xWC2[100]/D");
TBranch *b_yWC2 = outtree->Branch("yWC2",yWC2, "yWC2[100]/D");
TBranch *b_zWC2 = outtree->Branch("zWC2",zWC2, "zWC2[100]/D");
TBranch *b_tWC2 = outtree->Branch("tWC2",tWC2, "tWC2[100]/D");
TBranch *b_rtWC2 = outtree->Branch("rtWC2",rtWC2, "rtWC2[100]/D");
TBranch *b_PDGWC2 = outtree->Branch("PDGWC2",&PDGWC2, "PDGWC2[100]/D");
TBranch *b_PxWC2 = outtree->Branch("PxWC2",&PxWC2, "PxWC2[100]/D");
TBranch *b_PyWC2 = outtree->Branch("PyWC2",&PyWC2, "PyWC2[100]/D");
TBranch *b_PzWC2 = outtree->Branch("PzWC2",&PzWC2, "PzWC2[100]/D");
TBranch *b_ParentIDWC2 = outtree->Branch("ParentIDWC2",&ParentIDWC2, "ParentIDWC2[100]/D");
TBranch *b_TrackIDWC2 = outtree->Branch("TrackIDWC2",&TrackIDWC2, "TrackIDWC2[100]/D");
TBranch *b_EventIDWC2 = outtree->Branch("EventIDWC2",&EventIDWC2, "EventIDWC2[100]/D");


TBranch *b_xWC3 = outtree->Branch("xWC3",xWC3, "xWC3[100]/D");
TBranch *b_yWC3 = outtree->Branch("yWC3",yWC3, "yWC3[100]/D");
TBranch *b_zWC3 = outtree->Branch("zWC3",zWC3, "zWC3[100]/D");
TBranch *b_tWC3 = outtree->Branch("tWC3",tWC3, "tWC3[100]/D");
TBranch *b_rtWC3 = outtree->Branch("rtWC3",rtWC3, "rtWC3[100]/D");
TBranch *b_PDGWC3 = outtree->Branch("PDGWC3",&PDGWC3, "PDGWC3[100]/D");
TBranch *b_PxWC3 = outtree->Branch("PxWC3",&PxWC3, "PxWC3[100]/D");
TBranch *b_PyWC3 = outtree->Branch("PyWC3",&PyWC3, "PyWC3[100]/D");
TBranch *b_PzWC3 = outtree->Branch("PzWC3",&PzWC3, "PzWC3[100]/D");
TBranch *b_ParentIDWC3 = outtree->Branch("ParentIDWC3",&ParentIDWC3, "ParentIDWC3[100]/D");
TBranch *b_TrackIDWC3 = outtree->Branch("TrackIDWC3",&TrackIDWC3, "TrackIDWC3[100]/D");
TBranch *b_EventIDWC3 = outtree->Branch("EventIDWC3",&EventIDWC3, "EventIDWC3[100]/D");

TBranch *b_xWC4 = outtree->Branch("xWC4",xWC4, "xWC4[100]/D");
TBranch *b_yWC4 = outtree->Branch("yWC4",yWC4, "yWC4[100]/D");
TBranch *b_zWC4 = outtree->Branch("zWC4",zWC4, "zWC4[100]/D");
TBranch *b_tWC4 = outtree->Branch("tWC4",tWC4, "tWC4[100]/D");
TBranch *b_rtWC4 = outtree->Branch("rtWC4",rtWC4, "rtWC4[100]/D");
TBranch *b_PDGWC4 = outtree->Branch("PDGWC4",&PDGWC4, "PDGWC4[100]/D");
TBranch *b_PxWC4 = outtree->Branch("PxWC4",&PxWC4, "PxWC4[100]/D");
TBranch *b_PyWC4 = outtree->Branch("PyWC4",&PyWC4, "PyWC4[100]/D");
TBranch *b_PzWC4 = outtree->Branch("PzWC4",&PzWC4, "PzWC4[100]/D");
TBranch *b_ParentIDWC4 = outtree->Branch("ParentIDWC4",&ParentIDWC4, "ParentIDWC4[100]/D");
TBranch *b_TrackIDWC4 = outtree->Branch("TrackIDWC4",&TrackIDWC4, "TrackIDWC4[100]/D");  
TBranch *b_EventIDWC4 = outtree->Branch("EventIDWC4",&EventIDWC4, "EventIDWC4[100]/D");

TBranch *b_xDSTOF = outtree->Branch("xDSTOF",xDSTOF, "xDSTOF[100]/D");
TBranch *b_yDSTOF = outtree->Branch("yDSTOF",yDSTOF, "yDSTOF[100]/D");
TBranch *b_zDSTOF = outtree->Branch("zDSTOF",zDSTOF, "zDSTOF[100]/D");
TBranch *b_tDSTOF = outtree->Branch("tDSTOF",tDSTOF, "tDSTOF[100]/D");
TBranch *b_rtDSTOF = outtree->Branch("rtDSTOF",rtDSTOF, "rtDSTOF[100]/D");
TBranch *b_PDGDSTOF = outtree->Branch("PDGDSTOF",&PDGDSTOF, "PDGDSTOF[100]/D");
TBranch *b_PxDSTOF = outtree->Branch("PxDSTOF",&PxDSTOF, "PxDSTOF[100]/D");
TBranch *b_PyDSTOF = outtree->Branch("PyDSTOF",&PyDSTOF, "PyDSTOF[100]/D");
TBranch *b_PzDSTOF = outtree->Branch("PzDSTOF",&PzDSTOF, "PzDSTOF[100]/D");
TBranch *b_ParentIDDSTOF = outtree->Branch("ParentIDDSTOF",&ParentIDDSTOF, "ParentIDDSTOF[100]/D");
TBranch *b_TrackIDDSTOF = outtree->Branch("TrackIDDSTOF",&TrackIDDSTOF, "TrackIDDSTOF[100]/D");
TBranch *b_EventIDDSTOF = outtree->Branch("EventIDDSTOF",&EventIDDSTOF, "EventIDDSTOF[100]/D");

TBranch *b_beamtimesave = outtree->Branch("beam_timesave",&beam_timesave, "beam_timesave[350000]/D");  
TBranch *b_n45= outtree->Branch("n45",&n45,"n45/D");
TBranch *b_n345= outtree->Branch("n345",&n345,"n345/D");
TBranch *b_n2345= outtree->Branch("n2345",&n2345,"n2345/D");
TBranch *b_n12345= outtree->Branch("n12345",&n12345,"n12345/D");
TBranch *b_n012345= outtree->Branch("n012345",&n012345,"n012345/D");
TBranch *b_n1345= outtree->Branch("n1345",&n1345,"n1345/D");
TBranch *b_n01345= outtree->Branch("n01345",&n01345,"n01345/D");
TBranch *b_n245= outtree->Branch("n245",&n245,"n245/D");
TBranch *b_n1245= outtree->Branch("n1245",&n1245,"n1245/D");
TBranch *b_n01245= outtree->Branch("n01245",&n01245,"n01245/D");
TBranch *b_Spill= outtree->Branch("Spill",&spill,"spill/D");


  //Now that we've removed triggers that wouldnt actually happen because of dead time in the trigger readout, loop through the tree and save information for particles during the readout buffer for that trigger.
  for (int i=0; i<(int)triggertimes.size();i++)
  {
    ResetTriggerTree();
    for(int jentry=0; jentry<detNumHits[0];jentry++)
    {
      if ((detHits[0][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[0][jentry].spill_time-triggertimes.at(i)<WCbufferend) && detHits[0][jentry].PDGid!=22)
      {
	xUSTOF[USTOFiter]=detHits[0][jentry].x;
	yUSTOF[USTOFiter]=detHits[0][jentry].y;
	zUSTOF[USTOFiter]=detHits[0][jentry].z;
	tUSTOF[USTOFiter]=detHits[0][jentry].spill_time;
	rtUSTOF[USTOFiter]=detHits[0][jentry].t;
	PDGUSTOF[USTOFiter]=detHits[0][jentry].PDGid;
	PxUSTOF[USTOFiter]=detHits[0][jentry].Px; 
	PyUSTOF[USTOFiter]=detHits[0][jentry].Py; 
	PzUSTOF[USTOFiter]=detHits[0][jentry].Pz;
	ParentIDUSTOF[USTOFiter]=detHits[0][jentry].ParentID;
	TrackIDUSTOF[USTOFiter]=detHits[0][jentry].TrackID;
	EventIDUSTOF[USTOFiter]=detHits[0][jentry].EventID;
	USTOFiter++;
      }
    }

    for(int jentry=0; jentry<detNumHits[1];jentry++)
    {
      if ((detHits[1][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[1][jentry].spill_time-triggertimes.at(i)<WCbufferend) && detHits[1][jentry].PDGid!=22)
      {
	xWC1[wc1iter]=detHits[1][jentry].x;
	yWC1[wc1iter]=detHits[1][jentry].y;
	zWC1[wc1iter]=detHits[1][jentry].z;
	tWC1[wc1iter]=detHits[1][jentry].spill_time;
	rtWC1[wc1iter]=detHits[1][jentry].t;
	PDGWC1[wc1iter]=detHits[1][jentry].PDGid;
	PxWC1[wc1iter]=detHits[1][jentry].Px; 
	PyWC1[wc1iter]=detHits[1][jentry].Py; 
	PzWC1[wc1iter]=detHits[1][jentry].Pz;
	ParentIDWC1[wc1iter]=detHits[1][jentry].ParentID;
	TrackIDWC1[wc1iter]=detHits[1][jentry].TrackID;
	EventIDWC1[wc1iter]=detHits[1][jentry].EventID;
	wc1iter++;
      }
    }
    for(int jentry=0; jentry<detNumHits[2];jentry++)
    {    
      if ((detHits[2][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[2][jentry].spill_time-triggertimes.at(i)<WCbufferend) && detHits[2][jentry].PDGid!=22)
      {
	xWC2[wc2iter]=detHits[2][jentry].x;
	yWC2[wc2iter]=detHits[2][jentry].y;
	zWC2[wc2iter]=detHits[2][jentry].z;
	tWC2[wc2iter]=detHits[2][jentry].spill_time;
        rtWC2[wc2iter]=detHits[2][jentry].t;
	PDGWC2[wc2iter]=detHits[2][jentry].PDGid;
	PxWC2[wc2iter]=detHits[2][jentry].Px; 
	PyWC2[wc2iter]=detHits[2][jentry].Py; 
	PzWC2[wc2iter]=detHits[2][jentry].Pz;
	ParentIDWC2[wc2iter]=detHits[2][jentry].ParentID;
	TrackIDWC2[wc2iter]=detHits[2][jentry].TrackID;
	EventIDWC2[wc2iter]=detHits[2][jentry].EventID;
	wc2iter++;
      }
    }
    for(int jentry=0; jentry<detNumHits[3];jentry++)
    {      
      if ((detHits[3][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[3][jentry].spill_time-triggertimes.at(i)<WCbufferend)&& detHits[3][jentry].PDGid!=22)
      {
	xWC3[wc3iter]=detHits[3][jentry].x;
	yWC3[wc3iter]=detHits[3][jentry].y;
	zWC3[wc3iter]=detHits[3][jentry].z;
	tWC3[wc3iter]=detHits[3][jentry].spill_time;
        rtWC3[wc3iter]=detHits[3][jentry].t;
	PDGWC3[wc3iter]=detHits[3][jentry].PDGid;
	PxWC3[wc3iter]=detHits[3][jentry].Px; 
	PyWC3[wc3iter]=detHits[3][jentry].Py; 
	PzWC3[wc3iter]=detHits[3][jentry].Pz;
	ParentIDWC3[wc3iter]=detHits[3][jentry].ParentID;
	TrackIDWC3[wc3iter]=detHits[3][jentry].TrackID;
	EventIDWC3[wc3iter]=detHits[3][jentry].EventID;
	wc3iter++;
      }
    }  
    for(int jentry=0; jentry<detNumHits[4];jentry++)
    {      
      if ((detHits[4][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[4][jentry].spill_time-triggertimes.at(i)<WCbufferend)&& detHits[4][jentry].PDGid!=22)
      {
	xWC4[wc4iter]=detHits[4][jentry].x;
	yWC4[wc4iter]=detHits[4][jentry].y;
	zWC4[wc4iter]=detHits[4][jentry].z;
	tWC4[wc4iter]=detHits[4][jentry].spill_time;
	rtWC4[wc4iter]=detHits[4][jentry].t;
	PDGWC4[wc4iter]=detHits[4][jentry].PDGid;
	PxWC4[wc4iter]=detHits[4][jentry].Px; 
	PyWC4[wc4iter]=detHits[4][jentry].Py; 
	PzWC4[wc4iter]=detHits[4][jentry].Pz;
	ParentIDWC4[wc4iter]=detHits[4][jentry].ParentID;
	TrackIDWC4[wc4iter]=detHits[4][jentry].TrackID;
	EventIDWC4[wc4iter]=detHits[4][jentry].EventID;
	wc4iter++;
      }                  
    }
    
    for(int jentry=0; jentry<detNumHits[5];jentry++)
    {
      if ((detHits[5][jentry].spill_time-triggertimes.at(i)>WCbufferstart && detHits[5][jentry].spill_time-triggertimes.at(i)<WCbufferend) && detHits[5][jentry].PDGid!=22)
      {
	xDSTOF[DSTOFiter]=detHits[5][jentry].x;
	yDSTOF[DSTOFiter]=detHits[5][jentry].y;
	zDSTOF[DSTOFiter]=detHits[5][jentry].z;
	tDSTOF[DSTOFiter]=detHits[5][jentry].spill_time;
	rtDSTOF[DSTOFiter]=detHits[5][jentry].t;
	PDGDSTOF[DSTOFiter]=detHits[5][jentry].PDGid;
	PxDSTOF[DSTOFiter]=detHits[5][jentry].Px; 
	PyDSTOF[DSTOFiter]=detHits[5][jentry].Py; 
	PzDSTOF[DSTOFiter]=detHits[5][jentry].Pz;
	ParentIDDSTOF[DSTOFiter]=detHits[5][jentry].ParentID;
	TrackIDDSTOF[DSTOFiter]=detHits[5][jentry].TrackID;
	EventIDDSTOF[DSTOFiter]=detHits[5][jentry].EventID;
	DSTOFiter++;
      }
    }    
    if(i==0)
    {
        n45=trig45Num;
    n345=trig345Num;
    n2345=trig2345Num;
    n12345=trig12345Num;
    n012345=trig012345Num;
    n1345=trig1345Num;
    n01345=trig01345Num;
    n245=trig245Num;
    n1245=trig1245Num;
    n01245=trig01245Num;
    spill=RawSPILLNUM;
      for(int ihit=0;ihit<350000;++ihit)
      {
        beam_timesave[ihit]=-9999;
      }
    }
    //Fill the tree
    outtree->Fill();
    evt++; 
  }
 //Write and Close
 outtree->Write();
 outfile->Close(); 
 delete outfile;   
}


void lariat_spillSPILLNUM(){

/// Initialize beam variable
//Pushback 6 times, one for each detector

Int_t NPrimary  = 350000;
Int_t BatchPOr  = 7;
Int_t BucketPBa = 84;
Int_t BucketPOr = BatchPOr * BucketPBa;

Double_t BucketCentSpace = 18.8e-9; // in seconds
Double_t BucketWidth     = 2.2e-9;  // in seconds
Double_t BatchLen        = BucketCentSpace * (Double_t)BucketPBa;
Double_t OrbitLen        = BatchLen * (Double_t)BatchPOr;
Double_t SpillLen        = 4.2; //in seconds
Double_t OrbitsPSpill    = SpillLen / OrbitLen;

// Check if these make sense
Double_t SpillTime =   (BucketPBa - 1) * BucketCentSpace 
                     + (BatchPOr - 1) * BatchLen 
                     + (OrbitsPSpill - 1) * OrbitLen
                     + BucketWidth;

std::cout<<"SpillTime = "<<SpillTime<<std::endl; 

//  Calculate a random time within the spill for each event 
for(Int_t iran = 0; iran < NPrimary; iran++){
  TRandom *r3 = new TRandom3(0);
  Int_t BucketInBatch = r3->Integer(BucketPBa);
  Int_t BatchInOrbit  = r3->Integer(BatchPOr - 1);
  Int_t OrbitInSpill  = r3->Integer((Int_t) OrbitsPSpill);
  Double_t PlaceInBucket = r3->Gaus(0,BucketWidth/3);


  beam_time[iran]     =       BucketWidth/2 + PlaceInBucket
                            + BucketInBatch * BucketCentSpace 
                            + BatchInOrbit * BatchLen   
                            + OrbitInSpill * OrbitLen;
  beam_timesave[iran]=beam_time[iran];	
}

// Sort the times so that smaller event numbers have earlier times.
std:sort(beam_time,beam_time+NPrimary);

std::cout<<"BucketPOr = "<<BucketPOr<<", BatchLen = "<<BatchLen<<", OrbitLen  = "<<OrbitLen<<", OrbitsPSpill = "<<OrbitsPSpill<<std::endl;

// open the spill file
TFile *larg4beam = new TFile("indirbase");
larg4beam->cd("VirtualDetector");

// read the detector trees and calculate the spill time for each hit
// and time-order the hits.

// USTOF hits
TTree *TOFus = (TTree*)gDirectory->Get("TOFus");
ReadTree(0, TOFus);
std::sort(detHits[0], detHits[0]+detNumHits[0], timecomp);

// WC1 hits
TTree *WC1 = (TTree*)gDirectory->Get("Det1");
ReadTree(1, WC1);
std::sort(detHits[1], detHits[1]+detNumHits[1], timecomp);

// WC2 hits
TTree *WC2 = (TTree*)gDirectory->Get("Det2");
ReadTree(2, WC2);
std::sort(detHits[2], detHits[2]+detNumHits[2], timecomp);

// WC3 hits
TTree *WC3 = (TTree*)gDirectory->Get("Det3");
ReadTree(3, WC3);
std::sort(detHits[3], detHits[3]+detNumHits[3], timecomp);

// WC4 hits
TTree *WC4 = (TTree*)gDirectory->Get("Det4");
ReadTree(4, WC4);
std::sort(detHits[4], detHits[4]+detNumHits[4], timecomp);

// DSTOF hits
TTree *TOFds = (TTree*)gDirectory->Get("TOFds");
ReadTree(5, TOFds);
std::sort(detHits[5], detHits[5]+detNumHits[5], timecomp);

larg4beam->Close();



// Find hit coincidences which can form triggers
FindTriggers();



}



