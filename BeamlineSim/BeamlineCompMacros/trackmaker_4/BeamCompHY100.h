//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  6 11:30:31 2018 by ROOT version 6.08/06
// from TTree TriggerTree/TriggerTree
// found on file: TriggerTreeSpillSPILLNUM00.root
//////////////////////////////////////////////////////////

#ifndef BeamCompHY100SpillSPILLNUM_h
#define BeamCompHY100SpillSPILLNUM_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BeamCompHY100SpillSPILLNUM	 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Double_t        xWC1[100];
   Double_t        yWC1[100];
   Double_t        zWC1[100];
   Double_t        tWC1[100];
   Double_t        rtWC1[100];
   Double_t        PDGWC1[100];
   Double_t        PxWC1[100];
   Double_t        PyWC1[100];
   Double_t        PzWC1[100];
   Double_t        ParentIDWC1[100];
   Double_t        TrackIDWC1[100];
   Double_t        EventIDWC1[100];
   Double_t        xWC2[100];
   Double_t        yWC2[100];
   Double_t        zWC2[100];
   Double_t        tWC2[100];
   Double_t        rtWC2[100];
   Double_t        PDGWC2[100];
   Double_t        PxWC2[100];
   Double_t        PyWC2[100];
   Double_t        PzWC2[100];
   Double_t        ParentIDWC2[100];
   Double_t        TrackIDWC2[100];
   Double_t        EventIDWC2[100];
   Double_t        xWC3[100];
   Double_t        yWC3[100];
   Double_t        zWC3[100];
   Double_t        tWC3[100];
   Double_t        rtWC3[100];
   Double_t        PDGWC3[100];
   Double_t        PxWC3[100];
   Double_t        PyWC3[100];
   Double_t        PzWC3[100];
   Double_t        ParentIDWC3[100];
   Double_t        TrackIDWC3[100];
   Double_t        EventIDWC3[100];
   Double_t        xWC4[100];
   Double_t        yWC4[100];
   Double_t        zWC4[100];
   Double_t        tWC4[100];
   Double_t        rtWC4[100];
   Double_t        PDGWC4[100];
   Double_t        PxWC4[100];
   Double_t        PyWC4[100];
   Double_t        PzWC4[100];
   Double_t        ParentIDWC4[100];
   Double_t        TrackIDWC4[100];
   Double_t        EventIDWC4[100];
   
   Double_t        xUSTOF[100];
   Double_t        yUSTOF[100];
   Double_t        zUSTOF[100];
   Double_t        tUSTOF[100];
   Double_t        rtUSTOF[100];
   Double_t        PDGUSTOF[100];
   Double_t        PxUSTOF[100];
   Double_t        PyUSTOF[100];
   Double_t        PzUSTOF[100];
   Double_t        ParentIDUSTOF[100];
   Double_t        TrackIDUSTOF[100];
   Double_t        EventIDUSTOF[100];
   
   Double_t        xDSTOF[100];
   Double_t        yDSTOF[100];
   Double_t        zDSTOF[100];
   Double_t        tDSTOF[100];
   Double_t        rtDSTOF[100];
   Double_t        PDGDSTOF[100];
   Double_t        PxDSTOF[100];
   Double_t        PyDSTOF[100];
   Double_t        PzDSTOF[100];
   Double_t        ParentIDDSTOF[100];
   Double_t        TrackIDDSTOF[100];
   Double_t        EventIDDSTOF[100];
   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_xWC1;   //!
   TBranch        *b_yWC1;   //!
   TBranch        *b_zWC1;   //!
   TBranch        *b_tWC1;   //!
   TBranch        *b_rtWC1;   //!
   TBranch        *b_PDGWC1;   //!
   TBranch        *b_PxWC1;   //!
   TBranch        *b_PyWC1;   //!
   TBranch        *b_PzWC1;   //!
   TBranch        *b_ParentIDWC1;   //!
   TBranch        *b_TrackIDWC1;   //!
   TBranch        *b_EventIDWC1;   //!
   TBranch        *b_xWC2;   //!
   TBranch        *b_yWC2;   //!
   TBranch        *b_zWC2;   //!
   TBranch        *b_tWC2;   //!
   TBranch        *b_rtWC2;   //!
   TBranch        *b_PDGWC2;   //!
   TBranch        *b_PxWC2;   //!
   TBranch        *b_PyWC2;   //!
   TBranch        *b_PzWC2;   //!
   TBranch        *b_ParentIDWC2;   //!
   TBranch        *b_TrackIDWC2;   //!
   TBranch        *b_EventIDWC2;   //!
   TBranch        *b_xWC3;   //!
   TBranch        *b_yWC3;   //!
   TBranch        *b_zWC3;   //!
   TBranch        *b_tWC3;   //!
   TBranch        *b_rtWC3;   //!
   TBranch        *b_PDGWC3;   //!
   TBranch        *b_PxWC3;   //!
   TBranch        *b_PyWC3;   //!
   TBranch        *b_PzWC3;   //!
   TBranch        *b_ParentIDWC3;   //!
   TBranch        *b_TrackIDWC3;   //!
   TBranch        *b_EventIDWC3;   //!
   TBranch        *b_xWC4;   //!
   TBranch        *b_yWC4;   //!
   TBranch        *b_zWC4;   //!
   TBranch        *b_tWC4;   //!
   TBranch        *b_rtWC4;   //!
   TBranch        *b_PDGWC4;   //!
   TBranch        *b_PxWC4;   //!
   TBranch        *b_PyWC4;   //!
   TBranch        *b_PzWC4;   //!
   TBranch        *b_ParentIDWC4;   //!
   TBranch        *b_TrackIDWC4;   //!
   TBranch        *b_EventIDWC4;   //!
   TBranch        *b_xUSTOF;   //!
   TBranch        *b_yUSTOF;   //!
   TBranch        *b_zUSTOF;   //!
   TBranch        *b_tUSTOF;   //!
   TBranch        *b_rtUSTOF;   //!
   TBranch        *b_PDGUSTOF;   //!
   TBranch        *b_PxUSTOF;   //!
   TBranch        *b_PyUSTOF;   //!
   TBranch        *b_PzUSTOF;   //!
   TBranch        *b_ParentIDUSTOF;   //!
   TBranch        *b_TrackIDUSTOF;   //!
   TBranch        *b_EventIDUSTOF;   //!
   TBranch        *b_xDSTOF;   //!
   TBranch        *b_yDSTOF;   //!
   TBranch        *b_zDSTOF;   //!
   TBranch        *b_tDSTOF;   //!
   TBranch        *b_rtDSTOF;   //!
   TBranch        *b_PDGDSTOF;   //!
   TBranch        *b_PxDSTOF;   //!
   TBranch        *b_PyDSTOF;   //!
   TBranch        *b_PzDSTOF;   //!
   TBranch        *b_ParentIDDSTOF;   //!
   TBranch        *b_TrackIDDSTOF;   //!
   TBranch        *b_EventIDDSTOF;   //!
   BeamCompHY100SpillSPILLNUM(TTree *tree=0);
   virtual ~BeamCompHY100SpillSPILLNUM();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   double                BuildFourPointTrack(std::vector<TVector3> wc1hits, std::vector<TVector3> wc2hits, std::vector<TVector3> wc3hits, std::vector<TVector3> wc4hits);
   double                FindTOF(std::vector<double> ushits, std::vector<double> dshits);
   void                  BuildThreePointTrack(std::vector<TVector3> wc1hits, std::vector<TVector3> wc2hits, std::vector<TVector3> wc3hits, std::vector<TVector3> wc4hits);
   std::vector<float>    Regression(double (&y)[4], double (&z)[4], int & WCMissed);
   void                  ResetVars();
   double                CalculateFourPointMomentum(int (&wcindex)[4][2], float yslope);
   double                CalculateThreePointMomentum(int (&wcindex)[4][2], float yslope, int WCMissed);
   bool                  PassWCFilter(int (&wcindex)[4][2]);
   std::vector<double>   ProjToZ(std::vector<double> hit0, std::vector<double> hit1, double zpos);
   bool                  CheckUpstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2);
   bool                  CheckDownstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2);
   bool                  CheckDownstreamCollimatorAperture(std::vector<double> hit1, std::vector<double> hit2);
   bool                  PassWCFilter(int (&wcindex)[4][2], int WCMissed, std::vector<float> trackstats, double momentum);
   bool                  MPtoWC4(std::vector<double> hit1, std::vector<double> hit2, std::vector<double> hit3, std::vector<double> hit4);
   bool                  MidplaneMatch(std::vector<double> hit1, std::vector<double> hit2, std::vector<double> hit3, std::vector<double> hit4);
   bool                  PureFrankTrackStatus(int (&ID)[4][2]);
   bool                  PureFrankTrackStatus(int (&ID)[4][2], int WCMissed);
};

#endif

#ifdef BeamCompHY100SpillSPILLNUM_cxx
BeamCompHY100SpillSPILLNUM::BeamCompHY100SpillSPILLNUM(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("indirbaseSpillSPILLNUM");
      if (!f || !f->IsOpen()) {
         f = new TFile("indirbaseSpillSPILLNUM");
      }
      f->GetObject("TriggerTree",tree);

   }
   Init(tree);
   Loop();
}

BeamCompHY100SpillSPILLNUM::~BeamCompHY100SpillSPILLNUM()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BeamCompHY100SpillSPILLNUM::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BeamCompHY100SpillSPILLNUM::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BeamCompHY100SpillSPILLNUM::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_evt);
   fChain->SetBranchAddress("xWC1", xWC1, &b_xWC1);
   fChain->SetBranchAddress("yWC1", yWC1, &b_yWC1);
   fChain->SetBranchAddress("zWC1", zWC1, &b_zWC1);
   fChain->SetBranchAddress("tWC1", tWC1, &b_tWC1);
   fChain->SetBranchAddress("rtWC1", rtWC1, &b_rtWC1);
   fChain->SetBranchAddress("PDGWC1", PDGWC1, &b_PDGWC1);
   fChain->SetBranchAddress("PxWC1", PxWC1, &b_PxWC1);
   fChain->SetBranchAddress("PyWC1", PyWC1, &b_PyWC1);
   fChain->SetBranchAddress("PzWC1", PzWC1, &b_PzWC1);
   fChain->SetBranchAddress("ParentIDWC1", ParentIDWC1, &b_ParentIDWC1);
   fChain->SetBranchAddress("TrackIDWC1", TrackIDWC1, &b_TrackIDWC1);
   fChain->SetBranchAddress("EventIDWC1", EventIDWC1, &b_EventIDWC1);
   fChain->SetBranchAddress("xWC2", xWC2, &b_xWC2);
   fChain->SetBranchAddress("yWC2", yWC2, &b_yWC2);
   fChain->SetBranchAddress("zWC2", zWC2, &b_zWC2);
   fChain->SetBranchAddress("tWC2", tWC2, &b_tWC2);
   fChain->SetBranchAddress("rtWC2", rtWC2, &b_rtWC2);
   fChain->SetBranchAddress("PDGWC2", PDGWC2, &b_PDGWC2);
   fChain->SetBranchAddress("PxWC2", PxWC2, &b_PxWC2);
   fChain->SetBranchAddress("PyWC2", PyWC2, &b_PyWC2);
   fChain->SetBranchAddress("PzWC2", PzWC2, &b_PzWC2);
   fChain->SetBranchAddress("ParentIDWC2", ParentIDWC2, &b_ParentIDWC2);
   fChain->SetBranchAddress("TrackIDWC2", TrackIDWC2, &b_TrackIDWC2);
   fChain->SetBranchAddress("EventIDWC2", EventIDWC2, &b_EventIDWC2);
   fChain->SetBranchAddress("xWC3", xWC3, &b_xWC3);
   fChain->SetBranchAddress("yWC3", yWC3, &b_yWC3);
   fChain->SetBranchAddress("zWC3", zWC3, &b_zWC3);
   fChain->SetBranchAddress("tWC3", tWC3, &b_tWC3);
   fChain->SetBranchAddress("rtWC3", rtWC3, &b_rtWC3);
   fChain->SetBranchAddress("PDGWC3", PDGWC3, &b_PDGWC3);
   fChain->SetBranchAddress("PxWC3", PxWC3, &b_PxWC3);
   fChain->SetBranchAddress("PyWC3", PyWC3, &b_PyWC3);
   fChain->SetBranchAddress("PzWC3", PzWC3, &b_PzWC3);
   fChain->SetBranchAddress("ParentIDWC3", ParentIDWC3, &b_ParentIDWC3);
   fChain->SetBranchAddress("TrackIDWC3", TrackIDWC3, &b_TrackIDWC3);
   fChain->SetBranchAddress("EventIDWC3", EventIDWC3, &b_EventIDWC3);
   fChain->SetBranchAddress("xWC4", xWC4, &b_xWC4);
   fChain->SetBranchAddress("yWC4", yWC4, &b_yWC4);
   fChain->SetBranchAddress("zWC4", zWC4, &b_zWC4);
   fChain->SetBranchAddress("tWC4", tWC4, &b_tWC4);
   fChain->SetBranchAddress("rtWC4", rtWC4, &b_rtWC4);
   fChain->SetBranchAddress("PDGWC4", PDGWC4, &b_PDGWC4);
   fChain->SetBranchAddress("PxWC4", PxWC4, &b_PxWC4);
   fChain->SetBranchAddress("PyWC4", PyWC4, &b_PyWC4);
   fChain->SetBranchAddress("PzWC4", PzWC4, &b_PzWC4);
   fChain->SetBranchAddress("ParentIDWC4", ParentIDWC4, &b_ParentIDWC4);
   fChain->SetBranchAddress("TrackIDWC4", TrackIDWC4, &b_TrackIDWC4);
   fChain->SetBranchAddress("EventIDWC4", EventIDWC4, &b_EventIDWC4);
   fChain->SetBranchAddress("xUSTOF", xUSTOF, &b_xUSTOF);
   fChain->SetBranchAddress("yUSTOF", yUSTOF, &b_yUSTOF);
   fChain->SetBranchAddress("zUSTOF", zUSTOF, &b_zUSTOF);
   fChain->SetBranchAddress("tUSTOF", tUSTOF, &b_tUSTOF);
   fChain->SetBranchAddress("rtUSTOF", rtUSTOF, &b_rtUSTOF);
   fChain->SetBranchAddress("PDGUSTOF", PDGUSTOF, &b_PDGUSTOF);
   fChain->SetBranchAddress("PxUSTOF", PxUSTOF, &b_PxUSTOF);
   fChain->SetBranchAddress("PyUSTOF", PyUSTOF, &b_PyUSTOF);
   fChain->SetBranchAddress("PzUSTOF", PzUSTOF, &b_PzUSTOF);
   fChain->SetBranchAddress("ParentIDUSTOF", ParentIDUSTOF, &b_ParentIDUSTOF);
   fChain->SetBranchAddress("TrackIDUSTOF", TrackIDUSTOF, &b_TrackIDUSTOF);
   fChain->SetBranchAddress("EventIDUSTOF", EventIDUSTOF, &b_EventIDUSTOF);
      fChain->SetBranchAddress("xDSTOF", xDSTOF, &b_xDSTOF);
   fChain->SetBranchAddress("yDSTOF", yDSTOF, &b_yDSTOF);
   fChain->SetBranchAddress("zDSTOF", zDSTOF, &b_zDSTOF);
   fChain->SetBranchAddress("tDSTOF", tDSTOF, &b_tDSTOF);
   fChain->SetBranchAddress("rtDSTOF", rtDSTOF, &b_rtDSTOF);
   fChain->SetBranchAddress("PDGDSTOF", PDGDSTOF, &b_PDGDSTOF);
   fChain->SetBranchAddress("PxDSTOF", PxDSTOF, &b_PxDSTOF);
   fChain->SetBranchAddress("PyDSTOF", PyDSTOF, &b_PyDSTOF);
   fChain->SetBranchAddress("PzDSTOF", PzDSTOF, &b_PzDSTOF);
   fChain->SetBranchAddress("ParentIDDSTOF", ParentIDDSTOF, &b_ParentIDDSTOF);
   fChain->SetBranchAddress("TrackIDDSTOF", TrackIDDSTOF, &b_TrackIDDSTOF);
   fChain->SetBranchAddress("EventIDDSTOF", EventIDDSTOF, &b_EventIDDSTOF);

   Notify();
}

Bool_t BeamCompHY100SpillSPILLNUM::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BeamCompHY100SpillSPILLNUM::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BeamCompHY100SpillSPILLNUM::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BeamCompHY100SpillSPILLNUM_cxx
