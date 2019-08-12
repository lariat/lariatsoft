//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 11 15:26:01 2019 by ROOT version 6.10/04
// from TTree TrackTree/TrackTree
// found on file: CombinedNeg100Tracks.root
//////////////////////////////////////////////////////////

#ifndef BeamComp64GeVNeg100ARound2_h
#define BeamComp64GeVNeg100ARound2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BeamComp64GeVNeg100ARound2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        SpillNumber;
   Double_t        Momentum;
   Double_t        MomentumError[2];
   Double_t        TrueHitMom[4][2];
   Double_t        trackIDs[4][2];
   Double_t        MotherIDs[4][2];
   Double_t        trackpdgs[4][2];
   Bool_t          FourptTrack;
   Bool_t          WCFilt;
   Bool_t          S2Track;
   Bool_t          S3Track;
   Bool_t          PickyTrack;
   Bool_t          HYTrack;
   Bool_t          PureTrack;
   Bool_t          Frankentrack;
   Double_t        thetaus;
   Double_t        thetads;
   Double_t        yslope;
   Double_t        MinRes;
   Double_t        MidDiff;
   Double_t        MidDiffX;
   Double_t        MidDiffY;
   Double_t        WC4Diff;
   Double_t        WC4DiffX;
   Double_t        WC4DiffY;
   Double_t        WC1Mult;
   Double_t        WC2Mult;
   Double_t        WC3Mult;
   Double_t        WC4Mult;
   Double_t        WC1Miss;
   Double_t        XYZTrack[4][3];
   Double_t        TOF;
   Double_t        Mass;

   // List of branches
   TBranch        *b_SpillNumber;   //!
   TBranch        *b_momentum;   //!
   TBranch        *b_error;   //!
   TBranch        *b_truemom;   //!
   TBranch        *b_trackIDs;   //!
   TBranch        *b_MotherIDs;   //!
   TBranch        *b_trackpdgs;   //!
   TBranch        *b_is4Pt;   //!
   TBranch        *b_wcfilt;   //!
   TBranch        *b_isS2;   //!
   TBranch        *b_isS3;   //!
   TBranch        *b_isPicky;   //!
   TBranch        *b_isHY;   //!
   TBranch        *b_isPure;   //!
   TBranch        *b_isFrankentrack;   //!
   TBranch        *b_theta_x_us;   //!
   TBranch        *b_theta_x_ds;   //!
   TBranch        *b_yslope;   //!
   TBranch        *b_MinRes;   //!
   TBranch        *b_MidDiff;   //!
   TBranch        *b_MidDiffX;   //!
   TBranch        *b_MidDiffY;   //!
   TBranch        *b_WC4Diff;   //!
   TBranch        *b_WC4DiffX;   //!
   TBranch        *b_WC4DiffY;   //!
   TBranch        *b_WC1Mult;   //!
   TBranch        *b_WC2Mult;   //!
   TBranch        *b_WC3Mult;   //!
   TBranch        *b_WC4Mult;   //!
   TBranch        *b_WC1Miss;   //!
   TBranch        *b_XYZTrack;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_mass;   //!

   BeamComp64GeVNeg100ARound2(TTree *tree=0);
   virtual ~BeamComp64GeVNeg100ARound2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   TPaveText* MakeTextBox(float x, float y, float textSize, float numLines, float width);
   double FindPDG(double(&pdg)[4][2], double (&mother)[4][2], double (&id)[4][2]);
};

#endif

#ifdef BeamComp64GeVNeg100ARound2_cxx
BeamComp64GeVNeg100ARound2::BeamComp64GeVNeg100ARound2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tracks.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Tracks.root");
      }
      f->GetObject("TrackTree",tree);

   }
   Init(tree);
   Loop();
}

BeamComp64GeVNeg100ARound2::~BeamComp64GeVNeg100ARound2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BeamComp64GeVNeg100ARound2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BeamComp64GeVNeg100ARound2::LoadTree(Long64_t entry)
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

void BeamComp64GeVNeg100ARound2::Init(TTree *tree)
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

   fChain->SetBranchAddress("SpillNumber", &SpillNumber, &b_SpillNumber);
   fChain->SetBranchAddress("Momentum", &Momentum, &b_momentum);
   fChain->SetBranchAddress("MomentumError", &MomentumError, &b_error);
   fChain->SetBranchAddress("TrueHitMom", TrueHitMom, &b_truemom);
   fChain->SetBranchAddress("trackIDs", trackIDs, &b_trackIDs);
   fChain->SetBranchAddress("MotherIDs", MotherIDs, &b_MotherIDs);
   fChain->SetBranchAddress("trackpdgs", trackpdgs, &b_trackpdgs);
   fChain->SetBranchAddress("FourptTrack", &FourptTrack, &b_is4Pt);
   fChain->SetBranchAddress("WCFilt", &WCFilt, &b_wcfilt);
   fChain->SetBranchAddress("S2Track", &S2Track, &b_isS2);
   fChain->SetBranchAddress("S3Track", &S3Track, &b_isS3);
   fChain->SetBranchAddress("PickyTrack", &PickyTrack, &b_isPicky);
   fChain->SetBranchAddress("HYTrack", &HYTrack, &b_isHY);
   fChain->SetBranchAddress("PureTrack", &PureTrack, &b_isPure);
   fChain->SetBranchAddress("Frankentrack", &Frankentrack, &b_isFrankentrack);
   fChain->SetBranchAddress("thetaus", &thetaus, &b_theta_x_us);
   fChain->SetBranchAddress("thetads", &thetads, &b_theta_x_ds);
   fChain->SetBranchAddress("yslope", &yslope, &b_yslope);
   fChain->SetBranchAddress("MinRes", &MinRes, &b_MinRes);
   fChain->SetBranchAddress("MidDiff", &MidDiff, &b_MidDiff);
   fChain->SetBranchAddress("MidDiffX", &MidDiffX, &b_MidDiffX);
   fChain->SetBranchAddress("MidDiffY", &MidDiffY, &b_MidDiffY);
   fChain->SetBranchAddress("WC4Diff", &WC4Diff, &b_WC4Diff);
   fChain->SetBranchAddress("WC4DiffX", &WC4DiffX, &b_WC4DiffX);
   fChain->SetBranchAddress("WC4DiffY", &WC4DiffY, &b_WC4DiffY);
   fChain->SetBranchAddress("WC1Mult", &WC1Mult, &b_WC1Mult);
   fChain->SetBranchAddress("WC2Mult", &WC2Mult, &b_WC2Mult);
   fChain->SetBranchAddress("WC3Mult", &WC3Mult, &b_WC3Mult);
   fChain->SetBranchAddress("WC4Mult", &WC4Mult, &b_WC4Mult);
   fChain->SetBranchAddress("WC1Miss", &WC1Miss, &b_WC1Miss);
   fChain->SetBranchAddress("XYZTrack", XYZTrack, &b_XYZTrack);
   fChain->SetBranchAddress("TOF", &TOF, &b_tof);
   fChain->SetBranchAddress("Mass", &Mass, &b_mass);
   Notify();
}

Bool_t BeamComp64GeVNeg100ARound2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
TPaveText* BeamComp64GeVNeg100ARound2::MakeTextBox(float x, float y, float textSize, float numLines, float width)
{
  TPaveText *pt = new TPaveText(x, y-numLines*textSize, x+width, y, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextSize(textSize);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  return pt;
}
void BeamComp64GeVNeg100ARound2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BeamComp64GeVNeg100ARound2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BeamComp64GeVNeg100ARound2_cxx
