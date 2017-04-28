////////////////////////////////////////////////////////////////////////
// Class:       Lifetime
// Module Type: analyzer
// File:        Lifetime_module.cc
//
// Generated at Wed Mar 16 02:10:43 2016 by Roberto Acciarri using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "Utilities/DatabaseUtilityT1034.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
//#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
//#include "lardata/AnalysisAlg/CalorimetryAlg.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
#include "TError.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTimeStamp.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include <TStyle.h>
#include <TCanvas.h>

#include <TH2.h>
#include <TH3.h>
#include "TVectorD.h"
#include "TVector3.h"
#include "TProfile2D.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <ctime>



Double_t langaufun(Double_t *x,Double_t *par) {
    
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.
    
    // Numeric constants
    const Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    const Double_t mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    //      Double_t np = 18000.0;      // number of convolution steps
    //      Double_t sc =    8.0;      // convolution extends to +-sc Gaussian sigmas
    
    const Double_t np = 5000.0;      // number of convolution steps
    const Double_t sc =    5.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;
    
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0];
    //mpc = par[1];
    
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    
    step = (xupp-xlow) / np;
    
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    
    return (par[2] * step * sum * invsq2pi / par[3]);
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
Double_t  gaulangaufun( Double_t *x,Double_t *par)  // gauss + langau
{
    //     par[0] = Norm of Pedestal Gauss
    //     par[1] = mean of Pedestal Gauss
    //     par[2] = sigma of Pedestal Gauss
    //     par[3] = Width (scale) parameter of Landau density
    //     par[4] = Most Probable Value (MPV, location) parameter of Landau density
    //     par[5] = Total area (integral -inf to inf, normalization constant)
    //     par[6] = Width (sigma) of convoluted Gaussian function
    
    Double_t *p = &par[3];
    return par[0]*TMath::Gaus(x[0],par[1],par[2]) + langaufun(x,p);
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

Double_t gaulangaufun_area(Double_t *x, Double_t *par)// gauss + langau
{
    Double_t sqrt2PI = 9.8696044010894;
    Double_t *p = &par[3];
    par[1]=2*par[4];
    par[2]=par[3]+par[6];
    return par[0]*TMath::Gaus(x[0],par[1],par[2])/(sqrt2PI*par[2]) + langaufun(x,p);
}

TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   //Once again, here are the Landau * Gaussian parameters:
   //par[0] = Norm of Pedestal Gauss
   //par[1] = mean of Pedestal Gauss
   //par[2] = sigma of Pedestal Gauss
   //par[3] = Width (scale) parameter of Landau density
   //par[4] = Most Probable Value (MPV, location) parameter of Landau density
   //par[5] = Total area (integral -inf to inf, normalization constant)
   //par[6] = Width (sigma) of convoluted Gaussian function   //

   //Variables for langaufit call:
   //his             histogram to fit
   //fitrange[2]     lo and hi boundaries of fit range
   //startvalues[4]  reasonable start values for the fit
   //parlimitslo[4]  lower parameter limits
   //parlimitshi[4]  upper parameter limits
   //fitparams[4]    returns the final fit parameters
   //fiterrors[4]    returns the final fit errors
   //ChiSqr          returns the chi square
   //NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,gaulangaufun_area,fitrange[0],fitrange[1],7);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Norm2G","Mean2G","Sigma2G","Width","MPV","Area","GSigma");
   
   for (i=0; i<7; i++) 
   {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   ffit->FixParameter(1,100); // fix the mean of the second gaussian
   ffit->FixParameter(2,200); // fix sigma of the second gaussian
//   ffit->FixParameter(0,0); // fix the amplitude of the second gaussian
//   ffit->FixParameter(1,0); // fix the mean of the second gaussian
//   ffit->FixParameter(2,1); // fix sigma of the second gaussian

   //his->Fit(FunName,"RB0ME");   // fit within specified range, use ParLimits, do not plot
   his->Fit(FunName,"RLME");   // fit within specified range, use ParLimits, do not plot RLME

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<7; i++) 
   {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
}


const int kMaxHits       = 20000;    //maximum number of hits
namespace lariat 
{
  class Lifetime;
}

class lariat::Lifetime : public art::EDAnalyzer {
public:
   explicit Lifetime(fhicl::ParameterSet const & p);
   virtual ~Lifetime();
 
   void analyze(art::Event const & e) override;

   // Selected optional functions.
  // void beginJob();
   void endJob();

   void reconfigure(fhicl::ParameterSet const & p);

private:

   // Declare member data here.

   void ResetVars();

 //  FILE *values_file;

   bool   fVerbose;		/// Print - not print some output
   int    run;
   int    subrun;
   int    event;
   int    nhits;
   int    ntracks_reco;
  static const int kMaxTracks=260;
   int    ntrkhits[kMaxTracks];
   
   double trkx[kMaxTracks][1000];
   double trky[kMaxTracks][1000];
   double trkz[kMaxTracks][1000];
   int    hit_plane[kMaxHits];
   int    hit_wire[kMaxHits];
   int    hit_trkkey[kMaxHits];
   double hit_ph[kMaxHits];
   double hit_peakT[kMaxHits];
   double hit_charge[kMaxHits];
   double trkpitch[1000][3];
   float  hit_x[kMaxHits];        //<---hit x coordinate
   float  hit_y[kMaxHits];        //<---hit y coordinate
   float  hit_z[kMaxHits];        //<---hit z coordinate
  
   double trkstartdcosx[kMaxTracks];
   double trkstartdcosy[kMaxTracks];
   double trkstartdcosz[kMaxTracks];
   double trkenddcosx[kMaxTracks];
   double trkenddcosy[kMaxTracks];
   double trkenddcosz[kMaxTracks];
   double trkvtxx[kMaxTracks];
   double trkvtxy[kMaxTracks];
   double trkvtxz[kMaxTracks];
   double trkendx[kMaxTracks];
   double trkendy[kMaxTracks];
   double trkendz[kMaxTracks];


   int tothit;
   int minitr;
   int checkitr;
   int entries;
   int fpretrig;             // number of time ticks of pre trigger
   double avfittau;	    // average tau from fit
   double eavfittau;	    // error on average tau from fit
   double avtau;            // final average tau
   double errtau;           // error on final average tau
   double chisqr;
   double fchargecut;        // cut on hitcharge
   double fitvalue;
   double lik;
   double tau;
   double mpc;             // Corrected most probable value
   double dqdx;            // dQ0/dx after tau scaling
   double prob;
   double mintau;
   double minlik;
   double minerrtau;
   double maxerrtau;
   double sigminoserrlow[3];
   double sigminoserrhi[3];
   double lowtaulimit;     // lower limit on tau minimization
   double hitaulimit;      // upper limit on tau minimization
   double lowsigmalimit;   // lower limit on sigma minimization
   double hisigmalimit;    // upper limit on sigma minimization
   double lowdqdxolimit;   // lower limit on dqdxo minimization
   double hidqdxolimit;    // upper limit on dqdxo minimization
   double hitau;
   double lowtau;
   double distance;
   double expo0;
   double expo1;
   double landau1;
   double lowcut;
   double hicut;
   double fsampling;   // TPC fsampling rate for time ticks to microseconds conversion
   std::vector<double> fvariable;
   std::vector<double> fstep;
   std::vector<double> tothitvec;
   std::vector<double> wire;
   std::vector<double> ftime;
   std::vector<double> fwire;
   std::vector<double> fcharge;
   std::vector<double> distancevec;
   std::vector<double> tracktime;
   std::vector<double> etracktimelow;
   std::vector<double> etracktimehi;
   std::vector<double> minvec;
   std::vector<double> dqdxovec;
   std::vector<double> edqdxoveclow;
   std::vector<double> edqdxovechi;
   std::vector<double> tauvec;
   std::vector<double> etauveclow;
   std::vector<double> etauvechi;
   std::vector<double> sigmavec;
   std::vector<double> esigmaveclow;
   std::vector<double> esigmavechi;
   std::vector<double> taucheck;
   std::vector<double> likcheck;
   std::vector<double> vts;
   std::vector<double> vazi;
   std::vector<double> valti;

    int       besttime    = -1;
    int       ngoodtracks = 0;
    int       itime       = 0; 
    static const int tbinsize    = 256;             // 192, 384
    static const int ntbins      = 3072/tbinsize;   // being 3072 the time ticks of a drift window
    int       pretrig     = 189;             // pre trigger time ticks. 189= 24.2 mus of pre-sampling
    int       hitcut      = 5;               // minimum number of hits in track 
    double    aa          = 0.0;
    double    bb          = 0.0;
    double    chargecut   = 6000; //1400.; Run 1 value           // charge threshold to flag a hit as belonging to a delta
    double    sampling    = 0.128;           // TPC sampling rate in mus
    double    trkzenith   = 0.0;             // track angle
    double    tbangle     = 0.0;             // track angle respect to beam
    double    strght      = 0.0;             // Track straightness
    double    length      = 0.0;             // Track length
    double    chisqrc     = 0.0;             // ChiSquare for fits
    double    chisqri     = 0.0;             // ChiSquare for fits
    double    cfact[240];                    // Correction factor for wire charge non uniformity
    int       goodtrk; 
    int       nentry      = 0;
    double    pedestalvalue = 504;
   	                	   
    Bool_t Object_exists;

    double year  = 0.0;
    double month = 0.0;
    double day   = 0.0;
    double hour  = 0.0;
    double min   = 0.0;
    double sec   = 0.0;

    TH1D *hitcharge[9999][ntbins];
    TH1D *hittime[9999][ntbins];          // Rob -> I only need this 
    TFile *myfile[9999][ntbins]; 
    TH1D *hitcharge_read[9999][ntbins];

    TH1D *hzenith = new TH1D("hzenith","Track Zenith Angle",360,0,180);
    TH1D *htbangle = new TH1D("htbangle","Track Angle wrt beam",360,0,180);
    TH1D *hstraight = new TH1D("hstraight","Track straightness",400,0.95,1.05);
    TH1D *htrklen = new TH1D("htrklen","Track Length",260,-10,120);
    TH1D *hchi2_col = new TH1D("hchi2_col","ChiSquare Collection Plane",500,0,500);
    TH1D *hchi2_ind = new TH1D("hchi2_ind","ChiSquare Induction Plane",500,0,500);
    
    TH1D *hxstart = new TH1D("hxstart","Track X-start position",208,-2,50);
    TH1D *hxend = new TH1D("hxend","Track X-end position",208,-2,50);
    TH3D *htrkstart = new TH3D("htrkstart","Track Start position",416,-2,102,208,-2,50,176,-22,22);
    TH3D *htrkend = new TH3D("htrkend","Track End position",416,-2,102,208,-2,50,176,-22,22);
    TH3D *htrk = new TH3D("htrk","Track",416,-2,102,208,-2,50,176,-22,22);
    TH1D *chargetonoise_wire = new TH1D("chargetonoise_wire","Charge To Noise (Wire)", 100, 0.,100.);
    TH1D *chargetonoise_cathode = new TH1D("chargetonoise_cathode","Charge to Noise (Cathode)",100,0.,100.);

    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fClusterModuleLabel;

};


lariat::Lifetime::Lifetime(fhicl::ParameterSet const & pset)
  :EDAnalyzer(pset) 
{
   this->reconfigure(pset);
}

lariat::Lifetime::~Lifetime()
{
//ResetVars();
  // Clean up dynamic memory and other resources here.
}

void lariat::Lifetime::reconfigure(fhicl::ParameterSet const & pset)
{
   fHitsModuleLabel      	= pset.get< std::string >("HitsModuleLabel");
   fTrackModuleLabel		= pset.get< std::string >("TrackModuleLabel");
   fClusterModuleLabel          = pset.get< std::string >("ClusterModuleLabel");

   fchargecut                   = pset.get< double               >("ChargeCut");
   fsampling                    = pset.get< double               >("SamplingTime");
   fpretrig                     = pset.get< int                  >("PreSampling");
   fvariable                    = pset.get< std::vector<double>  >("Variable");
   fstep                        = pset.get< std::vector<double>  >("Step");
   fVerbose                     = pset.get< bool                 >("Verbose", false);

   return;
}
void lariat::Lifetime::ResetVars()
{

   run = -99999;
   subrun = -99999;
   event = -99999;
   ntracks_reco = -99999;
   nhits = -99999;
 //  evttime = 9999;		

   for (int i = 0; i < 1000; ++i)
   {
     if(i < kMaxTracks){
       trkvtxx[i] = -99999;
       trkendx[i] = -99999;
     }
     for (int j = 0; j<2; ++j) trkpitch[i][j] = -99999;
   }

   for (int i = 0; i<20000; ++i)
   {
      hit_plane[i] = -99999;
      hit_wire[i] = -99999;
      hit_peakT[i] = -99999;
      hit_charge[i] = -99999;
      hit_trkkey[i] = -99999;
   }

}


void lariat::Lifetime::analyze(art::Event const & evt)
{

   // ### Reset variables before we get started ###
   // #############################################
   ResetVars();



   // #############################################
  // a IMultiGenFunction type   
  // ROOT::Math::Functor f(&gaulangaufun,3); 
   // starting point

   // #######################################
   // ### Get potentially useful services ###
   // #######################################
   // === Geometry Service ===
   art::ServiceHandle<geo::Geometry> geom;
   // === Liquid Argon Properties Services ===
   //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
   // === Detector properties service ===
  // auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   // === BackTracker service ===
   art::ServiceHandle<cheat::BackTracker> bt;
  // const sim::ParticleList& plist = bt->ParticleList();

   // === Run Number ===
   run = evt.run();
   // === Sub-Run Number ===
   subrun = evt.subRun();
   // === Event Number ===
   event = evt.id().event();
   
   std::cout << std::endl
             << "=========================================" << std::endl
             << "Run = " << run << ", SubRun = " << subrun << ", Evt = " << event << std::endl
             << "=========================================" << std::endl 
             << std::endl;

   // #####################################
   // ### Getting the Track Information ###
   // #####################################
   art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
   // === Filling the tracklist from the tracklistHandle ===
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}

   // ###################################
   // ### Getting the Hit Information ###
   // ###################################
   art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define tracklist as a pointer to recob::tracks

   // === Filling the hitlist from the hitlistHandle ===
   if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      {art::fill_ptr_vector(hitlist, hitListHandle);}

   // ##########################################
   // ### Getting the 2D Cluster Information ###
   // ##########################################
   art::Handle< std::vector<recob::Cluster> > clusterListHandle; //<---Define clusterListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Cluster> > clusterlist; //<---Define cluster as a pointer to recob::cluster

   // === Filling the clusterlist from the clusterlistHandle ===
   if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      {art::fill_ptr_vector(clusterlist, clusterListHandle);}

   // ##########################################################
   // ### Grabbing associations for use later in the AnaTool ###
   // ##########################################################
   // === Associations between hits and raw digits ===
   art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
   // === Association between SpacePoints and Tracks ===
   art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
   // === Association between Tracks and 2d Hits ===
   art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
   // ==== Association between Clusters and Hits ===
   art::FindManyP<recob::Cluster>     fmc(hitListHandle,   evt, fClusterModuleLabel);
   // ==== Association between Tracks and Hits
   art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
   art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);

  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE EVENT INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------



   if(fVerbose) std::cout << "Tracklist size " << tracklist.size() << std::endl;

   ntracks_reco=tracklist.size();
   double larStart[3];
   double larEnd[3];
   std::vector<double> trackStart;
   std::vector<double> trackEnd;
   double    sampling    = 0.128;
   int       pretrig     = 189; 
   double    chargecut   = 6000;


   for(size_t i=0; i<tracklist.size();++i)
   {
      // ### Clearing the vectors for each track ###
      trackStart.clear();
      trackEnd.clear();
    
      // ### Setting the track information into memory ###
      auto extent = tracklist[i]-> Extent();
      trkvtxx[i]        = extent.first.X();
      trkvtxy[i]        = extent.first.Y();
      trkvtxz[i]        = extent.first.Z();

      trkendx[i]        = extent.second.X();
      trkendy[i]        = extent.second.Y();
      trkendz[i]        = extent.second.Z();

      // ### Recording the directional cosine at the start of the track ###
      memset(larStart, 0, 3);
      memset(larEnd, 0, 3);
      tracklist[i]->Direction(larStart,larEnd);
      trkstartdcosx[i]  = larStart[0];
      trkstartdcosy[i]  = larStart[1];
      trkstartdcosz[i]  = larStart[2];
    
      // ### Recording the directional cosine at the end of the track ###
      trkenddcosx[i]    = larEnd[0];
      trkenddcosy[i]    = larEnd[1];
      trkenddcosz[i]    = larEnd[2];

    
      // ### Grabbing the SpacePoints associated with this track ###
      ntrkhits[i] = fmsp.at(i).size();
      std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
      


    //############################################
    //###Retriveing 3D information for the hits###
    //############################################

        if (fmthm.isValid()){
           auto vhit = fmthm.at(i);
           auto vmeta = fmthm.data(i);
           for (size_t h = 0; h < vhit.size(); ++h){
            	if (vhit[h].key()<kMaxHits){
 hit_trkkey[vhit[h].key()] = tracklist[i].key();
 hit_x[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).X();
 hit_y[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Y();
 hit_z[vhit[h].key()] = tracklist[i]->LocationAtPoint(vmeta[h]->Index()).Z();

            }
     }//loop over all hits 
 }//fmthm is valid 



      // ###########################################
      // ### Calculating the pitch for the track ###
      // ###########################################
      // === Looping over our two planes (0 == Induction, 1 == Collection) ===
      for (int j = 0; j<2; ++j)
      {
         // ### Putting in a protective try in case we can't set the pitch ###
	 try
	 {
	    // ### If we are in the induction plane calculate the tracks pitch in that view ###
	    if (j==0) trkpitch[i][j] = lar::util::TrackPitchInView(*tracklist[i], geo::kU);
	    // ### If we are in the collection plane calculate the tracks pitch in that view ###
	    else if (j==1) trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kV);
	 }//<---End Try statement
	 catch( cet::exception &e)
	 {
            mf::LogWarning("Lifetime")<<"caught exeption "<<e<<"\n setting pitch to 0";
	    trkpitch[i][j] = 0;
 	 }//<---End catch statement
      }// <---End looping over planes (j)
   }// <---End looping over tracks   

   nhits = hitlist.size();

   for (size_t i = 0; i<hitlist.size() && int(i)< 20000; ++i)
   {
      cet::maybe_ref<raw::RawDigit const> rdref(ford.at(i));
      geo::WireID wireid = hitlist[i]->WireID();

      hit_plane[i]   = wireid.Plane;
      hit_wire[i]    = wireid.Wire;
      hit_peakT[i]   = hitlist[i]->PeakTime();
      hit_charge[i]  = hitlist[i]->Integral();
      hit_ph[i]      = hitlist[i]->PeakAmplitude();


      if (fmtk.isValid())
      {
         if (fmtk.at(i).size()!=0) hit_trkkey[i] = fmtk.at(i)[0].key();
      }

   }
    
  
    // ##################################################
    // #####Applying correction factors to the wires#####
    // ##################################################
    
    for (int i = 0; i<240; ++i)
    {
        cfact[i] = 1.0; //assign start value to correction factor. 1 = no correction applied
    }
    
    //   while (correctionfile.good() && !correctionfile.eof())
    //   {
    //      correctionfile >> aa >> bb;
    //      cfact[int(aa)]=1.000/bb;
    //   }

    
    //#########################################
    //###Checking if they are crossing muons###
    //#########################################
    
//    for (Long64_t jentry=0; jentry<nentries;jentry++)   // Loop on events -- Não preciso mais! Estamos rodando dentro de um analyser
        //   for (Long64_t jentry=1000; jentry<2000;jentry++)   // Loop on events
// { -- Era o parênteses do loop on the events
        
   std::vector<double> timecol;
   std::vector<double> wirecol;
   std::vector<double> timeind;
   std::vector<double> wireind;

        timecol.clear();
        wirecol.clear();
        timeind.clear();
        wireind.clear();
        chisqrc = 0.0;
        chisqri = 0.0;


	 for (int i = 0; i<ntbins; ++i)
         {

//	       hitcharge[itime][i] = new TH1D(Form("hitcharge_%d_%d",itime,i),Form("hitcharge_%d_%d",itime,i),80,0,4000); Run 1 settings
	       hitcharge[itime][i] = new TH1D(Form("hitcharge_%d_%d",itime,i),Form("hitcharge_%d_%d",itime,i),120,0,6000);
	       hitcharge[itime][i]->Sumw2();
               
	    }

         int numtrk = 0;
         int itrk = -1;
         
         for (int i = 0; i < ntracks_reco; ++i)
          {
          if (ntrkhits[i] >= hitcut)
          {
          ++numtrk;
          itrk = i;
          }
         }
        
        
        trkzenith = acos(trkstartdcosy[itrk])*180/TMath::Pi();
        tbangle = acos(trkstartdcosz[itrk])*180/TMath::Pi();
        strght = trkstartdcosx[itrk]*trkenddcosx[itrk]+
        trkstartdcosy[itrk]*trkenddcosy[itrk]+
        trkstartdcosz[itrk]*trkenddcosz[itrk];
        length = sqrt(pow(trkvtxx[itrk]-trkendx[itrk],2)+
                      pow(trkvtxy[itrk]-trkendy[itrk],2)+
                      pow(trkvtxz[itrk]-trkendz[itrk],2));
        
        std::map<int, int >hitmap;
        std::map<int, int >deltahit;
        int minw = 999999;
        int maxw = -999999;

	
        if (((trkvtxx[itrk] < 2.0   && trkvtxx[itrk] > -0.5)  || (trkvtxx[itrk] > 45.5 && trkvtxx[itrk] < 48.5)  ||
             (trkvtxy[itrk] < -18.0 && trkvtxy[itrk] > -20.5) || (trkvtxy[itrk] > 18.0 && trkvtxy[itrk] < 20.5)  ||
             (trkvtxz[itrk] < 2.0   && trkvtxz[itrk] > -0.5)  || (trkvtxz[itrk] > 88.0 && trkvtxz[itrk] < 90.5)) &&
            ((trkendx[itrk] < 2.0   && trkendx[itrk] > -0.5)  || (trkendx[itrk] > 45.5 && trkendx[itrk] < 48.5)  ||
             (trkendy[itrk] < -18.0 && trkendy[itrk] > -20.5) || (trkendy[itrk] > 18.0 && trkendy[itrk] < 20.5)  ||
             (trkendz[itrk] < 2.0   && trkendz[itrk] > -0.5)  || (trkendz[itrk] > 88.0 && trkendz[itrk] < 90.5)))    // Requiring crossing tracks
        {
          
            
            for (int i = 0; i < nhits; ++i)              // Loop over hits
            {
                // Fill vectors for induction and collection plane to check straightness
                if( hit_trkkey[i] == itrk                                                       &&
                   ((trkzenith > 20 && trkzenith < 80) || (trkzenith > 100 && trkzenith < 160)) && // select tracks not-vertical and not-beam alligned
                   ((tbangle > 20 && tbangle < 80) || (tbangle > 100 && tbangle < 160))         && // select tracks not-drift and not-beam alligned
                   strght > 0.99                                                                && // reject not-straight tracks
                   length > 10
                   )
                {
                    if (hit_plane[i] == 1)   // Fill vectors for collection plane
                    {
                        timecol.push_back((hit_peakT[i]-pretrig)*sampling);
                        wirecol.push_back(hit_wire[i]);
                    }
                    else if (hit_plane[i] == 0) // Fill vectors for induction plane
                    {

                        timeind.push_back((hit_peakT[i]-pretrig)*sampling);
                        wireind.push_back(hit_wire[i]);

                    }
                } // Fill vectors for induction and collection plane to check straightness
              
                if (hit_plane[i] == 1 && hit_trkkey[i] >= 0)   // Selecting collection plane
                { 
                    if (ntrkhits[hit_trkkey[i]] >= hitcut)        // Requiring minimum number of hits in the track
                    {
                        hitmap[hit_wire[i]]++;
                        if (hit_wire[i]>maxw) maxw = hit_wire[i];
                        if (hit_wire[i]<minw) minw = hit_wire[i];
                        if ((hit_charge[i]/trkpitch[hit_trkkey[i]][1]) > chargecut) 
				{
				deltahit[hit_wire[i]] = 1;  // If too much charge flag as delta
				}				

                    } //minimum number of hits in track requirement
                } // collection plane selection
            } //loop over hits

            // make straightness fit only if vectors are filled
            if (timecol.size() > 0)
            {
                // Make a fit of wire VS time both on collection and induction plane to find whether track is straight or not
                TGraph *wgrcol = new TGraph(timecol.size(),&timecol[0],&wirecol[0]);
                TF1 *wfunc = new TF1("wfunc","pol1",0,400);
                wgrcol->SetTitle(Form("Fit of Track # %d",nentry+1));
                wgrcol->GetXaxis()->SetTitle("Time (#mus)");
                wgrcol->GetYaxis()->SetTitle("Wire number");
                wgrcol->SetMarkerStyle(22);
                wgrcol->SetMarkerSize(1.5);
                //            wgrcol->Draw("AP");
                wgrcol->Fit("wfunc");
                
                chisqrc=wfunc->GetChisquare()/(wfunc->GetNDF());
                
                TGraph *wgrind = new TGraph(timeind.size(),&timeind[0],&wireind[0]);
                TF1 *wfuni = new TF1("wfuni","pol1",0,400);
                wgrind->SetTitle(Form("Fit of Track # %d",nentry+1));
                wgrind->GetXaxis()->SetTitle("Time (#mus)");
                wgrind->GetYaxis()->SetTitle("Wire number");
                wgrind->SetMarkerStyle(22);
                wgrind->SetMarkerSize(1.5);
                //            wgrind->Draw("AP");
                wgrind->Fit("wfuni");
                
                chisqri=wfuni->GetChisquare()/(wfuni->GetNDF());
		
		nentry++;
            }  // make straightness fit only if vectors are filled
        } // Requiring crossing tracks
        
        if (((trkvtxx[itrk] < 2.0   && trkvtxx[itrk] > -0.5)  || (trkvtxx[itrk] > 45.5 && trkvtxx[itrk] < 48.5)  ||
             (trkvtxy[itrk] < -18.0 && trkvtxy[itrk] > -20.5) || (trkvtxy[itrk] > 18.0 && trkvtxy[itrk] < 20.5)  ||
             (trkvtxz[itrk] < 2.0   && trkvtxz[itrk] > -0.5)  || (trkvtxz[itrk] > 88.0 && trkvtxz[itrk] < 90.5)) &&
            ((trkendx[itrk] < 2.0   && trkendx[itrk] > -0.5)  || (trkendx[itrk] > 45.5 && trkendx[itrk] < 48.5)  ||
             (trkendy[itrk] < -18.0 && trkendy[itrk] > -20.5) || (trkendy[itrk] > 18.0 && trkendy[itrk] < 20.5)  ||
             (trkendz[itrk] < 2.0   && trkendz[itrk] > -0.5)  || (trkendz[itrk] > 88.0 && trkendz[itrk] < 90.5)))    // Requiring crossing tracks
        {
            for (int i = 0; i < nhits; ++i)             // Loop over track hits
            {
                if (hit_plane[i] == 1 && hit_trkkey[i] >= 0)   // Selecting collection plane
                {
                    if (ntrkhits[hit_trkkey[i]] >= hitcut)      // Requiring minimum number of hits in the track
                    {
                        bool deltalike = deltahit[hit_wire[i]] == 1 && deltahit[hit_wire[i]-1] == 1 && deltahit[hit_wire[i]+1] == 1;  // if the wire considered, the one before and after have too much charge, flag as a delta

                       if (deltahit[hit_wire[i]] == 1 && deltahit[hit_wire[i]+1] == 1 && deltahit[hit_wire[i]+2] == 1) deltalike = true;
                       if (deltahit[hit_wire[i]] == 1 && deltahit[hit_wire[i]-1] == 1 && deltahit[hit_wire[i]-2] == 1) deltalike = true;

                        // Cuts
                        if (hitmap[hit_wire[i]] == 1                                                     && // select only collection hits
                            hit_wire[i] != minw                                                          && // reject first hit of the track
                            hit_wire[i] != maxw                                                          && // reject last hit of the track
                            //                      (hit_wire[i] < 90 || hit_wire[i] > 160)                                      && // reject weird middle part
                            deltalike != true                                                            && // reject deltas alligned to the track
                            ((trkzenith > 20 && trkzenith < 80) || (trkzenith > 100 && trkzenith < 160)) && // select tracks not-vertical and not-beam alligned
                            ((tbangle > 20 && tbangle < 80) || (tbangle > 100 && tbangle < 160))         && // select tracks not-drift and not-beam alligned
                            (strght > 0.99 && chisqrc < 10 && chisqri  < 10)                             && // reject not-straight tracks
                            //                      strght > 0.99                                                                && // reject not-straight tracks
                            length > 10) // Cuts                                                            // reject tracks less than 10 cm long
                        {
                            if (hit_ph[i]/trkpitch[hit_trkkey[i]][1] > 0)  // if hit height/trackpitch>0 also track charge is surely above zero
                            {
                                if (hit_charge[i]/trkpitch[hit_trkkey[i]][1] > 0)
                                {

                                   // myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)]    = new TFile(Form("hitcharge_run_%d_%d_%d.root",run, itime, int((hit_peakT[i]-pretrig)/tbinsize)), "UPDATE");
                                    myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)]    = new TFile(Form("hitcharge_run_%d.root",run), "UPDATE");

				    if ((Object_exists = myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)]-> GetListOfKeys()-> Contains(Form("hitcharge_%d_%d", itime, int((hit_peakT[i]-pretrig)/tbinsize)))) == true) //checking if the histograms already exists so they are not recreate
                                   {

                                    hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)] = (TH1D*)myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)] -> Get(Form("hitcharge_%d_%d", itime, int((hit_peakT[i]-pretrig)/tbinsize))); //using the information of the existing histogram

                                    hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)] -> Fill((hit_charge[i]/trkpitch[hit_trkkey[i]][1])*cfact[hit_wire[i]]); // Filling important histograms - the only important one

                                    hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)] -> Write(0,2,0);

		                    myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)]-> Close();
				    }

				   else //if the histogram already exists, only fill the old histogram
                                     {
					hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)] -> Fill((hit_charge[i]/trkpitch[hit_trkkey[i]][1])*cfact[hit_wire[i]]); // Filling important histograms - the only important one

                                        hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)]->SetXTitle("dQ/dx (ADC/cm)");
                                        hitcharge[itime][int((hit_peakT[i]-pretrig)/tbinsize)] -> Write(0,2,0);
		                        myfile[itime][int((hit_peakT[i]-pretrig)/tbinsize)]-> Close();
				    }
   
                            } // if hit height/trackpitch>0 also track charge is surely above zero
                        }  // end of cuts
                    }    // requiring minimum number of hits in the track     
                }       // selecting collection plane
            }          // Loop over track hits
        }             // selecting crossing track
    

  //  } // End loop on events -- Tirei o loop on events! Estamos num analyze e ele roda sozinho sobre os eventos
 
} // end Loop class


}

void lariat::Lifetime::endJob()
{

cout << "***************************************************************"<< endl
     << "                      Histograms created!                      "<< endl
     << "***************************************************************"<< endl;
}

DEFINE_ART_MODULE(lariat::Lifetime)

