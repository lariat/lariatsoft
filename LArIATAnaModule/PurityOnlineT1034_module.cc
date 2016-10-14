////////////////////////////////////////////////////////////////////////
// Class:       PurityOnlineT1034
// Module Type: analyzer
// File:        PurityOnlineT1034_module.cc
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
#include "lardata/RecoBaseArt/TrackUtils.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
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
#include "larreco/Calorimetry/CalorimetryAlg.h"

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
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

double likely(const double *par)
{

   // Numeric constants
   double mpshift  = -0.22278298;       // Landau maximum 

   int index = 0;
   const double tau=par[0];
   const double sigma=par[1];
   const double dqdxo=par[2];
   double mpc;             // Corrected most probable value
   double dqdx;            // charge after tau scaling
   double prob;            // Landau probability
   double lik;             // likelihood
   double aa;
   double bb;

   vector<double> charge;
   vector<double> dtime;

   ifstream file ("tmpdontcareaboutthis.txt");

   while (file.good() && !file.eof())
   {
      file >> aa >> bb;
      dtime.push_back(aa);
      charge.push_back(bb); 
   }
   index=dtime.size()-1;
   dtime.resize(index);
   charge.resize(index); 

   lik=0;
   for (int ii = 0; ii<index; ++ii)              // Loop over hits
   { 
      dqdx=dqdxo*TMath::Exp(-(dtime.at(ii)/tau));
      mpc = dqdx-mpshift*sigma;                      
      prob = TMath::Landau(charge.at(ii),mpc,sigma,kTRUE); 
      lik-=TMath::Log(prob);
   }     // loop over hits   
   
   file.close();
   dtime.clear();
   charge.clear();     
   return lik;
}

namespace lariat 
{
  class PurityOnlineT1034;
}


class lariat::PurityOnlineT1034 : public art::EDAnalyzer {
public:
   explicit PurityOnlineT1034(fhicl::ParameterSet const & p);
   virtual ~PurityOnlineT1034();
 
   void analyze(art::Event const & e) override;

   // Selected optional functions.
   void beginJob();
   void endJob();

   void reconfigure(fhicl::ParameterSet const & p);

private:

   // Declare member data here.

   void ResetVars();

   FILE *values_file;

   bool   fVerbose;		/// Print - not print some output
   int    run;
   int    subrun;
   int    event;
   int    nhits;
   int    ntracks_reco; 
   int    hit_plane[20000];
   int    hit_wire[20000];
   int    hit_trkkey[20000];
   double hit_peakT[20000];
   double hit_charge[20000];
   double trkpitch[1000][2];
   double trkvtxx[1000];
   double trkendx[1000];

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
   std::vector<double> time;
   std::vector<double> wire;
   std::vector<double> charge;
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

   TH1D *h97;
   TH1D *h971;
   TH1D *h98;
   TH1D *h981;
   TH1D *h99;  

   std::string fHitsModuleLabel;
   std::string fTrackModuleLabel;
   std::string fClusterModuleLabel;

   ROOT::Math::Minimizer* min = 
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

};


lariat::PurityOnlineT1034::PurityOnlineT1034(fhicl::ParameterSet const & pset)
  :EDAnalyzer(pset) 
{
   this->reconfigure(pset);
}

lariat::PurityOnlineT1034::~PurityOnlineT1034()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::PurityOnlineT1034::reconfigure(fhicl::ParameterSet const & pset)
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

void lariat::PurityOnlineT1034::analyze(art::Event const & evt)
{

   // ### Reset variables before we get started ###
   // #############################################
   ResetVars();

   // #############################################
  // a IMultiGenFunction type   
   ROOT::Math::Functor f(&likely,3); 
   // starting point

   // #######################################
   // ### Get potentially useful services ###
   // #######################################
   // === Geometry Service ===
   art::ServiceHandle<geo::Geometry> geom;
   // === Liquid Argon Properties Services ===
   //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
   // === Detector properties service ===
//   auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   // === BackTracker service ===
//   art::ServiceHandle<cheat::BackTracker> bt;
//   const sim::ParticleList& plist = bt->ParticleList();

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
   std::vector<double> trackStart;
   std::vector<double> trackEnd;

   for(size_t i=0; i<tracklist.size();++i)
   {
      // ### Clearing the vectors for each track ###
      trackStart.clear();
      trackEnd.clear();
    
      // ### Setting the track information into memory ###
      tracklist[i]->Extent(trackStart,trackEnd); 
      trkvtxx[i]        = trackStart[0];
      trkendx[i]        = trackEnd[0];    


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
            mf::LogWarning("PurityOnline")<<"caught exeption "<<e<<"\n setting pitch to 0";
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
      if (fmtk.isValid())
      {
         if (fmtk.at(i).size()!=0) hit_trkkey[i] = fmtk.at(i)[0].key();
      }
   }

  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							PURITY CODE
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

   ftime.clear();
   fwire.clear();
   fcharge.clear();
   distancevec.clear();
   distance=0.0;
   chisqr=0.0;
   lowcut = 0.0;
   hicut = 0.0;

   for (int ij = 0; ij<nhits; ++ij)              // Loop over hits
   { 
      if (hit_plane[ij]==1&&hit_trkkey[ij]==0&&(hit_charge[ij]/trkpitch[0][1])<fchargecut)   // Selecting collection plane, first track of the plane
      {
         h99->Fill(hit_charge[ij]/trkpitch[0][1]);
	 ftime.push_back((hit_peakT[ij]-fpretrig)*fsampling);
         fwire.push_back(hit_wire[ij]);
         fcharge.push_back(hit_charge[ij]/trkpitch[0][1]);
      }  // end selection collection plane and long tracks

   }     // loop over hits  

   if (ftime.size() < 80 || fwire.size() < 80 || fcharge.size() < 80)
   {
      if(fVerbose) std::cout << "Not enough hits, skipping..." << endl;
      return;
   }
       
   // Make a fit of wire VS time to find whether track is straight or not 
   TGraph *wgr = new TGraph(ftime.size(),&ftime[0],&fwire[0]);
   TF1 *wfun = new TF1("wfun","pol1",0,400);
   wgr->Fit("wfun");

   for(int n=0; n<2; n++)  // make fit twice to make sure to pick up every hit
   {
      tothit=0;
      time.clear();
      wire.clear();
      charge.clear();
      for(size_t i=0; i< ftime.size(); i++)  //select only hits close to the fit
      {
         fitvalue=wfun->GetParameter(0)+(ftime[i]*wfun->GetParameter(1));
         if(TMath::Abs((fwire[i]-fitvalue)/fwire[i])<0.02)
         {
            tothit++;
            time.push_back(ftime[i]);
            wire.push_back(fwire[i]);
            charge.push_back(fcharge[i]);
         }
      }  //select only hits close to the fit
    
      TGraph *wgr2 = new TGraph(time.size(),&time[0],&wire[0]);
      wgr2->Fit("wfun","R");//,"N");
   } 

   chisqr = wfun->GetChisquare()/wfun->GetNDF();  // obtain chi^2 divided by number of degrees of freedom
   if(chisqr!=chisqr) return;                   // skip the event if the chi square is a nan value
   h98->Fill(wfun->GetParameter(1)); 
   h97->Fill(chisqr);     
 
   ////////////////////////////////////////////////////////////////////////////////
   ///              Here I select the events for the minimization               ///
   ////////////////////////////////////////////////////////////////////////////////

   if (tothit>=80                   &&                   
       chisqr<2.0                   &&
       chisqr!=0.0                  &&       
       ((wfun->GetParameter(1) > 0.3&&wfun->GetParameter(1) < 2.1)||
        (wfun->GetParameter(1) > -2.1&&wfun->GetParameter(1) < -0.3))&&
       ((trkvtxx[0]<6&&trkendx[0]>41)||
        (trkvtxx[0]>41&&trkendx[0]<6)))
   {  
      tothitvec.push_back(tothit);
      h981->Fill(wfun->GetParameter(1)); 
      h971->Fill(chisqr);          
      // Make a fit of charge VS time to find charge at time=0, i.e., dqdxo 

      TGraph *gr = new TGraph(time.size(),&time[0],&charge[0]);
      TF1 *fun = new TF1("fun","expo",0,400);
      gr->Fit("fun"); //,"B");//,"N");

      expo0=fun->GetParameter(0);
      expo1=fun->GetParameter(1);


      TH1D *h1 = new TH1D("h1", "Charge distance", 200, -4000, +4000);
      for (size_t i = 0; i<charge.size(); ++i)  h1->Fill(charge[i]-TMath::Exp(expo0+(expo1*time[i])));
      TF1 *gfun = new TF1("gfun","landau",-2000,4000);
      h1->Fit("gfun");

      landau1=gfun->GetParameter(2);

      gr->Draw("AP");

      TF1 *gfun1;
      TF1 *gfun2;

      if((4*landau1)<1200) gfun1 = new TF1("gfun1","TMath::Exp([0]+([1]*x))-(4*[2])",0,400);    
      else gfun1 = new TF1("gfun1","TMath::Exp([0]+([1]*x))-1200",0,400);                       
      if((4*landau1)<1200) gfun2 = new TF1("gfun2","TMath::Exp([0]+([1]*x))+(4*[2])",0,400);   
      else gfun2 = new TF1("gfun2","TMath::Exp([0]+([1]*x))+1200",0,400);                      

      gfun1->SetParameter(0,expo0);
      gfun1->SetParameter(1,expo1);
      gfun1->SetParameter(2,landau1);
      gfun2->SetParameter(0,expo0);
      gfun2->SetParameter(1,expo1);
      gfun2->SetParameter(2,landau1);

      ////////////////////////////////////////////////////////////////////////////////
      ///  Now I have the cut, move on writing the file for the minimization       ///
      ////////////////////////////////////////////////////////////////////////////////

      if(expo1 < 0.0)  // I don't want an exponential fit of the charge going upward (i.e. more charge at the end of drift time rather than bottom)
      {

         values_file=fopen("tmpdontcareaboutthis.txt","w");

         for (size_t i = 0; i<charge.size(); i++)
         {

            if((4*landau1)<1200) lowcut=TMath::Exp(expo0+(expo1*time[i]))-(4*landau1);
            else lowcut=TMath::Exp(expo0+(expo1*time[i]))-1200;
            if((4*landau1)<1200) hicut=TMath::Exp(expo0+(expo1*time[i]))+(4*landau1);
            else hicut=TMath::Exp(expo0+(expo1*time[i]))+1200;

            if(charge[i] > lowcut && charge[i] < hicut) fprintf(values_file,"%6.2f %6.0f\n",time[i],charge[i]);
         }
            
         fclose(values_file);

      }

      ////////////////////////////////////////////////////////////////////////////////
      ///////////////////     Minimization starts HERE       /////////////////////////
      ////////////////////////////////////////////////////////////////////////////////

      min->SetFunction(f);
 
      // Set the free variables to be minimized!
      min->SetVariable(0,"Tau",fvariable[0], fstep[0]);
      min->SetVariable(1,"Sigma",fvariable[1], fstep[1]);
      min->SetVariable(2,"dqdx0",fvariable[2], fstep[2]);

      lowtaulimit=10;     // lower limit on tau minimization
      hitaulimit=10000;      // upper limit on tau minimization
      lowsigmalimit=0;   // lower limit on sigma minimization
      hisigmalimit=6000;    // upper limit on sigma minimization
      lowdqdxolimit=100;   // lower limit on dqdxo minimization
      hidqdxolimit=40000;    // upper limit on dqdxo minimization

      min->SetVariableLimits(0,lowtaulimit,hitaulimit);
      min->SetVariableLimits(1,lowsigmalimit,hisigmalimit);
      min->SetVariableLimits(2,lowdqdxolimit,hidqdxolimit);

      min->SetVariableValue(2,TMath::Exp(fun->GetParameter(0)));

      // do the minimization
      min->Minimize(); 
 
      const double *xs = min->X();

      if( xs[0]!=xs[0] || xs[1]!=xs[1] || xs[2]!=xs[2] ) 
      {
         if(fVerbose) std::cout << "MINIMIZATION GAVE NAN VALUE!!!! Skipping the event" << endl; 
         return;
      }

      // Get Errors from MiNOS
      min->SetErrorDef(0.5); // If you want the n-sigma level, you have to put here n*n/2. So for 1 sigma, you put 0.5. If you have likelihood. If you have chi2 fit, then you remove the /2. See also https://root.cern.ch/root/html/ROOT__Minuit2__FCNBase.html

      for (int l=0; l<3; l++)
      {
         sigminoserrlow[l] = 0;
         sigminoserrhi[l] = 0;
      }

      min->GetMinosError(0, sigminoserrlow[0], sigminoserrhi[0]); // Get MINOS errors accordingly.
      min->GetMinosError(1, sigminoserrlow[1], sigminoserrhi[1]); // Get MINOS errors accordingly.
      min->GetMinosError(2, sigminoserrlow[2], sigminoserrhi[2]); // Get MINOS errors accordingly.

      if (xs[0]>lowtaulimit+1   && xs[0]<hitaulimit-1   &&	// Fill vectors only if values are not at limit
          xs[1]>lowsigmalimit+1 && xs[1]<hisigmalimit-1 &&
          xs[2]>lowdqdxolimit+1 && xs[2]<hidqdxolimit-1) 
      {
         tauvec.push_back(xs[0]);
         etauveclow.push_back(-sigminoserrlow[0]);
         etauvechi.push_back(sigminoserrhi[0]);

         sigmavec.push_back(xs[1]);
         esigmaveclow.push_back(-sigminoserrlow[1]);
         esigmavechi.push_back(sigminoserrhi[1]);

         dqdxovec.push_back(xs[2]);
         edqdxoveclow.push_back(-sigminoserrlow[2]);
         edqdxovechi.push_back(sigminoserrhi[2]);

         minvec.push_back(min->MinValue());

         tracktime.push_back((run/1.)+(subrun/10000.));
      }   // End fill vectors only if values are not at limit

      min->Clear();
   } // End Event selection 
}

void lariat::PurityOnlineT1034::beginJob()
{
   art::ServiceHandle<art::TFileService> tfs;
   h97 = tfs->make<TH1D>("h97", "Chi Square", 100, 0, 10);
   h971 = tfs->make<TH1D>("h971", "Chi Square", 100, 0, 10);
   h98 = tfs->make<TH1D>("h98", "Track Angle", 100, -3, 3);
   h981 = tfs->make<TH1D>("h981", "Track Angle", 100, -3, 3);
   h99 = tfs->make<TH1D>("h99", "Hit Charge", 600, 0, 6000);  

   h971->SetLineColor(2);
   h981->SetLineColor(2);

   tothitvec.clear();
   minvec.clear();
   dqdxovec.clear();
   edqdxoveclow.clear();
   edqdxovechi.clear();
   tauvec.clear();
   etauveclow.clear();
   etauvechi.clear();
   sigmavec.clear();
   esigmaveclow.clear();
   esigmavechi.clear();
   taucheck.clear();
   likcheck.clear();

   // set tolerance , etc...
   min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000);  // for GSL 
   min->SetTolerance(0.0001);
   min->SetPrintLevel(1);

}

void lariat::PurityOnlineT1034::endJob()
{

   gSystem->Exec("rm tmpdontcareaboutthis.txt");

   if(fVerbose)
   {
      std::cout << endl 
           << endl
           << "////////////////////////////////////////////////////////////" << endl
           << "/////                  SUMMARY                         /////" << endl
           << "////////////////////////////////////////////////////////////" << endl
           << "# Tracks selected: " << tauvec.size() << endl
           << endl 
           << endl;
   }

   avtau=0.0;
   errtau=0.0;
   hitau=0.0;
   lowtau=0.0;

   for(size_t s=0; s<tauvec.size(); s++)
   {
      avtau+=tauvec[s];
      hitau+=tauvec[s]+etauvechi[s];
      lowtau+=tauvec[s]-etauveclow[s];
      if(fVerbose)
      {
         std::cout << " Number of hits passing selection: " << tothitvec[s] << endl
                   << " Tau: " << tauvec[s] << " + " << etauvechi[s] << " - " << etauveclow[s] << endl
                   << " Sigma: " << sigmavec[s] << " + " << esigmavechi[s] << " - " << esigmaveclow[s] << endl
                   << " dQ0/dx: " << dqdxovec[s] << " + " << edqdxovechi[s] << " - " << edqdxoveclow[s] << endl
                   << " Likelihood value @ minimum: " << minvec[s] << endl
                   << "  " << endl;
      }
   } // loop on selected tracks
   avtau/=tauvec.size();
   hitau/=tauvec.size();
   lowtau/=tauvec.size();
   
   entries=tauvec.size();

   std::cout << endl << endl
   << "***************************************************************" << endl
   << "***************************************************************" << endl
   << "******                SINGLE TRACK METHOD                ******" << endl
   << " Average value of tau over " <<  tauvec.size() << " tracks: " << endl
   << " Tau= " << avtau  << " + " << hitau-avtau << " - " << avtau - lowtau << " mus" << endl
   << "***************************************************************" << endl
   << "***************************************************************" << endl;

}

void lariat::PurityOnlineT1034::ResetVars()
{

   run = -99999;
   subrun = -99999;
   event = -99999;
   ntracks_reco = -99999;
   nhits = -99999;

   for (int i = 0; i < 1000; ++i) 
   {
      trkvtxx[i] = -99999;
      trkendx[i] = -99999;
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

DEFINE_ART_MODULE(lariat::PurityOnlineT1034)
