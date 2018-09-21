////////////////////////////////////////////////////////////////////////
// Class:       LifetimeX
// Module Type: analyzer
// File:        LifetimeX_module.cc
//
// Generated at Thu Mar 22 15:06:33 2018 by William Foreman using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom2.h"
#include "TVector2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TF1.h"

// C++ includes
#include <iostream>
#include <memory>
#include <random>
#include <algorithm>

#include "LArIATRecoAlg/TriggerFilterAlg.h"

float RAD_TO_DEG = 180. / 3.1415926535;

// ###################################################################
double LandauGaus(Double_t *x,Double_t *par) {
    
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
    const Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    const Double_t mpshift  = -0.22278298;       // Landau maximum location
    const Double_t np =   500.0;      // number of convolution steps
    const Double_t sc =    5.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    double xx, mpc, fland, xlow, xupp, step, i;
    double sum = 0.0;
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0];
    
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





class LifetimeX;

class LifetimeX : public art::EDAnalyzer {
public:
  explicit LifetimeX(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LifetimeX(LifetimeX const &) = delete;
  LifetimeX(LifetimeX &&) = delete;
  LifetimeX & operator = (LifetimeX const &) = delete;
  LifetimeX & operator = (LifetimeX &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void GetDetProperties();
  
  // fcl-configurable params
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fTrackCalModuleLabel;
  int         fNumDriftBins;
  float       fMarginX;
  float       fMindX;
  float       fMaxdX;
  float       fMinTrkLength;
  float       fMaxTrkPitch;
  int         fMinRun;
  int         fMaxRun;
  float       fMinElectronLifetimeFromDB;
  float       fMaxElectronLifetimeFromDB;
  

private:
  
  geo::GeometryCore     *fGeo;
  
  int         fRunNumber;
  int         fSubRunNumber;
  int         fEventNumber;
  int         fEventTime;

  float       fGeoRangeX[2];
  float       fDriftVelocity;
  float       fXTicksOffset[2];
  float       fSamplingPeriod;
  float       fElectronLifetimeFromDB;
  float       fBinWidth;

  TH1D*               h_EventSelection;
  TH1D*               h_runNumber;
  TH1D*               h_electronLifetimeFromDB;
  TH1D*               h_trkNodeX;
  TH1D*                h_trkdX;
  TH1D*               h_hitAmplitude;
  TH1D*               h_hitRMS;
  TH1D*               h_trkPitch;
  TH1D*               h_trkZenithAngle;
  std::vector<TH1D*>  h_dQdx;
  TH1D*               h_T_vs_dQdx;
  std::vector<float>  fX1;
  std::vector<float>  fX2;
  std::vector<float>  fT1;
  std::vector<float>  fT2;
  std::vector<float>  fTc;


};


//##############################################################################3
LifetimeX::LifetimeX(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
{
   this->reconfigure(pset);
}

//##############################################################################3
void LifetimeX::reconfigure(fhicl::ParameterSet const & pset)
{
   fHitsModuleLabel      	= pset.get< std::string > ("HitsModuleLabel","trajcluster");
   fTrackModuleLabel		= pset.get< std::string > ("TrackModuleLabel","tracktc");
   fTrackCalModuleLabel		= pset.get< std::string > ("TrackCalModuleLabel","calotc");
   fNumDriftBins                 = pset.get< int >         ("NumDriftBins",10);
   fMarginX                     = pset.get< float >       ("MarginX",3.75);
   fMindX                       = pset.get< float >       ("MindX",45.);
   fMaxdX                       = pset.get< float >       ("MaxdX",50.);
   fMinTrkLength                = pset.get< float >       ("MinTrkLength",40.);
   fMaxTrkPitch                 = pset.get< float >       ("MaxTrkPitch",3.);
   fMinRun                      = pset.get< int >         ("MinRun",-999);
   fMaxRun                      = pset.get< int >         ("MaxRun",-999);
   fMinElectronLifetimeFromDB   = pset.get< float >       ("MinElectronLifetimeFromDB",0);
   fMaxElectronLifetimeFromDB   = pset.get< float >       ("MaxElectronLifetimeFromDB",9999);
}



//##############################################################################3
void LifetimeX::beginJob()
{
  LOG_VERBATIM("LifetimeX")
  <<"================================================\n"
  <<"Configuring LifetimeX Module";
 
  if( fMinRun < 0 ) fMinRun = 0;
  if( fMaxRun < 0 ) fMaxRun = 999999;
  
  art::ServiceHandle<art::TFileService> tfs;
  h_EventSelection    = tfs->make<TH1D>("EventSelection","Event selection",4,0,4);
    h_EventSelection  ->SetOption("HIST TEXT");
    h_EventSelection  ->GetXaxis()->SetBinLabel(1,"Total evts"); 
    h_EventSelection  ->GetXaxis()->SetBinLabel(2,"Run num"); 
    h_EventSelection  ->GetXaxis()->SetBinLabel(3,"ACP track"); 
    h_EventSelection  ->GetXaxis()->SetBinLabel(4,"Pitch cut"); 
  h_runNumber         = tfs->make<TH1D>("RunNumber","Run numbers per event",10000,0,10000);
  h_electronLifetimeFromDB = tfs->make<TH1D>("ElectLifetimeFromDB","Electron lifetimes stored in DB",125,0,2500.); 
  h_trkNodeX          = tfs->make<TH1D>("TrkNodeX","Track node X;X [cm]",120,-5.,55.);
  h_trkdX          = tfs->make<TH1D>("TrkdX","Track dx;dx [cm]",160,0.,80.);
  h_hitAmplitude      = tfs->make<TH1D>("HitAmplitude","Hit amplitude",100,0.,200.);
  h_hitRMS            = tfs->make<TH1D>("HitRMS","Hit RMS",100,0.,30.);
  h_trkPitch          = tfs->make<TH1D>("TrkPitch","Pitch in collection plane for crossing tracks;Pitch [cm]",200,0.,20.);
  h_trkZenithAngle    = tfs->make<TH1D>("TrkZenithAngle","Zenith angle of crossing tracks;Angle [deg]",180,0.,180.);

  // Get a pointer to the geometry service provider
  fGeo          = &*(art::ServiceHandle<geo::Geometry>());
  fGeoRangeX[0] = 0.;
  fGeoRangeX[1] = fGeo->DetHalfWidth()*2.;
  
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fElectronLifetimeFromDB   = detprop->ElectronLifetime();  // microseconds
  fDriftVelocity      = detprop->DriftVelocity();
  fSamplingPeriod     = detprop->SamplingRate()*1e-3;
 
  // Set time bin ranges
  float dX        = (fGeoRangeX[1]-fGeoRangeX[0] - 2.*fMarginX) / fNumDriftBins;
  float dT        = dX / fDriftVelocity;
  LOG_VERBATIM("LifetimeX")
  <<"  Efield   = "<<detprop->Efield(0)<<"\n"
  <<"  X1       = "<<fGeoRangeX[0]+fMarginX<<"\n"
  <<"  X2       = "<<fGeoRangeX[1]-fMarginX<<"\n"
  <<"  marginX  = "<<fMarginX<<"\n"
  <<"  driftVel = "<<fDriftVelocity<<" cm/us\n"
  <<"  dX       = "<<dX<<" cm\n"
  <<"  dT       = "<<dT<<" us\n"
  <<"  samp per = "<<fSamplingPeriod<<"\n";
  for(int i=0; i<fNumDriftBins; i++){
    fX1.push_back(fMarginX+i*dX);
    fX2.push_back(fMarginX+(i+1)*dX);
    fT1.push_back(fX1[i] / fDriftVelocity);
    fT2.push_back(fX2[i] / fDriftVelocity);
    fTc.push_back(fT1[i] + dT/2.);
    h_dQdx.push_back(tfs->make<TH1D>(Form("dQdx_%d",i),Form("dQdx for hits in drift range T=(%3.0f,%3.0f) #mus;dQ/dx [ADC/cm]",fT1[i],fT2[i]),150,0,15e3));
    LOG_VERBATIM("LifetimeX")
    <<"  bin "<<i+1<<"  T="<<fT1[i]<<"-"<<fT2[i]<<"  X="<<fX1[i]<<"-"<<fX2[i];
  }
  
  std::cout<<"configuring h_tvsdQdx\n";
  std::cout<<" fT1[0] = "<<fT1[0]<<"\n";
  std::cout<<" fT1["<<fNumDriftBins-1<<"] = "<<fT2[fNumDriftBins-1]<<"\n";
  h_T_vs_dQdx = tfs->make<TH1D>("T_vs_dQdx","Mean dQ/dx vs. drift time;Drift time [#mus];Mean dQ/dx [ADC/cm]",fNumDriftBins,fT1[0],fT2[fNumDriftBins-1]);

  fBinWidth = dT;

  LOG_VERBATIM("LifetimeX")<<"\n"
  <<"================================================";
}



//##############################################################################3
void LifetimeX::analyze(art::Event const & evt)
{
  // ====================================
  // Get run, subrun and event number
  fRunNumber    = (int)evt.run();
  fSubRunNumber = (int)evt.subRun();
  fEventNumber  = (int)evt.event();
  fEventTime    = (int)evt.getSubRun().beginTime().value();
    
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fXTicksOffset[0]    = detprop->GetXTicksOffset(0,0,0);
  fXTicksOffset[1]    = detprop->GetXTicksOffset(1,0,0);

  h_EventSelection->Fill(0);

  if( fRunNumber < fMinRun || fRunNumber > fMaxRun ) return;
  if( detprop->ElectronLifetime() < fMinElectronLifetimeFromDB
      || detprop->ElectronLifetime() > fMaxElectronLifetimeFromDB ) {
    return;
  }
  
  h_runNumber->Fill(fRunNumber);
  h_electronLifetimeFromDB->Fill( detprop->ElectronLifetime() );
  
  h_EventSelection->Fill(1);


        
  /* 
  if( evt.isRealData() ) {
    // Implementation of required member function here.
    //Get the triggers from the event record
    art::Handle< std::vector<raw::Trigger> > triggerHandle;
    evt.getByLabel("daq",triggerHandle);
    // Get Trigger info
    if (triggerHandle.isValid() && triggerHandle->size() > 0){
      uint32_t triggerBits = triggerHandle->at(0).TriggerBits();
      std::bitset<32> triggerBitSet(triggerBits);
      std::cout<<"TRIGGER BIT 13 (MICHEL) = "<<triggerBitSet[13]<<"\n";
      std::cout<<"TRIGGER BIT 12 (LARRY) = "<<triggerBitSet[12]<<"\n";
      std::cout<<"TRIGGER BIT 11 (COSMIC) = "<<triggerBitSet[11]<<"\n";
    }
  }
  */
  
   
  LOG_VERBATIM("LifetimeX")<<"\n"
  <<"Beginning Run "<<fRunNumber<<", subrun "<<fSubRunNumber<<", event "<<fEventNumber;
  
  // ================================================================
  // Get the tracks
  art::Handle< std::vector< recob::Track >> trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}

  // =====================================
  // Getting the hit information
   art::Handle< std::vector<recob::Hit> > hitListHandle; 
   std::vector<art::Ptr<recob::Hit> > hitlist;
   if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      {art::fill_ptr_vector(hitlist, hitListHandle);}

  // ===================================== 
  // Association between Tracks and Hits
  art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
    
  // Get calo object
  
  LOG_VERBATIM("LifetimeX")
  <<"Found "<<tracklist.size()<<" tracks";
  
  // =====================================
  // Loop over tracks
  for(size_t iTrk=0; iTrk<tracklist.size();iTrk++)
  {
    // ------------------------------------------------------------
    // Get the recob::Track object and record its endpoint/vertex
    art::Ptr< recob::Track > TrkPtr(trackListHandle,iTrk);
    recob::Track Trk = *TrkPtr;
    float trkX1 = Trk.Vertex().X();
    float trkX2 = Trk.End().X();
    float dX = fabs(trkX1-trkX2);
    h_trkNodeX->Fill( trkX1 );
    h_trkNodeX->Fill( trkX2 );
    h_trkdX   ->Fill( dX );
    
    std::cout<<"Trk "<<iTrk<<" has dx = "<<dX<<"\n";
  
    // ----------------------------------------------------------------
    // Look for tracks that pass from cathode to anode
    if(
            (trkX1 < fGeoRangeX[0]+fMarginX || trkX2 < fGeoRangeX[0]+fMarginX)
        &&  (trkX1 > fGeoRangeX[1]-fMarginX || trkX2 > fGeoRangeX[1]-fMarginX) ) {
        //dX >= fMindX && dX <= fMaxdX &&  Trk.Length() > fMinTrkLength ){
      
      h_EventSelection->Fill(2);
      
      LOG_VERBATIM("LifetimeX")
      <<"Track passes from cathode-anode:  Xvert= "<<trkX1<<"   Xend= "<<trkX2;
      
      float trkZenithAngle  = -9.;
      float trkPitch        = -9.;

      // ------------------------------------------------- 
      // Calculate the track's X-angle and zenith angle
      TVector3 vert(0.,1.,0.);
//      TVector3 dir = (Trk.Vertex() - Trk.End());
      float dx = fabs(fGeoRangeX[0]-fGeoRangeX[1]);
      float dy = Trk.Vertex().Y()-Trk.End().Y();
      float dz = Trk.Vertex().Z()-Trk.End().Z();
      TVector3 dir(dx,dy,dz);
      dir.SetMag(1.);
      trkZenithAngle = dir.Angle(vert);
      if( trkZenithAngle > TMath::Pi()/2. ) trkZenithAngle = TMath::Pi() - trkZenithAngle;
      h_trkZenithAngle->Fill( trkZenithAngle*RAD_TO_DEG);
      
      // Track pitch on collection plane
      float wirePitch  = fGeo->WirePitch( geo::kV, 0, 0 );
      float wireVertAng = fGeo->WireAngleToVertical( geo::kV, 0, 0 ) - 0.5*TMath::Pi();
      float cosgamma = fabs( sin(wireVertAng)*dir.Y() + cos(wireVertAng)*dir.Z() );
      if( cosgamma > 0 && wirePitch > 0 ){
        trkPitch = wirePitch / cosgamma;   
        h_trkPitch->Fill(trkPitch);
      }

      if( trkPitch <= 0 || trkPitch > fMaxTrkPitch ) continue;
      
      h_EventSelection->Fill(3);

      // -------------------------------------------------------------
      // Create vector of hit keys associated with this track
      std::vector<int> trkHits;
      art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
      if (fmthm.isValid()){
        auto vhit = fmthm.at(iTrk);
        for (size_t h = 0; h < vhit.size(); ++h) trkHits.push_back(vhit[h].key());
      }

      // -----------------------------------------------------------
      

      // ------------------------------------------------------------
      // Loop through the hits associated with this track, and for each
      // hit, calculate its dQ/dx and save into the appropriate histogram
      for(size_t j=0; j<trkHits.size(); j++){
        
        size_t iHit = trkHits.at(j);
        int plane = hitlist[iHit]->WireID().Plane; 

        if( plane != 1 ) continue;

        if( hitlist[iHit]->PeakAmplitude() < 40 ) continue;
    
        h_hitAmplitude->Fill( hitlist[iHit]->PeakAmplitude() );
        h_hitRMS      ->Fill( hitlist[iHit]->RMS() );

        float hitTime     = fSamplingPeriod*(hitlist[iHit]->PeakTime()-fXTicksOffset[plane]);
        float hitdQdx     = hitlist[iHit]->Integral() / trkPitch;

        LOG_VERBATIM("LifetimeX")
        <<"   peak T= "<<hitlist[iHit]->PeakTime()<<"   Xoffset= "<<fXTicksOffset[plane]<<"     hit T= "<<hitTime<<"    dQdx= "<<hitdQdx;
        
        int driftBin       = -1;
        for(int i=0; i<fNumDriftBins; i++){
          if( hitTime >= fT1[i] && hitTime < fT2[i] ) {
            driftBin = i;
            break;
          }
        }
        if( driftBin < 0 ) continue;

        h_dQdx[driftBin]    ->Fill( hitdQdx );
        
      } // end loop over track hits







      break;
    }// end ACP selection

  }// end loop over tracks
}



//##############################################################################3
void LifetimeX::endJob()
{
  std::cout<<"Min run = "<<fMinRun<<", max run = "<<fMaxRun<<"\n";
  std::cout<<"Min LT = "<<fMinElectronLifetimeFromDB<<", max LT = "<<fMaxElectronLifetimeFromDB<<"\n";
  TGraphErrors gr;
  
  TF1 landau("LandauGaus",LandauGaus,0,20000,4);
  
  for(int i=0; i<fNumDriftBins; i++){
    if( h_dQdx[i]->GetEntries() < 50 ) {
      std::cout<<"DRIFT BIN "<<i+1<<" ONLY HAS "<<h_dQdx[i]->GetEntries()<<" ENTRIES --> not filling this bin!\n";
      continue;
    } 
  
    std::cout<<"Filling bin "<<i+1<<"\n";

    // Only use first and last 2 bins... 
    //if( i>1 && i < fNumDriftBins-2) {
    //  std::cout<<"(skipping)\n";
    //  continue;
    //}
    
    float mean  = h_dQdx[i]->GetMean();
    float rms   = h_dQdx[i]->GetRMS(); 

  landau.SetParameter(0,rms/10.); // width scale of landau
  landau.SetParLimits(0,1.,rms);
  landau.SetParameter(1,mean);  // mpv
  landau.SetParLimits(1,0.,mean*1.5);
  landau.SetParameter(2,h_dQdx[i]->GetEntries()*100); // total integral
  landau.SetParameter(3,rms/20.); // convoluted gaussian width
  landau.SetParLimits(3,0.,rms/5.); // convoluted gaussian width
  h_dQdx[i]->Fit("LandauGaus","R");
    float val = landau.GetParameter(1);
    float err = landau.GetParError(1);
  
    gr.SetPoint(gr.GetN(),fTc[i],val);
    gr.SetPointError(gr.GetN()-1, fBinWidth/sqrt(12), err);
    
    h_T_vs_dQdx->Fill( fTc[i], val );
    h_T_vs_dQdx->SetBinError(i+1, err );
    std::cout<<val<<"\n";
  }

  TF1 expFit("expFit","[0]*exp(-x/[1])");
  // limit range to exclude edge bins
  //expFit.SetRange(fT2[0],fT2[fNumDriftBins-2]);
  expFit.SetRange(fT1[0],fT2[fNumDriftBins-1]);
  expFit.SetParameter(0, h_T_vs_dQdx->GetBinContent(0));
  expFit.SetParameter(1, 1000 );
  h_T_vs_dQdx->Fit("expFit","R");
 

  LOG_VERBATIM("LifetimeX")
  <<"\n"
  <<"******************************************************************\n"
  <<"Lifetime from fit  = "<<expFit.GetParameter(1)<<" +/- "<<expFit.GetParError(1)<<" microseconds\n"
  <<"Chi^2 / NDF-1      = "<<expFit.GetChisquare()/(expFit.GetNDF()-1)<<"\n"
  <<"******************************************************************\n";

}

//########################################################################################
void LifetimeX::GetDetProperties(){
    //std::cout<<"XTicksOffset "<<fXTicksOffset[1]<<"\n";
}



DEFINE_ART_MODULE(LifetimeX)
