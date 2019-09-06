////////////////////////////////////////////////////////////////////////
// Class:       AveragePMTWfm
// Module Type: filter
// File:        AveragePMTWfm_module.cc
//
// Generated Fri June 7, 2019
// W. Foreman
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>

// ROOT includes
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"

//LAriatSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"

class AveragePMTWfm;

class AveragePMTWfm : public art::EDAnalyzer {
public:
  explicit AveragePMTWfm(fhicl::ParameterSet const & p);
  AveragePMTWfm(AveragePMTWfm const &) = delete;
  AveragePMTWfm(AveragePMTWfm &&) = delete;
  AveragePMTWfm & operator = (AveragePMTWfm const &) = delete;
  AveragePMTWfm & operator = (AveragePMTWfm &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void AddToAverageWfm( std::vector<float> &, short, short, TH1F*, int&, int polarity);

private:
  
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory waveformDir = tfs->mkdir("Waveforms");
  
  // Alg object
  OpHitBuilderAlg     fOpAlg;
  
  // Counters and things
  int    fRunNumber;
  int    fSubRunNumber;
  int    fEventNumber;
  int    fNumActiveEvents[10];
  int    fNumSavedWaveforms[10];
 
  // Timetamp-related variables  
  double fTimestamp;
  
  // Histograms
  TH1F* h_Timestamps;
  TH1F* h_Runs;
  TH1F* h_OpHitTimes[10];
  TH1F* h_AveWfm[10];

  short   fAvePulse_preWindow;
  short   fAvePulse_postWindow;
  int     fAvePulse_total;
  
  std::vector<raw::OpDetPulse> fPulses;

  // Tunable parameters defined by fcl
  size_t              fMaxSavedWaveforms;
  float               fADCThresh;
  float               fSavedWaveformADCThresh;
  std::string         fDAQModule;
  std::string         fInstanceName;

  
  std::vector<std::string> fSelectEventTypes; 
  std::vector<size_t> fSelectChannels;      // PMT channels to look at (0=HMM, 1=ETL)
  short               fTruncateWfm;         // truncate PMT wfms to this length
  short               fBaselineWindowLength;// pedestal region of wfm for baseline/RMS calculation
  std::vector<float>  fMaxWfmRMS;           // max wfm RMS (ADC)
  int                 fWfmSmoothingRange;   // local window for waveform smoothing (+/- val)
  std::string         fCorrectOvershootMode;// overshoot correciton mode ("int","amp") 
  bool                fMskBaselineSubtr;    // do masked baseline subtraction on wfms
  std::vector<float>  fGradHitThresh;       // gradient threshold for hit-finding (PMT-specific)
  std::vector<float>  fMinOpHitWidth;       // min ophit FWHM
  
  TCanvas * c;
  std::map<size_t,size_t> iCh;
  char	              buffer[200], histName[200], histTitle[200];

};


//#######################################################################
void AveragePMTWfm::reconfigure(fhicl::ParameterSet const & p)
{
  fAvePulse_preWindow     = p.get< short >        ("AvePulse_preWindow",700);
  fAvePulse_postWindow    = p.get< short >        ("AvePulse_postWindow",7000);
  fSelectChannels         = p.get< std::vector<size_t> >  ("SelectChannels", {1} );
  fTruncateWfm            = p.get< short >                ("TruncateWfm",18000);
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fMaxWfmRMS              = p.get< std::vector< float >>  ("MaxWfmRMS", {2.0} );
  fWfmSmoothingRange      = p.get< int  >                 ("WfmSmoothingRange",2);
  fCorrectOvershootMode   = p.get< std::string >          ("CorrectOvershootMode", "");
  fGradHitThresh          = p.get< std::vector< float >>  ("GradHitThresh",{-15.} );

  fMaxSavedWaveforms      = p.get< size_t > ("MaxSavedWaveforms",10000);
  fADCThresh              = p.get< float > ("ADCThresh",50);
  fSavedWaveformADCThresh = p.get< float > ("SavedWaveformADCThresh",100);

  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
}


//#######################################################################
AveragePMTWfm::AveragePMTWfm(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
{
  this                ->reconfigure(p);
  
  fAvePulse_total = int(fAvePulse_preWindow + fAvePulse_postWindow);

  fRunNumber = 0;
  fSubRunNumber = 0;
  fEventNumber = 0; 
  fPulses.reserve(10);
  
  for(size_t i=0; i< fSelectChannels.size();i++) {
    iCh[fSelectChannels[i]]    = i; // map from channel to index
    fNumSavedWaveforms[i]     = 0;
    fNumActiveEvents[i]       = 0;
  }

}


//#######################################################################
void AveragePMTWfm::beginJob()
{
  h_Timestamps  = tfs->make<TH1F>("Timestamps","Timestamp for all events;sec",240,0.,60.);
  h_Runs        = tfs->make<TH1F>("Runs","Runs analyzed;Run number",10000,4000,14000);
  
  for(size_t i=0; i<fSelectChannels.size(); i++){
    size_t  ch  = fSelectChannels.at(i);
    h_OpHitTimes[ch] = waveformDir.make<TH1F>(
      Form("%lu_OpHitTimes",ch),
      Form("Optical hit times for opdet %lu;Time [ns]",ch),
      5600, 0, 28000);
    h_AveWfm[ch] = waveformDir.make<TH1F>(
      Form("%lu_AveWfm",ch),
      Form("Averaged waveform for opdet %lu;Time [ns];Averaged signal [ADC]",ch),
      fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
  }
}


//#######################################################################
void AveragePMTWfm::analyze(art::Event const & e)
{
  // Get run, subrun and event number
  fRunNumber    = (int)e.run();
  fSubRunNumber = (int)e.subRun();
  fEventNumber  = (int)e.event();
  h_Runs->Fill(fRunNumber);
  
  std::cout<<"event "<<fEventNumber<<"\n";
  
  // Get the PMT data, saving the OpDetPulses specified in fhicl
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);
  if( (size_t)WaveformHandle->size() == 0 ){
    std::cout << "No optical detector data found; skipping event.\n";
    return;
  }
  for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){
    art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
    raw::OpDetPulse thePulse = *ThePulsePtr;
    for( size_t i=0; i<fSelectChannels.size(); i++){
      if( fSelectChannels[i] == thePulse.OpChannel() ) fPulses.push_back(thePulse);
    }
  }
  
  // Grab the first OpDetPulse in the handle just to check whether its
  // the right event type for use in this analysis
  fTimestamp   = (float(WaveformHandle->at(0).PMTFrame())*8.)/1.0e09;
  h_Timestamps ->Fill(fTimestamp);
   
  // Loop over each PMT separately
  for( size_t ipmt = 0; ipmt < fPulses.size(); ipmt++){

    size_t ch   = fPulses[ipmt].OpChannel(); 
    size_t ich  = iCh[ch];
    
    // Set proper gradient hit threshold for hit-finding
    fOpAlg.Reset();
    fOpAlg.SetGradHitThresh( fGradHitThresh[ich] );
      
    // Save the PMT waveform
    if( fTruncateWfm == -1 ) fTruncateWfm = fPulses[ipmt].Waveform().size();
    std::vector<short> PMT_wfm_raw(fPulses[ipmt].Waveform().begin(), fPulses[ipmt].Waveform().begin() + fTruncateWfm); 
    std::vector<float> PMT_wfm(PMT_wfm_raw.begin(), PMT_wfm_raw.end());
      
    // Smooth out high-freq noise
    fOpAlg.SmoothOutVector(PMT_wfm,fWfmSmoothingRange);
      
    // Get waveform pedestal/RMS, and subtract off pedestal
    fOpAlg  .CalcBaselineAndRMS( PMT_wfm, 0, fBaselineWindowLength );
    fOpAlg  .SubtractBaseline( PMT_wfm );

    // PMT overshoot correction (HMM ONLY, CH 0)
    if( fCorrectOvershootMode != "" && ch == 0 ){ 
      std::vector<short> hitTimesTmp = fOpAlg.GetHits(PMT_wfm, (size_t)fPulses[ipmt].FirstSample() );
      fOpAlg.CorrectWfmOvershoot(PMT_wfm, hitTimesTmp, fCorrectOvershootMode );
    }
    
    /*   
    // Masked baseline subtraction
    if( fMskBaselineSubtr ) {
        fvbs[ch].resize(fPMT_wfm[ch].size());
        fOpAlg.MaskedBaselineSubtraction(fPMT_wfm[ch], fvbs[ch]);
        for(size_t i=0; i<fPMT_wfm[ch].size(); i++) fPMT_wfm[ch].at(i) = fPMT_wfm[ch].at(i) - fvbs[ch].at(i);
      }
    */
      
    // determine if this event is "active"
    std::cout<<"   pmt "<<ch<<" --> checking if active...\n";
    for(size_t i=0; i<PMT_wfm.size();i++){
      if( fabs(PMT_wfm.at(i)) >= fADCThresh ) {
        std::cout<<"Active!\n";
        fNumActiveEvents[ich]++; 
        break;
      }
    }

    // hit finding and filtering
    std::vector<short> hitTimes = fOpAlg.GetHits(PMT_wfm, (size_t)fPulses[ipmt].FirstSample() );
     
    std::cout<<"We found "<<hitTimes.size()<<" hits\n"; 
    if( hitTimes.size() == 1 ) {  
        
      float hittime = hitTimes.at(0);
      float a = fOpAlg.GetLocalMaximum(PMT_wfm, hittime ); 
      int x1 = (int)hittime - (int)fAvePulse_preWindow;
      int x2 = (int)hittime + (int)fAvePulse_postWindow;

      if( x1 >= 0 && x2 < (int)PMT_wfm.size() && a >= fSavedWaveformADCThresh )
        AddToAverageWfm( PMT_wfm, x1, x2, h_AveWfm[ch], fNumSavedWaveforms[ich], -1);

    }
    
  }//<-- end first-pass loop over PMTs


}


//#######################################################################
void AveragePMTWfm::endJob()
{
  std::cout 
    << "\n============================================================\n"
    << "Ending AveragePMTWfm\n";
  
  for(size_t i=0; i<fSelectChannels.size(); i++){
    
    size_t ich = iCh[fSelectChannels[i]];
    
    std::cout<<"PMT "<<fSelectChannels[i]<<": "<<fNumActiveEvents[ich]<<" active events, "<<fNumSavedWaveforms[ich]<<" saved waveforms\n";
  
  }
  
  std::cout<< "\n============================================================\n";
  
}

//#######################################################################
void AveragePMTWfm::endRun(art::Run const & r){
}

//#######################################################################
void AveragePMTWfm::endSubRun(art::SubRun const & sr){
}

//########################################################################################
void AveragePMTWfm::AddToAverageWfm( std::vector<float> &wfm, short t1, short t2, TH1F* h, int &counter, int polarity){
  int nbins = h->GetNbinsX();
  short T = t2 - t1;
  if( nbins ==  T ) {
    float xmin = h->GetXaxis()->GetXmin();
    for(short i=0; i<T; i++) h->Fill(xmin+i, polarity*wfm.at(t1+i) );
    counter++;
  } else {
    LOG_VERBATIM("AveragePMTWfm")<<"ERROR: Number of bins ("<<nbins<<") does not match range given ("<<T<<")";
  }
}


DEFINE_ART_MODULE(AveragePMTWfm)
