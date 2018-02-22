////////////////////////////////////////////////////////////////////////
// Class:       OpDetPhelRate
// Module Type: analyzer
// File:        OpDetPhelRate_module.cc
//
// This module looks for single photoelectron (PE) candidates in the
// optical detector waveforms and makes a histogram of their integrals
// for calibration of the detector's single electron response (SER).
//
// Generated at Fri Mar  4 04:40:04 2016 by William Foreman using artmod
// from cetpkgsupport v1_08_06.
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

// ROOT includes
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TRandom2.h"
#include "TTree.h"

//LAriatSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"
#include "LArIATRecoAlg/TriggerFilterAlg.h" 

class OpDetPhelRate;

class OpDetPhelRate : public art::EDAnalyzer {
public:
  explicit OpDetPhelRate(fhicl::ParameterSet const & p);
  OpDetPhelRate(OpDetPhelRate const &) = delete;
  OpDetPhelRate(OpDetPhelRate &&) = delete;
  OpDetPhelRate & operator = (OpDetPhelRate const &) = delete;
  OpDetPhelRate & operator = (OpDetPhelRate &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  double IntegrateGaus(double sigma, double x1, double x2);
  std::vector<float> FitSER(TH1D* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit);
  std::vector<int> SimpleCluster(std::vector<float> hitlist, float window );
  void resetVariables();

private:

  double PrevTimestamp;
  double TimestampOffset;
  double ModTimestamp;
  
  

  // Tunable parameters defined by fcl
  std::vector<size_t> fSelectChannels;
  std::vector<size_t> fOpDetChannels;
  std::vector<float>  fOpDetPolarity;
  std::vector<float>  fWfmAbsRMSCut;
  std::vector<float>  fWfmAmpCut;
  std::vector<float>  fPrePE_RMSFactorCut;
  std::vector<float>  fPulseHitRMSThresh;
  std::vector<float>  fSetRMS;
  std::vector<float>  fPulseHitThreshHigh;
  std::vector<float>  fPedestalMaxWidth;
  std::vector<float>  fMean_set;
  std::vector<float>  fMean_lowerLim;
  std::vector<float>  fMean_upperLim;
  std::vector<float>  fSinglePE;
  std::vector<float>  fSinglePEWidth;
  std::vector<float>  fSinglePE_tolerance;
  std::vector<std::vector<float>>  SERWaveform;
  std::vector<short>  fPreIntWindow;
  std::vector<short>  fPostIntWindow;
  std::vector<bool>   fFixMean;
  std::vector<float>  fPEIntegralThresh;
  std::vector<std::vector<float>> PhelTimes_integralMethod;
  std::vector<std::vector<float>> PhelTimes_threshMethod;
  std::vector<float>              PhelTimes;
  
  // TTree and its branch variables
  TTree*              fTree;
  int                 tRun;
  int                 tSubRun;
  int                 tEvent;
  int                 tChannel;
  float               tSPE;
  float               tSPEWidth;
  float               tTimeInFrame_ns;
  float               tTimeInSpill_sec;
  float               tTime_sec;
  float               tIntegral;
  float               tAmplitude;

  int                 subrunCount;

  bool                fFixMeanToSinglePE;
  bool		      fSaveAvePhelWfm;
  bool                fLookUpSER;
  bool		      fSkipRunIfNotInTable;
  float               fSecondFitChi2Thresh;
  std::vector<float>  fSER_x1;
  std::vector<float>  fSER_x2;
  size_t              fMaxSavedWaveforms;
  size_t              fNSamples;
  bool                fSubtractRunningBaseline;
  size_t              fRunningBaselineLength;
  bool                fVerbose;
  std::string         fDAQModule;
  std::string         fInstanceName;
  float               fMvPerADC;
  float		      fGradientCut;
  short               fBaselineWindowLength;
  short               fT1;
  short               fT2;
  std::vector<std::string> fSelectEventTypes;
  float               fTimestamp_T1;
  float               fTimestamp_T2;
  float               fWindowExtensionThresh; 
  float               fMaxWindowFactor;
  short               fPrePEBaselineWindow;
  short               fThreshPersist;
  float               fThreshPersistFactor;
  short               fDeadTimeMin;
  short               fQuietTimeMin;
  short               SER_bins[10];
  float               Timestamp;
  char		      histName[100];
  char                histTitle[100];
  char		      buffer[200];

  std::map<size_t,size_t> iCh;

  

  // Alg objects
  OpHitBuilderAlg     fOpHitBuilderAlg;
  TriggerFilterAlg    fTriggerFilterAlg;
  std::string         fTriggerUtility;

  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory waveformDir = tfs->mkdir("Waveforms");
  art::TFileDirectory diagnosticDir = tfs->mkdir("Diagnostics");
  
  // Histograms
  TH1I* h_TotalEvents;
  TH1I* h_TotalSubruns;
  TH1D* h_TotalAveRate;
  TH1D* h_TotalRateVsTime;
  TH1D* h_RateVsSubrun[10];
  TH1D* h_RateVsSubrunRebin[10];
  TH1D* h_Timestamps[10];
  TH1D* h_LiveSamples[10];
  TH1D* h_LiveTime[10];
  TH1D* h_Rate[10];
  TH1D* h_RateStrict[10];
  TH1D* h_RatePedestal[10];
  TH1D* h_AveRate[10];
  TH1D* h_SER[10];
  TH2D* h_SERvTime[10];
  TH1D* h_SERWindowSize[10];
  TH1D* h_AvePEWfm[10];
  TH1I* h_AvePECount[10];
  TH1D* h_RawBaseline[10];
  TH1D* h_BaselineRMS[10];
  TH1D* h_PrePhelRMS[10];
  TH1D* h_PhelTime[10];
  TH1D* h_PhelTimeStrict[10];
  TH1D* h_PedestalTime[10];
  TH1I* h_NumberPhelCandidates_threshMethod[10];
  TH1I* h_NumberPhelCandidates_integralMethod[10];
  TH1D* h_dT[10];
  TH1D* h_dT_nn[10];
  TH1D* h_PhelSample[10];
  size_t fNumberActiveEvents[10];
  TH1D* h_Amplitudes[10];
  size_t fNumberSaturatedEvents[10];
  TH1D* h_PEsPerEvent[10];
  TH1D* h_PhelClusterCount[10];
  TH1D* h_PhelClusterCountBothPMTs;
  TH1D* h_PMTCorrelation;

  TH1I* h_TrigMult_cosmic;
  TH1I* h_TrigMult_beam;
  TH1I* h_TrigsPerSubrun_beam;
  TH1I* h_TrigsPerSubrun_cosmic;
  std::vector<int> TrigsPerSubrun_beam;
  std::vector<int> TrigsPerSubrun_cosmic;
  int   NTrigs_beam;
  int   NTrigs_cosmic;
  int	eventNum_subrun;
  double PE1_scaleFactor[10];
  
  TH1D* h_Waveform;
  TH1D* h_Waveform_raw;
  //TH1D* h_WaveformQuiet;
  //TH1D* h_WaveformDead;
  //TH1D* h_WaveformLive;
  TH1D* h_WaveformSER;
  TH1D* h_WaveformSERGood;
  
  // testing things, andrzej
 // TH1*  h_FFT_output;
//   TH1D* h_Waveform_filt;
//   TH1D* h_Waveform_run;
 // TH1D* h_Waveform_FFT_tots;
 // TCanvas* c2;
  // end testing things
  size_t WaveformCount[10];
  size_t LiveSamples[10];
  size_t LiveSamplesInSubrun[10];
  double NPEsInSubrun[10];
  size_t NPEs[10];
  TCanvas* c;
  //TCanvas*  c_EventNumVsTimestamp;
  //TGraph*   g_EventNumVsTimestamp;

  bool skipThisRun = false;
  

};


//#######################################################################
OpDetPhelRate::OpDetPhelRate(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg")),
fTriggerFilterAlg(p.get<fhicl::ParameterSet>("TriggerFilterAlg"))
{
  this                ->reconfigure(p);
  
  std::cout<<"OpDetPhelRate construction...\n";
  
  SERWaveform.resize(fOpDetChannels.size());
  PhelTimes_integralMethod.resize(fOpDetChannels.size());
  PhelTimes_threshMethod.resize(fOpDetChannels.size());
  PhelTimes.reserve(5000);
 
  // Create the channel index map
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    iCh[fOpDetChannels[i]] = i;
    SER_bins[i] = fPreIntWindow[i]+fPostIntWindow[i];
    SERWaveform[i].resize(SER_bins[i]);
    PhelTimes_integralMethod[i].reserve(1000);
    PhelTimes_threshMethod[i].reserve(1000);
    WaveformCount[i] = 0;
    LiveSamples[i] = 0;
    NPEs[i] = 0;
    fNumberSaturatedEvents[i] = 0;
    fNumberActiveEvents[i]    = 0;
  }

  TrigsPerSubrun_beam   .reserve(1000);
  TrigsPerSubrun_cosmic .reserve(1000);

  std::cout<<"Constructed.\n";
}

//#######################################################################
void OpDetPhelRate::analyze(art::Event const & e)
{
  resetVariables();
  int eventnr   = e.id().event();
  int runnr     = e.run();
  int subrunnr  = e.subRun();
  tEvent        = e.id().event();
  tRun          = e.run();
  tSubRun       = e.subRun();

  // Skip subruns with known strangely high rates
  if( runnr==10125 && ( ( subrunnr >= 46)&&(subrunnr <= 50) ) ) return;
  if( runnr==10129 && ( ( subrunnr <= 5 ) ) ) return;
  if( runnr==10134 && ( ( subrunnr >= 100 ) ) ) return;
  if( runnr==10053 && ( ( subrunnr == 114 )||(subrunnr>=35 && subrunnr <=40) ) ) return;
  if( runnr==10024 && ( ( subrunnr >= 76)&&(subrunnr <= 81) ) ) return;
  if( runnr==10052 && ( ( subrunnr == 112) ) ) return;
  if( runnr==10027 && ( ( subrunnr >= 77)&&(subrunnr <= 78))) return;
  if( runnr==10051 && ( ( subrunnr >= 121)&&(subrunnr <= 122))) return;
  if( runnr==10137 && ( ( subrunnr == 86 ) ) ) return;
  if( runnr==10054 && ( ( subrunnr >= 77)&&(subrunnr <= 84))) return;
  if( runnr==11278 && subrunnr <= 3 ) return;
  if( runnr==11281 && ( subrunnr >= 54 && subrunnr <= 55) ) return;
  if( runnr==11527 && subrunnr==1 ) return;
  if( runnr==11529 && (subrunnr==31 || subrunnr==42 || subrunnr==55) ) return;
  if( runnr==11530 && (subrunnr==22 || subrunnr==23 || (subrunnr >= 34 && subrunnr <= 38 ))) return;
  if( runnr==11531 && (subrunnr==18 || subrunnr==49 || subrunnr==59)) return;
  if( runnr==11532 && subrunnr >= 44 ) return; // beam contamination in last 15 subruns of this "no beam" run 
   
  // If this run's SER value does not appear in the 
  // SER database table, then skip it  
  if( fSkipRunIfNotInTable && skipThisRun ) {
    std::cout<<"Skipping this run...\n";
    return;
  }

  if(fVerbose) std::cout<<"\n";
  eventNum_subrun++;
 
  // reset vectors that will hold all our found phel times for later analysis
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    PhelTimes_integralMethod[i].clear();
    PhelTimes_threshMethod[i].clear();
  }
  PhelTimes.clear();

  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);
  if( (size_t)WaveformHandle->size() == 0 ){
    if(fVerbose) std::cout << "No optical detector data found; skipping event.\n";
    return;
  }

  // Grab the first OpDetPulse in the handle just to check whether its
  // the right event type for use in this analysis
  Timestamp   = (float(WaveformHandle->at(0).PMTFrame())*8.)/1.0e09;
  //g_EventNumVsTimestamp   ->SetPoint(g_EventNumVsTimestamp->GetN(),eventNum_subrun,Timestamp);
  // correct for funny bug in timestamp assignment when we collect beyond 17.18 sec
  if( (PrevTimestamp - Timestamp) > 1. ) TimestampOffset += 17.18;
  ModTimestamp = Timestamp + TimestampOffset;
  if(fVerbose) {
    std::cout << "\n"; 
    std::cout << "Run "<<runnr<<"   sr "<<subrunnr<<"   e "<<eventnr<<" ("<<eventNum_subrun<<" within sr)\n";
    std::cout << "Timestamp: "<<Timestamp<<"  (prev: "<<PrevTimestamp<<")   mod: "<<ModTimestamp<<"    offset = "<<TimestampOffset<<"\n";
  }
  PrevTimestamp = Timestamp; 
  
  // Save timestamp (for tree)
  tTimeInSpill_sec  = ModTimestamp;
 
  std::string eventType = fOpHitBuilderAlg.eventType(ModTimestamp); 
  if( eventType == "cosmic" ) NTrigs_cosmic++;
  if( eventType == "beam"   ) NTrigs_beam++;
  if (  !fOpHitBuilderAlg.eventTypeFilter(ModTimestamp,fSelectEventTypes)
      ||  !(ModTimestamp >= fTimestamp_T1 || ModTimestamp <= fTimestamp_T2)     ){
    if(fVerbose) std::cout<<"Timestamp ("<<ModTimestamp<<" sec) out of range --> skipping event\n";
    return;
  }
   
  
  // Finally, find the amplitudes of each channel.  If any amplitude exceeds
  // the cut defined in the fhicl, OR if the baseline is too noisy, skip event.
  bool skipThisEvent = false;
  for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){
    
    // Get the OpDetPulse from the handle
    art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
    raw::OpDetPulse ThePulse = *ThePulsePtr;
      
    // Check if this optical channel is among those selected
    // for SER analysis in the fhicl script
    size_t ch   = ThePulse.OpChannel();
    size_t ich  = iCh[ch];
    bool gotPMT = false;
    for( size_t i=0; i<fSelectChannels.size(); i++){
      if( fSelectChannels[i] == ch ) { gotPMT=true; break;}
    }
    if(!gotPMT) continue;
    
    // Convert this to a vector<float> for easier manipulation later on
    std::vector<short>  Wfm_raw   = ThePulse.Waveform();
    std::vector<float> Wfm( Wfm_raw.begin(), Wfm_raw.end() );
    size_t              NSamples  = Wfm_raw.size();

    // Get raw baseline
    std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
    float raw_baseline  = tmp[0];
    float raw_rms       = tmp[1];
    h_BaselineRMS[ich]    ->Fill(raw_rms*fMvPerADC);
    h_RawBaseline[ich]    ->Fill(raw_baseline*fMvPerADC);
    
    if( raw_rms*fMvPerADC >= fWfmAbsRMSCut[ich] ) {
      //std::cout<<"Baseline too noisy on channel "<<ich<<" ("<<rms*fMvPerADC<<" mV)... skipping this event.\n";
      skipThisEvent = true;
    } else {
    
      // Find amplitude
      float amp = -1.;
      for(size_t i=0; i<NSamples; i++) {
        float wfm_corrected = (raw_baseline - Wfm[i])*fMvPerADC;
        if( wfm_corrected > amp  ) amp = wfm_corrected;
      }
      
      // Skip this event if the pulse amplitude exceeds our cut
      if( amp >= fWfmAmpCut[ich]) {
        std::cout<<"Amplitude too large on channel "<<ich<<" ("<<amp<<" mV) -- skipping this event ("<<eventnr<<").\n";
        skipThisEvent = true;
      }  
      
    }
  }  
  if( skipThisEvent ) return;
  
  // Count total number of events containing PMT data that have
  // the right timestamp
  h_TotalEvents->Fill(0);

  // --------------------------------------------------------------------
  // Now do full analysis on each pulse
  for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){

    // Get the OpDetPulse from the handle
    art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
    raw::OpDetPulse ThePulse = *ThePulsePtr;
      
    // Check if this optical channel is among those selected
    // for SER analysis in the fhicl script
    size_t ch   = ThePulse.OpChannel();
    size_t ich  = iCh[ch];
    bool gotPMT = false;
    for( size_t i=0; i<fSelectChannels.size(); i++){
      if( fSelectChannels[i] == ch ) { gotPMT=true; break;}
    }
    if(!gotPMT) continue;

    // Set channel (for tree)
    tChannel          = ch;
    tSPE              = fSinglePE[ich];
    tSPEWidth         = fSinglePEWidth[ich];

    // Reset NPEs counter
    double NPEsPerEvent = 0.;
    double NPEsPerEventStrict = 0.;
    size_t liveSamples = 0;

    // Get the raw waveform
    std::vector<short>  Wfm_raw   = ThePulse.Waveform();
    size_t              NSamples  = Wfm_raw.size();
    
    // Convert this to a vector<float> for easier manipulation later on
    std::vector<float> Wfm( Wfm_raw.begin(), Wfm_raw.end() );

    
    // Subtract running baseline to remove oscillations (useful for random pulser runs where 
    // the signal is expected to remain approximately at baseline).
    if(fSubtractRunningBaseline) fOpHitBuilderAlg.SubtractRunningBaseline(Wfm, Wfm, NSamples, fRunningBaselineLength);
    
    // Find waveform baseline and RMS
    std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 300, fBaselineWindowLength );
    float baseline  = tmp[0];
    float rms       = tmp[1];
    
    if(fVerbose){ 
      std::cout << "OpDet "<<ch<<", index "<<ich<<" Timestamp "<<ModTimestamp<<" sec\n";
      std::cout << "Raw waveform baseline: " << baseline << " +/- " << rms <<" ADC\n";
    }

      
    bool areSavingWaveform = false; 
    if( WaveformCount[ich] < fMaxSavedWaveforms ){
	
      sprintf(histName,"Ch%lu_wfmRaw%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,subrunnr,eventnr);
      sprintf(histTitle,"OpDet%lu (raw): run %i, subrun %i, event %i;Sample [ns];Signal [ADC]",ch,runnr,subrunnr,eventnr); 
      h_Waveform_raw	= waveformDir.make<TH1D>(histName,histTitle,Wfm_raw.size(),0,Wfm_raw.size());
      h_Waveform_raw  ->SetOption("HIST");
        
      sprintf(histName,"Ch%lu_wfm%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,subrunnr,eventnr);
      sprintf(histTitle,"OpDet%lu: run %i, subrun %i, event %i;Sample [ns];Signal [mV]",ch,runnr,subrunnr,eventnr); 
      h_Waveform      = waveformDir.make<TH1D>(histName,histTitle,Wfm.size(),0,Wfm.size());
      h_Waveform      ->SetOption("HIST");
	
      sprintf(histName,"Ch%lu_SEROverlay_wfm%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,subrunnr,eventnr);
      sprintf(histTitle,"OpDet%lu: run %i, subrun %i, event %i;Sample [ns];Signal [mV]",ch,runnr,subrunnr,eventnr); 
      c               = waveformDir.make<TCanvas>(histName,"c",700,500);

      	
      //h_WaveformQuiet   = new TH1D("quiet",histTitle,Wfm.size(),0,Wfm.size());
      //h_WaveformQuiet   ->SetLineColor(kCyan); 
      //h_WaveformDead   = new TH1D("dead",histTitle,Wfm.size(),0,Wfm.size());
      //h_WaveformDead   ->SetLineColor(kViolet); 
      //h_WaveformLive   = new TH1D("livetime",histTitle,Wfm.size(),0,Wfm.size());
      //h_WaveformLive   ->SetLineColor(kYellow); 
      h_WaveformSER   = new TH1D("serwindow",histTitle,Wfm.size(),0,Wfm.size());
      h_WaveformSER   ->SetLineColor(kPink+9); 
      h_WaveformSERGood     = new TH1D("serwindow_good",histTitle,Wfm.size(),0,Wfm.size());
      h_WaveformSERGood     ->SetLineColor(kRed);
      h_WaveformSERGood     ->SetLineWidth(1);
      for(short ii=0; ii< short(Wfm.size()); ii++) {
        //h_WaveformQuiet    ->Fill(ii,-0.01);
        //h_WaveformDead   ->Fill(ii,-0.02);
        //h_WaveformLive    ->Fill(ii,-0.03);
        h_WaveformSERGood ->Fill(ii,0.00);
        h_WaveformSER     ->Fill(ii,0.00);
      } 
        
      areSavingWaveform = true;
    }
      
    // Baseline subtraction and signal inversion.
    std::vector<float> wfm_corrected(NSamples);
    float amp = -1.;
    bool isSaturated = false;
    for(size_t i=0; i<NSamples; i++) {
      wfm_corrected[i] = fOpDetPolarity[ich]*((float)Wfm[i]-baseline);
      if(areSavingWaveform) {
	h_Waveform_raw  ->Fill(i,(float)Wfm_raw[i]);
	h_Waveform	  ->Fill(i,wfm_corrected[i]*fMvPerADC);
      }
      if( wfm_corrected[i] > amp  ) amp = wfm_corrected[i];
      if( Wfm_raw[i] == 0	    ) isSaturated = true;
    }
    h_Amplitudes[ich]	->Fill(amp*fMvPerADC);
    if( amp > 20.*rms ) {
      fNumberActiveEvents[ich]++;
      if( isSaturated )	fNumberSaturatedEvents[ich]++;
    }
     
    // Skip this event if the pulse amplitude exceeds our cut
    if( amp*fMvPerADC >= fWfmAmpCut[ich]) {
      std::cout<<"Amplitude too large on channel "<<ich<<" ("<<amp*fMvPerADC<<" mV) -- skipping this event.\n";
      continue; 
    }

      
    // ------------------------------------------------------------ 
    // Now that all the event-level cuts have been passed, keep a count 
    // of the total events analyzed as well as their timestamps
    h_Timestamps[ich]->Fill(ModTimestamp);
    
     
      
    // Determine range of samples to scan over
    short t1 = std::max((short)ThePulse.FirstSample() + fT1, 100);
    short t2 = std::min((short)ThePulse.FirstSample() + fT2, (int)NSamples-100);
      
    // Declare and initialize variables
    float   integral = 0;
    bool    flag = false;
    int     windowsize = 0;
    int     counter = -(int)fPreIntWindow[ich];
    int     dipcounter = 0;
    short   quietTime = 0;
    short   deadTime = 0;
    float   prePE_baseline = 0;
    float   prePE_rms = 999;
    int     hit_time = -1;
    float   PE_amp  = -1;
   
    
    // Make gradient
    //std::vector<float> grad = fOpHitBuilderAlg.MakeGradient(Wfm);
    
    // Make empty vector in which to store waveform for each PE candidate
    // to be reset after each PE (or after added to SERWaveform)
    std::vector<float> tmp_wfm(SER_bins[ich]);
    int tmp_wfm_i = 0;

    // Now begin looking for single PEs
    if(fVerbose){
    std::cout
    << "OpDet "<<ch<<"\n"
    << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh[ich] 
    << ", threshPersist " << fThreshPersist <<"/"<<fThreshPersistFactor<<") \n";
    }

    
    //bool isGoodPE = false;

    if( fSetRMS[ich] > 0. ) rms = fSetRMS[ich]/fMvPerADC;

    // ------------------------------------------------------------
    // Scan waveform and search for hits within specified time window
    for(short i = t1; i < t2; i++){
        
      if(!flag) prePE_baseline = 0;
      float y               = wfm_corrected[i] - prePE_baseline;
      float y_mV            = y*fMvPerADC;
      bool  IsLive          = ((quietTime >= fQuietTimeMin)&&(deadTime >= fDeadTimeMin));
      bool  IsOverThresh    = (y  >= fPulseHitRMSThresh[ich]*rms);
      bool  IsOverThreshNeg = (y  <= -3.0*rms);
      if(IsOverThreshNeg)   { dipcounter++;}
      else{                   dipcounter=0;}
      bool  NegativeDipDetected = (dipcounter >= 2);
      bool  IsOverLimit     = (y_mV >= fPulseHitThreshHigh[ich]); 
      //bool  IsPECandidate   = ((IsOverThresh) && (!IsOverLimit) && ( fabs(grad[i]) >= fabs(fGradientCut)) );
      bool  IsPECandidate   = ((IsOverThresh) && (!IsOverLimit) );
     
//      std::cout<<i<<"  "<<Wfm.at(i)<<"     "<<y<<" pulseHitThresh x"<<fPulseHitRMSThresh[ich]<<"   rms "<<rms<<"   thresh = "<<fPulseHitRMSThresh[ich]*rms<<"   isLive "<<IsLive<<"\n";
      
      if(IsLive || flag) {
        //h_WaveformLive->SetBinContent(i,0.6*fPulseHitRMSThresh[ich]*rms*fMvPerADC);
        liveSamples++;
      }
        
      //printf("%5u %7.3f ADC \n",i,y);

      // If we're already in a PE window, increment the
      // counters and add to integral
      if( flag ) {
        
        if(fVerbose) {
        //std::cout
        //<< "  " << i << "  y_mV = " << y_mV << " mV (window size " <<windowsize<< "), "
        //<< " thresh " << fPulseHitRMSThresh[ich]*rms*fMvPerADC << " mV, g " << grad[i]<<", quietTime = "<<quietTime<<"   deadTime = "<<deadTime<<"\n";
        }

        counter++;
        windowsize++;
        integral += y;
       
        // Scan for amplitude
        if( y_mV > PE_amp ) PE_amp = y;

        // Store signal values
        if(tmp_wfm_i < SER_bins[ich]){
          tmp_wfm[tmp_wfm_i] = y;
          tmp_wfm_i++;
        }
          
        // If (a) a "dip" below baseline was detected, (b) a pulse extends above upper 
        // limit, or (c) if window length is exceeded due to multiple merges, abort mission.
        if( NegativeDipDetected || IsOverThreshNeg || IsOverLimit || (windowsize > fMaxWindowFactor*SER_bins[ich])){
          if(fVerbose) std::cout << "  abort!\n";

          // avoid infinite loops by setting 'i' a little bit ahead
          // of the hit that triggered the start of this integration.
          if( i < hit_time + 3) i = hit_time + 3;
            
          counter = -(int)fPreIntWindow[ich];
          hit_time = -1;
          PE_amp = -1;
          integral = 0.;
          flag = false;
          windowsize = 0;
          tmp_wfm.clear();
          tmp_wfm_i = 0;
          continue;
        }
          
	//----------------------------------
	// End integration if we've reached the minimum required window size and
	// the signal has returned to baseline for at least 2 consecutive samples
        if( counter >= (int)fPostIntWindow[ich] && deadTime >= 2 ){
            
	  if(fVerbose) printf("  * %f ADCs, windowsize %d \n",integral,windowsize);
          
          h_SER[ich]          ->Fill(integral);
          h_SERvTime[ich]     ->Fill(ModTimestamp,integral);
          h_SERWindowSize[ich]->Fill(windowsize);
            
          // Fill tree!
          tTimeInFrame_ns = (float)hit_time;
          tTime_sec       = (subrunCount-1) + tTimeInSpill_sec;
          tIntegral       = integral;
          tAmplitude      = PE_amp; 
          fTree           ->Fill(); 
          

          // if this is pedestal-like, add to histogram
          if(integral < fSinglePE[ich] - fSinglePEWidth[ich] ) h_PedestalTime[ich]->Fill(ModTimestamp);

          // If this exceeds cut, plot the time:
          if( integral > fPEIntegralThresh[ich] ) {
          

            if( fabs(integral-fSinglePE[ich]) < 0.5 * fSinglePEWidth[ich] ){  
              //isGoodPE = true;
              //PhelTimes_integralMethod[ich].push_back(hit_time);
	      h_PhelSample[ich] ->Fill(hit_time);
              NPEsPerEventStrict++;
                
              // Mark on waveform
              if(areSavingWaveform){
                short start = hit_time-fPreIntWindow[ich];
                short end   = start+windowsize;
                for(short ii=start; ii<end; ii++) h_WaveformSERGood ->SetBinContent(ii+1,fPulseHitRMSThresh[ich]*rms*fMvPerADC);
              }
            
            }

            // Divide up into 1pe, 2pe, 3pe
            if( (integral >= fSinglePE[ich] - fSinglePEWidth[ich])&&(integral < 1.5*fSinglePE[ich]) ) {
              PhelTimes_integralMethod[ich].push_back(hit_time);
              NPEs[ich]     += 1;
              NPEsPerEvent  += 1;
            } else 
            if( (integral >= 1.5*fSinglePE[ich])&&(integral < 2.5*fSinglePE[ich]) ) {
              PhelTimes_integralMethod[ich].push_back(hit_time);
              PhelTimes_integralMethod[ich].push_back(hit_time);
              NPEs[ich]     += 2;
              NPEsPerEvent  += 2;
            } else 
            if( (integral >= 2.5*fSinglePE[ich])&&(integral < 3.5*fSinglePE[ich]) ) {
              PhelTimes_integralMethod[ich].push_back(hit_time);
              PhelTimes_integralMethod[ich].push_back(hit_time);
              PhelTimes_integralMethod[ich].push_back(hit_time);
              NPEs[ich]     += 3;
              NPEsPerEvent  += 3;
            }

          }
            
           
          // Add to average waveform if it looks good
          if( fSaveAvePhelWfm && windowsize == SER_bins[ich] && fabs(integral - fSinglePE[ich])/fSinglePE[ich] <= fSinglePE_tolerance[ich] ){
            if(fVerbose) std::cout<<"    --> Adding to average PE waveform.\n";
            for(short ii=0; ii<SER_bins[ich]; ii++)
	      SERWaveform[ich].at(ii) += tmp_wfm[ii]*fMvPerADC;
            h_AvePECount[ich]->Fill(0);
          }

          // Mark on waveform
          if(areSavingWaveform){
            short start = hit_time-fPreIntWindow[ich];
            short end   = start+windowsize;
            for(short ii=start; ii<end; ii++) h_WaveformSER ->SetBinContent(ii+1,fPulseHitRMSThresh[ich]*rms*fMvPerADC);
          }
    
          integral = 0;
          hit_time = -1;
          PE_amp    = -1;
          windowsize = 0;
          counter = -(int)fPreIntWindow[ich];
          flag = false;
          tmp_wfm.clear();
          tmp_wfm_i=0;
          continue;
        }
	//----------------------------------
        
      } // endif flag
   
        
       
      // If we're not yet integrating a PE window, signal is within bounds,
      // and enough quietTime has elapsed, then we're in business
      if( !flag && IsPECandidate && IsLive ){
          
        // Find pre-PE baseline
        prePE_baseline  = 99;
        prePE_rms	  = 99;
        std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS(wfm_corrected,i-fPrePEBaselineWindow,i);
        prePE_baseline  = tmp[0];
        prePE_rms	  = tmp[1];
        h_PrePhelRMS[ich] ->Fill(prePE_rms*fMvPerADC);
        
        // When taking the pre-PE baseline into account, require that the 
        // candidate still passes the threshold:
        y             = wfm_corrected[i] - prePE_baseline;
        y_mV          = y*fMvPerADC;
        IsOverThresh  = (y    >= fPulseHitRMSThresh[ich]*rms);
        IsOverLimit   = (y_mV >= fPulseHitThreshHigh[ich]); 
        IsPECandidate = ((IsOverThresh) && (!IsOverLimit)); 
        
        // If "threshPersist" functionality being used, look ahead and make 
        // sure signal passes the threshold for requisite consecutive samples
        bool flagger = true;
        if( fThreshPersist > 0 ){
          for( short j=1; j <= fThreshPersist; j++ ){
            if ( wfm_corrected[i+j] - prePE_baseline < fThreshPersistFactor*fPulseHitRMSThresh[ich]*rms ) flagger = false;
          }
        }
          
        if( fVerbose ){
          std::cout
          << "  "<<i<<"   "<<y_mV<<" (thresh "<<fPulseHitRMSThresh[ich]*rms*fMvPerADC<<")   wfm rms ="<<rms*fMvPerADC<<"   prePE rms = "<<prePE_rms*fMvPerADC<<"\n";
        }

        // Require flat pre-PE region and that threshold persist requirement was met
        if( (IsPECandidate) && (prePE_rms <= fPrePE_RMSFactorCut[ich]*rms) && (flagger) ){
          
          // Found a "PE hit"
          flag      = true;
          hit_time  = i;
          PhelTimes_threshMethod[ich].push_back(hit_time);
         
          if(fVerbose) printf  ("  candidate at %i, %f mV after %i quietTime / %i deadtime\n",hit_time,y_mV,quietTime,deadTime);
          
          // Go back "preWindow" number of samples
          i -= fPreIntWindow[ich];
        
        } 
    
      } // <-- end if(PE cand)
        
      // "quietTime" will increment for every sample where the signal 
      // is below the limit (default 5mV) and where a sustained negative 
      // dip has not been seen.  If either of these conditions is met, 
      // the counter is reset.
      if( !IsOverLimit && !NegativeDipDetected ){ 
        quietTime++;  
        //h_WaveformQuiet->SetBinContent(i,0.5*fPulseHitRMSThresh[ich]*rms*fMvPerADC);
      } else { 
        quietTime = 0;  
      }

      // "deadTime" is similar, but with stricter conditions.
      if( !IsOverThresh && !IsOverThreshNeg ) { 
        deadTime++; 
        //h_WaveformDead->SetBinContent(i,0.4*fPulseHitRMSThresh[ich]*rms*fMvPerADC);
      } else { 
        deadTime = 0; 
      }


    } // <-- end scan over waveform 
    if(fVerbose) std::cout<<"Ending scan over waveform\n";
      
    if( areSavingWaveform ){
	h_Waveform_raw->SetOption("HIST");
	h_Waveform    ->SetOption("HIST");
        WaveformCount[ich]++;
        c->cd();
        h_Waveform->DrawCopy("HIST");
//	h_Waveform_filt->SetLineColor(kRed);
// 	h_Waveform_filt->DrawCopy("hist same");
// 	h_Waveform_run->SetLineColor(kBlack);
// 	h_Waveform_run->DrawCopy("hist same");
        //h_WaveformSER->Scale(std::max(0.5,h_Waveform->GetMaximum()/10.));
        //h_WaveformQuiet->Draw("hist same");
        //h_WaveformDead->Draw("hist same");
        //h_WaveformLive->Draw("hist same");
        h_WaveformSER->Draw("hist same");
        h_WaveformSERGood->Draw("hist same");
        c->Write(histName);
	
// 	c2->cd();
// 	h_Waveform->Rebin(128);
//         h_Waveform->Draw("hist");
// 	h_Waveform_filt->SetLineColor(kRed);
// 	h_Waveform_filt->Rebin(128);
// 	h_Waveform_filt->Draw("hist same");
//         h_Waveform_run->Rebin(128);
// 	h_Waveform_run->Draw("hist same");
//         //h_WaveformSER->Scale(std::max(0.5,h_Waveform->GetMaximum()/10.));
//         c2->Write(histNameFilt);
	
    }
      

    // Exclude events with abnormal amount of light:
    if( NPEsPerEvent >= 100 ) return;

    h_PEsPerEvent[ich]    ->Fill(NPEsPerEvent);
    for(size_t ii=0; ii<(size_t)NPEsPerEvent; ii++)         h_PhelTime[ich]   ->Fill(ModTimestamp);
    for(size_t ii=0; ii<(size_t)NPEsPerEventStrict; ii++)   h_PhelTimeStrict[ich]   ->Fill(ModTimestamp);
    for(size_t ii=0; ii<(size_t)liveSamples; ii++)          h_LiveSamples[ich]->Fill(ModTimestamp);
    LiveSamples[ich]         += liveSamples;
    NPEsInSubrun[ich]        += NPEsPerEvent;
    LiveSamplesInSubrun[ich] += liveSamples; 
      
    if(fVerbose) std::cout<<"Live samples for this event... "<<liveSamples<<"\n"; 
    
    h_NumberPhelCandidates_threshMethod[ich]->Fill(PhelTimes_threshMethod[ich].size());
    h_NumberPhelCandidates_integralMethod[ich]->Fill(PhelTimes_integralMethod[ich].size());
    

    // analyze for clusters
    std::vector<int> clusterCounts = SimpleCluster(PhelTimes_integralMethod[ich], 5000);
    for(size_t i=0; i<clusterCounts.size(); i++) h_PhelClusterCount[ich]->Fill(clusterCounts.at(i));

    // look at hit-time difference distribution for this PMT
    if( fVerbose) std::cout<<"Looping through found phel times...\n";
    for(size_t i=0; i<PhelTimes_integralMethod[ich].size(); i++){
      //std::cout<<"  "<<i<<": "<<ti<<"\n";
      //int tprev = PhelTimes_integralMethod[ich].at(i-1); 
      for(size_t j=i; j<PhelTimes_integralMethod[ich].size(); j++){
        if( i == j ) continue;
        float dt = fabs( PhelTimes_integralMethod[ich].at(i) - PhelTimes_integralMethod[ich].at(j) ); 
        if( dt > 250. ) h_dT[ich]->Fill(dt); 
      }
      // Nearest-neighbor dT
      if( i>0 ) h_dT_nn[ich]->Fill( PhelTimes_integralMethod[ich].at(i) - PhelTimes_integralMethod[ich].at(i-1) );
    }

    
  } // endLoop over OpDets
    
  // Combine both PMT phel time lists
  PhelTimes.insert( PhelTimes.end(), PhelTimes_integralMethod[0].begin(), PhelTimes_integralMethod[0].end() );
  PhelTimes.insert( PhelTimes.end(), PhelTimes_integralMethod[1].begin(), PhelTimes_integralMethod[1].end() );
  std::sort( PhelTimes.begin(), PhelTimes.end() );
  
  // look for clusters in combined list   
  std::vector<int> clusterCounts = SimpleCluster(PhelTimes, 5000);
  for(size_t i=0; i<clusterCounts.size(); i++) h_PhelClusterCountBothPMTs->Fill(clusterCounts.at(i)); 

  // Calculate time correlations both between the two PMTs
  for(size_t i=0; i<PhelTimes_integralMethod[0].size(); i++){
    for(size_t j=0; j<PhelTimes_integralMethod[1].size(); j++){
      h_PMTCorrelation    ->Fill(fabs( PhelTimes_integralMethod[0].at(i) - PhelTimes_integralMethod[1].at(j) ) ); 
    }
  }


}

//#######################################################################
void OpDetPhelRate::beginJob()
{
  subrunCount = 0;
  
  // Make the tree
  fTree = tfs->make<TTree>("pheltree","pheltree");
  fTree ->Branch("run",             &tRun);
  fTree ->Branch("subrun",          &tSubRun);
  fTree ->Branch("event",           &tEvent);
  fTree ->Branch("channel",         &tChannel);
  fTree ->Branch("spe",             &tSPE); 
  fTree ->Branch("spewidth",        &tSPEWidth); 
  fTree ->Branch("timeinframe_ns",  &tTimeInFrame_ns);
  fTree ->Branch("timeinspill_sec", &tTimeInSpill_sec);
  fTree ->Branch("time_sec",        &tTime_sec);
  fTree ->Branch("integral",        &tIntegral);
  fTree ->Branch("amplitude",       &tAmplitude);

  size_t  spillCycleBins = 22;
  float	  maxTime = 55.;
  
  PrevTimestamp = 0;
  TimestampOffset = 0;
  ModTimestamp = 0;

  h_TotalEvents           = tfs->make<TH1I>("TotalEvents","Total events",1,0,1); 
  h_TotalSubruns          = tfs->make<TH1I>("TotalSubruns","Total subruns",1,0,1);
  h_TotalAveRate	  = tfs->make<TH1D>("TotalAveRate","Total PE rate (combined);;Single photoelectron rate [kHz]",1,0,1);
  h_TotalAveRate	  ->SetOption("E1");
  h_TotalRateVsTime	  = tfs->make<TH1D>("TotalRateVsTime","Total PE rate;Time in spill cycle [sec];Single photoelectron rate [kHz]",spillCycleBins,0,maxTime);
  h_TotalRateVsTime	  ->SetOption("E1");
  h_TotalRateVsTime	  ->Sumw2(0);
  h_TotalRateVsTime	  ->Sumw2(1);
  h_TrigMult_cosmic       = diagnosticDir.make<TH1I>("TrigMult_cosmic","Distribution of cosmic triggers per subruns",150,0,150);
  h_TrigMult_beam         = diagnosticDir.make<TH1I>("TrigMult_beam","Distribution of beam triggers per subruns",150,0,150);

  //g_EventNumVsTimestamp	  = new TGraph;

  h_PMTCorrelation = tfs->make<TH1D>("PMTCorrelation","Phel candidate time differences between two PMTs;#DeltaT [ns]",112,0.,28000.);
  h_PMTCorrelation->Sumw2(); 
    
  h_PhelClusterCountBothPMTs = tfs->make<TH1D>("PhelClusterCountBothPMTs","Number phel candidates within cluster",100,0,100);

  //h_Waveform_FFT_tots     = tfs->make<TH1I>("FFT spectrum","FFT spectrum",1,0,1);
  //h_NTrigs          = tfs->make<TH1I>("NTrigs",  "Trigger bits fired during run",  32,0,32);
  
  // Create histograms
  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ch   = fSelectChannels[i];
    size_t ich  = iCh[ch]; 
  

    // Define SER scale factor for PE counting
    float offset = fPEIntegralThresh[ich] - fSinglePE[ich]; 
    float integral_full = IntegrateGaus(fSinglePEWidth[ich],-10.*fSinglePEWidth[ich],10.*fSinglePEWidth[ich]);
    float integral_part = IntegrateGaus(fSinglePEWidth[ich],offset, 10.*fSinglePEWidth[ich]); 
    PE1_scaleFactor[ich] = integral_full / integral_part;

    int x_bins = std::min((int)(fSER_x2[ich] - fSER_x1[ich]),500);
    
    sprintf(histName,"%lu_SER",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integrated ADC;Counts",ch);
    h_SER[ich]		  = tfs->make<TH1D>(histName,histTitle,x_bins,fSER_x1[ich],fSER_x2[ich]);
   
    sprintf(histName,"%lu_Timestamps",ch); 
    sprintf(histTitle,"Modified timestamp in spill [sec] for OpDet %lu",ch);
    h_Timestamps[ich]     = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    
    sprintf(histName,"%lu_LiveSamples",ch); 
    sprintf(histTitle,"SER live samples for OpDet %lu;Time in spill cycle [sec];Cumulative live samples",ch);
    h_LiveSamples[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_LiveSamples[ich]      ->SetOption("E1");
    
    sprintf(histName,"%lu_LiveTime",ch); 
    sprintf(histTitle,"SER livetime for OpDet %lu;Time in spill cycle [sec];Cumulative livetime [sec]",ch);
    h_LiveTime[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_LiveTime[ich]      ->SetOption("E1");
    
    sprintf(histName,"%lu_PhelTime",ch); 
    //sprintf(histTitle,"PE candidates (> %4.1f) for OpDet %lu;Time in spill cycle [sec];Number of candidates",fPEIntegralThresh[ich],ch);
    sprintf(histTitle,"PE candidates for OpDet %lu;Time in spill cycle [sec];Number of candidates",ch);
    h_PhelTime[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_PhelTime[ich]      ->SetOption("E1");
    
    sprintf(histName,"%lu_PhelTimeStrict",ch); 
    //sprintf(histTitle,"PE candidates (> %4.1f) for OpDet %lu;Time in spill cycle [sec];Number of candidates",fPEIntegralThresh[ich],ch);
    sprintf(histTitle,"PE candidates for OpDet %lu (strict selection);Time in spill cycle [sec];Number of candidates",ch);
    h_PhelTimeStrict[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_PhelTimeStrict[ich]      ->SetOption("E1");
    
    sprintf(histName,"%lu_PedestalTime",ch); 
    //sprintf(histTitle,"PE candidates (> %4.1f) for OpDet %lu;Time in spill cycle [sec];Number of candidates",fPEIntegralThresh[ich],ch);
    sprintf(histTitle,"Pedestal candidates for OpDet %lu;Time in spill cycle [sec];Number of candidates",ch);
    h_PedestalTime[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_PedestalTime[ich]      ->SetOption("E1");
    
    
    sprintf(histName,"%lu_AveRate",ch); 
    sprintf(histTitle,"PE rate for OpDet %lu averaged over spill cycle;;Single photoelectron rate [kHz]",ch);
    h_AveRate[ich]	  = diagnosticDir.make<TH1D>(histName,histTitle,1,0,1);
    
    sprintf(histName,"%lu_PhelSample",ch); 
    sprintf(histTitle,"PE candidates (> %4.1f) for OpDet %lu;Waveform sample;Number of candidates",fPEIntegralThresh[ich],ch);
    h_PhelSample[ich]	  = diagnosticDir.make<TH1D>(histName,histTitle,1000,0.,(float)fNSamples);
    
    sprintf(histName,"%lu_SERvTime",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integrated ADC;Timestamp",ch);
    h_SERvTime[ich]       = diagnosticDir.make<TH2D>(histName,histTitle,spillCycleBins,0.,maxTime,x_bins,fSER_x1[ich],fSER_x2[ich]);
    h_SERvTime[ich]       ->SetOption("colz");
    
    sprintf(histName,"%lu_SERWindowSize",ch); 
    sprintf(histTitle,"SER window sizes for OpDet %lu;Window size (samples);Counts",ch);
    h_SERWindowSize[ich]       = diagnosticDir.make<TH1D>(histName,histTitle,100,0.,fMaxWindowFactor*SER_bins[ich]+10);
  
    if(fSaveAvePhelWfm){ 
    sprintf(histName,"%lu_AvePEWfm",ch); 
    sprintf(histTitle,"Average PE waveform for OpDet %lu;ns;mV",ch);
    h_AvePEWfm[ich]  = tfs->make<TH1D>(histName,histTitle,SER_bins[ich],0.,(float)SER_bins[ich]);
    sprintf(histName,"%lu_AvePECount",ch); 
    sprintf(histTitle,"Total numbwer of averaged single-PE waveforms for OpDet %lu",ch);
    h_AvePECount[ich]  = diagnosticDir.make<TH1I>(histName,histTitle,1,0,1);
    }
    
    sprintf(histName,"%lu_RawBaseline",ch); 
    sprintf(histTitle,"Raw baselines for OpDet %lu;Baseline [ADC];Counts",ch);
    h_RawBaseline[ich]    = diagnosticDir.make<TH1D>(histName,histTitle,600,0.0,1200.);
    
    sprintf(histName,"%lu_BaselineRMS",ch); 
    sprintf(histTitle,"RMS for OpDet %lu;Baseline RMS [mV];Counts",ch);
    h_BaselineRMS[ich]    = diagnosticDir.make<TH1D>(histName,histTitle,200,0.0,1.0);
    
    sprintf(histName,"%lu_PrePhelRMS",ch); 
    sprintf(histTitle,"RMS of pre-PE candidate baselines, OpDet %lu;Pre-PE RMS [mV];Counts",ch);
    h_PrePhelRMS[ich]    = diagnosticDir.make<TH1D>(histName,histTitle,200,0.0,1.0);
    
    //sprintf(histName,"%lu_WfmRMS",ch); 
    //sprintf(histTitle,"Waveform RMS for OpDet %lu;Waveform RMS [mV];Counts",ch);
    //h_WfmRMS[ich]    = tfs->make<TH1D>(histName,histTitle,200,0.0,1.0);
    
    sprintf(histName,"%lu_Amplitudes",ch); 
    sprintf(histTitle,"Pulse amplitudes for OpDet %lu;Amplitude [mV];Counts",ch);
    h_Amplitudes[ich]    = diagnosticDir.make<TH1D>(histName,histTitle,2000.,0.,200.);
    
    sprintf(histName,"%lu_Rate",ch); 
    sprintf(histTitle,"PE rate for OpDet %lu;Time in spill cycle [sec];Single photoelectron rate [kHz]",ch);
    h_Rate[ich]		  = tfs->make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_Rate[ich]           ->SetMarkerStyle(7);
    
    sprintf(histName,"%lu_RateStrict",ch); 
    sprintf(histTitle,"PE rate for OpDet %lu (strict selection);Time in spill cycle [sec];Single photoelectron rate [kHz]",ch);
    h_RateStrict[ich]		  = tfs->make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_RateStrict[ich]           ->SetMarkerStyle(7);
    
    sprintf(histName,"%lu_RatePedestal",ch); 
    sprintf(histTitle,"Pedestal count rate for OpDet %lu;Time in spill cycle [sec];Rate [kHz]",ch);
    h_RatePedestal[ich]		  = tfs->make<TH1D>(histName,histTitle,spillCycleBins,0.,maxTime);
    h_RatePedestal[ich]           ->SetMarkerStyle(7);
    
    sprintf(histName,"%lu_RateVsSubrun",ich); 
    sprintf(histTitle,"OpDet %lu;Subrun;Rate [kHz]",ich);
    h_RateVsSubrun[ich]	  = tfs->make<TH1D>(histName,histTitle,200,0,200);
    //h_RateVsSubrun[ich]	  ->SetOption("E");
    
    sprintf(histName,"%lu_RateVsSubrunRebin",ich); 
    sprintf(histTitle,"OpDet %lu;Subrun;Rate [kHz]",ich);
    h_RateVsSubrunRebin[ich]	  = tfs->make<TH1D>(histName,histTitle,200,0,200);
    //h_RateVsSubrunRebin[ich]	  ->SetOption("E");
    
    
    sprintf(histName,"%lu_PEsPerEvent",ch); 
    sprintf(histTitle,"PE candidates per event for opdet%lu;Counts",ch);
    h_PEsPerEvent[ich]           = tfs->make<TH1D>(histName,histTitle,80,0,80);
    
    sprintf(histName,"%lu_NumberPhelCandidatesThresh",ch); 
    sprintf(histTitle,"Number threshold-based PE candidates for OpDet %lu;Number of candidates",ch);
    h_NumberPhelCandidates_threshMethod[ich]       = diagnosticDir.make<TH1I>(histName,histTitle,80,0,80);
    
    sprintf(histName,"%lu_NumberPhelCandidatesIntegral",ch); 
    sprintf(histTitle,"Number integral-based PE candidates for OpDet %lu;Number of candidates",ch);
    h_NumberPhelCandidates_integralMethod[ich]       = diagnosticDir.make<TH1I>(histName,histTitle,80,0,80);
    
    sprintf(histName,"%lu_PhelClusterCount",ch); 
    sprintf(histTitle,"OpDet %lu;Number phel candidates within cluster",ch);
    h_PhelClusterCount[ich] = tfs->make<TH1D>(histName,histTitle,50,0,50);
   
    sprintf(histName,"%lu_dT",ch); 
    sprintf(histTitle,"Distribution of time differences of phel candidates (integral method) for OpDet %lu;#DeltaT [ns]",ch);
    h_dT[ich]                                         = tfs->make<TH1D>(histName,histTitle,112,0.,28000.);
    h_dT[ich]->SetOption("E0 X0");
    h_dT[ich]->SetMarkerStyle(7);
    
    sprintf(histName,"%lu_dT_nn",ch); 
    sprintf(histTitle,"Distribution of nearest-neighbor time differences of phel candidates for OpDet %lu;#DeltaT [ns]",ch);
    h_dT_nn[ich]                                         = tfs->make<TH1D>(histName,histTitle,112,0.,28000.);
    h_dT_nn[ich]->SetOption("E0 X0");
    h_dT_nn[ich]->SetMarkerStyle(7);
    

  }
}

void OpDetPhelRate::reconfigure(fhicl::ParameterSet const & p)
{
  std::vector<size_t>     SelectChannelsDefault{0, 1};
  std::vector<size_t>     BaselineChannelDefault{};
  std::vector<float>      MeanSetDefault{60, 60};
  std::vector<float>      MeanLowerLimDefault{30, 30};
  std::vector<float>      MeanUpperLimDefault{90, 90};
  std::vector<float>      PrePeRMSFactorCutDefault{1.5, 1.5};
  std::vector<float>      PulseHitThreshHighDefault{8., 8.};
  std::vector<float>      PulseHitRMSThreshDefault{ 4., 4.};
  std::vector<float>      SinglePEDefault{ 60., 60.};
  std::vector<float>      SinglePEWidthDefault{ 15., 15.};
  std::vector<float>      SinglePEToleranceDefault{ 0.1, 0.1};
  std::vector<float>      WfmAbsRMSCutDefault{ 0.5, 0.5};
  std::vector<float>      PedestalMeanLowerLimDefault{ -10, -10};
  std::vector<float>      PedestalMeanUpperLimDefault{ 10, 10};
  std::vector<float>      x1Default{0,0};
  std::vector<float>      x2Default{350,350};
  std::vector<std::string> EventTypeDefaults{"beam"};
  std::vector<short>     PreWindowDefault{5,5};
  std::vector<short>     PostWindowDefault{20,20};
  std::vector<float>      PolarityDefault{-1.,-1.};
  std::vector<float>      PedestalMaxWidthDefault{30.,30.};
  std::vector<float>	  PEIntegralThreshDefault{30.,30.};
  std::vector<float>	  WfmAmpCutDefault{15,15};
  std::vector<float>      SetRMSDefault{0.,0.};

  fSER_x1                 = p.get< std::vector<float> >   ("SER_x1",x1Default);
  fSER_x2                 = p.get< std::vector<float> >   ("SER_x2",x2Default);
  fSelectChannels         = p.get< std::vector<size_t> >  ("SelectChannels",SelectChannelsDefault);
  fOpDetPolarity        = p.get< std::vector<float> >   ("OpDetPolarity",PolarityDefault);
  fOpDetChannels          = p.get< std::vector<size_t> >  ("OpDetChannels", SelectChannelsDefault);
  fMean_set               = p.get< std::vector<float> >   ("Mean_set",MeanSetDefault);
  fMean_lowerLim          = p.get< std::vector<float> >   ("Mean_lowerLim",MeanLowerLimDefault);
  fMean_upperLim          = p.get< std::vector<float> >   ("Mean_upperLim",MeanUpperLimDefault);
  fPrePE_RMSFactorCut     = p.get< std::vector<float> >   ("PrePE_RMSFactorCut",PrePeRMSFactorCutDefault);
  fGradientCut            = p.get< float >                ("GradientCut",0);
  fPulseHitThreshHigh     = p.get< std::vector<float> >   ("PulseHitThresh_high",PulseHitThreshHighDefault);
  fSetRMS                 = p.get< std::vector<float> >   ("SetRMS",SetRMSDefault);
  fPulseHitRMSThresh      = p.get< std::vector<float> >   ("PulseHitRMSThresh",PulseHitRMSThreshDefault);
  fSinglePE               = p.get< std::vector<float> >   ("SinglePE",SinglePEDefault);
  fSinglePEWidth          = p.get< std::vector<float> >   ("SinglePEWidth",SinglePEWidthDefault);
  fSinglePE_tolerance     = p.get< std::vector<float> >   ("SinglePE_tolerance",SinglePEToleranceDefault);
  fWfmAbsRMSCut           = p.get< std::vector<float> >   ("WfmAbsRMSCut", WfmAbsRMSCutDefault);
  fWfmAmpCut		  = p.get< std::vector<float> >   ("WfmAmpCut", WfmAmpCutDefault);
  fPedestalMaxWidth       = p.get< std::vector<float> >   ("PedestalMaxWidth",PedestalMaxWidthDefault);
  fPreIntWindow              = p.get< std::vector<short> >   ("PreIntWindow",PreWindowDefault);
  fPostIntWindow             = p.get< std::vector<short> >   ("PostIntWindow",PostWindowDefault);

  fFixMeanToSinglePE      = p.get< bool >         ("FixMeanToSinglePE",false);
  fSaveAvePhelWfm	  = p.get< bool >	  ("SaveAvePhelWfm",true);
  fSkipRunIfNotInTable	  = p.get< bool >	  ("SkipRunIfNotInTable",false);
  fLookUpSER              = p.get< bool >         ("LookUpSER",false);
  fSecondFitChi2Thresh    = p.get< float >        ("SecondFitChi2Thresh",1.2);
  fTriggerUtility         = p.get< std::string >  ("TriggerUtility","FragmentToDigit");
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fTimestamp_T1           = p.get< float >        ("Timestamp_T1",0);
  fTimestamp_T2           = p.get< float >        ("Timestamp_T2",60);
  fWindowExtensionThresh  = p.get< float >        ("WindowExtensionThresh",1.0);
  fMaxWindowFactor        = p.get< float >        ("MaxWindowFactor",2);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.1953);
  fPrePEBaselineWindow    = p.get< short >        ("PrePEBaselineWindow",100);
  fThreshPersist          = p.get< short >        ("ThreshPersist",0);
  fThreshPersistFactor    = p.get< float >        ("ThreshPersistFactor",1.0);
  fDeadTimeMin            = p.get< short >        ("DeadTimeMin",10);
  fQuietTimeMin           = p.get< short >        ("QuietTimeMin",200);
  fVerbose                = p.get< bool >         ("Verbose",true);
  fNSamples               = p.get< size_t >        ("NSamples",28672);
  fSelectEventTypes       = p.get< std::vector<std::string> > ("SelectEventTypes",EventTypeDefaults);
  fMaxSavedWaveforms      = p.get< size_t >       ("MaxSavedWaveforms",10);
  fSubtractRunningBaseline= p.get< bool >         ("SubtractRunningBaseline",false);
  fRunningBaselineLength  = p.get< size_t >         ("RunningBaselineLength",256);
  fPEIntegralThresh	  = p.get< std::vector<float >>	  ("PEIntegralThresh",PEIntegralThreshDefault);
}

void OpDetPhelRate::beginRun(art::Run const & r){
 
  skipThisRun = false; 
  //fTriggerFilterAlg.loadXMLDatabaseTable( r.run() );

    // Look up the SER for photodetectors
    if( fLookUpSER ){

    for(size_t i=0; i<fSelectChannels.size(); i++){
      size_t ch  = fSelectChannels[i];
      size_t ich = iCh[fSelectChannels[i]];
     
      sprintf(buffer,"LArIATPhotodetectorSER_ch%lu.txt",ch);
      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(buffer, fullname); 
   
      if (fullname.empty()) {
        //std::cout << "Input file " << buffer << " not found" << std::endl;
        throw cet::exception("File not found");
      } else {
        //std::cout << "Reading SER value from file = " <<buffer<< std::endl;
        //std::cout << "Run = "<<r.run()<<"\n";

        std::ifstream inFile(fullname, std::ios::in);
        std::string line;

        size_t lineNum = 0;
        size_t runNum = 0;
        bool foundRunInTable = false;
        while (std::getline(inFile,line)) {
          // skip first line (column headers)
          lineNum++;
          if(lineNum == 1 ) continue;
          size_t channel;
          float ser, serErr, serWidth;
          std::istringstream ss(line);
          ss >> runNum >> channel >> ser >> serErr >> serWidth;
          //std::cout<<"   "<<runNum<<"   "<<channel<<"   "<<ser<<"\n";
          if( runNum == (size_t)r.run() ){
            //std::cout<<"Found SER in database table: "<<ser<<"\n";
            fSinglePE[ich]	    = ser;
            fSinglePEWidth[ich] = serWidth;
            fMean_set[ich]      = ser;
            foundRunInTable = true;
            break;
          }
        }
   
        if(!foundRunInTable) {
          std::cout<<"!!!!! RUN NOT FOUND IN SER DATABASE TABLE   !!!!\n";
          skipThisRun = true;
        }
      }
    }
  }

}

void OpDetPhelRate::beginSubRun(art::SubRun const & sr){
  std::cout<<"beginSubRun\n";
  
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    LiveSamplesInSubrun[i] = 0;
    NPEsInSubrun[i] = 0;
  }
 
  subrunCount++; 
  eventNum_subrun = 0;
  NTrigs_cosmic   = 0;
  NTrigs_beam     = 0;
  TimestampOffset = 0;
  PrevTimestamp   = 0;
  ModTimestamp    = 0;
}

void OpDetPhelRate::endJob()
{

  std::cout<<"Ending job...\n";


  /*
  c_EventNumVsTimestamp	  = tfs->make<TCanvas>("EventNumVsTimestamp","c",700,500);
  c_EventNumVsTimestamp	  ->cd();
  g_EventNumVsTimestamp	  ->Draw("AP"); 
  g_EventNumVsTimestamp	  ->GetXaxis()->SetTitle("Event number (within subrun)");
  g_EventNumVsTimestamp	  ->GetYaxis()->SetTitle("Unmodified timestamp [sec]");
  c_EventNumVsTimestamp	  ->Update();
  c_EventNumVsTimestamp	  ->Write("EventNumVsTimestamp");
  */

  // Make histograms showing NTrigs for each subrun analyzed
  // (particularly useful if input list is sorted by creation time)
  //h_TrigsPerSubrun_beam   = tfs->make<TH1I>("TrigsPerSubrun_beam","Beam triggers per subrun",
  //  TrigsPerSubrun_beam.size(),0,TrigsPerSubrun_beam.size());
  //h_TrigsPerSubrun_cosmic = tfs->make<TH1I>("TrigsPerSubrun_cosmic","Cosmic triggers per subrun",
  //  TrigsPerSubrun_cosmic.size(),0,TrigsPerSubrun_cosmic.size());

  //for( size_t i = 0; i < h_TotalSubruns->GetEntries(); i++){
  //  h_TrigsPerSubrun_beam->Fill(i,TrigsPerSubrun_beam[i]);
  //  h_TrigsPerSubrun_cosmic->Fill(i,TrigsPerSubrun_cosmic[i]);
  //}


  h_PMTCorrelation->Scale(1./(double)h_TotalEvents->GetEntries());
    h_PhelClusterCountBothPMTs->Scale(1./double(h_TotalEvents->GetEntries()));


  // Print out general parameters used
  std::cout 
    << "\n============================================================\n"
    << "Ending SER program.\n";
  
  for(size_t i=0; i<fSelectChannels.size(); i++){
    
    size_t ich = iCh[fSelectChannels[i]];
  
    // Make rebinned/rescaled version of rat
    h_RateVsSubrunRebin[ich]->Add(h_RateVsSubrun[ich]);
    h_RateVsSubrunRebin[ich]->Rebin(5);
    //h_RateVsSubrunRebin[ich]->Sumw2();
    h_RateVsSubrunRebin[ich]->Scale(0.2);
    h_RateVsSubrunRebin[ich]->SetOption("E1");
    h_RateVsSubrun[ich]->SetOption("E1");

    // Calculate the total cumulative time the SER code was "live"
    float cumulativeTime = (float)LiveSamples[ich]*1.0e-9;

    // Scale live-time histogram to seconds and create
    // the rate histograms (vs. time in spill cycle)
    h_LiveTime[ich]   ->Sumw2();
    h_LiveSamples[ich]->Sumw2();
    h_LiveTime[ich]   ->Add(h_LiveSamples[ich]);
    h_LiveTime[ich]   ->Scale(1.0e-9);
    //h_PhelTime[ich]   ->Sumw2();
    h_Rate[ich]	      ->Sumw2(0);
    h_Rate[ich]	      ->Sumw2(1);
    h_Rate[ich]	      ->Add(h_PhelTime[ich]);
    h_Rate[ich]	      ->Divide(h_LiveTime[ich]);
    h_Rate[ich]	      ->Scale(1./1000.);
    h_Rate[ich]	      ->SetOption("E1 X0");
    // Zero-out bins for times < 5sec
    // behavior was wonky here, plus this is the 
    // beam window bin, which makes things confusing.
    h_Rate[ich]       ->SetBinContent(1,0.);
    h_Rate[ich]       ->SetBinError(1,0.);
    h_Rate[ich]       ->SetBinContent(2,0.);
    h_Rate[ich]       ->SetBinError(2,0.);
    
    h_TotalRateVsTime ->Add(h_Rate[ich]);
    h_TotalRateVsTime ->SetOption("E1");

    // "strict" rate (integral near SER peak)
    h_RateStrict[ich]	      ->Sumw2(0);
    h_RateStrict[ich]	      ->Sumw2(1);
    h_RateStrict[ich]	      ->Add(h_PhelTimeStrict[ich]);
    h_RateStrict[ich]	      ->Divide(h_LiveTime[ich]);
    h_RateStrict[ich]	      ->Scale(1./1000.);
    h_RateStrict[ich]	      ->SetOption("E1 X0");
    h_RateStrict[ich]       ->SetBinContent(1,0.);
    h_RateStrict[ich]       ->SetBinError(1,0.);
    h_RateStrict[ich]       ->SetBinContent(2,0.);
    h_RateStrict[ich]       ->SetBinError(2,0.);
    
    // "pedestal" rate
    h_RatePedestal[ich]	      ->Sumw2(0);
    h_RatePedestal[ich]	      ->Sumw2(1);
    h_RatePedestal[ich]	      ->Add(h_PedestalTime[ich]);
    h_RatePedestal[ich]	      ->Divide(h_LiveTime[ich]);
    h_RatePedestal[ich]	      ->Scale(1./1000.);
    h_RatePedestal[ich]	      ->SetOption("E1 X0");
    h_RatePedestal[ich]       ->SetBinContent(1,0.);
    h_RatePedestal[ich]       ->SetBinError(1,0.);
    h_RatePedestal[ich]       ->SetBinContent(2,0.);
    h_RatePedestal[ich]       ->SetBinError(2,0.);

   
    // Get light-like ratio
    double tt = 250. + 0.5*(28000. - 250.)*(2. - sqrt(2.));
    double b0 = h_dT[ich]->GetXaxis()->FindBin(0.);
    double b1 = h_dT[ich]->GetXaxis()->FindBin(tt);
    double b2 = h_dT[ich]->GetXaxis()->FindBin(28000.);
    double NA = h_dT[ich]->Integral(b0,b1);
    double NB = h_dT[ich]->Integral(b1,b2);
    double rr = NA / (NA + NB);
    double Pl = 2.*rr - 1.;
    std::cout<<"NA ( "<<b0<<" - "<<b1<<" ) = "<<NA<<"\n";
    std::cout<<"NB ( "<<b1<<" - "<<b2<<" ) = "<<NB<<"\n";
    std::cout<<"NA / total = "<<NA/(NA+NB)<<"\n";
    

    std::cout<<"Doing same for nearest neighbor dT...\n";
    NA = h_dT_nn[ich]->Integral(b0,b1);
    NB = h_dT_nn[ich]->Integral(b1,b2);
    double rrnn = NA / (NA + NB);
    double Plnn = 2.*rrnn - 1.;
    std::cout<<"NA ( "<<b0<<" - "<<b1<<" ) = "<<NA<<"\n";
    std::cout<<"NB ( "<<b1<<" - "<<b2<<" ) = "<<NB<<"\n";
    std::cout<<"NA / total = "<<NA/(NA+NB)<<"\n";
   
   
    h_dT[ich]->Sumw2(); 
    h_dT_nn[ich]->Sumw2(); 
    h_dT_nn[ich]   ->Scale(1./double(h_TotalEvents->GetEntries())); 
    h_dT[ich]   ->Scale(1./double(h_TotalEvents->GetEntries())); 
    h_PhelClusterCount[ich]->Scale(1./double(h_TotalEvents->GetEntries()));

    std::vector<float> params(17);
    params = FitSER(h_SER[ich],fSER_x1[ich], fSER_x2[ich], fMean_set[ich], fSinglePEWidth[ich], fPedestalMaxWidth[ich], true);

    float n0  = params[0];
    float dn0 = params[1];
    float mu0 = params[2];
    float dmu0 = params[3];
    float sig0 = params[4];
    float dsig0 = params[5];
    float n1  = params[6];
    float dn1 = params[7];
    float mu  = params[8];
    float dmu = params[9];
    float sig   = params[10];
    float dsig  = params[11];
    float n2   = params[12];
    float dn2  = params[13];
    float n3    = params[14];
    float dn3    = params[15];
    float chi  = params[16];
    
    double rate = 0., drate = 0.;
    double crate = 0., dcrate = 0.;
    double totalphels = n1 + 2.*n2 + 3.*n3;
    double dtotalphels = sqrt(pow(dn1,2) + pow(2.*dn2,2) + pow(3.*dn3,2));
    if( cumulativeTime > 0.) {
      rate = totalphels/cumulativeTime;
      drate= dtotalphels/cumulativeTime;
      crate = NPEs[ich] / cumulativeTime;
      dcrate = sqrt(double(NPEs[ich])) / cumulativeTime;
    }

    // Fit each of the limited-range SER plots
    for(int ibin = 0; ibin < h_SERvTime[ich]->GetNbinsX(); ibin++ ){
      params = FitSER(h_SERvTime[ich]->ProjectionY("_py",ibin,ibin),fSER_x1[ich], fSER_x2[ich], fMean_set[ich], fSinglePEWidth[ich], fPedestalMaxWidth[ich],false);
      h_PedestalTime[ich]->SetBinContent(ibin,params[0]);
      h_PedestalTime[ich]->SetBinError(ibin,params[1]);
    }
    
    //==========================================================================
    // 2016-12-12, Fit dT histograms and estimate rate based on parameterization.
    double p0   =  0.00332;
    double ep0  =  0.00067;
    double p1   = 18.900;
    double ep1  =  0.001; 

    sprintf(buffer,"sqrt( ((%f*y)/(2.*sqrt(x)))^2 + (%f)^2 + x*(%f)^2)",p1,ep0,ep1);
    TF2 ferr("ferr",buffer);

    double binwidth = h_dT[ich]->GetXaxis()->GetBinWidth(1);
    //sprintf(buffer,"[0]*(1.-x/28000.)*(1e-3)*(%f/100.) + [1]*exp(-x/[2])",binwidth);
    sprintf(buffer,"[0]*(1.-x/28000.)*(1e-3)*(%f/100.) + [1]*exp(-x/[2]) + [3]*exp(-x/[4])",binwidth);
    TF1 f_randrate("f_randrate",buffer,500.,28000.);
    f_randrate.SetParameter(0,1.);
    f_randrate.SetParameter(1,h_dT[ich]->GetMaximum());
    f_randrate.SetParLimits(1,0.,1.);
    f_randrate.SetParameter(2,1500.);
    f_randrate.SetParLimits(2,1200.,1800.);
    f_randrate.SetParameter(3,h_dT[ich]->GetMaximum()*0.1);
    f_randrate.SetParLimits(3,0.,1.);
    f_randrate.SetParameter(4,5000.);
    f_randrate.SetParLimits(4,3000.,20000.);
    h_dT[ich]->Fit("f_randrate","QR");
    std::cout<<"Resulting slope: "<<f_randrate.GetParameter(0)<<"\n";
    double randRate   = (1e3)*(p0 + p1*sqrt(f_randrate.GetParameter(0)));
    double drandRate  = (1e3)*ferr.Eval(f_randrate.GetParameter(0),f_randrate.GetParError(0));
    double lightrate  = (rate - randRate);
    double dlightrate = sqrt(pow(drate,2) + pow(drandRate,2));

    h_AveRate[ich]->Fill(0.,rate/1000.);
    h_AveRate[ich]->SetBinError(1,drate/1000.);
    h_TotalAveRate->Add(h_AveRate[ich]);

    // Print out photodetector-specific parameters
    std::cout
    << "\n-------------------------------------------\n"
    << "Channel	" << fOpDetChannels[ich] << "\n"
    << "  Pre-PE baseline length    " << fPrePEBaselineWindow << "\n"
    << "  Running BS subtraction    " << fSubtractRunningBaseline << " ("<<fRunningBaselineLength<<" samples)\n"
    << "  RMS thresh factor	    " << fPulseHitRMSThresh[ich] <<" x wfm RMS\n"
    << "  Pre-PE RMS cut	    " << fPrePE_RMSFactorCut[ich] << " x wfm RMS\n"
    << "  Threshhold persist	    " << fThreshPersist <<" (x "<<fThreshPersistFactor<<")\n"
    << "  Quiet time min            " << fQuietTimeMin <<"\n"
    << "  Dead time min             " << fDeadTimeMin << "\n"
    << "  Pre / Post window size    " << fPreIntWindow[ich] << "," << fPostIntWindow[ich] << "\n"
    << "  Abs Wfm RMS cut	    " << fWfmAbsRMSCut[ich] << " mV\n"
    << "  Phel integral cut	     >" << fPEIntegralThresh[ich] << " integrated ADC\n"
    << "  Set single PE             " << fSinglePE[ich] << " integrated ADC\n"
    << "  1PE count scale factor    " << PE1_scaleFactor[ich] << "\n"
    << "  Active events		    " << fNumberActiveEvents[ich] << "\n"
    << "  Saturation rate	    " << (float)fNumberSaturatedEvents[ich] / (float)fNumberActiveEvents[ich] << "\n"
    << "  Mean baseline RMS	    " << h_BaselineRMS[ich]->GetMean(1) << " mV\n"
    << "  --> PE candidates found   " << h_SER[ich]->GetEntries() << "\n"
    << "\n";
    
    std::cout 
    <<"*****************************************\n"
    <<"SER fit results for OpDet "<< fOpDetChannels[ich]<<":\n"
    <<"  Noise N         = "<<n0<<" +/- "<<dn0<<"\n"
    <<"  Noise mean      = "<<mu0<<" +/- "<<dmu0<<"\n"
    <<"  Noise sigma     = "<<sig0<<" +/- "<<dsig0<<"\n"
    <<"  1PE mean (SER)  = "<<mu<<" +/- "<<dmu<<"\n"
    <<"  1PE sigma       = "<<sig<<" +/- "<<dsig<<"\n"
    <<"  1PE N           = "<<n1<<" +/- "<<dn1<<"\n"
    <<"  2PE N           = "<<n2<<" +/- "<<dn2<<"\n"
    <<"  3PE N           = "<<n3<<" +/- "<<dn3<<"\n"
    <<"  chisquare / NDF = "<<chi<<"\n"
    <<"  1+2+3 PE intgrl = "<<totalphels<<" +/- "<<dtotalphels<<" PE\n"
    <<"  live time       = "<<cumulativeTime<<" sec\n\n"
    <<"  rate from final SER integral :  "<<rate/1000.<<" +/- "<<drate/1000.<<" kHz\n"
    <<"  rate from counting           :  "<<crate/1000.<<" +/- "<<dcrate/1000.<<" kHz\n"
    <<"    - random-like rate         :  "<<randRate/1000.<<" +/- "<<drandRate/1000.<<" kHz\n"
    <<"    - light-like rate          :  "<<lightrate/1000.<<" +/- "<<dlightrate/1000.<<" kHz\n" 
    <<"    - NA / total               :  "<<rr<<"\n"
    <<"    - light-likeness param     :  "<<Pl<<"\n"
    <<"    - NA / total (NN)          :  "<<rrnn<<"\n"
    <<"    - light-likeness (NN)      :  "<<Plnn<<"\n" 
    <<"*****************************************\n"
    <<"\n";
    


    // Normalize the summed PE waveform to get an average
    if(fSaveAvePhelWfm && h_AvePECount[ich]->GetEntries() > 0 ){
      //std::cout
      //<< "Ave PE waveform for OpDet "<<fOpDetChannels[ich]<<" ("<<h_AvePECount[ich]->GetEntries()<<" entries)\n"
      //<< "Selection criteria: "<<fSinglePE[ich]<<" ADC +/- "<<fSinglePE_tolerance[ich]*100.<<"%\n" 
      //<< "Format:  sample [ns]   ADC\n";
      float integral = 0;
      for( int i = 0; i < SER_bins[ich]; i++) {
        float w = SERWaveform[ich].at(i) / float(h_AvePECount[ich]->GetEntries()); 
        h_AvePEWfm[ich]->Fill(i,w);
        integral += w/fMvPerADC;
        //std::cout<<"  "<<i<<"     "<<w/fMvPerADC<<"\n";
      }
      h_AvePEWfm[ich]->SetOption("HIST");
      //std::cout<< "\n-------------------------------------------\n";
    }

    
    
  
  } // end loop over photodetectors
  
}

void OpDetPhelRate::endRun(art::Run const & r){
    
}


void OpDetPhelRate::endSubRun(art::SubRun const & sr){
  h_TotalSubruns	  ->Fill(0);
  h_TrigMult_beam         ->Fill(NTrigs_beam);
  h_TrigMult_cosmic       ->Fill(NTrigs_cosmic);
 
  // Record rate for this subrun
  for(size_t i=0; i< fSelectChannels.size();i++) {
  
    double totalRate      = 0.;
    double totalRateErr   = 0.;
    double liveTime       = (double)LiveSamplesInSubrun[i]*1.0e-9; 
    
    if( liveTime > 0. ) {
      totalRate	      = NPEsInSubrun[i]/liveTime;
      totalRateErr    = totalRate / sqrt(NPEsInSubrun[i]);
    }
    
    h_RateVsSubrun[i] ->SetBinContent(sr.subRun()+1, totalRate/1000.);
    h_RateVsSubrun[i] ->SetBinError  (sr.subRun()+1, totalRateErr/1000.);
  
   
//    if(fVerbose){
      std::cout<<"Done with SUBRUN "<<sr.subRun()<<"\n";
      std::cout<<"Live time for opdet "<<i<<"   "<<liveTime<<"\n";
      std::cout<<"NPEs "<<NPEsInSubrun[i]<<"\n";
      std::cout<<"Rate = "<<totalRate<<" +/- "<<totalRateErr<<"\n";
      //std::cout<<"PERatePerSubrun[ "<<i<<" ] = "<<PERatePerSubrun[i][PERatePerSubrun.size()-1]<<"\n";
  //  }
  
  }
}

// ============================================================================================
double OpDetPhelRate::IntegrateGaus(double sigma, double  x1, double x2 )
{
  double t1 = x1/sigma;
  double t2 = x2/sigma;
  return 0.5*(sigma*sqrt(2.*TMath::Pi())) * ( TMath::Erf(t1/sqrt(2.)) + TMath::Erf(t2/sqrt(2.)) );
}

// ============================================================================================
std::vector<float> OpDetPhelRate::FitSER(TH1D* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit = true)
{ 

  // Fit out the SER
  TF1 SER_fit("SER_fit","([0]/([2]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-[1])/[2],2))"
        "+ [3]*((1./(sqrt(1.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-1.*[4])/(sqrt(1)*[5]),2))" 
        "+ ([6]/(sqrt(2.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-2.*[4])/(sqrt(2)*[5]),2))" 
        "+ ([7]/(sqrt(3.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-3.*[4])/(sqrt(3)*[5]),2)))",x1,x2);
 
  float max = (float)h_SER->GetEntries();

    // Aug 3 2016: Bug in ROOT6 causes refitting to fail
    // in Fit Panel when parameters are named.
    SER_fit.SetParName(0,"N0");
    SER_fit.SetParName(1,"#mu0");
    SER_fit.SetParName(2,"#sigma0");
    SER_fit.SetParName(3,"N1");
    SER_fit.SetParName(4,"#mu1 (SER)");
    SER_fit.SetParName(5,"#sigma1");
    SER_fit.SetParName(6,"rN2");
    SER_fit.SetParName(7,"rN3");

    // "Noise" component (gaus)
    SER_fit.SetParameter(0,max);
    SER_fit.SetParLimits(0,0.,max);
    SER_fit.SetParameter(1,0);
    SER_fit.SetParLimits(1,-0.50*meanSet,0.50*meanSet);
    SER_fit.SetParameter(2,pedMaxWidth*0.5);
    SER_fit.SetParLimits(2,1.,pedMaxWidth);

    // 1PE (gaus)
    SER_fit.SetParameter(3,max);
    SER_fit.SetParLimits(3,0.,max);
    SER_fit.SetParameter(4,meanSet);
    SER_fit.SetParLimits(4,meanSet*0.8,meanSet*2.);
    SER_fit.SetParameter(5,width);
    SER_fit.SetParLimits(5,5.,width*3.);
    if(fFixMeanToSinglePE) SER_fit.FixParameter(4,meanSet);

    // N2 (ratio)
    SER_fit.SetParameter(6,0.01);
    SER_fit.SetParLimits(6,0.,0.3);

    // N3 (ratio)
    SER_fit.SetParameter(7,0.001);
    SER_fit.SetParLimits(7,0.,0.3);

    // Perform an initial fit.  If Chi2/NDF > 1.2, refit with restricted
    // range that excludes most of the pedestal peak
    h_SER->Fit("SER_fit","QR");
    float n0  = SER_fit.GetParameter(0);
    float dn0 = SER_fit.GetParError(0);
    if( refit && SER_fit.GetChisquare()/(SER_fit.GetNDF()-1) > fSecondFitChi2Thresh ) {
      SER_fit.SetRange(SER_fit.GetParameter(1)+1.5*SER_fit.GetParameter(2), x2);
      h_SER->GetFunction("SER_fit")->Delete();
      h_SER->Fit("SER_fit","QR");
    }

    float mu0 = SER_fit.GetParameter(1);
    float dmu0 = SER_fit.GetParError(1);
    float sig0 = SER_fit.GetParameter(2);
    float dsig0 = SER_fit.GetParError(2);
    float n1  = SER_fit.GetParameter(3);
    float dn1 = SER_fit.GetParError(3);
    float mu  = SER_fit.GetParameter(4);
    float dmu = SER_fit.GetParError(4);
    float sig   = SER_fit.GetParameter(5);
    float dsig  = SER_fit.GetParError(5);
    float n2    = n1*SER_fit.GetParameter(6);
    float dn2   = n1*SER_fit.GetParError(6);
    float n3    = n1*SER_fit.GetParameter(7);
    float dn3   = n1*SER_fit.GetParError(7);
    float chi   = SER_fit.GetChisquare()/(SER_fit.GetNDF()-1);

    std::vector<float> out(17);
    out[0] = n0;
    out[1] = dn0;
    out[2] = mu0;
    out[3] = dmu0;
    out[4] = sig0;
    out[5] = dsig0;
    out[6] = n1;
    out[7] = dn1;
    out[8] = mu;
    out[9] = dmu;
    out[10] = sig;
    out[11] = dsig;
    out[12] = n2;
    out[13] = dn2;
    out[14] = n3;
    out[15] = dn3;
    out[16] = chi;

    return out;

}

// =====================================================================
std::vector<int> OpDetPhelRate::SimpleCluster(std::vector<float> hitlist, float window ) {    

    std::vector<int> out;

    // sort list
    std::sort( hitlist.begin(), hitlist.end());

    int cluster_end = 0;
    int nphelcluster = 0;
    
    // Look at hit-time difference distribution for this PMT
    if( fVerbose) std::cout<<"Looping through found phel times...\n";
    for(size_t i=0; i<hitlist.size(); i++){
      int ti    = hitlist.at(i);

      // Look for clusters
      if( cluster_end > 0 ) {
        if( ti <= cluster_end ) {
          nphelcluster++;
          //std::cout<<"  Adding phel to cluster!\n";
        } else {
          out.push_back(nphelcluster);
          //std::cout<<"  Cluster ended... n = "<<nphelcluster<<"\n";
          nphelcluster = 0;
          cluster_end = 0;
        }
      }
      if( cluster_end == 0 ) {
        cluster_end = ti + window;
        nphelcluster = 1;
        //std::cout<<"  cluster end set to "<<cluster_end<<"\n";
      }   
    }

  return out;

}

// =====================================================================
void OpDetPhelRate::resetVariables(){
  tRun              = -999;
  tSubRun           = -999;
  tEvent            = -999;
  tChannel          = -999;
  tSPE              = -999.;
  tSPEWidth         = -999.;
  tTimeInFrame_ns   = -999.;
  tTimeInSpill_sec  = -999.;
  tTime_sec         = -999.;
  tIntegral         = -999.;
  tAmplitude        = -999.;
}

DEFINE_ART_MODULE(OpDetPhelRate)
