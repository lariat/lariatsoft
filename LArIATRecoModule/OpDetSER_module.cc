////////////////////////////////////////////////////////////////////////
// Class:       OpDetSER
// Module Type: analyzer
// File:        OpDetSER_module.cc
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
#include <cstdlib>

// ROOT includes
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"

//LAriatSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"

class OpDetSER;

class OpDetSER : public art::EDAnalyzer {
public:
  explicit OpDetSER(fhicl::ParameterSet const & p);
  OpDetSER(OpDetSER const &) = delete;
  OpDetSER(OpDetSER &&) = delete;
  OpDetSER & operator = (OpDetSER const &) = delete;
  OpDetSER & operator = (OpDetSER &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  std::vector<float> FitSER(TH1F* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit);
  void NormalizeWfmHistogram( TH1F* h, int count);

private:
  
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory waveformDir = tfs->mkdir("Waveforms");
  
  // Alg object
  OpHitBuilderAlg     fOpAlg;
  
  // Counters and things
  size_t fNumberActiveEvents[10];
  size_t fNumberSaturatedEvents[10];
  size_t LiveSamples[10];
  size_t WaveformCount[10];
  size_t AvePECount[10];
  int fRunNumber;
  int fSubRunNumber;
  int fEventNumber;
  
  // Timetamp-related variables  
  double fTimestamp;
  
  // Histograms
  TH1F* h_Timestamps;
  TH1F* h_Runs;
  TH1F* h_SER[10];
  TH1F* h_SERWindowSize[10];
  TH1F* h_AvePEWfm[10];
  TH1F* h_RawBaseline[10];
  TH1F* h_RawBaselineRMS[10];
  TH1F* h_BaselineRMS[10];
  TH1F* h_PrePhelRMS[10];
  TH1F* h_Amplitudes[10];
  TH1F* h_Waveform;
  TH1F* h_Waveform_raw;
  TH1F* h_WaveformSER;
  TH1F* h_Live;

  TH1F* h_Tau[10];

  short   fAvePulse_preWindow;
  short   fAvePulse_postWindow;
  int     fAvePulse_total;
  TH1F* h_AvePulse[10];
  TH1F* h_AvePulse_050_150[10];
  TH1F* h_AvePulse_150_250[10];
  TH1F* h_AvePulse_250_350[10];
  TH1F* h_AvePulse_350_450[10];
  TH1F* h_AvePulse_450_550[10];
  TH1F* h_AvePulse_550_650[10];
  TH1F* h_AvePulse_650_750[10];
  TH1F* h_AvePulse_750_850[10];
  TH1F* h_AvePulse_850_950[10];
  int fAvePulseCount[10];
  int fAvePulseCount_050_150[10];
  int fAvePulseCount_150_250[10];
  int fAvePulseCount_250_350[10];
  int fAvePulseCount_350_450[10];
  int fAvePulseCount_450_550[10];
  int fAvePulseCount_550_650[10];
  int fAvePulseCount_650_750[10];
  int fAvePulseCount_750_850[10];
  int fAvePulseCount_850_950[10];

  std::vector<short> fPEWfm_preWindow;
  std::vector<short> fPEWfm_postWindow;
  short fPEWfm_totalSamples[10];

  // Tunable parameters defined by fcl
  bool                fAnalyzePhelWfm;
  bool                fFindSinglePEs;
  std::vector<std::string> fSelectEventTypes;
  std::vector<size_t> fSelectChannels;
  std::vector<size_t> fOpDetChannels;
  std::vector<float>  fOpDetPolarity;
  std::vector<float>  fWfmAbsRMSCut;
  std::vector<float>  fPrePE_RMSFactorCut;
  std::vector<float>  fPulseHitRMSThresh;
  std::vector<float>  fPulseHitThreshHigh;
  std::vector<float>  fPedestalMaxWidth;
  std::vector<float>  fMean_set;
  std::vector<float>  fWidth_set;
  std::vector<float>  fMean_lowerLim;
  std::vector<float>  fMean_upperLim;
  std::vector<float>  fSinglePE;
  std::vector<float>  fSinglePE_tolerance;
  std::vector<short>  fPreIntWindow;
  std::vector<short>  fPostIntWindow;
  std::vector<float>  fSER_x1;
  std::vector<float>  fSER_x2;
  int                 fMinRunNumber;
  int                 fMaxRunNumber;
  bool		      fSaveAvePhelWfm;
  bool                fSaveAvePulses;
  bool                fDoSecondFit;
  size_t              fMaxSavedWaveforms;
  float               fSavedWaveformThresh;
  size_t              fNSamples;
  size_t              fRunningBaselineLength;
  float               fSecondFitChi2Thresh;
  bool                fSubtractRunningBaseline;
  bool                fSubtractMaskedBaseline;
  bool                fVerbose;
  int                 fSmoothingRange;
  short               fBaselineWindowLength;
  short               fT1;
  short               fT2;
  float               fTimestamp_T1;
  float               fTimestamp_T2;
  float               fMaxWindowFactor;
  short               fPrePEBaselineWindow;
  float               fPrePE_TruncP;
  int                 fPrePE_TruncSkew;
  short               fThreshPersist;
  float               fThreshPersistFactor;
  short               fDeadTimeMin;
  short               fQuietTimeMin;
  short               SER_bins[10];
  char		      histName[100];
  char                histTitle[100];
  char		      buffer[200];
  std::string         fDAQModule;
  std::string         fInstanceName;

  bool                fSecondFitDone;

  std::map<size_t,size_t> iCh;
  TCanvas * c;
  
  // Constants
  float               fMvPerADC;

};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
{
  this                ->reconfigure(p);
  fAvePulse_total = int(fAvePulse_preWindow + fAvePulse_postWindow);

  std::cout<<"Saving average pulses? "<<fSaveAvePulses<<"\n";
  fRunNumber = 0;
  fSubRunNumber = 0;
  fEventNumber = 0; 
  
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    iCh[fOpDetChannels[i]]    = i; // map from channel to index
    SER_bins[i]               = fPreIntWindow[i]+fPostIntWindow[i];
    WaveformCount[i]          = 0;
    LiveSamples[i]            = 0;
    fNumberSaturatedEvents[i] = 0;
    fNumberActiveEvents[i]    = 0;
    AvePECount[i]             = 0;
   
    fAvePulseCount[i]         = 0; 
    fAvePulseCount_050_150[i] = 0;
    fAvePulseCount_150_250[i] = 0;
    fAvePulseCount_250_350[i] = 0;
    fAvePulseCount_350_450[i] = 0;
    fAvePulseCount_450_550[i] = 0;
    fAvePulseCount_550_650[i] = 0;
    fAvePulseCount_650_750[i] = 0;
    fAvePulseCount_750_850[i] = 0;
    fAvePulseCount_850_950[i] = 0;
  
    fPEWfm_totalSamples[i] = fPEWfm_preWindow[i] + fPEWfm_postWindow[i];
  }

}

//#######################################################################
void OpDetSER::reconfigure(fhicl::ParameterSet const & p)
{
  std::vector<size_t>     SelectChannelsDefault{0, 1};
  std::vector<size_t>     BaselineChannelDefault{};
  std::vector<float>      MeanSetDefault{60, 60};
  std::vector<float>      WidthSetDefault{ 15., 15.};
  std::vector<float>      MeanLowerLimDefault{30, 30};
  std::vector<float>      MeanUpperLimDefault{90, 90};
  std::vector<float>      PrePeRMSFactorCutDefault{1.5, 1.5};
  std::vector<float>      PulseHitThreshHighDefault{8., 8.};
  std::vector<float>      PulseHitRMSThreshDefault{ 4., 4.};
  std::vector<float>      SinglePEDefault{ 60., 60.};
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
  std::vector<short>      PEWfm_preWindowDefault{20, 20};
  std::vector<short>      PEWfm_postWindowDefault{50, 50};

  fFindSinglePEs          = p.get< bool >                 ("FindSinglePEs",true);
  fAnalyzePhelWfm         = p.get< bool >                 ("AnalyzePhelWfm",true);
  fSER_x1                 = p.get< std::vector<float> >   ("SER_x1",x1Default);
  fSER_x2                 = p.get< std::vector<float> >   ("SER_x2",x2Default);
  fSelectChannels         = p.get< std::vector<size_t> >  ("SelectChannels",SelectChannelsDefault);
  fOpDetPolarity        = p.get< std::vector<float> >   ("OpDetPolarity",PolarityDefault);
  fOpDetChannels          = p.get< std::vector<size_t> >  ("OpDetChannels", SelectChannelsDefault);
  fMean_set               = p.get< std::vector<float> >   ("Mean_set",MeanSetDefault);
  fWidth_set              = p.get< std::vector<float> >   ("Width_set",WidthSetDefault);
  fMean_lowerLim          = p.get< std::vector<float> >   ("Mean_lowerLim",MeanLowerLimDefault);
  fMean_upperLim          = p.get< std::vector<float> >   ("Mean_upperLim",MeanUpperLimDefault);
  fPrePE_RMSFactorCut     = p.get< std::vector<float> >   ("PrePE_RMSFactorCut",PrePeRMSFactorCutDefault);
  fPulseHitThreshHigh     = p.get< std::vector<float> >   ("PulseHitThresh_high",PulseHitThreshHighDefault);
  fPulseHitRMSThresh      = p.get< std::vector<float> >   ("PulseHitRMSThresh",PulseHitRMSThreshDefault);
  fSinglePE               = p.get< std::vector<float> >   ("SinglePE",SinglePEDefault);
  fSinglePE_tolerance     = p.get< std::vector<float> >   ("SinglePE_tolerance",SinglePEToleranceDefault);
  fWfmAbsRMSCut           = p.get< std::vector<float> >   ("WfmAbsRMSCut", WfmAbsRMSCutDefault);
  fPedestalMaxWidth       = p.get< std::vector<float> >   ("PedestalMaxWidth",PedestalMaxWidthDefault);
  fPreIntWindow              = p.get< std::vector<short> >   ("PreIntWindow",PreWindowDefault);
  fPostIntWindow             = p.get< std::vector<short> >   ("PostIntWindow",PostWindowDefault);

  fMinRunNumber           = p.get< int >          ("MinRunNumber",-999);
  fMaxRunNumber           = p.get< int >          ("MaxRunNumber",-999);
  fSaveAvePhelWfm	  = p.get< bool >	  ("SaveAvePhelWfm",true);
  fPEWfm_preWindow        = p.get< std::vector<short> >          ("PEWfm_preWindow",PEWfm_preWindowDefault);
  fPEWfm_postWindow       = p.get< std::vector<short> >          ("PEWfm_postWindow",PEWfm_postWindowDefault);
  fSaveAvePulses	  = p.get< bool >	  ("SaveAvePulses",false);
  fDoSecondFit            = p.get< bool >         ("DoSecondFit",true);
  fSecondFitChi2Thresh    = p.get< float >        ("SecondFitChi2Thresh",1.2);
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fSmoothingRange         = p.get< int >          ("SmoothingRange",0);
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fTimestamp_T1           = p.get< float >        ("Timestamp_T1",0);
  fTimestamp_T2           = p.get< float >        ("Timestamp_T2",60);
  fMaxWindowFactor        = p.get< float >        ("MaxWindowFactor",2);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.1953);
  fPrePEBaselineWindow    = p.get< short >        ("PrePEBaselineWindow",100);
  fPrePE_TruncP           = p.get< float >        ("PrePE_TruncP",0.05);
  fPrePE_TruncSkew        = p.get< int  >         ("PrePE_TruncSkew",0);
  fThreshPersist          = p.get< short >        ("ThreshPersist",0);
  fThreshPersistFactor    = p.get< float >        ("ThreshPersistFactor",1.0);
  fDeadTimeMin            = p.get< short >        ("DeadTimeMin",10);
  fQuietTimeMin           = p.get< short >        ("QuietTimeMin",200);
  fVerbose                = p.get< bool >         ("Verbose",true);
  fNSamples               = p.get< size_t >        ("NSamples",28672);
  fSelectEventTypes       = p.get< std::vector<std::string> > ("SelectEventTypes",EventTypeDefaults);
  fMaxSavedWaveforms      = p.get< size_t >       ("MaxSavedWaveforms",10);
  fSavedWaveformThresh    = p.get< float >        ("SavedWaveformThresh",50);
  fSubtractRunningBaseline= p.get< bool >         ("SubtractRunningBaseline",false);
  fSubtractMaskedBaseline= p.get< bool >         ("SubtractMaskedBaseline",false);
  fRunningBaselineLength  = p.get< size_t >         ("RunningBaselineLength",256);
  fAvePulse_preWindow     = p.get< short >          ("AvePulse_preWindow",1000);
  fAvePulse_postWindow    = p.get< short >          ("AvePulse_postWindow",14000);
}

//#######################################################################
void OpDetSER::beginJob()
{
  h_Timestamps  = tfs->make<TH1F>("Timestamps","Timestamp for all events;sec",240,0.,60.);
  h_Runs        = tfs->make<TH1F>("Runs","Runs analyzed;Run number",10000,6000,16000);

  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ch   = fSelectChannels[i];
    size_t ich  = iCh[ch]; 
    int x_bins  = std::min((int)(fSER_x2[ich] - fSER_x1[ich]),500);
    
    sprintf(histName,"%lu_SER",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integrated ADC",ch);
    h_SER[ich]		  = tfs->make<TH1F>(histName,histTitle,x_bins,fSER_x1[ich],fSER_x2[ich]);
    
    sprintf(histName,"%lu_SERWindowSize",ch); 
    sprintf(histTitle,"SER window sizes for OpDet %lu;Window size (samples);Counts",ch);
    h_SERWindowSize[ich]       = tfs->make<TH1F>(histName,histTitle,100,0.,fMaxWindowFactor*SER_bins[ich]+10);
  
    sprintf(histName,"%lu_RawBaseline",ch); 
    sprintf(histTitle,"Raw baselines for OpDet %lu;Baseline [ADC];Counts",ch);
    h_RawBaseline[ich]    = tfs->make<TH1F>(histName,histTitle,600,0.0,1200.);
    
    sprintf(histName,"%lu_RawBaselineRMS",ch); 
    sprintf(histTitle,"Raw RMS for OpDet %lu;Baseline RMS [mV];Counts",ch);
    h_RawBaselineRMS[ich]    = tfs->make<TH1F>(histName,histTitle,300,0.0,1.5);
    
    sprintf(histName,"%lu_BaselineRMS",ch); 
    sprintf(histTitle,"RMS for OpDet %lu (after corrections);Baseline RMS [mV];Counts",ch);
    h_BaselineRMS[ich]    = tfs->make<TH1F>(histName,histTitle,300,0.0,1.5);
    
    sprintf(histName,"%lu_PrePhelRMS",ch); 
    sprintf(histTitle,"RMS of pre-PE candidate baselines, OpDet %lu;Pre-PE RMS [mV];Counts",ch);
    h_PrePhelRMS[ich]     = tfs->make<TH1F>(histName,histTitle,300,0.0,1.5);
    
    sprintf(histName,"%lu_Amplitudes",ch); 
    sprintf(histTitle,"Waveform amplitudes for OpDet %lu;Amplitude [mV];Counts",ch);
    h_Amplitudes[ich]    = tfs->make<TH1F>(histName,histTitle,2000.,0.,200.);
    
    if(fSaveAvePhelWfm){ 
      sprintf(histName,"%lu_AvePEWfm",ch); 
      sprintf(histTitle,"Average PE waveform for OpDet %lu (%5.1f ADC, %4.2f tolerance);ns;ADC",ch, fSinglePE[ch], fSinglePE_tolerance[ch]);
      //h_AvePEWfm[ich]  = tfs->make<TH1F>(histName,histTitle,SER_bins[ich],0.,(float)SER_bins[ich]);
      h_AvePEWfm[ich]  = tfs->make<TH1F>(histName,histTitle,fPEWfm_totalSamples[ich],-1.*fPEWfm_preWindow[ich],fPEWfm_postWindow[ich]);
      //h_AvePEWfm[ich]->SetBit(TH1::kIsAverage);
    }
  
  
    if( fSaveAvePulses ){
  
      art::TFileDirectory avePulseDir = tfs->mkdir("AvePulses");
      
      h_Tau[ich] = tfs->make<TH1F>(Form("%lu_Tau",ch),Form("Tau from fit to individual pulses;#tau [ns]"),150,0.,3000);

      sprintf(histName,"%lu_AvePulse",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu;ns;ADC",ch);
      h_AvePulse[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
       
      sprintf(histName,"%lu_AvePulse_050_150",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 50-150 ADC;ns;ADC",ch);
      h_AvePulse_050_150[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_150_250",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 150-250 ADC;ns;ADC",ch);
      h_AvePulse_150_250[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_250_350",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 250-350 ADC;ns;ADC",ch);
      h_AvePulse_250_350[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_350_450",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 350-450 ADC;ns;ADC",ch);
      h_AvePulse_350_450[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_450_550",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 450-550 ADC;ns;ADC",ch);
      h_AvePulse_450_550[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_550_650",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 550-650 ADC;ns;ADC",ch);
      h_AvePulse_550_650[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_650_750",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 650-750 ADC;ns;ADC",ch);
      h_AvePulse_650_750[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_750_850",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 750-850 ADC;ns;ADC",ch);
      h_AvePulse_750_850[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
      
      sprintf(histName,"%lu_AvePulse_850_950",ch); 
      sprintf(histTitle,"Average pulse, OpDet %lu, amplitude 850-950 ADC;ns;ADC",ch);
      h_AvePulse_850_950[ich]  = avePulseDir.make<TH1F>(histName,histTitle,fAvePulse_total, -1.*fAvePulse_preWindow, fAvePulse_postWindow);
    }
  }
}


//#######################################################################
void OpDetSER::analyze(art::Event const & e)
{
  fEventNumber   = e.id().event();
  int runnr     = e.run();
  fSubRunNumber  = e.subRun();

  if( 
      (fMinRunNumber > 0 && runnr < fMinRunNumber)
    ||(fMaxRunNumber > 0 && runnr > fMaxRunNumber) ) return;

  if( fVerbose ) std::cout<<"Analyzing run "<<runnr<<", subrun "<<fSubRunNumber<<", event "<<fEventNumber<<"\n";
 
  h_Runs->Fill(runnr);

  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);
  if( (size_t)WaveformHandle->size() == 0 ){
    if(fVerbose) std::cout << "No optical detector data found; skipping event.\n";
    return;
  }
  
  // Grab the first OpDetPulse in the handle just to check whether its
  // the right event type for use in this analysis
  fTimestamp   = (float(WaveformHandle->at(0).PMTFrame())*8.)/1.0e09;
  if ( (fTimestamp_T1 > 0 && fTimestamp < fTimestamp_T1) || (fTimestamp_T2 > 0 && fTimestamp > fTimestamp_T2) ) { 
    if(fVerbose) std::cout<<"Timestamp ("<<fTimestamp<<" sec) out of range --> skipping event\n";
    return;
  }


  h_Timestamps->Fill(fTimestamp);

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
    
    if( fVerbose ) std::cout<<"Analyzing opdet "<<ch<<"\n";

    // Reset live samples counter
    size_t liveSamples = 0;

    // Get the raw waveform
    std::vector<short>  Wfm_raw   = ThePulse.Waveform();
    size_t              NSamples  = Wfm_raw.size();
    
    // Convert this to a vector<float> for easier manipulation later on
    std::vector<float> Wfm( Wfm_raw.begin(), Wfm_raw.end() );

    // Get raw quantities
    fOpAlg.CalcBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
    h_RawBaseline[ich]      ->Fill(fOpAlg.GetBaseline());
    h_RawBaselineRMS[ich]   ->Fill(fOpAlg.GetRMS()*fMvPerADC);

    // --------------------------------------------------------------------
    // Smooth out waveform
    fOpAlg.SmoothOutVector(Wfm,fSmoothingRange);
    fOpAlg.CalcBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
    fOpAlg.SubtractBaseline( Wfm );
    
    // Masked baseline subtraction
    if( fSubtractMaskedBaseline ) {
      std::vector<float> vbs(Wfm.size(),0.);
      fOpAlg.MaskedBaselineSubtraction(Wfm, vbs);
      for(size_t i=0; i<Wfm.size(); i++) Wfm.at(i) = Wfm.at(i) - vbs.at(i);
    }
  
    // Subtract running baseline to remove oscillations (useful for random pulser runs where 
    // the signal is expected to remain approximately at baseline).
    if(fSubtractRunningBaseline) {
      if(fVerbose) std::cout<<"Performing running baseline subtraction...\n";
      std::vector<float> Wfm2(Wfm.size(), 0.);
      fOpAlg.SubtractRunningBaseline(Wfm, Wfm2, NSamples, fRunningBaselineLength);
      Wfm = Wfm2;
    }
    // --------------------------------------------------------------------
    
    // Find waveform baseline and RMS
    fOpAlg.CalcBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
    float baseline  = fOpAlg.GetBaseline();
    float rms       = fOpAlg.GetRMS();
    h_BaselineRMS[ich]      ->Fill(rms*fMvPerADC);
    
    // Skip if RMS is too noisy
    if( rms*fMvPerADC >= fWfmAbsRMSCut[ich] ) {
      if(fVerbose) std::cout<<"Baseline too noisy on channel "<<ich<<" ("<<rms*fMvPerADC<<" mV)... skipping this event.\n";
      continue;
    } 
    
    if(fVerbose){ 
      std::cout << "OpDet "<<ch<<", index "<<ich<<" Timestamp "<<fTimestamp<<" sec\n";
      std::cout << "Waveform raw baseline: " << baseline << " +/- " << rms <<" ADC\n";
    }
    
    // Baseline subtraction and signal inversion.
    std::vector<float> wfm_corrected(NSamples);
    float amp = -1.;
    bool isSaturated = false;
    for(size_t i=0; i<NSamples; i++) {
      wfm_corrected[i] = fOpDetPolarity[ich]*((float)Wfm[i]-baseline);
      if( wfm_corrected[i] > amp  ) amp = wfm_corrected[i];
      if( Wfm_raw[i] == 0         ) isSaturated = true;
    }
    h_Amplitudes[ich]	->Fill(amp*fMvPerADC);
    if( amp > 20.*rms ) {
      fNumberActiveEvents[ich]++;
      if( isSaturated )	fNumberSaturatedEvents[ich]++;
    }


      
    bool areSavingWaveform = false; 
    if( WaveformCount[ich] < fMaxSavedWaveforms && amp > fSavedWaveformThresh ){
	
      //sprintf(histName,"Ch%lu_wfmRaw%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,fSubRunNumber,fEventNumber);
      sprintf(histName,"Ch%lu_wfmRaw%lu",ch,WaveformCount[ich]);
      sprintf(histTitle,"OpDet%lu (raw): run %i, subrun %i, event %i;Sample [ns];Signal [ADC]",ch,runnr,fSubRunNumber,fEventNumber); 
      //h_Waveform_raw	= waveformDir.make<TH1F>(histName,histTitle,Wfm_raw.size(),0,Wfm_raw.size());
      //h_Waveform_raw  ->SetOption("HIST");
      
        
      //sprintf(histName,"Ch%lu_wfm%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,fSubRunNumber,fEventNumber);
      sprintf(histName,"Ch%lu_wfm%lu",ch,WaveformCount[ich]);
      sprintf(histTitle,"OpDet%lu: run %i, subrun %i, event %i;Sample [ns];Signal [ADC]",ch,runnr,fSubRunNumber,fEventNumber); 
      h_Waveform      = waveformDir.make<TH1F>(histName,histTitle,Wfm.size(),0,Wfm.size());
      h_Waveform      ->SetOption("HIST");
	
      sprintf(histName,"Ch%lu_SEROverlay_wfm%lu_r%i_sr%i_e%i",ch,WaveformCount[ich],runnr,fSubRunNumber,fEventNumber);
      sprintf(histTitle,"OpDet%lu: run %i, subrun %i, event %i;Sample [ns];Signal [mV]",ch,runnr,fSubRunNumber,fEventNumber); 
      c               = waveformDir.make<TCanvas>(histName,"c",700,500);

      h_WaveformSER   = new TH1F("serwindow","serwindow",Wfm.size(),0,Wfm.size());
      h_WaveformSER   ->SetLineColor(kRed); 
      h_Live   = new TH1F("livetime","livetime",Wfm.size(),0,Wfm.size());
      h_Live   ->SetLineColor(kGreen);
      
      // Save raw and corrected waveforms, and initialize
      // the livetime and SER histograms 
      for(short ii=0; ii< short(Wfm.size()); ii++) {
        h_Live            ->SetBinContent(ii+1,0.01);
        h_WaveformSER     ->SetBinContent(ii+1,0.01);
	//h_Waveform_raw    ->Fill(ii,(float)Wfm_raw[ii]);
	h_Waveform	  ->Fill(ii,wfm_corrected[ii]);
      }
        
      areSavingWaveform = true;
    }
      
   
    // ----------------------------------------------------
    // Saving average pulses
    if( fSaveAvePulses ) {
       
       std::cout<<"Saving average pulse\n"; 
      // Perform hit-finding
      std::vector<short> vHitTimesTmp = fOpAlg.GetHits(ThePulse);
     
      std::cout<<"Found "<<vHitTimesTmp.size()<<" hits\n"; 

      size_t nhits = vHitTimesTmp.size();

      for( size_t iHit = 0; iHit < nhits; iHit++){
      
        short hittime = vHitTimesTmp[iHit];

        // If hit saturates, skip it:
        //if( fOpAlg.GetLocalMinimum(Wfm, hittime) == 0 ) continue;
     
        // Make sure the hit is sufficiently isolated
        float t_prev = -9.;
        float t_next = -9.;
        if( iHit > 0 ) t_prev = vHitTimesTmp[iHit-1];
        if( iHit < nhits-1 ) t_next = vHitTimesTmp[iHit+1];
        if( t_prev > 0 && hittime-t_prev < 5000 ) continue;
        if( t_next > 0 && t_next-hittime < 2500 ) continue;

        // Require hit be at trigger point
        //std::cout<<"Is hit at trigger point? "<<hittime<<"   "<<ThePulse.FirstSample()<<"   "<<Wfm.size()*0.015<<"\n";
        //if( fabs( hittime - short(ThePulse.FirstSample()) ) >= Wfm.size()*0.05 ) continue; 
        //std::cout<<"yes\n";
        
        // Require flat pre-hit region
        //float preHitRMS = fOpAlg.GetBaselineAndRMS(wfm_corrected, hittime - 1000, hittime)[1];
//        std::cout<<"Is pre-hit region flat? "<<preHitRMS<<"   ( "<<rms<<")\n";
        //if( preHitRMS >= 1.5 * rms ) continue;

        float a = fOpAlg.GetLocalMaximum(wfm_corrected, hittime);
        int x1 = (int)hittime - (int)fAvePulse_preWindow;
        int x2 = (int)hittime + (int)fAvePulse_postWindow;
        
        if( x1 >= 0 && x2 < (int)wfm_corrected.size() && a >= 50. && a < 950. ) {
        
          TH1D *h_wfm = (TH1D*)h_AvePulse[ich]->Clone("tmp");
          h_wfm->Reset();
          for(int jj = 0; jj < fAvePulse_total; jj++) {
            int bin = h_wfm->GetXaxis()->FindBin(-1.*fAvePulse_preWindow+jj);
            h_wfm -> SetBinContent(bin, wfm_corrected[x1+jj]);
            h_wfm -> SetBinError( bin, rms);
            //h_wfm -> Fill(-1.*fAvePulse_preWindow+jj, wfm_corrected[x1 + jj]);
          }
          
          TF1 tau("tau","[0]*exp(-x/[1])",400,2000);
          tau.SetParameter(0,50);
          tau.SetParameter(1,1000);
          tau.SetParLimits(1,0,3000);
          h_wfm->Fit("tau","QR");
          if( tau.GetParameter(1) > 0 ) h_Tau[ich]->Fill(tau.GetParameter(1));
           
          std::cout<<"Adding pulse\n";
          h_AvePulse[ich]->Add( h_wfm ); fAvePulseCount[ich]++;
          if(       a >= 50.  && a < 150. ) { h_AvePulse_050_150[ich]->Add( h_wfm ); fAvePulseCount_050_150[ich]++; }
          else if(  a >= 150. && a < 250. ) { h_AvePulse_150_250[ich]->Add( h_wfm ); fAvePulseCount_150_250[ich]++; }
          else if(  a >= 250. && a < 350. ) { h_AvePulse_250_350[ich]->Add( h_wfm ); fAvePulseCount_250_350[ich]++; }
          else if(  a >= 350. && a < 450. ) { h_AvePulse_350_450[ich]->Add( h_wfm ); fAvePulseCount_350_450[ich]++; }
          else if(  a >= 450. && a < 550. ) { h_AvePulse_450_550[ich]->Add( h_wfm ); fAvePulseCount_450_550[ich]++; }
          else if(  a >= 550. && a < 650. ) { h_AvePulse_550_650[ich]->Add( h_wfm ); fAvePulseCount_550_650[ich]++; }
          else if(  a >= 650. && a < 750. ) { h_AvePulse_650_750[ich]->Add( h_wfm ); fAvePulseCount_650_750[ich]++; }
          else if(  a >= 750. && a < 850. ) { h_AvePulse_750_850[ich]->Add( h_wfm ); fAvePulseCount_750_850[ich]++; }
          else if(  a >= 850. && a < 950. ) { h_AvePulse_850_950[ich]->Add( h_wfm ); fAvePulseCount_850_950[ich]++; }
          
          delete h_wfm;

        }


      }


    } //<-- done saving ave pulses
    // -------------------------------------
     
     
    if( fFindSinglePEs ) {
       
    // Determine range of samples to scan over
    short t1 = std::max((short)ThePulse.FirstSample() + fT1, (int)fPrePEBaselineWindow);
    short t2 = std::min((short)ThePulse.FirstSample() + fT2, (int)NSamples-10);
      
    // Declare and initialize variables
    float   integral = 0;
    bool    flag = false;
    int     windowsize = 0;
    int     counter = -(int)fPreIntWindow[ich];
    int     dipcounter = 0;
    int     overthreshcounter = 0;
    short   quietTime = 0;
    short   deadTime = 0;
    float   prePE_baseline = 0;
    float   prePE_rms = 999;
    int     hit_time = -1;
    float   PE_amp  = -1;
   
    // Now begin looking for single PEs
    if(fVerbose){
    std::cout
    << "OpDet "<<ch<<"\n"
    << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh[ich] 
    << ", threshPersist " << fThreshPersist <<"/"<<fThreshPersistFactor<<") \n";
    }

    
    // ------------------------------------------------------------
    // Scan waveform and search for hits within specified time window
    for(short i = t1; i < t2; i++){
        
      if(!flag) prePE_baseline = 0.;
      float y               = wfm_corrected[i] - prePE_baseline;
      float ythresh         = fPulseHitRMSThresh[ich]*rms;
      float y_mV            = y*fMvPerADC;
      bool  IsLive          = ((quietTime >= fQuietTimeMin)&&(deadTime >= fDeadTimeMin));
      bool  IsOverThresh    = (y  >= ythresh);
      if(   IsOverThresh )  { overthreshcounter++;}
      else                  { overthreshcounter=0;}
      //bool  IsOverThreshNeg = (y  <= -3.0*rms);
      bool IsOverThreshNeg = (y <= -1.*fPulseHitRMSThresh[ich]*rms);
      if(IsOverThreshNeg)   { dipcounter++;}
      else{                   dipcounter=0;}
      bool  NegativeDipDetected = (dipcounter > 3);
      bool  IsOverLimit     = (y_mV >= fPulseHitThreshHigh[ich]); 
      bool  IsPECandidate   = ((IsOverThresh) && (!IsOverLimit) );
     
      //printf("%8d  y=%6.2f   thresh=%6.2f   overthreshcounter %5d   dipcounter %5d   quietTime %5d   deadTime %5d   live? %3d\n",i,y,ythresh,overthreshcounter,dipcounter,quietTime,deadTime,IsLive);
      
      if(IsLive || flag) {
        liveSamples++;
        if( areSavingWaveform ) h_Live->SetBinContent(i+1,ythresh);
      }
      
        
      // If we're already in a PE window, increment the
      // counters and add to integral
      if( flag ) {
        
        if(fVerbose) {
        std::cout
        << "  " << i << "  y = " << y << " prePE_bs = "<<prePE_baseline<<" (window size " <<windowsize+1<< "), "
        << " thresh " << ythresh << ", quietTime = "<<quietTime<<"   deadTime = "<<deadTime<<"\n";
        }

        counter++;
        windowsize++;
        integral += y;
       
        // Scan for amplitude
        if( y_mV > PE_amp ) PE_amp = y_mV;

        // If (a) a "dip" below baseline was detected, (b) a pulse extends above upper 
        // limit, or (c) if window length is exceeded due to multiple merges, abort mission.
        if( NegativeDipDetected || IsOverLimit || (windowsize > fMaxWindowFactor*SER_bins[ich])){
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
          continue;
        }
          
	//----------------------------------
	// End integration if we've reached the minimum required window size and
	// the signal has returned to baseline
        if( counter >= (int)fPostIntWindow[ich] && deadTime >= 10 ){
	 
          if(fVerbose) printf("  * %f ADCs, windowsize %d \n",integral,windowsize);
         
          h_SER[ich]          ->Fill(integral);
          h_SERWindowSize[ich]->Fill(windowsize);
         
          // Add to average waveform if it looks good
          if( fSaveAvePhelWfm && windowsize == SER_bins[ich] && fabs(integral - fSinglePE[ich])/fSinglePE[ich] <= fSinglePE_tolerance[ich] ){

            // Check ahead to make check for secondary pulses within the ave PE window
            short kk = hit_time + fPEWfm_postWindow[ich];
            bool goodWfm = true;
            for(short j=i; j<kk; j++){
              if( fabs(wfm_corrected.at(j)) > ythresh ) {
                goodWfm = false;
                break;
              }
            }

            if( goodWfm ) {
              if(fVerbose) std::cout<<"    --> Adding to average PE waveform (prePE baseline = "<<prePE_baseline<<")\n";
              for(short ii=0; ii<fPEWfm_totalSamples[ich]; ii++) h_AvePEWfm[ich]->Fill(-1.*fPEWfm_preWindow[ich]+ii, wfm_corrected[hit_time-1.*fPEWfm_preWindow[ich]+ii]-prePE_baseline);
              AvePECount[ich]++;
            }
          }

          // Mark on waveform
          if(areSavingWaveform){
            short start = hit_time-fPreIntWindow[ich];
            short end   = start+windowsize;
            //for(short ii=start; ii<end; ii++) h_WaveformSER ->SetBinContent(ii+1,ythresh*fMvPerADC);
            for(short ii=start; ii<end; ii++) h_WaveformSER ->SetBinContent(ii+1,ythresh);
          }
    
          integral = 0;
          hit_time = -1;
          PE_amp    = -1;
          windowsize = 0;
          counter = -(int)fPreIntWindow[ich];
          flag = false;
          continue;
        }
	//----------------------------------
        
      } // endif flag
   
        
       
      // If we're not yet integrating a PE window, signal is within bounds,
      // and enough quietTime has elapsed, then we're in business
      if( !flag && IsPECandidate && IsLive ){
        
        bool flagger = true;

        if( fVerbose ) std::cout<<"PE candidate found...\n";
          
        // Find pre-PE baseline
        prePE_baseline  = 99;
        prePE_rms	  = 99;

        short k1 = i-fPrePEBaselineWindow-fPreIntWindow[ich];
        short k2 = i-fPreIntWindow[ich];

        //std::vector<float> vals;
        //for(short k=k1; k<k2; k++) vals.push_back(wfm_corrected.at(k));
        //prePE_baseline = fOpAlg.CalcTruncatedMean(vals,fPrePE_TruncP,fPrePE_TruncSkew);
    
        std::vector<float> tmp = fOpAlg.GetBaselineAndRMS(wfm_corrected,k1,k2);
        prePE_baseline  = tmp[0];
        prePE_rms	  = tmp[1];
        h_PrePhelRMS[ich] ->Fill(prePE_rms*fMvPerADC);

        // Check shorter prePE baseline
        //k1 = i-fPreIntWindow[ich]*10;
        //k2 = i-fPreIntWindow[ich];
        //float bs2 = fOpAlg.GetBaselineAndRMS(wfm_corrected,k1,k2)[0];
        //if( fabs( bs2 - prePE_baseline ) > prePE_rms ) flagger = false; 


        //prePE_baseline = 0.;
        
        // When taking the pre-PE baseline into account, require that the 
        // candidate still passes the threshold:
        y             = wfm_corrected[i] - prePE_baseline;
        y_mV          = y*fMvPerADC;
        IsOverThresh  = (y    >= ythresh);
        IsOverLimit   = (y_mV >= fPulseHitThreshHigh[ich]); 
        IsPECandidate = ((IsOverThresh) && (!IsOverLimit)); 
       
        // If "threshPersist" functionality being used, look ahead and make 
        // sure signal passes the threshold for requisite consecutive samples
        short tmp_counter = 0;
        if( fThreshPersist > 0 ){
          for( short j=1; j <= fThreshPersist; j++ ){
            if ( wfm_corrected[i+j] - prePE_baseline > fThreshPersistFactor*ythresh ) tmp_counter++;
          }
          if (tmp_counter < fThreshPersist) flagger = false;
        }
        
         
        if( fVerbose ){
          std::cout
          << "  "<<i<<"   "<<y_mV<<" (thresh "<<ythresh*fMvPerADC<<")   wfm rms ="<<rms*fMvPerADC<<"   prePE rms = "<<prePE_rms*fMvPerADC<<"\n";
        }

        // Require flat pre-PE region and that threshold persist requirement was met
        if( (IsPECandidate) && (prePE_rms <= fPrePE_RMSFactorCut[ich]*rms) && (flagger) ){
          
          // Found a "PE hit"
          flag      = true;
          hit_time  = i;
         
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
      } else { 
        quietTime = 0;  
      }

      // "deadTime" is similar, but with stricter conditions.
      // Deadtime resets whenever a PE candidate is found,
      // or when the signal drops below neg threh
      if( !IsOverThresh && !IsOverThreshNeg ) { 
        deadTime++; 
      } else { 
        deadTime = 0; 
      }


    } // <-- end scan over waveform 
    if(fVerbose) std::cout<<"Ending scan over waveform\n";
      
    if( areSavingWaveform ){
	//h_Waveform_raw->SetOption("hist");
	h_Waveform    ->SetOption("hist");
        WaveformCount[ich]++;
        
        c->cd();
        h_Waveform->DrawCopy("hist");
        h_Live->Draw("hist same");
        h_WaveformSER->Draw("hist same");
        c->Write(histName);
        delete c;
        delete h_Waveform;
        delete h_WaveformSER;
        //delete h_Live; //<-- causes segfault?
    }

    LiveSamples[ich]         += liveSamples;
    if(fVerbose) std::cout<<"Live samples for this event... "<<liveSamples<<"\n"; 
 
  }// endif search for single pes
    
  } // endLoop over OpDets
    
}


//#######################################################################
void OpDetSER::beginSubRun(art::SubRun const & sr)
{
}

//#######################################################################
void OpDetSER::endJob()
{
  std::cout 
    << "\n============================================================\n"
    << "Ending OpDetSER\n";
  
  for(size_t i=0; i<fSelectChannels.size(); i++){
    
    size_t ich = iCh[fSelectChannels[i]];

    // Calculate the total cumulative time the SER code was "live"
    float cumulativeTime = (float)LiveSamples[ich]*1.0e-9;

    // Perform SER fit and save the parameters
    std::vector<float> params(17);
    params = FitSER(h_SER[ich],fSER_x1[ich], fSER_x2[ich], fMean_set[ich], fWidth_set[ich], fPedestalMaxWidth[ich], fDoSecondFit);
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
    
    // Print out photodetector-specific parameters and SER fit results
    std::cout
    << "\n**************************************\n"
    << "Channel	" << fOpDetChannels[ich] << "\n"
    << "  Wfm smoothing range       " << fSmoothingRange<<"\n"
    << "  Running BS subtr          " << fSubtractRunningBaseline<<" (r="<<fRunningBaselineLength<<")\n"
    << "  Masked BS subtr           " << fSubtractMaskedBaseline<<" (gradRmsF="<<fOpAlg.fMskBaselineSubtr_grmsfac<<", P="<<fOpAlg.fMskBaselineSubtr_P<<", r="<<fOpAlg.fMskBaselineSubtr_range<<")\n"
    << "  RMS thresh factor	    " << fPulseHitRMSThresh[ich] <<" x wfm RMS\n"
    << "  Threshhold persist	    " << fThreshPersist <<" (x "<<fThreshPersistFactor<<")\n"
    << "  Quiet time min            " << fQuietTimeMin <<"\n"
    << "  Dead time min             " << fDeadTimeMin << "\n"
    << "  Pre-PE baseline           " << fPrePEBaselineWindow << "\n"
    << "  Pre-PE trunc. (P,skew)    " << fPrePE_TruncP<<", "<<fPrePE_TruncSkew<<"\n"
    << "  Pre-PE RMS cut	    " << fPrePE_RMSFactorCut[ich] << " x wfm RMS\n"
    << "  Pre / Post window size    " << fPreIntWindow[ich] << "," << fPostIntWindow[ich] << "\n"
    << "  Active events	(>20*RMS)   " << fNumberActiveEvents[ich] << "\n"
    << "  Saturation rate	    " << (float)fNumberSaturatedEvents[ich] / (float)fNumberActiveEvents[ich] << "\n"
    << "  Mean baseline RMS	    " << h_BaselineRMS[ich]->GetMean(1) << " mV\n"
    << "  --> PE candidates found   " << h_SER[ich]->GetEntries() << "\n"
    <<"-----------------------------------------\n"
    <<"SER fit results for OpDet "<< fOpDetChannels[ich]<<":\n"
    <<"  2nd fit?        = "<<fSecondFitDone<<"\n"
    <<"  Noise N         = "<<n0<<" +/- "<<dn0<<"\n"
    <<"  Noise mean      = "<<mu0<<" +/- "<<dmu0<<"\n"
    <<"  Noise sigma     = "<<sig0<<" +/- "<<dsig0<<"\n"
    <<"  1PE mean (SER)  = "<<mu<<" +/- "<<dmu<<"\n"
    <<"  1PE sigma       = "<<sig<<" +/- "<<dsig<<"\n"
    <<"  1PE N           = "<<n1<<" +/- "<<dn1<<"\n"
    <<"  2PE N           = "<<n2<<" +/- "<<dn2<<" (r="<<n2/n1<<")\n"
    <<"  3PE N           = "<<n3<<" +/- "<<dn3<<" (r="<<n3/n1<<")\n"
    <<"  chisquare / NDF = "<<chi<<"\n"
    <<"  live time       = "<<cumulativeTime<<" sec\n";
    
   
    // Normalize the average pulse histograms
    if( fSaveAvePulses ) {
   
       
      std::cout<<"Ave pulse count             : "<<fAvePulseCount[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 50-150 ADC : "<<fAvePulseCount_050_150[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 150-250 ADC: "<<fAvePulseCount_150_250[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 250-350 ADC: "<<fAvePulseCount_250_350[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 350-450 ADC: "<<fAvePulseCount_350_450[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 450-550 ADC: "<<fAvePulseCount_450_550[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 550-650 ADC: "<<fAvePulseCount_550_650[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 650-750 ADC: "<<fAvePulseCount_650_750[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 750-850 ADC: "<<fAvePulseCount_750_850[ich]<<"\n"; 
      std::cout<<"Ave pulse count, 850-950 ADC: "<<fAvePulseCount_850_950[ich]<<"\n"; 

    
      if( fAvePulseCount[ich]         > 0 ){ NormalizeWfmHistogram(h_AvePulse[ich], fAvePulseCount[ich]);}
      if( fAvePulseCount_050_150[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_050_150[ich], fAvePulseCount_050_150[ich]);}
      if( fAvePulseCount_150_250[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_150_250[ich], fAvePulseCount_150_250[ich]);}
      if( fAvePulseCount_250_350[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_250_350[ich], fAvePulseCount_250_350[ich]);}
      if( fAvePulseCount_350_450[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_350_450[ich], fAvePulseCount_350_450[ich]);}
      if( fAvePulseCount_450_550[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_450_550[ich], fAvePulseCount_450_550[ich]);}
      if( fAvePulseCount_550_650[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_550_650[ich], fAvePulseCount_550_650[ich]);}
      if( fAvePulseCount_650_750[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_650_750[ich], fAvePulseCount_650_750[ich]);}
      if( fAvePulseCount_750_850[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_750_850[ich], fAvePulseCount_750_850[ich]);}
      if( fAvePulseCount_850_950[ich] > 0 ){ NormalizeWfmHistogram(h_AvePulse_850_950[ich], fAvePulseCount_850_950[ich]);}
  
      // do fit
      TF1 f_tail("f_tail","[0]*exp(-x/[1])",400,2000);
      f_tail.SetParameter(0,100);
      //f_tail.SetParLimits(0,1,1e5);
      f_tail.SetParameter(1,1000);
      f_tail.SetParLimits(1,300,1800);
      h_AvePulse[ich]->Fit("f_tail","RQ");
      std::cout
      <<"Tau fit to ave waveform (400ns - 2us) = "<<f_tail.GetParameter(1)<<" +/- "<<f_tail.GetParError(1)<<" ns\n"
      <<"Tau from event-by-event               = "<<h_Tau[ich]->GetMean()<<" +/- "<<h_Tau[ich]->GetRMS()/sqrt(h_Tau[ich]->Integral())<<" ns\n"
      <<"RMS = "<<h_Tau[ich]->GetRMS()<<" ns\n";
    }
    
    // Normalize the summed PE waveform to get an average
    if(fSaveAvePhelWfm && AvePECount[ich] > 0 ) NormalizeWfmHistogram(h_AvePEWfm[ich], AvePECount[ich] );
   
    if( fAnalyzePhelWfm ) {
      // Take a closer look at the averaged PE waveform
      std::cout
      <<"\n"
      <<"---------------------------------------\n"
      <<"Analyzing phel waveform\n";
      int N = fPEWfm_totalSamples[ich];
      std::vector<float> pe_vals(N,0);
      for(int i=1; i<=N; i++) pe_vals.at(i-1) = h_AvePEWfm[ich]->GetBinContent(i);
      int prek = fPEWfm_preWindow[ich];
      int k1 = 0;

      int k2 = prek - fPreIntWindow[ich];
      int k3 = prek + fPostIntWindow[ich];
      int k4 = prek + 100;
      if( k2-k1 >= 5 ) {
        float pe_bs = fOpAlg.GetBaselineAndRMS(pe_vals,k1,k2)[0];
        std::cout<<"Pre PE baseline = "<<pe_bs<<" (samples "<<k1<<"-"<<k2<<")\n";
      
        float integral_spe      = 0.;
        float integral_spe_bs   = 0.;
        float integral_full     = 0.;
        float integral_full_bs  = 0.;
        std::cout<<"Integrating SPE region from "<<k2-prek<<" to "<<k3-prek<<"\n";
        std::cout<<"Integrating full waveform from "<<k2-prek<<" to "<<k4-prek<<"\n\n";
        for(int i=k2; i<k4; i++){
          integral_full     += pe_vals.at(i);
          integral_full_bs  += pe_vals.at(i) - pe_bs;
          if( i >= k2 && i <= k3 ) {
            integral_spe += pe_vals.at(i);
            integral_spe_bs += pe_vals.at(i) - pe_bs;
          }
        }

        float factor = 0.;
        float factor2 = 0.;
        if( integral_full > 0 ) factor = integral_full / integral_spe;
        if( integral_full_bs > 0 ) factor2 = integral_full_bs / integral_spe_bs;

        std::cout
        <<"  SPE region integral   : "<<integral_spe<<"\n"
        <<"  Full integral         : "<<integral_full<<"\n"
        <<"  Corr. factor          = "<<factor<<"\n\n"
        <<"Accounting for wfm baseline...\n"
        <<"  SPE region integral   : "<<integral_spe_bs<<"\n"
        <<"  Full integral         : "<<integral_full_bs<<"\n"
        <<"  Corr. factor          = "<<factor2<<"\n\n"
        <<"Averaged corr. factor   = "<< 0.5*(factor + factor2)<<"\n\n";
      }
    }
   
    
  } // end loop over photodetectors
  std::cout<< "\n============================================================\n";
  
}

//#######################################################################
void OpDetSER::endRun(art::Run const & r){
}

//#######################################################################
void OpDetSER::endSubRun(art::SubRun const & sr){
}

//#######################################################################
std::vector<float> OpDetSER::FitSER(TH1F* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit = true)
{ 
  float max = (float)h_SER->GetEntries();
  std::vector<float> out(17);
  
  if( max == 0 ) {
    std::cout<<"FitSER ERROR: input histogram is empty!\n";
    return out;
  }
 
  // Fit out the SER
  TF1 SER_fit("SER_fit","([0]/([2]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-[1])/[2],2))"
        "+ [3]*((1./(sqrt(1.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-1.*[4])/(sqrt(1)*[5]),2))" 
        "+ ([6]/(sqrt(2.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-2.*[4])/(sqrt(2)*[5]),2))" 
        "+ ([7]/(sqrt(3.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-3.*[4])/(sqrt(3)*[5]),2)))",x1,x2);
  SER_fit.SetNpx(1000);

    SER_fit.SetParName(0,"N_{0}");
    SER_fit.SetParName(1,"#mu_{0}");
    SER_fit.SetParName(2,"#sigma_{0}");
    SER_fit.SetParName(3,"N_{1}");
    SER_fit.SetParName(4,"#mu_{1} (SPE)");
    SER_fit.SetParName(5,"#sigma_{1}");
    SER_fit.SetParName(6,"f_{2}");
    SER_fit.SetParName(7,"f_{3}");

    // "Noise" component (gaus)
    SER_fit.SetParameter(0,max);
    SER_fit.SetParLimits(0,0.,max);
    SER_fit.SetParameter(1,0.);
    SER_fit.SetParLimits(1,-0.80*meanSet,0.80*meanSet);
    SER_fit.SetParameter(2,pedMaxWidth*0.5);
    SER_fit.SetParLimits(2,1.,pedMaxWidth);

    // 1PE (gaus)
    SER_fit.SetParameter(3,max);
    SER_fit.SetParLimits(3,0.,max);
    SER_fit.SetParameter(4,meanSet);
    SER_fit.SetParLimits(4,meanSet*0.7,meanSet*2.);
    SER_fit.SetParameter(5,width);
    SER_fit.SetParLimits(5,5.,width*3.);

    // N2 (ratio)
    SER_fit.SetParameter(6,0.001);
    SER_fit.SetParLimits(6,0.,0.10);

    // N3 (ratio)
    SER_fit.SetParameter(7,0.001);
    SER_fit.SetParLimits(7,0.,0.05);

    // Perform an initial fit, but if goodness isn't good enough, fit with restricted
    // range that excludes range > 1.5 * SER
    h_SER->Fit("SER_fit","QR");
    float n0  = SER_fit.GetParameter(0);
    float dn0 = SER_fit.GetParError(0);
    fSecondFitDone = false;
    if( refit && SER_fit.GetChisquare()/(SER_fit.GetNDF()-1) > fSecondFitChi2Thresh ) {
      fSecondFitDone = true;
      //SER_fit.SetRange(SER_fit.GetParameter(1)+1.0*SER_fit.GetParameter(2), x2);
      //SER_fit.SetRange(SER_fit.GetParameter(1)-2.0*SER_fit.GetParameter(2), SER_fit.GetParameter(4)*2.5);
      SER_fit.SetRange(SER_fit.GetParameter(1)-1.0*SER_fit.GetParameter(2), SER_fit.GetParameter(4)*3.5);
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
    float rn2   = SER_fit.GetParameter(6);
    float drn2  = SER_fit.GetParError(6);
    float n2    = n1*rn2;
    float dn2   = n2*drn2;
    float rn3   = SER_fit.GetParameter(7);
    float drn3  = SER_fit.GetParError(7);
    float n3    = n1*rn3;
    float dn3   = n3*drn3;
    float chi   = SER_fit.GetChisquare()/(SER_fit.GetNDF()-1);

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


//#######################################################################
void OpDetSER::NormalizeWfmHistogram( TH1F* h, int count)
{
  if( count > 0 ) h->Scale(1./float(count));
  h->SetOption("hist");
}

DEFINE_ART_MODULE(OpDetSER)
