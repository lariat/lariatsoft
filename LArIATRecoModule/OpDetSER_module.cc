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

// ROOT includes
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
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
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) ;
  std::vector<float> FitSER(TH1D* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit);

private:
  
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory waveformDir = tfs->mkdir("Waveforms");
  
  // Alg object
  OpHitBuilderAlg     fOpHitBuilderAlg;
  
  // Counters and things
  size_t fNumberActiveEvents[10];
  size_t fNumberSaturatedEvents[10];
  size_t LiveSamples[10];
  size_t WaveformCount[10];
  size_t AvePECount[10];
  
  // Timetamp-related variables  
  double fTimestamp;
  double fPrevTimestamp;
  double fTimestampOffset;
  double fModTimestamp;
  
  // Histograms
  TH1D* h_Timestamps;
  TH1D* h_SER[10];
  TH1D* h_SERWindowSize[10];
  TH1D* h_AvePEWfm[10];
  TH1D* h_RawBaseline[10];
  TH1D* h_BaselineRMS[10];
  TH1D* h_PrePhelRMS[10];
  TH1D* h_Amplitudes[10];
  TH1D* h_Waveform;
  TH1D* h_Waveform_raw;
  TH1D* h_WaveformSER;

  // Tunable parameters defined by fcl
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
  std::vector<std::vector<float>>  SERWaveform;
  bool		      fSaveAvePhelWfm;
  bool                fDoSecondFit;
  size_t              fMaxSavedWaveforms;
  size_t              fNSamples;
  size_t              fRunningBaselineLength;
  float               fSecondFitChi2Thresh;
  bool                fSubtractRunningBaseline;
  bool                fVerbose;
  float               fMvPerADC;
  short               fBaselineWindowLength;
  short               fT1;
  short               fT2;
  float               fTimestamp_T1;
  float               fTimestamp_T2;
  float               fMaxWindowFactor;
  short               fPrePEBaselineWindow;
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

  std::map<size_t,size_t> iCh;
  TCanvas * c;

};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
{
  this                ->reconfigure(p);
  SERWaveform.resize(fOpDetChannels.size());
  
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    iCh[fOpDetChannels[i]]    = i; // map from channel to index
    SER_bins[i]               = fPreIntWindow[i]+fPostIntWindow[i];
    SERWaveform[i]            .resize(SER_bins[i]);
    WaveformCount[i]          = 0;
    LiveSamples[i]            = 0;
    fNumberSaturatedEvents[i] = 0;
    fNumberActiveEvents[i]    = 0;
  }
}

//#######################################################################
void OpDetSER::analyze(art::Event const & e)
{
  int eventnr   = e.id().event();
  int runnr     = e.run();
  int subrunnr  = e.subRun();
 
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
  // correct for funny bug in timestamp assignment when we collect beyond 17.18 sec
  if( (fPrevTimestamp - fTimestamp) > 1. ) fTimestampOffset += 17.18;
  fModTimestamp = fTimestamp + fTimestampOffset;
  fPrevTimestamp = fTimestamp; 
  if (  !fOpHitBuilderAlg.eventTypeFilter(fModTimestamp,fSelectEventTypes)
      ||  !(fModTimestamp >= fTimestamp_T1 || fModTimestamp <= fTimestamp_T2)     ){
    if(fVerbose) std::cout<<"Timestamp ("<<fModTimestamp<<" sec) out of range --> skipping event\n";
    return;
  }
  
  h_Timestamps->Fill(fModTimestamp);

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


    // Reset live samples counter
    size_t liveSamples = 0;

    // Get the raw waveform
    std::vector<short>  Wfm_raw   = ThePulse.Waveform();
    size_t              NSamples  = Wfm_raw.size();
    
    // Convert this to a vector<float> for easier manipulation later on
    std::vector<float> Wfm( Wfm_raw.begin(), Wfm_raw.end() );
   
    // Subtract running baseline to remove oscillations (useful for random pulser runs where 
    // the signal is expected to remain approximately at baseline).
    if(fSubtractRunningBaseline) {
      if(fVerbose) std::cout<<"Performing running baseline subtraction...\n";
      fOpHitBuilderAlg.SubtractRunningBaseline(Wfm, Wfm, NSamples, fRunningBaselineLength);
    }
    
    // Find waveform baseline and RMS
    std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
    float baseline  = tmp[0];
    float rms       = tmp[1];
    h_RawBaseline[ich]    ->Fill(baseline);
    h_BaselineRMS[ich]    ->Fill(rms*fMvPerADC);
    
    // Skip if RMS is too noisy
    if( rms*fMvPerADC >= fWfmAbsRMSCut[ich] ) {
      if(fVerbose) std::cout<<"Baseline too noisy on channel "<<ich<<" ("<<rms*fMvPerADC<<" mV)... skipping this event.\n";
      return;
    } 
    
    if(fVerbose){ 
      std::cout << "OpDet "<<ch<<", index "<<ich<<" Timestamp "<<fModTimestamp<<" sec\n";
      std::cout << "Waveform raw baseline: " << baseline << " +/- " << rms <<" ADC\n";
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

      h_WaveformSER   = new TH1D("serwindow",histTitle,Wfm.size(),0,Wfm.size());
      h_WaveformSER   ->SetLineColor(kRed); 
      for(short ii=0; ii< short(Wfm.size()); ii++) h_WaveformSER     ->Fill(ii,0.00);
        
      areSavingWaveform = true;
    }
      
    // Baseline subtraction and signal inversion.
    std::vector<float> wfm_corrected(NSamples);
    float amp = -1.;
    bool isSaturated = false;
    for(size_t i=0; i<NSamples; i++) {
      wfm_corrected[i] = fOpDetPolarity[ich]*((float)Wfm[i]-baseline);
      if(areSavingWaveform) {
	h_Waveform_raw    ->Fill(i,(float)Wfm_raw[i]);
	h_Waveform	  ->Fill(i,wfm_corrected[i]*fMvPerADC);
      }
      if( wfm_corrected[i] > amp  ) amp = wfm_corrected[i];
      if( Wfm_raw[i] == 0         ) isSaturated = true;
    }
    h_Amplitudes[ich]	->Fill(amp*fMvPerADC);
    if( amp > 20.*rms ) {
      fNumberActiveEvents[ich]++;
      if( isSaturated )	fNumberSaturatedEvents[ich]++;
    }
     
      
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
      bool  IsPECandidate   = ((IsOverThresh) && (!IsOverLimit) );
      
      if(IsLive || flag) {
        liveSamples++;
      }
        
      // If we're already in a PE window, increment the
      // counters and add to integral
      if( flag ) {
        
        if(fVerbose) {
        std::cout
        << "  " << i << "  y_mV = " << y_mV << " mV (window size " <<windowsize<< "), "
        << " thresh " << fPulseHitRMSThresh[ich]*rms*fMvPerADC << " mV, quietTime = "<<quietTime<<"   deadTime = "<<deadTime<<"\n";
        }

        counter++;
        windowsize++;
        integral += y;
       
        // Scan for amplitude
        if( y_mV > PE_amp ) PE_amp = y_mV;

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
          h_SERWindowSize[ich]->Fill(windowsize);
          
          // Add to average waveform if it looks good
          if( fSaveAvePhelWfm && windowsize == SER_bins[ich] && fabs(integral - fSinglePE[ich])/fSinglePE[ich] <= fSinglePE_tolerance[ich] ){
            if(fVerbose) std::cout<<"    --> Adding to average PE waveform.\n";
            for(short ii=0; ii<SER_bins[ich]; ii++) SERWaveform[ich].at(ii) += tmp_wfm[ii]*fMvPerADC;
            AvePECount[ich]++;
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
      if( !IsOverThresh && !IsOverThreshNeg ) { 
        deadTime++; 
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
        h_WaveformSER->Draw("hist same");
        c->Write(histName);
    }

    LiveSamples[ich]         += liveSamples;
    if(fVerbose) std::cout<<"Live samples for this event... "<<liveSamples<<"\n"; 
    
  } // endLoop over OpDets
    
}

//#######################################################################
void OpDetSER::beginJob()
{
  fPrevTimestamp   = 0;
  fTimestampOffset = 0;
  fModTimestamp    = 0;
  
  h_Timestamps  = tfs->make<TH1D>("Timestamps","Timestamp for all events;sec",240,0.,60.);

  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ch   = fSelectChannels[i];
    size_t ich  = iCh[ch]; 
    int x_bins  = std::min((int)(fSER_x2[ich] - fSER_x1[ich]),500);
    
    sprintf(histName,"%lu_SER",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integrated ADC;Counts",ch);
    h_SER[ich]		  = tfs->make<TH1D>(histName,histTitle,x_bins,fSER_x1[ich],fSER_x2[ich]);
    
    sprintf(histName,"%lu_SERWindowSize",ch); 
    sprintf(histTitle,"SER window sizes for OpDet %lu;Window size (samples);Counts",ch);
    h_SERWindowSize[ich]       = tfs->make<TH1D>(histName,histTitle,100,0.,fMaxWindowFactor*SER_bins[ich]+10);
  
    sprintf(histName,"%lu_RawBaseline",ch); 
    sprintf(histTitle,"Raw baselines for OpDet %lu;Baseline [ADC];Counts",ch);
    h_RawBaseline[ich]    = tfs->make<TH1D>(histName,histTitle,600,0.0,1200.);
    
    sprintf(histName,"%lu_BaselineRMS",ch); 
    sprintf(histTitle,"RMS for OpDet %lu;Baseline RMS [mV];Counts",ch);
    h_BaselineRMS[ich]    = tfs->make<TH1D>(histName,histTitle,300,0.0,1.5);
    
    sprintf(histName,"%lu_PrePhelRMS",ch); 
    sprintf(histTitle,"RMS of pre-PE candidate baselines, OpDet %lu;Pre-PE RMS [mV];Counts",ch);
    h_PrePhelRMS[ich]     = tfs->make<TH1D>(histName,histTitle,300,0.0,1.5);
    
    sprintf(histName,"%lu_Amplitudes",ch); 
    sprintf(histTitle,"Waveform amplitudes for OpDet %lu;Amplitude [mV];Counts",ch);
    h_Amplitudes[ich]    = tfs->make<TH1D>(histName,histTitle,2000.,0.,200.);
    
    if(fSaveAvePhelWfm){ 
      sprintf(histName,"%lu_AvePEWfm",ch); 
      sprintf(histTitle,"Average PE waveform for OpDet %lu;ns;mV",ch);
      h_AvePEWfm[ich]  = tfs->make<TH1D>(histName,histTitle,SER_bins[ich],0.,(float)SER_bins[ich]);
    }
    
  }
}

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

  fSaveAvePhelWfm	  = p.get< bool >	  ("SaveAvePhelWfm",true);
  fDoSecondFit            = p.get< bool >         ("DoSecondFit",true);
  fSecondFitChi2Thresh    = p.get< float >        ("SecondFitChi2Thresh",1.2);
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fTimestamp_T1           = p.get< float >        ("Timestamp_T1",0);
  fTimestamp_T2           = p.get< float >        ("Timestamp_T2",60);
  fMaxWindowFactor        = p.get< float >        ("MaxWindowFactor",2);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.2);
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
}

void OpDetSER::beginRun(art::Run const & r){
 
}

void OpDetSER::beginSubRun(art::SubRun const & sr){
  fTimestampOffset = 0;
  fPrevTimestamp   = 0;
  fModTimestamp    = 0;
}

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
    << "  Pre-PE baseline length    " << fPrePEBaselineWindow << "\n"
    << "  RMS thresh factor	    " << fPulseHitRMSThresh[ich] <<" x wfm RMS\n"
    << "  Pre-PE RMS cut	    " << fPrePE_RMSFactorCut[ich] << " x wfm RMS\n"
    << "  Threshhold persist	    " << fThreshPersist <<" (x "<<fThreshPersistFactor<<")\n"
    << "  Quiet time min            " << fQuietTimeMin <<"\n"
    << "  Dead time min             " << fDeadTimeMin << "\n"
    << "  Pre / Post window size    " << fPreIntWindow[ich] << "," << fPostIntWindow[ich] << "\n"
    << "  Set single PE             " << fSinglePE[ich] << " integrated ADC\n"
    << "  Active events	(>20*RMS)   " << fNumberActiveEvents[ich] << "\n"
    << "  Saturation rate	    " << (float)fNumberSaturatedEvents[ich] / (float)fNumberActiveEvents[ich] << "\n"
    << "  Mean baseline RMS	    " << h_BaselineRMS[ich]->GetMean(1) << " mV\n"
    << "  --> PE candidates found   " << h_SER[ich]->GetEntries() << "\n"
    <<"-----------------------------------------\n"
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
    <<"  live time       = "<<cumulativeTime<<" sec\n";
    
    // Normalize the summed PE waveform to get an average
    if(fSaveAvePhelWfm && AvePECount[ich] > 0 ){
      if( fVerbose ) {
        std::cout
        <<"-----------------------------------------\n"
        << "Ave PE waveform for OpDet "<<fOpDetChannels[ich]<<" ("<<AvePECount[ich]<<" entries)\n"
        << "Selection criteria: "<<fSinglePE[ich]<<" ADC +/- "<<fSinglePE_tolerance[ich]*100.<<"%\n" 
        << "Format:  sample [ns]   ADC\n";
      }
      for( int i = 0; i < SER_bins[ich]; i++) {
        float w = SERWaveform[ich].at(i) / float(AvePECount[ich]); 
        h_AvePEWfm[ich]->Fill(i,w);
        if( fVerbose) std::cout<<"  "<<i<<"     "<<w/fMvPerADC<<"\n";
      }
      h_AvePEWfm[ich]->SetOption("HIST");
    }
    
  } // end loop over photodetectors
  
}

//#######################################################################
void OpDetSER::endRun(art::Run const & r){
}

//#######################################################################
void OpDetSER::endSubRun(art::SubRun const & sr){
}

//#######################################################################
std::vector<float> OpDetSER::FitSER(TH1D* h_SER, float x1, float x2, float meanSet, float width, float pedMaxWidth, bool refit = true)
{ 

  // Fit out the SER
  TF1 SER_fit("SER_fit","([0]/([2]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-[1])/[2],2))"
        "+ [3]*((1./(sqrt(1.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-1.*[4])/(sqrt(1)*[5]),2))" 
        "+ ([6]/(sqrt(2.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-2.*[4])/(sqrt(2)*[5]),2))" 
        "+ ([7]/(sqrt(3.)*[5]*sqrt(2.*TMath::Pi())))*exp(-0.5*pow((x-3.*[4])/(sqrt(3)*[5]),2)))",x1,x2);
 
  float max = (float)h_SER->GetEntries();

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
      SER_fit.SetRange(SER_fit.GetParameter(1)+1.0*SER_fit.GetParameter(2), x2);
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

DEFINE_ART_MODULE(OpDetSER)
