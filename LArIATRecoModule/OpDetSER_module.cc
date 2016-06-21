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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <fstream>

// ROOT includes
#include <TF1.h>
#include <TH1.h>
#include <TH2F.h>

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
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Tunable parameters defined by fcl
  std::vector<size_t> fSelectChannels;
  std::vector<size_t> fOpDetChannels;
  std::vector<size_t> fBaselineChannel;
  std::vector<float>  fWfmAbsRMSCut;
  std::vector<float>  fPrePE_RMSFactorCut;
  std::vector<float>  fPulseHitRMSThresh;
  std::vector<float>  fPulseHitThreshHigh;
  std::vector<float>  fPedestalMean_lowerLim;
  std::vector<float>  fPedestalMean_upperLim;
  std::vector<float>  fMean_set;
  std::vector<float>  fMean_lowerLim;
  std::vector<float>  fMean_upperLim;
  std::vector<float>  fSinglePE;
  std::vector<float>  fSinglePE_tolerance;
  std::vector<std::vector<float>>  SERWaveform;
  std::vector<size_t> SERWaveform_count;

  bool                fVerbose;
  std::string         fDAQModule;
  std::string         fInstanceName;
  float               fMvPerADC;
  bool                fAttemptFit;
  bool                fWfmBaselineCorrection;
  float		      fGradientCut;
  short               fBaselineWindowLength;
  short               fT1;
  short               fT2;
  float               fTimestamp_T1;
  float               fTimestamp_T2; 
  short               fPostWindow;
  short               fPreWindow;
  float               fMaxWindowFactor;
  float               fUsePrePEBaseline;
  short               fPrePEBaselineWindow;
  short               fThreshPersist;
  short               fDeadTimeMin;
  short               fNSamples;

  short               SER_bins;
  bool                gotPMT;
  size_t              NSamples;
  std::vector<short>  Wfm;
  std::vector<short>  WfmBaseline;
  float               Timestamp;
  short               TrigSample;

  char		      histName[100];
  char                histTitle[100];
  std::map<size_t,size_t> iCh;

  // Alg objects
  OpHitBuilderAlg     fOpHitBuilderAlg;
  
  // Histograms
  TH1I* h_TotalEvents;
  TH1F* h_SER[10];
  TH1F* h_AvePEWfm[10];
  TH1F* h_WfmRMS[10];
  TH2F* h_GradVsSER[10]; 
  TH1I* h_PhelTime[10];
  TH1I* h_LiveTime[10];

  TH1F* h_HitTimes[10];
  TH1F* h_HitAmplitudes[10];

};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg"))
{
  this                ->reconfigure(p);
  
  // Create the channel index map
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    iCh[fOpDetChannels[i]] = i;
  }

  SER_bins            = fPreWindow + fPostWindow;
  SERWaveform         .resize(10);
  SERWaveform_count   .resize(10);
  for(size_t i=0; i< 10; i++){
    SERWaveform[i]        .resize(SER_bins);
    SERWaveform_count[i]  = 0;
  }

  if(fBaselineChannel.size() == 0 ) fWfmBaselineCorrection = false;

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

  // Skip event if no waveforms found
  if( (size_t)WaveformHandle->size() == 0 ){
    std::cout << "No optical detector data found; skipping event.\n";
    return;
 
  // If waveforms found, look for the specified PMTs
  } else {

    h_TotalEvents->Fill(0);

    /*
    // First, look for the "baseline" input if it exists
    WfmBaseline.clear();
    if( fBaselineChannel.size() > 0 ){
      for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){
        art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
        raw::OpDetPulse ThePulse = *ThePulsePtr;
        if(ThePulse.OpChannel() == fBaselineChannel[0]) WfmBaseline = ThePulse.Waveform();
      }
    }
    */

    // Now do full analysis on each pulse, subtracting the baseline pulse if specified
    for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){

      art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
      raw::OpDetPulse ThePulse = *ThePulsePtr;
      
      size_t ch = ThePulse.OpChannel();
      Timestamp   = (float(ThePulse.PMTFrame())*8.)/1.0e09;
      if( Timestamp < fTimestamp_T1 || Timestamp > fTimestamp_T2 ) continue;

      // Check if this optical channel is among those selected
      // for SER analysis in the fhicl script
      gotPMT = false;
      for( size_t i=0; i<fSelectChannels.size(); i++){
        if( fSelectChannels[i] == ch ) { gotPMT=true; break;}
      }
      if(!gotPMT) continue;

      Wfm       = ThePulse.Waveform();
      NSamples    = Wfm.size();
      //if(fWfmBaselineCorrection && ch != fBaselineChannel[0]){ for(size_t i=0; i<NSamples; i++) Wfm[i] = Wfm[i] - WfmBaseline[i];} // Correct for shared baseline noise
      TrigSample  = (short)ThePulse.FirstSample(); 
      size_t ich  = iCh[ch];
       
      std::cout << "\n"; 
      std::cout << runnr<<"/"<<subrunnr<<"/"<<eventnr<<"  OpDet "<<ch<<", index "<<ich<<": " << NSamples << " samples, T0 "<<TrigSample<<", timestamp "<<Timestamp<<" sec\n";
      
      // Find waveform baseline and RMS 
      std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
      float baseline  = tmp[0];
      float rms       = tmp[1];
      h_WfmRMS[ich]   ->Fill(rms*fMvPerADC);
     
      // Run standard hit-finder over the waveform:
      std::vector<short> hit_times = fOpHitBuilderAlg.GetHits( ThePulse );
      for(size_t i=0; i<hit_times.size(); i++){
        h_HitTimes[ich]      ->Fill(hit_times[i]);
        h_HitAmplitudes[ich] ->Fill((baseline-fOpHitBuilderAlg.GetLocalMinimum(Wfm,hit_times[i]))*fMvPerADC);
      }
      
      short t1 = std::max(TrigSample + fT1,0);
      short t2 = std::min(TrigSample + fT2, (int)NSamples);
      
      std::cout << "Waveform raw baseline: " << baseline << " +/- " << rms <<" ADC ("<<baseline*fMvPerADC<<" +/- "<<rms*fMvPerADC<<" mV)\n";
      if( rms*fMvPerADC >= fWfmAbsRMSCut[ich] ) {
        std::cout<<"Baseline too noisy... moving on.\n";
        continue;
      }

      float   integral = 0;
      bool    flag = false;
      int     windowsize = 0;
      int     counter = 0;
      short   quietTime = 0;
      float   prePE_baseline = -99;
      float   prePE_rms = 99;
  
      // Now begin looking for single PEs
      std::cout << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh[ich] << ", threshPersist " << fThreshPersist <<") \n";

      // Make gradient
      std::vector<float> g = fOpHitBuilderAlg.MakeGradient(Wfm);
      float hit_grad = 0;
      int hit_time = -1;

      // Baseline subtraction and signal inversion
      std::vector<float> wfm_corrected(NSamples);
      for(size_t i=0; i<NSamples; i++) wfm_corrected[i] = baseline - (float)Wfm[i];

      // Make empty vector in which to store waveform for each PE candidate
      // to be reset after each PE (or after added to SERWaveform)
      std::vector<float> tmp_wfm(SER_bins);
      int tmp_wfm_i = 0;

      
      // ------------------------------------------------------------
      // Scan waveform and search for hits within specified time window
      for(short i = t1; i < t2; i++){
       
        float y  = wfm_corrected[i];
        float yy = y*fMvPerADC;
        bool  IsLive        = (quietTime >= fDeadTimeMin);
        bool  IsOverThresh  = (y  >= fPulseHitRMSThresh[ich]*rms);
        bool  IsOverLimit   = (yy >= fPulseHitThreshHigh[ich]); 
        bool  IsPECandidate = ((IsOverThresh) && (!IsOverLimit) && ( fabs(g[i]) >= fabs(fGradientCut)) );
        
        h_LiveTime[ich]->Fill((int)(IsLive || flag));
    
        // If we're already in a PE window, increment the
        // counters and add to integral
        if( flag ) {
         
          //std::cout
          //<< "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
          //<< " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << " mV, g " << g[i]<<", deadTime = "<<quietTime<<"\n";
          
          counter++;
          windowsize++;
          float yc = wfm_corrected[i] - fUsePrePEBaseline*prePE_baseline;
          integral += yc;
    
          if(tmp_wfm_i < SER_bins){
            tmp_wfm[tmp_wfm_i] = yc;
            tmp_wfm_i++;
          }
          
          // If another PE is detected after at least 5 ns, extend the window by resetting counter
          if( (counter >= 5 + fPreWindow)  && IsPECandidate ){
            //std::cout << "  Secondary hit, extending window\n";
            counter = 0;
          }
          
          // If pulse extends above upper limit or if window length
          // is exceeded due to multiple merges, abort mission.
          if( IsOverLimit || (windowsize > fPreWindow + fMaxWindowFactor*fPostWindow) ){
            std::cout << "  abort!\n";
            counter = 0;
            hit_grad = 0;
            hit_time = -1;
            integral = 0;
            flag = false;
            windowsize = 0;
            quietTime = 0;
            tmp_wfm.clear();
            tmp_wfm_i = 0;
            continue;
          }
          
          // If we reached the end of the allotted window, add
          // integral to vector and reset everything
          if( counter == fPostWindow ){
            
            std::cout
            << "Finished PE window of size "<<windowsize<<", tmp_wfm_i = "<<tmp_wfm_i << "  "
            << integral << " ADCs";
            
            h_SER[ich]        ->Fill(integral);
            h_GradVsSER[ich]  ->Fill(hit_grad,integral);
            
            // If this is > 20 ADC, plot the time:
            if(integral > 20 ) h_PhelTime[ich]   ->Fill(hit_time);
            
            // Add to average waveform if it looks good
            if( (windowsize == SER_bins) && fabs(integral - fSinglePE[ich])<=fSinglePE_tolerance[ich] ){
              std::cout << "Add to average PE wfm.\n";
              for(short ii=0; ii<SER_bins; ii++) SERWaveform[ich].at(ii) += tmp_wfm[ii]*fMvPerADC;
              SERWaveform_count[ich]++;
            }
    
            integral = 0;
            hit_grad = 0;
            hit_time = -1;
            windowsize = 0;
            counter = 0;
            quietTime = 0;
            flag = false;
            tmp_wfm.clear();
            tmp_wfm_i=0;
            continue;
          }
        
        } // <-- end if(flag)
   
        
        
        if( !IsOverThresh ) quietTime++;
    
        
        
        // If we're not yet integrating a PE window, signal is within bounds,
        // and enough dead-time (quietTime) has elapsed, then we're in business
        if( !flag && IsPECandidate && IsLive ){
          
          // Find pre-PE baseline
          prePE_baseline = 99;
          prePE_rms      = 99;
          //std::vector<float> tmp = fOpHitBuilderAlg.GetPedestalAndRMS(wfm_corrected,i-buffer_min,i);
          std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS(wfm_corrected,i-fPrePEBaselineWindow,i);
          prePE_baseline = tmp[0];
          prePE_rms      = tmp[1];
         
          //std::cout
          //<< "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
          //<< " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag << "\n"
          //<< "  Potential PE!  preBS = "<< prePE_baseline*fMvPerADC << " mV, preRMS = "<< prePE_rms*fMvPerADC <<" mV, quietTime count "<<quietTime<<"\n";
     
          // Look a few samples ahead and make sure the signal
          // remains above threshold for some time...
          bool flagger = true;
          if( fThreshPersist > 0 ){
            //std::cout<< "  Peeking at next few samples....\n";
            for( short j=0; j<fThreshPersist; j++ ){
              //std::cout<<"    "<<j<<"  "<<wfm_corrected[i+j]*fMvPerADC << " mV \n";
              if ( wfm_corrected[i+j] < fPulseHitRMSThresh[ich]*rms ) flagger = false;
            }
          }
          // Require flat pre-PE region and that threshold persist requirement was met
          if( (prePE_rms <= fPrePE_RMSFactorCut[ich]*rms) && (flagger)){
            
            // Found a "PE hit", so start integral by
            // adding preceeding prepulse region
            flag      = true;
            hit_grad  = g[i];
            hit_time  = i;

    
            // Add up previous "prewindow" samples
            for(short ii=0; ii <= fPreWindow; ii++)
            {
              integral += wfm_corrected[i-fPreWindow+ii] - fUsePrePEBaseline*prePE_baseline;
              tmp_wfm[tmp_wfm_i] = wfm_corrected[i-fPreWindow+ii] - prePE_baseline;
              tmp_wfm_i++;
              windowsize++;
              counter = 1;
            }
          
            std::cout << "** Candidate at "<<hit_time<<", "<<wfm_corrected[i]*fMvPerADC<<" mV, g = "<<hit_grad<<".  Beginning integration...\n"; 
          } else {
            //std::cout << "  Doesn't pass quality cuts, moving on...\n"; 
          }
    
        } // <-- end if(PE cand)
        
        if( IsOverThresh ) quietTime = 0;
    
      } // <-- end scan over waveform 


    } // endLoop over OpDets
  } // endIf Waveformhandle is valid
  


}

//#######################################################################
void OpDetSER::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
 
  h_TotalEvents   = tfs->make<TH1I>("TotalEvents","Total events",1,0,1);
  
  // Create histograms
  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ch   = fSelectChannels[i];
    size_t ich  = iCh[ch]; 
    
    sprintf(histName,"%lu_SER",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integraded ADC;Counts",ch);
    h_SER[ich]       = tfs->make<TH1F>(histName,histTitle,400,-50.,350.);
   
    sprintf(histName,"%lu_AvePEWfm",ch); 
    sprintf(histTitle,"Average PE waveform for OpDet %lu;ns;mV",ch);
    h_AvePEWfm[ich]  = tfs->make<TH1F>(histName,histTitle,SER_bins,0.,(float)SER_bins);
    
    sprintf(histName,"%lu_RMS",ch); 
    sprintf(histTitle,"RMS for OpDet %lu;Baseline RMS [mV];Counts",ch);
    h_WfmRMS[ich]    = tfs->make<TH1F>(histName,histTitle,100,0.0,0.5);
    
    sprintf(histName,"%lu_GradVsSER",ch); 
    sprintf(histTitle,"Grad vs. SER for OpDet %lu;Gradient;SER [ADC]",ch);
    h_GradVsSER[ich]    = tfs->make<TH2F>(histName,histTitle,30,-10,5,  100,-50,350);
    h_GradVsSER[ich]    ->SetOption("colz");

    sprintf(histName,"%lu_PhelTime",ch); 
    sprintf(histTitle,"PE times for OpDet %lu;Sample [ns];Counts",ch);
    h_PhelTime[ich]    = tfs->make<TH1I>(histName,histTitle,(int)(fNSamples/100.),0,fNSamples);
    
    sprintf(histName,"%lu_LiveTime",ch); 
    sprintf(histTitle,"SER-finder live-time (1) vs. dead-time (0) for OpDet %lu;;total samples",ch);
    h_LiveTime[ich]    = tfs->make<TH1I>(histName,histTitle,2,0,2);
    
    sprintf(histName,"%lu_HitTimes",ch); 
    sprintf(histTitle,"Found hit times for OpDet %lu;Samples [ns];Counts",ch);
    h_HitTimes[ich]    = tfs->make<TH1F>(histName,histTitle,(int)(fNSamples/100.),0,fNSamples);
    
    sprintf(histName,"%lu_HitAmplitudes",ch); 
    sprintf(histTitle,"Found hit amplitudes for OpDet %lu;Amplitude [mV];Counts",ch);
    h_HitAmplitudes[ich]    = tfs->make<TH1F>(histName,histTitle,100,0,10);
    

  }
}

void OpDetSER::reconfigure(fhicl::ParameterSet const & p)
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
  std::vector<float>      SinglePEToleranceDefault{ 5, 5};
  std::vector<float>      WfmAbsRMSCutDefault{ 0.5, 0.5};
  std::vector<float>      PedestalMeanLowerLimDefault{ -10, -10};
  std::vector<float>      PedestalMeanUpperLimDefault{ 10, 10};

  fSelectChannels         = p.get< std::vector<size_t> >  ("SelectChannels",SelectChannelsDefault);
  fOpDetChannels          = p.get< std::vector<size_t> >  ("OpDetChannels", SelectChannelsDefault);
  fBaselineChannel        = p.get< std::vector<size_t> >  ("BaselineChannel",BaselineChannelDefault);
  fMean_set               = p.get< std::vector<float> >   ("Mean_set",MeanSetDefault);
  fMean_lowerLim          = p.get< std::vector<float> >   ("Mean_lowerLim",MeanLowerLimDefault);
  fMean_upperLim          = p.get< std::vector<float> >   ("Mean_upperLim",MeanUpperLimDefault);
  fPrePE_RMSFactorCut     = p.get< std::vector<float> >   ("PrePE_RMSFactorCut",PrePeRMSFactorCutDefault);
  fGradientCut            = p.get< float >   ("GradientCut",0);
  fPulseHitThreshHigh     = p.get< std::vector<float> >   ("PulseHitThresh_high",PulseHitThreshHighDefault);
  fPulseHitRMSThresh      = p.get< std::vector<float> >   ("PulseHitRMSThresh",PulseHitRMSThreshDefault);
  fSinglePE               = p.get< std::vector<float> >   ("SinglePE",SinglePEDefault);
  fSinglePE_tolerance     = p.get< std::vector<float> >   ("SinglePE_tolerance",SinglePEToleranceDefault);
  fWfmAbsRMSCut           = p.get< std::vector<float> >   ("WfmAbsRMSCut", WfmAbsRMSCutDefault);
  fPedestalMean_lowerLim  = p.get< std::vector<float> >   ("PedestalMean_lowerLim",PedestalMeanLowerLimDefault);
  fPedestalMean_upperLim  = p.get< std::vector<float> >   ("PedestalMean_upperLim",PedestalMeanUpperLimDefault);
  
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fAttemptFit             = p.get< bool >         ("AttemptFit","true");
  fWfmBaselineCorrection  = p.get< bool >         ("WfmBaselineCorrection","false");
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fTimestamp_T1           = p.get< float >        ("Timestamp_T1",0);
  fTimestamp_T2           = p.get< float >        ("Timestamp_T2",60);
  fPreWindow              = p.get< short >       ("PreWindow",5);
  fPostWindow             = p.get< short >       ("PostWindow",45);
  fMaxWindowFactor        = p.get< float >        ("MaxWindowFactor",2);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.2);
  fUsePrePEBaseline       = p.get< float >        ("UsePrePEBaseline",1);
  fPrePEBaselineWindow    = p.get< short >        ("PrePEBaselineWindow",100);
  fThreshPersist          = p.get< short >        ("ThreshPersist",3);
  fDeadTimeMin            = p.get< short >        ("DeadTimeMin",100);
  fVerbose                = p.get< bool >         ("Verbose",true);
  fNSamples               = p.get< short >        ("NSamples",28672);
}

void OpDetSER::beginRun(art::Run const & r){}

void OpDetSER::beginSubRun(art::SubRun const & sr){}

void OpDetSER::endJob()
{

  std::cout 
    << "============================================================\n"
    << "Ending SER program.\n"
    << "  Use pre-PE BS in int?   " << fUsePrePEBaseline << "\n"
    << "  Pre-PE baseline length  " << fPrePEBaselineWindow << "\n"
    << "  Threshhold persist      " << fThreshPersist <<"\n" 
    << "  Graident cut            " << fGradientCut << "\n"
    << "  Pre / Post Window       " << fPreWindow << "," << fPostWindow << "\n\n";

  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ich = iCh[fSelectChannels[i]];
    std::cout
    << "-------------------------------------------\n"
    << "  channel                 " << fOpDetChannels[ich] << "\n"
    << "  Abs Wfm RMS cut         " << fWfmAbsRMSCut[ich] << " mV\n"
    << "  Pre-PE RMS cut          " << fPrePE_RMSFactorCut[ich] << " x wfm RMS\n"
    << "  RMS thresh factor       " << fPulseHitRMSThresh[ich] <<" x wfm RMS\n"
    << "  Mean baseline RMS	  " << h_WfmRMS[ich]->GetMean(1) << "\n"
    << "  --> PE candidates found " << h_SER[ich]->GetEntries() << "\n\n";
    if( SERWaveform_count[ich] > 0 ){
      std::cout
      << "  Ave PE waveform (sample [ns], ADC):\n";
      
      float integral = 0;
      for( int i = 0; i < SER_bins; i++) {
        float w = SERWaveform[ich].at(i) / float(SERWaveform_count[ich]); 
        h_AvePEWfm[ich]->Fill(i,w);
        integral += w/fMvPerADC;
        std::cout<<"    "<<i<<"     "<<w<<"\n";
      }
    
      std::cout
      << "  Integral: "<<integral<< " ADC ("<<SERWaveform_count[ich]<<" waveforms averaged)\n\n";
    
    }
  
  
 
 
    if( fAttemptFit ){
      std::cout
      <<"Attempting SER fit...\n";
       
      // Fit out the SER
      TF1 f_gaus("f_gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))",-1000,1000);
      //TF1 SER_fit("SER_fit","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/[5],2)) + [6]*exp(-0.5*pow((x-2.*[4])/[7],2)) + [8]*exp(-0.5*pow((x-3.*[4])/[9],2))",-50,400);
      TF1 SER_fit("SER_fit","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/(sqrt(1)*[5]),2)) + [6]*exp(-0.5*pow((x-2.*[4])/(sqrt(2)*[5]),2)) + [7]*exp(-0.5*pow((x-3.*[4])/(sqrt(3)*[5]),2))",-50,400);
  
      float max = (float)h_SER[ich]->GetMaximum();
  
      SER_fit.SetParName(0,"Noise norm");
      SER_fit.SetParName(1,"Noise mean ADC");
      SER_fit.SetParName(2,"Noise sigma");
      SER_fit.SetParName(3,"1PE norm");
      SER_fit.SetParName(4,"1PE mean ADC");
      SER_fit.SetParName(5,"1PE sigma");
      SER_fit.SetParName(6,"2PE norm");
      SER_fit.SetParName(7,"3PE norm");
  
      // "Noise" component (gaus)
      SER_fit.SetParameter(0,max);
      SER_fit.SetParLimits(0,0.,1.2*max);
      SER_fit.SetParameter(1,0);
      SER_fit.SetParLimits(1,fPedestalMean_lowerLim[ich],fPedestalMean_upperLim[ich]);
      SER_fit.SetParameter(2,20);
      SER_fit.SetParLimits(2,0.,30.);

      // 1PE (gaus)
      SER_fit.SetParameter(3,max);
      SER_fit.SetParLimits(3,0.,1.2*max);
      SER_fit.SetParameter(4,fMean_set[ich]);
      SER_fit.SetParLimits(4,fMean_lowerLim[ich],fMean_upperLim[ich]);
      SER_fit.SetParameter(5,35);
      SER_fit.SetParLimits(5,10.,60.);

      // 2PE (gaus)
      SER_fit.SetParameter(6,0.1*max);
      SER_fit.SetParLimits(6,0.,0.3*max);

      // 3PE (gaus)
      SER_fit.SetParameter(7,1.);
      SER_fit.SetParLimits(7,0.,0.3*max);

      h_SER[ich]->Fit("SER_fit","R");

      /*
      ofstream file;
      file.open("SER_result.txt");
      file <<"=========================================================\n";
      file <<"SER fit results:\n";
      file <<"  1PE = "<<SER_fit.GetParameter(4)<<" +/- "<<SER_fit.GetParError(4)<<"\n";
      file <<"  sigma = "<<SER_fit.GetParameter(5)<<"\n";
      file <<"  red Chi2 = "<<SER_fit.GetChisquare()/(SER_fit.GetNDF()-1)<<"\n";
      file <<"=========================================================\n";
      file.close(); 
      */

      std::cout<<"SER fit results for OpDet "<< fOpDetChannels[ich]<<":\n";
      std::cout<<"  1PE = "<<SER_fit.GetParameter(4)<<" +/- "<<SER_fit.GetParError(4)<<"\n";
      std::cout<<"  sigma = "<<SER_fit.GetParameter(5)<<"\n";
      std::cout<<"  red Chi2 = "<<SER_fit.GetChisquare()/(SER_fit.GetNDF()-1)<<"\n";

      // Draw the components
      //TF1 f_noise("f_noise","[0]*exp(-0.5*pow((x-[1])/[2],2))",-20,350);
      //f_noise.SetParameter(0,SER_fit.GetParameter(0));
      //f_noise.SetParameter(1,SER_fit.GetParameter(1));
      //f_noise.SetParameter(2,SER_fit.GetParameter(2));
      //f_noise.SetLineColor(kBlack);

      //std::cout<<"f_noise norm: "<<f_noise.GetParameter(0);
    }
  
    std::cout<<"=====================================================\n";
 
  }
  
}

void OpDetSER::endRun(art::Run const & r){}

void OpDetSER::endSubRun(art::SubRun const & sr){}

DEFINE_ART_MODULE(OpDetSER)
