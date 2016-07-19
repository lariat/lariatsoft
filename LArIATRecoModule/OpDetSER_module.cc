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
#include "LArIATRecoAlg/TriggerFilterAlg.h" 

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

  bool                fVerbose;
  std::string         fDAQModule;
  std::string         fInstanceName;
  float               fMvPerADC;
  bool                fAttemptFit;
  float		      fGradientCut;
  short               fBaselineWindowLength;
  short               fT1;
  short               fT2;
  float               fTimestamp_T1;
  float               fTimestamp_T2; 
  short               fPostWindow;
  short               fPreWindow;
  float               fMaxWindowFactor;
  short               fPrePEBaselineWindow;
  short               fThreshPersist;
  short               fDeadTimeMin;
  short               fQuietTimeMin;
  short               fNSamples;

  short               SER_bins;
  float               Timestamp;

  char		      histName[100];
  char                histTitle[100];
  std::map<size_t,size_t> iCh;

  // Alg objects
  OpHitBuilderAlg     fOpHitBuilderAlg;
  TriggerFilterAlg    fTriggerFilterAlg;
  std::string         fTriggerUtility;

  
  // Histograms
  TH1I* h_TotalEvents;
  TH1F* h_SER[10];
  TH1F* h_SERAmp[10];
  TH1F* h_AvePEWfm[10];
  TH1I* h_AvePECount[10];
  TH1F* h_WfmRMS[10];
  TH1F* h_PrePhelRMS[10];
  TH1I* h_PhelTime[10];

  int   NTrigs_Pulser;
  int   NTrigs_Beam;
  int   NTrigs_Michel;
  int   currentSubrun;
  TH1I* h_NTrigs;
  TH1I* h_NTrigs_Pulser;
  TH1I* h_NTrigs_Beam;
  TH1I* h_NTrigs_Michel;
  TH1I* h_TriggersSeen;
  TH2F* h_Baselines;

};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p.get<fhicl::ParameterSet>("OpHitBuilderAlg")),
fTriggerFilterAlg(p.get<fhicl::ParameterSet>("TriggerFilterAlg"))
{
  this                ->reconfigure(p);
  
  // Create the channel index map
  for(size_t i=0; i< fOpDetChannels.size();i++) {
    iCh[fOpDetChannels[i]] = i;
  }

  SER_bins            = fPreWindow + fPostWindow;
  SERWaveform         .resize(10);
  for(size_t i=0; i< 10; i++) SERWaveform[i].resize(SER_bins);

  currentSubrun = -1;

}

//#######################################################################
void OpDetSER::analyze(art::Event const & e)
{
  std::cout<<"\n";
  int eventnr   = e.id().event();
  int runnr     = e.run();
  int subrunnr  = e.subRun();
  

  // This section temporarily used to investigate raw::Trigger objects. 
  // Commented out for now.  Will remove eventually...
  /*
  if( subrunnr != currentSubrun ){
    
    if( subrunnr != -1 ){
      h_NTrigs_Pulser   ->Fill(NTrigs_Pulser);
      h_NTrigs_Beam     ->Fill(NTrigs_Beam);
      h_NTrigs_Michel   ->Fill(NTrigs_Michel);
    }
    
    currentSubrun   = subrunnr;
    NTrigs_Pulser = 0;
    NTrigs_Beam     = 0;
    NTrigs_Michel   = 0;
  }

  art::Handle< std::vector< raw::Trigger >> triggerHandle;
  e.getByLabel(fDAQModule,fInstanceName,triggerHandle);
  if( triggerHandle->size() != 0 ){
    std::cout<<"Trigger handle not empty! "<<triggerHandle->size()<<"\n";
    //art::Ptr< raw::Trigger > TheTriggerPtr(triggerHandle,0);
    //raw::Trigger thisTrigger = *TheTriggerPtr;
    raw::Trigger thisTrigger = triggerHandle->at(0);
    NTrigs_Pulser     += thisTrigger.Triggered(9);
    NTrigs_Beam  += thisTrigger.Triggered(3);
    NTrigs_Michel  += thisTrigger.Triggered(13);

    std::cout<<"thisTrigger.Triggered(9)    = "<<thisTrigger.Triggered(9)<<"\n";
    std::cout<<"thisTrigger.Triggered(3)    = "<<thisTrigger.Triggered(3)<<"\n";
    std::cout<<"thisTrigger.Triggered(13)   = "<<thisTrigger.Triggered(13)<<"\n";
  
    int thereAreTriggers = 0; 
    for(int i=0; i<32; i++ ){
      //std::cout<<"  bit "<<i<<"   "<<thisTrigger.Triggered(i)<<"\n";
      h_NTrigs->Fill(i,thisTrigger.Triggered(i));
      if(thisTrigger.Triggered(i)==1) thereAreTriggers = 1;
    }
    std::cout<<"were there triggers?  "<<thereAreTriggers<<"\n";
    h_TriggersSeen->Fill(thereAreTriggers);

    //for(int i=0; i<32; i++) std::cout<<"  "<<i<<"  "<<thisTrigger.Triggered(i)<<"\n";
    //std::string beam  ="+USTOF";
    //std::string michel="+MICHEL+COSMICON";
    //std::string ped   ="+PULSER+PEDESTALON-SCINTGATE-BUSYA-BUSYC";
    //std::cout<<"Is pedestal event? "<<fTriggerFilterAlg.doesTriggerPassFilter(thisTrigger,ped)<<"\n";
    //std::cout<<"Is beam event? "<<fTriggerFilterAlg.doesTriggerPassFilter(thisTrigger,beam)<<"\n";
    //std::cout<<"Is michel event? "<<fTriggerFilterAlg.doesTriggerPassFilter(thisTrigger,michel)<<"\n";
  }
  */
 


  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);

  // Skip event if no waveforms found
  if( (size_t)WaveformHandle->size() == 0 ){
    std::cout << "No optical detector data found; skipping event.\n";
    return;
 
  // If waveforms found, look for the specified PMTs
  } else {


    // Grab the first OpDetPulse in the handle just to check the timestamp
    Timestamp   = (float(WaveformHandle->at(0).PMTFrame())*8.)/1.0e09;
    if( Timestamp < fTimestamp_T1 || Timestamp > fTimestamp_T2 ) {
      std::cout<<"Timestamp out of range; skipping event\n";
      return;
    }


    // Keep a count of the total events analyzed
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

      // Get the waveform
      std::vector<short>  Wfm       = ThePulse.Waveform();
      size_t              NSamples  = Wfm.size();
      
      std::cout << "\n"; 
      std::cout << runnr<<"/"<<subrunnr<<"/"<<eventnr<<"  OpDet "<<ch<<", index "<<ich<<": " << NSamples << " samples, timestamp "<<Timestamp<<" sec\n";

      // Find waveform baseline and RMS using two different methods, and compare in a 2D histogram
      std::vector<float> tmp = fOpHitBuilderAlg.GetPedestalAndRMS( Wfm, 0, fBaselineWindowLength );
      std::vector<float> tmp2= fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
      float baseline  = tmp[0];
      float rms       = tmp[1];
      h_WfmRMS[ich]   ->Fill(rms*fMvPerADC);
      h_Baselines     ->Fill(tmp2[0],tmp[0]);
      
      // Determine range of samples to scan over
      short t1 = std::max((short)ThePulse.FirstSample() + fT1, 100);
      short t2 = std::min((short)ThePulse.FirstSample() + fT2, (int)NSamples-100);
      
      std::cout << "Waveform raw baseline: " << baseline << " +/- " << rms <<" ADC ("<<baseline*fMvPerADC<<" +/- "<<rms*fMvPerADC<<" mV)\n";
      std::cout<<"Searing samples "<<t1<<" - "<<t2<<"\n";
      
      // If the waveform baseline is noisy, skip
      if( rms*fMvPerADC >= fWfmAbsRMSCut[ich] ) {
        std::cout<<"Baseline too noisy... moving on.\n";
        continue;
      }

      // Declare and initialize variables
      float   integral = 0;
      bool    flag = false;
      int     windowsize = 0;
      int     counter = -fPreWindow;
      int     dipcounter = 0;
      short   quietTime = 0;
      short   deadTime = 0;
      float   prePE_baseline = 0;
      float   prePE_rms = 999;
      float   hit_grad = 0;
      int     hit_time = -1;
      float   PE_amp  = -1;
  
      // Now begin looking for single PEs
      std::cout << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh[ich] << ", threshPersist " << fThreshPersist <<") \n";

      // Make gradient
      std::vector<float> g = fOpHitBuilderAlg.MakeGradient(Wfm);

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
      
        if(!flag) prePE_baseline = 0;
        float y               = wfm_corrected[i] - prePE_baseline;
        float y_mV            = y*fMvPerADC;
        bool  IsLive          = ((quietTime >= fQuietTimeMin)&&(deadTime >= fDeadTimeMin));
        bool  IsOverThresh    = (y  >= fPulseHitRMSThresh[ich]*rms);
        bool  IsOverThreshNeg = (y  <= -fPulseHitRMSThresh[ich]*rms);
        if(IsOverThreshNeg)   { dipcounter++;}
        else{                   dipcounter=0;}
        bool  NegativeDipDetected = (dipcounter >= 2);
        bool  IsOverLimit     = (y_mV >= fPulseHitThreshHigh[ich]); 
        bool  IsPECandidate   = ((IsOverThresh) && (!IsOverLimit) && ( fabs(g[i]) >= fabs(fGradientCut)) );

        // If we're already in a PE window, increment the
        // counters and add to integral
        if( flag ) {
         
          std::cout
          << "  " << i << "  y_mV = " << y_mV << " mV (window size " <<windowsize<< "), "
          << " thresh " << fPulseHitRMSThresh[ich]*rms*fMvPerADC << " mV, g " << g[i]<<", quietTime = "<<quietTime<<"   deadTime = "<<deadTime<<"\n";
          
          counter++;
          windowsize++;
          integral += y;
         
          // Scan for amplitude within 5ns of hit
          if( counter < 5 ){
            if( y_mV > PE_amp ) PE_amp = y_mV;
            std::cout<<"Pe amp "<<PE_amp<<"\n";
          }

          // Store signal values
          if(tmp_wfm_i < SER_bins){
            tmp_wfm[tmp_wfm_i] = y;
            tmp_wfm_i++;
          }
          
          // If another PE is detected after at least 5 ns, extend window by resetting counter
          if( counter >= 5 && IsPECandidate ){
            std::cout << "  Secondary hit, extending window\n";
            counter = 0;
          }
          
          // If (a) a "dip" below baseline was detected, (b) a pulse extends above upper 
          // limit, or (c) if window length is exceeded due to multiple merges, abort mission.
          if( NegativeDipDetected || IsOverLimit || (windowsize > fMaxWindowFactor*SER_bins)){
            std::cout << "  abort!\n\n";

            // avoid infinite loops by setting 'i' a little bit ahead
            // of the hit that triggered the start of this integration.
            if( i < hit_time + 5 ) i = hit_time + 5;
            
            counter = -fPreWindow;
            hit_grad = 0;
            hit_time = -1;
            PE_amp = -1;
            integral = 0;
            flag = false;
            windowsize = 0;
            tmp_wfm.clear();
            tmp_wfm_i = 0;
            continue;
          }
          
          // If we reached the end of the allotted window, add
          // integral to vector and reset everything
          if( counter == fPostWindow ){
            
            std::cout
            << "   Finished PE window of size "<<windowsize<<", tmp_wfm_i = "<<tmp_wfm_i << "  "
            << integral << " ADCs\n\n";
            
            h_SER[ich]        ->Fill(integral);
            h_SERAmp[ich]     ->Fill(PE_amp);
            
            // If this is > 20 ADC, plot the time:
            if(integral > 20 ) h_PhelTime[ich]   ->Fill(hit_time);
            
            // Add to average waveform if it looks good
            if( (windowsize == SER_bins) && fabs(integral - fSinglePE[ich])<=fSinglePE_tolerance[ich] ){
              for(short ii=0; ii<SER_bins; ii++) SERWaveform[ich].at(ii) += tmp_wfm[ii]*fMvPerADC;
              h_AvePECount[ich]->Fill(0);
            }
    
            integral = 0;
            hit_grad = 0;
            hit_time = -1;
            PE_amp    = -1;
            windowsize = 0;
            counter = -fPreWindow;
            flag = false;
            tmp_wfm.clear();
            tmp_wfm_i=0;
            continue;
          }
        
        } // endif flag
   
        
       
        // If we're not yet integrating a PE window, signal is within bounds,
        // and enough quietTime has elapsed, then we're in business
        if( !flag && IsPECandidate && IsLive ){
          
          // Find pre-PE baseline
          prePE_baseline = 99;
          prePE_rms      = 99;
          std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS(wfm_corrected,i-fPrePEBaselineWindow-fPreWindow,i-fPreWindow);
          prePE_baseline = tmp[0];
          prePE_rms      = tmp[1];
          h_PrePhelRMS[ich]->Fill(prePE_rms*fMvPerADC);
          
          // If "threshPersist" functionality being used, look ahead and make 
          // sure signal passes the threshold for requisite consecutive samples
          bool flagger = true;
          if( fThreshPersist > 0 ){
            for( short j=0; j<fThreshPersist; j++ ){
              if ( wfm_corrected[i+j] - prePE_baseline < fPulseHitRMSThresh[ich]*rms ) flagger = false;
            }
          }

          // When taking the pre-PE baseline into account, require that the 
          // candidate still passes the threshold:
          y             = wfm_corrected[i] - prePE_baseline;
          y_mV          = y*fMvPerADC;
          IsOverThresh  = (y    >= fPulseHitRMSThresh[ich]*rms);
          IsOverLimit   = (y_mV >= fPulseHitThreshHigh[ich]); 
          IsPECandidate = ((IsOverThresh) && (!IsOverLimit) && ( fabs(g[i]) >= fabs(fGradientCut)) ); 

          // Require flat pre-PE region and that threshold persist requirement was met
          if( (IsPECandidate) && (prePE_rms <= fPrePE_RMSFactorCut[ich]*rms) && (flagger) ){
            
            // Found a "PE hit"
            flag      = true;
            hit_grad  = g[i];
            hit_time  = i;
           
            std::cout 
            << "** OpDet "<<ch<<": Candidate at "<<hit_time<<", "<<y_mV<<" mV, g = "<<hit_grad
            << ", prePE_baseline  = "<<prePE_baseline*fMvPerADC<<" +/- "<<prePE_rms*fMvPerADC
            <<"\n";

            // Go back "preWindow" number of samples
            i -= fPreWindow;
          
          } 
    
        } // <-- end if(PE cand)
        
        // "quietTime" will increment for every sample where the signal 
        // is below the limit (default 5mV) and where a sustained negative 
        // dip has not been seen.  If either of these conditions is met, 
        // the counter is reset.
        if( !IsOverLimit && !NegativeDipDetected ){ quietTime++;    }
        else                                      { quietTime = 0;  }

        // "deadTime" is similar, but with stricter conditions.
        if( !IsOverThresh && !IsOverThreshNeg ) { deadTime++; }
        else                                    { deadTime = 0; }


      } // <-- end scan over waveform 


    } // endLoop over OpDets
  
  } // endIf Waveformhandle is valid

}

//#######################################################################
void OpDetSER::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
 
  h_TotalEvents     = tfs->make<TH1I>("TotalEvents","Total events",1,0,1);
  //h_TriggersSeen    = tfs->make<TH1I>("TriggersSeen","0 = NO trigger bits found for event",2,0,2);
  //h_NTrigs          = tfs->make<TH1I>("NTrigs",  "Trigger bits fired during run",  32,0,32);
  //h_NTrigs_Pulser   = tfs->make<TH1I>("NTrigs_Pulser",  "Pulser triggers per subrun",  15,0,15);
  //h_NTrigs_Beam     = tfs->make<TH1I>("NTrigs_Beam",      "Beam triggers per subrun",    100,0,100);
  //h_NTrigs_Michel   = tfs->make<TH1I>("NTrigs_Michel",  "Michel triggers per subrun",  100,0,100);
  h_Baselines       = tfs->make<TH2F>("Baselines","Baseline comparisons;Baseline [ADC];Pedestal [ADC]",200,915,935,200,915,935);
  h_Baselines       ->SetOption("colz");
  
  // Create histograms
  for(size_t i=0; i<fSelectChannels.size(); i++){

    size_t ch   = fSelectChannels[i];
    size_t ich  = iCh[ch]; 
    
    sprintf(histName,"%lu_SER",ch); 
    sprintf(histTitle,"SER for OpDet %lu;Integrated ADC;Counts",ch);
    h_SER[ich]       = tfs->make<TH1F>(histName,histTitle,400,-50.,350.);
    
    sprintf(histName,"%lu_SERAmp",ch); 
    sprintf(histTitle,"PE candidate amplitudes for OpDet %lu;mV;Counts",ch);
    h_SERAmp[ich]       = tfs->make<TH1F>(histName,histTitle,150,0.,3.);
   
    sprintf(histName,"%lu_AvePEWfm",ch); 
    sprintf(histTitle,"Average PE waveform for OpDet %lu;ns;mV",ch);
    h_AvePEWfm[ich]  = tfs->make<TH1F>(histName,histTitle,SER_bins,0.,(float)SER_bins);
    
    sprintf(histName,"%lu_AvePECount",ch); 
    sprintf(histTitle,"Total numbwer of averaged single-PE waveforms for OpDet %lu",ch);
    h_AvePECount[ich]  = tfs->make<TH1I>(histName,histTitle,1,0,1);
    
    sprintf(histName,"%lu_RMS",ch); 
    sprintf(histTitle,"RMS for OpDet %lu;Baseline RMS [mV];Counts",ch);
    h_WfmRMS[ich]    = tfs->make<TH1F>(histName,histTitle,100,0.0,0.5);
    
    sprintf(histName,"%lu_PrePhelRMS",ch);
    sprintf(histTitle,"Pre-PE candidate RMS for OpDet %lu;Pre-PE RMS [mV];Counts",ch);
    h_PrePhelRMS[ich]    = tfs->make<TH1F>(histName,histTitle,100,0.0,0.5);
    
    sprintf(histName,"%lu_PhelTime",ch); 
    sprintf(histTitle,"PE times for OpDet %lu;Sample [ns];Counts/500ns",ch);
    h_PhelTime[ich]    = tfs->make<TH1I>(histName,histTitle,(int)(fNSamples/500.),0,fNSamples);
   
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
  fGradientCut            = p.get< float >                ("GradientCut",0);
  fPulseHitThreshHigh     = p.get< std::vector<float> >   ("PulseHitThresh_high",PulseHitThreshHighDefault);
  fPulseHitRMSThresh      = p.get< std::vector<float> >   ("PulseHitRMSThresh",PulseHitRMSThreshDefault);
  fSinglePE               = p.get< std::vector<float> >   ("SinglePE",SinglePEDefault);
  fSinglePE_tolerance     = p.get< std::vector<float> >   ("SinglePE_tolerance",SinglePEToleranceDefault);
  fWfmAbsRMSCut           = p.get< std::vector<float> >   ("WfmAbsRMSCut", WfmAbsRMSCutDefault);
  fPedestalMean_lowerLim  = p.get< std::vector<float> >   ("PedestalMean_lowerLim",PedestalMeanLowerLimDefault);
  fPedestalMean_upperLim  = p.get< std::vector<float> >   ("PedestalMean_upperLim",PedestalMeanUpperLimDefault);

  fTriggerUtility         = p.get< std::string >  ("TriggerUtility","FragmentToDigit");
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fAttemptFit             = p.get< bool >         ("AttemptFit","true");
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fTimestamp_T1           = p.get< float >        ("Timestamp_T1",0);
  fTimestamp_T2           = p.get< float >        ("Timestamp_T2",60);
  fPreWindow              = p.get< short >       ("PreWindow",5);
  fPostWindow             = p.get< short >       ("PostWindow",45);
  fMaxWindowFactor        = p.get< float >        ("MaxWindowFactor",2);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.2);
  fPrePEBaselineWindow    = p.get< short >        ("PrePEBaselineWindow",100);
  fThreshPersist          = p.get< short >        ("ThreshPersist",3);
  fDeadTimeMin            = p.get< short >        ("DeadTimeMin",10);
  fQuietTimeMin           = p.get< short >        ("QuietTimeMin",200);
  fVerbose                = p.get< bool >         ("Verbose",true);
  fNSamples               = p.get< short >        ("NSamples",28672);
}

void OpDetSER::beginRun(art::Run const & r){
  fTriggerFilterAlg.loadXMLDatabaseTable( r.run() );
}

void OpDetSER::beginSubRun(art::SubRun const & sr){}

void OpDetSER::endJob()
{

  // Fit line to baseline comparison plot
  TF1 f1("f_line1","[0] + [1]*x",915,935);
  f1.SetParName(0,"Offset");
  f1.SetParName(1,"Slope");
  f1.SetParameter(0,0);
  f1.SetParameter(1,1);
  h_Baselines->Fit("f_line1","R");

  // Print out general parameters used
  std::cout 
    << "============================================================\n"
    << "Ending SER program.\n"
    << "  Pre-PE baseline length  " << fPrePEBaselineWindow << "\n"
    << "  Threshhold persist      " << fThreshPersist <<"\n" 
    << "  Graident cut            " << fGradientCut << "\n"
    << "  Pre / Post Window       " << fPreWindow << "," << fPostWindow << "\n\n";

  
  for(size_t i=0; i<fSelectChannels.size(); i++){

    // Print out photodetector-specific parameters
    size_t ich = iCh[fSelectChannels[i]];
    std::cout
    << "-------------------------------------------\n"
    << "  channel                 " << fOpDetChannels[ich] << "\n"
    << "  Abs Wfm RMS cut         " << fWfmAbsRMSCut[ich] << " mV\n"
    << "  Pre-PE RMS cut          " << fPrePE_RMSFactorCut[ich] << " x wfm RMS\n"
    << "  RMS thresh factor       " << fPulseHitRMSThresh[ich] <<" x wfm RMS\n"
    << "  Mean baseline RMS	  " << h_WfmRMS[ich]->GetMean(1) << "\n"
    << "  --> PE candidates found " << h_SER[ich]->GetEntries() << "\n\n";
    
    // Normalize the summed PE waveform to get an average
    if(h_AvePECount[ich]->GetEntries() > 0 ){
      std::cout
      << "  Ave PE waveform (sample [ns], ADC):\n";
      
      float integral = 0;
      for( int i = 0; i < SER_bins; i++) {
        float w = SERWaveform[ich].at(i) / float(h_AvePECount[ich]->GetEntries()); 
        h_AvePEWfm[ich]->Fill(i,w);
        integral += w/fMvPerADC;
        std::cout<<"    "<<i<<"     "<<w/fMvPerADC<<"\n";
      }
    
      std::cout
      << "  Integral: "<<integral<< " ADC ("<<h_AvePECount[ich]->GetEntries()<<" waveforms averaged)\n\n";
    
    }
 
 
    if( fAttemptFit ){
      std::cout
      <<"Attempting SER fit...\n";
       
      // Fit out the SER
      TF1 f_gaus("f_gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))",-1000,1000);
      TF1 SER_fit("SER_fit","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/(sqrt(1)*[5]),2)) + [6]*exp(-0.5*pow((x-2.*[4])/(sqrt(2)*[5]),2)) + [7]*exp(-0.5*pow((x-3.*[4])/(sqrt(3)*[5]),2))",-50,400);
 
      SER_fit.SetNpx(1000);
  
      float max = (float)h_SER[ich]->GetMaximum();
  
      SER_fit.SetParName(0,"Noise norm");
      SER_fit.SetParName(1,"Noise mean");
      SER_fit.SetParName(2,"Noise sigma");
      SER_fit.SetParName(3,"1PE norm");
      SER_fit.SetParName(4,"1PE mean (SER)");
      SER_fit.SetParName(5,"1PE sigma");
      SER_fit.SetParName(6,"2PE norm");
      SER_fit.SetParName(7,"3PE norm");
  
      // "Noise" component (gaus)
      SER_fit.SetParameter(0,max);
      SER_fit.SetParLimits(0,0.,max);
      SER_fit.SetParameter(1,0);
      SER_fit.SetParLimits(1,fPedestalMean_lowerLim[ich],fPedestalMean_upperLim[ich]);
      SER_fit.SetParameter(2,10);
      SER_fit.SetParLimits(2,3.,30.);

      // 1PE (gaus)
      SER_fit.SetParameter(3,max);
      SER_fit.SetParLimits(3,0.,max);
      SER_fit.SetParameter(4,fMean_set[ich]);
      SER_fit.SetParLimits(4,fMean_lowerLim[ich],fMean_upperLim[ich]);
      SER_fit.SetParameter(5,20);
      SER_fit.SetParLimits(5,10.,60.);

      // 2PE (gaus)
      SER_fit.SetParameter(6,0.1*max);
      SER_fit.SetParLimits(6,0.,0.5*max);

      // 3PE (gaus)
      SER_fit.SetParameter(7,1.);
      SER_fit.SetParLimits(7,0.,0.5*max);

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

    }
  
    std::cout<<"=====================================================\n";
 
  }
  
}

void OpDetSER::endRun(art::Run const & r){}

void OpDetSER::endSubRun(art::SubRun const & sr){}

DEFINE_ART_MODULE(OpDetSER)
