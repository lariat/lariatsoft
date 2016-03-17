////////////////////////////////////////////////////////////////////////
// Class:       OpDetSER
// Module Type: analyzer
// File:        OpDetSER_module.cc
//
// Generated at Fri Mar  4 04:40:04 2016 by William Foreman using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

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

// ROOT includes
#include <TF1.h>
#include <TH1F.h>

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
  std::string         fDAQModule;
  std::string         fInstanceName;
  size_t              fOpDetChannel;
  float               fMvPerADC;
  bool                fAttemptFit;
  short               fBaselineWindowLength;
  float               fMean_set;
  float               fMean_lowerLim;
  float               fMean_upperLim;
  short               fT1;
  short               fT2;
  float               fPulseHitThreshHigh;
  float               fPulseHitRMSThresh;
  float               fGradientCut;
  short               fPostWindow;
  short               fPreWindow;
  float               fSinglePE;
  float               fSinglePE_tolerance;
  float               fPrePE_RMSCut;
  float               fUsePrePEBaseline;
  short               fThreshPersist;
  float               fWfmAbsRMSCut;
  short               fDeadTimeMin;

  short               SER_bins;
  bool                gotPMT;
  size_t              NSamples;
  raw::OpDetPulse     OpDetPulse;
  std::vector<short>  Wfm;
  float               Timestamp;
  short               TrigSample;
  std::vector<float>  SERWaveform;
  int                 SERWaveform_count;

  // Alg objects
  OpHitBuilderAlg     fOpHitBuilderAlg;
  
  // Histograms
  TH1F* h_SER;
  TH1F* h_SER_g;
  TH1F* h_AvePEWfm;
};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p)
{
  this                ->reconfigure(p);
  SER_bins            = fPreWindow + fPostWindow;
  SERWaveform         .resize(SER_bins);
  SERWaveform_count   = 0;
}

//#######################################################################
void OpDetSER::analyze(art::Event const & e)
{
  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);

  // Skip event if no waveforms found
  if( (size_t)WaveformHandle->size() == 0 ){
    std::cout << "No optical detector data found; skipping event.\n";
    return;
 
  // If waveforms found, look for the specified PMT and save it
  } else {
    gotPMT = false;
    for( size_t ipulse = 0; ipulse < (size_t)WaveformHandle->size(); ipulse++){
      art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
      raw::OpDetPulse ThePulse = *ThePulsePtr;
      if( ThePulse.OpChannel() == fOpDetChannel ) {
        gotPMT = true;
        OpDetPulse  = ThePulse; 
        Wfm         = OpDetPulse.Waveform();
        NSamples    = (size_t)Wfm.size();
        Timestamp   = (float(OpDetPulse.PMTFrame())*8.)/1.0e09;
        TrigSample  = (short)ThePulse.FirstSample(); 
        
        std::cout << "PMT pulse recorded (" << NSamples << " samples, trigger at "<<TrigSample<<")\n"
        << "Timestamp " << Timestamp << " sec \n";

      }
    }
  }
  // If we somehow after all this we still don't have
  // the PMT we want, skip this weird event.
  if( !gotPMT ) return;
 
   

  short t1 = std::max(TrigSample + fT1,0);
  short t2 = std::min(TrigSample + fT2,(int)NSamples);

  // Find waveform baseline and RMS 
  std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
  float baseline  = tmp[0];
  float rms       = tmp[1]; 

  std::cout << "Waveform raw baseline: " << baseline << " +/- " << rms <<" ADC\n";
  if( rms*fMvPerADC >= fWfmAbsRMSCut ) return;

  float   integral = 0;
  bool    flag = false;
  int     windowsize = 0;
  int     counter = 0;
  //int     flat_samples_count = 0;
  short   quietTime = 0;
  float   prePE_baseline = -99;
  float   prePE_rms = 99;
  
  // Now begin looking for single PEs
  std::cout << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh << ", threshPersist " << fThreshPersist <<") \n";

  // Make gradient
  std::vector<float> g = fOpHitBuilderAlg.MakeGradient(Wfm);
  float hit_grad = 0;

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
    bool  IsOverThresh  = (y  >= fPulseHitRMSThresh*rms);
    bool  IsOverLimit   = (yy >= fPulseHitThreshHigh); 
    bool  IsPECandidate = ((IsOverThresh) && (!IsOverLimit) && (g[i] <= fGradientCut) );

    // If we're already in a PE window, increment the
    // counters and add to integral
    if( flag ) {
     
      std::cout
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << " mV, g " << g[i]<<", deadTime = "<<quietTime<<"\n";
      
      counter++;
      windowsize++;
      float yc = wfm_corrected[i] - fUsePrePEBaseline*prePE_baseline;
      integral += yc;

      if(tmp_wfm_i < SER_bins){
        tmp_wfm[tmp_wfm_i] = yc;
        tmp_wfm_i++;
      }
      
      // If another PE is detected after at least 5 ns, extend the window by resetting counter
      if( counter >=10 && IsPECandidate ){
        std::cout << "  Secondary hit, extending window\n";
        counter = 0;
      }
      
      // If pulse extends above upper limit or if window length
      // is exceeded due to multiple merges, abort mission.
      if( IsOverLimit || (windowsize > fPreWindow + 2*fPostWindow) ){
        std::cout << "  abort!\n";
        counter = 0;
        hit_grad = 0;
        integral = 0;
        flag = false;
        windowsize = 0;
        //flat_samples_count = 0;
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
        << integral << " ADCs, g = " << hit_grad << "\n";
        
        h_SER   ->Fill(integral);
        h_SER_g ->Fill(hit_grad);
        
        // Add to average waveform if it looks good
        if( (windowsize == SER_bins) && fabs(integral - fSinglePE)<=fSinglePE_tolerance ){
          
          std::cout << "Add to average PE wfm.\n";
          
          for(short ii=0; ii<SER_bins; ii++){
            SERWaveform.at(ii) += tmp_wfm[ii]*fMvPerADC;
          }
          SERWaveform_count++;
        }

        integral = 0;
        hit_grad = 0;
        windowsize = 0;
        counter = 0;
        //flat_samples_count = 0;
        quietTime = 0;
        flag = false;
        tmp_wfm.clear();
        tmp_wfm_i=0;
        continue;
      }
    
    } // <-- end if(flag)

    int buffer_min = 100;

    if( !IsOverThresh ) quietTime++;

    // If we're not yet integrating a PE window, signal is within bounds,
    // and enough dead-time (quietTime) has elapsed, then we're in business
    if( !flag && IsPECandidate && quietTime >= fDeadTimeMin ){
      
      // Find pre-PE baseline
      prePE_baseline = 99;
      prePE_rms      = 99;
      //std::vector<float> tmp = fOpHitBuilderAlg.GetPedestalAndRMS(wfm_corrected,i-buffer_min,i);
      std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS(wfm_corrected,i-buffer_min,i);
      prePE_baseline = tmp[0];
      prePE_rms      = tmp[1];
     
      std::cout
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag << "\n"
      << "  Potential PE!  preBS = "<< prePE_baseline*fMvPerADC << " mV, preRMS = "<< prePE_rms*fMvPerADC <<" mV, quietTime count "<<quietTime<<"\n";
 
      // Look a few samples ahead and make sure the signal
      // remains above threshold for some time...
      bool flagger = true;
      if( fThreshPersist > 0 ){
        std::cout<< "  Peeking at next few samples....\n";
        for( short j=0; j<fThreshPersist; j++ ){
          std::cout<<"    "<<j<<"  "<<wfm_corrected[i+j]*fMvPerADC << " mV \n";
          if ( wfm_corrected[i+j] < fPulseHitRMSThresh*rms ) flagger = false;
        }
      }
      // Require flat pre-PE region and that signal stays over thresh
      // for at least 5 samples
      if( (prePE_rms <= fPrePE_RMSCut*rms) && (flagger)){
        
        // Found a "PE hit", so start integral by
        // adding preceeding prepulse region
        flag = true;
        hit_grad = g[i];

        // Add up previous "prewindow" samples
        for(short ii=0; ii <= fPreWindow; ii++){
          integral += wfm_corrected[i-fPreWindow+ii] - fUsePrePEBaseline*prePE_baseline;
          tmp_wfm[tmp_wfm_i] = wfm_corrected[i-fPreWindow+ii] - prePE_baseline;
          tmp_wfm_i++;
          windowsize++;
          counter = 1;
        }
      
        std::cout << "  Looks good!  Beginning integration...\n"; 
      } else {
        std::cout << "  Doesn't pass quality cuts, moving on...\n"; 
      }

    } // <-- end if(PE cand)
    
    if( IsOverThresh ) quietTime = 0;

  } // <-- end scan over waveform 
}

//#######################################################################
void OpDetSER::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
  
  // Create histograms
  h_SER       = tfs->make<TH1F>("SER","SER;Integrated ADC;Counts",500,-50.,450.);
  h_SER_g     = tfs->make<TH1F>("SER_g","SER_g;Gradient at SER candidate hit;Counts",64,-8,8.);
  h_AvePEWfm  = tfs->make<TH1F>("AvePEWfm","Average PE waveform",SER_bins,0.,(float)SER_bins);
  h_AvePEWfm  ->GetXaxis()->SetTitle("ns");
  h_AvePEWfm  ->GetYaxis()->SetTitle("mV");
}

void OpDetSER::reconfigure(fhicl::ParameterSet const & p)
{
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fOpDetChannel           = p.get< size_t >       ("OpDetChannel",1);
  fBaselineWindowLength   = p.get< short >        ("BaselineWindowLength",1000);
  fMean_set               = p.get< float >        ("Mean_set",50);
  fMean_lowerLim          = p.get< float >        ("Mean_lowerLim",30);
  fMean_upperLim          = p.get< float >        ("Mean_upperLim",80);
  fPrePE_RMSCut           = p.get< float >        ("PrePE_RMSCut",1.5);
  fGradientCut            = p.get< float >        ("GradientCut",-4);
  fPulseHitThreshHigh     = p.get< float >        ("PulseHitThresh_high",5);
  fPulseHitRMSThresh      = p.get< float >        ("PulseHitRMSThresh",3.);
  fSinglePE               = p.get< float >        ("SinglePE",85);
  fSinglePE_tolerance     = p.get< float >        ("SinglePE_tolerance",5);
  fAttemptFit             = p.get< bool >         ("AttemptFit","true");
  fT1                     = p.get< short >       ("T1",3000);
  fT2                     = p.get< short >       ("T2",19000);
  fPreWindow              = p.get< short >       ("PreWindow",5);
  fPostWindow             = p.get< short >       ("PostWindow",45);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.2);
  fUsePrePEBaseline       = p.get< float >        ("UsePrePEBaseline",1);
  fThreshPersist          = p.get< short >        ("ThreshPersist",3);
  fWfmAbsRMSCut           = p.get< float >        ("WfmAbsRMSCut",0.5);
  fDeadTimeMin            = p.get< short >        ("DeadTimeMin",100);
}

void OpDetSER::beginRun(art::Run const & r){}

void OpDetSER::beginSubRun(art::SubRun const & sr){}

void OpDetSER::endJob()
{

  std::cout 
    << "============================================================\n"
    << "Ending SER program.\n"
    << "Parameters: \n"
    << "  Optical channel         " << fOpDetChannel << "\n"
    << "  Abs Wfm RMS cut         " << fWfmAbsRMSCut << " mV\n"
    << "  Pre-PE RMS cut          " << fPrePE_RMSCut << " x wfm RMS\n"
    << "  Use pre-PE BS in int?   " << fUsePrePEBaseline << "\n"
    << "  RMS thresh factor       " << fPulseHitRMSThresh <<" x wfm RMS\n"
    << "  Threshhold persist      " << fThreshPersist <<"\n" 
    << "  Graident cut            " << fGradientCut << "\n"
    << "  Pre / Post Window       " << fPreWindow << "," << fPostWindow << "\n\n";
 
  std::cout << "Found "<< h_SER->GetEntries() << " single PE candidates.";
  
  // SER average waveform
  float integral = 0;
  
  std::cout
    << "------------------------------------------------------\n"
    << "PE waveform (sample [ns], ADC):\n";
    
    if( SERWaveform_count > 0 ){
      
      for( int i = 0; i < SER_bins; i++) {
        float w = SERWaveform.at(i) / float(SERWaveform_count); 
        h_AvePEWfm->Fill(i,w);
        integral += w/fMvPerADC;
        std::cout<<i<<"     "<<w<<"\n";
      }
    }
    std::cout 
      << "-------------------------------------------------------\n"
      << "Ave PE wfm integral: "<<integral<< " ADC ("<<SERWaveform_count<<" waveforms averaged)\n";
 
    if( fAttemptFit ){
       
      // Fit out the SER
      TF1 f_gaus("f_gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))",-1000,1000);
      TF1 SER_fit("SER_fit","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/[5],2)) + [6]*exp(-0.5*pow((x-2.*[4])/[7],2)) + [8]*exp(-0.5*pow((x-3.*[4])/[9],2))",-50,400);
  
    float max = (float)h_SER->GetMaximum();
    std::cout<<"histogram max: "<<max<<"\n";
  
    SER_fit.SetParName(0,"Noise norm");
    SER_fit.SetParName(1,"Noise mean ADC");
    SER_fit.SetParName(2,"Noise sigma");
    SER_fit.SetParName(3,"1PE norm");
    SER_fit.SetParName(4,"1PE mean ADC");
    SER_fit.SetParName(5,"1PE sigma");
    SER_fit.SetParName(6,"2PE norm");
    SER_fit.SetParName(7,"2PE sigma");
    SER_fit.SetParName(8,"3PE norm");
    SER_fit.SetParName(9,"3PE sigma");
  
    // "Noise" component (gaus)
    SER_fit.SetParameter(0,max);
    SER_fit.SetParLimits(0,0.,1.2*max);
    SER_fit.SetParameter(1,0);
    SER_fit.SetParLimits(1,-10,10);
    SER_fit.SetParameter(2,20);
    SER_fit.SetParLimits(2,5,30.);

    // 1PE (gaus)
    SER_fit.SetParameter(3,max);
    SER_fit.SetParLimits(3,0.,1.2*max);
    SER_fit.SetParameter(4,fMean_set);
    SER_fit.SetParLimits(4,fMean_lowerLim,fMean_upperLim);
    SER_fit.SetParameter(5,35);
    SER_fit.SetParLimits(5,20.,60.);

    // 2PE (gaus)
    SER_fit.SetParameter(6,0.1*max);
    SER_fit.SetParLimits(6,0.,0.3*max);
    SER_fit.SetParameter(7,50);
    SER_fit.SetParLimits(7,20.,100.);

    // 3PE (gaus)
    SER_fit.SetParameter(8,1.);
    SER_fit.SetParLimits(8,0.,0.3*max);
    SER_fit.SetParameter(9,80);
    SER_fit.SetParLimits(9,40.,150.);

    h_SER->Fit("SER_fit","R");
    std::cout<<"SER fit results:\n";
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

}

void OpDetSER::endRun(art::Run const & r){}

void OpDetSER::endSubRun(art::SubRun const & sr){}

DEFINE_ART_MODULE(OpDetSER)
