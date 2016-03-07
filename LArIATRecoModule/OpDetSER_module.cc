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

//LAriatSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"
//#include "LArIATRecoAlg/TriggerFilterAlg.h"
//#include "Utilities/DatabaseUtilityT1034.h"

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

  // Custom functions
  //std::vector<float> GetBaselineAndRMS( std::vector<short>, short, short);
  //std::vector<float> GetBaselineAndRMS( std::vector<float>, short, short);
  //std::vector<float> MakeGradient( std::vector<short>);

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
  float               fPrePE_RMSCut;

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

};


//#######################################################################
OpDetSER::OpDetSER(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
fOpHitBuilderAlg(p)
{
  this->reconfigure(p);
  SER_bins      = fPreWindow + fPostWindow;
  SERWaveform   .resize(SER_bins);
}

//#######################################################################
void OpDetSER::analyze(art::Event const & e)
{
  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);

  // Skip event if no waveforms found
  if( (size_t)WaveformHandle->size() == 0 ){
    LOG_VERBATIM("OpDetSER") << "No optical detector data found; skipping event.";
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
        
        LOG_VERBATIM("OpDetSER")
        << "PMT pulse recorded (" << NSamples << " samples, trigger at "<<TrigSample<<")\n"
        << "Timestamp " << Timestamp << " sec";

      }
    }
  }
  // If we somehow after all this we still don't have
  // the PMT we want, skip this weird event.
  if( !gotPMT ) return;
   
   
  // Now begin looking for single PEs
  LOG_VERBATIM("OpDetSER")
  << "Beginning search for single PE candidates (RMS thresh x " << fPulseHitRMSThresh << ")";

  short t1 = std::max(TrigSample + fT1,0);
  short t2 = std::min(TrigSample + fT2,(int)NSamples);

  // Find waveform baseline and RMS 
  std::vector<float> tmp = fOpHitBuilderAlg.GetBaselineAndRMS( Wfm, 0, fBaselineWindowLength );
  float baseline  = tmp[0];
  float rms       = tmp[1]; 

  float   integral = 0;
  bool    flag = false;
  int     windowsize = 0;
  int     counter = 0;
  int     flat_samples_count = 0;
  float   prePE_baseline = -99;
  float   prePE_rms = 99;

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
     
      LOG_VERBATIM("OpDetSER")
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag;
      
      counter++;
      windowsize++;
      float yc = wfm_corrected[i] - prePE_baseline;
      integral += yc;

      if(tmp_wfm_i < SER_bins){
        tmp_wfm[tmp_wfm_i] = yc;
        tmp_wfm_i++;
      }
      
      // If another PE is detected after at least 5 ns, extend the window by resetting counter
      if( counter >=5 && IsPECandidate ){
        LOG_VERBATIM("OpDetSER") << "  Secondary hit, extending window";
        counter = 0;
      }
      
      // If pulse extends above upper limit or if window length
      // is exceeded due to multiple merges, abort mission.
      if( IsOverLimit || (windowsize > 2*fPostWindow) ){
        LOG_VERBATIM("OpDetSER") << "  abort!";
        counter = 0;
        hit_grad = 0;
        integral = 0;
        flag = false;
        windowsize = 0;
        flat_samples_count = 0;
        tmp_wfm.clear();
        tmp_wfm_i = 0;
        continue;
      }
      
      // If we reached the end of the allotted window, add
      // integral to vector and reset everything
      if( counter == fPostWindow ){
       
        LOG_VERBATIM("OpDetSER") 
        << "Finished PE window of size "<<windowsize<<", "
        << integral << " ADCs, g = " << hit_grad;
        
        h_SER   ->Fill(integral);
        h_SER_g ->Fill(hit_grad);
        
        // Add to average waveform if it looks good
        if( (windowsize+fPreWindow == SER_bins) && fabs(integral - fSinglePE)<=5 ){
          
          LOG_VERBATIM("OpDetSER") << "Add to average PE wfm.";
          
          for(short ii=0; ii<SER_bins; ii++){
            SERWaveform.at(ii) += tmp_wfm[ii]*fMvPerADC;
          }
          SERWaveform_count++;
        }

        integral = 0;
        hit_grad = 0;
        windowsize = 0;
        counter = 0;
        flat_samples_count = 0;
        flag = false;
        tmp_wfm.clear();
        tmp_wfm_i=0;
        continue;
      }
    
    } // <-- end if(flag)

    int buffer_min = 100;

    if( !IsOverThresh && flat_samples_count <= 50 ) flat_samples_count++;

    // If we're not yet integrating a PE window, signal is within bounds,
    // and the previous 10 samples were below threshold, then we're in business!
    if( !flag && IsPECandidate ){
      
      // Find pre-PE baseline
      prePE_baseline = 99;
      prePE_rms      = 99;
      std::vector<float> tmp = fOpHitBuilderAlg.GetPedestalAndRMS(wfm_corrected,i-buffer_min,i);
      prePE_baseline = tmp[0];
      prePE_rms      = tmp[1];
     
      LOG_VERBATIM("OpDetSER")
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag << "\n"
      << "  Potential PE!  preBS/RMS " << prePE_baseline << ", "<< prePE_rms;
 
      // Require flat pre-PE region
      if( prePE_rms < fPrePE_RMSCut){ 
        // Found a "PE hit", so start integral by
        // adding preceeding prepulse region
        flag = true;
        hit_grad = g[i];

        // Add up previous "prewindow" samples
        for(short ii=0; ii<=fPreWindow; ii++){
          integral += wfm_corrected[i-fPreWindow+ii] - prePE_baseline;
          tmp_wfm[tmp_wfm_i] = wfm_corrected[i-fPreWindow+ii] - prePE_baseline;
          tmp_wfm_i++;
        }
   
        LOG_VERBATIM("OpDetSER") << "  Looks good!  Beginning integration..."; 
      } else {
        LOG_VERBATIM("OpDetSER") << "  Doesn't pass preBS cut, moving on..."; 
      }

    } // <-- end if(PE cand)
    
    if( IsOverThresh ) flat_samples_count = 0;

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
}

void OpDetSER::reconfigure(fhicl::ParameterSet const & p)
{
  fDAQModule              = p.get< std::string >  ("DAQModule","daq");
  fInstanceName           = p.get< std::string >  ("InstanceName","");
  fOpDetChannel           = p.get< size_t >       ("OpDetChannel",1);
  fBaselineWindowLength   = p.get< short >       ("BaselineWindowLength",1000);
  fMean_set               = p.get< float >        ("Mean_set",50);
  fMean_lowerLim          = p.get< float >        ("Mean_lowerLim",30);
  fMean_upperLim          = p.get< float >        ("Mean_upperLim",80);
  fPrePE_RMSCut           = p.get< float >        ("PrePE_RMSCut",2.5);
  fSinglePE               = p.get< float >        ("SinglePE",85);
  fAttemptFit             = p.get< bool >         ("AttemptFit","true");
  fT1                     = p.get< short >       ("fT1",1000);
  fT2                     = p.get< short >       ("fT2",19000);
  fPulseHitThreshHigh     = p.get< float >        ("PulseHitThresh_high",10);
  fPulseHitRMSThresh      = p.get< float >        ("PulseHitRMSThresh",3.);
  fGradientCut            = p.get< float >        ("GradientCut",-4);
  fPreWindow              = p.get< short >       ("PreWindow",5);
  fPostWindow             = p.get< short >       ("PostWindow",45);
  fMvPerADC               = p.get< float >        ("MvPerADC",0.2);

}

void OpDetSER::beginRun(art::Run const & r){}

void OpDetSER::beginSubRun(art::SubRun const & sr){}

void OpDetSER::endJob(){}

void OpDetSER::endRun(art::Run const & r){}

void OpDetSER::endSubRun(art::SubRun const & sr){}

/*
//--------------------------------------------------------------
// Get (baseline,rms) of a segment of an std::vector<short>
std::vector<float> GetBaselineAndRMS( std::vector<short> wfm, short x1, short x2 )
{
  // Convert to vector<float> and send into
  // other version of function
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return GetBaselineAndRMS(wfm_float,x1,x2);
}

// Get (baseline,rms) of a segment of an std::vector<float>
std::vector<float> GetBaselineAndRMS( std::vector<float> wfm, short x1, short x2 )
{
  float mean = 0;
  float sumSquares = 0;
  short N = x2 - x1;
  for ( short i = x1; i < x2; i++ ) mean += wfm[i]/N;
  for ( short i = x1; i < x2; i++ ) sumSquares += pow(wfm[i]-mean,2);
  std::vector<float> out(2);
  out[0] = mean;
  out[1] = sqrt(sumSquares/N);
  return out;
}

//-------------------------------------------------------------- 
// MakeGradient 
std::vector<float> MakeGradient( std::vector<short> wfm ) 
{
  std::vector<float> g(wfm.size());
  g[0]=0;
  g[1]=0;
  for(size_t i=2; i<wfm.size(); i++) g[i] = float(wfm[i] - wfm[i-2])*0.5;
  return g;
}
*/

DEFINE_ART_MODULE(OpDetSER)
