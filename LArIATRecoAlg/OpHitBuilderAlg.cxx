////////////////////////////////////////////////////////////////////
//                                                                  
// These are functions used to process optical detector information
// (particulary PMTs) and find/integrate hits.  They've been used
// mostly to analyze events in the stopping/decaying cosmic muon 
// samples, and to aid in single photoelectron calibration, but will 
// can adapted for more general use.
// 
// Authors: William Foreman, wforeman@uchicago.edu
//
/////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/TriggerData.h"

// LArIATSoft includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"

// C++ includes
#include <iostream>
#include <cmath>
#include <cstdlib>

// ROOT includes
#include <TF1.h>
#include <TGraph.h>



//####################################################################
OpHitBuilderAlg::OpHitBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;
 
  // Set size of vector to hold summed waveforms 
  AveWfmBins = fPrePulseDisplay + fFullWindowLength;
  AverageWaveform.resize(AveWfmBins);
  AverageWaveform_count = 0; 

  // This specifies whether or not waveforms are to be summed
  // to eventually produce average waveform.  Intended to be
  // changed externally, ie:
 
  // Initialize some publically-accessible values
  prepulse_baseline = 0.;
  prepulse_rms = 0.;
  fit_ZeroPoint = 0.;
  fit_SlowNorm = 0.;
  fit_SlowTau = 0.;
  fit_ReducedChi2 = 0.;
  prepulse_x1 = 0;

  std::cout
  <<"Initializing OpHitBuilderAlg:\n"
  <<"  - make histograms     : "<<fMakeHistograms<<"\n"
  <<"  - use prepulse fit?   : "<<fUsePrepulseFit<<"\n"
  <<"  - prepulse tau limits : "<<fPrePulseTau1<<"-"<<fPrePulseTau2<<"\n"
  <<"  - hit finding mode    : "<<fHitFindingMode<<"\n"
  <<"  - grad hit thresh     : "<<fGradHitThresh<<"\n"
  <<"  - grad RMS thresh     : "<<fGradRMSThresh<<"\n";
  
  // Create directory for histograms
  if( fMakeHistograms ) {
    art::TFileDirectory histDir = tfs->mkdir("ophitbuilderalg");
    hGradient       = histDir.make<TH1D>("Gradient","Gradient",600,-30,30);
  }

}


//--------------------------------------------------------------  
//Destructor
OpHitBuilderAlg::~OpHitBuilderAlg()
{
}

//--------------------------------------------------------------
void OpHitBuilderAlg::reconfigure( fhicl::ParameterSet const& pset ){

  fTau                  = pset.get<float> ("Tau",-1300.);

  fReturnRawIntegrals   = pset.get< bool > ("ReturnRawIntegrals", false);
  fDAQModule            = pset.get< std::string >("DAQModule","daq");
  fInstanceName         = pset.get< std::string >("InstanceName","");
  fGradHitThresh        = pset.get< float >("GradHitThresh",-10);
  fSignalHitThresh      = pset.get< float >("SignalHitThresh",4);
  fPulseHitThreshLow    = pset.get< float >("PulseHitThreshLow",0.);
  fPulseHitThreshHigh   = pset.get< float >("PulseHitThreshHigh",999);
  fGradRMSThresh        = pset.get< float >("GradRMSThresh",5); 
  fMinHitSeparation     = pset.get< short >("MinHitSeparation",20);
  fBaselineWindowLength = pset.get< short >("BaselineWindowLength",1000);
  fPrePulseBaselineFit  = pset.get< short >("PrePulseBaselineFit",100);
  fPrePulseDisplay      = pset.get< short >("PrePulseDisplay",500);
  fPrePulseTau1         = pset.get< float >("PrePulseTau1",800.);
  fPrePulseTau2         = pset.get< float >("PrePulseTau2",1600.);
  fFullWindowLength     = pset.get< short >("FullWindowLength",7000);
  fIntegrationWindows   = pset.get< std::vector<short> >("IntegrationWindows",{100,7000});
  fMvPerADC             = pset.get< float >("MvPerADC",0.1953);
  fUsePrepulseFit       = pset.get< bool  >("UsePrepulseFit","true");
  fHitTimeCutoffLow     = pset.get< int   >("HitTimeCutoffLow",-100000);
  fHitTimeCutoffHigh    = pset.get< int   >("HitTimeCutoffHigh",100000);
  fHitFindingMode       = pset.get< std::string >("HitFindingMode","grad");
  fAddHitsToAverageWaveform = pset.get< bool >("AddHitsToAverageWaveform",false);
  fMakeHistograms       = pset.get< bool > ("MakeHistograms",true);

  fMskBaselineSubtr_waitprd = pset.get<int>   ("MskBaselineSubtr_waitprd",50);
  fMskBaselineSubtr_minseg  = pset.get<size_t>("MskBaselineSubtr_minseg",50);
  fMskBaselineSubtr_range   = pset.get<size_t>("MskBaselineSubtr_range",100);
  fMskBaselineSubtr_P       = pset.get<float> ("MskBaselineSubtr_P",0.15);
  fMskBaselineSubtr_grmsfac = pset.get<float> ("MskBaselineSubtr_grmsfac",3);
  fMskBaselineSubtr_adcthresh = pset.get<float>("MskBaselineSubtr_adcthresh",5);

}

//-------------------------------------------------------------------------
// Return specific OpDetPulse object from event
raw::OpDetPulse OpHitBuilderAlg::GetPulse( const art::Event& e, int opchannel){
  
  raw::OpDetPulse out;

  // Get the OpDetPulses of the photodetectors
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);
  
  for( size_t i = 0; i < WaveformHandle->size(); i++){
    art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,i);
    raw::OpDetPulse ThePulse = *ThePulsePtr;
    if( ThePulse.OpChannel() == opchannel ) {
      out = ThePulse;
      break;
    }
  }
 
 return out; 
}

//-------------------------------------------------------------------------
// Resets all of the member data
void OpHitBuilderAlg::Reset() {
  fBaseline = -999.;
  fRMS      = -999.;
}

//--------------------------------------------------------------------------
//GetHits:  The 'meat & potatoes' of OpHitFinding!  This function takes a whole
//          opDetPulse object and return a vector of hit times. Hits are found using 
//          a gradient-threshold method, and each is required to exceed some 
//          multiple of the local RMS of the gradient (calculated in a window 
//          preceding the hit) to reduce the frequency of fake hits.
//          
//          Parameters to set in your FCL.  The defaults have been chosen for
//          best effectiveness.  It's recommended these remain unchanged.
//              - GradHitThresh (default -10)
//              - PulseHitThresh (default 0 [mV])
//              - GradRMSThresh (default 5)
//              - MinHitSeparation (default 20 [ns])
//
std::vector<short> OpHitBuilderAlg::GetHits( raw::OpDetPulse &opdetpulse ) {
  std::vector<short> wfm = opdetpulse.Waveform();
  std::vector<float> wfm_float(wfm.begin(), wfm.end()); 
  size_t TriggerTime = (size_t)opdetpulse.FirstSample();
  return GetHits(wfm_float, TriggerTime);
}
std::vector<short> OpHitBuilderAlg::GetHits( std::vector<float>& wfm, size_t TriggerTime ) 
{

  // Get baseline and RMS of the waveform.  First, we want to check
  // that there are no pulses in the first part of the waveform where
  // the baseline will be calculated.  If there are, we should mask
  // these out.
  
  // First make the gradient
  std::vector<float> g = MakeGradient(wfm);
  
  // Create an empty vector to be filled with select
  // waveform values within the baseline region.  
  std::vector<short> BaselineWindow;
  BaselineWindow.reserve(fBaselineWindowLength);

  // After every gradient hit is found, skip this many
  // of the following samples
  size_t mask_interval = 100;
   
  // Look for hits in the gradient within this baseline window.
  // If a hit is found in the gradient, exclude next N samples
  // (where N is the defined by mask_interval above)
  for(size_t i = 0; i < (size_t)fBaselineWindowLength; i++){
    if ( g[i] <= fGradHitThresh ){
      i = std::min( (size_t)fBaselineWindowLength, i + mask_interval );
    }
    else {
      BaselineWindow.push_back(wfm[i]);
    }
  }
 
  // Now get baseline using this masked region as input 
  float baseline = GetBaselineAndRMS(BaselineWindow,0,BaselineWindow.size())[0];

  // Hit finding limits (+/- trigger time) set in fcl
  size_t t1 = std::max((int)TriggerTime + fHitTimeCutoffLow,0);
  size_t t2 = std::min((int)TriggerTime + fHitTimeCutoffHigh,(int)wfm.size());

  // Make vector hits to be filled
  std::vector<short> hits;
  std::vector<short> hits_filtered;

  // Begin hit-finding, discriminating on either the gradient or the signal
  if( fHitFindingMode=="grad" )
  {
    
    LOG_DEBUG("OpHitBuilder")
    << "Scanning for hits using gradient-threshold method: \n"
    << "(1st-pass thresh " << fGradHitThresh << " ADC/ns, RMS thresh x"<<fGradRMSThresh<<")";
  
    // Set rising_edge to false before starting hit finding loop
    bool rising_edge = false;
    
    // Scan over gradient
    for(size_t i = 0; i < g.size(); i++){
      
      // If gradient is under threshold, and we are within the hit finding
      // time limits set in the fhicl, we have a hit candidate
      if (g[i] <= fGradHitThresh && rising_edge == false && i > t1 && i < t2){
        rising_edge = true;
        hits.insert(hits.end(),i);
      }
         
      // When we leave the hit candidate, reset rising_edge to false
      if (g[i] > fGradHitThresh && rising_edge == true) rising_edge = false;
  
    } // <-- end scan over gradient
 
    LOG_DEBUG("OpHitBuilder")
    << "Found " << hits.size()
    << " hits in first gradient threshold pass"; 
 
    if( fGradRMSThresh > 0 ) { 
    // So now we have a list of hits, but clusters of these could be
    // fakes due to a "noisy" gradient. Filter out the hits that
    // are not high enough above the gradient's local RMS.
    
    // Loop over the hits
    for( size_t i = 0; i < hits.size(); i++ ) {
  
      // Find the window limits to use in calculating local RMS.  
      // For all hits after the first one, the window's start point 
      // will be limited to fMinHitSeparation after the first hit point
      short rms_window_size   = 100;
      short rms_window_end    = hits[i] - 5;
      short rms_window_start  = rms_window_end - rms_window_size;
      if( rms_window_start < 0 ) 
        rms_window_start = 0;
      if( i >= 1 && hits[0] + fMinHitSeparation > rms_window_start ) 
        rms_window_start = hits[0] + fMinHitSeparation;
      rms_window_size = rms_window_end - rms_window_start;
  
      // Find gradient's local RMS for the hit using this window
      std::vector<float> tmp = GetBaselineAndRMS( g, rms_window_start, rms_window_end );
      float g_mean = tmp[0];
      float g_rms  = tmp[1];
  
      // Find amplitude of gradient in neighborhood of this hit
      float g_amp = g_mean - GetLocalMinimum(g,hits[i]);
  
      // Calculate quick amplitude of pulse to use in discrimination
      float hit_amp   = (baseline-GetLocalMinimum(wfm,hits[i]))*fMvPerADC;
      
      LOG_DEBUG("OpHitBuilder")
      << "  - " << hits[i] << "(pre-pulse " << rms_window_start << "-" << rms_window_end << ")"
      << " gRMS " << g_rms << " (thresh " << g_rms*fGradRMSThresh
      << ", grad amp " << g_amp 
      << ", pulse amp " << hit_amp;
  
      // First check if hit's gradient and pulse amplitude exceed the set thresholds
      //if( (g_amp > g_rms*fGradRMSThresh) && (hit_amp > fPulseHitThreshLow) ){
      if( (g_amp > g_rms*fGradRMSThresh)){
        hits_filtered.push_back(hits[i]);
        LOG_DEBUG("OpHitBuilder") << "  --> hit passes.";
      }
    
    } // <-- end loop over hits
    } else {
      hits_filtered = hits; 
    }
    LOG_DEBUG("OpHitBuilder") << hits_filtered.size() << " hits pass filter.";
  
  } // endif GRAD mode
  
  else if ( fHitFindingMode == "signal" ){
    
    LOG_DEBUG("OpHitBuilder") 
    << "Scanning for hits using signal/pulse-threshold method: \n"
    << "(threshold " << fSignalHitThresh << ")";
    
    // Set rising_edge to false before starting hit finding loop
    bool rising_edge = false;
  
    for(size_t i = 0; i < wfm.size(); i++){
      
      // If signal is under threshold, and we are within the hit finding
      // time limits set in the fhicl, we have a hit candidate
      if ( (baseline - wfm[i])*fMvPerADC <= fSignalHitThresh && rising_edge == false && i > (size_t)t1 && i < (size_t)t2){
        rising_edge = true;
        hits_filtered.insert(hits_filtered.end(),i);
        LOG_DEBUG("OpHitBuilder")
        << " - " << i;
      }
         
      // When we leave the hit candidate, reset rising_edge to false
      if ( (baseline - wfm[i])*fMvPerADC > fSignalHitThresh && rising_edge == true){
        rising_edge = false;
      }
  
    } // <-- end scan over waveform
    LOG_DEBUG("OpHitBuilder") << "Found " << hits_filtered.size() << " hits.";
  
  } // Endif SIGNAL mode
  // End hitfinding
  
  // Now merge all remaining hits using shorter spacing
  std::vector<short> hits_merged = HitMerger(hits_filtered,fMinHitSeparation,1);
  LOG_DEBUG("OpHitBuilder") << "Post-merging: " << hits_merged.size() << " hits.";

  return hits_merged;
  
}



//--------------------------------------------------------------
// MakeGradient
std::vector<float> OpHitBuilderAlg::MakeGradient( const std::vector<short>& wfm )
{
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return MakeGradient(wfm_float);
}
std::vector<float> OpHitBuilderAlg::MakeGradient( const std::vector<float>& wfm )
{
  // Using "pseudotime" from time of flight alg 
  std::vector<float> g(wfm.size(),0.);
  for(size_t i=2; i<wfm.size()-2; i++){
    g[i] = (-wfm[i+2]+8.*wfm[i+1]-8*wfm[i-1]+wfm[i-2])/12.;
    if( fMakeHistograms ) hGradient->Fill(g[i]); 
  }
  return g;
}




//-------------------------------------------------------------
// Merge hits
std::vector<short> OpHitBuilderAlg::HitMerger( std::vector<short>& hits, short spacing, int option)
{
  LOG_DEBUG("OpHitBuilder")
  << "Merging hits... (" << hits.size() << ")";

  std::vector<short> hits_merged;
  for(size_t i=0; i<hits.size(); i++){
    // Always add the first hit
    if(i==0) {
      hits_merged.push_back(hits[i]);
    } else {    
      switch (option){
        // option 0: check only against first hit
        case 0:
          if( (hits[i]-hits[0]) >= spacing ) hits_merged.push_back(hits[i]);
          break;
        // option 1: check each hit against previous hit
        case 1:
          if( (hits[i]-hits_merged[hits_merged.size()-1]) >= spacing ) hits_merged.push_back(hits[i]);
          break;
      }
    }
  }
  return hits_merged;
}


//--------------------------------------------------------------
std::vector<float> OpHitBuilderAlg::GetBaselineAndRMS( const std::vector<short>& wfm, short x1, short x2 )
{
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return GetBaselineAndRMS(wfm_float,x1,x2);
}

std::vector<float> OpHitBuilderAlg::GetBaselineAndRMS( const std::vector<float>& wfm, short x1, short x2 )
{
  CalcBaselineAndRMS(wfm,x1,x2);
  std::vector<float> out(2);
  out[0] = fBaseline;
  out[1] = fRMS;
  return out;
}

//-------------------------------------------------------------
void OpHitBuilderAlg::CalcBaselineAndRMS( const std::vector<short>& wfm, short x1, short x2 ) { 
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  CalcBaselineAndRMS(wfm_float,x1,x2);
}
void OpHitBuilderAlg::CalcBaselineAndRMS( const std::vector<float>& wfm, short x1, short x2 ) {  
  float mean = 0;
  float sumSquares = 0;
  short N = x2 - x1;
  for ( short i = x1; i < x2; i++ ) mean += wfm[i]/N;
  for ( short i = x1; i < x2; i++ ) sumSquares += pow(wfm[i]-mean,2);
  fBaseline = mean;
  fRMS      = sqrt(sumSquares/N);
}


//------------------------------------------------------------
void OpHitBuilderAlg::SubtractBaseline( std::vector<float>& wfm ) {
  SubtractBaseline( wfm, fBaseline );
}
void OpHitBuilderAlg::SubtractBaseline( std::vector<float>& wfm, float bs ) {
  for(size_t i=0; i<wfm.size();i++) wfm.at(i) -= bs;
}

//-----------------------------------------------------------
float OpHitBuilderAlg::MeanInRange( const std::vector<float>& v, int x1, int x2){
  float mean = 0;
  int N = x2 - x1;
  for(int i=x1; i<x2; i++) mean += v.at(i) / N;
  return mean;
}

//-----------------------------------------------------------
void OpHitBuilderAlg::SmoothOutVector( std::vector<float>& v, int range ) {
  if( range == 0 ) return;
  std::vector<float> v_orig(v.begin(), v.end() );
  int N = v.size();
  for(int i=0; i<N; i++ ) {
    int k1 = std::max(0, i - range);
    int k2 = std::min(N, i + range);
    v.at(i) = MeanInRange(v_orig, k1, k2);
  }
}

//-----------------------------------------------------------
void OpHitBuilderAlg::RebinVector( std::vector<float>& v, int range ) {
  if( range == 0 ) return;
  int N = v.size();
  if( N % range != 0 ) return; 
  for(int i=0; i<N; i+=range){
    float mean = MeanInRange(v,i,i+range);
    for(int j=i; j<i+range; j++) v.at(j) = mean;
  }
}


// -------------------------------------------------------------
// Returns a vector<float> containing (1) the hit's amplitude, and 
// all the specified integrals
std::vector<float> OpHitBuilderAlg::GetHitInfo( const std::vector<float>& wfm, short hit, short prev_hit, std::vector<short>& windows){
  return GetHitInfo(wfm, hit, prev_hit, windows, fUsePrepulseFit);
}

std::vector<float> OpHitBuilderAlg::GetHitInfo( const std::vector<float>& wfm, short hit, short prev_hit, std::vector<short>& windows, bool usePrepulseFit)
{

  LOG_DEBUG("OpHitBuilder")
  << "GETHITINFO: Processing hit at sample " << hit << "(prev hit @ " << prev_hit << ")";

  int polarity = -1;

  PulseWidth = -99.;


  // Save number of samples in the waveform
  size_t nSamples = wfm.size();

  // Create vector to be returned (amplitude, integralWindow1, integralWindow2, ...)
  std::vector<float> hit_info(windows.size()+1, -9999.);
  
  // If the hit is too early to reliably calculate a baseline, OR if the 
  // previous hit is too close, stop now and return defaults
  if( (hit < 110)||(hit-prev_hit < 250 ) || ( hit > (short)nSamples-100) ) {
    LOG_DEBUG("OpHitBuilder") << "!!! Hit is too early or close to prev hit -- abort.";
    return hit_info;
  }

  // Get baseline
  size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hit-10));
  CalcBaselineAndRMS(wfm, 0, baseline_win_size );
  float baseline  = fBaseline;
  float rms       = fRMS;

  // Determine bounds for fit and integration 
  short x1,x1b,x2,x3;
  x1  = std::max(int(hit - fPrePulseBaselineFit),prev_hit+200);
  x1b = std::max(int(hit - fPrePulseDisplay),int(x1));
  x2  = hit - 5;
  x3  = std::min(int(hit + windows.at(windows.size()-1)),int(nSamples)-1);
  const int prepulse_bins = int(x2) - int(x1);
  const int total_bins    = int(x3) - int(x1);
  LOG_DEBUG("OpHitBuilder") 
  << "  x1, x2, x3 = " << x1 << "  " << x2 << "  " << x3;
  
  prepulse_baseline = 0.;

  // Find overall prepulse baseline/RMS
  CalcBaselineAndRMS(wfm,x1b,x2);
  prepulse_baseline = fBaseline;
  prepulse_rms      = fRMS;
  float diff       = prepulse_baseline - baseline;
  
  LOG_DEBUG("OpHitBuilder") 
  << "  Prepulse baseline (x1-x2) " << prepulse_baseline << ", rms " << prepulse_rms << "\n"
  << "  Waveform baseline         " << baseline <<", rms " << rms << "\n"
  << "  Difference                " << diff << " (" << diff/rms << " * rms)";

  // Get mean of first few samples of the prepulse region
  int  nn = 10;
  float x1_baseline = GetBaselineAndRMS(wfm,x1,x1+nn)[0];
  float norm_set    = (x1_baseline - baseline);
  if (norm_set > 0.) norm_set = 0.;
  
  // Define exponential function and initialize parameters
  TF1 pfit("pfit","[0] + [1]*exp(-(x-[2])/[3])",0.,30000.);
//  if(prepulse_rms <= rms*1.5) { pfit.FixParameter(0,prepulse_baseline);
//  } else {                      pfit.FixParameter(0,baseline);}
  pfit.FixParameter(0,baseline); 
  pfit.FixParameter(2,x1+nn/2); 
  pfit.SetParameter(1,norm_set);
  pfit.SetParLimits(1,-1200.,0.);
  if( fTau > 0 ) {
    pfit.FixParameter(3,fTau);
  } else {
    pfit.SetParameter(3,0.5*(fPrePulseTau1+fPrePulseTau2));
    pfit.SetParLimits(3,fPrePulseTau1,fPrePulseTau2);
  }

  if(usePrepulseFit){ 
    
    // Fill x,y arrays to be used in the TGraph fit
    std::vector<int> x(prepulse_bins,0);
    std::vector<int> y(prepulse_bins,0);
    for( int i = 0; i < prepulse_bins; i++) {
      x[i] = int(x1) + i;
      y[i] = int(wfm[x[i]]);
    }
      
    // Fit prepulse region with above function
    TGraph graph(prepulse_bins,x.data(),y.data());
    graph.Fit("pfit","QN0");

  } else {

    // Prepulse fit turned OFF; just use baseline
    // (zero out exponential component)
    pfit.FixParameter(1,0.); 
  
  }

  // Save fit info into publicly-accessible members   
  fit_SlowNorm        = pfit.GetParameter(1);
  fit_ZeroPoint       = pfit.GetParameter(2);
  fit_SlowTau         = pfit.GetParameter(3);
  int NDF             = pfit.GetNDF();
  fit_ReducedChi2     = pfit.GetChisquare()/float(NDF); 
  prepulse_x1         = x1;
 
  LOG_DEBUG("OpHitBuilder")
  << "  Resulting fit parameters \n"
  << "    norm            : " << fit_SlowNorm     << "\n"
  << "    tau             : " << fit_SlowTau      << "\n"
  << "    Chi2/NDF        : " << fit_ReducedChi2  << "\n"  
  << "  Fit function values (baseline subtracted)\n"
  << "    f(x1)  = " << pfit.Eval(x1) - baseline    << "\n"
  << "    f(x2)  = " << pfit.Eval(x2) - baseline    << "\n"
  << "    f(hit) = " << pfit.Eval(hit) - baseline   << "\n"
  << "    actual hit value = " << wfm[hit] - baseline           << "\n"
  << "    f(x3)  = " << pfit.Eval(x3) - baseline;

  // Save average waveform during integration?
  bool flag_ave = ( (fAddHitsToAverageWaveform)&&(total_bins>=AveWfmBins));
  if(flag_ave) AverageWaveform_count++;

  // Integrate using the fitted function as running baseline
  float integral  = 0.;
  float amplitude = -9999.;
  short iamp      = 0;
  int   iwin      = 0;
  
  std::vector<float> wfm_corr(wfm.size(),0.);
  for(size_t i=0; i<wfm.size(); i++){
    wfm_corr.at(i) = polarity*( wfm.at(i) - pfit.Eval(i) ); 
  }

  LOG_DEBUG("OpHitBuilder") 
  << "  Starting integration from "<<x2<<" to "<<x3;

  for( short xx = x1; xx < x3; xx++){
   
//    float yy   = CorrectWfmTF1( wfm, xx, pfit, -1);
    float yy = wfm_corr.at(xx);
     
    // Once past x2, start integration
    if( xx >= x2 ) {
      if( fReturnRawIntegrals ) 
        integral += polarity*( wfm.at(xx) - baseline ); 
      else 
        integral += yy;
    }

    // Save integral at appropriate window sizes
    //if( xx - x2 == windows[iwin]-1) {
    if( xx - hit == windows[iwin]-1) {
      hit_info[iwin+1] = integral;
      iwin++;
    }

    // Scan for amplitude when near the hit
    if ( ((hit-xx > -5)||(hit-xx < 100)) && (yy > amplitude) ) {
      amplitude   = yy;
      iamp        = xx;
    }

    // Add to average waveform if specified
    if( (flag_ave) && ( xx >= x1b)) AverageWaveform.at(xx-x1b) += yy*fMvPerADC;
  
  }
  
  // Save amplitude as first element of hit info vector
  hit_info[0] = amplitude;

  // Calculate pulse width, ie FWHM.
  // Start at the peak and walk backwards to find where signal
  // drops below half the amplitude, extrapolating between
  // neighboring samples for better accuracy.  Then, starting at
  // the peak again, walk forwards and find the corresponding point
  // on the high-side.
//  LOG_DEBUG("OpHitBuilderAlg") << "Determining pulse width, iamp = "<<iamp;
  short k1 = std::max(0, iamp - 50);
  short k2 = std::min((short)nSamples-2, iamp + 50);
  float p1 = 0.;
  float p2 = 0.;
  float thresh = amplitude*0.5;
  for(short i = iamp; i > k1; i--){
    float y1 = wfm_corr.at(i);
    if( y1 <= thresh ) {
      float y2 = wfm_corr.at(i+1);
      p1 = FindIntersection(i, y1, i+1, y2, thresh); 
      break;
    }
  }
//  std::cout<<"p1 = "<<p1<<"\n";
//  std::cout<<"Now working our way up to i = "<<k2<<"\n";
  for(short i = iamp; i < k2; i++){
    float y2   = wfm_corr.at(i);
    if( y2 <= thresh ) {
      float y1 = wfm_corr.at(i-1);
      p2 = FindIntersection(i-1, y1, i, y2, thresh); 
      break;
    }
  }
  PulseWidth = p2 - p1;

//  std::cout<<"Done\n";
   
  return hit_info;

}

std::vector<float> OpHitBuilderAlg::GetHitInfo( const std::vector<short>& wfm, short hit, short prev_hit, std::vector<short>& windows)
{
  return GetHitInfo(wfm, hit, prev_hit, windows, fUsePrepulseFit);
}
std::vector<float> OpHitBuilderAlg::GetHitInfo( const std::vector<short>& wfm, short hit, short prev_hit, std::vector<short>& windows, bool usePrepulseFit)
{
  std::vector<float> wfm2(wfm.begin(),wfm.end());
  return GetHitInfo(wfm2, hit, prev_hit, windows, usePrepulseFit);
}


//-------------------------------------------------------------
// Scans segment of vector<float> around a designated sample point
// and returns the local minimum
float OpHitBuilderAlg::GetLocalMinimum(const std::vector<float>& v, short hit)
{
  size_t r1 = 5;
  size_t rmin = 10;
  size_t r2 = 45; 
  float  low  = 9999;
  size_t x1 = hit-r1;
  size_t x2 = hit+r2;
  
  int counter = 0;

  // Go at LEAST "rmin" samples forward of hit, but no more than
  // r2.  If the value fails to exceed the amplitude for 5 
  // consecutive samples following the first 5 samples, or if
  // we've reached r2, end the scan.
  for( size_t j = x1; j<x2; j++){
    if( v[j] <= low ) {
      low = v[j];
      counter=0;
    }
    if( j > hit+rmin && v[j] > low ) counter++;

    if( counter==5 ) break;
  }

  return low;
}

short OpHitBuilderAlg::GetLocalMinimum(const std::vector<short>& v, short hit)
{
  std::vector<float> vd(v.begin(),v.end());
  return (short)GetLocalMinimum(vd,hit);
}


//-------------------------------------------------------------
// Scans segment of vector<float> around a designated sample point
// and returns the local minimum
float OpHitBuilderAlg::GetLocalMaximum(const std::vector<float>& v, short hit)
{
  size_t r1 = 5;
  size_t rmin = 10;
  size_t r2 = 45; 
  float  high  = -9999;
  size_t x1 = hit-r1;
  size_t x2 = hit+r2;
  
  int counter = 0;

  // Go at LEAST "rmin" samples forward of hit, but no more than
  // r2.  If the value fails to exceed the amplitude for 5 
  // consecutive samples following the first 5 samples, or if
  // we've reached r2, end the scan.
  for( size_t j = x1; j<x2; j++){
    if( v[j] >= high ) {
      high = v[j];
      counter=0;
    }
    if( j > hit+rmin && v[j] < high ) counter++;

    if( counter==5 ) break;
  }

  return high;
}

short OpHitBuilderAlg::GetLocalMaximum(const std::vector<short>& v, short hit)
{
  std::vector<float> vd(v.begin(),v.end());
  return (short)GetLocalMaximum(vd,hit);
}


//--------------------------------------------------------------
// Get only amplitude of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant amplitude information).
float OpHitBuilderAlg::GetHitAmplitude(const std::vector<short>& wfm, short hit,short prev_hit)
{
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return GetHitAmplitude(wfm_float, hit, prev_hit);
}

float OpHitBuilderAlg::GetHitAmplitude(const std::vector<float>& wfm, short hit,short prev_hit)
{
  std::vector<short> intWindows{100};
  return GetHitInfo(wfm,hit,prev_hit,intWindows)[0];
}


//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitPromptIntegral(const std::vector<short>& wfm, short hit, short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit,fIntegrationWindows)[1];
}



//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitFullIntegral(const std::vector<short>& wfm, short hit, short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit,fIntegrationWindows)[2];
}


//--------------------------------------------------------------
// Decide if a given OpDetPulse qualifies as a "clean beam waveform."
// ie, (1) Timestamp within beam window, and (2) exactly one optical hit
// that occurs at the trigger time (+/- 1% tolerance)
bool OpHitBuilderAlg::IsCleanBeamWaveform( raw::OpDetPulse &opdetpulse )
{
  std::vector<short> hits = GetHits(opdetpulse);
  std::vector<short> wfm = opdetpulse.Waveform();
  int NSamples = wfm.size();

  if( (hits.size() == 1) && ( fabs(int(hits[0]) - opdetpulse.FirstSample()) <= 0.01*(float)NSamples) ){
    return true;
  } else {
    return false;
  }
}



//--------------------------------------------------------------------------
// Get pedestal.  This algorithm first calls the much-faster GetBaselineAndRMS 
// function in order to set the limits of the histogram to be used in proper 
// pedestal calculation.  Fitting to the histogram is more computationally taxing, 
// so unless the baseline needs to be known very precisely, it is more efficient 
// to use the standard GetBaselineAndRMS().
std::vector<float> OpHitBuilderAlg::GetPedestalAndRMS( const std::vector<float>& wfm, short x1, short x2)
{
  std::vector<float> out(2);
  std::vector<float> tmp(2); 
  tmp = GetBaselineAndRMS( wfm, x1, x2);
  float baseline  = tmp[0];
  float rms      = tmp[1]; 
  
  TH1F h("h","h",100,baseline-5.*rms,baseline+5.*rms);
  for(int i=x1; i<x2; i++) h.Fill(wfm[i]); 

  TF1 f("f","gaus");
  f.SetParameter(0,h.GetMaximum());
  f.SetParameter(1,baseline);
  f.SetParameter(2,rms);
  
  h.Fit("f","QN0");
  float ped   = f.GetParameter(1);
  float sigma = f.GetParameter(2);
  //float rchi2 = f.GetChisquare()/(float)(f.GetNDF()-1);
  
  out[0] = ped;
  out[1] = sigma;
  return out;
}

// Get pedestal (for vector<short> input)
std::vector<float> OpHitBuilderAlg::GetPedestalAndRMS( const std::vector<short>& v, short x1, short x2)
{
  std::vector<float> wfm(v.begin(),v.end());
  return GetPedestalAndRMS(wfm,x1,x2);
}


//--------------------------------------------------------------------------
// Performs a running baseline subtraction
void OpHitBuilderAlg::SubtractRunningBaseline(const std::vector<float>& wvform, std::vector<float> &wvformout, const size_t nsamples, const size_t range)
{
  if( wvformout.size() != wvform.size() ) return; 
  for(size_t i=0; i<nsamples; i++) {
    int kstart  = (((int)i-(int)range)/2 > 0 ) ? i-range/2 : 0;
    int kend	  = (i+range/2 >= nsamples ) ? nsamples : i+range/2;
    double kaverage=0;
    for(int k=kstart;k<kend;k++) kaverage += wvform[k];
    kaverage	  /=(double)(kend-kstart);
    wvformout[i]  = wvform[i]-kaverage;
  }
}

void OpHitBuilderAlg::SubtractRunningBaseline(const std::vector<short>& wvform, std::vector<float> &wvformout, const size_t nsamples, const size_t range)
{  
  std::vector<float> wvform_float(wvform.begin(), wvform.end());
  SubtractRunningBaseline(wvform_float, wvformout, nsamples, range);
}

//--------------------------------------------------------------------------
// Performs a "masked" baseline subtraction using the gradient information to 
// zero out uninteresting portions of the waveform in an attempt to remove
// oscillation noise.
void OpHitBuilderAlg::MaskedBaselineSubtraction(const std::vector<float>& wvform, std::vector<float> &wvformout){
 
  // Params
  size_t nsamples = wvform.size(); 
  float gradRmsFactor = fMskBaselineSubtr_grmsfac;
  float adcThresh = fMskBaselineSubtr_adcthresh;
  if ( gradRmsFactor < 0 )  gradRmsFactor = 1e6;
  if ( adcThresh < 0 )      adcThresh = 1e6;
  
  // Create "baseline" vector
  std::vector<float> bs(nsamples, 0.);
  std::vector<bool>  isBs(nsamples, false);
  size_t nQuietSamples = 0;
 
  // Get gradient and find its baseline / RMS 
  std::vector<float> g = MakeGradient(wvform);
  CalcBaselineAndRMS( g, 0, fBaselineWindowLength );
  float gradMean  = fBaseline;
  float gradRMS   = fRMS;
  LOG_DEBUG("OpHitBuilderAlg")
  <<"MASKED BASELINE SUBTRACTION, grad mean/rms = "<<gradMean<<", "<<gradRMS;
  
  int counter = 0;
  for(size_t i=0; i<wvform.size(); i++){
    
    // whenever a jump in the gradient is detected, or an actual
    // hit is found, update the waiting period and negate the
    // last 10 baseline samples
    if( fabs(g[i]) > gradRmsFactor*gradRMS || fabs(wvform.at(i)) > adcThresh ){
      if( counter < fMskBaselineSubtr_waitprd ) counter = fMskBaselineSubtr_waitprd;
      size_t k1 = std::max((size_t)0, i-10);
      for(size_t ii=k1; ii<i; ii++) isBs.at(ii) = false;
    }
    
    // whenever the counter is zero (ie, we are in a "quiet" part
    // of the waveform), add to baseline vector
    if( counter == 0 ) {
      nQuietSamples++;
      isBs.at(i)  = true;
    }
  
    if( counter > 0 ) counter--;
  
  }
  
  // remove baseline regions that are small in length
  if( nQuietSamples != nsamples ) {
    LOG_DEBUG("OpHitBuilderAlg")<<"Removing short baseline segments...";
    for(size_t i=0; i<nsamples; i++){
      if( isBs.at(i) ) {
        for(size_t j=i; j<nsamples; j++){
          if( !isBs.at(j) || j == nsamples - 1 ) {
            if( j-i < fMskBaselineSubtr_minseg ) {
              for(size_t k=i; k<j; k++) isBs.at(k) = false;
            }
            i=j;
            break;
          }
        }
      }
    }
  }
  
  // start and end of waveform must be "baseline" to 
  // prevent nasty edge effect crashes later on
  for(size_t i=0; i<fMskBaselineSubtr_minseg; i++){
    isBs.at(i) = true;
    isBs.at(nsamples-1-i) = true;
  }

  LOG_DEBUG("OpHitBuilderAlg")
  <<"Done creating masked baseline -- "<<nQuietSamples<<" out of "<<nsamples<<" total samples included";
  
  // Smooth out the quiet areas with a running baseline subtraction,
  // using a truncated mean approach.
  for(size_t i=0; i<nsamples; i++){
    if( !isBs.at(i) ) continue;
    // create the vector of values to use
    int kstart  = (((int)i-(int)fMskBaselineSubtr_range)/2 > 0 ) ? i-fMskBaselineSubtr_range/2 : 0;
    int kend	  = (i+fMskBaselineSubtr_range/2 >= nsamples ) ? nsamples : i+fMskBaselineSubtr_range/2;
    std::vector<float> vals;
    for(int k=kstart; k<kend; k++){
      if( isBs.at(k) ) vals.push_back(wvform[k]);
    }
    bs.at(i) = CalcTruncatedMean(vals, fMskBaselineSubtr_P);
    //bs.at(i) = wvform[i];
  }

  LOG_DEBUG("OpHitBuilderAlg")
  <<"Done smoothing out baseline regions using truncated mean";

  // Extrapolate between gaps
  for(size_t i=0; i<nsamples; i++){ 
    if( !isBs.at(i) ) {
      float x1 = float(i-1);
      float x2 = 0.; 
      float y1 = bs.at(i-1);
      float y2 = 0.;
      
      // move forward until we find the next baseline region
      for(size_t ii=i; ii<nsamples; ii++){
        if( isBs.at(ii) ) {
          x2 = float(ii);
          y2 = bs.at(ii);
            
          /*  
            float m = (y2 - y1) / (x2 - x1);
            for(size_t j=i; j<ii; j++){
              bs.at(j) = y1 + m*(j-x1);
            }
            */
          
          // apply linear function between;
          // if gap is longer than 300, just
          // assume BS = 0
          if( x2-x1 < 500 ) {
            float m = (y2 - y1) / (x2 - x1);
            for(size_t j=i; j<ii; j++){
              bs.at(j) = y1 + m*(j-x1);
            }
          } else {
            // avoid artificially creating sudden jumps in waveform that could 
            // trick gradient-based hit finder; create 20-sample-long "ramp" 
            // on both ends
            size_t nramp = 20;
            for(size_t j=0; j<nramp; j++){
              bs.at(i+j)    = y1 - (j+1)*(y1/nramp);
              bs.at(ii-j-1) = y2 - (j+1)*(y2/nramp);
            }
          }
          
          i = ii;
          break;
        
        }
      }
    }
  }

  for(size_t i=0; i<nsamples; i++){
    wvformout.at(i) = bs.at(i);
   
    /* 
    if( isBs.at(i) ) {
      wvformout.at(i) = wvform.at(i);
    } else {
      wvformout.at(i) = -99999.;
    }
    */

  }
  
}

 
//-------------------------------------------------------------------------
float OpHitBuilderAlg::CalcTruncatedMean( std::vector<float>& v, float p){
  return CalcTruncatedMean(v, p, 0);
}
float OpHitBuilderAlg::CalcTruncatedMean( std::vector<float>& v, float p, int skew = 0){

 // Sort the vector
 if( p > 0 ) std::sort( v.begin(), v.end() );

  int N = v.size();

  // Calculate truncated mean
  int ntrim = p*N;
  if( N < 3 ) ntrim = 0;
  float sum = 0.;
  int NN = 0;
  int low = ntrim;
  int high = N - ntrim;
  if( skew == -1 ) low  = 0;
  if( skew ==  1 ) high = N;
  for(int i=0; i < N; i++){
    if( i >= low && i < high ){
      sum += v.at(i);
      NN++;
    }
  }

  if(NN > 0){ 
    return sum / NN;
  } else {
    return 0.;
  }

}

//--------------------------------------------------------------------------
// Experimental overshoot correction for Hamamatsu PMT in Run IIb+.
// A falling exponential is applied after each found "hit", whose parameters
// are based on the pulse amplitude.  Parameterization found by Jose Ignacio 
// Cevallos Aleman.
//
//    f(t) = A*exp(-t/B)
//    B = 11496 ns
//
//    Using pulse amplitude:
//      A(amp) = 0.3324 + 0.005689*amp
//
//    Using pulse integral,
//    50 samples, -5 to 45 relative to hit time:
//      A(int) = 0.2777 + 0.0002466*int
//
//void OpHitBuilderAlg::CorrectWfmOvershoot(const std::vector<float>& wfm_in, std::vector<float>& wfm_out, std::vector<short> hits, std::string mode="int") 
void OpHitBuilderAlg::CorrectWfmOvershoot(std::vector<float>& wfm_in, std::vector<short> hits, std::string mode="int") 
{
//  if( wfm_out.size() != wfm_in.size() ) return;

  // Define hit number limit
  const size_t maxHits = 50;
  size_t kmax = std::min(hits.size(), maxHits); 

  // Time constant of overshoot
  const float tau = 11496.;

  // Linear function for calculating "A"
  TF1 fA("fA","[0]+[1]*x",0.,1500.);
  if( mode == "amp" ) {
    fA.SetRange(0., 1500.);
    fA.FixParameter(0, 0.3324);
    fA.FixParameter(1, 0.005689);
  } else 
  if ( mode == "int" ) {
    fA.SetRange(0., 25000.);
    fA.FixParameter(0, 0.2777);
    fA.FixParameter(1, 0.0002466);
  }
  
  LOG_DEBUG("OpHitBuilderAlg")<<"Correcting waveform overshoot -- found "<<kmax<<" hits:";

//  std::cout<<"Hits.size = "<<hits.size()<<"\n";

  // Calculate hit info and save into vector
  std::vector<float> vx;
  std::vector<float> vt;
  std::vector<short> intWindows{50};
  for(size_t i=0; i<kmax; i++){
  // std::cout<<"  "<<i<<"   "<<hits.at(i)<<"\n";
    std::vector<float> hit_info(intWindows.size()+1,-999.); 
    if( i == 0 )  hit_info = GetHitInfo(wfm_in,hits.at(i),0,intWindows,false);
    else          hit_info = GetHitInfo(wfm_in,hits.at(i),hits.at(i-1),intWindows,true);
 //   std::cout<<"Done getting hit info, "<<hit_info[0]<<", "<<hit_info[1]<<"\n";
    if( mode == "amp" && hit_info[0] > 0. ) {vx.push_back(hit_info[0]); vt.push_back(hits.at(i));}
    if( mode == "int" && hit_info[1] > 0. ) {vx.push_back(hit_info[1]);  vt.push_back(hits.at(i));}
  }

  // redefine kmax
  kmax = vt.size();
  
  // Create array of TF1s
  TF1* f[maxHits];
  char buffer[200];
  for(size_t i=0; i < kmax; i++){
    sprintf(buffer,"f%lu",i);
    f[i] = new TF1(buffer,"[0]*exp(-(x-[1])/[2])",0.,float(wfm_in.size()));
    f[i]->SetParameter(0, fA.Eval(vx.at(i)) );
    f[i]->SetParameter(1, vt.at(i) );
    f[i]->SetParameter(2, tau);
    //LOG_VERBATIM("OpHitBuilderAlg")<<"  f"<<i<<": A = "<<f[i]->GetParameter(0);
  }

  // Loop through input waveform and create new corrected waveforms
  for(size_t i=0; i<wfm_in.size(); i++) {
    //wfm_out.at(i) = wfm_in.at(i);
    for(size_t iHit=0; iHit<kmax; iHit++){
      if( (short)i >= vt.at(iHit) ) wfm_in.at(i) -= f[iHit]->Eval( i );
    }
  }
}


//--------------------------------------------------------------------------
float OpHitBuilderAlg::CorrectWfmTF1(const std::vector<float>& wfm, short i, const TF1& f, int polarity){
  return polarity*( wfm[i] - f.Eval(i) ); 
}

//-------------------------------------------------------------------------
float OpHitBuilderAlg::FindIntersection(float x1, float y1, float x2, float y2, float thresh){
    float dx = x2 - x1;
    if( dx <= 0 ) return 0.;
    float m = (y2-y1) / dx;
    float b = y1 - m*x1;
    return (thresh - b)/m;
}


//--------------------------------------------------------------------------
// eventType: Classify event based on its timestamp.
std::string OpHitBuilderAlg::eventType(float T){
  if      ( T >= 0.   &&  T < 1.2 ) {return "pedestal";}
  else if ( T >= 1.2  &&  T < 5.5 ) {return "beam";}
  else if ( T >= 5.5 ) {return "cosmic";}
  else                              {return "unknown";} 
}

//--------------------------------------------------------------------------
// eventType: filter events.  Returns TRUE if Timestamp matches any of
// the input categories 
bool OpHitBuilderAlg::eventTypeFilter(float T, std::vector<std::string> categories){
  if( categories.size() == 0 ) {  return true; }
  else {
    bool out = false;
    for(size_t i = 0; i < categories.size(); i++){
      if(   eventType(T)  == categories[i]
        ||  categories[i] == "all"          ) out = true; }
    return out;
  }
}

