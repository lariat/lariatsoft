//////////////////////////////////////////////////////////////////////
//                                                                  
// These are functions used to process optical detector information
// (particulary PMTs) and find/integrate hits.  They've been used
// mostly to analyze events in the stopping/decaying cosmic muon 
// samples, but will be adapted for more general use.
//
// Oct 2, 2015:
// The code is currently being tested and optimized, but feel free to
// to try out the algorithms if you'd like.  Comments prior to each
// function specify its use.  More complete documentation coming soon.
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

// LArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "RawData/TriggerData.h"

//LAriatSoft Includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"

// C++ includes
#include <iostream>
#include <cmath>
#include <cstdlib>

// ROOT includes
#include <TF1.h>
#include <TGraph.h>


//--------------------------------------------------------------
//Constructor
OpHitBuilderAlg::OpHitBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;
 
  // Set size of vector to hold summed waveforms 
  AverageWaveform.resize(fPrePulseDisplay + fFullWindowLength);
  AverageWaveform_count = 0; 
  AddHitToAverageWaveform = 0;
  
  // Initialize some values
  prepulse_baseline = 0.;
  prepulse_rms = 0.;
  fit_FastNorm = 0.;
  fit_FastTau = 0.;
  fit_SlowNorm = 0.;
  fit_SlowTau = 0.;
  fit_ReducedChi2 = 0.;
  
}


//--------------------------------------------------------------  
//Destructor
OpHitBuilderAlg::~OpHitBuilderAlg()
{
}


//--------------------------------------------------------------
void OpHitBuilderAlg::reconfigure( fhicl::ParameterSet const& pset ){
  fGradHitThresh       = pset.get< float >("GradHitThresh",-10);
  fPulseHitThreshLow       = pset.get< float >("PulseHitThreshLow",0.);
  fPulseHitThreshHigh      = pset.get< float >("PulseHitThreshHigh",999);
  fPulseHitRMSThresh        = pset.get< float>("PulseHitRMSThresh",2.5);
  fGradientRMSFilterThreshold = pset.get< float >("GradientRMSFilter",5); 
  fMinHitSeparation           = pset.get< short >("MinHitSeparation",20);
  fFirstHitSeparation         = pset.get< short >("FirstHitSeparation",250);
  fBaselineWindowLength       = pset.get< short >("BaselineWindowLength",1000);
  fPrePulseBaselineFit        = pset.get< short >("PrePulseBaselineFit",2000);
  fPrePulseDisplay        = pset.get< short >("PrePulseDisplay",500);
  fPrePulseTau1               = pset.get< float>("PrePulseTau1",1400.);
  fPrePulseTau2               = pset.get< float>("PrePulseTau1",1600.);
  fPromptWindowLength         = pset.get< short >("PromptWindowLength",100);
  fFullWindowLength           = pset.get< short >("FullWindowLength",7000);
  fMvPerADC                   = pset.get< float>("MvPerADC",0.2);
  fUsePrepulseFit             = pset.get< bool >("UsePrepulseFit","true");
  fTimestampCut               = pset.get< float>("TimestampCut",5.25);
  bVerbose                    = pset.get< bool >("Verbosity","true");
  fHitTimeCutoffLow       = pset.get< int>("HitTimeCutoffLow",-100000);
  fHitTimeCutoffHigh      = pset.get< int>("HitTimeCutoffHigh",100000);
  fSER_PreWindow          = pset.get<short>("SER_PreWindow",5);
  fSER_PostWindow         = pset.get<short>("SER_PostWindow",25);
}


//--------------------------------------------------------------------------
//GetHits:  The 'meat & potatoes' of OpHitFinding!  This function takes a 
//          waveform and returns a vector of hit times. Hits are found using 
//          a gradient-threshold method, and each is required to exceed some 
//          multiple of the local RMS of the gradient (calculated in a window 
//          preceding the hit) to reduce the frequency of fake hits.
//          
//          Parameters to set in your FCL.  The defaults have been chosen for
//          best effectiveness.  It's recommended these remain unchanged.
//              - GradientHitThreshold (default -10)
//              - PulseHitThreshold (default 0 [mV])
//              - GradientRMSFilterThreshold (default 5)
//              - MinHitSeparation (default 50 [ns])
//
std::vector<short> OpHitBuilderAlg::GetHits( raw::OpDetPulse opdetpulse) 
{
  std::vector<short> wfm = opdetpulse.Waveform();
  int TriggerTime = (int)opdetpulse.FirstSample();
  int t1 = std::max(TriggerTime + fHitTimeCutoffLow,0);
  int t2 = std::min(TriggerTime + fHitTimeCutoffHigh,(int)wfm.size());

  // Make vector hits to be filled
  std::vector<short> hits;

  // Make the gradient
  std::vector<float> g = MakeGradient(wfm);
  
  // Set rising_edge to false before starting hit finding loop
  bool rising_edge = false;

  for(size_t i = 0; i < g.size(); i++){
    
    // If gradient is under threshold, and we are within the hit finding
    // time limits set in the fhicl, we have a hit candidate
    if (g[i] <= fGradHitThresh && rising_edge == false && i > (size_t)t1 && i < (size_t)t2){
      rising_edge = true;
      hits.insert(hits.end(),i);
    }
       
    // When we leave the hit candidate, reset rising_edge to false
    if (g[i] > fGradHitThresh && rising_edge == true){
      rising_edge = false;
    }

  } // <-- end scan over gradient

  if(bVerbose) std::cout<<"Found "<<hits.size()<<" hits in first gradient theshold pass\n";

  // So now we have a list of hits, but clusters of these could be
  // fakes due to a "noisy" gradient. Filter out the hits that
  // are not high enough above the gradient's local RMS.
  
  std::vector<short> hits_filtered;
  
  for( size_t i = 0; i < hits.size(); i++ ) {

    short     rms_window_size = 100;
    short     rms_window_start = 0;
    short     rms_window_end  = rms_window_size;

    // Find the window limits to use in calculating local RMS.  
    // For all hits after the first one, the window's start point 
    // will be limited to 20ns after the first hit point (for now)
    rms_window_end    = hits[i]-5;
    rms_window_start  = rms_window_end - rms_window_size;
    if( rms_window_start < 0 ) rms_window_start = 0;
    if( i >= 1 ){
      if( hits[0] +10  > rms_window_start) rms_window_start = hits[0]+10;
    }
    rms_window_size = rms_window_end - rms_window_start;

    // Find gradient's local RMS for the hit using this window
    std::vector<float> tmp = GetBaselineAndRMS( g, rms_window_start, rms_window_end );
    float    g_mean = tmp[0];
    float    g_rms = tmp[1];


    // Find amplitude of gradient in neighborhood of this hit
    float g_amp = g_mean - GetLocalMinimum(g,hits[i]);

    // Calculate quick amplitude of pulse to use in discrimination
    size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hits[i]-10));
    tmp = GetBaselineAndRMS(wfm,0,baseline_win_size);
    float baseline = tmp[0];
    float hit_amp = (baseline-GetLocalMinimum(wfm,hits[i]))*fMvPerADC;
    
    if(bVerbose){
      std::cout<<"     "<<hits[i]<<": window: ("<<rms_window_start<<","<<rms_window_end<<")  rms: "<<g_rms
      <<"   rms_thresh: "<<g_rms*fGradientRMSFilterThreshold<<"   g_amp: "<<g_amp<<"  amp "<<hit_amp<<"\n";
    }

    // First check if hit's gradient and pulse amplitude exceed the set thresholds
    if( (g_amp > g_rms*fGradientRMSFilterThreshold) && (hit_amp > fPulseHitThreshLow) && (hit_amp < fPulseHitThreshHigh)){
      hits_filtered.push_back(hits[i]);
      if(bVerbose) std::cout<<"Hit passes!\n";
    }
  
  } // <-- end loop over hits

  if(bVerbose) std::cout<<"Post-filter: "<<hits_filtered.size()<<" hits\n";
  // Merge float-hits within 250 ns of first pulse
  std::vector<short> hits_merged = HitMerger(hits_filtered,fFirstHitSeparation,0);

  // Now merge all remaining hits using shorter spacing
  std::vector<short> hits_remerged = HitMerger(hits_merged,fMinHitSeparation,1);

  if(bVerbose) std::cout<<"Post-merging: "<<hits_remerged.size()<<" hits\n";
  return hits_remerged;
  
}


//--------------------------------------------------------------
// MakeGradient
std::vector<float> OpHitBuilderAlg::MakeGradient( std::vector<short> wfm )
{
  std::vector<float> g(wfm.size());
  g[0]=0; 
  g[1]=0; 
  for(size_t i=2; i<wfm.size(); i++) g[i] = float(wfm[i] - wfm[i-2])*0.5;
  return g;
}


//-------------------------------------------------------------
// Merge hits
std::vector<short> OpHitBuilderAlg::HitMerger( std::vector<short> hits, short spacing, int option)
{
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
          if( (hits[i]-hits[i-1]) >= spacing ) hits_merged.push_back(hits[i]);
          break;
      
      }
    }
  }

  return hits_merged;
}


//--------------------------------------------------------------
// Get (baseline,rms) of a segment of an std::vector<short>
std::vector<float> OpHitBuilderAlg::GetBaselineAndRMS( std::vector<short> wfm, short x1, short x2 )
{
  // Convert to vector<float> and send into
  // other version of function
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return GetBaselineAndRMS(wfm_float,x1,x2);
}


//--------------------------------------------------------------
// Get (baseline,rms) of a segment of an std::vector<float>
std::vector<float> OpHitBuilderAlg::GetBaselineAndRMS( std::vector<float> wfm, short x1, short x2 )
{
  float mean = 0;
  float sumSquares = 0;
  int N = x2 - x1;

  for ( int i = x1; i < x2; i++ ) mean += wfm[i]/N;
  for ( int i = x1; i < x2; i++ ) sumSquares += pow(wfm[i]-mean,2);
  
  std::vector<float> out(2);
  out[0] = mean;
  out[1] = sqrt(sumSquares/N);

  return out;
}


// -------------------------------------------------------------
// Returns a vector<float> containing (1) the hit's amplitude,
// (2) its prompt integral, and (3) its full integral
std::vector<float> OpHitBuilderAlg::GetHitInfo( std::vector<short> wfm, short hit, short prev_hit)
{

  //if(bVerbose) std::cout<<"GetHitInfo(wfm,"<<hit<<")"<<std::endl;

  // Create vector to be returned 
  std::vector<float> hit_info {-9.,-99.,-999.};
  
  // If the hit is too early to reliably calculate 
  // a baseline, OR if the previous hit is too close,
  // stop now and return zeros
  if( (hit < 100)||(hit-prev_hit < 100 ) ) return hit_info;

  // Get baseline
  size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hit-10));
  std::vector<float> tmp = GetBaselineAndRMS(wfm,0,baseline_win_size);
  float baseline = tmp[0];
  float rms      = tmp[1];

  // Determine bounds for fit and integration 
  short x1,x1b,x2,x3;
  x1  = std::max(int(hit - fPrePulseBaselineFit),prev_hit+50);
  x1b = std::max(int(hit - fPrePulseDisplay),int(x1));
  x2  = hit - 10;
  x3  = std::min(int(hit + fFullWindowLength),int(wfm.size()));
  const int prepulse_bins = int(x2) - int(x1);
  const int total_bins    = int(x3) - int(x1);
  if(bVerbose) {
    std::cout<<"Integration limits: x1,x2,x3 "<<x1<<"   "<<x2<<"  "<<x3<<std::endl;
    std::cout<<"Prev hit: "<<prev_hit<<"\n";
  }
  
  // Fill x,y arrays to be used in the TGraph
  // during the fit later
  int x[prepulse_bins],y[prepulse_bins];
  prepulse_baseline = 0.;
  for( int i = 0; i < prepulse_bins; i++) {
    x[i] = int(x1) + i;
    y[i] = int(wfm[x[i]]);
  }

  // Find overall prepulse baseline/RMS
  std::vector<float> prepulse_info = GetBaselineAndRMS(wfm,x1b,x2);
  prepulse_baseline = prepulse_info[0];
  prepulse_rms      = prepulse_info[1];
  float diff       = prepulse_baseline - baseline;
  
  if(bVerbose){
    std::cout<<"Prepulse baseline "<<prepulse_baseline<<"   rms "<<prepulse_rms<<std::endl;
    std::cout<<"Waveform baseline "<<baseline<<"   rms "<<rms<<std::endl;
    std::cout<<"Difference        "<<diff<<"   ("<<diff/rms<<" * rms)"<<std::endl;
  }

  // Get mean of first few samples of the prepulse region
  int  nn = 10;
  float x1_baseline = GetBaselineAndRMS(wfm,x1,x1+nn)[0];
  //bool PrepulseIsFlat = 0;
  // Check if prepulse baseline is flat.  If so,
  // use this as the "local" baseline and force
  // p1=0 in fit
  //if( prepulse_rms <= 2. ){ 
  //  std::cout<<"Prepulse looks flat... "<<prepulse_rms<<std::endl;
  //  baseline = prepulse_baseline;
  //  PrepulseIsFlat = "true";
  //}

  float norm_set    = (x1_baseline - baseline);
  if(bVerbose) std::cout<<"Norm set: "<<norm_set<<std::endl;
  if (norm_set > 0.) norm_set = 0.;
  
  // Two-component version of fit (for future use)
  // TF1 prepulse_exp_fit("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])+[4]*exp(-(x-[2])/[5])",0.,30000.);
  // prepulse_exp_fit.SetParameter(4,-0.1);
  // prepulse_exp_fit.SetParLimits(4,-200.,0.);
  // prepulse_exp_fit.SetParameter(5,7.);

  TF1 prepulse_exp_fit("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])",0.,30000.);
  
  prepulse_exp_fit.FixParameter(0,baseline);
  prepulse_exp_fit.FixParameter(2,x1+nn/2); 
  prepulse_exp_fit.SetParameter(1,norm_set);
  prepulse_exp_fit.SetParLimits(1,-200.,0.);
  prepulse_exp_fit.SetParameter(3,0.5*(fPrePulseTau1+fPrePulseTau2));
  prepulse_exp_fit.SetParLimits(3,fPrePulseTau1,fPrePulseTau2);

  if(fUsePrepulseFit){ 
      TGraph graph(prepulse_bins,x,y);
      graph.Fit("prepulse_exp_fit","QN0");
  } else {
    // Prepulse fit turned OFF; just use baseline
    prepulse_exp_fit.FixParameter(1,0.); // Zero out exponential component
  }

  // Save fit info into publicly-accessible members   
  fit_SlowNorm        = prepulse_exp_fit.GetParameter(1);
  fit_SlowTau         = prepulse_exp_fit.GetParameter(3);
  //fit_FastNorm      = prepulse_exp_fit.GetParameter(4);
  //fit_FastTau       = prepulse_exp_fit.GetParameter(5);
  int NDF             = prepulse_exp_fit.GetNDF();
  fit_ReducedChi2     = prepulse_exp_fit.GetChisquare()/float(NDF); 
  
  if(bVerbose){
    std::cout<<"Parameters of fitted prepulse region: \n";
    std::cout<<"    Norm_slow           : "<<fit_SlowNorm<<"\n";
    std::cout<<"    Tau_slow            : "<<fit_SlowTau<<"\n";
    //std::cout<<"    Norm_fast           : "<<fit_FastNorm<<"\n";
    //std::cout<<"    Tau_fast            : "<<fit_FastTau<<"\n";
    std::cout<<"    waveform baseline   : "<<baseline<<"\n";
    std::cout<<"    BS at x1            : "<<x1_baseline<<"\n";
    std::cout<<"    actual hit value    : "<<wfm[hit]<<"\n";
    std::cout<<"    f(x1)               : "<<prepulse_exp_fit.Eval(x1)<<" ("<<prepulse_exp_fit.Eval(x1)-baseline<<" from BS)\n";
    std::cout<<"    f(x2)               : "<<prepulse_exp_fit.Eval(x2)<<" ("<<prepulse_exp_fit.Eval(x2)-baseline<<" from BS)\n";
    std::cout<<"    f(x3)               : "<<prepulse_exp_fit.Eval(x3)<<" ("<<prepulse_exp_fit.Eval(x2)-baseline<<" from BS)\n";
    std::cout<<"    ChiSquared          : "<<fit_ReducedChi2<<"\n";
  }
  
  // Save average waveform during integration?
  bool flag_ave = ( (AddHitToAverageWaveform)&&(total_bins>=fPrePulseDisplay+fFullWindowLength));
  if(flag_ave) AverageWaveform_count++;
  
  // Integrate using the fitted function as running baseline
  float integral = 0, integral_prompt = 0;
  float amplitude = 0.;
  for( int i = 0; i < total_bins; i++){
    
    short xx    = x1 + short(i);
    float yy   = prepulse_exp_fit.Eval(xx) - (float)wfm[xx];  
    
    // Once past x2, start integration
    if( xx >= x2 ) integral += yy;

    // Save prompt integral
    if( xx - hit == fPromptWindowLength) integral_prompt = integral;

    // Scan for amplitude when near the hit
    if ( (abs(hit-xx) < 30) && (yy > amplitude) ) amplitude = yy;

    // Add to average waveform if specified
    if( (flag_ave) && ( xx >= x1b)) AverageWaveform.at(xx-x1b) += yy*fMvPerADC;
 
    //std::cout<<"   i "<<xx<<"  v:"<<yy*fMvPerADC
    //<<"  integrating? "<<(xx>=x2)
    //<<"  "<<integral<<"   "<<integral_prompt<<std::endl;
  
  }

  hit_info[0] = amplitude*fMvPerADC;
  hit_info[1] = integral_prompt;
  hit_info[2] = integral;

  return hit_info;

}


//-------------------------------------------------------------
// Performs prepulse fit on a hit and returns baseline -corrected 
// vector<float> for easier integration or for making ave waveform
//std::vector<float> OpHitBuilderAlg::HitBaselineCorrection( std::vector<short> wfm, short hit, int pre, int post){
//
//}

//-------------------------------------------------------------
// (Overloaded function, version 1/2)
// Scans segment of vector<float> around a designated sample point
// and returns the local minimum
float OpHitBuilderAlg::GetLocalMinimum(std::vector<float> v, short hit)
{
  size_t r1 = 5;
  size_t rmin = 5;
  size_t r2 = 25; 
  float  low                 = 9999;
  size_t x1 = hit-r1;
  size_t x2 = hit+r2;
  
  int counter = 0;

  // Go at LEAST "rmin" samples forward of hit, but no more than
  // r2.  If the value fails to exceed the amplitude for 5 
  // consecutive samples following the first 5 samples, or if
  // we've reached r2, end the scan.
  for( size_t j = x1; j<x2; j++){
    if( v[j] < low ) {
      low = v[j];
      counter=0;
    }
    if( j > hit+rmin && v[j] > low ) counter++;

    if( counter==5 ) break;
  }

  return low;
}


//-------------------------------------------------------------
// (Overloaded function, version 2/2)
// Scans segment of vector<short> around a designated sample point
// and returns the local minimum
short OpHitBuilderAlg::GetLocalMinimum(std::vector<short> v, short hit)
{
  // Convert vector<short> to vector<float> and pass
  // to the previous function
  std::vector<float> vd(v.begin(),v.end());
  return (short)GetLocalMinimum(vd,hit);
}


//--------------------------------------------------------------
// Get only amplitude of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant amplitude information).
float OpHitBuilderAlg::GetHitAmplitude(std::vector<short> wfm, short hit,short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit)[0];
  //std::vector<float> tmp = GetHitInfo(wfm,hit,prev_hit);
  //return tmp[0]; 
}


//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitPromptIntegral(std::vector<short> wfm, short hit, short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit)[1];
  //std::vector<float> tmp = GetHitInfo(wfm,hit,prev_hit);
  //return tmp[1];
}


//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitFullIntegral(std::vector<short> wfm, short hit, short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit)[2];
  //std::vector<float> tmp = GetHitInfo(wfm,hit,prev_hit);
  //return tmp[2];
}


//--------------------------------------------------------------
// Decide if a given OpDetPulse qualifies as a "clean beam waveform."
// ie, (1) Timestamp within beam window, and (2) exactly one optical hit
// that occurs at the trigger time (+/- 1% tolerance)
bool OpHitBuilderAlg::IsCleanBeamWaveform( raw::OpDetPulse opdetpulse )
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


//------------------------------------------------------------
// Single PE finder
std::vector<float> OpHitBuilderAlg::GetSinglePEs( raw::OpDetPulse opdetpulse )
{
  std::vector<float> integrals;
  
  short TriggerTime       = opdetpulse.FirstSample();
  std::vector<short> wfm  = opdetpulse.Waveform();
  int NSamples            = wfm.size(); 
  int t1 = std::max(TriggerTime + fHitTimeCutoffLow,0);
  int t2 = std::min(TriggerTime + fHitTimeCutoffHigh,(int)wfm.size());

  // Find waveform baseline and RMS 
  std::vector<float> tmp = GetBaselineAndRMS( wfm, 0, fBaselineWindowLength );
  float baseline  = tmp[0];
  float rms       = tmp[1]; 

  float integral = 0;
  bool flag = false;
  int windowsize = 0;
  int counter = 0;

  // Baseline subtraction and signal inversion
  std::vector<short> wfm_corrected(NSamples);
  for(int i=0; i<NSamples; i++) wfm_corrected[i] = baseline - wfm[i];

  // Scan waveform and search for hits within specified time window
  for(int i=t1; i<t2; i++){
   
    float y  = wfm_corrected[i];
    float yy = y*fMvPerADC;
    bool  IsOverThresh  = (y  >= fPulseHitRMSThresh*rms);
    bool  IsOverLimit   = (yy >= fPulseHitThreshHigh); 

    if( flag ) {

      counter++;
      windowsize++;
      integral += wfm_corrected[i];
      
      // If we're already integrating a PE window and the pulse
      // again rises above thresh, extend the window by resetting counter
      if( IsOverThresh && !IsOverLimit ){
        counter = 0;
      }
      
      // If pulse extends above upper limit or if window length
      // is exceeded due to multiple merges, abort mission.
      if( IsOverLimit || (windowsize > 2*fSER_PostWindow) ){
        counter = 0;
        integral = 0;
        flag = false;
        windowsize = 0;
        continue;
      }
      
      // If we reached the end of the allotted window, add
      // integral to vector and reset everything
      if( counter == fSER_PostWindow ){
        if(bVerbose) std::cout<<"finished PE window of size "<<windowsize<<": "<<integral<<" ADCs\n";
        integrals.push_back(integral);
        integral = 0;
        windowsize = 0;
        counter = 0;
        flag = false;
        

        continue;
      }
    
    }

    if( !flag && IsOverThresh && !IsOverLimit ){
      // Found a "PE hit", so start integral by
      // adding preceeding prepulse region
      flag = true;
      for(int ii=0; ii<fSER_PreWindow; ii++) integral += wfm_corrected[i-fSER_PreWindow+ii];
      integral += wfm_corrected[i];
    }
  } 

  return integrals;

}

