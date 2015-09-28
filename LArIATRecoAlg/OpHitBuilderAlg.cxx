//////////////////////////////////////////////////////////////////////
//                                                                  
// These are functions used for the optical detector hit finder    
// algorithm, which is used to build hits from cryo-PMT waveforms   
// in LArIAT.  Each hit will have a time, amplitude, charge (in     
// integrated ADC), and a flag marking whether or not it is         
// saturated (exceeds the ADC range of the V1751 digitizer).        
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

// LArIAT includes
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
  AverageWaveform.resize(fPrePulseBaselineFit + fFullWindowLength);
  AverageWaveform_count = 0; 
  
  // Initialize some variables
  prepulse_baseline = 0.;
  prepulse_rms      = 0.;
  fit_norm          = 0.;
  fit_tau           = 0.;
  
}

//--------------------------------------------------------------  
//Destructor
OpHitBuilderAlg::~OpHitBuilderAlg()
{
}

//--------------------------------------------------------------
void OpHitBuilderAlg::reconfigure( fhicl::ParameterSet const& pset ){
  fGradientHitThreshold       = pset.get< Double_t >("GradientHitThreshold",-10);
  fGradientRMSFilterThreshold = pset.get< Double_t >("GradientRMSFilterThreshold",5); 
  fMinHitSeparation           = pset.get< Double_t >("MinHitSeparation",20);
  fBaselineWindowLength       = pset.get< short >("BaselineWindowLength",1000);
  fPrePulseBaselineFit        = pset.get< short >("PrePulseBaselineFit",300);
  fPromptWindowLength         = pset.get< short >("PromptWindowLength",100);
  fFullWindowLength           = pset.get< short >("FullWindowLength",7000);
  fMvPerADC                   = pset.get< double>("MvPerADC",0.2);
  fUsePrepulseFit             = pset.get< bool >("UsePrepulseFit","true");
}


//--------------------------------------------------------------------------
//GetHits:  The 'meat & potatoes' of OpHitFinding!  This function takes a 
//          waveform and returns a vector of hit times. Hits are found using 
//          a gradient-threshold method, and each is required to exceed some 
//          multiple of the local RMS of the gradient (calculated in a window 
//          preceding the hit) to reduce the frequency of fake hits.
//          
//          Parameters to set in your FCL if you want to fiddle with things:
//              - GradientHitThreshold (default -10)
//              - GradientRMSFilterThreshold (default 5)
//              - MinHitSeparation (default 100 [ns])
//
std::vector<short> OpHitBuilderAlg::GetHits( std::vector<short>& wfm) 
{
  std::vector<short> hits;

  // Make the gradient
  std::vector<double> g = MakeGradient(wfm);
  
  // Set rising_edge to false before starting hit finding loop
  bool rising_edge = false;

  for(size_t i = 0; i < g.size(); i++){
    
    // If gradient is under threshold, we have a hit candidate
    if (g[i] <= fGradientHitThreshold && rising_edge == false){

      rising_edge = true;
      
      // If there was already a hit before this one,
      // ensure we are well-separated from it:
      if (  (hits.size() >= 1) && (i - hits[hits.size()-1] >= fMinHitSeparation)  ) { 
        hits.insert(hits.end(),i);
      }

      // If first hit of the waveform, go ahead and add 
      // it to the hit list without condition
      if ( hits.size() == 0 ){
        hits.insert(hits.end(),i);
      }

    }
       
    // When we leave the hit candidate, reset rising_edge to false
    if (g[i] > fGradientHitThreshold && rising_edge == true){
      rising_edge = false;
    }

  } // <-- endloop over the gradient


  // So now we have a list of hits, but clusters of these could be
  // fakes due to a "noisy" gradient. Filter out the hits that
  // are not high enough above the gradient's local RMS.
  
  std::vector<short> hits_filtered;
  
  for( size_t i = 0; i < hits.size(); i++ ) {

    Double_t  rms = 0;
    short     rms_window_size = 100;
    short     rms_window_start = 0;
    short     rms_window_end  = rms_window_size;

    // Find the window limits to use in calculating local RMS.  
    // For all hits after the first one, the window's start point 
    // will be limited to the first hit point (for now)
    rms_window_end    = hits[i]-5;
    rms_window_start  = rms_window_end - rms_window_size;
    if( rms_window_start < 0 ) rms_window_start = 0;
    if( i >= 1 ){
      if( hits[0] > rms_window_start ) rms_window_start = hits[0];
    }
    rms_window_size = rms_window_end - rms_window_start;

    // Find gradient's local RMS for the hit using this window
    rms = GetLocalRMSOfGradient( g, rms_window_start, rms_window_end );

    // Find lowest point of gradient in neighborhood of each hit
    Double_t low = 9999.;
    for( size_t j = 0; j<20; j++){
      size_t index = hits[i]-5+j;
      if( g[index] < low ) low = g[index];
    }
    
    //std::cout<<"     "<<hits[i]<<": window: ("<<rms_window_start<<","<<rms_window_end<<")  rms: "<<rms<<"   g: "<<low<<"\n";

    // Only if the hit exceeds the filter threshold,
    // save it into filtered hits
    if( low <= -1.*rms*fGradientRMSFilterThreshold ){
      //std::cout<<"  Hit passes.  rms: "<<rms<<"   low point: "<<low<<"\n";
      hits_filtered.insert(hits_filtered.end(),hits[i]);
    }
  
  } // <-- end loop over hits

  return hits_filtered;
  
}

//--------------------------------------------------------------
// MakeGradient
std::vector<double> OpHitBuilderAlg::MakeGradient( std::vector<short> wfm )
{
  std::vector<double> g(wfm.size());
  g[0]=0;
  g[1]=0;
  for(size_t i=2; i<wfm.size(); i++) g[i] = double(wfm[i] - wfm[i-2])*0.5;
  return g;
}

//--------------------------------------------------------------
// Get (baseline,rms) of a segment of an std::vector<short>
std::vector<double> OpHitBuilderAlg::GetBaselineAndRMS( std::vector<short> wfm, short x1, short x2 )
{
  double mean = 0;
  double sumSquares = 0;
  int N = x2 - x1;

  for ( int i = x1; i < x2; i++ ) mean += double(wfm[i])/N;
  for ( int i = x1; i < x2; i++ ) sumSquares += pow(wfm[i]-mean,2);
  
  std::vector<double> out(2);
  out[0] = mean;
  out[1] = sqrt(sumSquares/N);
  
  return out;
}


//--------------------------------------------------------------
// GetLocalRMS of a segment of std::vector<double>
Double_t OpHitBuilderAlg::GetLocalRMSOfGradient( std::vector<Double_t> wfm, short x1, short x2 )
{
  Double_t mean = 0;
  Double_t sumSquares = 0;
  int N = x2 - x1;

  // find mean over interval
  for ( int i = x1; i < x2; i++ ) mean += wfm[i]/N;

  // now get sum of squares
  for ( int i = x1; i < x2; i++ ) sumSquares += pow(wfm[i]-mean,2);
  
  // return the RMS
  return sqrt(sumSquares/N);
}

// ----------------------
// Prompt pulse amplitude and integral
std::vector<double> OpHitBuilderAlg::IntegrateHit( std::vector<short> wfm, short hit, short n)
{
  // Create vector to be returned 
  std::vector<double> hit_info {0.,0.};
  short xx;
  double yy;
  
  // If the hit is too early to reliably calculate 
  // a baseline, stop now and return zero
  if( hit < 100 ) return hit_info;

  // Get baseline
  size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hit-10));
  std::vector<double> tmp = GetBaselineAndRMS(wfm,0,baseline_win_size);
  double baseline = tmp[0];
  double rms      = tmp[1];

  //std::cout<<"Baseline: "<<baseline<<"\n";
 
  // Determine bounds for fit and integration 
  short x1,x2,x3;
  x2 = hit - 5;
  x1 = std::max(int(hit - fPrePulseBaselineFit),0);
  x3 = std::min(int(hit + n),int(wfm.size()));
  const int prepulse_bins = int(x2) - int(x1);
  const int integral_bins = int(x3) - int(x2);
  const int total_bins    = int(x3) - int(x1);
  //std::cout<<"Integration limits: x1,x2,x3 "<<x1<<"   "<<x2<<"  "<<x3<<std::endl;

  // Fill x,y arrays
  int x[prepulse_bins];
  int y[prepulse_bins];
  for( int i = 0; i < prepulse_bins; i++) {
    x[i] = int(x1)+i;
    y[i] = int(wfm[x1+short(i)]);
  }

  // Find mean of first 20ns -- use this as a limit
  // later in the fit
  double mean = 0;
  for( int i = 0; i < 20; i++){
    mean += (y[i] - baseline)/20.;
  }

  // Find overall prepulse baseline/RMS
  std::vector<double> prepulse_info = GetBaselineAndRMS(wfm,x1,x2);
  prepulse_baseline  = prepulse_info[0];
  prepulse_rms       = prepulse_info[1];
  double diff        = prepulse_baseline - baseline;

  std::cout<<"Prepulse baseline "<<prepulse_baseline<<"   rms "<<prepulse_rms<<std::endl;
  std::cout<<"Waveform baseline "<<baseline<<"   rms "<<rms<<std::endl;
  std::cout<<"Difference        "<<diff<<"   ("<<diff/rms<<" * rms)"<<std::endl;

  // Function for exponential pre-pulse baseline fit
  TF1 prepulse_exp_fit("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])",0.,10000.);

  if(fUsePrepulseFit){ 
    // If the prepulse baseline is consistent with
    // overall waveform baseline, don't do fit
    if( fabs(prepulse_baseline-baseline) < 0.5*rms ){
      std::cout<<"Prepulse region consistent with baseline -- no fit used.\n";
      prepulse_exp_fit.FixParameter(0,baseline);
      prepulse_exp_fit.FixParameter(1,0.);
      prepulse_exp_fit.FixParameter(2,x2);
      prepulse_exp_fit.FixParameter(3,1600.);
    } else {
      // Otherwise, fit it
      TGraph graph(prepulse_bins,x,y);
      prepulse_exp_fit.FixParameter(0,baseline);
      prepulse_exp_fit.SetParameter(1,wfm[x2]-baseline);
      prepulse_exp_fit.FixParameter(2,x2);
      prepulse_exp_fit.FixParameter(3,1600);
      graph.Fit("prepulse_exp_fit","QN0");
    }
  } else {
    // Prepulse fit turned OFF; just use baseline
    prepulse_exp_fit.FixParameter(0,baseline);
    prepulse_exp_fit.FixParameter(1,0.);
    prepulse_exp_fit.FixParameter(2,x2);
    prepulse_exp_fit.FixParameter(3,1600.);
  }

  // Save fit info into publicly-accessible members   
  fit_norm  = prepulse_exp_fit.GetParameter(1);
  fit_tau   = prepulse_exp_fit.GetParameter(3);
   
  std::cout<<"Parameters of fitted prepulse region: \n";
  std::cout<<"    Norm                : "<<fit_norm<<"\n";
  std::cout<<"    Tau                 : "<<fit_tau<<"\n";
  std::cout<<"    fit baseline at hit : "<<prepulse_exp_fit.Eval(double(hit))<<"\n";
  std::cout<<"    waveform baseline   : "<<baseline<<"\n";
  std::cout<<"    actual hit value    : "<<wfm[hit]<<"\n";
  std::cout<<"    mean first 20ns     : "<<mean<<"\n";
  std::cout<<"    f(x1)               : "<<prepulse_exp_fit.Eval(x1)<<"\n";


  // Add to average waveform vector
  std::vector<double> wfm_corrected(total_bins);
  if( n==fFullWindowLength ){

    // Subtract fitted baseline
    for( int i = 0; i < total_bins; i++){
      xx  = x1 + short(i);
      yy  = prepulse_exp_fit.Eval((double)xx) - (double)wfm[xx];
      AverageWaveform.at(i) += yy*double(fMvPerADC);
    }
    AverageWaveform_count++;
  }
  
  
  // Integrate using the fitted function as running baseline
  double integral = 0.;
  double amplitude = 0.;
  for( int i = 0; i < integral_bins; i++){
    xx    = x2 + short(i);
    yy    = prepulse_exp_fit.Eval((double)xx) - (double)wfm[xx];  
    integral += yy;
    if ( yy > amplitude ) amplitude = yy;
  }

  hit_info[0] = amplitude*fMvPerADC;
  hit_info[1] = integral;

  return hit_info;

}

// --------------------------------------------
// Return vector of single p.e. window integrals
std::vector<double> OpHitBuilderAlg::GetSinglePEs( std::vector<short> wfm )
{  
  std::vector<double> integrals;

  return integrals;
}
