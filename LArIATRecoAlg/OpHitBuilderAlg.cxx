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
  
  // Function to fit prepulse regions
  prepulse_exp_fit = new TF1("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])",0.,10000.);
  prepulse_exp_fit->SetParLimits(3,0.,1600.);
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
  fMinHitSeparation           = pset.get< Double_t >("MinHitSeparation",100);
  fBaselineWindowLength       = pset.get< short >("BaselineWindowLength",1000);
  fPrePulseBaselineFit        = pset.get< short >("PrePulseBaselineFit",300);
  fPromptWindowLength         = pset.get< short >("PromptWindowLength",100);
  fMvPerADC                   = pset.get< double>("MvPerADC",0.2);
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

  //std::cout<<"   found "<<hits.size()<<" hits on first pass.\n";

  // So now we have a list of hits, but clusters of these could be
  // fakes due to a "noisy" gradient. Filter out the hits that
  // are not high enough above the gradient's local RMS.
  
  std::vector<short> hits_filtered;
  
  for( size_t i = 0; i < hits.size(); i++ ) {

    Double_t  rms = 0;
    short     rms_window_size = 1000;
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

    // Find local RMS for the hit using this window
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

  // find mean over interval
  for ( int i = x1; i < x2; i++ ) mean += wfm[i]/N;

  // now get sum of squares
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
// Prompt pulse integral
std::vector<double> OpHitBuilderAlg::IntegrateHit( std::vector<short> wfm, short hit )
{
  // create the vector to fill in subsequent steps
  std::vector<double> hit_info(2);

  size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hit-10));
  Double_t baseline = 0;
  for( size_t i = 0; i < baseline_win_size; i++){
    baseline += wfm[i];
  }
  baseline = baseline/double(baseline_win_size);

  std::cout << "Baseline = " << baseline << " (win size: " << baseline_win_size << ")\n";
  
  short x1,x2,x3;
  x1 = std::max(int(hit - fPrePulseBaselineFit),0);
  x2 = std::min(int(x1) + int(fPromptWindowLength),int(wfm.size()));
  x3 = std::max(int(hit)-5,0);
  const int prepulse_bins = int(x3) - int(x1);
  const int integral_bins = int(x2) - int(x3);

  int x[prepulse_bins];
  int y[prepulse_bins];
  for( int i = 0; i < prepulse_bins; i++) {
    x[i] = int(x1)+i;
    y[i] = int(wfm[x1+short(i)]);
  }

  graph = new TGraph(prepulse_bins,x,y);

  prepulse_exp_fit->FixParameter(0,baseline);
  prepulse_exp_fit->FixParameter(2,hit);
  prepulse_exp_fit->SetParLimits(3,1400,1600);
  graph->Fit(prepulse_exp_fit,"N0");

  std::cout<<"Hit was fit.  Parameters: \n";
  std::cout<<"    Norm            : "<<prepulse_exp_fit->GetParameter(1)<<"\n";
  std::cout<<"    Tau             : "<<prepulse_exp_fit->GetParameter(3)<<"\n";
  std::cout<<"    value at hit    : "<<prepulse_exp_fit->Eval(int(hit))<<"\n";
  std::cout<<"    actual hit value: "<<wfm[hit]<<"\n";

  double integral = 0;
  double amplitude = 0;
  // now integrate
  for( int i = 0; i < integral_bins; i++){
    short xx    = hit - 5 + i;
    short dInt  = prepulse_exp_fit->Eval(xx) - wfm[xx];  
    integral += dInt;
    if ( dInt > amplitude ) amplitude = dInt;
  }

  hit_info[0] = amplitude*fMvPerADC;
  hit_info[1] = integral;

  return hit_info;

}


