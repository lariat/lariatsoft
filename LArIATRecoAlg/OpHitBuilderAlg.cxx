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


//--------------------------------------------------------------
//Constructor
OpHitBuilderAlg::OpHitBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;
  
  // Function to fit prepulse regions
  prepulse_exp_fit = new TF1("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])",0.,10000.);
  prepulse_exp_fit->SetParLimits(3,1200,2000);
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
}


//--------------------------------------------------------------------------
//GetHits:  The 'meat & potatoes' of OpHitFinding!  This function takes a 
//          waveform and returns a vector of hit times. Hits are found using 
//          a gradient-threshold method, and each is required to exceed some 
//          multiple of the local RMS of the gradient (calculated in a window 
//          preceding the hit) to reduce the frequency of fake hits.
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
    rms = GetLocalRMS( g, rms_window_start, rms_window_end );

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
// GetLocalRMS
Double_t OpHitBuilderAlg::GetLocalRMS( std::vector<Double_t> wfm, short x1, short x2 )
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

//--------------------------------------------------------------
// CalcHitInfo: returns a vector<vector<double>>, where each sub-vector stores
// (a) the hit amplitude,
// (b) the hit's prompt integral (100ns nominally) 
// UNDER CONSTRUCTION!!!
std::vector<std::vector<double>> OpHitBuilderAlg::CalcHitInfo( std::vector<short> wfm, std::vector<short> hit_times )
{
  // initialize supervector to store all the waveform info
  std::vector<std::vector<double>> hitamps( hit_times.size(), std::vector<double> (2,0));

  // if hit vector is empty, end here
  if ( hit_times.size() == 0 ) {
    return hitamps;
  } else {

      // otherwise, proceed...
      // first find baseline
      size_t baseline_win_size = std::min(1000,int(hit_times[0]));
      Double_t baseline = 0;
      for ( size_t i = 0; i < baseline_win_size; i++ ){
        baseline += wfm[i];
      }
      baseline = baseline/double(baseline_win_size);
      
      std::cout << "Baseline == " << baseline << "  (win size: "<<baseline_win_size<<")\n";

      // loop through each of the hits; for each, define an integration window 
      // including some pre-hit region and 100ns of post-hit
      for ( size_t hit_index = 0; hit_index < hit_times.size(); hit_index++ ){
      
      }

      // fit exponential fall to pre-pulse region and extrapolate
      // forward in the integration window; use this as running baseline
      prepulse_exp_fit->FixParameter(0,baseline);
      


      return hitamps;
  }
  
}


