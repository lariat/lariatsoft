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

#include <iostream>
#include <cmath>
#include <cstdlib>

//--------------------------------------------------------------
//Constructor
OpHitBuilderAlg::OpHitBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;

}

//--------------------------------------------------------------  
//Destructor
OpHitBuilderAlg::~OpHitBuilderAlg()
{

}

//--------------------------------------------------------------
void OpHitBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{

}

std::vector<short> OpHitBuilderAlg::GetHits( std::vector<short>& wfm, Double_t g_threshold) 
{
  std::vector<short> hits;
 // const short Nsamples = wfm.size(); 

  // Threshold, min hit separation
  Double_t min_hit_separation = 20;     // ns
  //Double_t g_RMSFilterThreshold = 5; 

  // Make the gradient
  const std::vector<short> g = MakeGradient(wfm);
  std::cout<<"   made the gradient!\n";

  bool rising_edge = false;

  for(size_t i = 0; i < g.size(); i++){
    
    
    // If gradient is under threshold, we have a hit candidate
    if (g[i] <= g_threshold && rising_edge == false){

      rising_edge = true;
      
      // if there was already a hit before this one,
      // ensure we are well-separated from it:
      if (  (hits.size() >= 1) && (i - hits[hits.size()-1] >= min_hit_separation)  ) { 
        hits.insert(hits.end(),i);
      }

      // if this is the first hit of the waveform, 
      // go ahead and add it to the hit list.
      if ( hits.size() == 0 ){
        hits.insert(hits.end(),i);
      }

    }
       
    // When we leave the hit candidate, reset rising_edge to false
    if (g[i] > g_threshold && rising_edge == true){
      rising_edge = false;
    }

  }

  // So now we have a list of hits, but some of these could be
  // fakes due to a "noisy" gradient. Filter out the hits that
  // are not high enough above the gradient's local RMS.
  


  return hits;
  
}

std::vector<short> OpHitBuilderAlg::MakeGradient( std::vector<short> wfm )
{
  std::vector<short> g(wfm.size());
  g[0]=0;
  g[1]=0;
  for(size_t i=2; i<wfm.size(); i++){
    g[i] = (wfm[i] - wfm[i-2])*0.5;
  }
  return g;
}
 
