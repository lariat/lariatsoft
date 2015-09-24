////////////////////////////////////////////////////////////////
//
// Class definition for OpHitBuilderAlg, a collection of tools
// to be used in creating OpHit objects from PMT waveforms.
//
// Authors: William Foreman, wforeman@uchicago.edu
//
////////////////////////////////////////////////////////////////


#ifndef OPHITBUILDERALG_H
#define OPHITBUILDERALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/AuxDetDigit.h"

//ROOT
#include <TH1F.h>
#include <TGraph.h>


//--------------------------------------------
class OpHitBuilderAlg{
 public:
  
  //Constructor/destructor
  OpHitBuilderAlg( fhicl::ParameterSet const& pset );
  ~OpHitBuilderAlg();

  void reconfigure( fhicl::ParameterSet const& pset );

  std::vector<short>                GetHits( std::vector<short>&);
  std::vector<Double_t>             MakeGradient( std::vector<short> );
  std::vector<Double_t>             GetBaselineAndRMS( std::vector<short>, short, short);
  Double_t                          GetLocalRMSOfGradient( std::vector<Double_t>, short, short);
  std::vector<Double_t>               IntegrateHit( std::vector<short>, short );
  
 private:
  
  Double_t  fGradientHitThreshold;
  Double_t  fGradientRMSFilterThreshold;
  Double_t  fMinHitSeparation;  
  short     fBaselineWindowLength;
  short     fPrePulseBaselineFit;
  short     fPromptWindowLength;
  double    fMvPerADC;
 
  // ROOT objects
  TF1 *prepulse_exp_fit;
  TGraph *graph;
  
};


#endif
