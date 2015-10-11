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
#include "RawData/OpDetPulse.h"

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

  std::vector<short>              GetHits( std::vector<short>);
  std::vector<short>              HitMerger( std::vector<short>, short, int);
  std::vector<double>             MakeGradient( std::vector<short> );
  std::vector<double>             GetBaselineAndRMS( std::vector<short>, short, short);
  std::vector<double>               GetBaselineAndRMS( std::vector<double>, short, short);
  std::vector<double>             GetHitInfo( std::vector<short>, short);
  double                          GetHitAmplitude( std::vector<short>, short);
  double                          GetHitPromptIntegral( std::vector<short>, short);
  double                          GetHitFullIntegral(   std::vector<short>, short);  
  bool                              IsCleanBeamWaveform( raw::OpDetPulse );
  short                             GetLocalMinimum( std::vector<short>, short);
  double                          GetLocalMinimum( std::vector<double>, short);
  
  // Average waveform vector
  std::vector<double>   AverageWaveform;
  int                   AverageWaveform_count;
  int                   AddHitToAverageWaveform;

  // Fit parameters
  double prepulse_baseline;
  double prepulse_rms;
  double fit_FastNorm;
  double fit_FastTau;
  double fit_SlowNorm;
  double fit_SlowTau;
  double fit_ReducedChi2;

  bool bVerbose;
  bool      fUsePrepulseFit;

 private:
  
  double  fGradientHitThreshold;
  double  fPulseHitThreshold;
  double  fGradientRMSFilterThreshold;
  double  fMinHitSeparation;  
  short     fBaselineWindowLength;
  short     fPrePulseBaselineFit;
  short     fPromptWindowLength;
  short     fFullWindowLength;
  double    fMvPerADC;
  double    fTimestampCut;
  double    fPrePulseTau1;
  double    fPrePulseTau2;
  
};


#endif
