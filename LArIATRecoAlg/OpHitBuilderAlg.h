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

  std::vector<short>              GetHits( raw::OpDetPulse );
  std::vector<short>              HitMerger( std::vector<short>, short, int);
  std::vector<float>             MakeGradient( std::vector<short> );
  std::vector<float>             GetBaselineAndRMS( std::vector<short>, short, short);
  std::vector<float>               GetBaselineAndRMS( std::vector<float>, short, short);
  std::vector<float>             GetHitInfo( std::vector<short>, short, short);
  float                          GetHitAmplitude( std::vector<short>, short, short);
  float                          GetHitPromptIntegral( std::vector<short>, short, short);
  float                          GetHitFullIntegral(   std::vector<short>, short, short);  
  bool                              IsCleanBeamWaveform( raw::OpDetPulse );
  short                             GetLocalMinimum( std::vector<short>, short);
  float                          GetLocalMinimum( std::vector<float>, short);
  std::vector<float>             GetSinglePEs( raw::OpDetPulse );
  
  // Average waveform vector
  std::vector<float>   AverageWaveform;
  int                   AverageWaveform_count;
  int                   AddHitToAverageWaveform;

  // Fit parameters
  float prepulse_baseline;
  float prepulse_rms;
  float fit_FastNorm;
  float fit_FastTau;
  float fit_SlowNorm;
  float fit_SlowTau;
  float fit_ReducedChi2;

  bool  bVerbose;
  float fPulseHitRMSThresh;
  bool  fUsePrepulseFit;
  float fGradHitThresh;
  float fPulseHitThreshLow;
  float fPulseHitThreshHigh;
  float fGradientRMSFilterThreshold;
  short fMinHitSeparation; 
  short fFirstHitSeparation; 
  short fBaselineWindowLength;
  short fPrePulseBaselineFit;
  short fPrePulseDisplay;
  short fPromptWindowLength;
  short fFullWindowLength;
  float fMvPerADC;
  float fTimestampCut;
  float fPrePulseTau1;
  float fPrePulseTau2;
  int   fHitTimeCutoffLow;
  int   fHitTimeCutoffHigh;
  short fSER_PreWindow;
  short fSER_PostWindow;

 private:
  
  
};


#endif
