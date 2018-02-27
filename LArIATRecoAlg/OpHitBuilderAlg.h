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
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "lardataobj/RawData/OpDetPulse.h"

//ROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>

//--------------------------------------------
class OpHitBuilderAlg{
 
 public:
  
  //Constructor/destructor
  OpHitBuilderAlg( fhicl::ParameterSet const& pset );
  ~OpHitBuilderAlg();

  void reconfigure( fhicl::ParameterSet const& pset );

  raw::OpDetPulse               GetPulse( const art::Event&, int);
  std::vector<short>            GetHits( raw::OpDetPulse& );
  std::vector<short>            GetHits( std::vector<float>&, size_t );
  bool                          IsCleanBeamWaveform( raw::OpDetPulse& );
 
  float                         CalcTruncatedMean( std::vector<float>&, float);  
  float                         CalcTruncatedMean( std::vector<float>&, float, int);  
  
  float                         MeanInRange( const std::vector<float>&, int, int);
  void                          SmoothOutVector( std::vector<float>&, int);
  void                          RebinVector( std::vector<float>&, int);
  
  std::vector<short>            HitMerger( std::vector<short>, short, int);
  std::vector<float>            MakeGradient( const std::vector<short>& );
  std::vector<float>            MakeGradient( const std::vector<float>& );
  std::vector<float>            GetBaselineAndRMS( const std::vector<short>&, short, short);
  std::vector<float>            GetBaselineAndRMS( const std::vector<float>&, short, short);

  void                          CalcBaselineAndRMS( const std::vector<short>&, short, short);
  void                          CalcBaselineAndRMS( const std::vector<float>&, short, short);
  void                          SubtractBaseline( std::vector<float>& );
  void                          SubtractBaseline( std::vector<float>&, float );

  std::vector<float>            GetHitInfo( const std::vector<short>&, short, short, std::vector<short>);
  std::vector<float>            GetHitInfo( const std::vector<short>&, short, short, std::vector<short>, bool);
  std::vector<float>            GetHitInfo( const std::vector<float>&, short, short, std::vector<short>);
  std::vector<float>            GetHitInfo( const std::vector<float>&, short, short, std::vector<short>, bool);
  float                         GetHitAmplitude( const std::vector<short>&, short, short);
  float                         GetHitAmplitude( const std::vector<float>&, short, short);
  float                         GetHitPromptIntegral( const std::vector<short>&, short, short);
  float                         GetHitFullIntegral( const std::vector<short>&, short, short);  
  short                         GetLocalMinimum(const std::vector<short>&, short);
  float                         GetLocalMinimum(const std::vector<float>&, short);
  short                         GetLocalMaximum(const std::vector<short>&, short);
  float                         GetLocalMaximum(const std::vector<float>&, short);
  std::vector<float>		GetPedestalAndRMS( const std::vector<float>&, short, short);
  std::vector<float>		GetPedestalAndRMS( const std::vector<short>&, short, short);
  void                          MaskedBaselineSubtraction(const std::vector<float>&, std::vector<float>&);
  void				SubtractRunningBaseline(const std::vector<short>&, std::vector<float>&, const size_t, const size_t);
  void				SubtractRunningBaseline(const std::vector<float>&, std::vector<float>&, const size_t, const size_t);
//  void				CorrectWfmOvershoot(const std::vector<float>&, std::vector<float>&, std::vector<short>, std::string);
  void				CorrectWfmOvershoot(std::vector<float>&, std::vector<short>, std::string);
  float                         CorrectWfmTF1(const std::vector<float>&, short, const TF1&, int);
  float                         FindIntersection(float x1, float y1, float x2, float y2, float thresh);
  std::string                   eventType(float);
  bool                          eventTypeFilter(float,std::vector<std::string>);
//  void				SubtractRunningBaseline_v2(const std::vector<short>, std::vector<float>&, const size_t, float );
//  void				SubtractRunningBaseline_v2(const std::vector<float>, std::vector<float>&, const size_t, float );
 
  void                          Reset();

  void                          SetGradHitThresh(float th) { fGradHitThresh = th; } 
  void                          SetTau(float tau) { fTau = tau; } 
  
  float                         GetBaseline() const { return fBaseline; }
  float                         GetRMS() const { return fRMS; }
  
  // Fhicl parameters
  std::string fHitFindingMode;
  std::string fDAQModule;
  std::string fInstanceName;
  bool  fUsePrepulseFit;
  bool  fMakeHistograms; 
  float fSignalHitThresh;
  float fPulseHitThreshLow;
  float fPulseHitThreshHigh;
  float fGradRMSThresh;
  short fMinHitSeparation; 
  short fBaselineWindowLength;
  short fPrePulseBaselineFit;
  short fPrePulseDisplay;
  short fFullWindowLength;
  float fMvPerADC;
  float fPrePulseTau1;
  float fPrePulseTau2;
  int   fHitTimeCutoffLow;
  int   fHitTimeCutoffHigh;
  float fTau;
  std::vector<short> fIntegrationWindows;
  
  // Average waveform vector
  std::vector<float>    AverageWaveform;
  int                   AverageWaveform_count;
  int                   fAddHitsToAverageWaveform;
  int                   AveWfmBins;

  // Fit parameters
  float prepulse_baseline;
  float prepulse_rms;
  short prepulse_x1;
  float fit_ZeroPoint;
  float fit_SlowNorm;
  float fit_SlowTau;
  float fit_ReducedChi2;

  // Pulse width (FWHM)
  float PulseWidth;
  
  float fGradHitThresh;
  float fBaseline;
  float fRMS;
 
  int     fMskBaselineSubtr_waitprd;
  size_t  fMskBaselineSubtr_minseg;
  size_t  fMskBaselineSubtr_range;
  float   fMskBaselineSubtr_P;
  float   fMskBaselineSubtr_grmsfac;
  float   fMskBaselineSubtr_adcthresh;
 
 private:
  
  bool fReturnRawIntegrals;

  TH1D* hGradient;
  
};


#endif
