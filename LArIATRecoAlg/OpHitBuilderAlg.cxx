//////////////////////////////////////////////////////////////////////
//                                                                  
// These are functions used to process optical detector information
// (particulary PMTs) and find/integrate hits.  They've been used
// mostly to analyze events in the stopping/decaying cosmic muon 
// samples, but will be adapted for more general use.
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

  // This specifies whether or not waveforms are to be summed
  // to eventually produce average waveform.  Intended to be
  // changed externally, ie:
  //   fOpHitBuilderAlgInstance.AddHitToAverageWaveform=1
  AddHitToAverageWaveform = 0;
 
  // SER waveform vector to sum waveforms into
  SERWaveform.resize(fSER_PreWindow + fSER_PostWindow);
  SERWaveform_count++;
  
  // Initialize some publically-accessible values
  prepulse_baseline = 0.;
  prepulse_rms = 0.;
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
  fGradHitThresh        = pset.get< float >("GradHitThresh",-10);
  fSignalHitThresh      = pset.get< float >("SignalHitThresh",4);
  fPulseHitThreshLow    = pset.get< float >("PulseHitThreshLow",0.);
  fPulseHitThreshHigh   = pset.get< float >("PulseHitThreshHigh",999);
  fPulseHitRMSThresh    = pset.get< float >("PulseHitRMSThresh",2.5);
  fGradRMSThresh        = pset.get< float >("GradRMSThresh",5); 
  fMinHitSeparation     = pset.get< short >("MinHitSeparation",20);
  fFirstHitSeparation   = pset.get< short >("FirstHitSeparation",250);
  fBaselineWindowLength = pset.get< short >("BaselineWindowLength",1000);
  fPrePulseBaselineFit  = pset.get< short >("PrePulseBaselineFit",500);
  fPrePulseDisplay      = pset.get< short >("PrePulseDisplay",500);
  fPrePulseTau1         = pset.get< float >("PrePulseTau1",1400.);
  fPrePulseTau2         = pset.get< float >("PrePulseTau1",1600.);
  fPromptWindowLength   = pset.get< short >("PromptWindowLength",100);
  fFullWindowLength     = pset.get< short >("FullWindowLength",7000);
  fIntegrationWindows   = pset.get< std::vector<short> >("IntegrationWindows");
  fMvPerADC             = pset.get< float >("MvPerADC",0.2);
  fUsePrepulseFit       = pset.get< bool  >("UsePrepulseFit","true");
  fTimestampCut         = pset.get< float >("TimestampCut",5.25);
  fHitTimeCutoffLow     = pset.get< int   >("HitTimeCutoffLow",-100000);
  fHitTimeCutoffHigh    = pset.get< int   >("HitTimeCutoffHigh",100000);
  fSER_PreWindow        = pset.get< short >("SER_PreWindow",5);
  fSER_PostWindow       = pset.get< short >("SER_PostWindow",25);
  fSER_PrePE_RMS_cut    = pset.get< float >("SER_PrePE_RMS_cut",2.5);
  fSER_Grad_cut         = pset.get< float >("SER_Grad_cut",-2.5);
  fSinglePE             = pset.get< float >("SinglePE",85);
  fHitFindingMode       = pset.get< std::string >("HitFindingMode","grad");
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
std::vector<short> OpHitBuilderAlg::GetHits( raw::OpDetPulse &opdetpulse) 
{
  // Extract relevant information from the OpDetPulse object
  std::vector<short> wfm = opdetpulse.Waveform();
  int TriggerTime = (int)opdetpulse.FirstSample();

  // Hit finding limits (+/- trigger time) set in fcl
  int t1 = std::max(TriggerTime + fHitTimeCutoffLow,0);
  int t2 = std::min(TriggerTime + fHitTimeCutoffHigh,(int)wfm.size());

  // Make vector hits to be filled
  std::vector<short> hits;
  std::vector<short> hits_filtered;

  // Begin hit-finding, discriminating on either the gradient or the signal
  if( fHitFindingMode=="grad" )
  {
    
    LOG_VERBATIM("OpHitBuilder") 
    << "Scanning for hits using gradient-threshold method: \n"
    << "(1st-pass thresh " << fGradHitThresh << " ADC/ns, RMS thresh x"<<fGradRMSThresh<<")";
  
    // Set rising_edge to false before starting hit finding loop
    bool rising_edge = false;

    // Make the gradient
    std::vector<float> g = MakeGradient(wfm);
    
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
 
    LOG_VERBATIM("OpHitBuilder")
    << "Found " << hits.size()
    << " hits in first gradient threshold pass"; 
  
    // So now we have a list of hits, but clusters of these could be
    // fakes due to a "noisy" gradient. Filter out the hits that
    // are not high enough above the gradient's local RMS.
    
    // std::vector<short> hits_filtered;
    
    for( size_t i = 0; i < hits.size(); i++ ) {
  
      short     rms_window_size = 100;
      short     rms_window_start = 0;
      short     rms_window_end  = rms_window_size;
  
      // Find the window limits to use in calculating local RMS.  
      // For all hits after the first one, the window's start point 
      // will be limited to fMinHitSeparation after the first hit point
      rms_window_end    = hits[i]-5;
      rms_window_start  = rms_window_end - rms_window_size;
      if( rms_window_start < 0 ) rms_window_start = 0;
      if( (i >= 1) && (hits[0] + fMinHitSeparation > rms_window_start) ) rms_window_start = hits[0] + fMinHitSeparation;
      rms_window_size = rms_window_end - rms_window_start;
  
      // Find gradient's local RMS for the hit using this window
      std::vector<float> tmp = GetBaselineAndRMS( g, rms_window_start, rms_window_end );
      float g_mean = tmp[0];
      float g_rms  = tmp[1];
  
      // Find amplitude of gradient in neighborhood of this hit
      float g_amp = g_mean - GetLocalMinimum(g,hits[i]);
  
      // Calculate quick amplitude of pulse to use in discrimination
      size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hits[i]-10));
      tmp = GetBaselineAndRMS(wfm,0,baseline_win_size);
      float baseline  = tmp[0];
      float hit_amp   = (baseline-GetLocalMinimum(wfm,hits[i]))*fMvPerADC;
      

      LOG_VERBATIM("OpHitBuilder")
      << "  - " << hits[i] << "(pre-pulse " << rms_window_start << "-" << rms_window_end << ")"
      << " gRMS " << g_rms << " (thresh " << g_rms*fGradRMSThresh
      << ", grad amp " << g_amp 
      << ", pulse amp " << hit_amp;
  
      // First check if hit's gradient and pulse amplitude exceed the set thresholds
      if( (g_amp > g_rms*fGradRMSThresh) && (hit_amp > fPulseHitThreshLow) ){
        hits_filtered.push_back(hits[i]);
        LOG_VERBATIM("OpHitBuilder") << "  --> hit passes.";
      }
    
    } // <-- end loop over hits
    LOG_VERBATIM("OpHitBuilder") << hits_filtered.size() << " hits pass filter.";
  
  } // endif GRAD mode
  
  else if ( fHitFindingMode == "signal" ){
    
    LOG_VERBATIM("OpHitBuilder") 
    << "Scanning for hits using signal/pulse-threshold method: \n"
    << "(threshold " << fSignalHitThresh << ")";

    // Get wfm baseline
    float baseline = GetPedestalAndRMS(wfm,0,fBaselineWindowLength)[0];
    
    // Set rising_edge to false before starting hit finding loop
    bool rising_edge = false;
  
    for(size_t i = 0; i < wfm.size(); i++){
      
      // If signal is under threshold, and we are within the hit finding
      // time limits set in the fhicl, we have a hit candidate
      if ( (baseline - wfm[i])*fMvPerADC <= fSignalHitThresh && rising_edge == false && i > (size_t)t1 && i < (size_t)t2){
        rising_edge = true;
        hits_filtered.insert(hits_filtered.end(),i);
        LOG_VERBATIM("OpHitBuilder")
        << " - " << i;
      }
         
      // When we leave the hit candidate, reset rising_edge to false
      if ( (baseline - wfm[i])*fMvPerADC < fSignalHitThresh && rising_edge == true){
        rising_edge = false;
      }
  
    } // <-- end scan over waveform
    LOG_VERBATIM("OpHitBuilder") << "Found " << hits_filtered.size() << " hits.";

  
  } // Endif SIGNAL mode
  // End hitfinding
  
  // Now merge all remaining hits using shorter spacing
  std::vector<short> hits_merged = HitMerger(hits_filtered,fMinHitSeparation,1);
  LOG_VERBATIM("OpHitBuilder") << "Post-merging: " << hits_merged.size() << " hits.";

  return hits_merged;
  
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
  LOG_VERBATIM("OpHitBuilder")
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
// Get (baseline,rms) of a segment of an std::vector<short>
std::vector<float> OpHitBuilderAlg::GetBaselineAndRMS( std::vector<short> wfm, short x1, short x2 )
{
  // Convert to vector<float> and send into
  // other version of function
  std::vector<float> wfm_float(wfm.begin(),wfm.end());
  return GetBaselineAndRMS(wfm_float,x1,x2);
}

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
// Returns a vector<float> containing (1) the hit's amplitude, and 
// all the specified integrals
std::vector<float> OpHitBuilderAlg::GetHitInfo( std::vector<short> wfm, short hit, short prev_hit, std::vector<short> windows)
{

  LOG_VERBATIM("OpHitBuilder")
  << "GETHITINFO: Processing hit at sample " << hit << "(prev hit @ " << prev_hit << ")";

  // Create vector to be returned (amplitude, integralWindow1, integralWindow2, ...)
  // Initialize dummy values.
  std::vector<float> hit_info(windows.size()+1);
  for(size_t i=0; i<windows.size()+1; i++) hit_info[i] = -9999.;
  
  // If the hit is too early to reliably calculate a baseline, OR if the 
  // previous hit is too close, stop now and return defaults
  if( (hit < 100)||(hit-prev_hit < 200 ) ) {
    LOG_VERBATIM("OpHitBuilder") << "!!! Hit is too early or close to prev hit -- abort.";
    return hit_info;
  }

  // Get baseline
  size_t baseline_win_size = std::min(int(fBaselineWindowLength),int(hit-10));
  std::vector<float> tmp = GetBaselineAndRMS(wfm,0,baseline_win_size);
  float baseline = tmp[0];
  float rms      = tmp[1];

  // Determine bounds for fit and integration 
  short x1,x1b,x2,x3;
  x1  = std::max(int(hit - fPrePulseBaselineFit),prev_hit+100);
  x1b = std::max(int(hit - fPrePulseDisplay),int(x1));
  x2  = hit - 10;
  x3  = std::min(int(hit + fFullWindowLength),int(wfm.size()));
  const int prepulse_bins = int(x2) - int(x1);
  const int total_bins    = int(x3) - int(x1);
  LOG_VERBATIM("OpHitBuilder") 
  << "  x1, x2, x3 = " << x1 << "  " << x2 << "  " << x3;
  
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
  
  LOG_VERBATIM("OpHitBuilder") 
  << "  Prepulse baseline (x1-x2) " << prepulse_baseline << ", rms " << prepulse_rms << "\n"
  << "  Waveform baseline         " << baseline <<", rms " << rms << "\n"
  << "  Difference                " << diff << " (" << diff/rms << " * rms)";


  // Get mean of first few samples of the prepulse region
  int  nn = 10;
  float x1_baseline = GetBaselineAndRMS(wfm,x1,x1+nn)[0];
  float norm_set    = (x1_baseline - baseline);
  if (norm_set > 0.) norm_set = 0.;
  
  // Define exponential function and initialize parameters
  TF1 prepulse_exp_fit("prepulse_exp_fit","[0] + [1]*exp(-(x-[2])/[3])",0.,30000.);
  prepulse_exp_fit.FixParameter(0,baseline);
  prepulse_exp_fit.FixParameter(2,x1+nn/2); 
  prepulse_exp_fit.SetParameter(1,norm_set);
  prepulse_exp_fit.SetParLimits(1,-500.,0.);
  prepulse_exp_fit.SetParameter(3,0.5*(fPrePulseTau1+fPrePulseTau2));
  prepulse_exp_fit.SetParLimits(3,fPrePulseTau1,fPrePulseTau2);

  if(fUsePrepulseFit){ 
      
      // Fit prepulse region with above function
      TGraph graph(prepulse_bins,x,y);
      graph.Fit("prepulse_exp_fit","QN0");

  } else {

    // Prepulse fit turned OFF; just use baseline
    // (zero out exponential component)
    prepulse_exp_fit.FixParameter(1,0.);  
  
  }

  // Save fit info into publicly-accessible members   
  fit_SlowNorm        = prepulse_exp_fit.GetParameter(1);
  fit_SlowTau         = prepulse_exp_fit.GetParameter(3);
  int NDF             = prepulse_exp_fit.GetNDF();
  fit_ReducedChi2     = prepulse_exp_fit.GetChisquare()/float(NDF); 
 
  LOG_VERBATIM("OpHitBuilder")
  << "  Resulting fit parameters \n"
  << "    norm            : " << fit_SlowNorm     << "\n"
  << "    tau             : " << fit_SlowTau      << "\n"
  << "    Chi2/NDF        : " << fit_ReducedChi2  << "\n"  
  << "  Function values (baseline subtracted)\n"
  << "    f(x1)  = " << prepulse_exp_fit.Eval(x1) - baseline    << "\n"
  << "    f(x2)  = " << prepulse_exp_fit.Eval(x2) - baseline    << "\n"
  << "    f(hit) = " << prepulse_exp_fit.Eval(hit) - baseline   << "\n"
  << "    actual hit value = " << wfm[hit] - baseline           << "\n"
  << "    f(x3)  = " << prepulse_exp_fit.Eval(x3) - baseline;
   
  // Save average waveform during integration?
  bool flag_ave = ( (AddHitToAverageWaveform)&&(total_bins>=fPrePulseDisplay+fFullWindowLength));
  if(flag_ave) AverageWaveform_count++;

  // Integrate using the fitted function as running baseline
  float integral  = 0 ;
  float amplitude = 0.;
  int   iWindow   = 0 ;

  LOG_VERBATIM("OpHitBuilder") 
  << "  Starting integration from "<<x2<<" to "<<x3;

  for( int i = 0; i < total_bins; i++){
   
    short xx    = x1 + short(i);
    float yy   = prepulse_exp_fit.Eval(xx) - (float)wfm[xx];  
    
    // Once past x2, start integration
    if( xx >= x2 ) integral += yy;

    // Save integral at appropriate window sizes
    if( xx - hit == windows[iWindow]-1) {

      LOG_VERBATIM("OpHitBuilder")
      << "    window "<< iWindow << " (size " << windows[iWindow] << ")"
      << "  = " << integral << " ADC ";
      
      hit_info[iWindow+1] = integral;
      iWindow++;
    }

    // Scan for amplitude when near the hit
    if ( (abs(hit-xx) < 30) && (yy > amplitude) ) amplitude = yy;

    // Add to average waveform if specified
    if( (flag_ave) && ( xx >= x1b)) AverageWaveform.at(xx-x1b) += yy*fMvPerADC;
  
  }

  hit_info[0] = amplitude*fMvPerADC;
  LOG_VERBATIM("OpHitBuilder") << "  amplitude = " << hit_info[0];

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
  return GetHitInfo(wfm,hit,prev_hit,fIntegrationWindows)[0];
}


//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitPromptIntegral(std::vector<short> wfm, short hit, short prev_hit)
{
  return GetHitInfo(wfm,hit,prev_hit,fIntegrationWindows)[1];
}


//--------------------------------------------------------------
// Get only integral of the specified hit (this function calls the broader
// "GetHitInfo" function, but returns only the relevant integral information).
float OpHitBuilderAlg::GetHitFullIntegral(std::vector<short> wfm, short hit, short prev_hit)
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




//------------------------------------------------------------
// Single PE finder
std::vector<std::pair<float,float>> OpHitBuilderAlg::GetSinglePEs( raw::OpDetPulse &opdetpulse )
{

  LOG_VERBATIM("OpHitBuilder") 
  << "Searching for single PE candidates (RMS thresh x " << fPulseHitRMSThresh << ")";

  std::vector<std::pair<float,float>> out;
  
  short TriggerTime       = opdetpulse.FirstSample();
  std::vector<short> wfm  = opdetpulse.Waveform();
  int NSamples            = wfm.size(); 
  int t1 = std::max(TriggerTime + fHitTimeCutoffLow,0);
  int t2 = std::min(TriggerTime + fHitTimeCutoffHigh,(int)wfm.size());

  // Find waveform baseline and RMS 
  std::vector<float> tmp = GetBaselineAndRMS( wfm, 0, fBaselineWindowLength );
  float baseline  = tmp[0];
  float rms       = tmp[1]; 

  float   integral = 0;
  bool    flag = false;
  int     windowsize = 0;
  int     counter = 0;
  int     flat_samples_count = 0;
  float   prePE_baseline = -99;
  float   prePE_rms = 99;

  // Make gradient
  std::vector<float> g = MakeGradient(wfm);
  float hit_grad = 0;

  // Baseline subtraction and signal inversion
  std::vector<float> wfm_corrected(NSamples);
  for(int i=0; i<NSamples; i++) wfm_corrected[i] = baseline - (float)wfm[i];

  // Make empty vector in which to store waveform for each PE candidate
  // to be reset after each PE (or after added to average waveform)
  int SER_bins = fSER_PreWindow + fSER_PostWindow;
  std::vector<float> tmp_wfm(SER_bins);
  int tmp_wfm_i = 0;

  // ------------------------------------------------------------
  // Scan waveform and search for hits within specified time window
  for(int i=t1; i<t2; i++){
   
    float y  = wfm_corrected[i];
    float yy = y*fMvPerADC;
    bool  IsOverThresh  = (y  >= fPulseHitRMSThresh*rms);
    bool  IsOverLimit   = (yy >= fPulseHitThreshHigh); 
    bool  IsPECandidate = ((IsOverThresh) && (!IsOverLimit) && (g[i] <= fSER_Grad_cut) );

    // If we're already in a PE window, increment the
    // counters and add to integral
    if( flag ) {
     
      LOG_VERBATIM("OpHitBuilder")
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag;
      
      counter++;
      windowsize++;
      float yc = wfm_corrected[i] - prePE_baseline;
      integral += yc;

      if(tmp_wfm_i < SER_bins){
        tmp_wfm[tmp_wfm_i] = yc;
        tmp_wfm_i++;
      }
      
      // If another PE is detected after at least 5 ns, extend the window by resetting counter
      if( counter >=5 && IsPECandidate ){
        LOG_VERBATIM("OpHitBuilder") << "  Secondary hit, extending window";
        counter = 0;
      }
      
      // If pulse extends above upper limit or if window length
      // is exceeded due to multiple merges, abort mission.
      if( IsOverLimit || (windowsize > 2*fSER_PostWindow) ){
        LOG_VERBATIM("OpHitBuilder") << "  abort!";
        counter = 0;
        hit_grad = 0;
        integral = 0;
        flag = false;
        windowsize = 0;
        flat_samples_count = 0;
        tmp_wfm.clear();
        tmp_wfm_i = 0;
        continue;
      }
      
      // If we reached the end of the allotted window, add
      // integral to vector and reset everything
      if( counter == fSER_PostWindow ){
       
        LOG_VERBATIM("OpHitBuilder") 
        << "Finished PE window of size "<<windowsize<<", "
        << integral << " ADCs, g = " << hit_grad;
        out.push_back(std::make_pair(integral,hit_grad));
        
        // Add to average waveform if it looks good
        if( (windowsize+fSER_PreWindow == SER_bins) && fabs(integral - fSinglePE)<=5 ){
          
          LOG_VERBATIM("OpHitBuilder") << "Add to average PE wfm.";
          
          for(int ii=0; ii<SER_bins; ii++){
            SERWaveform.at(ii) += tmp_wfm[ii]*fMvPerADC;
          }
          SERWaveform_count++;
        }

        integral = 0;
        hit_grad = 0;
        windowsize = 0;
        counter = 0;
        flat_samples_count = 0;
        flag = false;
        tmp_wfm.clear();
        tmp_wfm_i=0;
        continue;
      }
    
    } // <-- end if(flag)

    int buffer_min = 100;

    if( !IsOverThresh && flat_samples_count <= 50 ) flat_samples_count++;

    // If we're not yet integrating a PE window, signal is within bounds,
    // and the previous 10 samples were below threshold, then we're in business!
    if( !flag && IsPECandidate ){
      
      
      // Find pre-PE baseline
      prePE_baseline = 99;
      prePE_rms      = 99;
      std::vector<float> tmp = GetPedestalAndRMS(wfm_corrected,i-buffer_min,i);
      prePE_baseline = tmp[0];
      prePE_rms      = tmp[1];
     
      LOG_VERBATIM("OpHitBuilder")
      << "  " << i << "  yy = " << yy << " mV (wfm RMS " << rms*fMvPerADC << " mV), "
      << " thresh " << fPulseHitRMSThresh*rms*fMvPerADC << ", g " << g[i] << ", flag " << flag << "\n"
      << "  Potential PE!  preBS/RMS " << prePE_baseline << ", "<< prePE_rms;
 
      // Require flat pre-PE region
      if( prePE_rms < fSER_PrePE_RMS_cut){ 
        // Found a "PE hit", so start integral by
        // adding preceeding prepulse region
        flag = true;
        hit_grad = g[i];

        // Add up previous "prewindow" samples
        for(int ii=0; ii<=fSER_PreWindow; ii++){
          integral += wfm_corrected[i-fSER_PreWindow+ii] - prePE_baseline;
          tmp_wfm[tmp_wfm_i] = wfm_corrected[i-fSER_PreWindow+ii] - prePE_baseline;
          tmp_wfm_i++;
        }
   
        LOG_VERBATIM("OpHitBuilder") << "  Looks good!  Beginning integration..."; 
      } else {
        LOG_VERBATIM("OpHitBuilder") << "  Doesn't pass preBS cut, moving on..."; 
      }

    } // <-- end if(PE cand)
    
    if( IsOverThresh ) flat_samples_count = 0;

  } // <-- end scan over waveform 

  return out;

}



//--------------------------------------------------------------------------
// Get pedestal.  This algorithm first calls the much-faster GetBaselineAndRMS 
// function in order to set the limits of the histogram to be used in proper 
// pedestal calculation.  Fitting to the histogram is more computationally taxing, 
// so unless the baseline needs to be known very precisely, it is more efficient 
// to use the standard GetBaselineAndRMS().
std::vector<float> OpHitBuilderAlg::GetPedestalAndRMS( std::vector<float> wfm, short x1, short x2)
{
  std::vector<float> out(2);
  std::vector<float> tmp(2); 
  tmp = GetBaselineAndRMS( wfm, x1, x2);
  float baseline  = tmp[0];
  float rms1      = tmp[1]; 
  
  TH1F h("h","h",100,baseline-3.*rms1,baseline+3.*rms1);
  for(int i=x1; i<x2; i++) h.Fill(wfm[i]); 

  TF1 f("f","gaus",baseline-3.*rms1,baseline+3.*rms1);
  f.SetParameter(0,1.);
  f.SetParameter(1,baseline);
  f.SetParameter(2,rms1);
  
  h.Fit("f","QN0");
  float ped = f.GetParameter(1);
  float rms2 = f.GetParameter(2);
  float rchi2 = f.GetChisquare()/(float)(f.GetNDF()-1);

  LOG_VERBATIM("OpHitBuilder")
  << "Calculating pedestal        ... " << ped << " (traditional BS " << baseline << ")\n"
  << "Calculating gaussian spread ... " << rms2 << "(traditional RMS " << rms1 << ")\n"
  << "Reduced Chi2 of fit         ... " << rchi2;


  // If resulting pedestal (from fit) is trustworthy,
  // then return it; otherwise, return traditional baseline/rms
  //if( fabs(ped-baseline) <= 1.*rms1){
  if( rchi2 <= 1.2 ){
    out[0] = ped;
    out[1] = rms2;
  } else {
    out[0] = baseline;
    out[1] = rms1;
  }

  return out;
}

// Get pedestal (for vector<short> input)
std::vector<float> OpHitBuilderAlg::GetPedestalAndRMS( std::vector<short> v, short x1, short x2)
{
  std::vector<float> wfm(v.begin(),v.end());
  return GetPedestalAndRMS(wfm,x1,x2);
}
