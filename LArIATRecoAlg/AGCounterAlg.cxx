////////////////////////////////////////////////////////////////
//                                                            //
// Reconstruction Algorithm for the Aerogel Cherenkov Counter //
//                                                            //
// Authors: UT Austin Karol Lang Group			      //
//							      //
//         *Dung Phan (brianp.dung@gmail.com)		      //
//	    Will Flanagan (will.flanagan@utexas.edu)	      //
//	    Brandon Soubasis (brandon.soubasis@gmail.com      //
//                                                            //
////////////////////////////////////////////////////////////////


// ART
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArIAT SOFT
#include "LArIATRecoAlg/AGCounterAlg.h"

// C++ STDLIB
#include <iostream>
#include <cmath>
#include <cstdlib>

//------------------------------------------------------------------------------
AGCounterAlg::AGCounterAlg(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);

  art::ServiceHandle<art::TFileService> tfs;

  fDer1p06_1 = tfs->make<TH1F>("Der1p06_1","Der1p06_1",400,-200.0,200.0);
  fDer1p06_2 = tfs->make<TH1F>("Der1p06_2","Der1p06_2",400,-200.0,200.0);
  fDer1p10_1 = tfs->make<TH1F>("Der1p10_1","Der1p10_1",400,-200.0,200.0);
  fDer1p10_2 = tfs->make<TH1F>("Der1p10_2","Der1p10_2",400,-200.0,200.0);

  fPed1p06_1 = tfs->make<TH1F>("Ped1p06_1","1p06_1 Ped",1000,0.0,2000.0);
  fPed1p06_2 = tfs->make<TH1F>("Ped1p06_2","1p06_2 Ped",1000,0.0,2000.0);
  fPed1p10_1 = tfs->make<TH1F>("Ped1p10_1","1p10_1 Ped",1000,0.0,2000.0);
  fPed1p10_2 = tfs->make<TH1F>("Ped1p10_2","1p10_2 Ped",1000,0.0,2000.0);

}

//------------------------------------------------------------------------------
AGCounterAlg::~AGCounterAlg() {}

//------------------------------------------------------------------------------
void AGCounterAlg::reconfigure( fhicl::ParameterSet const& pset ) {}

//------------------------------------------------------------------------------
void AGCounterAlg::ImportWaveform(std::string                   const& AuxDetName,
                                  std::vector<raw::AuxDetDigit> const& Digits) {

  if (AuxDetName == "AG1p10_1") { AG1p10_1Digits = Digits; }   
  else if (AuxDetName == "AG1p10_2") { AG1p10_2Digits = Digits; }   
  else if (AuxDetName == "AG1p06_1") { AG1p06_1Digits = Digits; }   
  else if (AuxDetName == "AG1p06_2") { AG1p06_2Digits = Digits; }   
  
  return;
}

//------------------------------------------------------------------------------
// Used algorithm here is based on TOF HitsFinder developed by Dan Smith
std::vector<int> AGCounterAlg::HitsFinder(std::string const& AuxDetName, raw::AuxDetDigit WaveformDigit) {

  // This threshold was determined from the gradient histograms
  //   created by this same function call
  float Threshold = -8.25;
  
  std::vector<int> Hits;
  Hits.clear();
  
  bool RisingEdge = false;
  int wait = 0;

  // Requires the waveform is sufficiently long to prevent segfaults
  if(WaveformDigit.NADC() < 200) {
    if (Hits.size() == 0) { Hits.push_back(0); }    
    return Hits;
  }

  // Loop over all the ADC values
  //   within a bounds of 5 to NADC-50 to prevent seg faults
  //   in the differentiation equation
  for (size_t i = 5; (int)i < (int)WaveformDigit.NADC() - 50; ++i) {

    // Numerical differentiation
    double gradient = double(8.0*WaveformDigit.ADC(i+1) - 8.0*WaveformDigit.ADC(i-1) +
                             WaveformDigit.ADC(i-2) - WaveformDigit.ADC(i+2))/12.0;                  
    // Keep track of the pedistals and the gradient
    if(AuxDetName == "AG1p10_1") {
      fPed1p10_1->Fill(WaveformDigit.ADC(i));
      fDer1p10_1->Fill(gradient);
    } else if(AuxDetName == "AG1p10_2") {
      fPed1p10_2->Fill(WaveformDigit.ADC(i));
      fDer1p10_2->Fill(gradient);
    } else if(AuxDetName == "AG1p06_1") {
      fPed1p06_1->Fill(WaveformDigit.ADC(i));
      fDer1p06_2->Fill(gradient);
    } else if(AuxDetName == "AG1p06_2") {
      fPed1p06_2->Fill(WaveformDigit.ADC(i));
      fDer1p06_2->Fill(gradient);
    }
    
    // The hit finding algorithm:
    //
    // If the gradient is below a threshold, raise a flag
    //     to say that the loop is currently 'in' a hit
    // The alg then begins a counter
    // If the gradient rises above the threshold
    //     or if the counter gets too high (30 ticks)
    //     then the hit is recorded
    // If the hit wasn't one tick long, it isn't recorded
    //     to prevent exceptionally high ADC fluctuation
    //

    if(gradient < Threshold && RisingEdge == false) {
      RisingEdge = true;
      wait = 0;
    }

    if(RisingEdge) { wait++; }
    
    if((gradient > Threshold && RisingEdge == true) || wait > 30) {
      if(wait != 0) { Hits.push_back(i-wait); }
      RisingEdge = false;
    }

  }
  
  // Gives a hit of 0 if nothing else was found to prevent seg faults
  if (Hits.size() == 0) { Hits.push_back(0); }
  
  return Hits;
}


//------------------------------------------------------------------------------
std::vector<std::vector<int> > AGCounterAlg::HitsTimeMatching(std::vector<int> const&AG1p10_1Hits,
                                                              std::vector<int> const&AG1p10_2Hits,
                                                              std::vector<int> const&AG1p06_1Hits,
                                                              std::vector<int> const&AG1p06_2Hits)
{
  
  std::vector<std::vector<int> > AllHitTimePairing;
  std::vector<int> HitTimePairing;

  // I believe the two AG detectors should be independent
  //   changed it to be coded that way. Might be wrong.

  // AG1p10 hit pairing
  for (size_t i = 0; i < AG1p10_1Hits.size(); i++) {
    for (size_t j = 0; j < AG1p10_2Hits.size(); j++) {

      HitTimePairing.push_back(AG1p10_1Hits.at(i));
      HitTimePairing.push_back(AG1p10_2Hits.at(j));

      // To stand in for the 1p06 pairings
      // that is now done seperately
      // HitsWrapper should be updated 
      // so this isn't required
      HitTimePairing.push_back(0);
      HitTimePairing.push_back(0);

      if(std::abs(AG1p10_1Hits.at(i) - AG1p10_2Hits.at(j)) < 15.) { 
	AllHitTimePairing.push_back(HitTimePairing);
      }
	  
      HitTimePairing.clear();
    }
  }

  // AG1p06 hit pairing
  for (size_t i = 0; i < AG1p06_1Hits.size(); i++) {
    for (size_t j = 0; j < AG1p06_2Hits.size(); j++) {

      // To stand in for the 1p10 pairings
      // that is now done seperately
      // HitsWrapper should be updated 
      // so this isn't required
      HitTimePairing.push_back(0);
      HitTimePairing.push_back(0);

      HitTimePairing.push_back(AG1p06_1Hits.at(i));
      HitTimePairing.push_back(AG1p06_2Hits.at(j));
	  
      if(std::abs(AG1p06_1Hits.at(i) - AG1p06_2Hits.at(j)) < 15.) { 
	AllHitTimePairing.push_back(HitTimePairing);
      }

      HitTimePairing.clear();
    }
  }
  
  return AllHitTimePairing;

}

//------------------------------------------------------------------------------
// This function isn't used anymore after my change in the hit pairing
// Again, the updated hit finding might be a bad idea
// Worked well for the aerogel study for Run I

bool AGCounterAlg::CheckMatched(std::vector<int> const& HitTimePairing)
{
  short TimeThreshold = 15; ///\todo: Check Hardcoded! Hit Time Seperation
  bool FormedHitGroup;
  std::vector<int> TimeMetrics;
  
  for (size_t i = 0; i < HitTimePairing.size(); i++) {
    if (HitTimePairing[i] != 0) {
      TimeMetrics.push_back(HitTimePairing[i]);
    }
  }
  
  FormedHitGroup = true;
  for (size_t i = 0; i < TimeMetrics.size(); i++) {
    for (size_t j = i+1; j < TimeMetrics.size(); j++) {
      if (std::abs(TimeMetrics[i]-TimeMetrics.at(j)) > TimeThreshold) {
        FormedHitGroup = false;
        j = TimeMetrics.size();
        i = j;
      }
    }
  }
  
  return FormedHitGroup;
}



//------------------------------------------------------------------------------
float AGCounterAlg::PulseAreaFinder(raw::AuxDetDigit  const& Digit,
                                    long unsigned int const& HitTime)
{
  
  float Ped = 0.;
  for (size_t i = 0; i < 100; ++i) {
    Ped = Ped + (float)Digit.ADC(i);
  }
  Ped = Ped/100;
  
  float PulseArea = 0.;
  for (size_t i = 0; i < 60; ++i) {
    PulseArea = PulseArea + (float)Digit.ADC(HitTime + i) - Ped;
  }
  
  return PulseArea;
}

//------------------------------------------------------------------------------
std::vector<std::vector<ldp::AGCHits> > AGCounterAlg::AGCHitsWrapper()
{
  std::vector<ldp::AGCHits> AllHitsGroup;
  ldp::AGCHits SingleHit;  

  for (size_t i = 0; i < AG1p10_1Digits.size(); ++i) {

    std::vector<int> AG1p10_1Hits = HitsFinder("AG1p10_1", AG1p10_1Digits[i]); 
    std::vector<int> AG1p10_2Hits = HitsFinder("AG1p10_2", AG1p10_2Digits[i]); 
    std::vector<int> AG1p06_1Hits = HitsFinder("AG1p06_1", AG1p06_1Digits[i]);
    std::vector<int> AG1p06_2Hits = HitsFinder("AG1p06_2", AG1p06_2Digits[i]);

    std::vector<std::vector<int> > AllHitsTimePairing = HitsTimeMatching(AG1p10_1Hits, AG1p10_2Hits, AG1p06_1Hits, AG1p06_2Hits);

    for (size_t j = 0; j < AllHitsTimePairing.size(); j++) {

      SingleHit.TriggerTimeStamp = AG1p10_1Digits[i].TimeStamp();
      
      SingleHit.HitTimeStamp1p10_1 = AllHitsTimePairing.at(j).at(0);
      SingleHit.HitTimeStamp1p10_2 = AllHitsTimePairing.at(j).at(1);
      SingleHit.HitTimeStamp1p06_1 = AllHitsTimePairing.at(j).at(2);
      SingleHit.HitTimeStamp1p06_2 = AllHitsTimePairing.at(j).at(3);
						
      SingleHit.HitExist1p10_1 = true;
      SingleHit.HitExist1p10_2 = true;
      SingleHit.HitExist1p06_1 = true;
      SingleHit.HitExist1p06_2 = true;     

      if (AllHitsTimePairing.at(j).at(0) == 0) {
        SingleHit.HitExist1p10_1 = false;
        SingleHit.HitPulseArea1p10_1 = 0.;
      }
      if (AllHitsTimePairing.at(j).at(1) == 0) {
        SingleHit.HitExist1p10_2 = false;
        SingleHit.HitPulseArea1p10_2 = 0.;
      }
      if (AllHitsTimePairing.at(j).at(2) == 0) {
        SingleHit.HitExist1p06_1 = false;
        SingleHit.HitPulseArea1p06_1 = 0.;
      }
      if (AllHitsTimePairing.at(j).at(3) == 0) {
        SingleHit.HitExist1p06_2 = false;
        SingleHit.HitPulseArea1p06_2 = 0.;
      }
      
      if (SingleHit.HitExist1p10_1 == true) {
        SingleHit.HitPulseArea1p10_1 = PulseAreaFinder(AG1p10_1Digits[i], SingleHit.HitTimeStamp1p10_1);
      }
      if (SingleHit.HitExist1p10_2 == true) {
        SingleHit.HitPulseArea1p10_2 = PulseAreaFinder(AG1p10_2Digits[i], SingleHit.HitTimeStamp1p10_2);
      }
      if (SingleHit.HitExist1p06_1 == true) {
        SingleHit.HitPulseArea1p06_1 = PulseAreaFinder(AG1p06_1Digits[i], SingleHit.HitTimeStamp1p06_1);
      }
      if (SingleHit.HitExist1p06_2 == true) {
        SingleHit.HitPulseArea1p06_2 = PulseAreaFinder(AG1p06_2Digits[i], SingleHit.HitTimeStamp1p06_2);
      }
      
      AllHitsGroup.push_back(SingleHit);
    }
    
    fAllHitsInEvent.push_back(AllHitsGroup);
  }
  
  return fAllHitsInEvent;

}

//------------------------------------------------------------------------------
void AGCounterAlg::clear_aerogel()
{
  fAllHitsInEvent.clear();
}
