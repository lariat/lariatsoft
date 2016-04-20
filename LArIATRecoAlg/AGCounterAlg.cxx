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
}

//------------------------------------------------------------------------------
AGCounterAlg::~AGCounterAlg() {}

//------------------------------------------------------------------------------
void AGCounterAlg::reconfigure( fhicl::ParameterSet const& pset ) {}

//------------------------------------------------------------------------------
void AGCounterAlg::ImportWaveform(std::string                   const& AuxDetName,
                                  std::vector<raw::AuxDetDigit> const& Digits) {
  if (AuxDetName == "USE") {
    fUSEDigits = Digits;
  } else if (AuxDetName == "USW") {
    fUSWDigits = Digits;
  } else if (AuxDetName == "DS1") {
    fDS1Digits = Digits;
  } else if (AuxDetName == "DS2") {
    fDS2Digits = Digits;
  }
  
  return;
}

//------------------------------------------------------------------------------
// Used algorithm here is based on TOF HitsFinder developed by Elena
std::vector<int> AGCounterAlg::HitsFinder(raw::AuxDetDigit WaveformDigit) {
  float Threshold = -40;	///\todo: Check Hardcoded! Hit Threshold
  
  std::vector<int> Hits;
  Hits.clear();
  
  bool RisingEdge = false;
  
  for (size_t i = 2; i < WaveformDigit.NADC()-2; ++i) {
      // Five-point stencil derivative
    float Derivative = float(8*(WaveformDigit.ADC(i+1)) -
                             8*(WaveformDigit.ADC(i-1)) +
                             WaveformDigit.ADC(i-2)     -
                             WaveformDigit.ADC(i+2))/12;
    
    if ((Derivative < Threshold)&&(RisingEdge == false)) {
      Hits.push_back(i);
      RisingEdge = true;
    }
    
    if ((Derivative > std::abs(Threshold))&&(RisingEdge == true)) {
      RisingEdge = false;
    }
  }
  
  // Gives a hit of 0 if nothing else was found to prevent seg faults
  if (Hits.size() == 0) { Hits.push_back(0); }
  
  return Hits;
}


//------------------------------------------------------------------------------
std::vector<std::vector<int> > AGCounterAlg::HitsTimeMatching(std::vector<int> const& USEHits,
                                                              std::vector<int> const& USWHits,
                                                              std::vector<int> const& DS1Hits,
                                                              std::vector<int> const& DS2Hits)
{
  
  std::vector<std::vector<int> > 	AllHitTimePairing;
  std::vector<int>		HitTimePairing;
  for (size_t i = 0; i < USEHits.size(); i++) {
    for (size_t j = 0; j < USWHits.size(); j++) {
      for (size_t k = 0; k < DS1Hits.size(); k++) {
        for (size_t l = 0; l < DS2Hits.size(); l++) {
          HitTimePairing.push_back(USEHits[i]);
          HitTimePairing.push_back(USWHits.at(j));
          HitTimePairing.push_back(DS1Hits.at(k));
          HitTimePairing.push_back(DS2Hits.at(l));
          if (CheckMatched(HitTimePairing)) {
            AllHitTimePairing.push_back(HitTimePairing);
          }
          HitTimePairing.clear();
        }
      }
    }
  }
  
  return AllHitTimePairing;
}

//------------------------------------------------------------------------------
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
  ldp::AGCHits				      SingleHit;
  
  for (size_t i = 0; i < fUSEDigits.size(); ++i) {
    std::vector<int> USEHits = HitsFinder(fUSEDigits[i]);
    std::vector<int> USWHits = HitsFinder(fUSWDigits[i]);
    std::vector<int> DS1Hits = HitsFinder(fDS1Digits[i]);
    std::vector<int> DS2Hits = HitsFinder(fDS2Digits[i]);
    
    std::vector<std::vector<int> > AllHitsTimePairing = HitsTimeMatching(USEHits, USWHits, DS1Hits, DS2Hits);
    
    for (size_t j = 0; j < AllHitsTimePairing.size(); j++) {
      SingleHit.TriggerTimeStamp = fUSEDigits[i].TimeStamp();
      
      SingleHit.HitTimeStampUSE = AllHitsTimePairing.at(j).at(0);
      SingleHit.HitTimeStampUSW = AllHitsTimePairing.at(j).at(1);
      SingleHit.HitTimeStampDS1 = AllHitsTimePairing.at(j).at(2);
      SingleHit.HitTimeStampDS2 = AllHitsTimePairing.at(j).at(3);
						
      SingleHit.HitExistUSE = true;
      SingleHit.HitExistUSW = true;
      SingleHit.HitExistDS1 = true;
      SingleHit.HitExistDS2 = true;
      
      if (AllHitsTimePairing.at(j).at(0) == 0) {
        SingleHit.HitExistUSE = false;
        SingleHit.HitPulseAreaUSE = 0.;
      }
      if (AllHitsTimePairing.at(j).at(1) == 0) {
        SingleHit.HitExistUSW = false;
        SingleHit.HitPulseAreaUSW = 0.;
      }
      if (AllHitsTimePairing.at(j).at(2) == 0) {
        SingleHit.HitExistDS1 = false;
        SingleHit.HitPulseAreaDS1 = 0.;
      }
      if (AllHitsTimePairing.at(j).at(3) == 0) {
        SingleHit.HitExistDS2 = false;
        SingleHit.HitPulseAreaDS2 = 0.;
      }
      
      if (SingleHit.HitExistUSE == true) {
        SingleHit.HitPulseAreaUSE = PulseAreaFinder(fUSEDigits[i], SingleHit.HitTimeStampUSE);
      }
      if (SingleHit.HitExistUSW == true) {
        SingleHit.HitPulseAreaUSW = PulseAreaFinder(fUSWDigits[i], SingleHit.HitTimeStampUSW);
      }
      if (SingleHit.HitExistDS1 == true) {
        SingleHit.HitPulseAreaDS1 = PulseAreaFinder(fDS1Digits[i], SingleHit.HitTimeStampDS1);
      }
      if (SingleHit.HitExistDS2 == true) {
        SingleHit.HitPulseAreaDS2 = PulseAreaFinder(fDS2Digits[i], SingleHit.HitTimeStampDS2);
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
