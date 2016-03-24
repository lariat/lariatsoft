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


#ifndef AGCOUNTERALG_H
#define AGCOUNTERALG_H

// C++ STDLIB
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>

// ART
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LAR/LARIAT-SOFT
#include "larcore/Geometry/Geometry.h"
#include "lardata/RawData/AuxDetDigit.h"

#include "LArIATDataProducts/AGCounter.h"

// ROOT
#include <TH1F.h>

/*
struct AGCHits {
	long unsigned int 	TriggerTimeStamp;
	
	long unsigned int 	HitTimeStampUSE;
	long unsigned int 	HitTimeStampUSW;
	long unsigned int 	HitTimeStampDS1;
	long unsigned int 	HitTimeStampDS2;
	
	float		HitPulseAreaUSE;
	float		HitPulseAreaUSW;
	float		HitPulseAreaDS1;
	float		HitPulseAreaDS2;
	
	bool 			HitExistUSE;
	bool 			HitExistUSW;
	bool 			HitExistDS1;
	bool 			HitExistDS2;
};
*/


//--------------------------------------------
class AGCounterAlg{
public:
  
  AGCounterAlg(fhicl::ParameterSet const& pset);
  ~AGCounterAlg();
  
  void reconfigure( fhicl::ParameterSet const& pset );
  
  void 						ImportWaveform(std::string, std::vector<raw::AuxDetDigit>);  
  std::vector<std::vector<AGCHits> >		AGCHitsWrapper();

  void clear_aerogel();

private:
  std::vector<raw::AuxDetDigit> 		fUSEDigits;
  std::vector<raw::AuxDetDigit> 		fUSWDigits;
  std::vector<raw::AuxDetDigit>		 	fDS1Digits;
  std::vector<raw::AuxDetDigit> 		fDS2Digits;
  
  std::vector<std::vector<AGCHits> >		fAllHitsInEvent;
  
  std::vector<int>			HitsFinder(raw::AuxDetDigit);
  std::vector<std::vector<int> >	HitsTimeMatching(std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>);
  float					PulseAreaFinder(raw::AuxDetDigit, long unsigned int);
  bool					CheckMatched(std::vector<int>);
};

#endif
