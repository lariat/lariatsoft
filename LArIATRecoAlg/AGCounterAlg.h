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


//--------------------------------------------
class AGCounterAlg{
public:
  
  AGCounterAlg(fhicl::ParameterSet const& pset);
  ~AGCounterAlg();
  
  void reconfigure( fhicl::ParameterSet const& pset );
  
  void ImportWaveform(std::string                   const& string,
                      std::vector<raw::AuxDetDigit> const& digits);
  std::vector<std::vector<ldp::AGCHits> >	AGCHitsWrapper();

  void clear_aerogel();

private:

  TH1F* fDer1p06_1;
  TH1F* fDer1p06_2;
  TH1F* fDer1p10_1;
  TH1F* fDer1p10_2;

  TH1F* fPed1p06_1;
  TH1F* fPed1p06_2;
  TH1F* fPed1p10_1;
  TH1F* fPed1p10_2;

  std::vector<int> HitsFinder(std::string const& AuxDetName, raw::AuxDetDigit);
  std::vector<std::vector<int> > HitsTimeMatching(std::vector<int> const&,
                                                  std::vector<int> const&,
                                                  std::vector<int> const&,
                                                  std::vector<int> const&);
  float PulseAreaFinder(raw::AuxDetDigit  const& digit,
                        long unsigned int const&);
  bool	CheckMatched(std::vector<int> const&);

  std::vector<raw::AuxDetDigit> AG1p10_1Digits;
  std::vector<raw::AuxDetDigit> AG1p10_2Digits; 
  std::vector<raw::AuxDetDigit> AG1p06_1Digits; 
  std::vector<raw::AuxDetDigit> AG1p06_2Digits; 
  
  std::vector<std::vector<ldp::AGCHits> > fAllHitsInEvent;
  
};

#endif
