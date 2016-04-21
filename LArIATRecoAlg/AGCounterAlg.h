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
  std::vector<int> HitsFinder(raw::AuxDetDigit);
  std::vector<std::vector<int> > HitsTimeMatching(std::vector<int> const&,
                                                  std::vector<int> const&,
                                                  std::vector<int> const&,
                                                  std::vector<int> const&);
  float PulseAreaFinder(raw::AuxDetDigit  const& digit,
                        long unsigned int const&);
  bool	CheckMatched(std::vector<int> const&);

  std::vector<raw::AuxDetDigit> 		fUSEDigits;
  std::vector<raw::AuxDetDigit> 		fUSWDigits;
  std::vector<raw::AuxDetDigit>		 	fDS1Digits;
  std::vector<raw::AuxDetDigit> 		fDS2Digits;
  
  std::vector<std::vector<ldp::AGCHits> > fAllHitsInEvent;
  
};

#endif
