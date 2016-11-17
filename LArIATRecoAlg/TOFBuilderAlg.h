////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for the time of flight          //
// builder algorithm, used to reconstruct the time of flight  //
// of test-beam particles passing through LArIAT's            //
// upstream and downstream time of flight                     //
//                                                            //
// Authors: Elena Gramellini elena.gramellini@yale.edu        //
//                                                            //
////////////////////////////////////////////////////////////////


#ifndef TOFBUILDERALG_H
#define TOFBUILDERALG_H

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

//ROOT
#include <TH1F.h>



//--------------------------------------------
class TOFBuilderAlg{
 public:
  
  //Constructor/destructor
  TOFBuilderAlg( fhicl::ParameterSet const& pset );
  ~TOFBuilderAlg();

  void reconfigure( fhicl::ParameterSet const& pset );
  
  std::vector<float> find_hits(std::vector<float> wv);
  std::vector<float> find_hits(std::vector<float> wv, std::string tag);
  std::vector<float> match_hits(std::vector<float> hits1, std::vector<float> hits2);
  std::vector<float> match_hits(std::vector<float> hits1, std::vector<float> hits2, std::string tag);

  std::pair <std::vector<float>, std::vector<long> > get_TOF_and_TimeStamp (std::vector<const raw::AuxDetDigit*> ust_wv,
									    std::vector<const raw::AuxDetDigit*> dst_wv);
  std::vector<float> get_tof(std::vector<const raw::AuxDetDigit*> ust_wv,
			     std::vector<const raw::AuxDetDigit*> dst_wv);

  std::vector<long> get_timeStampDst(std::vector<const raw::AuxDetDigit*> ust_wv,
				     std::vector<const raw::AuxDetDigit*> dst_wv);

  void clear_tof_and_timeStampDst();

 private:
  
  float  fLinear;
  float  fMultiple;
  double fHitThreshold;
  double fHitWait;
  double fHitDiffMeanUS;
  double fHitDiffMeanDS;
  double fHitMatchThreshold;   

  TH1F*  fdeltaHit;
  TH1F*  fdeltaHitUS;
  TH1F*  fdeltaHitDS;
  TH1F*  fderHit;
  TH1F*  fLenHit;

  TH1F*  fDerUSA;
  TH1F*  fDerUSB;
  TH1F*  fDerDSA;
  TH1F*  fDerDSB;
  TH1F*  fampHitUSA;
  TH1F*  fampHitUSB;
  TH1F*  fampHitDSA;
  TH1F*  fampHitDSB;

  std::vector<float> TOF;	  
  std::vector<long> Dst_Timestamp;
  
};


#endif
