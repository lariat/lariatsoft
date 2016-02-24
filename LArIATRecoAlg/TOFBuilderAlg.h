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
#include "lardata/RawData/AuxDetDigit.h"

//ROOT
#include <TH1F.h>



//--------------------------------------------
class TOFBuilderAlg{
 public:
  
  //Constructor/destructor
  TOFBuilderAlg( fhicl::ParameterSet const& pset );
  ~TOFBuilderAlg();

  void reconfigure( fhicl::ParameterSet const& pset );
  
  std::vector<short> find_hits(std::vector<short> wv);
  std::vector<short> match_hits(std::vector<short> hits1, std::vector<short> hits2);

  std::pair <std::vector<short>, std::vector<long> > get_TOF_and_TimeStamp (std::vector<const raw::AuxDetDigit*> ust_wv,
									    std::vector<const raw::AuxDetDigit*> dst_wv);
  std::vector<short> get_tof(std::vector<const raw::AuxDetDigit*> ust_wv,
			     std::vector<const raw::AuxDetDigit*> dst_wv);

  std::vector<long> get_timeStampDst(std::vector<const raw::AuxDetDigit*> ust_wv,
				     std::vector<const raw::AuxDetDigit*> dst_wv);

  void clear_tof_and_timeStampDst();

 private:
  
  float  fLinear;
  float  fMultiple;
  double fHitThreshold;
  double fHitWait;
  double fHitMatchThreshold;   

  std::vector<short> TOF;	  
  std::vector<long> Dst_Timestamp;
  
};


#endif
