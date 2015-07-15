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

//ROOT
#include <TH1F.h>



//--------------------------------------------
class OpHitBuilderAlg{
 public:
  
  //Constructor/destructor
  OpHitBuilderAlg( fhicl::ParameterSet const& pset );
  ~OpHitBuilderAlg();

  void reconfigure( fhicl::ParameterSet const& pset );
  
 private:

  // ROOT histograms go here
  // TH1F* ustof_histo;
  // TH1F* timestamp_histo;
  
};


#endif
