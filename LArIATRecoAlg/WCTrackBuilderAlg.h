////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for the wire chamber track      //
// builder algorithm, used to reconstruct momentum and other  //
// geometrical properties of test-beam particles passing      //
// through LArIAT's four wire chambers                        //
//                                                            //
// Authors: Ryan Linehan, rlinehan@stanford.edu               //  
//          Johnny Ho, johnnyho@uchicago.edu                  //
//          Jason St. John, stjohn@fnal.gov                   //
//                                                            //
////////////////////////////////////////////////////////////////


#ifndef WCTRACKBUILDERALG_H
#define WCTRACKBUILDERALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"


//--------------------------------------------
class WCTrackBuilderAlg{
 public:
  
  //Constructor/destructor
  WCTrackBuilderAlg( fhicl::ParameterSet const& pset );
  ~WCTrackBuilderAlg();
  
  
  void reconfigure( fhicl::ParameterSet const& pset );
  
  void firstFunction();
  
  
 private:
  
};


#endif
