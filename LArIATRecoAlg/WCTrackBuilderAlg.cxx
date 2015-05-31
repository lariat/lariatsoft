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

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "WCTrackBuilderAlg.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

//--------------------------------------------------------------
//Constructor
WCTrackBuilderAlg::WCTrackBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
}

//--------------------------------------------------------------  
//Destructor
WCTrackBuilderAlg::~WCTrackBuilderAlg()
{

}

//--------------------------------------------------------------
void WCTrackBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  
}

//--------------------------------------------------------------
void WCTrackBuilderAlg::firstFunction()
{

  std::cout << "First function called." << std::endl;

}
