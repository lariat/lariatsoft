////////////////////////////////////////////////////////////////
//                                                            //
// This is the class definition for an algorithm made to      //
// filter through LArIAT triggers and select those that are   //
// desirable based on entring an input pattern string.        //
//                                                            //
// Author: Ryan Linehan, rlinehan@stanford.edu                //
// July 10, 2015                                              //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef TRIGGERFILTERALG_H
#define TRIGGERFILTERALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//LArSoft includes
#include "RawData/TriggerData.h"

//--------------------------------------------
class TriggerFilterAlg{
 public:

  //Constructor/Destructor
  TriggerFilterAlg(fhicl::ParameterSet const& pset);
  ~TriggerFilterAlg();

  void reconfigure( fhicl::ParameterSet const& pset );

  bool doesTriggerPassFilter( raw::Trigger theTrigger, std::string filterPattern );
  void parseFilterPattern( std::string filterPattern,
			   std::vector<std::string> & triggerInputs,
			   std::vector<std::string> & vetoInputs );
  bool comparePatternToTrigger( std::vector<std::string> triggerInputs,
				std::vector<std::string> vetoInputs,
				raw::Trigger theTrigger );
  void initializeBitsToStrings( std::map<std::string,bool> & whatIsTriggered,
				raw::Trigger theTrigger );
  


 private:
  bool fVerbose;

};

#endif
