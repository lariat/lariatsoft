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
#include <string>
#include <map>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"


//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//LArSoft includes
#include "RawData/TriggerData.h"

//LArIATSoft includes
#include "Utilities/DatabaseUtilityT1034.h"

//--------------------------------------------
class TriggerFilterAlg{
 public:

  //Constructor/Destructor
  TriggerFilterAlg(fhicl::ParameterSet const& pset);
  ~TriggerFilterAlg();

  void reconfigure( fhicl::ParameterSet const& pset );



  bool doesTriggerPassFilter( raw::Trigger theTrigger, std::string filterPattern );
  void parseANDPatterns( std::string filterPattern,
			 std::vector<std::string> & ANDGroups );
  void parseFilterPattern( std::string filterPattern,
			   std::vector<std::string> & triggerInputs,
			   std::vector<std::string> & vetoInputs );
  bool comparePatternToTrigger( std::vector<std::string> triggerInputs,
				std::vector<std::string> vetoInputs,
				raw::Trigger theTrigger );
  void initializeBitsToStrings( std::map<std::string,bool> & whatIsTriggered,
				raw::Trigger theTrigger );
  void loadXMLDatabaseTable( int run );

  bool verifyInputPattern( std::vector<std::string> filterTriggerInputs,
			   std::vector<std::string> filterVetoInputs );
    

 private:
  bool fVerbose;
  size_t fNumTrigInputs;
  art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;
  int fRun;
  std::vector<std::string> fTriggerInputConfigParams;
  std::map<std::string,std::string> fTriggerInputConfigValues;



};

#endif
