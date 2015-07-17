////////////////////////////////////////////////////////////////
//                                                            //
// This is the class implementation for an algorithm made to  //
// filter through LArIAT triggers and select those that are   //
// desirable based on entring an input pattern string.        //
//                                                            //
// Author: Ryan Linehan, rlinehan@stanford.edu                //
// July 10, 2015                                              //
//                                                            //
////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>
#include <bitset>

//LArSoft/LArIATSoft includes
#include "TriggerFilterAlg.h"
#include "Utilities/DatabaseUtilityT1034.h"

//--------------------------------------------------------------
//Constructor
TriggerFilterAlg::TriggerFilterAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  fNumTrigInputs = 16;
  
}

//--------------------------------------------------------------
//Destructor
TriggerFilterAlg::~TriggerFilterAlg()
{

}

//--------------------------------------------------------------
void TriggerFilterAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  fVerbose          = pset.get<bool   >("Verbosity",         0         );
  return;
}

//--------------------------------------------------------------
bool TriggerFilterAlg::doesTriggerPassFilter(raw::Trigger theTrigger, std::string filterPattern )
{
  //First convert the string format of the filter pattern to a 32-bit value
  //to be compared to the TriggerBits member of the Trigger class

  //Creating the TRIGGER and VETO choice vectors
  std::vector<std::string> filterTriggerInputs;
  std::vector<std::string> filterVetoInputs;
  parseFilterPattern( filterPattern, filterTriggerInputs, filterVetoInputs );

  //Checking to make sure that the filter patterns are all included in current 16 options
  bool isPatternValid = verifyInputPattern( filterTriggerInputs, filterVetoInputs );
  if( isPatternValid == false ){
    std::cout << "You have input a filter pattern specifying a device that is not in the 16 trigger inputs for run "
	      << fRun << ". Please correct. Returning false for filter abort." << std::endl;
    return isPatternValid;
  }

  //If the input patterns are all valid, we finally compare them to the triggerBits
  bool didTriggerPassFilter = comparePatternToTrigger( filterTriggerInputs, filterVetoInputs, theTrigger );
  if( fVerbose ){ 
    if( didTriggerPassFilter ){ std::cout << "Trigger passed filter!" << std::endl; }
    else{ std::cout << "Trigger failed filter." << std::endl; }
  }
  return didTriggerPassFilter;
  
}

//============================================================================================
void TriggerFilterAlg::parseFilterPattern( std::string filterPattern,
					   std::vector<std::string> & triggerInputs,
					   std::vector<std::string> & vetoInputs )
{
  //Break up the string into TRIGGER and VETO choices
  bool readingTriggerInput = false;                                  //Specifically, input that is NOT a veto
  std::string newInput;
  for( size_t iChar = 0; iChar < filterPattern.size() ; ++iChar ){
    
    //If there is a +, then add the last word to the triggerInputs vector
    if( filterPattern.at(iChar) == '+' ){
      if( iChar == 0 ){
	readingTriggerInput = true;
	continue;
      }
      else{
	if( readingTriggerInput == true ) triggerInputs.push_back(newInput);
	if( readingTriggerInput == false ) vetoInputs.push_back(newInput);
	newInput.clear();
	readingTriggerInput = true;
      }
    }

    //If there is a -, then add the last word to the vetoInputs vector 
    else if( filterPattern.at(iChar) == '-' ){
      if( iChar == 0 ){
	readingTriggerInput = false;
	continue;
      }
      else{
	if( readingTriggerInput == true ) triggerInputs.push_back(newInput);
	if( readingTriggerInput == false ) vetoInputs.push_back(newInput);
	newInput.clear();
	readingTriggerInput = false;
      }
    }

    //Special case: the last character in the string
    else if( iChar == filterPattern.size()-1 ){
      newInput.push_back(filterPattern.at(iChar));
      if( readingTriggerInput == true ) triggerInputs.push_back(newInput);
      if( readingTriggerInput == false ) vetoInputs.push_back(newInput);
      newInput.clear();
    }

    //Anything else? Extend the word by this character.
    else newInput.push_back(filterPattern.at(iChar));
  }

  if( fVerbose ){  
    std::cout << "--------->TRIGGER PATTERNS" << std::endl;
    for( size_t iWord = 0; iWord < triggerInputs.size(); ++iWord ){
      std::cout << triggerInputs.at(iWord) << std::endl;
    }
    
    std::cout << "--------->VETO PATTERNS" << std::endl;
    for( size_t iWord = 0; iWord < vetoInputs.size(); ++iWord ){
      std::cout << vetoInputs.at(iWord) << std::endl;
    }
  }
}

//============================================================================================
bool TriggerFilterAlg::comparePatternToTrigger( std::vector<std::string> triggerInputs,
						std::vector<std::string> vetoInputs,
						raw::Trigger theTrigger )
{
  //Initializing what inputs are triggered for this trigger
  std::map<std::string,bool> whatIsTriggered;
  initializeBitsToStrings( whatIsTriggered, theTrigger );
  
  //Loop through triggerInputs and see if corresponding bits are 1.
  //If one of them is not, return false.
  for( size_t iTrig = 0; iTrig < triggerInputs.size() ; ++iTrig )
    if( !whatIsTriggered.at(triggerInputs.at(iTrig)) ) return false;
  
  //Loop through the veto inputs and see if corresponding bits are 0
  //If any of them is not, return false.
  for( size_t iVeto = 0; iVeto < vetoInputs.size() ; ++iVeto )
    if( whatIsTriggered.at(vetoInputs.at(iVeto)) ) return false;

  //Else, return true. The trigger has passed.
  return true;
}

//============================================================================================
//This is the function that gets modified when the database method comes along
void TriggerFilterAlg::initializeBitsToStrings( std::map<std::string,bool> & whatIsTriggered,
						raw::Trigger theTrigger )
{
  //Now using the database-retrieval method for getting trigger input values
  for( size_t iConfig = 0; iConfig < fNumTrigInputs ; ++iConfig ){
    std::string xmlTriggerInput = fTriggerInputConfigValues.at(fTriggerInputConfigParams.at(iConfig));
    whatIsTriggered.emplace(xmlTriggerInput,theTrigger.Triggered(iConfig));
  }
  /*
  whatIsTriggered.emplace("WC1",theTrigger.Triggered(0));
  whatIsTriggered.emplace("WC2",theTrigger.Triggered(1));
  whatIsTriggered.emplace("WC3",theTrigger.Triggered(2));
  whatIsTriggered.emplace("WC4",theTrigger.Triggered(3));
  whatIsTriggered.emplace("BEAMON",theTrigger.Triggered(4));
  whatIsTriggered.emplace("USTOF",theTrigger.Triggered(5));
  whatIsTriggered.emplace("DSTOF",theTrigger.Triggered(6));
  whatIsTriggered.emplace("PUNCH",theTrigger.Triggered(7));
  whatIsTriggered.emplace("HALO",theTrigger.Triggered(8));
  whatIsTriggered.emplace("PULSER",theTrigger.Triggered(9));
  whatIsTriggered.emplace("COSMICON",theTrigger.Triggered(10));
  whatIsTriggered.emplace("COSMIC",theTrigger.Triggered(11));
  whatIsTriggered.emplace("SC1 CFD",theTrigger.Triggered(12));
  whatIsTriggered.emplace("SC2 CFD",theTrigger.Triggered(13));
  whatIsTriggered.emplace("SC3 CFD",theTrigger.Triggered(14));
  whatIsTriggered.emplace("MuRS",theTrigger.Triggered(15));
  */
  
  if( fVerbose ){
    std::cout << "---------> TRIGGER INPUT CONFIG" << std::endl;
    for( size_t iTrig = 0; iTrig < fNumTrigInputs; ++iTrig ){
      std::cout << "Bit: " << iTrig << ", Device: " << fTriggerInputConfigValues.at(fTriggerInputConfigParams.at(iTrig)) << std::endl;
    }
  }


  if( fVerbose ){
    //Temporary testing
    std::bitset<32> bits (theTrigger.TriggerBits());
    std::cout << "Activated trigger inputs: ";
    for( size_t iBit = 0; iBit < bits.size() ; ++iBit ){
      std::cout << bits[iBit];
    }
    std::cout << std::endl;
  }
}

  
//============================================================================================
//This is the function that is called to load correct row of the lariat_xml_database table for a run. This must be called within the beginRun() method of your analysis module
void TriggerFilterAlg::loadXMLDatabaseTable( int run )
{
  fRun = run;

  
  if( fVerbose ){
    std::cout << "**************************************************" << std::endl;
    std::cout << " NEW RUN: " << fRun << std::endl;
    std::cout << "**************************************************" << std::endl;
  }
  


  //Specifying the trigger input config labels
  for( size_t iConfig = 0; iConfig < fNumTrigInputs; ++iConfig ){
    char name[40];
    sprintf(name,"v1495_config_v1495_in%d_name",int(iConfig));
    fTriggerInputConfigParams.push_back(name);
  }
  
  //Loading the corresponding string values into a map with the labels
  fTriggerInputConfigValues = fDatabaseUtility->GetConfigValues(fTriggerInputConfigParams,fRun);
}

//============================================================================================
//This is used to determine if there are user-input trigger patterns 
bool TriggerFilterAlg::verifyInputPattern( std::vector<std::string> filterTriggerInputs,
					   std::vector<std::string> filterVetoInputs )
{
  /*
  for( size_t iV = 0; iV < filterVetoInputs.size() ; ++iV )
    std::cout << "FilterVetoInputsVerify: " << filterVetoInputs.at(iV) << std::endl;
  for( size_t iT = 0; iT < filterTriggerInputs.size(); ++iT )
    std::cout << "FilterTriggerInputsVerify: " << filterTriggerInputs.at(iT) << std::endl;
  */
  
  //Checking to make sure that each user defined filter trigger string matches a database trigger string
  for( size_t iTrig = 0; iTrig < filterTriggerInputs.size(); ++iTrig ){
    bool hasMatch = false;
    std::map< std::string, std::string >::iterator iter;
    for( iter = fTriggerInputConfigValues.begin(); iter != fTriggerInputConfigValues.end(); ++iter ){
      if( iter->second == filterTriggerInputs.at(iTrig) ) hasMatch = true;
    }
    if( hasMatch == false ) return false;
  }

  //Checking to make sure that each user defined filter veto string matches a database veto string
  for( size_t iVeto = 0; iVeto < filterVetoInputs.size(); ++iVeto ){
    bool hasMatch = false;
    std::map< std::string, std::string >::iterator iter;
    for( iter = fTriggerInputConfigValues.begin(); iter != fTriggerInputConfigValues.end(); ++iter ){
      if( iter->second == filterVetoInputs.at(iVeto) ) hasMatch = true;
    }
    if( hasMatch == false ) return false;
  }

  return true;
}
