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

//--------------------------------------------------------------
//Constructor
TriggerFilterAlg::TriggerFilterAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
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
  std::vector<std::string> triggerInputs;
  std::vector<std::string> vetoInputs;
  parseFilterPattern( filterPattern, triggerInputs, vetoInputs );
  bool didTriggerPassFilter = comparePatternToTrigger( triggerInputs, vetoInputs, theTrigger );
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
    std::cout << "*********TRIGGER PATTERNS**********" << std::endl;
    for( size_t iWord = 0; iWord < triggerInputs.size(); ++iWord ){
      std::cout << triggerInputs.at(iWord) << std::endl;
    }
    
    std::cout << "*********VETO PATTERNS*************" << std::endl;
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

  
