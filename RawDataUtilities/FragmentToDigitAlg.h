////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for fragment-to-digit           //
// conversion.                                                //
//                                                            //
// Author: Brian Rebel, brebel@fnal.gov                       //
//                                                            //
//                                                            //
////////////////////////////////////////////////////////////////


#ifndef FRAGMENTTODIGIT_H
#define FRAGMENTTODIGIT_H

//C++ includes
#include <vector>
#include <set>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Provenance/EventID.h"

namespace raw{
  class RawDigit;
  class AuxDetDigit;
  class OpDetPulse;
  class TriggerData;
}

#include "LArIATFragments/TDCFragment.h"
#include "Utilities/DatabaseUtilityT1034.h"

class CAENFragment;

typedef std::vector<TDCFragment::TdcEventData> TDCDataBlock; 

//=====================================================================
class FragmentToDigitAlg{
 public:

  //Constructor/destructor
  FragmentToDigitAlg( fhicl::ParameterSet const& pset );
  ~FragmentToDigitAlg();
   
  void reconfigure( fhicl::ParameterSet const& pset );
  
  
  std::string TimestampToString(std::time_t const& Timestamp); 
  void InitializeRun(const art::Run& run, art::RunNumber_t runNumber);
  

  void makeTheDigits(std::vector<CAENFragment> caenFrags,
                     std::vector<TDCDataBlock> tdcDataBlocks,
                     std::vector<raw::AuxDetDigit> & auxDigits,
                     std::vector<raw::RawDigit> & rawDigits,
                     std::vector<raw::OpDetPulse> & opPulses );
  uint32_t triggerBits(std::vector<CAENFragment> const& caenFrags);
  void makeTPCRawDigits(std::vector<CAENFragment> const& caenFrags,
                        std::vector<raw::RawDigit>     & tpcDigits);
  float findPedestal(const std::vector<short> & adcVec);
  void makeOpDetPulses(std::vector<CAENFragment>    const& caenFrags,
                       std::vector<raw::OpDetPulse>      & opDetPulse);
  void caenFragmentToAuxDetDigits(std::vector<CAENFragment>     const& caenFrags,
                                  std::vector<raw::AuxDetDigit>      & auxDetDigits,
                                  uint32_t                      const& boardId,
                                  std::set<uint32_t>            const& boardChans,
                                  uint32_t                      const& chanOffset,
                                  std::string                   const& detName);


//<<-------------------------------THESE NEED TO BE CHANGED TO ADD CALL TO DATABASE----------------------------------------->>
  void makeMuonRangeDigits(std::vector<CAENFragment>     const& caenFrags,
                           std::vector<raw::AuxDetDigit>      & mrAuxDigits);
  void makeTOFDigits(std::vector<CAENFragment>     const& caenFrags,
                     std::vector<raw::AuxDetDigit>      & tofAuxDigits);
  void makeAeroGelDigits(std::vector<CAENFragment>     const& caenFrags,
                         std::vector<raw::AuxDetDigit>      & agAuxDigits);
  void makeHaloDigits(std::vector<CAENFragment>     const& caenFrags,
                      std::vector<raw::AuxDetDigit>      & hAuxDigits);
  void makeTriggerDigits(std::vector<CAENFragment>     const& caenFrags,
                         std::vector<raw::AuxDetDigit>      & trAuxDigits);
  void InitializeMWPCContainers();
  void CleanUpMWPCContainers();
  void makeMWPCDigits(std::vector<TDCFragment::TdcEventData> const& tdcEventData,
                      std::vector<raw::AuxDetDigit>               & mwpcAuxDigits);

  std::vector<raw::Trigger> makeTheTriggers(art::EventNumber_t                                    const& EventNumber,
                                            std::vector<CAENFragment>                             const& caenFrags,
                                            std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcDataBlocks);

  

 private:

  int                                        		fRunNumber;               	///< current run number
  int                                        		fRun;               	///< current run number
  std::string                                		fRawFragmentLabel;        	///< label for module producing artdaq fragments
  std::string                                		fRawFragmentInstance;     	///< instance label for artdaq fragments        
  size_t                                     		fMaxNumberFitIterations;  	///< number of fit iterations before stopping
  std::map<uint32_t, std::set<uint32_t> >    		fOpticalDetChannels;      	///< key is the board ID, set are channels on that board
  std::map< int, std::vector<CAENFragment> > 		fTriggerToCAENDataBlocks; 	///< map trigger ID to vector of CAEN blocks
  std::map< int, std::vector<TDCDataBlock> > 		fTriggerToTDCDataBlocks;  	///< map trigger ID to vector of TDC blocks
  std::map<size_t, size_t>                   		fTDCToStartWire;          	///< map TDCs to first wire attached to TDC
  std::map<size_t, size_t>                   		fTDCToChamber;            	///< map TDCs to the chamber they are attached
  std::vector<std::string>                   		fMWPCNames;               	///< vector to hold detector names of the MWPCs
  size_t                                     		fTriggerDecisionTick;     	///< tick at which to expect the trigger decision
  float                                      		fTrigger1740Pedestal;     	///< pedestal value for the 1740 readout of the triggers
  float                                      		fTrigger1740Threshold;    	///< 1740 readout must go below the pedestal this much to trigger
  float                                      		fV1751PostPercent;        	///< 1751 PostPercent setting (ranges 0-100)
  art::ServiceHandle<util::DatabaseUtilityT1034> 	fDatabaseUtility;     		///< handle to the DatabaseUtility1034
  std::map< std::string, std::string >       		fConfigValues;            	///< (key, value) pair for the database query result
  std::map< std::string, std::string >       		fHardwareConnections;          	///< (key, value) pair for the hardware database query result----------------jess lines
  std::vector<std::string>                   		fHardwareParams;            	///< vector of parameter names to be queried---------------------------------jess lines
  std::vector<std::string>                   		fConfigParams;            	///< vector of parameter names to be queried
  std::string		                   		fRunDateTime;  	           	///< string of Date/Time to use to querry HardwareConnectionsTable-----------jess lines
  std::uint32_t		                   		fRunTimestamp; 	           	///<Timestamp from runNumber to use to querry HardwareConnectionsTable-------jess lines

};


#endif
