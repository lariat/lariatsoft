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
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// artdaq
#include "artdaq-core/Data/Fragments.hh"
#include "artdaq-core/Data/Fragment.hh"

// lardata
#include "RawData/RawDigit.h"

// LArIAT
#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/WUTFragment.h"
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/TDCFragment.h"
#include "LArIATFragments/V1495Fragment.h"
#include "SimpleTypesAndConstants/RawTypes.h"

#include "RawData/RawDigit.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/TriggerData.h"
#include "SummaryData/RunData.h"
#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"




typedef std::vector<TDCFragment::TdcEventData> TDCDataBlock; 

//=====================================================================
class FragmentToDigitAlg{
 public:

  //Constructor/destructor
  FragmentToDigitAlg( fhicl::ParameterSet const& pset );
  ~FragmentToDigitAlg();
   
  void reconfigure( fhicl::ParameterSet const& pset );
  
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

  void InitializeRun( art::RunNumber_t runNumber );
  

 private:

  art::ServiceHandle<art::TFileService>      tfs;                      ///< handle to the TFileService
  std::string                                fRawFragmentLabel;        ///< label for module producing artdaq fragments
  std::string                                fRawFragmentInstance;     ///< instance label for artdaq fragments        
  size_t                                     fMaxNumberFitIterations;  ///< number of fit iterations before stopping
  std::map<uint32_t, std::set<uint32_t> >    fOpticalDetChannels;      ///< key is the board ID, set are channels on that board
  std::map< int, std::vector<CAENFragment> > fTriggerToCAENDataBlocks; ///< map trigger ID to vector of CAEN blocks
  std::map< int, std::vector<TDCDataBlock> > fTriggerToTDCDataBlocks;  ///< map trigger ID to vector of TDC blocks
  std::map<size_t, size_t>                   fTDCToStartWire;          ///< map TDCs to first wire attached to TDC
  std::map<size_t, size_t>                   fTDCToChamber;            ///< map TDCs to the chamber they are attached
  std::vector<std::string>                   fMWPCNames;               ///< vector to hold detector names of the MWPCs
  int                                        fRunNumber;               ///< current run number
  size_t                                     fTriggerDecisionTick;     ///< tick at which to expect the trigger decision
  float                                      fTrigger1740Pedestal;     ///< pedestal value for the 1740 readout of the triggers
  float                                      fTrigger1740Threshold;    ///< 1740 readout must go below the pedestal this much to trigger


};


#endif
