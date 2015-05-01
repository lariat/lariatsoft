////////////////////////////////////////////////////////////////////////
// Class:       FragmentToDigit
// Module Type: producer
// File:        FragmentToDigit_module.cc
//
// Generated at Mon Dec  1 11:28:13 2014 by Will Flanagan using artmod
// from cetpkgsupport v1_07_01.
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// TODO
//////////////////////////////////////////////////////////////
// [x] Add CAEN V1751 board 1
// [x] Add CAEN V1751 board 2
// [x] Add channels 32 to 64 of CAEN V1740 board 8
// [x] Add WUT
// [x] Add MWPCs
// [ ] Add trigger associations (Brian)
// [ ] Improve the matching algorithm (Johnny)
// [ ] Add CAEN V1740 channels for boards 1-7 and channels 1 
//     to 32 of board 8
// [ ] Put Pawel's OpDetPulse modifications back in (Pawel)
// [ ] Add SpillTrailer fragments
// [ ] Add helpful comments throughout code. This may never be
//     checked off.
//////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "artdaq-core/Data/Fragment.hh"

#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/WUTFragment.h"
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/TDCFragment.h"

#include "SimpleTypesAndConstants/RawTypes.h"
#include "RawData/RawDigit.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "TTree.h"

#include <memory>
#include <functional>
#include <vector>
#include <string>

enum {
  V1740_N_CHANNELS = 64,
  V1740_N_SAMPLES = 1536,
  V1751_N_CHANNELS = 8,
  V1751_N_SAMPLES = 1792,
  WUT_N_TDC_CHANNELS = 16,
  WUT_MAX_HITS = 128,
};

class FragmentToDigit;

class FragmentToDigit : public art::EDProducer {
public:
  explicit FragmentToDigit(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FragmentToDigit(FragmentToDigit const &) = delete;
  FragmentToDigit(FragmentToDigit &&) = delete;
  FragmentToDigit & operator = (FragmentToDigit const &) = delete;
  FragmentToDigit & operator = (FragmentToDigit &&) = delete;

  // Required functions.
  void produce(art::Event & evt) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

  void matchFragments(uint32_t & Ntriggers,
                      std::vector<size_t> & v1751InTrigger,
                      std::vector<size_t> & v1740InTrigger,
                      std::vector<size_t> & TDCInTrigger,
                      LariatFragment * data);
  void makeTPCDigits (LariatFragment *data,
		      std::unique_ptr< std::vector<raw::RawDigit> > & tpcDigits);
  void make1751Digits(int i, LariatFragment * data,
                      std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1751Board1Vec,
                      std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1751Board2Vec);
  void make1740Digits(int i, LariatFragment * data,
                      std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1740Board8Vec);
  void makeTDCDigits(int i, LariatFragment * data,
                     std::vector< std::reference_wrapper< std::unique_ptr< std::vector<raw::AuxDetDigit> > > > & mwpcTdcVecs,
                     std::string mwpcTdcLabels[16]);
  void makeWUTDigits(LariatFragment * data,
                     std::unique_ptr< std::vector<raw::AuxDetDigit> > & wutVec);

private:

  std::string fRawFragmentLabel;     ///< label for module producing artdaq fragments
  std::string fRawFragmentInstance;  ///< instance label for artdaq fragments        
  std::string fCaenV1740Board8Label;
  std::string fCaenV1751Board1Label;
  std::string fCaenV1751Board2Label;
  std::string fCaenOpLabel1;
  std::string fCaenOpLabel2;
  std::string fWutLabel;
  std::string fMwpcTdc01Label;
  std::string fMwpcTdc02Label;
  std::string fMwpcTdc03Label;
  std::string fMwpcTdc04Label;
  std::string fMwpcTdc05Label;
  std::string fMwpcTdc06Label;
  std::string fMwpcTdc07Label;
  std::string fMwpcTdc08Label;
  std::string fMwpcTdc09Label;
  std::string fMwpcTdc10Label;
  std::string fMwpcTdc11Label;
  std::string fMwpcTdc12Label;
  std::string fMwpcTdc13Label;
  std::string fMwpcTdc14Label;
  std::string fMwpcTdc15Label;
  std::string fMwpcTdc16Label;

  // variables from the SpillTrailer fragments
  uint32_t runNumber;
  uint32_t spillNumber;
  uint32_t timeStamp;

  uint32_t Ntriggers;
  std::vector<size_t> v1751InTrigger;
  std::vector<size_t> v1740InTrigger;
  std::vector<size_t> TDCInTrigger;

  std::vector<std::vector<int>> fOpDetChID;
};

//------------------------------------------------------------------------------
FragmentToDigit::FragmentToDigit(fhicl::ParameterSet const & p)
//  : EDProducer(p)
{
  this->reconfigure(p);
  
  produces< std::vector<raw::RawDigit> >();

  produces< std::vector<raw::AuxDetDigit> >(fCaenV1740Board8Label);
  produces< std::vector<raw::AuxDetDigit> >(fCaenV1751Board1Label);
  produces< std::vector<raw::AuxDetDigit> >(fCaenV1751Board2Label);
  produces< std::vector<raw::AuxDetDigit> >(fWutLabel);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc01Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc02Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc03Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc04Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc05Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc06Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc07Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc08Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc09Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc10Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc11Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc12Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc13Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc14Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc15Label);
  produces< std::vector<raw::AuxDetDigit> >(fMwpcTdc16Label);
  produces< std::vector<raw::OpDetPulse> >(fCaenOpLabel1);
  produces< std::vector<raw::OpDetPulse> >(fCaenOpLabel2);
}

//------------------------------------------------------------------------------
void FragmentToDigit::reconfigure(fhicl::ParameterSet const & p)
{
  fRawFragmentLabel = p.get< std::string >("RawFragmentLabel", "daq");
  fRawFragmentInstance = p.get< std::string >("RawFragmentInstance", "SPILL");
  fCaenV1740Board8Label = p.get< std::string >("CaenV1740Board8Label",
                                               "CaenV1740Board8");
  fCaenV1751Board1Label = p.get< std::string >("CaenV1751Board1Label",
                                               "CaenV1751Board1");
  fCaenV1751Board2Label = p.get< std::string >("CaenV1751Board2Label",
                                               "CaenV1751Board2");
  fWutLabel = p.get< std::string >("WutLabel", "Wut");
  fMwpcTdc01Label = p.get< std::string >("MwpcTdc01Label", "MwpcTdc01");
  fMwpcTdc02Label = p.get< std::string >("MwpcTdc02Label", "MwpcTdc02");
  fMwpcTdc03Label = p.get< std::string >("MwpcTdc03Label", "MwpcTdc03");
  fMwpcTdc04Label = p.get< std::string >("MwpcTdc04Label", "MwpcTdc04");
  fMwpcTdc05Label = p.get< std::string >("MwpcTdc05Label", "MwpcTdc05");
  fMwpcTdc06Label = p.get< std::string >("MwpcTdc06Label", "MwpcTdc06");
  fMwpcTdc07Label = p.get< std::string >("MwpcTdc07Label", "MwpcTdc07");
  fMwpcTdc08Label = p.get< std::string >("MwpcTdc08Label", "MwpcTdc08");
  fMwpcTdc09Label = p.get< std::string >("MwpcTdc09Label", "MwpcTdc09");
  fMwpcTdc10Label = p.get< std::string >("MwpcTdc10Label", "MwpcTdc10");
  fMwpcTdc11Label = p.get< std::string >("MwpcTdc11Label", "MwpcTdc11");
  fMwpcTdc12Label = p.get< std::string >("MwpcTdc12Label", "MwpcTdc12");
  fMwpcTdc13Label = p.get< std::string >("MwpcTdc13Label", "MwpcTdc13");
  fMwpcTdc14Label = p.get< std::string >("MwpcTdc14Label", "MwpcTdc14");
  fMwpcTdc15Label = p.get< std::string >("MwpcTdc15Label", "MwpcTdc15");
  fMwpcTdc16Label = p.get< std::string >("MwpcTdc16Label", "MwpcTdc16");
  fOpDetChID = p.get< std::vector<std::vector<int>> >("pmt_channel_ids");
  fCaenOpLabel1 = p.get< std::string >("OpDetBoardLabel1", "Caenv1751Optical1");
  fCaenOpLabel2 = p.get< std::string >("OpDetBoardLabel2", "Caenv1751Optical2");
}

//------------------------------------------------------------------------------
void FragmentToDigit::beginJob()
{
  return;
}

//------------------------------------------------------------------------------
void FragmentToDigit::produce(art::Event & evt)
{

  art::Handle< std::vector<artdaq::Fragment> > fragments;
  evt.getByLabel(fRawFragmentLabel, fRawFragmentInstance, fragments);

  std::unique_ptr< std::vector<raw::RawDigit> > tpcDigitVec(new std::vector<raw::RawDigit>);

  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1740Board8Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1751Board1Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1751Board2Vec
      (new std::vector<raw::AuxDetDigit>);

  std::unique_ptr< std::vector<raw::AuxDetDigit> > wutVec
      (new std::vector<raw::AuxDetDigit>);

  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc01Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc02Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc03Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc04Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc05Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc06Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc07Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc08Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc09Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc10Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc11Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc12Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc13Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc14Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc15Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mwpcTdc16Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr<std::vector< raw::OpDetPulse > >  OpDetVec1
      (new std::vector<raw::OpDetPulse>);
  std::unique_ptr<std::vector< raw::OpDetPulse > >  OpDetVec2
      (new std::vector<raw::OpDetPulse>);

  // the following line is way too long
  std::vector< std::reference_wrapper< std::unique_ptr< std::vector<raw::AuxDetDigit> > > > mwpcTdcVecs = {
      mwpcTdc01Vec, mwpcTdc02Vec, mwpcTdc03Vec, mwpcTdc04Vec,
      mwpcTdc05Vec, mwpcTdc06Vec, mwpcTdc07Vec, mwpcTdc08Vec,
      mwpcTdc09Vec, mwpcTdc10Vec, mwpcTdc11Vec, mwpcTdc12Vec,
      mwpcTdc13Vec, mwpcTdc14Vec, mwpcTdc15Vec, mwpcTdc16Vec
      };

  std::string mwpcTdcLabels[16] = {
      fMwpcTdc01Label, fMwpcTdc02Label, fMwpcTdc03Label, fMwpcTdc04Label,
      fMwpcTdc05Label, fMwpcTdc06Label, fMwpcTdc07Label, fMwpcTdc08Label,
      fMwpcTdc09Label, fMwpcTdc10Label, fMwpcTdc11Label, fMwpcTdc12Label,
      fMwpcTdc13Label, fMwpcTdc14Label, fMwpcTdc15Label, fMwpcTdc16Label
      };

  if ( !fragments.isValid() )
      throw cet::exception("FragmentToDigit")
      << "artdaq::Fragment handle is not valid, bail";
  if ( fragments->size() != 1 )
      throw cet::exception("FragmentToDigit")
      << "artdaq::Fragment handle contains more than one fragment, bail";

  // get the fragments we are interested in
  const auto& frag((*fragments)[0]);

  const char * bytePtr = reinterpret_cast<const char *> (&*frag.dataBegin());
  LariatFragment * data = new LariatFragment((char *) bytePtr,
      frag.dataSize() * sizeof(unsigned long long));
  mf::LogInfo("FragmentToDigit")
      << "Have data fragment "
      << frag.dataSize() * sizeof(unsigned long long);
  data->print();
  data->printSpillTrailer();

  LariatFragment::SpillTrailer & spillTrailer = data->spillTrailer;
  runNumber = spillTrailer.runNumber;
  spillNumber = spillTrailer.spillNumber;
  timeStamp = spillTrailer.timeStamp;

  mf::LogInfo("FragmentToDigit")
      << "evt.run(): " << evt.run()
      << "; evt.subRun(): " << evt.subRun()
      << "; evt.event(): " << evt.event()
      << "; evt.time().timeLow(): " << evt.time().timeLow()
      << "; evt.time().timeHigh(): " << evt.time().timeHigh();

  mf::LogInfo("FragmentToDigit")
      << "runNumber: " << runNumber << "; spillNumber: " << spillNumber
      << "; timeStamp: " << timeStamp;

  this->makeTPCDigits(data, tpcDigitVec);

  FragmentToDigit::matchFragments(Ntriggers, v1751InTrigger, v1740InTrigger, TDCInTrigger, data);
  std::cout<<"Ntriggers is: "<<Ntriggers<<std::endl;
  std::cout<<"The size of v1751InTrigger is: "<<v1751InTrigger.size()<<std::endl;
  std::cout<<"The size of v1740InTrigger is: "<<v1740InTrigger.size()<<std::endl;
  std::cout<<"The size of TDCInTrigger is: "<<TDCInTrigger.size()<<std::endl;

  for (size_t i = 0; i < Ntriggers; ++i) {

    std::cout<<"Trigger "<<i<<" has a V1751 Fragment with index "<<v1751InTrigger[i]
                            <<" and a V1740 Fragment with index "<<v1740InTrigger[i]
                            <<", and a TDC Fragment with index "<<TDCInTrigger[i]<<std::endl;

    if(v1751InTrigger[i]){FragmentToDigit::make1751Digits(v1751InTrigger[i], data, caenV1751Board1Vec, caenV1751Board2Vec);}
    if(v1740InTrigger[i]){FragmentToDigit::make1740Digits(v1740InTrigger[i], data, caenV1740Board8Vec);}
    if(TDCInTrigger[i]){FragmentToDigit::makeTDCDigits(TDCInTrigger[i], data, mwpcTdcVecs, mwpcTdcLabels);}

  }

  FragmentToDigit::makeWUTDigits(data, wutVec);

  evt.put(std::move(tpcDigitVec));
  evt.put(std::move(caenV1740Board8Vec), fCaenV1740Board8Label);
  evt.put(std::move(caenV1751Board1Vec), fCaenV1751Board1Label);
  evt.put(std::move(caenV1751Board2Vec), fCaenV1751Board2Label);

  evt.put(std::move(wutVec), fWutLabel);

  for (size_t i = 0; i < TDCFragment::MAX_TDCS; ++i) {
    evt.put(std::move(mwpcTdcVecs[i].get()), mwpcTdcLabels[i]);
  }

  return;  
}

// Matching v1751, v1740, and TDC fragments
// CAEN triggerTimeTag must be multiplied by 0.008 to get microseconds since the start of spill
// TDC tdcTimeStamp must be divided by 106.208 to get microseconds since the start of spill
// I am currently calling triggers within 200 microseconds a match
void FragmentToDigit::matchFragments(uint32_t & Ntriggers,
                                     std::vector<size_t> & v1751InTrigger,
                                     std::vector<size_t> & v1740InTrigger,
                                     std::vector<size_t> & TDCInTrigger,
                                     LariatFragment * data)
{

  size_t v1751FragNumber = 1;
  size_t numberOf1740Matches = 0;
  size_t numberOfTDCMatches = 0;

  std::cout<<"Lets do some fragment timestamp matching so that we can group triggers..."<<std::endl;

  const size_t numberCaenFrags = data->caenFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberCaenFrags << " CAEN fragments";

  const int numberTdcFrags = data->tdcFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberTdcFrags << " TDC fragments";

  if (numberCaenFrags > 0) {
    mf::LogInfo("FragmentToDigit") << "Begin trigger matching against V1751 CAEN fragments...";
  }

  for (size_t i = 0; i < numberCaenFrags; ++i) {

    numberOf1740Matches = 0;
    numberOfTDCMatches = 0;
    CAENFragment & caenFrag = data->caenFrags[i];

    if (caenFrag.header.boardId == 8) {

      v1751InTrigger.push_back(i);
      std::cout<<"v1751 fragment number "<<v1751FragNumber<<" at t="<<caenFrag.header.triggerTimeTag*0.008<<"ns"<<std::endl;

      for (size_t k = 0; k < numberCaenFrags; ++k) {

        CAENFragment & secondCaenFrag = data->caenFrags[k];

        if (secondCaenFrag.header.boardId == 7 && abs(secondCaenFrag.header.triggerTimeTag*0.008 - caenFrag.header.triggerTimeTag*0.008)<200) {

          numberOf1740Matches++;
          if(numberOf1740Matches==1){v1740InTrigger.push_back(k);}

        }

      }

      if(numberOf1740Matches==0){v1740InTrigger.push_back(0);}
      std::cout<<" has "<<numberOf1740Matches<<" matching v1740 fragments "<<std::endl;


      // There is always one TDC fragment, so I am no longer looping over TDC fragments.
      // Different TDC triggers can be found by looping over TDC events

      if( data->tdcFrags.size() < 1) continue;

      TDCFragment & tdcFrag = data->tdcFrags[0];

      for (size_t l = 0; l < tdcFrag.tdcEvents.size(); ++l) {

        if (abs(tdcFrag.tdcEvents[l].at(0).tdcEventHeader.tdcTimeStamp/106.208 - caenFrag.header.triggerTimeTag*0.008)<200) {

          numberOfTDCMatches++;
          if(numberOfTDCMatches==1){TDCInTrigger.push_back(l);}

        }

      }

      if(numberOfTDCMatches==0){TDCInTrigger.push_back(0);}
      std::cout<<" and "<<numberOfTDCMatches<<" matching TDC fragments "<<std::endl;

      v1751FragNumber++;

    }

  }

  Ntriggers=v1751FragNumber-1;

}

void FragmentToDigit::makeTPCDigits(LariatFragment *data,
				    std::unique_ptr< std::vector<raw::RawDigit> > & tpcDigits)
{

  std::vector<CAENFragment> const& caenFrags = data->caenFrags;

  raw::ChannelID_t tpcChan = 0;
  size_t maxChan = 64;

  for(auto const& frag : caenFrags){
    
    // the TPC mapping has the readout going to boards 0-7 of
    // the CAEN 1751, channels 0-63 of the boards 0-6, channels 0-31 of board 7
    if(frag.header.boardId > 7) continue;
    else{
      if(frag.header.boardId < 7) maxChan = 64;
      else maxChan = 32;
      for(size_t chan = 0; chan < maxChan; ++chan){ 
	if(chan > frag.waveForms.size() )
	  throw cet::exception("FragmentToDigit") << "attempting to access channel "
						  << chan << " from 1751 fragment with only "
						  << frag.waveForms.size() << " channels";

	tpcChan = (frag.header.boardId * 64) + chan;
	std::vector<short> const adc(frag.waveForms[chan].data.begin(), frag.waveForms[chan].data.end());
	raw::RawDigit rd(tpcChan, adc.size(), adc);
	tpcDigits->push_back(rd);
      } // end loop to fill channels from this board
    }// end if it is a TPC board      
  }// end loop over caen fragments

  return;
  
}

void FragmentToDigit::make1751Digits(int i, LariatFragment * data,
                                     std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1751Board1Vec,
                                     std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1751Board2Vec)
{

  CAENFragment & caenFrag = data->caenFrags[i];

  uint32_t boardId = caenFrag.header.boardId;
  uint32_t triggerTimeTag = caenFrag.header.triggerTimeTag;

  if (boardId == 8) {
    for (size_t j = 0; j < V1751_N_CHANNELS; ++j) {
      std::vector<short> caenFragWaveForm
          (caenFrag.waveForms[j].data.begin(),
           caenFrag.waveForms[j].data.end());
      caenV1751Board1Vec->push_back(
          raw::AuxDetDigit(
              static_cast <unsigned short> (j),
              caenFragWaveForm,
              fCaenV1751Board1Label,
              static_cast <unsigned long long> (triggerTimeTag)
              )
          );
    }

  }

  else if (boardId == 9) {
    for (size_t j = 0; j < V1751_N_CHANNELS; ++j) {
      std::vector<short> caenFragWaveForm
          (caenFrag.waveForms[j].data.begin(),
           caenFrag.waveForms[j].data.end());
      caenV1751Board2Vec->push_back(
          raw::AuxDetDigit(
              static_cast <unsigned short> (j),
              caenFragWaveForm,
              fCaenV1751Board2Label,
              static_cast <unsigned long long> (triggerTimeTag)
              )
          );
    }

  }

}

void FragmentToDigit::make1740Digits(int i, LariatFragment * data,
                                     std::unique_ptr< std::vector<raw::AuxDetDigit> > & caenV1740Board8Vec)
{

  CAENFragment & caenFrag = data->caenFrags[i];

  uint32_t boardId = caenFrag.header.boardId;
  uint32_t triggerTimeTag = caenFrag.header.triggerTimeTag;

  if (boardId == 7) {
    for (size_t j = 31; j < V1740_N_CHANNELS; ++j) {
      std::vector<short> caenFragWaveForm
          (caenFrag.waveForms[j].data.begin(),
           caenFrag.waveForms[j].data.end());
      caenV1740Board8Vec->push_back(
          raw::AuxDetDigit(
              static_cast <unsigned short> (j),
              caenFragWaveForm,
              fCaenV1740Board8Label,
              static_cast <unsigned long long> (triggerTimeTag)
              )
          );
    }

  }

}

void FragmentToDigit::makeTDCDigits(int i, LariatFragment * data,
                                    std::vector< std::reference_wrapper< std::unique_ptr< std::vector<raw::AuxDetDigit> > > > & mwpcTdcVecs,
                                    std::string mwpcTdcLabels[16])
{

  TDCFragment & tdcFrag = data->tdcFrags[0];

  std::vector< std::vector<TDCFragment::TdcEventData> > &
      tdcEvents = tdcFrag.tdcEvents;

  if (tdcFrag.controllerHeader.nTDCs != tdcEvents[i].size()) {
    mf::LogError("FragmentToDigit")
        << "*** Fatal nTDCs mismatch: " << tdcEvents[i].size()
        << " != " << tdcFrag.controllerHeader.nTDCs;
  }

  for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS;
       ++tdc_index) {
    TDCFragment::TdcEventData tdcEventData = tdcEvents[i].at(tdc_index);

    std::vector< std::vector<short> > hitsInChannel;

    hitsInChannel.resize(TDCFragment::N_CHANNELS);
    for (size_t channel = 0; channel < TDCFragment::N_CHANNELS;
         ++channel) {
      hitsInChannel[channel].reserve(TDCFragment::MAX_HITS);
    }

    uint32_t tdcTimeStamp = tdcEventData.tdcEventHeader.tdcTimeStamp;

    for (size_t hit_index = 0; hit_index < tdcEventData.tdcHits.size();
         ++hit_index) {
      TDCFragment::TdcHit & hit = tdcEventData.tdcHits[hit_index];
      hitsInChannel[size_t (hit.channel)].push_back(short (hit.timeBin));
    }

    for (size_t channel = 0; channel < TDCFragment::N_CHANNELS;
         ++channel) {

      mwpcTdcVecs[tdc_index].get()->push_back(
          raw::AuxDetDigit(
              static_cast <unsigned short> (channel),
              hitsInChannel[channel],
              mwpcTdcLabels[tdc_index],
              static_cast <unsigned long long> (tdcTimeStamp)
              )
          );

    }

  }

}

void FragmentToDigit::makeWUTDigits(LariatFragment * data,
                                    std::unique_ptr< std::vector<raw::AuxDetDigit> > & wutVec)
{

  const int numberWutFrags = data->wutFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberWutFrags << " WUT fragments";

  if (numberWutFrags > 0) {
    mf::LogInfo("FragmentToDigit") << "Looking at WUT fragments...";
  }

  for (int i = 0; i < numberWutFrags; ++i) {
    WUTFragment & wutFrag = data->wutFrags[i];
    uint32_t numberHits = wutFrag.header.nHits;
    std::vector<WUTFragment::WutHit> & hits = wutFrag.hits;

    std::vector< std::vector<short> > hitsInChannel;

    hitsInChannel.resize(WUT_N_TDC_CHANNELS);
    for (size_t channel = 0; channel < WUT_N_TDC_CHANNELS;
         ++channel) {
      hitsInChannel[channel].reserve(WUT_MAX_HITS);
    }

    uint32_t timeHeader = wutFrag.header.timeHeader;

    for (size_t j = 0; j < numberHits; ++j) {
      WUTFragment::WutHit & hit = hits[j];
      hitsInChannel[size_t (hit.channel)].push_back(short (hit.timeBin));
    }

    for (size_t channel = 0; channel < WUT_N_TDC_CHANNELS;
         ++channel) {

      wutVec->push_back(
          raw::AuxDetDigit(
              static_cast <unsigned short> (channel),
              hitsInChannel[channel],
              fWutLabel,
              static_cast <unsigned long long> (timeHeader)
              )
          );

    }
  }

}

DEFINE_ART_MODULE(FragmentToDigit)
