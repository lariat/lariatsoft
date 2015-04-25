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
std::vector<std::vector<int>> fOpDetChID;
};

//------------------------------------------------------------------------------
FragmentToDigit::FragmentToDigit(fhicl::ParameterSet const & p)
//  : EDProducer(p)
{
  this->reconfigure(p);
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

  //art::EventNumber_t spillNumber = evt.event();

  mf::LogInfo("FragmentToDigit")
      << "evt.run(): " << evt.run()
      << "; evt.subRun(): " << evt.subRun()
      << "; evt.event(): " << evt.event()
      << "; evt.time().timeLow(): " << evt.time().timeLow()
      << "; evt.time().timeHigh(): " << evt.time().timeHigh();

  mf::LogInfo("FragmentToDigit")
      << "runNumber: " << runNumber << "; spillNumber: " << spillNumber
      << "; timeStamp: " << timeStamp;

  const size_t numberCaenFrags = data->caenFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberCaenFrags << " CAEN fragments";

  if (numberCaenFrags > 0) {
    mf::LogInfo("FragmentToDigit") << "Looking at CAEN fragments...";
  }

  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment & caenFrag = data->caenFrags[i];

    uint32_t boardId = caenFrag.header.boardId;
    uint32_t triggerTimeTag = caenFrag.header.triggerTimeTag;

    //caenFrag.print();
    //LOG_DEBUG("FragmentToDigit")
    //    << "CAEN event counter: "
    //    << caenFrag.header.eventCounter;

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

    else if (boardId == 8) {
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

//OpDetPulse modification - channels to be chosen accordingly in fcl
	if(fOpDetChID[boardId].size()!=0){
			if(int(fOpDetChID[boardId].size())>int(j)){
				if(fOpDetChID[boardId][j]==int(j)){
        				OpDetVec1->push_back(
            				raw::OpDetPulse(
                				static_cast <unsigned short> (j),
                				caenFragWaveForm,
                				0,
                				static_cast <unsigned int> (triggerTimeTag)
                			)
				);
			}
		}
            
	}

//OpDetPulse modification - channels to be chosen accordingly in fcl
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

//OpDetPulse modification - channels to be chosen accordingly in fcl
	if(fOpDetChID[boardId].size()!=0){
			if(int(fOpDetChID[boardId].size())>int(j)){
				if(fOpDetChID[boardId][j]==int(j)){
        				OpDetVec1->push_back(
            				raw::OpDetPulse(
                				static_cast <unsigned short> (j),
                				caenFragWaveForm,
                				0,
                				static_cast <unsigned int> (triggerTimeTag)
                			)
				);
			}
		}
            
	}

//OpDetPulse modification - channels to be chosen accordingly in fcl
      }
    }

  }

  const int numberWutFrags = data->wutFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberWutFrags << " WUT fragments";

  if (numberWutFrags > 0) {
    mf::LogInfo("FragmentToDigit") << "Looking at WUT fragments...";
  }

  for (int i = 0; i < numberWutFrags; ++i) {
    WUTFragment & wutFrag = data->wutFrags[i];
    //wutFrag.print();
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

  const int numberTdcFrags = data->tdcFrags.size();
  mf::LogInfo("FragmentToDigit")
      << "Found " << numberTdcFrags << " TDC fragments";

  if (numberTdcFrags > 0) {
    mf::LogInfo("FragmentToDigit") << "Looking at TDC fragments...";
  }

  for (int i = 0; i < numberTdcFrags; ++i) {

    TDCFragment & tdcFrag = data->tdcFrags[i];

    std::vector< std::vector<TDCFragment::TdcEventData> > &
        tdcEvents = tdcFrag.tdcEvents;

    //LOG_DEBUG("FragmentToDigit")
    //    << "tdcEvents.size(): " << tdcEvents.size();
    //tdcFrag.print();

    for (size_t j = 0; j < tdcEvents.size(); ++j) {

      if (tdcFrag.controllerHeader.nTDCs != tdcEvents[j].size()) {
        mf::LogError("FragmentToDigit")
            << "*** Fatal nTDCs mismatch: " << tdcEvents[j].size()
            << " != " << tdcFrag.controllerHeader.nTDCs;
      }

      for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS;
           ++tdc_index) {
        TDCFragment::TdcEventData tdcEventData = tdcEvents[j].at(tdc_index);

        std::vector< std::vector<short> > hitsInChannel;

        hitsInChannel.resize(TDCFragment::N_CHANNELS);
        for (size_t channel = 0; channel < TDCFragment::N_CHANNELS;
             ++channel) {
          hitsInChannel[channel].reserve(TDCFragment::MAX_HITS);
        }

        //uint32_t triggerCounter = tdcEventData.tdcEventHeader.triggerCounter;
        //uint16_t controllerTimeStamp =
        //    tdcEventData.tdcEventHeader.controllerTimeStamp;

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
  }

  evt.put(std::move(caenV1740Board8Vec), fCaenV1740Board8Label);
  evt.put(std::move(caenV1751Board1Vec), fCaenV1751Board1Label);
  evt.put(std::move(caenV1751Board2Vec), fCaenV1751Board2Label);
evt.put(std::move(OpDetVec1), fCaenOpLabel1);
evt.put(std::move(OpDetVec2), fCaenOpLabel2);
  evt.put(std::move(wutVec), fWutLabel);

  for (size_t i = 0; i < TDCFragment::MAX_TDCS; ++i) {
    evt.put(std::move(mwpcTdcVecs[i].get()), mwpcTdcLabels[i]);
  }

  return;  
}

DEFINE_ART_MODULE(FragmentToDigit)
