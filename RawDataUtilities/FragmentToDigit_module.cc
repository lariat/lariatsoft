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
// [ ] Add WUT
// [ ] Add MWPCs
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
#include "LariatFragment.h"
//#include "WUTFragment.h"
#include "CAENFragment.h"
//#include "TDCFragment.h"

#include "RawData/AuxDetDigit.h"

#include "TTree.h"

#include <memory>
#include <vector>
#include <string>

enum {
  V1740_N_CHANNELS = 64,
  V1740_N_SAMPLES = 1536,
  V1751_N_CHANNELS = 8,
  V1751_N_SAMPLES = 1792,
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

  std::string           fRawFragmentLabel;    ///< label for module producing artdaq fragments
  std::string  		    fRawFragmentInstance; ///< instance label for artdaq fragments        
  std::string           fCaenV1740Board8Label;
  std::string           fCaenV1751Board1Label;
  std::string           fCaenV1751Board2Label;

};

//------------------------------------------------------------------------------
FragmentToDigit::FragmentToDigit(fhicl::ParameterSet const & p)
//  : EDProducer(p)
{
  this->reconfigure(p);
  produces< std::vector<raw::AuxDetDigit> >("a");
  produces< std::vector<raw::AuxDetDigit> >("b");
  produces< std::vector<raw::AuxDetDigit> >(fCaenV1740Board8Label);
  produces< std::vector<raw::AuxDetDigit> >(fCaenV1751Board1Label);
  produces< std::vector<raw::AuxDetDigit> >(fCaenV1751Board2Label);
}

//------------------------------------------------------------------------------
void FragmentToDigit::reconfigure(fhicl::ParameterSet const & p)
{
  fRawFragmentLabel    = p.get< std::string >("RawFragmentLabel", "daq");
  fRawFragmentInstance = p.get< std::string >("RawFragmentInstance", "SPILL");
  fCaenV1740Board8Label = p.get< std::string >("CaenV1740Board8Label",
                                               "CaenV1740Board8");
  fCaenV1751Board1Label = p.get< std::string >("CaenV1751Board1Label",
                                               "CaenV1751Board1");
  fCaenV1751Board2Label = p.get< std::string >("CaenV1751Board2Label",
                                               "CaenV1751Board2");
}

//------------------------------------------------------------------------------
void FragmentToDigit::beginJob()
{
  return;
}

//------------------------------------------------------------------------------
void FragmentToDigit::produce(art::Event & evt)
{

  ////////////////////////////////////////////////////////////
  // Begin dummies
  ////////////////////////////////////////////////////////////

  std::unique_ptr< std::vector<raw::AuxDetDigit> > partCol (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > partCol2 (new std::vector<raw::AuxDetDigit>);

  std::vector<short> ADCarray (3,1);
  std::cout<<"ADCarray[1]: "<<ADCarray[1]<<std::endl;
  unsigned short TheChannel = 10;
  std::cout<<"TheChannel: "<<TheChannel<<std::endl;
  std::string TheDetector ("DetectorA");
  std::string TheDetector2 ("DetectorB");
  std::cout<<"TheDetector: "<<TheDetector<<std::endl;

  raw::AuxDetDigit Name;
  raw::AuxDetDigit Name2;
  Name = raw::AuxDetDigit(TheChannel,ADCarray,TheDetector);
  Name2 = raw::AuxDetDigit(TheChannel,ADCarray,TheDetector2);
  std::cout<<"Name.NADC(): "<<Name.NADC()<<std::endl;

  partCol->push_back(Name);
  partCol2->push_back(Name2);

  evt.put(std::move(partCol), "a");
  evt.put(std::move(partCol2), "b");

  ////////////////////////////////////////////////////////////
  // End dummies
  ////////////////////////////////////////////////////////////

  art::Handle< std::vector<artdaq::Fragment> > fragments;
  evt.getByLabel(fRawFragmentLabel, fRawFragmentInstance, fragments);

  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1740Board8Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1751Board1Vec
      (new std::vector<raw::AuxDetDigit>);
  std::unique_ptr< std::vector<raw::AuxDetDigit> > caenV1751Board2Vec
      (new std::vector<raw::AuxDetDigit>);

  if ( !fragments.isValid() )
      throw cet::exception("FragmentToDigit")
      << "artdaq::Fragment handle is not valid, bail";
  if ( fragments->size() != 1 )
      throw cet::exception("FragmentToDigit")
      << "artdaq::Fragment handle contains more than one fragment, bail";

  art::EventNumber_t spillNumber = evt.event();

  // get the fragments we are interested in
  const auto& frag((*fragments)[0]);

  const char * bytePtr = reinterpret_cast<const char *> (&*frag.dataBegin());
  LariatFragment * data = new LariatFragment((char *) bytePtr,
      frag.dataSize() * sizeof(unsigned long long));
  std::cout << "Have data fragment "
            << frag.dataSize() * sizeof(unsigned long long)
            << std::endl;
  data->print();

  std::cout << "Run: " << evt.run() << "; subrun: " << evt.subRun()
            << "; spill: " << spillNumber << std::endl;

  const size_t numberCaenFrags = data->caenFrags.size();
  std::cout << "Found " << numberCaenFrags << " CAEN fragments" << std::endl;

  if (numberCaenFrags > 0) {
    std::cout << "Looking at CAEN fragments..." << std::endl;
  }

  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment & caenFrag = data->caenFrags[i];

    uint32_t boardId = caenFrag.header.boardId;
    uint32_t triggerTimeTag = caenFrag.header.triggerTimeTag;

    //caenFrag.print();
    //std::cout << "CAEN event counter: " << caenFrag.header.eventCounter
    //          << std::endl;

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

  evt.put(std::move(caenV1740Board8Vec), fCaenV1740Board8Label);
  evt.put(std::move(caenV1751Board1Vec), fCaenV1751Board1Label);
  evt.put(std::move(caenV1751Board2Vec), fCaenV1751Board2Label);

  return;  
}

DEFINE_ART_MODULE(FragmentToDigit)
