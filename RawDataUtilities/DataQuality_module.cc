////////////////////////////////////////////////////////////////////////
// Class:       DataQuality
// Module Type: analyzer
// File:        DataQuality_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "WUTFragment.h"
#include "CAENFragment.h"
#include "TDCFragment.h"

#include "TTree.h"

#include <vector>
#include <string>

enum {
  V1740_N_CHANNELS = 64,
  V1740_N_SAMPLES = 1536,
  V1751_N_CHANNELS = 8,
  V1751_N_SAMPLES = 1792,
  WUT_MAX_HITS = 128,
};

class DataQuality;

class DataQuality : public art::EDAnalyzer {
public:
  explicit DataQuality(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DataQuality(DataQuality const &) = delete;
  DataQuality(DataQuality &&) = delete;
  DataQuality & operator = (DataQuality const &) = delete;
  DataQuality & operator = (DataQuality &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  TTree *     fSpillTrailerTree;     ///< Tree holding the data from the SpillTrailer fragments
  TTree *     fCaenV1740DataTree;    ///< Tree holding the data from the CAEN V1740 fragments 
  TTree *     fCaenV1751DataTree;    ///< Tree holding the data from the CAEN V1751 fragments 
  TTree *     fWutDataTree;          ///< Tree holding the data from the Wave Union TDC fragments 
  TTree *     fMwpcTdcDataTree;      ///< Tree holding the data from the MWPC TDC fragments 
  std::string fRawFragmentLabel;     ///< label for module producing artdaq fragments
  std::string fRawFragmentInstance;  ///< instance label for artdaq fragments        

  // variables that will go into fSpillTrailerTree
  uint32_t runNumber;
  uint32_t spillNumber;
  uint32_t timeStamp;

  // spill number that goes into fCaenV1740DataTree, fCaenV1751DataTree, fMwpcTdcDataTree, fWutDataTree
  uint32_t spill;

  // variables that will go into fCaenV1740DataTree and/or fCaenV1751DataTree
  uint32_t caen_fragment;
  uint32_t caen_board_id;
  uint32_t caen_event_counter;
  uint32_t caen_trigger_time_tag;  // Each count in the V1751 trigger time tag is 8 ns
  std::vector< std::vector<uint16_t> > caen_v1751_waveform;
  std::vector< std::vector<uint16_t> > caen_v1740_waveform;

  // variables that will go into fMwpcTdcDataTree
  uint32_t mwpc_trigger_counter;
  uint16_t mwpc_controller_time_stamp;
  uint32_t mwpc_tdc_time_stamp;
  uint32_t mwpc_number_hits;
  std::vector<uint16_t> mwpc_tdc_number;
  std::vector<uint16_t> mwpc_hit_channel;
  std::vector<uint16_t> mwpc_hit_time_bin;

  // variables that will go into fWutDataTree
  uint32_t wut_time_header;  // Each count in the time header is 16 us
  uint32_t wut_number_hits;
  std::vector<uint16_t> wut_hit_channel;
  std::vector<uint32_t> wut_hit_time_bin;

};

//------------------------------------------------------------------------------
DataQuality::DataQuality(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);
}

//------------------------------------------------------------------------------
void DataQuality::reconfigure(fhicl::ParameterSet const & p)
{
  fRawFragmentLabel    = p.get< std::string >("RawFragmentLabel", "daq");
  fRawFragmentInstance = p.get< std::string >("RawFragmentInstance", "SPILL");
}

//------------------------------------------------------------------------------
void DataQuality::beginJob()
{

  caen_v1740_waveform.resize(V1740_N_CHANNELS);
  for (size_t i = 0; i < V1740_N_CHANNELS; ++i) {
    caen_v1740_waveform[i].reserve(V1740_N_SAMPLES);
  }

  caen_v1751_waveform.resize(V1751_N_CHANNELS);
  for (size_t i = 0; i < V1751_N_CHANNELS; ++i) {
    caen_v1751_waveform[i].reserve(V1751_N_SAMPLES);
  }

  mwpc_tdc_number.reserve(TDCFragment::MAX_TDCS * TDCFragment::MAX_HITS);
  mwpc_hit_channel.reserve(TDCFragment::MAX_TDCS * TDCFragment::MAX_HITS);
  mwpc_hit_time_bin.reserve(TDCFragment::MAX_TDCS * TDCFragment::MAX_HITS);

  wut_hit_channel.reserve(4 * WUT_MAX_HITS);
  wut_hit_time_bin.reserve(4 * WUT_MAX_HITS);

  art::ServiceHandle<art::TFileService> tfs;

  fSpillTrailerTree = tfs->make<TTree>("spillTrailer", "spillTrailer");
  fSpillTrailerTree->Branch("runNumber", &runNumber, "runNumber/i");
  fSpillTrailerTree->Branch("spillNumber", &spillNumber, "spillNumber/i");
  fSpillTrailerTree->Branch("timeStamp", &timeStamp, "timeStamp/i");

  fCaenV1740DataTree = tfs->make<TTree>("v1740", "v1740");
  fCaenV1740DataTree->Branch("spill", &spill, "spill/i");
  fCaenV1740DataTree->Branch("fragment", &caen_fragment, "fragment/i");
  fCaenV1740DataTree->Branch("event_counter", &caen_event_counter,
                             "event_counter/i");
  fCaenV1740DataTree->Branch("board_id", &caen_board_id, "board_id/i");
  fCaenV1740DataTree->Branch("trigger_time_tag", &caen_trigger_time_tag,
                             "trigger_time_tag/i");

  for (size_t i = 0; i < V1740_N_CHANNELS; ++i) {
    std::string branch_name = "channel_" + std::to_string(i);
    std::string leaf_list = "channel_" + std::to_string(i) + "[" +
                            std::to_string(V1740_N_SAMPLES) + "]/s";
    //std::cout << "branch_name: " << branch_name << std::endl;
    //std::cout << "leaf_list: " << leaf_list << std::endl;
    fCaenV1740DataTree->Branch(branch_name.c_str(),
                               caen_v1740_waveform[i].data(),
                               leaf_list.c_str());
  }

  fCaenV1751DataTree = tfs->make<TTree>("v1751", "v1751");
  fCaenV1751DataTree->Branch("spill", &spill, "spill/i");
  fCaenV1751DataTree->Branch("fragment", &caen_fragment, "fragment/i");
  fCaenV1751DataTree->Branch("event_counter", &caen_event_counter,
                             "event_counter/i");
  fCaenV1751DataTree->Branch("board_id", &caen_board_id, "board_id/i");
  fCaenV1751DataTree->Branch("trigger_time_tag", &caen_trigger_time_tag,
                             "trigger_time_tag/i");

  for (size_t i = 0; i < V1751_N_CHANNELS; ++i) {
    std::string branch_name = "channel_" + std::to_string(i);
    std::string leaf_list = "channel_" + std::to_string(i) + "[" +
                            std::to_string(V1751_N_SAMPLES) + "]/s";
    //std::cout << "branch_name: " << branch_name << std::endl;
    //std::cout << "leaf_list: " << leaf_list << std::endl;
    fCaenV1751DataTree->Branch(branch_name.c_str(),
                               caen_v1751_waveform[i].data(),
                               leaf_list.c_str());
  }

  fMwpcTdcDataTree = tfs->make<TTree>("mwpc", "mwpc");
  fMwpcTdcDataTree->Branch("spill", &spill, "spill/i");
  fMwpcTdcDataTree->Branch("trigger_counter", &mwpc_trigger_counter,
                           "trigger_counter/i");
  fMwpcTdcDataTree->Branch("controller_time_stamp",
                           &mwpc_controller_time_stamp,
                           "controller_time_stamp/s");
  fMwpcTdcDataTree->Branch("tdc_time_stamp", &mwpc_tdc_time_stamp,
                           "tdc_time_stamp/i");
  fMwpcTdcDataTree->Branch("number_hits", &mwpc_number_hits, "number_hits/i");
  fMwpcTdcDataTree->Branch("tdc_number", mwpc_tdc_number.data(),
                           "tdc_number[number_hits]/s");
  fMwpcTdcDataTree->Branch("hit_channel", mwpc_hit_channel.data(),
                           "hit_channel[number_hits]/s");
  fMwpcTdcDataTree->Branch("hit_time_bin", mwpc_hit_time_bin.data(),
                           "hit_time_bin[number_hits]/s");

  fWutDataTree = tfs->make<TTree>("wut", "wut");
  fWutDataTree->Branch("spill", &spill, "spill/i");
  fWutDataTree->Branch("time_header", &wut_time_header, "time_header/i");
  fWutDataTree->Branch("number_hits", &wut_number_hits, "number_hits/i");
  fWutDataTree->Branch("hit_channel", wut_hit_channel.data(),
                       "hit_channel[number_hits]/s");
  fWutDataTree->Branch("hit_time_bin", wut_hit_time_bin.data(),
                       "hit_time_bin[number_hits]/i");

  return;
}

//------------------------------------------------------------------------------
void DataQuality::analyze(art::Event const & evt)
{
  art::Handle< std::vector<artdaq::Fragment> > fragments;
  evt.getByLabel(fRawFragmentLabel, fRawFragmentInstance, fragments);

  if ( !fragments.isValid() )
      throw cet::exception("LARIATFragementReader")
      << "artdaq::Fragment handle is not valid, bail";
  if ( fragments->size() != 1 )
      throw cet::exception("LARIATFragementReader")
      << "artdaq::Fragment handle contains more than one fragment, bail";

  // get the fragments we are interested in
  const auto& frag((*fragments)[0]);

  const char * bytePtr = reinterpret_cast<const char *> (&*frag.dataBegin());
  LariatFragment * data = new LariatFragment((char *) bytePtr,
      frag.dataSize() * sizeof(unsigned long long));
  std::cout << "Have data fragment "
            << frag.dataSize() * sizeof(unsigned long long)
            << std::endl;
  data->print();
  data->printSpillTrailer();

  LariatFragment::SpillTrailer & spillTrailer = data->spillTrailer;
  runNumber = spillTrailer.runNumber;
  spillNumber = spillTrailer.spillNumber;
  timeStamp = spillTrailer.timeStamp;

  spill = spillNumber;

  std::cout << "evt.run(): " << evt.run() << "; evt.subRun(): " << evt.subRun()
            << "; evt.event(): " << evt.event() << std::endl;
  std::cout << "runNumber: " << runNumber << "; spillNumber: " << spillNumber
            << "; timeStamp: " << timeStamp << std::endl;

  const size_t numberCaenFrags = data->caenFrags.size();
  std::cout << "Found " << numberCaenFrags << " CAEN fragments" << std::endl;

  if (numberCaenFrags > 0) {
    std::cout << "Looking at CAEN fragments..." << std::endl;
  }

  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment & caenFrag = data->caenFrags[i];

    uint32_t board = caenFrag.header.boardId;

    if (board == 0 or board == 1 or board == 2 or
        board == 3 or board == 4 or board == 5 or
        board == 6 or board == 7) {
      //caenFrag.print();

      //std::cout << "CAEN event counter: "
      //          << caenFrag.header.eventCounter
      //          << std::endl;

      caen_fragment = i;
      caen_event_counter = caenFrag.header.eventCounter;
      caen_board_id = caenFrag.header.boardId;
      caen_trigger_time_tag = caenFrag.header.triggerTimeTag;

      for (size_t j = 0; j < V1740_N_CHANNELS; ++j) {
        caen_v1740_waveform[j].clear();
      }

      for (size_t sample = 0; sample < caenFrag.header.nSamples; ++sample) {
        for (size_t j = 0; j < V1740_N_CHANNELS; ++j) {
          caen_v1740_waveform[j].push_back(caenFrag.waveForms[j].data[sample]);
        }
      }

      fCaenV1740DataTree->Fill();
    }
    else if (board == 8 or board == 9) {
      //caenFrag.print();

      //std::cout << "CAEN event counter: "
      //          << caenFrag.header.eventCounter
      //          << std::endl;

      caen_fragment = i;
      caen_event_counter = caenFrag.header.eventCounter;
      caen_board_id = caenFrag.header.boardId;
      caen_trigger_time_tag = caenFrag.header.triggerTimeTag;

      for (size_t j = 0; j < V1751_N_CHANNELS; ++j) {
        caen_v1751_waveform[j].clear();
      }

      for (size_t sample = 0; sample < caenFrag.header.nSamples; ++sample) {
        for (size_t j = 0; j < V1751_N_CHANNELS; ++j) {
          caen_v1751_waveform[j].push_back(caenFrag.waveForms[j].data[sample]);
        }
      }

      fCaenV1751DataTree->Fill();
    }
  }

  const int numberTdcFrags = data->tdcFrags.size();
  std::cout << "Found " << numberTdcFrags << " TDC fragments" << std::endl;

  if (numberTdcFrags > 0) {
    std::cout << "Looking at TDC fragments..." << std::endl;
  }

  for (int i = 0; i < numberTdcFrags; ++i) {

    TDCFragment & tdcFrag = data->tdcFrags[i];

    std::vector< std::vector<TDCFragment::TdcEventData> > &
        tdcEvents = tdcFrag.tdcEvents;

    //std::cout << "tdcEvents.size(): " << tdcEvents.size() << std::endl;
    //tdcFrag.print();

    for (size_t j = 0; j < tdcEvents.size(); ++j) {

      if (tdcFrag.controllerHeader.nTDCs != tdcEvents[j].size()) {
        std::cout << "*** Fatal nTDCs mismatch: " << tdcEvents[j].size()
                  << " != " << tdcFrag.controllerHeader.nTDCs << std::endl;
      }

      mwpc_tdc_number.clear();
      mwpc_hit_channel.clear();
      mwpc_hit_time_bin.clear();

      mwpc_number_hits = 0;

      for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS;
           ++tdc_index) {
        TDCFragment::TdcEventData tdcEventData = tdcEvents[j].at(tdc_index);

        mwpc_trigger_counter = tdcEventData.tdcEventHeader.triggerCounter;

        //std::cout << "triggerCounter: " << mwpc_trigger_counter << std::endl;

        mwpc_controller_time_stamp = tdcEventData.tdcEventHeader.controllerTimeStamp;
        mwpc_tdc_time_stamp = tdcEventData.tdcEventHeader.tdcTimeStamp;

        mwpc_number_hits += tdcEventData.tdcEventHeader.nHits;

        for (size_t hit_index = 0; hit_index < tdcEventData.tdcHits.size();
             ++hit_index) {
          TDCFragment::TdcHit & hit = tdcEventData.tdcHits[hit_index];
          mwpc_hit_channel.push_back(uint16_t (hit.channel));
          mwpc_hit_time_bin.push_back(hit.timeBin);
          mwpc_tdc_number.push_back(tdc_index + 1);
        }
      }

      if (mwpc_number_hits == 0) continue;

      fMwpcTdcDataTree->Fill();

    }
  }

  const int numberWutFrags = data->wutFrags.size();
  std::cout << "Found " << numberWutFrags << " WUT fragments" << std::endl;

  if (numberWutFrags > 0) {
    std::cout << "Looking at WUT fragments..." << std::endl;
  }

  for (int i = 0; i < numberWutFrags; ++i) {
    WUTFragment & wutFrag = data->wutFrags[i];
    //wutFrag.print();
    uint32_t numberHits = wutFrag.header.nHits;
    std::vector<WUTFragment::WutHit> & hits = wutFrag.hits;

    wut_number_hits = wutFrag.header.nHits;

    wut_hit_channel.clear();
    wut_hit_time_bin.clear();

    wut_time_header = wutFrag.header.timeHeader;

    for (size_t j = 0; j < numberHits; ++j) {
      WUTFragment::WutHit & hit = hits[j];

      wut_hit_channel.push_back(uint16_t (hit.channel));
      wut_hit_time_bin.push_back(hit.timeBin);

    }

    fWutDataTree->Fill();

  }

  fSpillTrailerTree->Fill();

  return;  
}

DEFINE_ART_MODULE(DataQuality)
