////////////////////////////////////////////////////////////////////////
// Class:       LArIATFragmentReader
// Module Type: analyzer
// File:        LArIATFragmentReader_module.cc
//
// Generated at Tue Oct 28 14:15:01 2014 by Brian_Rebel using artmod
// from cetpkgsupport v1_07_01.
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

#include "daq/include/LariatFragment.h"
#include "daq/include/WUTFragment.h"
#include "daq/include/CAENFragment.h"
#include "daq/lariat-artdaq/lariat-artdaq/Overlays/GenericFragment.hh"
#include "daq/lariat-artdaq/lariat-artdaq/Overlays/WUTWrapperFragment.hh"
#include "daq/lariat-artdaq/lariat-artdaq/Overlays/CAENWrapperFragment.hh"
#include "daq/lariat-artdaq/lariat-artdaq/ArtModules/ModuleUtils.hh"

#include "TTree.h"

#include <vector>
#include <string>

namespace rdu {
  class LArIATFragmentReader;

  class WUTData {

  public:
    uint32_t time_header;  // Each count in the time header is 16 us
    std::vector<uint16_t> hit_channel;
    std::vector<uint32_t> hit_time_bin;
    std::vector<uint64_t> hit_time;
  };
  
  class CAENData {

  public:
    uint32_t trigger_time_tag;  // Each count in the trigger time tag is 8 ns
    std::vector<uint16_t> ustof1_logic;
    std::vector<uint16_t> ustof2_logic;
    std::vector<uint16_t> ustof3_logic;
    std::vector<uint16_t> ustof4_logic;
    std::vector<uint16_t> dstof1_logic;
    std::vector<uint16_t> dstof2_logic;
  };
  
  typedef struct{
    uint16_t spill;
    uint16_t fragment_id;
  } SpillInfo;
}


class rdu::LArIATFragmentReader : public art::EDAnalyzer {
public:
  explicit LArIATFragmentReader(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArIATFragmentReader(LArIATFragmentReader const &) = delete;
  LArIATFragmentReader(LArIATFragmentReader &&) = delete;
  LArIATFragmentReader & operator = (LArIATFragmentReader const &) = delete;
  LArIATFragmentReader & operator = (LArIATFragmentReader &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  void FillCAENInfo(std::vector<const uint8_t*>& caenFragPtrs);
  void FillWUTInfo(std::vector<const uint8_t*>&  wutFragPtrs);

  TTree                *fDataTree;            ///< Tree holding the data from the various fragments 
  CAENData              fCAEN;    	      ///< data from CAEN V1751				 
  WUTData    	        fWUT;      	      ///< data from WUT				       	 
  SpillInfo   	        fSpill;    	      ///< data from the event/spill                       
  std::vector<CAENData> fCAENs;               ///< collection of all CAEN fragments
  std::vector<WUTData>  fWUTs;                ///< collection of all CAEN fragments
  std::string           fRawFragmentLabel;    ///< label for module producing artdaq fragments
  std::string  		fRawFragmentInstance; ///< instance label for artdaq fragments        
};


//------------------------------------------------------------------------------
rdu::LArIATFragmentReader::LArIATFragmentReader(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  this->reconfigure(p);
}

//------------------------------------------------------------------------------
void rdu::LArIATFragmentReader::reconfigure(fhicl::ParameterSet const & p)
{
  fRawFragmentLabel    = p.get< std::string >("RawFragmentLabel", "daq");
  fRawFragmentInstance = p.get< std::string >("RawFragmentInstance", "SPILL");
}

//------------------------------------------------------------------------------
void rdu::LArIATFragmentReader::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fDataTree = tfs->make<TTree>("LArIATData", "LArIATData");
  fDataTree->Branch("event", &fSpill, "spill/i:fragment_id/i");
  fDataTree->Branch("caen", &fCAEN);
  fDataTree->Branch("wut",   &fWUT);

  return;
}

//------------------------------------------------------------------------------
void rdu::LArIATFragmentReader::analyze(art::Event const & e)
{
  art::Handle< std::vector<artdaq::Fragment> > fragments;
  e.getByLabel(fRawFragmentLabel, fRawFragmentInstance, fragments);

  if( !fragments.isValid() )
    throw cet::exception("LARIATFragementReader") << "artdaq::Fragment handle is not valid, bail";
  if( fragments->size() != 1 )
    throw cet::exception("LARIATFragementReader") << "artdaq::Fragment handle contains more than one fragment, bail";

  // get the fragments we are interested in
  const auto& frag((*fragments)[0]);

  std::vector<const uint8_t*> wutFragPtrs;
  std::vector<const uint8_t*> caenFragPtrs;

  // first the WUT fragments
  lariat::getLariatFragments(frag, 
			     LariatFragment::LariatFragmentType::FRAGMENT_TYPE_WUT_HIT_HEADER,
			     wutFragPtrs);
  this->FillWUTInfo(wutFragPtrs);

  // now get the CAEN fragments
  lariat::getLariatFragments(frag, 
			     LariatFragment::LariatFragmentType::FRAGMENT_TYPE_CAEN_ADC,
			     caenFragPtrs);
  this->FillCAENInfo(caenFragPtrs);


  // fill the trees
  if( fWUTs.size() != fCAENs.size() )
    throw cet::exception("LArIATFragmentReader") << "Different number of CAEN and WUT fragments - is that expected?";

  fSpill.spill = e.id().event();

  for(size_t s = 0; s < fCAENs.size(); ++s){
    fWUT  = fWUTs[s];
    fCAEN = fCAENs[s];

    fDataTree->Fill();
  }
  
  return;  
}

//------------------------------------------------------------------------------
void rdu::LArIATFragmentReader::FillCAENInfo(std::vector<const uint8_t*>& caenFragPtrs)
{
  const size_t numberCaenFrags = caenFragPtrs.size();
  LOG_VERBATIM("LArIATFragmentReader") << "Found " << numberCaenFrags << " CAEN fragments";
  
  LOG_VERBATIM("LArIATFragmentReader") << "Looking at CAEN fragments...";

  fCAENs.clear();
  CAENData data;
  for(size_t i = 0; i < numberCaenFrags; ++i) {
    lariat::CAENWrapperFragment caenWrapper(caenFragPtrs[i]);
    CAENFragment const& frag = *(caenWrapper.GetCAENFragment());
    
    if (frag.header.boardId == 8) {

      fSpill.fragment_id = i;
      data.trigger_time_tag = frag.header.triggerTimeTag;
      
      LOG_VERBATIM("LArIATFragmentReader") << "///////////////////////////////////////"
					   << "\nFragment number: " << fSpill.fragment_id 
					   << "\n///////////////////////////////////////"
					   << "\nBoard ID: " << frag.header.boardId
					   << "\nNumber of samples: " << frag.header.nSamples;
      data.ustof1_logic.clear();
      data.ustof2_logic.clear();
      data.ustof3_logic.clear();
      data.ustof4_logic.clear();
      data.dstof1_logic.clear();
      data.dstof2_logic.clear();
      
      for(size_t time = 0; time < frag.header.nSamples; ++time){
      	data.ustof1_logic.push_back(frag.waveForms[2].data[time]);
      	data.ustof2_logic.push_back(frag.waveForms[3].data[time]);
      	data.ustof3_logic.push_back(frag.waveForms[4].data[time]);
      	data.ustof4_logic.push_back(frag.waveForms[5].data[time]);
      	data.dstof1_logic.push_back(frag.waveForms[6].data[time]);
      	data.dstof2_logic.push_back(frag.waveForms[7].data[time]);
      }
      
      fCAENs.push_back(data);
      
    }// end if board ID is 8

  }// end loop over fragment pointers

  return;
}

//------------------------------------------------------------------------------
void rdu::LArIATFragmentReader::FillWUTInfo(std::vector<const uint8_t*>& wutFragPtrs)
{

  const size_t numberWutFrags = wutFragPtrs.size();
  LOG_VERBATIM("LArIATFragmentReader") << "Found " << numberWutFrags << " WUT fragments";

  fWUTs.clear();
  WUTData data;

  for(size_t i = 0; i < numberWutFrags; ++i){
    lariat::WUTWrapperFragment wutWrapper(wutFragPtrs[i]);

    WUTFragment const& frag = *(wutWrapper.GetWUTFragment());

    size_t numberHits = frag.header.nHits;

    LOG_VERBATIM("LArIATFragmentReader") << "  WUT fragment ID: " << i 
					 << "\n  Number of WUT hits: " << numberHits;

    data.hit_channel.clear();
    data.hit_time_bin.clear();
    data.hit_time.clear();
    
    data.time_header = frag.header.timeHeader;

    for (size_t j = 0; j < numberHits; ++j) {
      WUTFragment::WutHit const& hit = frag.hits[j];
      
      // hit time since beginning of spill
      uint64_t hitTime = ((uint64_t) data.time_header << 20) | ((uint64_t) hit.timeBin);
      
      data.hit_channel.push_back((uint16_t) hit.channel);
      data.hit_time_bin.push_back(hit.timeBin);
      data.hit_time.push_back(hitTime);

      LOG_VERBATIM("LArIATFragmentReader") << "    Hit ID: " << j
					   << "\n   Channel: " << (uint16_t) hit.channel
					   << "\n  Time bin: " << hit.timeBin 
					   << "\n  Time since BOS: " << hitTime
					   << "\n  Time since BOS (s): " << hitTime * 15.625e-12 ;

    }// end loop over WUT hits

    fWUTs.push_back(data);
  }// end loop over WUT fragment pointers

  return;
}

DEFINE_ART_MODULE(rdu::LArIATFragmentReader)
