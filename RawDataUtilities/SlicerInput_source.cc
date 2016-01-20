//////////////////////////////////////////////////////////////
// Name:      SlicerInput_source.cc
// Date:      8 September 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef SlicerInput_source
#define SlicerInput_source

// art includes
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Root/rootNames.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// artdaq includes
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/Fragments.hh"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/RawDigit.h"
#include "RawData/TriggerData.h"
#include "SimpleTypesAndConstants/RawTypes.h"
#include "SummaryData/RunData.h"
#include "Utilities/AssociationUtil.h"

// LArIATFragment includes
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/TDCFragment.h"
#include "LArIATFragments/V1495Fragment.h"
#include "LArIATFragments/WUTFragment.h"

// LArIATSoft includes
#include "LArIATDataProducts/ConditionsSummary.h"
#include "RawDataUtilities/SlicerAlg.h"
#include "RawDataUtilities/FragmentToDigitAlg.h"
#include "Utilities/DatabaseUtilityT1034.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>
#include <utility>
#include <initializer_list>

//-------------------------------------------------------------------------
// unnamed namespace
namespace {

  // Retrieves branch name (a la art convention) where object resides
  template <typename PROD>
  const char * getBranchName(art::InputTag const& tag)
  {
    std::ostringstream oss;
    oss << art::TypeID(typeid(PROD)).friendlyClassName()
          << '_'
          << tag.label()
          << '_'
          << tag.instance()
          << '_'
          << tag.process()
          << ".obj";
    return oss.str().data();
  }

  artdaq::Fragments * getFragments(TBranch * br, unsigned entry)
  {
    br->GetEntry(entry);
    return reinterpret_cast <artdaq::Fragments *> (br->GetAddress());
  }

  // Assumed file format is
  //
  //     "lariat_r[digits]_sr[digits]_other_stuff.root"
  //
  // This regex object is used to extract the run and subrun numbers.
  // The '()' groupings correspond to regex matches that can be
  // extracted using the std::regex_match facility.

  std::regex const filename_format(".*\\lariat_r(\\d+)_sr(\\d+).*\\.root");

}

//-------------------------------------------------------------------------
namespace raw {

  // Enable 'pset.get<raw::Compress_t>("compression")'
  void decode(boost::any const & a, Compress_t & result) {
    unsigned tmp;
    fhicl::detail::decode(a, tmp);
    result = static_cast <Compress_t> (tmp);
  }
}

//-------------------------------------------------------------------------
namespace rdu
{
  // The class Slicer is to be used as the template parameter for
  // art::Source<T>. It understands how to read art's ROOT data files,
  // to extract artdaq::Fragments from them, how to convert the
  // artdaq::Fragments into vectors of raw::RawDigit objects, and
  // finally how to re-combine those raw::RawDigit objects to present
  // the user of the art::Source<Slicer> with a sequence of
  // art::Events that have the desired event structure, different from
  // the event structure of the data file(s) being read.

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class Slicer
  {

   public:

    // Constructor and destructor.
    explicit Slicer(fhicl::ParameterSet        const& pset,
                    art::ProductRegistryHelper      & prhelper,
                    art::SourceHelper               & shelper);
    virtual ~Slicer();

    ///////////////////////////////////////////////////////////////////
    // See art/Framework/IO/Sources/Source.h for a description of each
    // of the public member functions of Slicer.
    ///////////////////////////////////////////////////////////////////

    // Open the file of the given name, returning a new FileBlock.
    // If readFile is unable to return a valid FileBlock it should
    // throw.
    bool readFile(std::string const& filename, art::FileBlock * & fileblock);

    // Read the next part of the current file. Return false if nothing
    // was read; return true and set the appropriate 'out' arguments if
    // something was read.
    bool readNext(art::RunPrincipal    * const& inRun,
                  art::SubRunPrincipal * const& inSubRun,
                  art::RunPrincipal    *      & outRun,
                  art::SubRunPrincipal *      & outSubRun,
                  art::EventPrincipal  *      & outEvent);

    // Close the current file.
    void closeCurrentFile();

    ///////////////////////////////////////////////////////////////////

    // Read in any parameters from the .fcl files.
    void reconfigure(fhicl::ParameterSet const& pset);

   private:

    // Private member functions are appended with an underscore.

    void loadDigits_(LariatFragment * & LArIATFragment);

    void makeEventAndPutDigits_(art::EventPrincipal * & outEvent);

    // Data members are prepended with the letter f.

    std::string            fSourceName;
    std::string            fLastFileName;
    std::unique_ptr<TFile> fFile;
    bool                   fDoneWithFile;
    art::InputTag          fInputTag;
    art::SourceHelper      fSourceHelper;
    TBranch *              fFragmentsBranch;
    TBranch *              fEventAuxBranch;
    size_t                 fNumberInputEvents;
    size_t                 fTreeIndex;
    art::RunNumber_t       fRunNumber;
    art::SubRunNumber_t    fSubRunNumber;
    art::EventNumber_t     fEventNumber;
    art::RunNumber_t       fCachedRunNumber;
    art::SubRunNumber_t    fCachedSubRunNumber;

    LariatFragment * fLariatFragment;

    // collection of data blocks
    std::vector< rdu::DataBlockCollection > fCollections;
    size_t fCollectionIndex;

    // DatabaseUtilityT1034 service handle
    art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

    // slicer algorithm
    rdu::SlicerAlg fSlicerAlg;

    // fragment-to-digit algorithm (?)
    FragmentToDigitAlg fFragmentToDigitAlg;

    // load up the RunPrincipal
    void commenceRun(art::RunPrincipal * & outRun);

    // load up the SubRunPrincipal
    void commenceSubRun(art::SubRunPrincipal * & outSubRun);
    
    // get parameters from lariat_prd database
    void getDatabaseParameters_(art::RunNumber_t const& RunNumber);

    // cast from string to size_t
    size_t castToSizeT_(std::string const& String);

    // vector of parameter names to be queried from the
    // lariat_xml_database table for a specified run
    std::vector<std::string> fConfigParams;
    // (key, value) pair for the database query result
    std::map< std::string, std::string > fConfigValues;

    // used to filter out events with number of TPC readouts outside the range of
    // fMinTPCReadoutsPerEvent <= NumberTPCReadouts <= fMaxTPCReadoutsPerEvent
    unsigned int fMinTPCReadoutsPerEvent;
    unsigned int fMaxTPCReadoutsPerEvent;

    // V1495 delay for the delayed trigger
    size_t fV1495DelayTicks;
    double fV1495Delay;

    // sample reduction of the CAEN V1740 digitizers
    size_t fV1740SampleReduction;
    size_t fV1740BSampleReduction;

    // post percent of CAEN digitizers
    double fV1740PostPercent;
    double fV1740BPostPercent;
    double fV1751PostPercent;

    // record lengths (number of time ticks) of the CAEN
    // V1740 and V1751 digitizers
    size_t fV1740RecordLength;
    size_t fV1740BRecordLength;
    size_t fV1751RecordLength;

    // sampling rate in MHz
    double fV1740SamplingRate;
    double fV1740BSamplingRate;
    double fV1751SamplingRate;

    // readout window lengths in microseconds
    double fV1740ReadoutWindow;
    double fV1740BReadoutWindow;
    double fV1751ReadoutWindow;
    double fTDCReadoutWindow;

    // TDC parameters
    double fTDCPipelineDelay;
    double fTDCGateWidth;

    // pre-/post-acquisition and acquisition windows
    double fV1740PreAcquisitionWindow;
    double fV1740PostAcquisitionWindow;
    double fV1740AcquisitionWindow;
    double fV1740BPreAcquisitionWindow;
    double fV1740BPostAcquisitionWindow;
    double fV1740BAcquisitionWindow;
    double fV1751PreAcquisitionWindow;
    double fV1751PostAcquisitionWindow;
    double fV1751AcquisitionWindow;
    double fTDCPreAcquisitionWindow;
    double fTDCPostAcquisitionWindow;
    double fTDCAcquisitionWindow;

    // timestamp from SpillTrailer
    std::uint64_t fTimestamp;

  }; // class Slicer

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  Slicer::Slicer(fhicl::ParameterSet        const& pset,
                 art::ProductRegistryHelper      & prhelper,
                 art::SourceHelper               & shelper)
    : fSourceName("daq")
    , fLastFileName(pset.get< std::vector< std::string > >("fileNames", {}).back())
    , fFile()
    , fDoneWithFile(false)
    , fInputTag("daq:SPILL:DAQ")
    , fSourceHelper(shelper)
    , fFragmentsBranch(nullptr)
    , fNumberInputEvents()
    , fRunNumber(1)       // Defaults in case input filename does not
    , fSubRunNumber(0)    // follow assumed filename_format above.
    , fEventNumber()
    , fCachedRunNumber(-1)
    , fCachedSubRunNumber(-1)
    , fSlicerAlg(pset.get<fhicl::ParameterSet>("SlicerAlg"))
    , fFragmentToDigitAlg(pset.get<fhicl::ParameterSet>("FragmentToDigitAlg"))
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);

    prhelper.reconstitutes<std::vector<raw::AuxDetDigit>, art::InEvent >(fSourceName);
    prhelper.reconstitutes<std::vector<raw::RawDigit>,    art::InEvent >(fSourceName);
    prhelper.reconstitutes<std::vector<raw::OpDetPulse>,  art::InEvent >(fSourceName);
    prhelper.reconstitutes<std::vector<raw::Trigger>,     art::InEvent >(fSourceName);
    prhelper.reconstitutes<sumdata::RunData,              art::InRun   >(fSourceName);
    prhelper.reconstitutes<ldp::ConditionsSummary,        art::InSubRun>(fSourceName);

    // set config parameters to get from the lariat_prd database
    fConfigParams.clear();
    fConfigParams.push_back("v1495_config_v1495_delay_ticks");
    fConfigParams.push_back("v1740_config_caen_postpercent");
    fConfigParams.push_back("v1740_config_caen_recordlength");
    fConfigParams.push_back("v1740b_config_caen_postpercent");
    fConfigParams.push_back("v1740b_config_caen_recordlength");
    fConfigParams.push_back("v1751_config_caen_postpercent");
    fConfigParams.push_back("v1751_config_caen_recordlength");
    fConfigParams.push_back("v1740_config_caen_v1740_samplereduction");
    fConfigParams.push_back("v1740b_config_caen_v1740_samplereduction");
    fConfigParams.push_back("tdc_config_tdc_pipelinedelay");
    fConfigParams.push_back("tdc_config_tdc_gatewidth");
  }

  //-----------------------------------------------------------------------
  // destructor
  Slicer::~Slicer()
  {}

  //-----------------------------------------------------------------------
  void Slicer::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSlicerAlg.reconfigure(pset.get<fhicl::ParameterSet>("SlicerAlg"));

    fSourceName = pset.get< std::string >("SourceName", "daq");

    fMinTPCReadoutsPerEvent = pset.get< unsigned int >("MinTPCReadoutsPerEvent", 0);
    fMaxTPCReadoutsPerEvent = pset.get< unsigned int >("MaxTPCReadoutsPerEvent", 100);

    fV1740PreAcquisitionWindow   = pset.get< double >("V1740PreAcquisitionWindow",   -1);
    fV1740PostAcquisitionWindow  = pset.get< double >("V1740PostAcquisitionWindow",  -1);
    fV1740AcquisitionWindow      = pset.get< double >("V1740AcquisitionWindow",      -1);
    fV1740BPreAcquisitionWindow  = pset.get< double >("V1740BPreAcquisitionWindow",  -1);
    fV1740BPostAcquisitionWindow = pset.get< double >("V1740BPostAcquisitionWindow", -1);
    fV1740BAcquisitionWindow     = pset.get< double >("V1740BAcquisitionWindow",     -1);
    fV1751PreAcquisitionWindow   = pset.get< double >("V1751PreAcquisitionWindow",   -1);
    fV1751PostAcquisitionWindow  = pset.get< double >("V1751PostAcquisitionWindow",  -1);
    fV1751AcquisitionWindow      = pset.get< double >("V1751AcquisitionWindow",      -1);
    fTDCPreAcquisitionWindow     = pset.get< double >("TDCPreAcquisitionWindow",     -1);
    fTDCPostAcquisitionWindow    = pset.get< double >("TDCPostAcquisitionWindow",    -1);
    fTDCAcquisitionWindow        = pset.get< double >("TDCAcquisitionWindow",        -1);

    return;
  }

  //-----------------------------------------------------------------------
  bool Slicer::readFile(std::string const& filename, art::FileBlock * & fileblock)
  {
    // Run numbers determined based on file name... see comment in
    // unnamed namespace above.
    std::smatch matches;
    if (std::regex_match(filename, matches, filename_format)) {
      fRunNumber    = std::stoul(matches[1]);
      fSubRunNumber = std::stoul(matches[2]);
    }

    LOG_VERBATIM("SlicerInput")
    << "\n////////////////////////////////////"
    << "\nfRunNumber:       " << fRunNumber
    << "\nfSubRunNumber:    " << fSubRunNumber
    << "\nfCachedRunNumber: " << fCachedRunNumber
    << "\n////////////////////////////////////\n";

    // get database parameters
    // TODO: Check to see if the current run number is the same as the
    //       previous run number. If it is, don't bother querying the
    //       database again since it will return the same results as
    //       before.
    this->getDatabaseParameters_(fRunNumber);

    // configure the slicer algorithm
    fSlicerAlg.Configure(fV1740PreAcquisitionWindow,
                         fV1740PostAcquisitionWindow,
                         fV1740AcquisitionWindow,
                         fV1740BPreAcquisitionWindow,
                         fV1740BPostAcquisitionWindow,
                         fV1740BAcquisitionWindow,
                         fV1751PreAcquisitionWindow,
                         fV1751PostAcquisitionWindow,
                         fV1751AcquisitionWindow,
                         fTDCPreAcquisitionWindow,
                         fTDCPostAcquisitionWindow,
                         fTDCAcquisitionWindow);

    // get artdaq::Fragments branch
    fFile.reset(new TFile(filename.data()));
    TTree * eventTree  = reinterpret_cast <TTree *> (fFile->Get(art::rootNames::eventTreeName().c_str()));
    fFragmentsBranch   = eventTree->GetBranch(getBranchName<artdaq::Fragments>(fInputTag)); // get branch for specific input tag
    //fEventAuxBranch    = eventTree->GetBranch("EventAuxiliary");
    fNumberInputEvents = static_cast <size_t> (fFragmentsBranch->GetEntries()); // Number of fragment-containing events to read in from input file
    fTreeIndex         = 0ul;

    // new fileblock
    fileblock = new art::FileBlock(art::FileFormatVersion(), filename);
    if (fileblock == nullptr) {
      throw art::Exception(art::errors::FileOpenError)
        << "Unable to open file " << filename << ".\n";
    }

    // reset flag
    fDoneWithFile = false;

    // load the data blocks
    this->loadDigits_(fLariatFragment);

    // group data blocks into collections
    fCollectionIndex = 0;
    fCollections.clear();
    fCollections = fSlicerAlg.Slice(fLariatFragment);

    // we are done with this file if there are no data blocks
    if (fCollections.size() < 1) fDoneWithFile = true;

    delete fLariatFragment;

    return true;
  }

  //-----------------------------------------------------------------------
  bool Slicer::readNext(art::RunPrincipal    * const& inRun,
                        art::SubRunPrincipal * const& inSubRun,
                        art::RunPrincipal    *      & outRun,
                        art::SubRunPrincipal *      & outSubRun,
                        art::EventPrincipal  *      & outEvent)
  {
    if (fDoneWithFile) return false;

    art::Timestamp timestamp = fTimestamp;

    if (fRunNumber != fCachedRunNumber) {
      outRun = fSourceHelper.makeRunPrincipal(fRunNumber, timestamp);
      fCachedRunNumber = fRunNumber;
      fEventNumber = 0ul;

      this->commenceRun(outRun);
    }

    if (fSubRunNumber != fCachedSubRunNumber) {
      outSubRun = fSourceHelper.makeSubRunPrincipal(fRunNumber, fSubRunNumber, timestamp);
      fCachedSubRunNumber = fSubRunNumber;
      // fEventNumber = 0ul;
      
      this->commenceSubRun(outSubRun);
    }

    LOG_VERBATIM("SlicerInput") << "fCollections.size(): " << fCollections.size();

    this->makeEventAndPutDigits_(outEvent);

    return true;
  }

  //-----------------------------------------------------------------------
  void Slicer::closeCurrentFile()
  {
    fFile.reset(nullptr);
  }

  //-----------------------------------------------------------------------
  void Slicer::loadDigits_(LariatFragment * & LArIATFragment)
  {

    if (fTreeIndex != fNumberInputEvents) {
      artdaq::Fragments * fragments = getFragments(fFragmentsBranch, fTreeIndex++);

      if ((*fragments).size() > 1)
        throw cet::exception("Slicer") << "artdaq::Fragment vector contains more than one fragment.";

      artdaq::Fragment frag = fragments->at(0);
      const char * bytePtr = reinterpret_cast <const char *> (&*frag.dataBegin());
      LArIATFragment = new LariatFragment((char *) bytePtr, frag.dataSize() * sizeof(unsigned long long));

      // get SpillTrailer
      LariatFragment::SpillTrailer const& spillTrailer = fLariatFragment->spillTrailer;

      // get timestamp from SpillTrailer, cast as uint64_t
      fTimestamp = (static_cast <std::uint64_t> (spillTrailer.timeStamp));
    }

  }

  //-----------------------------------------------------------------------
  void Slicer::makeEventAndPutDigits_(art::EventPrincipal * & outEvent)
  {

    ++fEventNumber;

    outEvent = fSourceHelper.makeEventPrincipal(fRunNumber, fSubRunNumber, fEventNumber, art::Timestamp());

    std::vector<raw::AuxDetDigit> auxDigits;
    std::vector<raw::RawDigit>    rawDigits;
    std::vector<raw::OpDetPulse>  opDetPulses;

    // will deal with this later

    //bool processed = false;

    //while (!processed) {

    //  if (fCollectionIndex >= fCollections.size()) break;

    //  rdu::DataBlockCollection const& Collection = fCollections[fCollectionIndex];

    //  size_t const& NumberTPCReadouts = Collection.numberTPCReadouts;

    //  std::cout << "fCollectionIndex: " << fCollectionIndex << std::endl;
    //  std::cout << "NumberTPCReadouts: " << NumberTPCReadouts << std::endl;
    //  std::cout << "Collection.numberTPCReadouts: " << Collection.numberTPCReadouts << std::endl;
    //  std::cout << "Collection.caenBlocks.size(): " << Collection.caenBlocks.size() << std::endl;
    //  std::cout << "Collection.tdcBlocks.size(): " << Collection.tdcBlocks.size() << std::endl;

    //  if ((fMinTPCReadoutsPerEvent <= NumberTPCReadouts) and (fMaxTPCReadoutsPerEvent >= NumberTPCReadouts)) {

    //    std::cout << "Adding digits to run " << fRunNumber
    //              << ", sub-run "            << fSubRunNumber
    //              <<  ", event "             << fEventNumber
    //              << "..."                   << std::endl;

    //    std::vector< CAENFragment > caenDataBlocks = Collection.caenBlocks;
    //    std::vector< std::vector<TDCFragment::TdcEventData> > tdcDataBlocks = Collection.tdcBlocks;

    //    // will deal with this later
    //    //std::vector<raw::Trigger> triggerVector = makeTheTriggers(caenDataBlocks,
    //    //                                                          tdcDataBlocks);

    //    fFragmentToDigitAlg.makeTheDigits(caenDataBlocks,
    //                                      tdcDataBlocks,
    //                                      auxDigits,
    //                                      rawDigits,
    //                                      opDetPulses);

    //    art::put_product_in_principal(std::make_unique< std::vector<raw::AuxDetDigit> > (auxDigits),
    //                                  * outEvent,
    //                                  fSourceName);
    //    art::put_product_in_principal(std::make_unique< std::vector<raw::RawDigit> > (rawDigits),
    //                                  * outEvent,
    //                                  fSourceName);
    //    art::put_product_in_principal(std::make_unique< std::vector<raw::OpDetPulse> > (opDetPulses),
    //                                  * outEvent,
    //                                  fSourceName);

    //    processed = true;
    //  }

    //  ++fCollectionIndex;

    //}

    rdu::DataBlockCollection const& Collection = fCollections[fCollectionIndex];

    size_t const& NumberTPCReadouts = Collection.numberTPCReadouts;

    LOG_VERBATIM("SlicerInput")
    << "fCollectionIndex: " << fCollectionIndex
    << "\nNumberTPCReadouts: " << NumberTPCReadouts
    << "\nCollection.numberTPCReadouts: " << Collection.numberTPCReadouts
    << "\nCollection.caenBlocks.size(): " << Collection.caenBlocks.size()
    << "\nCollection.tdcBlocks.size(): " << Collection.tdcBlocks.size()
    << "\n\nAdding digits to run " << fRunNumber
    << ", sub-run "                << fSubRunNumber
    <<  ", event "                 << fEventNumber
    << "...";

    std::vector< CAENFragment > caenDataBlocks = Collection.caenBlocks;
    std::vector< std::vector<TDCFragment::TdcEventData> > tdcDataBlocks = Collection.tdcBlocks;

    // will deal with this strange function later...
    std::vector<raw::Trigger> triggerVector = fFragmentToDigitAlg.makeTheTriggers(fEventNumber,
                                                                                  caenDataBlocks,
                                                                                  tdcDataBlocks);

    fFragmentToDigitAlg.makeTheDigits(caenDataBlocks,
                                      tdcDataBlocks,
                                      auxDigits,
                                      rawDigits,
                                      opDetPulses);

    // if there are no TPC readouts, clear the rawDigits vector
    // since we do not want any partial TPC readouts
    if (NumberTPCReadouts < 1) rawDigits.clear();

    art::put_product_in_principal(std::make_unique<std::vector<raw::Trigger> > (triggerVector),
                                  * outEvent,
                                  fSourceName);
    art::put_product_in_principal(std::make_unique< std::vector<raw::AuxDetDigit> > (auxDigits),
                                  * outEvent,
                                  fSourceName);
    art::put_product_in_principal(std::make_unique< std::vector<raw::RawDigit> > (rawDigits),
                                  * outEvent,
                                  fSourceName);
    art::put_product_in_principal(std::make_unique< std::vector<raw::OpDetPulse> > (opDetPulses),
                                  * outEvent,
                                  fSourceName);

    ++fCollectionIndex;

    if (fCollectionIndex >= fCollections.size()) {
      fDoneWithFile = true;
    }
    else {
      fDoneWithFile = false;
    }

    return;
  }

  //-----------------------------------------------------------------------
  void Slicer::getDatabaseParameters_(art::RunNumber_t const& RunNumber)
  {
    fConfigValues.clear();

    fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams,
                                                      static_cast <int> (RunNumber));

    LOG_VERBATIM("SlicerInput")
    << "//////////////////////////////////////////////"
    << "V1495DelayTicks:       " << fConfigValues["v1495_config_v1495_delay_ticks"]           
    << "V1740PostPercent:      " << fConfigValues["v1740_config_caen_postpercent"]            
    << "V1740BPostPercent:     " << fConfigValues["v1740b_config_caen_postpercent"]           
    << "V1751PostPercent:      " << fConfigValues["v1751_config_caen_postpercent"]            
    << "V1740RecordLength:     " << fConfigValues["v1740_config_caen_recordlength"]           
    << "V1740BRecordLength:    " << fConfigValues["v1740b_config_caen_recordlength"]          
    << "V1751RecordLength:     " << fConfigValues["v1751_config_caen_recordlength"]           
    << "V1740SampleReduction:  " << fConfigValues["v1740_config_caen_v1740_samplereduction"]  
    << "V1740BSampleReduction: " << fConfigValues["v1740b_config_caen_v1740_samplereduction"] 
    << "fTDCPipelineDelay:     " << fConfigValues["tdc_config_tdc_pipelinedelay"]             
    << "fTDCGateWidth:         " << fConfigValues["tdc_config_tdc_gatewidth"]                 
    << "//////////////////////////////////////////////";

    // cast from string to size_t
    fV1495DelayTicks       = this->castToSizeT_(fConfigValues["v1495_config_v1495_delay_ticks"]);
    fV1740PostPercent      = this->castToSizeT_(fConfigValues["v1740_config_caen_postpercent"]);
    fV1740BPostPercent     = this->castToSizeT_(fConfigValues["v1740b_config_caen_postpercent"]);
    fV1751PostPercent      = this->castToSizeT_(fConfigValues["v1751_config_caen_postpercent"]);
    fV1740RecordLength     = this->castToSizeT_(fConfigValues["v1740_config_caen_recordlength"]);
    fV1740BRecordLength    = this->castToSizeT_(fConfigValues["v1740b_config_caen_recordlength"]);
    fV1751RecordLength     = this->castToSizeT_(fConfigValues["v1751_config_caen_recordlength"]);
    fV1740SampleReduction  = this->castToSizeT_(fConfigValues["v1740_config_caen_v1740_samplereduction"]);
    fV1740BSampleReduction = this->castToSizeT_(fConfigValues["v1740b_config_caen_v1740_samplereduction"]);
    fTDCPipelineDelay      = this->castToSizeT_(fConfigValues["tdc_config_tdc_pipelinedelay"]);
    fTDCGateWidth          = this->castToSizeT_(fConfigValues["tdc_config_tdc_gatewidth"]);

    // V1495 delay
    fV1495Delay = static_cast <double> (fV1495DelayTicks) * 0.01;

    // CAEN V1740 and V1751 documentations:
    //   https://cdcvs.fnal.gov/redmine/projects/lariat-online/wiki/LArIAT_DAQ_Commercial_Products
    // TDC documentation:
    //   https://cdcvs.fnal.gov/redmine/projects/lariat-online/wiki/TDC_Readout_Documentation

    // sampling rate in MHz
    fV1740SamplingRate  = 62.5 / fV1740SampleReduction;
    fV1740BSamplingRate = 62.5 / fV1740BSampleReduction;
    fV1751SamplingRate  = 1e3;

    // readout window length in microseconds
    fV1740ReadoutWindow  = fV1740RecordLength  / fV1740SamplingRate;
    fV1740BReadoutWindow = fV1740BRecordLength / fV1740BSamplingRate;
    fV1751ReadoutWindow  = fV1751RecordLength  / fV1751SamplingRate;
    fTDCReadoutWindow    = fTDCGateWidth * 8 * 0.001177;

    // pre-/post-acquisition and acquisition windows
    if (fV1740PreAcquisitionWindow   < 0) fV1740PreAcquisitionWindow   = fV1740ReadoutWindow;
    if (fV1740PostAcquisitionWindow  < 0) fV1740PostAcquisitionWindow  = 0.128;
    if (fV1740AcquisitionWindow      < 0) fV1740AcquisitionWindow      = fV1740ReadoutWindow;
    if (fV1740BPreAcquisitionWindow  < 0) fV1740BPreAcquisitionWindow  = fV1740BReadoutWindow;
    if (fV1740BPostAcquisitionWindow < 0) fV1740BPostAcquisitionWindow = 0.128;
    if (fV1740BAcquisitionWindow     < 0) fV1740BAcquisitionWindow     = fV1740BReadoutWindow;
    if (fV1751PreAcquisitionWindow   < 0) fV1751PreAcquisitionWindow   = 0.128;
    if (fV1751PostAcquisitionWindow  < 0) fV1751PostAcquisitionWindow  = 0.128;
    if (fV1751AcquisitionWindow      < 0) fV1751AcquisitionWindow      = fV1751ReadoutWindow;
    if (fTDCPreAcquisitionWindow     < 0) fTDCPreAcquisitionWindow     = fTDCReadoutWindow;
    if (fTDCPostAcquisitionWindow    < 0) fTDCPostAcquisitionWindow    = 0;
    if (fTDCAcquisitionWindow        < 0) fTDCAcquisitionWindow        = fTDCReadoutWindow;

    LOG_VERBATIM("SlicerInput")
    << "//////////////////////////////////////////////"
    << "\nV1495DelayTicks:             " << fV1495DelayTicks
    << "\nV1495Delay:                  " << fV1495Delay
    << "\nV1740PreAcquisitionWindow:   " << fV1740PreAcquisitionWindow
    << "\nV1740PostAcquisitionWindow:  " << fV1740PostAcquisitionWindow
    << "\nV1740AcquisitionWindow:      " << fV1740AcquisitionWindow
    << "\nV1740BPreAcquisitionWindow:  " << fV1740BPreAcquisitionWindow
    << "\nV1740BPostAcquisitionWindow: " << fV1740BPostAcquisitionWindow
    << "\nV1740BAcquisitionWindow:     " << fV1740BAcquisitionWindow
    << "\nV1751PreAcquisitionWindow:   " << fV1751PreAcquisitionWindow
    << "\nV1751PostAcquisitionWindow:  " << fV1751PostAcquisitionWindow
    << "\nV1751AcquisitionWindow:      " << fV1751AcquisitionWindow
    << "\nTDCPreAcquisitionWindow:     " << fTDCPreAcquisitionWindow
    << "\nTDCPostAcquisitionWindow:    " << fTDCPostAcquisitionWindow
    << "\nTDCAcquisitionWindow:        " << fTDCAcquisitionWindow
    << "//////////////////////////////////////////////";

    return;
  }

  //-----------------------------------------------------------------------
  size_t Slicer::castToSizeT_(std::string const& String)
  {
    size_t SizeT;

    if (!String.empty()) {
        SizeT = static_cast <size_t> (std::stoi(String));
    }
    else {
      SizeT = 0;
    }

    return SizeT;
  }

  //-----------------------------------------------------------------------
  void Slicer::commenceRun(art::RunPrincipal * & outRun)  // wtf is this? idk.
  {
    
    fFragmentToDigitAlg.InitializeRun(outRun, fRunNumber, fTimestamp);

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    // std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    sumdata::RunData runcol = sumdata::RunData(geo->DetectorName());
    art::put_product_in_principal(std::make_unique<sumdata::RunData>(runcol),
                                  * outRun,
                                  fSourceName);

    return;
  }
  
  //-----------------------------------------------------------------------
  void Slicer::commenceSubRun(art::SubRunPrincipal * & outSubRun)
  {
//    std::initializer_list<std::string> sam_params = {
//      "detector.cathode_voltage",
//      "detector.collection_voltage",
//      "detector.induction_voltage",
//      "detector.pmt_etl",
//      "detector.pmt_ham",
//      "detector.shield_voltage",
//      "detector.sipm_ham",
//      "detector.sipm_sensl",
//      "secondary.intensity",
//      "secondary.momentum",
//      "secondary.polarity",
//      "tertiary.beam_counters",
//      "tertiary.cherenkov1",
//      "tertiary.cherenkov2",
//      "tertiary.cosmic_counters",
//      "tertiary.DSTOF",
//      "tertiary.USTOF",
//      "tertiary.halo_paddle",
//      "tertiary.magnet_current",
//      "tertiary.magnet_polarity",
//      "tertiary.muon_range_stack",
//      "tertiary.MWPC1",
//      "tertiary.MWPC2",
//      "tertiary.MWPC3",
//      "tertiary.MWPC4",
//      "tertiary.number_MuRS",
//      "tertiary.punch_through",
//      "file_format"
//    };
    
    std::initializer_list<std::string> run_params = {
      "v1751_config_caen_enablereadout",
      "larasic_config_larasic_collection_filter",
      "larasic_config_larasic_collection_gain",
      "larasic_config_larasic_enablereadout",
      "larasic_config_larasic_pulseron",
      "larasic_config_larasic_channelscan",
      "v1740_config_caen_recordlength",
    };

    std::initializer_list<std::string> subrun_params = {
      "end_f_mc7sc1"
    };

    std::vector<std::string> runparams(run_params);
    std::vector<std::string> subrunparams(subrun_params);
    
    // get the relevant information from the database
    auto runValues    = fDatabaseUtility->GetConfigValues(runparams,
                                                          static_cast <int> (fRunNumber));
    auto subrunValues = fDatabaseUtility->GetIFBeamValues(subrunparams,
                                                          static_cast <int> (fRunNumber),
                                                          static_cast <int> (fSubRunNumber));

//    for(auto itr : runValues) 
//      LOG_VERBATIM("SlicerInput") << itr.first << " " << itr.second;
//
//    for(auto itr : subrunValues) 
//      LOG_VERBATIM("SlicerInput") << itr.first << " " << itr.second;

    // create the ConditionsSummary object and put it into the subrun
    // for the time being, several parameters are accessible from the
    // SAM database, and we don't yet know how to get those, so comment them out
    
//    std::vector<bool> mwpc(4);
//    mwpc[0] = (strcmp(runValues["tertiary.MWPC1"].c_str(), "on") == 0) ? true : false;
//    mwpc[1] = (strcmp(runValues["tertiary.MWPC1"].c_str(), "on") == 0) ? true : false;
//    mwpc[2] = (strcmp(runValues["tertiary.MWPC1"].c_str(), "on") == 0) ? true : false;
//    mwpc[3] = (strcmp(runValues["tertiary.MWPC1"].c_str(), "on") == 0) ? true : false;
//    
//    size_t              secondaryIntensity     = this->castToSizeT_(runValues["secondary.intensity"]);
//    size_t              secondaryMomentum      = this->castToSizeT_(runValues["secondary.momentum"]);
//    size_t              secondaryPolarity      = (strcmp(runValues["secondary.polarity"].c_str(), "Negative") == 0) ? 0 : 1;
//    size_t              magnetCurrent          = this->castToSizeT_(runValues["tertiary.magnet_current"]);
//    size_t              magnetPolarity         = (strcmp(runValues["tertiary.magnet_polarity"].c_str(), "Negative") == 0) ? 0 : 1;
//    size_t              tpcCathodeHV           = this->castToSizeT_(runValues["detector.cathode_voltage"]);
//    size_t              tpcCollectionV         = this->castToSizeT_(runValues["detector.collection_voltage"]);
//    size_t              tpcInductionV          = this->castToSizeT_(runValues["detector.induction_voltage"]);
//    size_t              tpcShieldV             = this->castToSizeT_(runValues["detector.shield_voltage"]);
//    size_t              etlPMTHV               = 0;
//    size_t              hamamatsuPMTHV         = 0;
//    size_t              hamamatsuSiPMHV        = 0;
//    size_t              senslSiPMHV            = 0;
//    size_t              tertiaryBeamCounters   = 0;
//    size_t              tertiaryCherenkov1     = 0;
//    size_t              tertiaryCherenkov2     = 0;
//    size_t              tertiaryCosmicCounters = 0;
//    size_t              dsTOF                  = 0;
//    size_t              usTOF                  = 0;
//    size_t              haloPaddle             = 0;
//    size_t              muonRangeStack         = 0;
//    size_t              numberMuRS             = this->castToSizeT_(runValues["tertiary.number_MuRS"]);
//    size_t              punchThrough           = 0;
//    bool                correctFileFormat      = false;

    LOG_VERBATIM("SlicerInput") << subrunValues["end_f_mc7sc1"];

    size_t              endMC7SC1              = this->castToSizeT_(subrunValues["end_f_mc7sc1"]);
    bool                v1751CaenEnableReadout = this->castToSizeT_(runValues["v1751_config_caen_enablereadout"]);
    size_t              asicCollectionFilter   = this->castToSizeT_(runValues["larasic_config_larasic_collection_filter"]);
    size_t              asicCollectionGain     = this->castToSizeT_(runValues["larasic_config_larasic_collection_gain"]);
    bool                asicEnableReadout      = this->castToSizeT_(runValues["larasic_config_larasic_enablereadout"]);
    bool                asicPulserOn           = this->castToSizeT_(runValues["larasic_config_larasic_pulseron"]);
    bool                asicChannelScan        = this->castToSizeT_(runValues["larasic_config_larasic_channelscan"]);
    size_t              v1740RecordLength      = this->castToSizeT_(runValues["v1740_config_caen_recordlength"]);
    
    

    std::unique_ptr<ldp::ConditionsSummary> conditionsSummary(new ldp::ConditionsSummary(/*secondaryIntensity,
                                                                                         secondaryMomentum,
                                                                                         secondaryPolarity,
                                                                                         magnetCurrent,
                                                                                         magnetPolarity,
                                                                                         tpcCathodeHV,
                                                                                         tpcCollectionV,
                                                                                         tpcInductionV,
                                                                                         tpcShieldV,
                                                                                         etlPMTHV,
                                                                                         hamamatsuPMTHV,
                                                                                         hamamatsuSiPMHV,
                                                                                         senslSiPMHV,
                                                                                         tertiaryBeamCounters,
                                                                                         tertiaryCherenkov1,
                                                                                         tertiaryCherenkov2,
                                                                                         tertiaryCosmicCounters,
                                                                                         dsTOF,
                                                                                         usTOF,
                                                                                         haloPaddle,
                                                                                         muonRangeStack,
                                                                                         numberMuRS,
                                                                                         punchThrough,,
                                                                                         correctFileFormat,
                                                                                         mwpc,*/
                                                                                         endMC7SC1,
                                                                                         v1751CaenEnableReadout,
                                                                                         asicCollectionFilter,
                                                                                         asicCollectionGain,
                                                                                         asicEnableReadout,
                                                                                         asicPulserOn,
                                                                                         asicChannelScan,
                                                                                         v1740RecordLength)
                                                              );
    art::put_product_in_principal(std::move(conditionsSummary), *outSubRun, fSourceName);

    return;
  }


  DEFINE_ART_INPUT_SOURCE(art::Source<rdu::Slicer>)

} // namespace rdu

#endif // SlicerInput_source
