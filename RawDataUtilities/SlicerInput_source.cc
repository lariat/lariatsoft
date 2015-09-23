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
#include "art/Framework/Services/Optional/TFileService.h"
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

    void commenceRun(art::RunPrincipal * & outRun);

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

    prhelper.reconstitutes<std::vector<raw::AuxDetDigit>, art::InEvent>(fSourceName);
    prhelper.reconstitutes<std::vector<raw::RawDigit>,    art::InEvent>(fSourceName);
    prhelper.reconstitutes<std::vector<raw::OpDetPulse>,  art::InEvent>(fSourceName);
    prhelper.reconstitutes<std::vector<raw::Trigger>,     art::InEvent>(fSourceName);
    prhelper.reconstitutes<sumdata::RunData,              art::InRun  >(fSourceName);

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

    std::cout << "\n////////////////////////////////////" << std::endl;
    std::cout << "fRunNumber:       " << fRunNumber       << std::endl;
    std::cout << "fSubRunNumber:    " << fSubRunNumber    << std::endl;
    std::cout << "fCachedRunNumber: " << fCachedRunNumber << std::endl;
    std::cout << "////////////////////////////////////\n" << std::endl;

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

    art::Timestamp timestamp; // LBNE should decide how to initialize this

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
    }

    std::cout << "fCollections.size(): " << fCollections.size() << std::endl;

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

    std::cout << "fCollectionIndex: " << fCollectionIndex << std::endl;
    std::cout << "NumberTPCReadouts: " << NumberTPCReadouts << std::endl;
    std::cout << "Collection.numberTPCReadouts: " << Collection.numberTPCReadouts << std::endl;
    std::cout << "Collection.caenBlocks.size(): " << Collection.caenBlocks.size() << std::endl;
    std::cout << "Collection.tdcBlocks.size(): " << Collection.tdcBlocks.size() << std::endl;

    std::cout << "Adding digits to run " << fRunNumber
              << ", sub-run "            << fSubRunNumber
              <<  ", event "             << fEventNumber
              << "..."                   << std::endl;

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

    std::cout << "//////////////////////////////////////////////"                                       << std::endl;
    std::cout << "V1495DelayTicks:       " << fConfigValues["v1495_config_v1495_delay_ticks"]           << std::endl;
    std::cout << "V1740PostPercent:      " << fConfigValues["v1740_config_caen_postpercent"]            << std::endl;
    std::cout << "V1740BPostPercent:     " << fConfigValues["v1740b_config_caen_postpercent"]           << std::endl;
    std::cout << "V1751PostPercent:      " << fConfigValues["v1751_config_caen_postpercent"]            << std::endl;
    std::cout << "V1740RecordLength:     " << fConfigValues["v1740_config_caen_recordlength"]           << std::endl;
    std::cout << "V1740BRecordLength:    " << fConfigValues["v1740b_config_caen_recordlength"]          << std::endl;
    std::cout << "V1751RecordLength:     " << fConfigValues["v1751_config_caen_recordlength"]           << std::endl;
    std::cout << "V1740SampleReduction:  " << fConfigValues["v1740_config_caen_v1740_samplereduction"]  << std::endl;
    std::cout << "V1740BSampleReduction: " << fConfigValues["v1740b_config_caen_v1740_samplereduction"] << std::endl;
    std::cout << "fTDCPipelineDelay:     " << fConfigValues["tdc_config_tdc_pipelinedelay"]             << std::endl;
    std::cout << "fTDCGateWidth:         " << fConfigValues["tdc_config_tdc_gatewidth"]                 << std::endl;
    std::cout << "//////////////////////////////////////////////"                                       << std::endl;

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

    std::cout << "//////////////////////////////////////////////"                << std::endl;
    std::cout << "V1495DelayTicks:             " << fV1495DelayTicks             << std::endl;
    std::cout << "V1495Delay:                  " << fV1495Delay                  << std::endl;
    std::cout << "V1740PreAcquisitionWindow:   " << fV1740PreAcquisitionWindow   << std::endl;
    std::cout << "V1740PostAcquisitionWindow:  " << fV1740PostAcquisitionWindow  << std::endl;
    std::cout << "V1740AcquisitionWindow:      " << fV1740AcquisitionWindow      << std::endl;
    std::cout << "V1740BPreAcquisitionWindow:  " << fV1740BPreAcquisitionWindow  << std::endl;
    std::cout << "V1740BPostAcquisitionWindow: " << fV1740BPostAcquisitionWindow << std::endl;
    std::cout << "V1740BAcquisitionWindow:     " << fV1740BAcquisitionWindow     << std::endl;
    std::cout << "V1751PreAcquisitionWindow:   " << fV1751PreAcquisitionWindow   << std::endl;
    std::cout << "V1751PostAcquisitionWindow:  " << fV1751PostAcquisitionWindow  << std::endl;
    std::cout << "V1751AcquisitionWindow:      " << fV1751AcquisitionWindow      << std::endl;
    std::cout << "TDCPreAcquisitionWindow:     " << fTDCPreAcquisitionWindow     << std::endl;
    std::cout << "TDCPostAcquisitionWindow:    " << fTDCPostAcquisitionWindow    << std::endl;
    std::cout << "TDCAcquisitionWindow:        " << fTDCAcquisitionWindow        << std::endl;
    std::cout << "//////////////////////////////////////////////"                << std::endl;

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
    fFragmentToDigitAlg.InitializeRun(fRunNumber);

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    // std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    sumdata::RunData runcol = sumdata::RunData(geo->DetectorName());
    art::put_product_in_principal(std::make_unique<sumdata::RunData>(runcol),
                                  * outRun,
                                  fSourceName);

    return;
  }

  DEFINE_ART_INPUT_SOURCE(art::Source<rdu::Slicer>)

} // namespace rdu

#endif // SlicerInput_source
