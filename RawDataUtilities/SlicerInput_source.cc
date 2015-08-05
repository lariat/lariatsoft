// art
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Root/rootNames.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// artdaq
#include "artdaq-core/Data/Fragments.hh"
#include "artdaq-core/Data/Fragment.hh"

// lardata
#include "RawData/RawDigit.h"

// C++
#include <functional>
#include <memory>
#include <regex>
#include <string>
#include <vector>
#include <map>

// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"

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

#include "FragmentToDigitAlg.h"

using raw::RawDigit;
using std::vector;
using std::string;

typedef std::map< int, std::map< int, std::vector< std::map< unsigned int, std::vector<unsigned int> > > > > match_maps;
typedef std::map< int, std::map< int, std::vector< std::pair<double, double> > > > fit_params_maps;
typedef std::vector<TDCFragment::TdcEventData> TDCDataBlock; 


//==========================================================================

namespace {


  // Retrieves branch name (a la art convention) where object resides
  template <typename PROD>
  const char* getBranchName( art::InputTag const & tag )
  {
    std::ostringstream pat_s;
    pat_s << art::TypeID(typeid(PROD)).friendlyClassName()
          << '_'
          << tag.label()
          << '_'
          << tag.instance()
          << '_'
          << tag.process()
          << ".obj";
    return pat_s.str().data();
  }

  artdaq::Fragments*
  getFragments( TBranch* br, unsigned entry )
  {
    br->GetEntry( entry );
    return reinterpret_cast<artdaq::Fragments*>( br->GetAddress() );
  }

  /*  
  art::EventAuxiliary*
  getEvtAuxiliary( TBranch * br, unsigned entry )
  {
    std::cout << "Pre GetEntry" << std::endl;
    br->GetEntry( entry );
    std::cout << "Post GetEntry" << std::endl;
    return reinterpret_cast<art::EventAuxiliary*>( br->GetAddress() );
  }
  */
  /*
  void getEventTimes( TTree* evtree, int & timehigh, int & timelow )
  {
    evtree->SetMakeClass(1);
    //    TBranch * branch = evtree->GetBranch("EventAuxiliary.time_.timeHigh_");
    // branch->SetAddress(&timehigh);
    //TBranch * branch1 = evtree->GetBranch("EventAuxiliary.time_.timeLow_");
    //branch->Print();
    //branch1->Print();
    //  evtree->Print();
    //branch1->SetAddress(&timelow);
    evtree->SetBranchAddress("EventAuxiliary.time_.timeHigh_",&timehigh);
    evtree->SetBranchAddress("EventAuxiliary.time_.timeLow_",&timelow);
  
    if( evtree->GetEntries() > 0 ){
      evtree->GetEntry( 0 );
    }
    else{ std::cout << "Event tree has no entries." << std::endl; }
    std::cout << "Event timeHigh: " << timehigh << ", timeLow: " << timelow << std::endl;
  }    
*/
  // Assumed file format is
  //
  //     "lbne_r[digits]_sr[digits]_other_stuff.root"
  //
  // This regex object is used to extract the run and subrun numbers.
  // The '()' groupings correspond to regex matches that can be
  // extracted using the std::regex_match facility.
  std::regex const filename_format( ".*\\lariat_r(\\d+)_sr(\\d+).*\\.root" );

}

//==========================================================================

namespace raw {

  // Enable 'pset.get<raw::Compress_t>("compression")'
  void decode (boost::any const & a, Compress_t & result ){
    unsigned tmp;
    fhicl::detail::decode(a,tmp);
    result = static_cast<Compress_t>(tmp);
  }
}

//==========================================================================

namespace DAQToOffline
{
  // The class SlicerInput is to be used as the template parameter for
  // art::Source<T>. It understands how to read art's ROOT data files,
  // to extract artdaq::Fragments from them, how to convert the
  // artdaq::Fragments into vectors of raw::RawDigit objects, and
  // finally how to re-combine those raw::RawDigit objects to present
  // the user of the art::Source<SlicerInput> with a sequence of
  // art::Events that have the desired event structure, different from
  // the event structure of the data file(s) being read.

  class SlicerInput
  {
  public:
    SlicerInput(fhicl::ParameterSet const& ps,
		art::ProductRegistryHelper& prh,
		art::SourceHelper& sh);

    // See art/Framework/IO/Sources/Source.h for a description of each
    // of the public member functions of SlicerInput.
    bool readFile(string const& filename, art::FileBlock*& fb);

    bool readNext(art::RunPrincipal* const& inR,
                  art::SubRunPrincipal* const& inSR,
                  art::RunPrincipal*& outR,
                  art::SubRunPrincipal*& outSR,
                  art::EventPrincipal*& outE);

    void closeCurrentFile();

  private:

    using rawDigits_t = vector<RawDigit>;


    art::ServiceHandle<art::TFileService>      tfs;                      ///< handle to the TFileService
    string                                     fSourceName;
    string                                     fLastFileName;
    std::unique_ptr<TFile>                     fFile;
    bool                                       fDoneWithFiles;
    art::InputTag                              fInputTag;
    art::SourceHelper                          fSH;
    TBranch*                                   fFragmentsBranch;
    TBranch*                                   fEvtAuxBranch;
    size_t                                     fNInputEvts;
    size_t                                     fTreeIndex;
    art::RunNumber_t                           fRunNumber;
    art::SubRunNumber_t                        fSubRunNumber;
    art::EventNumber_t                         fEventNumber;
    art::RunNumber_t                           fCachedRunNumber;
    art::SubRunNumber_t                        fCachedSubRunNumber;
    size_t                                     fMaxNumberFitIterations;  ///< number of fit iterations before stopping
    bool                                       fVerbose;

    LariatFragment                           fLArIATFragment; 
    std::map<int,std::vector<CAENFragment> > fTriggerToCAENDataBlocks;
    std::map<int,std::vector<TDCDataBlock> > fTriggerToTDCDataBlocks;
    bool                                     fDoneWithThisFile;
    FragmentToDigitAlg                       fFrag2DigAlg;
    size_t                                   fTriggerDecisionTick;     ///< tick at which to expect the trigger decision
    float                                    fTrigger1740Pedestal;     ///< pedestal value for the 1740 readout of the triggers
    float                                    fTrigger1740Threshold;    ///< 1740 readout must go below the pedestal this much to trigger


    bool miniFragmentUtility();
    void makeEventAndPutFragments( art::EventPrincipal*& outE );

    //////// MATCHING ALG STUFF ////////////
    void matchDataBlocks(const LariatFragment * data);
    
    void coarseMatch(int   const& deviceAID,
		     int   const& deviceBID,
		     double       range[2],
		     std::map< int, std::map<unsigned int, double> > timeStamps,
		     match_maps & matchMaps);
    
    void fineMatch(int      const& deviceAID,
		   int      const& deviceBID,
		   double          range[2],    
		   fit_params_maps fitParametersMaps,
		   std::map< int, std::map<unsigned int, double> > timeStamps,
		   match_maps    & matchMaps);
    
    void printMatchMap(int   const& deviceAID,
		       int   const& deviceBID,
		       match_maps & matchMaps);
    
    double line(std::pair<double, double> const& parameters, 
		double                    const& x);
    
    double clockDriftCorr(std::pair<double, double> const& parameters,
			  double                    const& x);
    
    void fitClockDrift(int         const& deviceAID,
		       int         const& deviceBID,
		       std::map< int, std::map<unsigned int, double> > timeStamps,
		       match_maps       & matchMaps,
		       fit_params_maps  & fitParametersMaps,
		       std::string const& graphNamePrefix);
    
    void matchFitIter(int        const& deviceAID,
		      int        const& deviceBID,
		      std::map< int, std::map<unsigned int, double> > timeStamps,
		      match_maps      & matchMaps,
		      fit_params_maps & fitParametersMaps,
		      double            coarseRange[2],
		      double            fineEps[2]);
    ////////////////////////////////////////

    void commenceRun( art::RunPrincipal*& outR );
    std::vector<raw::Trigger> makeTheTriggers( std::vector<CAENFragment> theCAENDataBlocks,
					       std::vector<TDCDataBlock> theTDCDataBlocks );
    uint32_t triggerBits               (std::vector<CAENFragment>     const& caenFrags);   				  

  };
}

//=======================================================================================
DAQToOffline::SlicerInput::SlicerInput(fhicl::ParameterSet const& ps,
                                 art::ProductRegistryHelper& prh,
                                 art::SourceHelper& sh) :
  fSourceName("SlicerInput"),
  fLastFileName(ps.get<vector<string>>("fileNames",{}).back()), //REL maybe for identifying when subrun files are done being split?
  fFile(),
  fDoneWithFiles(false),
  fInputTag("daq:SPILL:DAQ"), // "moduleLabel:instance:processName"
  fSH(sh),
  fFragmentsBranch(nullptr),                                    //REL not sure if to branch of raw data keeping fragments, or something else...
  fNInputEvts(),    
  fRunNumber(1),      // Defaults in case input filename does not
  fSubRunNumber(0),   // follow assumed filename_format above.
  fEventNumber(),
  fCachedRunNumber(-1),
  fCachedSubRunNumber(-1),
  //  fBufferedDigits(),
  fMaxNumberFitIterations(ps.get<size_t>("maxNumberFitIterations")),
  fVerbose(true),
  fFrag2DigAlg(ps.get< fhicl::ParameterSet >("FragmentToDigitAlg")),
  fTriggerDecisionTick(size_t(137)),
  fTrigger1740Pedestal(2000),
  fTrigger1740Threshold(0)

{
  // Will use same instance name for the outgoing products as for the
  // incoming ones.
  prh.reconstitutes<std::vector<raw::AuxDetDigit>,art::InEvent>( fSourceName );
  prh.reconstitutes<std::vector<raw::RawDigit>,art::InEvent>( fSourceName );
  prh.reconstitutes<std::vector<raw::OpDetPulse>,art::InEvent>( fSourceName );
  prh.reconstitutes<sumdata::RunData,art::InRun>( fSourceName );
  prh.reconstitutes<std::vector<raw::Trigger>,art::InEvent>( fSourceName );
}

//=======================================================================================
bool
DAQToOffline::SlicerInput::readFile(string const& filename, art::FileBlock*& fb)
{
  // Run numbers determined based on file name...see comment in
  // anon. namespace above.
  std::smatch matches;
  if ( std::regex_match( filename, matches, filename_format ) ) {
    fRunNumber       = std::stoul( matches[1] );
    fSubRunNumber    = std::stoul( matches[2] );
  }
  
  // Get fragments branch
  fFile.reset( new TFile(filename.data()) );
  TTree* evtree    = reinterpret_cast<TTree*>(fFile->Get(art::rootNames::eventTreeName().c_str()));
  fFragmentsBranch = evtree->GetBranch( getBranchName<artdaq::Fragments>( fInputTag ) ); // get branch for specific input tag
  //  std::cout << "Pre GetBranch" << std::endl;
  //  fEvtAuxBranch    = evtree->GetBranch( "EventAuxiliary" );
  fNInputEvts      = static_cast<size_t>( fFragmentsBranch->GetEntries() );              //Number of fragment-containing events to read in from input file 
  fTreeIndex       = 0ul;
  
  // New fileblock
  fb = new art::FileBlock(art::FileFormatVersion(),filename);
  if ( fb == nullptr ) {
    throw art::Exception(art::errors::FileOpenError)
      << "Unable to open file " << filename << ".\n";
  }
  
  //Resetting transient variables
  fDoneWithThisFile = false;

  /*  
  //LEAVE THIS FOR NOW - IT'S FOR GETTING CORRECT EVENT TIME FROM ROOT FILE - STILL UNDER CONSTRUCTION
  size_t nEvtAuxEvents = fEvtAuxBranch->GetEntries();
  std::cout << "Number of entries: " << nEvtAuxEvents << std::endl;
  art::EventAuxiliary* evtAux = getEvtAuxiliary( fEvtAuxBranch, 0 );
  std::cout << "Pre search time" << std::endl;
  std::cout << "EventAuxiliary time: " << evtAux->time().timeLow() << std::endl;
  std::cout << "EventAuxiliary time High: " << evtAux->time().timeHigh() << std::endl;
  std::cout << "EventAuxiliary run: " << evtAux->id().run() << std::endl;

  //Getting event time
  //  int timehigh = 0 ;
  //  int timelow = 0;
  //  getEventTimes( evtree, timehigh, timelow );
  */
  
  //Stores processed LArIATFragment in member variable
  bool fragProcessedOkay = miniFragmentUtility();
  if( !fragProcessedOkay ) std::cout << "Fragment not processed well." << std::endl;
					
  //Accesses member variables and creates maps of ints to fragment vects
  //(matching data blocks with triggers)
  matchDataBlocks( &fLArIATFragment );
					
					

  return true;
}

//=======================================================================================
bool
DAQToOffline::SlicerInput::readNext(art::RunPrincipal*    const& inR,
                                 art::SubRunPrincipal* const& inSR,
                                 art::RunPrincipal*    & outR,
                                 art::SubRunPrincipal* & outSR,
                                 art::EventPrincipal*  & outE)
{
  if ( fDoneWithThisFile ) {
    return false;
  }

  //Check to see if we're on our last file
  std::cout << "fFile->GetName: " << fFile->GetName() << ", fLastFileName: " << fLastFileName << std::endl;
  // if( fFile->GetName() == fLastFileName )
  //  fDoneWithFiles = true;

  

  art::Timestamp ts; // LBNE should decide how to initialize this
  //REL if this is a new run number, make a new run principal - This should only call once for a given file
  if ( fRunNumber != fCachedRunNumber ){
    outR = fSH.makeRunPrincipal(fRunNumber,ts);
    fCachedRunNumber = fRunNumber;
    fEventNumber = 0ul;
    
    //Add run information to run and initialize some
    //trigger-to-digit algorithm stuff
    commenceRun(outR);
    
  }

  std::cout << "Run: " << fRunNumber << std::endl;

  //REL if this is a new subrun number, make a new subrun principal
  if ( fSubRunNumber != fCachedSubRunNumber ) {
    outSR = fSH.makeSubRunPrincipal(fRunNumber,fSubRunNumber,ts);
    fCachedSubRunNumber = fSubRunNumber;
    fEventNumber = 0ul;
  }

  std::cout << "Subrun: " << fSubRunNumber << std::endl;
  

  //Pseudocode
  // - Convert the artdaq::fragments to artdaq::fragment to LArIATFragment using a miniature version of FragmentUtility
  // - Group the lariat fragments together into triggers using johnny's algorithm
  // - for each trigger, make event and put fragments

  if( fVerbose )
    std::cout << "Number of CAEN Data Blocks, pre-matching: " << (&fLArIATFragment)->caenFrags.size() << ", number of TDC Data Blocks, pre-matching: " << (&fLArIATFragment)->tdcFrags.size() << std::endl;
  
  //Make event corresponding to one entry in the two maps
  makeEventAndPutFragments( outE ); 

  return true;

  /*
  // eventIsFull_ is what LBNE should modify based on its needs
  while ( !eventIsFull_( fBufferedDigits ) ) {
    if ( loadedDigits_.empty() && !loadDigits_() ) {
      if ( fFile->GetName() != fLastFileName )
        return false;
      else {
        fDoneWithFiles = true;
        break;
      }
    }
    fBufferedDigits.emplace_back( loadedDigits_.next() );
  }
  makeEventAndPutDigits_( outE );
  */

  //Temp
  //  return true;

}


//=======================================================================================
void
DAQToOffline::SlicerInput::closeCurrentFile()
{ 
  fFile.reset(nullptr);
}

//=======================================================================================
bool DAQToOffline::SlicerInput::miniFragmentUtility()
{
  if( fVerbose ){
    std::cout << "Number of input events: " << fNInputEvts << std::endl;
    std::cout << "Value of treeIndex: " << fTreeIndex << std::endl;
  }
  if( fTreeIndex != fNInputEvts ){
    //Read in the artdaq fragments from the root file (should only be one fragment per file)
    artdaq::Fragments* fragments = getFragments( fFragmentsBranch, fTreeIndex++ );
    if( (*fragments).size() > 1 )
      throw cet::exception("SlicerInput") << "artdaq::Fragment vector contains more than one fragment. Bail.";
    
    //Convert the artdaq fragment to a lariat fragment
    artdaq::Fragment frag = fragments->at(0);
    const char* bytePtr = reinterpret_cast<const char*> (&*frag.dataBegin());
    LariatFragment * LArIATFragment = new LariatFragment((char *) bytePtr, frag.dataSize() * sizeof(unsigned long long));
    fLArIATFragment = *LArIATFragment;
    if( fVerbose )
      std::cout << "MiniFragmentUtility executed well." << std::endl;
    return true;
  }
  return false;
}

//=======================================================================================
void
DAQToOffline::SlicerInput::makeEventAndPutFragments(art::EventPrincipal*& outE){
  ++fEventNumber;
  if( fVerbose )
    std::cout << "Event number: " << fEventNumber << std::endl;
  outE = fSH.makeEventPrincipal( fRunNumber, fSubRunNumber, fEventNumber, art::Timestamp() );
 
  //Debugging
  if( fVerbose )
    std::cout << "Number of triggers with CAEN Data Blocks: " << fTriggerToCAENDataBlocks.size() << ", number of triggers with TDC Data Blocks: " << fTriggerToTDCDataBlocks.size() << std::endl;

  //Find the smallest trigNum in the two maps
  int smallestTrigNum = 9999;
  std::map< int, std::vector<CAENFragment> >::iterator iterCAEN;
  std::map< int, std::vector<TDCDataBlock> >::iterator iterTDC;
  for( iterCAEN = fTriggerToCAENDataBlocks.begin() ; iterCAEN != fTriggerToCAENDataBlocks.end() ; ++iterCAEN )
    if( iterCAEN->first < smallestTrigNum ) smallestTrigNum = iterCAEN->first;
  for( iterTDC = fTriggerToTDCDataBlocks.begin() ; iterTDC != fTriggerToTDCDataBlocks.end() ; ++iterTDC )
    if( iterTDC->first < smallestTrigNum ) smallestTrigNum = iterTDC->first;

  //Retrieve the last vector of fragments in the CAEN and TDC maps
  //Make empty vectors in case there are events being created where there aren't data blocks for one or the other
  std::vector<CAENFragment> theCAENDataBlocks;
  std::vector<TDCDataBlock> theTDCDataBlocks;
  if( fTriggerToCAENDataBlocks.count(smallestTrigNum) > 0 ){
    theCAENDataBlocks = fTriggerToCAENDataBlocks.at(smallestTrigNum);
    fTriggerToCAENDataBlocks.erase(smallestTrigNum);
  }
  //Now do this for the TDC Datablocks (entering them into the event)
  if( fTriggerToTDCDataBlocks.count(smallestTrigNum) > 0 ){
    theTDCDataBlocks = fTriggerToTDCDataBlocks.at(smallestTrigNum);
    fTriggerToTDCDataBlocks.erase(smallestTrigNum);
  }
  

  /*
  //GOOD STUFF, BUT ACCIDENTALLY REVERSES THE EVENT ORDERING
  //Find the largest trigNum in the two maps
  int largestTrigNum = 0;
  std::map< int, std::vector<CAENFragment> >::iterator iterCAEN;
  std::map< int, std::vector<TDCDataBlock> >::iterator iterTDC;
  for( iterCAEN = fTriggerToCAENDataBlocks.begin() ; iterCAEN != fTriggerToCAENDataBlocks.end() ; ++iterCAEN )
    if( iterCAEN->first > largestTrigNum ) largestTrigNum = iterCAEN->first;
  for( iterTDC = fTriggerToTDCDataBlocks.begin() ; iterTDC != fTriggerToTDCDataBlocks.end() ; ++iterTDC )
    if( iterTDC->first > largestTrigNum ) largestTrigNum = iterTDC->first;
    
  //Retrieve the last vector of fragments in the CAEN and TDC maps
  //Make empty vectors in case there are events being created where there aren't data blocks for one or the other
  std::vector<CAENFragment> theCAENDataBlocks;
  std::vector<TDCDataBlock> theTDCDataBlocks;
  if( fTriggerToCAENDataBlocks.count(largestTrigNum) > 0 ){
    theCAENDataBlocks = fTriggerToCAENDataBlocks.at(largestTrigNum);
    fTriggerToCAENDataBlocks.erase(largestTrigNum);
  }
  //Now do this for the TDC Datablocks (entering them into the event)
  if( fTriggerToTDCDataBlocks.count(largestTrigNum) > 0 ){
    theTDCDataBlocks = fTriggerToTDCDataBlocks.at(largestTrigNum);
    fTriggerToTDCDataBlocks.erase(largestTrigNum);
  }
  */
  //Creating digits for insertion into the event
  std::vector<raw::AuxDetDigit>          auxDigits;	 
  std::vector<raw::RawDigit>    	 rawDigits;	 
  std::vector<raw::OpDetPulse>   	 opPulses;	 

  //Make triggers from the fragments
  
  std::vector<raw::Trigger> trigVect = makeTheTriggers( theCAENDataBlocks,
							theTDCDataBlocks );
				

  //Take fragments and convert them to the created digits
  fFrag2DigAlg.makeTheDigits( theCAENDataBlocks,
			      theTDCDataBlocks,
			      auxDigits,
			      rawDigits,
			      opPulses );

  //Now we have vectors of auxDigits, rawDigits, and opPulses that we feed into the event.
  art::put_product_in_principal( std::make_unique<std::vector<raw::Trigger> >(trigVect),
				 *outE,
				 fSourceName);
  art::put_product_in_principal( std::make_unique<std::vector<raw::AuxDetDigit> >(auxDigits),
				 *outE,
				 fSourceName);
  art::put_product_in_principal( std::make_unique<std::vector<raw::RawDigit> >(rawDigits),
				 *outE,
				 fSourceName);
  art::put_product_in_principal( std::make_unique<std::vector<raw::OpDetPulse> >(opPulses),
				 *outE,
				 fSourceName);
  
  //Check to see if we're done with this file yet
  if( fTriggerToCAENDataBlocks.size() == 0 && fTriggerToTDCDataBlocks.size() == 0 )
    fDoneWithThisFile = true;
  else{ fDoneWithThisFile = false; }
  
  mf::LogDebug("FragmentsTest") << "Producing event: " << outE->id();

  return;
  
}

//=======================================================================================
std::vector<raw::Trigger> DAQToOffline::SlicerInput::makeTheTriggers( std::vector<CAENFragment> theCAENDataBlocks,
								      std::vector<TDCDataBlock> theTDCDataBlocks )
{
  //Hardcoded for now until I can find out how to 
  //extract this from the raw event root file
  float eventTime = 0;
  
  bool caenDataPresent = false;		     
  std::vector<raw::Trigger> trigVect;

  //If there are caen fragments present, set the trigger info
  //based on the caen data block information
  if( theCAENDataBlocks.size() > 0 ){
    trigVect.push_back(raw::Trigger(fEventNumber, theCAENDataBlocks.front().header.triggerTimeTag, eventTime, this->triggerBits(theCAENDataBlocks) ) );
    caenDataPresent = true;  
  }
  else
    LOG_WARNING("FragmentToDigit") << "There are no CAEN Fragments for event " << fEventNumber
				   << " that may be OK, so continue";
  
  //If there are no caen fragments present, then set the trigger info
  //based on the TDCDataBlock information
  if( theTDCDataBlocks.size() > 0 ){    
    if(!caenDataPresent){
      trigVect.push_back(raw::Trigger(fEventNumber, theTDCDataBlocks.front().front().tdcEventHeader.tdcTimeStamp, 
				      eventTime, this->triggerBits(theCAENDataBlocks)));
    }
  }
  else
    LOG_WARNING("FragmentToDigit") << "There are no TDC Fragments for event " << fEventNumber
				   << " that may be OK, so continue";

  return trigVect;

}

//=======================================================================================
uint32_t DAQToOffline::SlicerInput::triggerBits(std::vector<CAENFragment> const& caenFrags)
{

  // the trigger bits are piped into the V1740 board in slot 7, inputs 48 to 63
  // after run 6154 the bits were piped into a V1740 in slot 24, inputs 48 to 63
  // these are example connections as of May 08, 2015
  // 0   WC1      | OR of 2 X view TDCs ANDed with OR of 2 Y
  // 1   WC2      | "                                      " 
  // 2   WC3      | "                                      " 
  // 3   WC4      | "                                      " 
  // 4   BEAMON   | Spill gate : STARTs on $21, STOPs on $36 (cable says $26 but Bill says $36)
  // 5   USTOF    | OR of 4 PMTs
  // 6   DSTOF    | OR of 2 PMTs
  // 7   PUNCH    | OR of 2 X view paddles ANDed with OR of 2 Y
  // 8   HALO     | OR of 2 PMTs
  // 9   PULSER   |
  // 10  COSMICON | Cosmic gate : STARTs on $36, STOPs on $00 (not optimal, would like to stop before $00)
  // 11  COSMIC   | the trigger signal from the cosmic rack
  // 12  PILEUP   | Coincidence of any later LARSCINT with a delayed gate initiated by itself. Higher discrimination thresh. 
  // 13  MICHEL   | Coincidence of two light flashes in TPC (LARSCINT) occurring within a 5us time window
  // 14  LARSCINT | Coincidence of Hamamatsu and ETL PMTs (discriminated)
  // 15  MuRS     | Any coincidence of two planes.  Each plane is the OR of the discriminated pulses of 4 paddles. 

  // Each waveform corresponds to a single trigger channel.  If the (pedestal subtracted?) value of any ADC
  // in a waveform is less than 0, then the trigger for that channel fired

  // Need database eventually to set this correctly for different data-taking periods.
  // would set the fTriggerDecisionTick, fTrigger1740Pedestal, fTrigger1740Threshold values
  // in the beginRun method

  std::bitset<16> triggerBits;

  size_t minChan  = 48;
  size_t maxChan  = 64;

  for(auto const& frag : caenFrags){

    if     (frag.header.boardId != 7  && fRunNumber < 6155) continue;
    else if(frag.header.boardId != 24 && fRunNumber > 6154) continue;

    for(size_t chan = minChan; chan < maxChan; ++chan){ 
      if(chan > frag.waveForms.size() )
	throw cet::exception("FragmentToDigit") << "attempting to access channel "
						<< chan << " from 1740 fragment with only "
						<< frag.waveForms.size() << " channels";
      
      // only look at the specific tick of the waveform where the trigger decision is taken
      if(frag.waveForms[chan].data.size() > fTriggerDecisionTick - 1)
	// the trigger waveform goes below the pedestal (low) if the trigger is on
	//Hard-coded temporarily to hit the right window of ticks for trigge (it's not just one tick)
	//Studies about a precise window for this first trigger tick are underway
	if(fTrigger1740Pedestal - frag.waveForms[chan].data[135] > fTrigger1740Threshold ||
	   fTrigger1740Pedestal - frag.waveForms[chan].data[136] > fTrigger1740Threshold ||
	   fTrigger1740Pedestal - frag.waveForms[chan].data[137] > fTrigger1740Threshold ||
	   fTrigger1740Pedestal - frag.waveForms[chan].data[138] > fTrigger1740Threshold ) 
	   triggerBits.set(chan - minChan);

    } // end loop over channels on the board
  } // end loop over caen fragments

  return triggerBits.to_ulong();
}  

//=======================================================================================
void DAQToOffline::SlicerInput::matchDataBlocks(const LariatFragment * data) 
{

  // maps for matching fragments
  std::map< int, std::map<unsigned int, double> > dataBlockTimeStamps;
  match_maps matchMaps;

  //////////////////////////////////////////////////////////////////////
  // Some notes
  //////////////////////////////////////////////////////////////////////
  //
  // There are multiple CAEN fragments and only one TDC fragment for
  // each spill. Each TDC fragment holds multiple TDC events.
  //
  // Since there are multiple CAEN fragments and only one TDC fragment
  // per spill, I will use the term "data block" to refer to a single
  // CAEN fragment or a single TDC event.
  //
  // dataBlockTimeStamps will store the data block time stamps from all
  // devices (except the WUT) to help with the matching.
  //
  // Here is some pseudocode on how to access the time stamps:
  //
  // for deviceID in devices:
  //     for dataBlockIndex in dataBlocks:
  //         dataBlockTimeStamp = dataBlockTimeStamps[deviceID][dataBlockIndex];
  //
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Device ID key
  //////////////////////////////////////////////////////////////////////
  //  Device ID  Device
  //  ---------  ------
  //  0          CAEN boardId 0, V1740 board 0
  //  1          CAEN boardId 1, V1740 board 1
  //  2          CAEN boardId 2, V1740 board 2
  //  3          CAEN boardId 3, V1740 board 3
  //  4          CAEN boardId 4, V1740 board 4
  //  5          CAEN boardId 5, V1740 board 5
  //  6          CAEN boardId 6, V1740 board 6
  //  7          CAEN boardId 7, V1740 board 7
  //  8          CAEN boardId 8, V1751 board 0
  //  9          CAEN boardId 9, V1751 board 1
  //  24         CAEN boardId 24, V1740 board 24
  //  32         Multi-wire proportional chambers, 16 TDCs
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //@\ BEGIN: loop through the data blocks to get their time stamps
  //          into dataBlockTimeStamps
  //////////////////////////////////////////////////////////////////////

  size_t numberCaenDataBlocks[32] = {};
  size_t numberMwpcDataBlocks = 0;

  const size_t numberCaenFrags = data->caenFrags.size();
  //REL LOG_VERBATIM("FragmentToDigit") << "Found " << numberCaenFrags << " CAEN fragments";

  if (numberCaenFrags > 0) 
    //REL LOG_VERBATIM("FragmentToDigit") << "Looking at CAEN fragments...";

  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment const& caenFrag = data->caenFrags[i];
    unsigned int boardId = static_cast <unsigned int> (caenFrag.header.boardId);
    unsigned int index = numberCaenDataBlocks[boardId];
    int deviceID = boardId;
    // each CAEN Trigger Time Tag count is 8 ns
    double timeStamp = caenFrag.header.triggerTimeTag * 0.008;  // convert to microseconds
    dataBlockTimeStamps[deviceID][index] = timeStamp;
    numberCaenDataBlocks[boardId] += 1;
  }

  const size_t numberTdcFrags = data->tdcFrags.size();
  //REL  LOG_VERBATIM("FragmentToDigit") << "Found " << numberTdcFrags << " TDC fragments";

  if (numberTdcFrags > 0) 
    //REL  LOG_VERBATIM("FragmentToDigit") << "Looking at TDC fragments...";

  for (size_t i = 0; i < numberTdcFrags; ++i) {

    TDCFragment const& tdcFrag = data->tdcFrags[i];

    std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcEvents = tdcFrag.tdcEvents;

    //LOG_DEBUG("FragmentToDigit")
    //    << "tdcEvents.size(): " << tdcEvents.size();
    //tdcFrag.print();

    numberMwpcDataBlocks = tdcEvents.size();

    for (size_t j = 0; j < tdcEvents.size(); ++j) {

      if (tdcFrag.controllerHeader.nTDCs != tdcEvents[j].size()) {
        mf::LogError("FragmentToDigit") << "*** Fatal nTDCs mismatch: " << tdcEvents[j].size()
            << " != " << tdcFrag.controllerHeader.nTDCs<< " "<< j;
        continue;
      }

      LOG_DEBUG("FragmentToDigit") << "TDC event: " << j;

      for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS; ++tdc_index) {
        TDCFragment::TdcEventData tdcEventData = tdcEvents[j].at(tdc_index);

        uint16_t controllerTimeStamp = tdcEventData.tdcEventHeader.controllerTimeStamp;
        uint32_t tdcTimeStamp = tdcEventData.tdcEventHeader.tdcTimeStamp;
        // each TDC Time Stamp count is 1/106.208 microseconds
        double timeStamp = tdcEventData.tdcEventHeader.tdcTimeStamp / 106.208;  // convert to microseconds

        int deviceID = 32;

        LOG_DEBUG("FragmentToDigit") << "  TDC index: " << tdc_index;
        LOG_DEBUG("FragmentToDigit") << "  TDC time stamp: " << tdcTimeStamp;
        LOG_DEBUG("FragmentToDigit") << "  Controller time stamp: " << controllerTimeStamp;

        dataBlockTimeStamps[deviceID][j] = timeStamp;

      }
    }
  }

  // let's forget the WUT for now
  //const int numberWutFrags = data->wutFrags.size();
  //LOG_DEBUG("FragmentToDigit")
  //    << "Found " << numberWutFrags << " WUT fragments";

  //////////////////////////////////////////////////////////////////////
  //@\ END: loop through the data blocks to get their time stamps
  //        into dataBlockTimeStamps
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //@\ BEGIN: loop through dataBlockTimeStamps to check that they are there
  //////////////////////////////////////////////////////////////////////

  for (size_t i = 0; i < 10; ++i) {
    int deviceID = i;
    unsigned int indices = numberCaenDataBlocks[i];
    LOG_DEBUG("FragmentToDigit") << "Board ID: " << i << ", number of data blocks: " << numberCaenDataBlocks[i];
    LOG_DEBUG("FragmentToDigit") << "Device ID: " << deviceID;
    for (size_t j = 0; j < indices; ++j) {
      LOG_DEBUG("FragmentToDigit") << "  Index: " << j;
      LOG_DEBUG("FragmentToDigit") << "  Time stamp: " << dataBlockTimeStamps[deviceID][j];
    }
  }

  for (size_t i = 0; i < TDCFragment::MAX_TDCS; ++i) {
    int deviceID = 32;
    LOG_DEBUG("FragmentToDigit") << "TDC index: " << i << ", number of data blocks: " << numberMwpcDataBlocks;
    LOG_DEBUG("FragmentToDigit") << "Device ID: " << deviceID;
    for (size_t j = 0; j < numberMwpcDataBlocks; ++j) {
      LOG_DEBUG("FragmentToDigit") << "  Index: " << j;
      LOG_DEBUG("FragmentToDigit") << "  Time stamp: " << dataBlockTimeStamps[deviceID][j];
    }
  }

  //////////////////////////////////////////////////////////////////////
  //@\ END: loop through dataBlockTimeStamps to check that they are there
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //@\ BEGIN: intra V1740/V1751 matching
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // NOTE: The master clocks of the V1740s is boardId 0; the master
  //       clock of the V1751s is boardId 8.
  //////////////////////////////////////////////////////////////////////

  // the difference in time stamps should be in between this range for a match; microseconds
  double v1740IntraRange[2] = { -0.032, 0.032 };
  double v1751IntraRange[2] = { -0.032, 0.032 };

  this->coarseMatch(0, 1, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 2, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 3, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 4, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 5, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 6, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(0, 7, v1740IntraRange, dataBlockTimeStamps, matchMaps);
  this->coarseMatch(8, 9, v1751IntraRange, dataBlockTimeStamps, matchMaps);

  //printMatchMap(0, 1, matchMaps);
  //printMatchMap(0, 2, matchMaps);
  //printMatchMap(0, 3, matchMaps);
  //printMatchMap(0, 4, matchMaps);
  //printMatchMap(0, 5, matchMaps);
  //printMatchMap(0, 6, matchMaps);
  //printMatchMap(0, 7, matchMaps);
  //printMatchMap(8, 9, matchMaps);

  //////////////////////////////////////////////////////////////////////
  //@\ END: intra V1740/V1751 matching
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //@\ BEGIN: matching
  //////////////////////////////////////////////////////////////////////

  // the difference in time stamps should be in between this range for a match; microseconds
  // these variables should probably be placed in the .fcl file
  double v1751v1740InterRange[2] = { -500, 500 };
  double v1751MwpcInterRange[2]  = {    0, 160 };
  double v1740MwpcInterRange[2]  = {    0, 160 };

  // acceptance range from fitted line to allow a match; microseconds
  // these variables should probably be placed in the .fcl file
  // y - eps[0] <= y_test <= y + eps[1]
  double v1751v1740InterEps[2] = { 1, 1 };
  double v1751MwpcInterEps[2]  = { 1, 1 };
  double v1740MwpcInterEps[2]  = { 1, 1 };

  fit_params_maps fitParamsMaps;

  this->matchFitIter(8,  0, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  1, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  2, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  3, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  4, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  5, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  6, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8,  7, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8, 24, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751v1740InterRange, v1751v1740InterEps);
  this->matchFitIter(8, 32, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1751MwpcInterRange,  v1751MwpcInterEps);
  this->matchFitIter(0, 32, dataBlockTimeStamps, matchMaps, fitParamsMaps, v1740MwpcInterRange,  v1740MwpcInterEps);

  //////////////////////////////////////////////////////////////////////
  //@\ END: matching
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  //@\ BEGIN: add to trigger ID to data block map
  //////////////////////////////////////////////////////////////////////

  int triggerID = 0;

  size_t numberMatchedCaenDataBlocks[32] = {};
  size_t numberMatchedMwpcDataBlocks = 0;

  fitParamsMaps[8][8].push_back(std::make_pair<double, double>(0, 0));
  fitParamsMaps[8][9].push_back(std::make_pair<double, double>(0, 0));

  // this shouldn't be hard-coded, but there is no way to get
  // the decimation factor of the sample rate from the CAENFragment
  double v1740SampleTime = 0.128;  // microseconds
  double v1751SampleTime = 0.001;  // microseconds
  size_t numberV1740Samples = 0;
  size_t numberV1751Samples = 0;

  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment const& caenFrag = data->caenFrags[i];
    unsigned int boardId = static_cast <unsigned int> (caenFrag.header.boardId);
    if (boardId == 0) {
      numberV1740Samples = static_cast <size_t> (caenFrag.header.nSamples);
    }
    if (boardId == 8) {
      numberV1751Samples = static_cast <size_t> (caenFrag.header.nSamples);
    }
    if (numberV1740Samples > 0 and numberV1751Samples > 0) break;
  }

  // vector of corrected timestamps to data block indices:
  //[
  //  timestamp0, [ caenBoardID, caenFragmentIndex0 ],
  //  timestamp1, [ caenBoardID, caenFragmentIndex1 ],
  //  timestamp3, [      mwpcID,     tdcEventIndex0 ],
  //  timestamp4, [ caenBoardID, caenFragmentIndex2 ],
  //  timestamp5, [      mwpcID,     tdcEventIndex1 ],
  //  ...
  //]
  std::vector< std::pair<double, std::pair<int, int> > > timeStampToDataBlockIndices;

  // correct CAENFragment timestamp and add corrected timestamp and index to vector
  for (size_t i = 0; i < numberCaenFrags; ++i) {
    CAENFragment const& caenFrag = data->caenFrags[i];
    unsigned int boardId = static_cast <unsigned int> (caenFrag.header.boardId);
    int deviceID = boardId;

    // each CAEN Trigger Time Tag count is 8 ns
    double timeStamp = caenFrag.header.triggerTimeTag * 0.008;  // convert to microseconds

    std::pair<double, double> fitParams(0, 0);

    fitParams = fitParamsMaps[8][deviceID].back();
    double corrTimeStamp = this->clockDriftCorr(fitParams, timeStamp);

    /*REL
    LOG_VERBATIM("FragmentToDigit") << "\n  deviceID:      " << deviceID
                                    << "\n  timeStamp:     " << timeStamp
                                    << "\n  corrTimeStamp: " << corrTimeStamp;
    */
    std::pair<int, int> dataBlockIndex(boardId, i);

    timeStampToDataBlockIndices.push_back(std::make_pair(corrTimeStamp, dataBlockIndex));

  }

  // correct TDCEvent timestamp and add corrected timestamp and index to vector
  for (size_t i = 0; i < numberTdcFrags; ++i) {

    TDCFragment const& tdcFrag = data->tdcFrags[i];

    std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcEvents = tdcFrag.tdcEvents;

    //LOG_DEBUG("FragmentToDigit")
    //    << "tdcEvents.size(): " << tdcEvents.size();
    //tdcFrag.print();

    int deviceID = 32;
    std::pair<double, double> fitParams(0, 0);
    fitParams = fitParamsMaps[8][deviceID].back();

    for (size_t j = 0; j < tdcEvents.size(); ++j) {

      if (tdcFrag.controllerHeader.nTDCs != tdcEvents[j].size()) {
        mf::LogError("FragmentToDigit") << "*** Fatal nTDCs mismatch: " << tdcEvents[j].size()
            << " != " << tdcFrag.controllerHeader.nTDCs<< " "<< j;
        continue;
      }

      //LOG_DEBUG("FragmentToDigit") << "TDC event: " << j;

      // Loop over TDCs to count the number of TDC timestamps that fall within range.
      // If that number is greater than 0, we consider this TDC event a match. We are
      // doing this to work around the mismatches of the TDC time bits.

      std::map<double, unsigned int> timeStampCounts;

      for (size_t tdc_index = 0; tdc_index < TDCFragment::MAX_TDCS; ++tdc_index) {
        TDCFragment::TdcEventData tdcEventData = tdcEvents[j].at(tdc_index);

        // each TDC Time Stamp count is 1/106.208 microseconds
        double timeStamp = tdcEventData.tdcEventHeader.tdcTimeStamp / 106.208;  // microseconds
        double corrTimeStamp = this->clockDriftCorr(fitParams, timeStamp);

        timeStampCounts[corrTimeStamp] += 1;

      } // end loop over TDCs

      unsigned int counts = 0;
      double corrTimeStamp = 0;

      for (auto const& k : timeStampCounts) {
        LOG_DEBUG("FragmentToDigit") << "timestamp: " << k.first
                                        << "\ncounts: " << k.second;
        if (k.second > counts) {
          corrTimeStamp = k.first;
          counts = k.second;
        }
      }

      std::pair<int, int> dataBlockIndex(32, j);

      timeStampToDataBlockIndices.push_back(std::make_pair(corrTimeStamp, dataBlockIndex));

    } // end loop over TDCEvents
  } // end loop over TDCFragments

  //LOG_VERBATIM("FragmentToDigit") << "Unsorted timestamps:";
  //for (auto const& i : timeStampToDataBlockIndices) {
  //  LOG_VERBATIM("FragmentToDigit") << "  timestamp: " << i.first;
  //}

  //std::sort(timeStampToDataBlockIndices.begin(), timeStampToDataBlockIndices.end());

  //LOG_VERBATIM("FragmentToDigit") << "Sorted timestamps:";

  //for (auto const& i : timeStampToDataBlockIndices) {
  //  LOG_VERBATIM("FragmentToDigit") << "  timestamp: " << i.first;
  //  //LOG_VERBATIM("FragmentToDigit") << "    index: (" << i.second.first << ", " << i.second.second << ")";
  //}

  std::vector<int> tstdbi_;  // master bookkeeper of already matched data blocks

  //REL  LOG_VERBATIM("FragmentToDigit") << "Sorted timestamps:";
  for (size_t i = 0; i < timeStampToDataBlockIndices.size(); ++i) {

    if (std::find(tstdbi_.begin(), tstdbi_.end(), i) != tstdbi_.end()) continue;

    double baseTimeStamp = timeStampToDataBlockIndices[i].first;
    //REL std::pair<int, int> baseIndexPair = timeStampToDataBlockIndices[i].second;

    //REL  LOG_VERBATIM("FragmentToDigit") << "  baseTimeStamp:    " << baseTimeStamp;
    //REL LOG_VERBATIM("FragmentToDigit") << "    baseIndexPair: (" << baseIndexPair.first << ", " << baseIndexPair.second << ")";

    double timeThresholdLow = baseTimeStamp - 0.032;
    double timeThresholdHigh = baseTimeStamp + 0.032;

    std::vector<int> selectedIndices;  // local bookkeeper  of selected data blocks

    for (size_t j = 0; j < timeStampToDataBlockIndices.size(); ++j) {

      if (std::find(tstdbi_.begin(), tstdbi_.end(), j) != tstdbi_.end()) continue;

      double timeStamp = timeStampToDataBlockIndices[j].first;
      //REL std::pair<int, int> indexPair = timeStampToDataBlockIndices[j].second;

      if (timeThresholdLow <= timeStamp and timeStamp <= timeThresholdHigh) {
	//REL LOG_VERBATIM("FragmentToDigit") << "    timeStamp:   " << timeStamp;
        //REL LOG_VERBATIM("FragmentToDigit") << "      indexPair: (" << indexPair.first << ", " << indexPair.second << ")";

        selectedIndices.push_back(j);  // add to local bookkeeper of selected data blocks
        tstdbi_.push_back(j);          // add to master bookkeeper of already matched data blocks
      } // if timestamp falls within range

    } // for each data block

    // now loop over the selected data blocks for this trigger ID
    // and see if we can add more data blocks that occur during
    // the readout windows of the CAEN V1740s and V1751s
    for (size_t j = 0; j < selectedIndices.size(); ++j) {
      int index = selectedIndices[j];
      std::pair<int, int> indexPair = timeStampToDataBlockIndices[index].second;
      int deviceID = indexPair.first;

      //REL LOG_VERBATIM("FragmentToDigit") << "  index: " << index;

      if (deviceID == 0) { // look for CAEN V1740 data blocks

        double timeStamp_ = timeStampToDataBlockIndices[index].first;
        double recordLength = numberV1740Samples * v1740SampleTime;  // microseconds

        double timeThresholdLow = timeStamp_ - recordLength;                // microseconds
        double timeThresholdHigh = timeStamp_ + recordLength;        // microseconds

        for (size_t k = 0; k < timeStampToDataBlockIndices.size(); ++k) {

          if (std::find(tstdbi_.begin(), tstdbi_.end(), k) != tstdbi_.end()) continue;

          double timeStamp = timeStampToDataBlockIndices[k].first;

          if (timeThresholdLow <= timeStamp and timeStamp <= timeThresholdHigh) {
            selectedIndices.push_back(k);  // add to local bookkeeper of selected data blocks
            tstdbi_.push_back(k);          // add to master bookkeeper of already matched data blocks
          } // if timestamp falls within range

        } // for each data block

      } // if CAEN board id is 0

      if (deviceID == 8) { // look for CAEN V1751 data blocks

        double timeStamp_ = timeStampToDataBlockIndices[index].first;
        double recordLength = numberV1751Samples * v1751SampleTime;  // microseconds

        double timeThresholdLow = timeStamp_;                // microseconds
        double timeThresholdHigh = timeStamp_ + recordLength;        // microseconds

        for (size_t k = 0; k < timeStampToDataBlockIndices.size(); ++k) {

          if (std::find(tstdbi_.begin(), tstdbi_.end(), k) != tstdbi_.end()) continue;

          double timeStamp = timeStampToDataBlockIndices[k].first;

          if (timeThresholdLow <= timeStamp and timeStamp <= timeThresholdHigh) {
            selectedIndices.push_back(k);  // add to local bookkeeper of selected data blocks
            tstdbi_.push_back(k);          // add to master bookkeeper of already matched data blocks
          } // if timestamp falls within range

        } // for each data block

      } // if CAEN board id is 8

    } // end loop over selected data blocks

    int v1740DataBlockCount = 0;
    int v1751DataBlockCount = 0;
    int WChamDataBlockCount = 0;

    // now that we have our final list of selected data blocks for this
    // trigger ID, let's add them to the trigger ID to data block maps
    for (size_t j = 0; j < selectedIndices.size(); ++j) {
      int index = selectedIndices[j];
      std::pair<int, int> indexPair = timeStampToDataBlockIndices[index].second;
      int deviceID = indexPair.first;
      int idx = indexPair.second;

      if (deviceID < 10 || deviceID == 24) {
        CAENFragment const& caenFrag = data->caenFrags[idx];
        unsigned int boardId = static_cast <unsigned int> (caenFrag.header.boardId);
        fTriggerToCAENDataBlocks[triggerID].push_back(caenFrag);
        numberMatchedCaenDataBlocks[boardId] += 1;
        if (deviceID == 0) ++v1740DataBlockCount; //Also need to check for the other v1740s somehow
        if (deviceID == 8) ++v1751DataBlockCount; //Also need to check for the other v17451 somehow
      }

      else if (deviceID == 32) {
        if (numberTdcFrags > 0) {
          TDCFragment const& tdcFrag = data->tdcFrags[0]; //The first and only spill's worth of TDC data.
          std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcEvents = tdcFrag.tdcEvents;
          fTriggerToTDCDataBlocks[triggerID].push_back(tdcEvents[idx]);
          numberMatchedMwpcDataBlocks += 1;
          ++WChamDataBlockCount;
        }//if numberTdcFrags > 0
      }//if deviceID == 32
    } // end loop over selected data blocks

    triggerID += 1;

    //REL LOG_VERBATIM("FragmentToDigit") << "  Trig: " << triggerID << " T/P/W:  " << v1740DataBlockCount <<"/"<<v1751DataBlockCount<<"/"<<WChamDataBlockCount;

    /*
    if (v1740DataBlockCount == 0)
      FragCountsSameTrigger_1751vsTDC_NoTPC    ->Fill(WChamDataBlockCount, v1751DataBlockCount);
    else if (v1740DataBlockCount == 1)
      FragCountsSameTrigger_1751vsTDC_WithTPC  ->Fill(WChamDataBlockCount, v1751DataBlockCount);
    else if (v1740DataBlockCount > 1)
      FragCountsSameTrigger_1751vsTDC_ExtraTPC ->Fill(WChamDataBlockCount, v1751DataBlockCount);
    */    
  } // for each data block

  //REL  LOG_VERBATIM("FragmentToDigit") << "  timeStampToDataBlockIndices.size(): " << timeStampToDataBlockIndices.size();
  //REL LOG_VERBATIM("FragmentToDigit") << "  tstdbi_.size(): " << tstdbi_.size();

  //////////////////////////////////////////////////////////////////////
  //@\ END: add to trigger ID to data block map
  //////////////////////////////////////////////////////////////////////

  // print matching summary

  //REL  LOG_VERBATIM("FragmentToDigit") << "\n";

  /*REL
  for (size_t i = 0; i < 32; ++i) {
    LOG_VERBATIM("FragmentToDigit") << "    boardId " << i << " matches: "
                                    << numberMatchedCaenDataBlocks[i] << " / "
                                    << numberCaenDataBlocks[i];
  }
  LOG_VERBATIM("FragmentToDigit") << "\n    MWPC matches: "
                                  << numberMatchedMwpcDataBlocks << " / "
                                  << numberMwpcDataBlocks;
  LOG_VERBATIM("FragmentToDigit") << "\n";
  */  

  return;
}

//-----------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::coarseMatch(int   const& deviceAID,
                                  int   const& deviceBID,
                                  double       range[2],      
                                  std::map< int, std::map<unsigned int, double> > timeStamps,
                                  match_maps & matchMaps) 
{

  std::map< unsigned int, std::vector<unsigned int> > matchAB;

  size_t numberADataBlocks = timeStamps[deviceAID].size();
  size_t numberBDataBlocks = timeStamps[deviceBID].size();

  for (size_t a = 0; a < numberADataBlocks; ++a) {
    double timeStampA = timeStamps[deviceAID][a];
    for (size_t b = 0; b < numberBDataBlocks; ++b) {
      double timeStampB = timeStamps[deviceBID][b];
      double difference = timeStampA - timeStampB;
      if (range[0] <= difference and difference <= range[1]) {
        matchAB[a].push_back(b);
      }
    }
  }

  matchMaps[deviceAID][deviceBID].push_back(matchAB);

  return;
}

//-----------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::fineMatch(int      const& deviceAID,
                                int      const& deviceBID,
                                double             eps[2],      
                                fit_params_maps fitParametersMaps,
                                std::map< int, std::map<unsigned int, double> > timeStamps,
                                match_maps    & matchMaps) 
{

  std::pair<double, double> fitParameters = fitParametersMaps[deviceAID][deviceBID].back();

  std::map< unsigned int, std::vector<unsigned int> > matchAB;

  size_t numberADataBlocks = timeStamps[deviceAID].size();
  size_t numberBDataBlocks = timeStamps[deviceBID].size();

  for (size_t a = 0; a < numberADataBlocks; ++a) {
    double timeStampA = timeStamps[deviceAID][a];
    for (size_t b = 0; b < numberBDataBlocks; ++b) {
      double timeStampB = timeStamps[deviceBID][b];
      double difference = timeStampA - timeStampB;
      double y = this->line(fitParameters, timeStampA);
      double yLow = y - eps[0];
      double yHigh = y + eps[1];
      if (yLow <= difference and difference <= yHigh) {
        matchAB[a].push_back(b);
      }
    }
  }

  matchMaps[deviceAID][deviceBID].push_back(matchAB);

  return;
}

//-----------------------------------------------------------------------------------
double DAQToOffline::SlicerInput::line(std::pair<double, double> const& parameters, 
                             double                    const& x) 
{
  double intercept = parameters.first;
  double slope = parameters.second;
  return intercept + slope * x;
}

//-----------------------------------------------------------------------------------
double DAQToOffline::SlicerInput::clockDriftCorr(std::pair<double, double> const& parameters, 
                                       double                    const& x) 
{
  double intercept = parameters.first;
  double slope = parameters.second;
  return (intercept + x) / (1 - slope);
}

//-----------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::fitClockDrift(int         const& deviceAID,
                                    int         const& deviceBID,
                                    std::map< int, std::map<unsigned int, double> > timeStamps,
                                    match_maps       & matchMaps,
                                    fit_params_maps  & fitParametersMaps,
                                    std::string const& graphNamePrefix) 
{

  std::vector<double> x;
  std::vector<double> y;

  std::vector<double> x_;
  std::vector<double> y_;

  std::map< unsigned int, std::vector<unsigned int> > matchAB = matchMaps[deviceAID][deviceBID].back();

  LOG_DEBUG("FragmentToDigit") << "Matches between " << deviceAID << " and " << deviceBID;
  LOG_DEBUG("FragmentToDigit") << "  Matching index pair format: " << "(" << deviceAID << ", " << deviceBID << ")";

  for (auto const & itA : matchAB) {
    LOG_DEBUG("FragmentToDigit") << "  Index: " << itA.first << "; number of matches: " << itA.second.size();
    for (auto const & itB : itA.second) {
      double timeStampA = timeStamps[deviceAID][itA.first];
      double timeStampB = timeStamps[deviceBID][itB];
      double difference = timeStampA - timeStampB;
      double reference = timeStampA;
      LOG_DEBUG("FragmentToDigit") << "    (" << itA.first << ", " << itB << ")"
                                   << "; difference: " << difference << " usec"
                                   << "; reference: " << reference << " usec";
      if (itA.second.size() == 1) {
        x.push_back(reference);
        y.push_back(difference);
      }
      else {
        x_.push_back(reference);
        y_.push_back(difference);
      }
    }
  }

  std::pair<double, double> intSlp(0., 0.);

  // try {
  //   this->LinFitUnweighted(x, y, intSlp.second, intSlp.first);
  // }
  // catch (cet::exception &e) {
  //   LOG_WARNING("FragmentToDigit") << "caught exception:\n" << e
  //                                  << "\n returning intercept = 0 and slope = 0";
  // }

  // return intSlp;

  std::string graphName = graphNamePrefix + "_drift_" + std::to_string(deviceAID) + "_" + std::to_string(deviceBID);
  std::string graphTitles = ("; Time since beginning of spill (using deviceID " +
                             std::to_string(deviceAID) +
  			                 " clock) [#mus]; #Delta t between device ID " +
                             std::to_string(deviceAID) +
  			                 " and device ID " +
                             std::to_string(deviceBID) +
                             " [#mus]");

  TGraph * graph = tfs->make<TGraph>(x.size(), &x[0], &y[0]);

  try {
    TF1 * f = tfs->make<TF1>("f", "pol1", 0, 30e6);
    graph->Fit("f", "Q");
    LOG_DEBUG("FragmentToDigit") << "Fit parameters: intercept, "
                                 << f->GetParameter(0) << " usec; slope, "
                                 << f->GetParameter(1) << " usec/usec";
    intSlp.first = f->GetParameter(0);
    intSlp.second = f->GetParameter(1);
  }
  catch (cet::exception &e) {
    LOG_WARNING("FragmentToDigit") << "caught exception:\n" << e
                                   << "\nTLinearFitter failed"
                                   << "\n returning intercept = 0 and slope = 0";
  }

  graph->SetMarkerStyle(20);
  graph->SetTitle(graphTitles.c_str());
  graph->Write(graphName.c_str());
  // return std::make_pair<double, double>(f->GetParameter(0), f->GetParameter(1));

  fitParametersMaps[deviceAID][deviceBID].push_back(intSlp);

  return;
}

//------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::matchFitIter(int        const& deviceAID,
                                   int        const& deviceBID,
                                   std::map< int, std::map<unsigned int, double> > timeStamps,
                                   match_maps      & matchMaps,
                                   fit_params_maps & fitParametersMaps,
                                   double            coarseRange[2],
                                   double            fineEps[2])
{

  LOG_DEBUG("FragmentToDigit") << "matchFitIter -- "
                               << "deviceA: " << deviceAID << "; deviceB: " << deviceBID;

  this->coarseMatch(deviceAID, deviceBID, coarseRange, timeStamps, matchMaps);
  this->fitClockDrift(deviceAID, deviceBID, timeStamps, matchMaps, fitParametersMaps, "coarse_match");

  std::pair<double, double> fitParameters = fitParametersMaps[deviceAID][deviceBID].back();
  LOG_DEBUG("FragmentToDigit") << "  intercept: " << fitParameters.first << "; slope: " << fitParameters.second;

  for (size_t i = 0; i < fMaxNumberFitIterations; ++i) {

    this->fineMatch(deviceAID, deviceBID, fineEps, fitParametersMaps, timeStamps, matchMaps);
    this->fitClockDrift(deviceAID, deviceBID, timeStamps, matchMaps, fitParametersMaps, "fine_match_" + std::to_string(i));

    fitParameters = fitParametersMaps[deviceAID][deviceBID].back();
    LOG_DEBUG("FragmentToDigit") << "  intercept: " << fitParameters.first << "; slope: " << fitParameters.second;

    std::pair<double, double> prevFitParameters = fitParametersMaps[deviceAID][deviceBID].end()[-2];

    if (fitParameters == prevFitParameters) {
      LOG_DEBUG("FragmentToDigit") << "  Fit parameters did not change!";
      break;
    }

  }

  return;
}

//------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::printMatchMap(int   const& deviceAID,
                                    int   const& deviceBID,
                                    match_maps & matchMaps) 
{

  std::map< unsigned int, std::vector<unsigned int> > matchAB = matchMaps[deviceAID][deviceBID].back();

  LOG_DEBUG("FragmentToDigit") << "Matches between " << deviceAID << " and " << deviceBID
                               << "\n Matching index pair format: "
                               << "(deviceA " << deviceAID
                               << ", deviceB " << deviceBID << ")";

  for (auto const & itA : matchAB) {
    for (auto const & itB : itA.second) {
      LOG_DEBUG("FragmentToDigit") << "  (" << itA.first << ", " << itB << ")";
    }
  }

  return;
}

//------------------------------------------------------------------------------
void DAQToOffline::SlicerInput::commenceRun( art::RunPrincipal*& outR )
{
  fFrag2DigAlg.InitializeRun(fRunNumber);

  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  //  std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
  sumdata::RunData runcol = sumdata::RunData(geo->DetectorName());
  art::put_product_in_principal( std::make_unique<sumdata::RunData>(runcol),
				 *outR,
				 fSourceName);
				 
  return;
}


//=======================================================================================
DEFINE_ART_INPUT_SOURCE(art::Source<DAQToOffline::SlicerInput>)
//=======================================================================================
