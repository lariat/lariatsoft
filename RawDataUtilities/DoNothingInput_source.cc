//-------------------------------------------------------------------------
//
// Name:   DoNothingInput_source.cc
// Date:   10 October 2016
// Author: Everybody is an author!
//
//-------------------------------------------------------------------------
//
// The purpose of this source is to test the order in which the specified
// input files are opened by the readFile function.
//
//-------------------------------------------------------------------------

#ifndef DoNothingInput_source
#define DoNothingInput_source

// art includes
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <string>
#include <vector>

//-------------------------------------------------------------------------
namespace rdu
{
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class DoNothingInput
  {

   public:

    // constructor and destructor
    explicit DoNothingInput(fhicl::ParameterSet        const& pset,
                            art::ProductRegistryHelper      & prhelper,
                            art::SourceHelper               & shelper);
    virtual ~DoNothingInput();

    bool readFile(std::string const& filename, art::FileBlock * & fileblock);

    bool readNext(art::RunPrincipal    * const& inRun,
                  art::SubRunPrincipal * const& inSubRun,
                  art::RunPrincipal    *      & outRun,
                  art::SubRunPrincipal *      & outSubRun,
                  art::EventPrincipal  *      & outEvent);

    void closeCurrentFile();

    void reconfigure(fhicl::ParameterSet const& pset);

   private:

    std::vector< std::string > fFileNames;

  }; // class DoNothingInput

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  DoNothingInput::DoNothingInput(fhicl::ParameterSet        const& pset,
                                 art::ProductRegistryHelper      & prhelper,
                                 art::SourceHelper               & shelper)
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);

    mf::LogVerbatim("DoNothingInput")
      << "You have specified the following " << fFileNames.size()
      << " file(s):";

    for (auto const& fileName : fFileNames)
    {
      mf::LogVerbatim("DoNothingInput") << "  " << fileName << std::endl;
    }
  }

  //-----------------------------------------------------------------------
  // destructor
  DoNothingInput::~DoNothingInput()
  {}

  //-----------------------------------------------------------------------
  void DoNothingInput::reconfigure(fhicl::ParameterSet const& pset)
  {
    fFileNames = pset.get< std::vector< std::string > >("fileNames", {});
    return;
  }

  //-----------------------------------------------------------------------
  bool DoNothingInput::readFile(std::string const& filename, art::FileBlock * & fileblock)
  {
    // new fileblock
    fileblock = new art::FileBlock(art::FileFormatVersion(), filename);
    if (fileblock == nullptr)
    {
      throw art::Exception(art::errors::FileOpenError)
        << "Unable to open file " << filename << "\n";
    }

    mf::LogInfo("DoNothingInput") << "Opened input file " << filename;

    return true;
  }

  //-----------------------------------------------------------------------
  bool DoNothingInput::readNext(art::RunPrincipal    * const& inRun,
                                art::SubRunPrincipal * const& inSubRun,
                                art::RunPrincipal    *      & outRun,
                                art::SubRunPrincipal *      & outSubRun,
                                art::EventPrincipal  *      & outEvent)
  {
    return false;
  }

  //-----------------------------------------------------------------------
  void DoNothingInput::closeCurrentFile()
  {
  }

  DEFINE_ART_INPUT_SOURCE(art::Source<rdu::DoNothingInput>)

} // namespace rdu

#endif // DoNothingInput_source
