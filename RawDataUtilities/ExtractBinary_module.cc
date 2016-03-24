//////////////////////////////////////////////////////////////
// Name:      ExtractBinary_module.cc
// Date:      21 March 2016
// Author:    JCF
//////////////////////////////////////////////////////////////

// JCF, 10/24/14

// "ExtractBinary" is a simple utility Art module which will take the
// entire payload of the artdaq::Fragment object and save it to disk
// (using the filename "<output_filename_base>_r<run>_sr<subrun>.dat")
// as a binary file. The main purpose for this module is to be able to
// cross check that the payload is the same as a pre-existing binary
// which was processed by the SpillFragmentReader fragment generator
// (modulo the 8-byte-wide padding at the end of the payload expected
// in artdaq)

#ifndef ExtractBinary_Module
#define ExtractBinary_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// artdaq includes
#include "artdaq-core/Data/Fragments.hh"

// C++ includes
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

namespace ExtractBinary {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class ExtractBinary : public art::EDAnalyzer 
  {

   public:

    // Standard constructor and destructor for an ART module.
    explicit ExtractBinary(fhicl::ParameterSet const& pset);
    virtual ~ExtractBinary();

    // This method is called once, at the start of the job.
    void beginJob();

    // This method is called once, at the start of each run.
    void beginRun(const art::Run& run);

    // This method is called once, at the start of each sub-run
    void beginSubRun(const art::SubRun& subrun);

    // This method reads in any parameters from the .fcl files.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze(const art::Event& event); 

   private:

    std::string raw_fragment_label_;    ///< label for module producing artdaq fragments
    std::string raw_fragment_instance_; ///< instance label for artdaq fragments

    std::string output_filename_base_;
    std::stringstream output_filename_stream_;
    std::string output_filename_;

    std::ofstream binary_file_;

  }; // class ExtractBinary


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // constructor
  ExtractBinary::ExtractBinary(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , output_filename_base_(pset.get<std::string>("output_filename_base", "extract_binary"))
  {
    // read in the parameters from the .fcl file
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // destructor
  ExtractBinary::~ExtractBinary() 
  {}

  //-----------------------------------------------------------------------
  void ExtractBinary::beginJob()
  {}

  //-----------------------------------------------------------------------
  void ExtractBinary::beginRun(const art::Run& run)
  {}

  //-----------------------------------------------------------------------
  void ExtractBinary::beginSubRun(const art::SubRun& subrun)
  {

    output_filename_stream_ << output_filename_base_ << "_r" << std::setfill('0') << std::setw(6) << subrun.run() << "_sr" << std::setfill('0') << std::setw(4) << subrun.subRun() << ".dat";
    output_filename_ = output_filename_stream_.str();

    if (binary_file_.is_open()) {
      binary_file_.close();
    }

    binary_file_.open(output_filename_, std::ios::out);

    if (! binary_file_.is_open()) {
      throw cet::exception("ExtractBinary") << "Error in ExtractBinary: unable to open " <<
    output_filename_ ;
    }

  }

  //-----------------------------------------------------------------------
  void ExtractBinary::reconfigure(fhicl::ParameterSet const& pset)
  {
    raw_fragment_label_    = pset.get< std::string >("raw_fragment_label",    "daq"  );
    raw_fragment_instance_ = pset.get< std::string >("raw_fragment_instance", "SPILL");

    return;
  }

  //-----------------------------------------------------------------------
  void ExtractBinary::analyze(const art::Event& event) 
  {

    art::Handle<artdaq::Fragments> raw;

    event.getByLabel(raw_fragment_label_, raw_fragment_instance_, raw);

    if ( !raw.isValid()) {

      throw cet::exception("Error in ExtractBinary_module: expected at least one valid spill fragment");

    } else {

      std::cout << "----------------------------------------------------------------------" << std::endl;
      for (size_t i_f = 0; i_f < raw->size(); ++i_f) {

        const auto& spillfrag((*raw)[i_f]);

        std::cout << "From " << reinterpret_cast<const void*>( &*spillfrag.dataBegin() ) <<
      ", will write " << spillfrag.dataSize()*sizeof(artdaq::RawDataType) <<
      " bytes into " << output_filename_ << std::endl;

        binary_file_.write( reinterpret_cast<const char*>( &*spillfrag.dataBegin() ), spillfrag.dataSize()*sizeof(artdaq::RawDataType) );
      }
      std::cout << "----------------------------------------------------------------------" << std::endl;
    }


    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file.
  DEFINE_ART_MODULE(ExtractBinary)

} // namespace ExtractBinary

#endif // ExtractBinary_Module
