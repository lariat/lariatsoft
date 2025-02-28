////////////////////////////////////////////////////////////////////////
// Class:       WCEfficiencies
// Module Type: analyzer
// File:        WCEfficiencies_module.cc
//
// Generated at Tue Aug  9 14:35:14 2016 by Lucas Mendes Santos using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

// LArSoft Libraries
#include "lardataobj/RawData/AuxDetDigit.h"

//ROOT Libraries
#include "TTree.h"

class WCEfficiencies;

class WCEfficiencies : public art::EDAnalyzer {
public:
  explicit WCEfficiencies(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCEfficiencies(WCEfficiencies const &) = delete;
  WCEfficiencies(WCEfficiencies &&) = delete;
  WCEfficiencies & operator = (WCEfficiencies const &) = delete;
  WCEfficiencies & operator = (WCEfficiencies &&) = delete;
  void beginJob() override;
  void  reconfigure(fhicl::ParameterSet const & p) ;
  
  // Required functions.
  void analyze(art::Event const & e) override;


private:

  // Declare member data here.
  std::string fAuxDetDigit;
  
  TTree *fTree;

};


WCEfficiencies::WCEfficiencies(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ this->reconfigure(p);}

void WCEfficiencies::reconfigure(fhicl::ParameterSet const & p)
{
   fAuxDetDigit = p.get<std::string>("DAQ", "daq");
}

void WCEfficiencies::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("wceff","WC Efficiencies");
}
void WCEfficiencies::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  
  
}

DEFINE_ART_MODULE(WCEfficiencies)
