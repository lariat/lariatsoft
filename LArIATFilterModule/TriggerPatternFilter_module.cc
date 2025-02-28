////////////////////////////////////////////////////////////////////////
// Class:       TriggerPatternFilter
// Module Type: filter
// File:        TriggerPatternFilter_module.cc
//
// Generated at Thu Aug  6 14:14:55 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArIAT includes
#include "LArIATRecoAlg/TriggerFilterAlg.h"

#include <memory>

class TriggerPatternFilter;

class TriggerPatternFilter : public art::EDFilter {
public:
  explicit TriggerPatternFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerPatternFilter(TriggerPatternFilter const &) = delete;
  TriggerPatternFilter(TriggerPatternFilter &&) = delete;
  TriggerPatternFilter & operator = (TriggerPatternFilter const &) = delete;
  TriggerPatternFilter & operator = (TriggerPatternFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  bool beginRun(art::Run & r) override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const  &fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  std::string fSlicerSourceLabel;
  std::string fFilterPattern;
  TriggerFilterAlg fTriggerFilterAlg;
};


TriggerPatternFilter::TriggerPatternFilter(fhicl::ParameterSet const & p)
: EDFilter(p),
  fTriggerFilterAlg(p.get< fhicl::ParameterSet > ("TriggerFilterAlg"))
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool TriggerPatternFilter::filter(art::Event & e)
{
  // Implementation of required member function here.
  //Get the triggers from the event record
  art::Handle< std::vector<raw::Trigger> > triggerHandle;
  e.getByLabel(fSlicerSourceLabel,triggerHandle);

  std::cout << "In body of filter module." << std::endl;

  //Kill an event from consideration in the event record if it does not pass the filter pattern
  std::string myFilterPattern = fFilterPattern;
  bool filterPassed = fTriggerFilterAlg.doesTriggerPassFilter( triggerHandle->at(0), fFilterPattern );
  std::cout << "Returning: " << filterPassed << std::endl;
  return filterPassed;

  

}

void TriggerPatternFilter::beginJob()
{
  // Implementation of optional member function here.

}

bool TriggerPatternFilter::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
  fTriggerFilterAlg.loadXMLDatabaseTable( r.run() );
  return true;
}

void TriggerPatternFilter::endJob()
{
  // Implementation of optional member function here.
}

void TriggerPatternFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fSlicerSourceLabel          = p.get<std::string>("SourceLabel", "daq");
  fFilterPattern              = p.get<std::string>("TriggerPattern");
}

void TriggerPatternFilter::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void TriggerPatternFilter::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void TriggerPatternFilter::respondToOpenInputFile(art::FileBlock const  &fb)
{
  // Implementation of optional member function here.
}

void TriggerPatternFilter::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerPatternFilter)
