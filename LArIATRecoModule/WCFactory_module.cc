////////////////////////////////////////////////////////////////////////
// Class:       WCFactory
// Module Type: producer
// File:        WCFactory_module.cc
//
// Generated at Wed Oct 21 15:46:13 2015 by Greg Pulliam using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"

//ROOT Things
#include <TH1F.h>

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
//#include "LArIATRecoAlg/WCTrackBuilderAlgBase.h"
#include "LArIATRecoAlg/WCHitFinderAlg.h"
#include "LArIATDataProducts/WCTrack.h"
#include "Utilities/DatabaseUtilityT1034.h"
#include "LArIATRecoAlg/WCTrackBuilderFactoryTest.h"



#include <utility>
#include <string>
#include <memory>

class WCFactory;

class WCFactory : public art::EDProducer {
public:
  explicit WCFactory(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCFactory(WCFactory const &) = delete;
  WCFactory(WCFactory &&) = delete;
  WCFactory & operator = (WCFactory const &) = delete;
  WCFactory & operator = (WCFactory &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  std::string WCTrackBuilderAlgName;
 // WCTrackAlgBase *pAlg;
};


WCFactory::WCFactory(fhicl::ParameterSet const & p)  
// :
// Initialize member data here.
{
   WCTrackBuilderAlgName = p.get<std::string>("FactoryAlgLabel"); //The name of the algorithm to 

  // Call appropriate produces<>() functions here.
}

void WCFactory::produce(art::Event & e)
{

  // Implementation of required member function here.
}

void WCFactory::beginJob()
{
      WCTrackAlgBase *pAlg = NULL;
      pAlg = WCTrackBuilderFactory::Get()->CreateAlg(WCTrackBuilderAlgName);
      
      std::cout << "////////////////////////////////////////////////" << std::endl;
      if (pAlg){
      std::cout << "It works?" << std::endl;
      pAlg->Hello();
      }
      else{
      std::cout<<"We dun goofd"<<std::endl;
      }
      std::cout << "////////////////////////////////////////////////" << std::endl;
		
      
  // Implementation of optional member function here.
}

void WCFactory::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WCFactory::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WCFactory::endJob()
{
  // Implementation of optional member function here.
}

void WCFactory::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WCFactory::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WCFactory::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

void WCFactory::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCFactory::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCFactory::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCFactory::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(WCFactory)
