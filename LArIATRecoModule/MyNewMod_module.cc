////////////////////////////////////////////////////////////////////////
// Class:       MyNewMod
// Module Type: producer
// File:        MyNewMod_module.cc
//
// Generated at Fri May 29 10:13:48 2015 by Jonathan Asaadi using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>
#include <TH1F.h>
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
// ### LArIAT Things ###
#include "RawDataUtilities/TriggerDigitUtility.h"

#include <memory>

namespace lrm {
  class MyNewMod;
}

class lrm::MyNewMod : public art::EDProducer {
public:
  explicit MyNewMod(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyNewMod(MyNewMod const &) = delete;
  MyNewMod(MyNewMod &&) = delete;
  MyNewMod & operator = (MyNewMod const &) = delete;
  MyNewMod & operator = (MyNewMod &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) ;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;
    

private:
  
  std::string fTriggerUtility; //<---Label for the module producing the triggers
  TH1F* fWC1Size;
  TH1F* fWC1XPlaneADC;
  
  // Declare member data here.

};


lrm::MyNewMod::MyNewMod(fhicl::ParameterSet const & p) //:fTrigFiltAlg(p.get< fhicl::ParameterSet > ("TriggerFilterAlg"))
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

void lrm::MyNewMod::produce(art::Event & e)
{
  // Implementation of required member function here.
  
  // ---- Example of how to make a collection of new object you are putting on the event ---
  //std::unique_ptr<std::vector<recob::WireChamberTrack> > WCtrackcol(new std::vector<recob::WireChamberTrack>);
  
  
  // ---- Example of how to make an association between your new object and the trigger object -----
  //std::unique_ptr<art::Assns<raw::Trigger, recob::WireChamberTrack> > TrigWCtrackAssn(new art::Assns<raw::Trigger, recob::WireChamberTrack>);
  
  // ###########################################
  // ### Grab the trigger data utility (tdu) ###
  // ###########################################
  
  // ### This is a crap way to pass a string.....FIX ME!!!! ####
  fTriggerUtility = "FragmentToDigit";
  rdu::TriggerDigitUtility tdu(e, fTriggerUtility);  //input event and name of whatever created our digits
    
    
  // ##############################
  // ### Loop over the triggers ###
  // ##############################
  for(size_t trig = 0; trig < tdu.NTriggers(); ++trig)
     {
     
       //### Global Trigger that you will use to make associations ###
       //       art::Ptr<raw::Trigger> trigger = tdu.EventTriggersPtr()[trig];

       // ### Getting a vector of Wire Chamber 1 Digits ###
       std::vector<const raw::AuxDetDigit*> WireChamber1Digi = tdu.TriggerMWPC1Digits(trig);
       // art::PtrVector<raw::AuxDetDigit> WireChamber1Digi = tdu.TriggerWMPC1DigitsPtr(trig); //if you want to use more correct pointer instead
       
       std::cout<<"Size of WC1Digi = "<<WireChamber1Digi.size()<<std::endl;
       fWC1Size->Fill(WireChamber1Digi.size());
       
       // ### Looping over Wire Chamber #1's Digits ###
       for(size_t wc1 = 0; wc1 < WireChamber1Digi.size(); ++wc1)
	 {
	   //std::cout<<"WireChamber1Digi Channel Number = "<<WireChamber1Digi.at(wc1)->Channel()<<std::endl;
	   ///std::cout<<"WireChamber1Digi TDC = "<<WireChamber1Digi.at(wc1)->NADC()<<std::endl;
	   auto wcDigit = WireChamber1Digi[wc1];
	   
	   // ### Looping over all the TDC hits (annoyingly called ADC's) ###
	   for (size_t i =0; i < wcDigit->NADC(); ++i)
	     {
	       
	       // ### Skipping any TDC hits (which returns the time tick) that is zero ###
	       if(wcDigit->ADC(i) == 0){continue;}
	       
	       
	       //### Clustering over X Plane ###
	       if(wcDigit->Channel() < 128) //<---(Channels 0 - 127 are the X Plane)
		 {
		   // Get wire and time pairs, cluster with dBScan (or whatever), save the good hits
		   fWC1XPlaneADC->Fill(wcDigit->ADC(i));
		 }//<---End X Plane
	       
	       // ### Clustering over the Y Plane ###
	       if(wcDigit->Channel() > 127) //<---(Channels 128 - 255 are the Y Plane)
		 {
		   
		   // Get wire and time pairs, cluster with dBScan (or whatever), save the good hits
		   
		 }//<---End Y Plane
	       
	       
	     }// <--- End i loop
	   
	 }//<--End wc1
       
     }//<---End trig loop
  
  
}

void lrm::MyNewMod::beginJob()
{
  // Implementation of optional member function here.
  
  art::ServiceHandle<art::TFileService> tfs;
  fWC1Size = tfs->make<TH1F>("WC1Size","WC1Size", 100, 0, 100);
  fWC1XPlaneADC  = tfs->make<TH1F>("WC1XPlaneADC", "WC1ADC", 1024, 0, 1024);
}

void lrm::MyNewMod::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::beginSubRun(art::SubRun & sr)
{
   

}

void lrm::MyNewMod::endJob()
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::reconfigure(fhicl::ParameterSet const & p)
{

   /// Here lies some pseudo-code for making objects and their associations 
   
   //produces< std::vector<recob::WireChamberTrack>>(); (DON'T USE ME EXPLICITLY....FOR ILLUSTRATION ONLY)
   //produces< art::Assns<raw::Trigger, recob::WireChamberTrack>>(); (DON'T USE ME EXPLICITLY....FOR ILLUSTRATION ONLY)
   
   // implementing the ability to pass the name of the TriggerUtililty 
   fTriggerUtility   = p.get< std::string >("TriggerUtility", "FragmentToDigit");

  // Implementation of optional member function here.
}

void lrm::MyNewMod::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::MyNewMod::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(lrm::MyNewMod)
