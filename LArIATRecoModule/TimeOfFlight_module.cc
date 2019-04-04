/*

LArIATSoft 

Time Of Flight producer module

Module: TimeOfFlight_module.cc 
FHiCL: TimeOfFlight.fcl
Header: LArIATDataProducts/TOF.h
Class: LArIATDataProducts/TOF.cxx
Dictionary: LArIATDataProducts/classes.h and classes_def.xml

Authores:
Elena Gramellini - elena.gramellini@yale.edu
Irene Nutini     - irene.nutini@stud.unifi.it
Daniel Smith     - dsmith@fnal.gov

*/

// Art Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"

//C++ Includes
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>

//ROOT Includes
#include <TH1F.h>

// LArSoft Includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/TriggerData.h"

//LAriatSoft Includes
#include "LArIATDataProducts/TOF.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/TOFBuilderAlg.h"

namespace lrm {
  class TimeOfFlight;
}

class lrm::TimeOfFlight : public art::EDProducer {
public:
  explicit TimeOfFlight(fhicl::ParameterSet const & p);

  TimeOfFlight(TimeOfFlight const &) = delete;
  TimeOfFlight(TimeOfFlight &&) = delete;
  TimeOfFlight & operator = (TimeOfFlight const &) = delete;
  TimeOfFlight & operator = (TimeOfFlight &&) = delete;

  void produce(art::Event & e) override;


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

  ///< Label for the module producing the triggers
  std::string fTriggerUtility; 
  TOFBuilderAlg fTOFAlg;
};

//----------------------------------------------------------------
lrm::TimeOfFlight::TimeOfFlight(fhicl::ParameterSet const & p)
  : fTOFAlg(p)
{

  // Configures the ROOT histograms and the 
  this->reconfigure(p);  

  // Produces the LArSoft object to be ultimately outputted
  produces<std::vector<ldp::TOF> >();
  produces<art::Assns<raw::Trigger, ldp::TOF> >();

}

void lrm::TimeOfFlight::produce(art::Event & e)
{
 
  // Setting up to begin looping over the Triggers in the inputted ROOT file

  rdu::TriggerDigitUtility tdu(e, fTriggerUtility);    
  art::PtrVector<raw::Trigger> const& EventTriggersPtr = tdu.EventTriggersPtr();

  // Associating the trigger and the TOF to be used in outputting the information to a ROOT file
  std::unique_ptr<art::Assns<raw::Trigger, ldp::TOF> > TriggerTOFAssn(new art::Assns<raw::Trigger, ldp::TOF>);
  
  // Creating the object for the TOF, the object that will be outputted
  std::unique_ptr<std::vector<ldp::TOF> > TOFCol(new std::vector<ldp::TOF> );

  // Loop over the triggers
  for(size_t trig = 0; trig < tdu.NTriggers(); ++trig) {

    // Getting a the current trigger
    art::Ptr<raw::Trigger> theTrigger = (EventTriggersPtr[trig]);

    // Retrieve the digits for the upstream and downstream paddles
    std::vector<const raw::AuxDetDigit*> ust_wv = tdu.TriggerUpStreamTOFDigits(trig);
    std::vector<const raw::AuxDetDigit*> dst_wv = tdu.TriggerDownStreamTOFDigits(trig);

    
    // Tests if the event has 2 PMTs, the amount needed for analysis
    if(ust_wv.size() == 2) { 
      
      // Time For some LArSoft magic. 
      // Creates the object, pushes it, and Associates it to the trigger.
      // Variables of our object
      
      std::pair<std::vector<float>, std::vector<long> > pair;
      pair = fTOFAlg.get_TOF_and_TimeStamp(ust_wv,dst_wv);

      ldp::TOF TOFObject(pair.first, pair.second);
      (*TOFCol).push_back( TOFObject );
      util::CreateAssn(*this, e, *TOFCol, theTrigger, *TriggerTOFAssn);
    
    }
  }

  // Move LArSoft magic to save the information into the final ROOT file 
  e.put(std::move(TriggerTOFAssn));
  e.put(std::move(TOFCol));

  return;
}





void lrm::TimeOfFlight::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
}

void lrm::TimeOfFlight::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::endJob()
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::reconfigure(fhicl::ParameterSet const & p)
{
   
   // Implementing the ability to pass the name of the TriggerUtililty 
   fTriggerUtility   = p.get< std::string >("TriggerUtility", "FragmentToDigit");

}

void lrm::TimeOfFlight::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(lrm::TimeOfFlight)
