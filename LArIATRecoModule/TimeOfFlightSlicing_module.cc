/*

LArIATSoft 

Time Of Flight producer module

Module: TimeOfFlightSlicing_module.cc 
FHiCL: TimeOfFlightSlicing.fcl
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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

//C++ Includes
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>

//ROOT Includes
#include <TH1F.h>

// LArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "RawData/TriggerData.h"

//LAriatSoft Includes
#include "LArIATDataProducts/TOF.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/TOFBuilderAlg.h"

namespace lrm {
  class TimeOfFlightSlicing;
}

class lrm::TimeOfFlightSlicing : public art::EDProducer {
public:
  explicit TimeOfFlightSlicing(fhicl::ParameterSet const & p);

  TimeOfFlightSlicing(TimeOfFlightSlicing const &) = delete;
  TimeOfFlightSlicing(TimeOfFlightSlicing &&) = delete;
  TimeOfFlightSlicing & operator = (TimeOfFlightSlicing const &) = delete;
  TimeOfFlightSlicing & operator = (TimeOfFlightSlicing &&) = delete;

  void produce(art::Event & e) override;


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

  ///< Label for the module producing the triggers
  std::string fTriggerUtility; 
  TOFBuilderAlg fTOFAlg;
  TH1F*         fTOFHisto;
  bool          fMakeHistograms;
  std::string   fSlicerSourceLabel;
};

//----------------------------------------------------------------
lrm::TimeOfFlightSlicing::TimeOfFlightSlicing(fhicl::ParameterSet const & p)
  : fTOFAlg(p)
{

  // Configures the ROOT histograms and the 
  this->reconfigure(p);  

  // Produces the LArSoft object to be ultimately outputted
  produces<std::vector<ldp::TOF> >();
  produces<art::Assns<raw::Trigger, ldp::TOF> >();

}

void lrm::TimeOfFlightSlicing::produce(art::Event & e)
{

  //Clearing the persistent members so we don't have repetition
  fTOFAlg.clear_tof_and_timeStampDst();
 
  // Creating the object for the TOF, the object that will be outputted
  std::unique_ptr<std::vector<ldp::TOF> > TOFCol(new std::vector<ldp::TOF> );

  //Retrieving the digits for the upstream and downstream paddles from the sliced event
  art::Handle< std::vector<raw::AuxDetDigit> > AuxDetDigitHandle;
  e.getByLabel(fSlicerSourceLabel,"SPILL",AuxDetDigitHandle);

  //Determine which digits are the USTOF and DSTOF ones
  std::vector<const raw::AuxDetDigit*> USTOF;
  std::vector<const raw::AuxDetDigit*> DSTOF;
  for( size_t iDig = 0; iDig < AuxDetDigitHandle->size() ; ++iDig ){
    if( AuxDetDigitHandle->at(iDig).AuxDetName() == "TOFUS" ) USTOF.push_back(&(AuxDetDigitHandle->at(iDig)));
    if( AuxDetDigitHandle->at(iDig).AuxDetName() == "TOFDS" ) DSTOF.push_back(&(AuxDetDigitHandle->at(iDig)));
  }
  
  // Tests if the event has 2 PMTs, the amount needed for analysis
  if(USTOF.size() == 2) { 
    
    // Time For some LArSoft magic. 
    // Creates the object, pushes it, and Associates it to the trigger.
    // Variables of our object
    
    std::pair<std::vector<short>, std::vector<long> > pair;
    pair = fTOFAlg.get_TOF_and_TimeStamp(USTOF,DSTOF);
    
    ldp::TOF TOFObject(pair.first, pair.second);
    (*TOFCol).push_back( TOFObject );
    
    //Fill histos if necessary
    if( fMakeHistograms ){
      std::cout << "NTOF: " << TOFObject.NTOF() << std::endl;
      for( size_t iTOF = 0; iTOF < TOFObject.NTOF(); ++iTOF ){
	fTOFHisto->Fill(TOFObject.SingleTOF(iTOF));
      }
    }
    
  }
  
  // Move LArSoft magic to save the information into the final ROOT file 
  e.put(std::move(TOFCol));
  
  return;
}





void lrm::TimeOfFlightSlicing::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;
  if( fMakeHistograms ){
    fTOFHisto = tfs->make<TH1F>("TOF","All TOF possibilities",100,0,100 );
    fTOFHisto->SetTitle("All TOF Possibilites");
    fTOFHisto->GetXaxis()->SetTitle("TOF (ns)");
    fTOFHisto->GetYaxis()->SetTitle("TOF hits per 1 ns");
  }
}

void lrm::TimeOfFlightSlicing::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::endJob()
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::reconfigure(fhicl::ParameterSet const & p)
{
   
   // Implementing the ability to pass the name of the TriggerUtililty 
   fTriggerUtility     = p.get< std::string >("TriggerUtility", "FragmentToDigit");
   fMakeHistograms     = p.get< bool >("MakeHistograms", true);
   fSlicerSourceLabel  = p.get< std::string >("SourceLabel","SlicerInput");
}

void lrm::TimeOfFlightSlicing::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlightSlicing::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(lrm::TimeOfFlightSlicing)
