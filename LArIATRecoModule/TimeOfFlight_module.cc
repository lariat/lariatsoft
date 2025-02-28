/*

LArIATSoft 

Time Of Flight producer module

Module: TimeOfFlight_module.cc 
FHiCL: timeOfFlight.fcl
Header: LArIATDataProducts/TOF.h
Class: LArIATDataProducts/TOF.cxx
Dictionary: LArIATDataProducts/classes.h and classes_def.xml

Authores:
Daniel Smith     - dsmith@fnal.gov
Elena Gramellini - elena.gramellini@yale.edu
Irene Nutini     - irene.nutini@stud.unifi.it

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
  TOFBuilderAlg fTOFAlg;

  TH1F*         fTOFHisto;
  TH1F*         fTimeStampHisto;
  TH1F*         fNTOFHisto;

  bool          fMakeHistograms;
  std::string   fSlicerSourceLabel;

  int fSubRun;
  int fRun;

};

//----------------------------------------------------------------
lrm::TimeOfFlight::TimeOfFlight(fhicl::ParameterSet const & p)
  : EDProducer(p),
    fTOFAlg(p)
{

  // Configures the ROOT histograms and the 
  this->reconfigure(p);  

  // Produces the LArSoft object to be ultimately outputted
  produces<std::vector<ldp::TOF> >();
}

void lrm::TimeOfFlight::produce(art::Event & e)
{

  //Clearing the persistent members so we don't have repetition
  fTOFAlg.clear_tof_and_timeStampDst();
 
  // Creating the object for the TOF, the object that will be outputted
  std::unique_ptr<std::vector<ldp::TOF> > TOFCol(new std::vector<ldp::TOF> );

  //Retrieving the digits for the upstream and downstream paddles from the sliced event
  art::Handle< std::vector<raw::AuxDetDigit> > AuxDetDigitHandle;
  e.getByLabel(fSlicerSourceLabel,AuxDetDigitHandle);

  //Get the digits for USTOF and DSTOF
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
    
    std::pair<std::vector<float>, std::vector<long> > pair;
    pair = fTOFAlg.get_TOF_and_TimeStamp(USTOF, DSTOF);

    ldp::TOF TOFObject(pair.first, pair.second);
    (*TOFCol).push_back(TOFObject);

    //Fill histos if necessary
    if( fMakeHistograms ){
      fNTOFHisto->Fill(TOFObject.NTOF());
      
      for( size_t iTOF = 0; iTOF < TOFObject.NTOF(); ++iTOF ){
	fTOFHisto->Fill(TOFObject.SingleTOF(iTOF));
	fTimeStampHisto->Fill(TOFObject.TimeStamp(iTOF));
      }
    }
    
  }
  
  // Move LArSoft magic to save the information into the final ROOT file 
  e.put(std::move(TOFCol));  
  return;
}

void lrm::TimeOfFlight::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;

  if( fMakeHistograms ){
    fTOFHisto = tfs->make<TH1F>("TOF","All TOF possibilities",500,0.0,100.0);
    fTOFHisto->SetTitle("All TOF Possibilites");
    fTOFHisto->GetXaxis()->SetTitle("TOF (ns)");
    fTOFHisto->GetYaxis()->SetTitle("TOF hits per 1 ns");

    fTimeStampHisto = tfs->make<TH1F>("TimeStamp","All Timestamps for TOFs",1000,0.0,9000000000.0);
    fNTOFHisto = tfs->make<TH1F>("NTOF","Number of TOF in an event",10,0.0,10.0);
  }
}

void lrm::TimeOfFlight::beginRun(art::Run & r)
{
  fRun = r.run();
}

void lrm::TimeOfFlight::beginSubRun(art::SubRun & sr)
{
  fSubRun = sr.subRun();
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
   fMakeHistograms     = p.get< bool >("MakeHistograms", true);
   fSlicerSourceLabel  = p.get< std::string >("SourceLabel");
   std::cout<<"TOF Label: "<<fSlicerSourceLabel<<std::endl;
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
