////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdentification
// Module Type: producer
// File:        ParticleIdentification_module.cc
//
// Generated at Tue Jul 21 10:59:53 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

//Lariatsoft includes
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "RawDataUtilities/TriggerDigitUtility.h"

//C++ includes
#include <vector>
#include <memory>
#include <iostream>

//ROOT includes
#include <TH2F.h>
#include <TH1F.h>


class ParticleIdentification;

class ParticleIdentification : public art::EDProducer {
public:
  explicit ParticleIdentification(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleIdentification(ParticleIdentification const &) = delete;
  ParticleIdentification(ParticleIdentification &&) = delete;
  ParticleIdentification & operator = (ParticleIdentification const &) = delete;
  ParticleIdentification & operator = (ParticleIdentification &&) = delete;

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

  // Declare member data here.
  std::string       fWCTrackModuleLabel;
  std::string       fTOFModuleLabel;
  bool              fVerbose;
  bool              fPlotHistograms;

  TH2F*             fPVsTOF;
  TH1F*             fNTOF;
  TH1F*             fP;
  TH1F*             fY_Kink;
  TH1F*             fX_Dist;
  TH1F*             fY_Dist;
  TH1F*             fZ_Dist;
  TH1F*             fX_Face_Dist;
  TH1F*             fY_Face_Dist;
  TH1F*             fTheta_Dist;
  TH1F*             fPhi_Dist;
  TH1F*             fTOF;
  
};


ParticleIdentification::ParticleIdentification(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);

  
}

void ParticleIdentification::produce(art::Event & e)
{
  // Implementation of required member function here.
  


  //Get the collection of WCTracks produced by the WCTrackBuilder module
  art::Handle< std::vector<ldp::WCTrack> > WCTrackColHandle;
  e.getByLabel(fWCTrackModuleLabel,WCTrackColHandle);
  
  //Get the collection of TOF objects produced by the TOF module
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  e.getByLabel(fTOFModuleLabel,TOFColHandle);

  std::cout << "Produce has commenced." << std::endl;

  //Getting the associations between triggers-WCTracks and triggers-TOF
  art::FindManyP<raw::Trigger> findTheWCTrackTrig(WCTrackColHandle,e,fWCTrackModuleLabel);
  art::FindManyP<raw::Trigger> findTheTOFTrig(TOFColHandle,e,fTOFModuleLabel);
  
  //Extract information about the triggers and compare the WCTrackTrigs against the TOFTrigs
  if( fVerbose ){
    std::cout << "Size of WCTrackColHandle: " << WCTrackColHandle->size() << std::endl;
    for( size_t iWCTrack = 0; iWCTrack < WCTrackColHandle->size(); ++iWCTrack ){
      std::vector<art::Ptr<raw::Trigger> > trigPtrVect = findTheWCTrackTrig.at(iWCTrack);
      std::cout << "iWCTrack: " << iWCTrack << ", Size of trigPtrVect: " << trigPtrVect.size() << ", TriggerID: " << trigPtrVect.at(0)->TriggerNumber() << std::endl;
    }    
    std::cout << "Size of TOFColHandle: " << TOFColHandle->size() << std::endl;
    for( size_t iTOF = 0; iTOF < TOFColHandle->size(); ++iTOF ){
      std::vector<art::Ptr<raw::Trigger> > trigPtrVect = findTheTOFTrig.at(iTOF);
      std::cout << "iTOF: " << iTOF << ", Size of trigPtrVect: " << trigPtrVect.size() << ", TriggerID: " << trigPtrVect.at(0)->TriggerNumber() << std::endl;
    }
  }

  //Quick-and-dirty method for matching TOF and WCTracks: good WCTracks only have one track per trigger.
  //I think (?) that only one TOF comes out of the TimeOfFlight module, but just in case that's not true, we'll require
  //that only one TOF can correspond to one trigger
  
  //Loop through the WCTracks and get the trigger ID
  for( size_t iWCTrack = 0; iWCTrack < WCTrackColHandle->size(); ++iWCTrack ){
    std::vector<art::Ptr<raw::Trigger> > trigPtrVectWCTrack = findTheWCTrackTrig.at(iWCTrack);
    unsigned int triggerIDWCTrack = trigPtrVectWCTrack.at(0)->TriggerNumber();
    //Loop through the TOF and get the trigger ID
    int matchCounter = 0;
    bool uniqueMatch = false;
    ldp::WCTrack theWCTrack = WCTrackColHandle->at(iWCTrack);
    ldp::TOF theTOF;
    for( size_t iTOF = 0; iTOF < TOFColHandle->size(); ++iTOF ){
      std::vector<art::Ptr<raw::Trigger> > trigPtrVectTOF = findTheTOFTrig.at(iTOF);
      unsigned int triggerIDTOF = trigPtrVectTOF.at(0)->TriggerNumber();
      if( triggerIDWCTrack == triggerIDTOF ){
	matchCounter++;
	uniqueMatch = true;
	theTOF = TOFColHandle->at(iTOF);
      }
      if( matchCounter > 1 ){
	uniqueMatch = false;
	break;
      }
    }
    if( uniqueMatch == true ){
      //Send to histo for cut analysis
      if( fVerbose )
	std::cout << "Number of hits in TOF: " << theTOF.NTOF() << std::endl;
      if( fPlotHistograms ){
	fNTOF->Fill(theTOF.NTOF());
	fP->Fill(theWCTrack.Momentum());
	fY_Kink->Fill(theWCTrack.YKink());
	fX_Dist->Fill(theWCTrack.DeltaDist(0));
	fY_Dist->Fill(theWCTrack.DeltaDist(1));
	fZ_Dist->Fill(theWCTrack.DeltaDist(2));
	fX_Face_Dist->Fill(theWCTrack.XYFace(0));
	fY_Face_Dist->Fill(theWCTrack.XYFace(1));
	fTheta_Dist->Fill(theWCTrack.Theta());
	fPhi_Dist->Fill(theWCTrack.Phi());
	fTOF->Fill(theTOF.SingleTOF(0));
	short timeOfFlight = theTOF.SingleTOF(0);
	float p = theWCTrack.Momentum();
	fPVsTOF->Fill(p,timeOfFlight);
      }

      //Send to PIDAlg for identification
    }
  }
  
}

void ParticleIdentification::beginJob()
{
  // Implementation of optional member function here.
  if( fPlotHistograms ){
    art::ServiceHandle<art::TFileService> tfs;
    fPVsTOF = tfs->make<TH2F>("PVsTOF","P vs. Time of Flight",160,0,1600,60,20,80);
    fNTOF = tfs->make<TH1F>("NTOF","Number of TOF values per TOF object",10,0,10);
    fP = tfs->make<TH1F>("Reco_P","Reconstructed momentum",180,0,1800);
    fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",100,-5*3.1415926/180,5*3.141592654/180);
    fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",120,-60,60);
    fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",120,-60,60);
    fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",120,-60,60);
    fX_Face_Dist = tfs->make<TH1F>("X_Face","X Location of Track's TPC Entry (mm)",800,-200,600);
    fY_Face_Dist = tfs->make<TH1F>("Y_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);
    fTheta_Dist = tfs->make<TH1F>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",100,0,0.2);
    fPhi_Dist = tfs->make<TH1F>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",100,0,6.28318);
    
    fP->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
    fP->GetYaxis()->SetTitle("Tracks per 10 MeV/c");
    fY_Kink->GetXaxis()->SetTitle("Reconstructed y_kink (radians)");
    fY_Kink->GetYaxis()->SetTitle("Tracks per 0.000872 radians");
    fX_Dist->GetXaxis()->SetTitle("X distance between US and DS track ends");
    fX_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Dist->GetXaxis()->SetTitle("Y distance between US and DS track ends");
    fY_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fZ_Dist->GetXaxis()->SetTitle("Z distance between US and DS track ends");
    fZ_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fX_Face_Dist->GetXaxis()->SetTitle("X (mm)");
    fX_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Face_Dist->GetXaxis()->SetTitle("Y (mm)");
    fY_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fTheta_Dist->GetXaxis()->SetTitle("Theta (radians)");
    fTheta_Dist->GetYaxis()->SetTitle("Tracks per .002 radians");
    fPhi_Dist->GetXaxis()->SetTitle("Phi (radians)");
    fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.0628 radians");

    fTOF = tfs->make<TH1F>("Reco_TOF","Reconstructed Time of Flight",70,20,90);
    


  }
}

void ParticleIdentification::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::endJob()
{
  // Implementation of optional member function here.
}

void ParticleIdentification::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fWCTrackModuleLabel           =p.get< std::string >("WCTrackModuleLabel");
  fTOFModuleLabel               =p.get< std::string >("TOFModuleLabel");
  fVerbose                      =p.get< bool >("Verbose");
  fPlotHistograms               =p.get< bool >("PlotHistograms");
}

void ParticleIdentification::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentification::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleIdentification)
