////////////////////////////////////////////////////////////////////////
// Class:       DoubleFlashFilter
// Module Type: filter
// File:        DoubleFlashFilter_module.cc
//
// This filter looks for delayed double pulses in the optical PMT 
// waveforms in beam events and passes event only if all of these apply:
//   (1) PMT contains 2 hits
//   (2) First hit occurs at the appropriate trigger location
//   (3) Second hit is spaced between T1 and T2 from the first
// 
// This could be used in conjunction with other filters, like for example,
// a stopping beam-particle filter, to look for muons that stop and decay.
// Could also be used to create a more complicated filter to tag pi -> mu 
// -> e. 
//
// Generated at Wed Apr  5 16:59:00 2017 by William Foreman using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// ROOT includes
#include <TH1F.h>
#include <TH2F.h>

// C++ includes
#include <iostream>
#include <memory>

// LArIATSoft includes
#include "LArIATRecoAlg/OpHitBuilderAlg.h"

class DoubleFlashFilter;

class DoubleFlashFilter : public art::EDFilter {
public:
  explicit DoubleFlashFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DoubleFlashFilter(DoubleFlashFilter const &) = delete;
  DoubleFlashFilter(DoubleFlashFilter &&) = delete;
  DoubleFlashFilter & operator = (DoubleFlashFilter const &) = delete;
  DoubleFlashFilter & operator = (DoubleFlashFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  // Optical detector alg
  OpHitBuilderAlg     fOpHitBuilderAlg;
  
  // Declare member data here.
  std::string       fDAQModuleLabel;
  std::string       fDAQModuleInstanceName;
  size_t            fOpDetChannel;

  short             fT1;
  short             fT2;
  short             fFirstPulseTolerance;

  TH1F* hTotalBeamEvents;
  TH1F* hNumOpHits;
  TH1F* hTimeDifference;
  TH1F* hTimeDifference_pass;

};


DoubleFlashFilter::DoubleFlashFilter(fhicl::ParameterSet const & p)
:
fOpHitBuilderAlg(p)
{
  this->reconfigure(p);
  art::ServiceHandle<art::TFileService> tfs;
  hTotalBeamEvents          = tfs->make<TH1F>("TotalBeamEvents","",1,0,1);
  hNumOpHits            = tfs->make<TH1F>("NumOpHits",";Number of optical hits found per event",20,0,20);
  hTimeDifference       = tfs->make<TH1F>("TimeDifference",";#Delta T for events w/ 2 optical hits;#Delta T [ns]",200,0.,20000);
  hTimeDifference_pass  = tfs->make<TH1F>("TimeDifference_pass",";#Delta T for events w/ 2 optical hits (passing filter);#Delta T [ns]",200,0.,20000);
}

bool DoubleFlashFilter::filter(art::Event & e)
{

  std::vector<short> waveform;
  raw::OpDetPulse opdetpulse;
  float timestamp=0.;

  // Get OpDetPulses
  art::Handle< std::vector< raw::OpDetPulse >> opdetHandle;
  e.getByLabel(fDAQModuleLabel, fDAQModuleInstanceName, opdetHandle);
  
  // Check size of OpDetPulse handle  
  if( (size_t)opdetHandle->size() == 0 ){
    std::cout << "No optical detector data found; skipping event.\n";
    return false;
  }

  // Get the waveform for the channel specificed by fOpDetChannel (and timestamp)
  for (size_t i_pulse = 0; i_pulse < opdetHandle->size(); i_pulse ++ ){
    art::Ptr< raw::OpDetPulse > ThePulsePtr(opdetHandle,i_pulse); 
    raw::OpDetPulse pulse = *ThePulsePtr;
    if( pulse.OpChannel() == fOpDetChannel ) {
      waveform = pulse.Waveform();    
      timestamp = ((float)pulse.PMTFrame()*8.)/1.0e09;
      opdetpulse = pulse;
      break;
    }
  }

  // If this is not a beam event, skip!
  if( timestamp < 1.2 || timestamp >= 5.2 ) return false;
  
  hTotalBeamEvents->Fill(0);
  
  std::vector<short> hitTimes = fOpHitBuilderAlg.GetHits( opdetpulse );

  size_t nHits = hitTimes.size();
  hNumOpHits->Fill(nHits);

  std::cout<<"Found number optical hits: "<<nHits<<"\n";

  if( nHits != 2 ) return false;


  short hit1 = hitTimes[0];
  short hit2 = hitTimes[1];
  short dt    = hit2 - hit1;
  short firstSample = opdetpulse.FirstSample();
  hTimeDifference->Fill(dt);
  
  std::cout<<"  "<<hit1<<"   "<<hit2<<"  ("<<dt<<")    first sample = "<<firstSample<<"\n";

  // Require first hit occurs near the "trigger time"
  // for c2: call to abs is ambiguous
  //if( abs(hit1 - opdetpulse.FirstSample()) >= 200 ) return false;
  if( util::absDiff(hit1, (short)opdetpulse.FirstSample()) >= 200 ) return false;
  
  // Require second hit be delayed by appropriate amount
  if( dt >= fT1 && dt <= fT2 ) {
    std::cout<<"Passes cut!\n";
    hTimeDifference_pass->Fill(dt);
    return true;
  }

  return false;
}

void DoubleFlashFilter::beginJob()
{
}

void DoubleFlashFilter::endJob()
{
}

void DoubleFlashFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fDAQModuleLabel               = p.get< std::string >  ("DAQModule","daq");
  fDAQModuleInstanceName        = p.get< std::string >  ("DAQInstanceName","");
  fOpDetChannel                 = p.get< size_t >       ("OpDetChannel",1);
  fT1                           = p.get< short >        ("T1",500);
  fT2                           = p.get< short >        ("T2",10000);
  fFirstPulseTolerance          = p.get< short >        ("FirstPulseTolerance",200);
  
  
}


DEFINE_ART_MODULE(DoubleFlashFilter)
