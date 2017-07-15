////////////////////////////////////////////////////////////////////////
// Class:       NNShowerFilter
// Module Type: filter
// File:        NNShowerFilter_module.cc
//
// Daniel Smith
// Boston University
// dansmith@bu.edu
//
// Module to filter out an EM shower that is matched to the projected WC
//  using a NN classifier that classifies hits as,
//  1) Track 2) Shower 3) Nothing in a data event. 
//
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"
#include "lardata/ArtDataHelper/MVAWriter.h"

#include "tbb/parallel_for.h"
#include "tbb/tbb.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larevt/Filters/ChannelFilter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "RawDataUtilities/TriggerDigitUtility.h"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/AuxDetParticleID.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/CoreUtils/ServiceUtil.h" 

// LArSoft includes
#include "larcorealg/Geometry/ChannelMapAlg.h" 
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT & C++
#include <memory>
#include "TCanvas.h"
#include <TH2F.h>
#include <TH1F.h>

#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include <sys/stat.h>
#include <iostream>

using namespace keras;

class NNShowerFilter;

class NNShowerFilter : public art::EDFilter {
public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg { Name("PointIdAlg") };

    fhicl::Atom<art::InputTag> WireLabel { Name("WireLabel"),      
	Comment("tag of deconvoluted ADC on wires (recob::Wire)") };

    fhicl::Atom<art::InputTag> HitModuleLabel { Name("HitModuleLabel"),      
	Comment("tag of hits to be EM/track tagged") };	

    fhicl::Atom<art::InputTag> TrackModuleLabel { Name("TrackModuleLabel"),      
	Comment("tag of 3D tracks to be EM/track tagged using accumulated results from hits in the best 2D projection") };

    fhicl::Atom<art::InputTag> WC2TPCModuleLabel { Name("WC2TPCModuleLabel"),      
	Comment("WC to TPC track Matching") };

    fhicl::Atom<art::InputTag> WCTrackLabel { Name("WCTrackLabel"),      
	Comment("wc building alg") };

    fhicl::Atom<float> ShowerThresh { Name("ShowerThresh"),      
	Comment("Threshold used to determine if a shower or not. ") };
  
  };
  
  using Parameters = art::EDProducer::Table<Config>;

  explicit NNShowerFilter(Parameters const & p);

  //  explicit NNShowerFilter(fhicl::ParameterSet const & p);

  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NNShowerFilter(NNShowerFilter const &) = delete;
  NNShowerFilter(NNShowerFilter &&) = delete;
  NNShowerFilter & operator = (NNShowerFilter const &) = delete;
  NNShowerFilter & operator = (NNShowerFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TH2F* fMapTrack;
  TH2F* fMapShower;
  TH2F* fMapNothing;

  TH1F* fProb0;
  TH1F* fProb1;
  TH1F* fProb2;

  TH1F* fProbHitTrack;
  TH1F* fProbHitShower;
  TH1F* fProbHitNothing;

  nnet::PointIdAlg fPointIdAlg;
  art::InputTag fWireProducerLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fTrackModuleLabel;
  art::InputTag fWC2TPCModuleLabel;
  art::InputTag fWCTrackLabel;
  float fShowerThresh;

};


//NNShowerFilter::NNShowerFilter(fhicl::ParameterSet const & p)
NNShowerFilter::NNShowerFilter(NNShowerFilter::Parameters const & config)
  : fPointIdAlg(config().PointIdAlg()),
    fWireProducerLabel(config().WireLabel()),
    fHitModuleLabel(config().HitModuleLabel()),
    fTrackModuleLabel(config().TrackModuleLabel()),
    fWC2TPCModuleLabel(config().WC2TPCModuleLabel()),
    fWCTrackLabel(config().WCTrackLabel()),
    fShowerThresh(config().ShowerThresh())
{
  // Call appropriate produces<>() functions here.
}

bool NNShowerFilter::filter(art::Event & event) { 

  // WC info
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
   
  if(event.getByLabel(fWCTrackLabel, wctrackHandle))
    art::fill_ptr_vector(wctrack, wctrackHandle);

  // Wire info
  art::Handle< std::vector<recob::Wire> > WireHandle;
  std::vector< art::Ptr<recob::Wire> > Wirelist;

  if(event.getByLabel(fWireProducerLabel, WireHandle))
    art::fill_ptr_vector(Wirelist, WireHandle);

  // Track info
  art::Handle< std::vector<recob::Track> > TrackHandle;
  std::vector< art::Ptr<recob::Track> > Tracklist;

  if(event.getByLabel(fTrackModuleLabel, TrackHandle))
    art::fill_ptr_vector(Tracklist, TrackHandle);

  // Hit info
  art::Handle< std::vector<recob::Hit> > HitHandle;
  std::vector< art::Ptr<recob::Hit> > Hitlist;

  if(event.getByLabel(fHitModuleLabel, HitHandle))
    art::fill_ptr_vector(Hitlist, HitHandle);

  // Association between Tracks and 2d Hits
  art::FindManyP<recob::Track> ass_trk_hits(HitHandle,   event, fTrackModuleLabel);

  // Association between Tracks and WC event
  art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, event, fWC2TPCModuleLabel);

  // Requiring one WC-TPC match
  if(fWC2TPC.size() != 1) { return false; }

  std::vector< float > probs (3, 0.);
  float numHits = 0.0;

  for(size_t i = 0; i < fWC2TPC.size(); i++) {

    // Get the TPC track
    cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(i));
    if(!trackWC2TPC) { return false; }  // Reject if no matches
    recob::Track const& matchedTrack(trackWC2TPC.ref()); 

    // Set up the classifier
    fPointIdAlg.setWireDriftData(*WireHandle, 1, 0, 0);
    
    for(size_t iHit = 0; iHit < Hitlist.size(); ++iHit) {

      if(ass_trk_hits.at(iHit).size() == 0) { continue; } // Make sure hit has an assocated track
      if(ass_trk_hits.at(iHit)[0]->ID() != matchedTrack.ID()) { continue; } // Make sure Hit is from WCtrack track
      if(Hitlist[iHit]->View() != 1) { continue; } // All training is from collection view

      int wireID = Hitlist[iHit]->WireID().Wire;
      int hitTime = Hitlist[iHit]->PeakTime();

      std::vector<float> results = fPointIdAlg.predictIdVector(wireID, hitTime);
   
      numHits += 1.0;
      probs[0] += results[0];
      probs[1] += results[1];
      probs[2] += results[2];

      fProbHitTrack->Fill(results[0]);
      fProbHitShower->Fill(results[1]);
      fProbHitNothing->Fill(results[2]);

      fMapTrack->Fill(wireID, hitTime, results[0]);
      fMapShower->Fill(wireID, hitTime, results[1]);
      fMapNothing->Fill(wireID, hitTime, results[2]);

    } // End of hitlist.size

  }

  // Plotting average NN outputs for each hit in a track
  /*
  std::cout << "results for matched track:" << std::endl;
  for(size_t j = 0; j < probs.size(); j++) { std::cout << probs[j] / float(numHits) << " "; }
  std::cout << "\n\n";
  */

  fProb0->Fill(probs[0] / numHits);
  fProb1->Fill(probs[1] / numHits);
  fProb2->Fill(probs[2] / numHits);
  
  return ((probs[1] / numHits) > fShowerThresh);

}

void NNShowerFilter::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fProb0 = tfs->make<TH1F>("fProb0","fProb0",50, 0.0, 1.0);
  fProb1 = tfs->make<TH1F>("fProb1","fProb1",50, 0.0, 1.0);
  fProb2 = tfs->make<TH1F>("fProb2","fProb2",50, 0.0, 1.0);

  fProbHitTrack = tfs->make<TH1F>("fProbHitTrack","fProbHitTrack",50, 0.0, 1.0);
  fProbHitShower = tfs->make<TH1F>("fProbHitShower","fProbHitShower",50, 0.0, 1.0);
  fProbHitNothing = tfs->make<TH1F>("fProbHitNothing","fProbHitNothing",50, 0.0, 1.0);

  fMapTrack = tfs->make<TH2F>("fMapTrack","fMapTrack",500, 0.0, 500.0, 300, 0.0, 3000.); 
  fMapShower = tfs->make<TH2F>("fMapShower","fMapShower",500, 0.0, 500.0, 300, 0.0, 3000.); 
  fMapNothing = tfs->make<TH2F>("fMapNothing","fMapNothing",500, 0.0, 500.0, 300, 0.0, 3000.); 

  fMapTrack->GetZaxis()->SetRangeUser(0.0, 1.0);
  fMapShower->GetZaxis()->SetRangeUser(0.0, 1.0);
  fMapNothing->GetZaxis()->SetRangeUser(0.0, 1.0);

}

void NNShowerFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NNShowerFilter)
