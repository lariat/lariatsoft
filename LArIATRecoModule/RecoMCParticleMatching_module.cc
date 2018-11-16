//////////////////////////////////////////////////////////////
// Name:      RecoMCParticleMatching
// Date:      5 October 2017
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef RecoMCParticleMatching_Module
#define RecoMCParticleMatching_Module

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// NuTools includes
#include "nusimdata/SimulationBase/MCParticle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"

// LArIATSoft includes
#include "NeutralPionAnalysis/SimulationUtil.h"

// ROOT includes
#include "TTree.h"

// C++ includes
//#include <algorithm>
#include <cmath>
//#include <iterator>
//#include <limits>
//#include <map>
#include <memory>
#include <string>
#include <vector>

struct metric_t {
  double score        = -1;
  double cleanliness  = -1;
  double completeness = -1;
  int    track_id     = -1;
};

bool score_sort(metric_t const& x, metric_t const& y)
{
  return x.score > y.score;
}

namespace npa
{
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class RecoMCParticleMatching : public art::EDProducer
  {

   public:

    // standard constructor and destructor for an art module
    explicit RecoMCParticleMatching(fhicl::ParameterSet const& pset);
    virtual ~RecoMCParticleMatching();

    // this method is called once, at the start of the job
    void beginJob() override;

    // this method is called once, at the start of each run
    void beginRun(art::Run & run) override;

    // this method is called once, at the start of each subrun
    void beginSubRun(art::SubRun & subrun) override;

    // this method is called once, at the end of the job
    void endJob() override;

    // this method is called once, at the end of each run
    void endRun(art::Run & run) override;

    // this method is called once, at the end of each subrun
    void endSubRun(art::SubRun & subrun) override;

    // this method reads in any parameters from the .fcl files
    void reconfigure(fhicl::ParameterSet const& pset) override;

    // the produce routine, called once per event
    void produce(art::Event & event) override;

   private:

    // parameters read from FHiCL (.fcl) file
    std::string simulation_producer_label_;
    std::string mctrack_producer_label_;
    std::string mcshower_producer_label_;
    std::string hit_producer_label_;
    std::string cluster_producer_label_;
    std::string track_producer_label_;
    std::string shower_producer_label_;
    bool        match_tracks_;
    bool        match_showers_;
    bool        debug_;

    // pointer to geometry provider
    geo::GeometryCore const* geometry_;

    // pointer to detector properties
    detinfo::DetectorProperties const* detector_;

    // back tracker service
    art::ServiceHandle< cheat::BackTrackerService > backtracker_;

    // simulation util
    simutil::SimulationUtil simutil_;

    // pointers to TTree object
    TTree * tree_;

    // variables that will go into the TTree objects
    int event_;     // number of the event being processed
    int run_;       // number of the run being processed
    int subrun_;    // number of the sub-run being processed

    std::vector< double > score_;
    std::vector< double > cleanliness_;
    std::vector< double > completeness_;
    std::vector< int >    pdg_;
    std::vector< int >    number_hits_;

    // fill metrics vector
    void fill_metrics_vector_(
        art::ValidHandle< std::vector< simb::MCParticle > > const& particle_handle,
        art::ValidHandle< std::vector< sim::MCTrack > >     const& mctrack_handle,
        art::ValidHandle< std::vector< sim::MCShower > >    const& mcshower_handle,
        std::set< int >                                     const& g4_trk_ids_hits,
        std::vector< art::Ptr< recob::Hit > >               const& hits,
        std::vector< art::Ptr< recob::Hit > >               const& hit_vector,
        simutil::ParticleMap                                const& particle_map,
        simutil::ShowerParticleMaps                         const& shower_particle_maps,
        std::vector< metric_t >                                  & metrics);

    // reset per event
    void reset_();
  };

  //-----------------------------------------------------------------------
  // constructor
  RecoMCParticleMatching::RecoMCParticleMatching(fhicl::ParameterSet const& pset)
  //  : EDProducer(pset)
      : simutil_(pset.get< fhicl::ParameterSet >("SimulationUtil"))
  {
    // reconfigure parameters
    this->reconfigure(pset);

    // get a pointer to the geometry service provider
    geometry_ = &*(art::ServiceHandle<geo::Geometry>());

    // get a pointer to the detector properties provider
    detector_ = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // produce associations
    //produces< art::Assns< recob::Cluster, simb::MCParticle, anab::BackTrackerMatchingData > > ();
    produces< art::Assns< recob::Track, simb::MCParticle, anab::BackTrackerMatchingData > > ();
    produces< art::Assns< recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData > > ();
  }

  //-----------------------------------------------------------------------
  // destructor
  RecoMCParticleMatching::~RecoMCParticleMatching() {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::beginJob()
  {
    // access art's TFileService
    art::ServiceHandle< art::TFileService > tfs;

    tree_ = tfs->make<TTree>(
        "RecoMCParticleMatching", "RecoMCParticleMatching");

    tree_->Branch("event",  &event_,  "event/I");
    tree_->Branch("run",    &run_,    "run/I");
    tree_->Branch("subrun", &subrun_, "subrun/I");

    tree_->Branch("score",        &score_);
    tree_->Branch("cleanliness",  &cleanliness_);
    tree_->Branch("completeness", &completeness_);
    tree_->Branch("pdg",          &pdg_);
    tree_->Branch("number_hits",  &number_hits_);
  }

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::beginRun(art::Run & /*run*/)
  {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::beginSubRun(art::SubRun & /*subrun*/)
  {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::endJob()
  {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::endRun(art::Run & /*run*/)
  {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::endSubRun(art::SubRun & /*subrun*/)
  {}

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::reconfigure(fhicl::ParameterSet const& pset)
  {
    // read parameters from the .fcl file
    simulation_producer_label_ = pset.get< std::string >("simulation_label", "largeant");
    mctrack_producer_label_    = pset.get< std::string >("mctrack_label",    "mcreco");
    mcshower_producer_label_   = pset.get< std::string >("mcshower_label",   "mcreco");
    hit_producer_label_        = pset.get< std::string >("hit_label");
    cluster_producer_label_    = pset.get< std::string >("cluster_label");
    track_producer_label_      = pset.get< std::string >("track_label");
    shower_producer_label_     = pset.get< std::string >("shower_label");
    match_tracks_              = pset.get< bool        >("match_tracks");
    match_showers_             = pset.get< bool        >("match_showers");
    debug_                     = pset.get< bool        >("debug");
  }

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::produce(art::Event & event)
  {
    //-------------------------------------------------------------------
    // associations
    //-------------------------------------------------------------------
    //std::unique_ptr< art::Assns< recob::Cluster, simb::MCParticle, anab::BackTrackerMatchingData > >
    //    cluster_mcparticle_assn(new art::Assns< recob::Cluster, simb::MCParticle, anab::BackTrackerMatchingData >);

    std::unique_ptr< art::Assns< recob::Track, simb::MCParticle, anab::BackTrackerMatchingData > >
        track_mcparticle_assn(new art::Assns< recob::Track, simb::MCParticle, anab::BackTrackerMatchingData >);

    std::unique_ptr< art::Assns< recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData > >
        shower_mcparticle_assn(new art::Assns< recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData >);

    //-------------------------------------------------------------------
    // reset once per event
    //-------------------------------------------------------------------

    this->reset_();

    //-------------------------------------------------------------------
    // get event, run, and subrun numbers
    //-------------------------------------------------------------------

    event_  = event.id().event();
    run_    = event.run();
    subrun_ = event.subRun();

    //-------------------------------------------------------------------
    // get simulated particles and simulated channels
    //-------------------------------------------------------------------

    // get all the simulated particles in the event
    // art::ValidHandle< std::vector< simb::MCParticle > >
    auto particle_handle = event.getValidHandle< std::vector< simb::MCParticle > >
        (simulation_producer_label_);

    // get all the simulated channels in the event
    // art::ValidHandle< std::vector< sim::SimChannel > >
    auto sim_channel_handle = event.getValidHandle< std::vector< sim::SimChannel > >
        (simulation_producer_label_);

    // get all the MCTracks in the event
    // art::ValidHandle< std::vector< sim::MCTrack > >
    auto mctrack_handle = event.getValidHandle< std::vector< sim::MCTrack > >
        (mctrack_producer_label_);

    // get all the MCShowers in the event
    // art::ValidHandle< std::vector< sim::MCShower > >
    auto mcshower_handle = event.getValidHandle< std::vector< sim::MCShower > >
        (mcshower_producer_label_);

    // fill vector of simulated particles
    std::vector< art::Ptr< simb::MCParticle > > particle_vector;
    art::fill_ptr_vector(particle_vector, particle_handle);

    // fill vector of simulated channels
    std::vector< art::Ptr< sim::SimChannel > > sim_channel_vector;
    art::fill_ptr_vector(sim_channel_vector, sim_channel_handle);

    // fill vector of MCTracks
    std::vector< art::Ptr< sim::MCTrack > > mctrack_vector;
    art::fill_ptr_vector(mctrack_vector, mctrack_handle);

    // fill vector of MCShowers
    std::vector< art::Ptr< sim::MCShower > > mcshower_vector;
    art::fill_ptr_vector(mcshower_vector, mcshower_handle);

    //-------------------------------------------------------------------
    // get hits, clusters, tracks, and vertices
    //-------------------------------------------------------------------

    // get all the hits in the event
    // art::ValidHandle< std::vector< recob::Hit > >
    auto hit_handle = event.getValidHandle< std::vector< recob::Hit > >
        (hit_producer_label_);

    // get all the clusters in the event
    // art::ValidHandle< std::vector< recob::Cluster > >
    auto cluster_handle = event.getValidHandle< std::vector< recob::Cluster > >
        (cluster_producer_label_);

    // get all the reconstructed tracks in the event
    // art::ValidHandle< std::vector< recob::Track > >
    auto track_handle = event.getValidHandle< std::vector< recob::Track > >
        (track_producer_label_);

    // fill vector of hits
    std::vector< art::Ptr< recob::Hit > > hit_vector;
    art::fill_ptr_vector(hit_vector, hit_handle);

    // fill vector of tracks
    std::vector< art::Ptr< recob::Track > > track_vector;
    art::fill_ptr_vector(track_vector, track_handle);

    // find many hits from tracks
    const art::FindManyP< recob::Hit >
        find_many_hits_from_trks(track_handle, event, track_producer_label_);

    //-------------------------------------------------------------------
    // get particle maps
    //-------------------------------------------------------------------
    auto const particle_map = simutil_.particle_map(particle_vector);
    auto const shower_particle_maps = simutil_.shower_particle_maps(
        particle_map, mcshower_vector);

    //-------------------------------------------------------------------
    // loop over reconstructed tracks
    //-------------------------------------------------------------------
    for (size_t trk_idx = 0; trk_idx < track_vector.size(); ++trk_idx)
    {
      if (!match_tracks_) break;

      // get reco. track ID
      //int const reco_trk_id = track_vector[trk_idx]->ID();

      // get hits from reco. track
      std::vector< art::Ptr< recob::Hit > > const& hits =
          find_many_hits_from_trks.at(trk_idx);

      // get set of track IDs from hits
      std::set< int > g4_trk_ids_hits = backtracker_->GetSetOfTrackIds(hits);

      // vector for matching metrics
      std::vector< metric_t > metrics;

      this->fill_metrics_vector_(
          particle_handle, mctrack_handle, mcshower_handle, g4_trk_ids_hits,
          hits, hit_vector, particle_map, shower_particle_maps, metrics);

      //// set of G4 track IDs that have been used
      //std::set< int > g4_trk_id_set;

      //// loop through MCTrack objects
      //for (auto const& mctrack : *mctrack_handle)
      //{
      //  if (!mctrack.size()) continue;

      //  // get G4 track ID
      //  int const g4_trk_id = mctrack.TrackID();

      //  // insert G4 track ID into set for BackTracker
      //  std::set< int > id;
      //  id.insert(g4_trk_id);
      //  g4_trk_id_set.insert(g4_trk_id);

      //  // continue if G4 track ID is not found in hits
      //  if (g4_trk_ids_hits.find(g4_trk_id) == g4_trk_ids_hits.end()) continue;

      //  // get charge-weighted cleanliness and charge-weighted
      //  // completeness using the BackTracker
      //  double const cleanliness =
      //      backtracker_->HitChargeCollectionPurity(id, hits);
      //  double const completeness =
      //      backtracker_->HitChargeCollectionEfficiency(
      //      id, hits, hit_vector, geo::k3D);

      //  // compute score
      //  double const score = cleanliness * completeness;

      //  // metric
      //  metric_t metric;
      //  metric.score = score;
      //  metric.cleanliness  = cleanliness;
      //  metric.completeness = completeness;
      //  metric.track_id = g4_trk_id;

      //  metrics.push_back(metric);
      //}

      //// loop through MCShower objects
      //for (auto const& mcshower : *mcshower_handle)
      //{
      //  // get G4 track ID
      //  int g4_trk_id = mcshower.TrackID();
      //  g4_trk_id_set.insert(g4_trk_id);

      //  // continue if track ID is not found in map
      //  if (shower_particle_maps.count(g4_trk_id) < 1) continue;

      //  // continue if no energy was deposited
      //  // FIXME: i don't think this works...
      //  //if (mcshower.DetProfile().E() < 1e-5) continue;

      //  // insert G4 track IDs into set for BackTracker
      //  std::set< int > id;

      //  // number of particles found in hits
      //  unsigned int number_particles = 0;

      //  for (auto const& particle_iter : shower_particle_maps.at(g4_trk_id))
      //  {
      //    id.insert(particle_iter.first);
      //    g4_trk_id_set.insert(g4_trk_id);

      //    // increase counter if G4 track ID is found in hits
      //    if (g4_trk_ids_hits.find(g4_trk_id) != g4_trk_ids_hits.end())
      //        ++number_particles;
      //  }

      //  // continue if G4 track IDs are not found in hits
      //  if (number_particles < 1) continue;

      //  // get charge-weighted cleanliness and charge-weighted
      //  // completeness using the BackTracker
      //  double const cleanliness =
      //      backtracker_->HitChargeCollectionPurity(id, hits);
      //  double const completeness =
      //      backtracker_->HitChargeCollectionEfficiency(
      //      id, hits, hit_vector, geo::k3D);

      //  // compute score
      //  double const score = cleanliness * completeness;

      //  // metric
      //  metric_t metric;
      //  metric.score = score;
      //  metric.cleanliness  = cleanliness;
      //  metric.completeness = completeness;
      //  metric.track_id = g4_trk_id;

      //  metrics.push_back(metric);
      //}

      //// loop through MCParticles
      //for (auto const& particle : *particle_handle)
      //{
      //  // get G4 track ID of particle
      //  int const g4_trk_id = particle.TrackId();

      //  // continue if G4 track ID has already been used
      //  if (g4_trk_id_set.find(g4_trk_id) != g4_trk_id_set.end()) continue;

      //  // continue of G4 track ID is not found in hits
      //  if (g4_trk_ids_hits.find(g4_trk_id) == g4_trk_ids_hits.end()) continue;

      //  // insert G4 track ID into set for BackTracker
      //  std::set< int > id;
      //  id.insert(g4_trk_id);
      //  g4_trk_id_set.insert(g4_trk_id);

      //  // get charge-weighted cleanliness and charge-weighted
      //  // completeness using the BackTracker
      //  double const cleanliness =
      //      backtracker_->HitChargeCollectionPurity(id, hits);
      //  double const completeness =
      //      backtracker_->HitChargeCollectionEfficiency(
      //      id, hits, hit_vector, geo::k3D);

      //  // compute score
      //  double const score = cleanliness * completeness;

      //  // metric
      //  metric_t metric;
      //  metric.score = score;
      //  metric.cleanliness  = cleanliness;
      //  metric.completeness = completeness;
      //  metric.track_id = g4_trk_id;

      //  metrics.push_back(metric);
      //}

      // if no matches, continue to the next reco. track
      if (!metrics.size()) continue;

      // sort metrics by score
      std::sort(metrics.begin(), metrics.end(), score_sort);

      // get best score
      double score = metrics.front().score;
      double cleanliness = metrics.front().cleanliness;
      double completeness = metrics.front().completeness;

      // continue if score is below threshold
      //if (score < 0.1) continue;

      // get best matched particle
      int trk_id = metrics.front().track_id;

      //if (score < 0.1)
      //{
      //  score = 0;
      //  cleanliness = 0;
      //  completeness = 0;
      //  //trk_id = -1;
      //  //continue;
      //}

      //const simb::MCParticle *tmp_particle =
      //    backtracker_->TrackIDToParticle(trk_id);

      //// continue to next reco. track if backtracker cannot find 
      //if (!tmp_particle) continue;

      art::Ptr< simb::MCParticle > particle = particle_map.at(trk_id);

      // set meta data
      anab::BackTrackerMatchingData bt_data;
      bt_data.cleanliness = cleanliness;
      bt_data.completeness = completeness;

      // add association between reco. track and simulated particle
      track_mcparticle_assn->addSingle(track_vector[trk_idx], particle, bt_data);

      score_.push_back(score);
      cleanliness_.push_back(cleanliness);
      completeness_.push_back(completeness);
      pdg_.push_back(particle->PdgCode());
      number_hits_.push_back(hits.size());
    }

    if (match_showers_)
    {
      // get all the reconstructed showers in the event
      // art::ValidHandle< std::vector< recob::Shower > >
      //auto shower_handle = event.getValidHandle< std::vector< recob::Shower > >
      //    (shower_producer_label_);

      art::Handle< std::vector< recob::Shower > > shower_handle;

      // fill vector of shower
      std::vector< art::Ptr< recob::Shower > > shower_vector;
      if (event.getByLabel(shower_producer_label_, shower_handle))
      {
        art::fill_ptr_vector(shower_vector, shower_handle);
      }

      // find many hits from showers
      const art::FindManyP< recob::Hit >
          find_many_hits_from_shwrs(shower_handle, event, shower_producer_label_);

      //-------------------------------------------------------------------
      // loop over reconstructed showers
      //-------------------------------------------------------------------
      for (size_t shwr_idx = 0; shwr_idx < shower_vector.size(); ++shwr_idx)
      {
        // get hits from reco. shower
        std::vector< art::Ptr< recob::Hit > > const& hits =
            find_many_hits_from_shwrs.at(shwr_idx);

        // get set of track IDs from hits
        std::set< int > g4_trk_ids_hits = backtracker_->GetSetOfTrackIds(hits);

        // vector for matching metrics
        std::vector< metric_t > metrics;

        this->fill_metrics_vector_(
            particle_handle, mctrack_handle, mcshower_handle, g4_trk_ids_hits,
            hits, hit_vector, particle_map, shower_particle_maps, metrics);

        // if no matches, continue to the next reco. track
        if (!metrics.size()) continue;

        // sort metrics by score
        std::sort(metrics.begin(), metrics.end(), score_sort);

        // get best score
        double score = metrics.front().score;
        double cleanliness = metrics.front().cleanliness;
        double completeness = metrics.front().completeness;

        // continue if score is below threshold
        //if (score < 0.1) continue;

        // get best matched particle
        int trk_id = metrics.front().track_id;

        //if (score < 0.1)
        //{
        //  score = 0;
        //  cleanliness = 0;
        //  completeness = 0;
        //  //trk_id = -1;
        //  //continue;
        //}

        //const simb::MCParticle *tmp_particle =
        //    backtracker_->TrackIDToParticle(trk_id);

        //// continue to next reco. track if backtracker cannot find 
        //if (!tmp_particle) continue;

        art::Ptr< simb::MCParticle > particle = particle_map.at(trk_id);

        // set meta data
        anab::BackTrackerMatchingData bt_data;
        bt_data.cleanliness = cleanliness;
        bt_data.completeness = completeness;

        // add association between reco. shower and simulated particle
        shower_mcparticle_assn->addSingle(shower_vector[shwr_idx], particle, bt_data);

        score_.push_back(score);
        cleanliness_.push_back(cleanliness);
        completeness_.push_back(completeness);
        pdg_.push_back(particle->PdgCode());
        number_hits_.push_back(hits.size());
      }
    }

    tree_->Fill();

    // put association into event
    //event.put(std::move(cluster_mcparticle_assn));
    event.put(std::move(track_mcparticle_assn));
    event.put(std::move(shower_mcparticle_assn));
  }

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::fill_metrics_vector_(
      art::ValidHandle< std::vector< simb::MCParticle > > const& particle_handle,
      art::ValidHandle< std::vector< sim::MCTrack > >     const& mctrack_handle,
      art::ValidHandle< std::vector< sim::MCShower > >    const& mcshower_handle,
      std::set< int >                                     const& g4_trk_ids_hits,
      std::vector< art::Ptr< recob::Hit > >               const& hits,
      std::vector< art::Ptr< recob::Hit > >               const& hit_vector,
      simutil::ParticleMap                                const& particle_map,
      simutil::ShowerParticleMaps                         const& shower_particle_maps,
      std::vector< metric_t >                                  & metrics)
  {
    // set of G4 track IDs that have been used
    std::set< int > g4_trk_id_set;

    // loop through MCTrack objects
    for (auto const& mctrack : *mctrack_handle)
    {
      if (!mctrack.size()) continue;

      // get G4 track ID
      int const g4_trk_id = mctrack.TrackID();

      // insert G4 track ID into set for BackTracker
      std::set< int > id;
      id.insert(g4_trk_id);
      g4_trk_id_set.insert(g4_trk_id);

      // continue if G4 track ID is not found in hits
      if (g4_trk_ids_hits.find(g4_trk_id) == g4_trk_ids_hits.end()) continue;

      // get charge-weighted cleanliness and charge-weighted
      // completeness using the BackTracker
      double const cleanliness =
          backtracker_->HitChargeCollectionPurity(id, hits);
      double const completeness =
          backtracker_->HitChargeCollectionEfficiency(
          id, hits, hit_vector, geo::k3D);

      // compute score
      double const score = cleanliness * completeness;

      // metric
      metric_t metric;
      metric.score = score;
      metric.cleanliness  = cleanliness;
      metric.completeness = completeness;
      metric.track_id = g4_trk_id;

      metrics.push_back(metric);
    }

    // loop through MCShower objects
    for (auto const& mcshower : *mcshower_handle)
    {
      // get G4 track ID
      int g4_trk_id = mcshower.TrackID();
      g4_trk_id_set.insert(g4_trk_id);

      // continue if track ID is not found in map
      if (shower_particle_maps.count(g4_trk_id) < 1) continue;

      // continue if no energy was deposited
      // FIXME: i don't think this works...
      //if (mcshower.DetProfile().E() < 1e-5) continue;

      // insert G4 track IDs into set for BackTracker
      std::set< int > id;

      // number of particles found in hits
      unsigned int number_particles = 0;

      for (auto const& particle_iter : shower_particle_maps.at(g4_trk_id))
      {
        int daughter_g4_trk_id = particle_iter.first;

        id.insert(daughter_g4_trk_id);
        g4_trk_id_set.insert(daughter_g4_trk_id);

        // increase counter if G4 track ID is found in hits
        if (g4_trk_ids_hits.find(daughter_g4_trk_id) != g4_trk_ids_hits.end())
            ++number_particles;
      }

      // continue if G4 track IDs are not found in hits
      if (number_particles < 1) continue;

      // get charge-weighted cleanliness and charge-weighted
      // completeness using the BackTracker
      double const cleanliness =
          backtracker_->HitChargeCollectionPurity(id, hits);
      double const completeness =
          backtracker_->HitChargeCollectionEfficiency(
          id, hits, hit_vector, geo::k3D);

      // compute score
      double const score = cleanliness * completeness;

      // metric
      metric_t metric;
      metric.score = score;
      metric.cleanliness  = cleanliness;
      metric.completeness = completeness;
      metric.track_id = g4_trk_id;

      metrics.push_back(metric);
    }

    // loop through MCParticles
    for (auto const& particle : *particle_handle)
    {
      // get G4 track ID of particle
      int const g4_trk_id = particle.TrackId();

      // continue if G4 track ID has already been used
      if (g4_trk_id_set.find(g4_trk_id) != g4_trk_id_set.end()) continue;

      // continue of G4 track ID is not found in hits
      if (g4_trk_ids_hits.find(g4_trk_id) == g4_trk_ids_hits.end()) continue;

      // insert G4 track ID into set for BackTracker
      std::set< int > id;
      id.insert(g4_trk_id);
      g4_trk_id_set.insert(g4_trk_id);

      // get charge-weighted cleanliness and charge-weighted
      // completeness using the BackTracker
      double const cleanliness =
          backtracker_->HitChargeCollectionPurity(id, hits);
      double const completeness =
          backtracker_->HitChargeCollectionEfficiency(
          id, hits, hit_vector, geo::k3D);

      // compute score
      double const score = cleanliness * completeness;

      // metric
      metric_t metric;
      metric.score = score;
      metric.cleanliness  = cleanliness;
      metric.completeness = completeness;
      metric.track_id = g4_trk_id;

      metrics.push_back(metric);
    }
  }

  //-----------------------------------------------------------------------
  void RecoMCParticleMatching::reset_()
  {
    score_.clear();
    cleanliness_.clear();
    completeness_.clear();
    pdg_.clear();
    number_hits_.clear();
  }

  DEFINE_ART_MODULE(RecoMCParticleMatching)

} // namespace npa

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // RecoMCParticleMatching_module
