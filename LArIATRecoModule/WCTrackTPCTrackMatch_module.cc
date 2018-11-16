//////////////////////////////////////////////////////////////
// Name:      WCTrackTPCTrackMatch
// Date:      16 July 2018
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef WCTrackTPCTrackMatch_Module
#define WCTrackTPCTrackMatch_Module

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
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
//#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
//#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

// LArIATSoft includes
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/WCTrack.h"

// ROOT includes
#include "TTree.h"
#include "TVector.h"

// C++ includes
//#include <algorithm>
#include <cmath>
//#include <iterator>
//#include <limits>
//#include <map>
#include <memory>
#include <string>
#include <vector>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class WCTrackTPCTrackMatch : public art::EDProducer
{

 public:

  // standard constructor and destructor for an art module
  explicit WCTrackTPCTrackMatch(fhicl::ParameterSet const& pset);
  virtual ~WCTrackTPCTrackMatch();

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
  //std::string hit_producer_label_;
  //std::string cluster_producer_label_;
  std::string track_producer_label_;
  std::string tof_producer_label_;
  std::string wctrack_producer_label_;
  double      alpha_cut_;
  double      min_track_length_proj_z_;
  double      min_upstream_z_;
  double      max_upstream_z_;
  double      circular_cut_x_center_;
  double      circular_cut_y_center_;
  double      circular_cut_radius_;
  bool        isThisMC_;

  // vector for primary vtx
  TVector3 primary_vtx_;

  // pointer to geometry provider
  geo::GeometryCore const* geometry_;

  // pointer to detector properties
  detinfo::DetectorProperties const* detector_;

  // pointers to TTree object
  TTree * ttree_;

  // variables that will go into the TTree objects
  int event_;     // number of the event being processed
  int run_;       // number of the run being processed
  int subrun_;    // number of the sub-run being processed

  int number_tofs_;
  int number_wctracks_;

  std::vector< double > reco_tof_;
  std::vector< double > reco_wctrack_momentum_;

  std::vector< int > wctrack_missed_;
  std::vector< int > wctrack_picky_;

  std::vector< double > wctrack_x_;
  std::vector< double > wctrack_y_;
  std::vector< double > wctrack_theta_;
  std::vector< double > wctrack_phi_;

  std::vector< int >    track_key_;
  std::vector< double > delta_x_;
  std::vector< double > delta_y_;
  std::vector< double > delta_z_;
  std::vector< double > delta_r_;
  std::vector< double > alpha_;

  std::vector< double > preselected_delta_x_;
  std::vector< double > preselected_delta_y_;
  std::vector< double > preselected_delta_z_;
  std::vector< double > preselected_delta_r_;
  std::vector< double > preselected_alpha_;
  std::vector< double > preselected_upstream_x_;
  std::vector< double > preselected_upstream_y_;
  std::vector< double > preselected_upstream_z_;
  std::vector< double > preselected_downstream_x_;
  std::vector< double > preselected_downstream_y_;
  std::vector< double > preselected_downstream_z_;

  std::vector< double > selected_delta_x_;
  std::vector< double > selected_delta_y_;
  std::vector< double > selected_delta_z_;
  std::vector< double > selected_delta_r_;
  std::vector< double > selected_alpha_;
  std::vector< double > selected_upstream_x_;
  std::vector< double > selected_upstream_y_;
  std::vector< double > selected_upstream_z_;
  std::vector< double > selected_downstream_x_;
  std::vector< double > selected_downstream_y_;
  std::vector< double > selected_downstream_z_;

  std::vector< int > wctrack_proj_missed_;
  std::vector< int > wctrack_proj_picky_;

  std::vector< double > wctrack_proj_x_;
  std::vector< double > wctrack_proj_y_;
  std::vector< double > wctrack_proj_z_;
  std::vector< double > wctrack_proj_theta_;
  std::vector< double > wctrack_proj_phi_;

  std::vector< int >    track_proj_key_;
  std::vector< double > delta_proj_x_;
  std::vector< double > delta_proj_y_;
  std::vector< double > delta_proj_z_;
  std::vector< double > delta_proj_r_;
  std::vector< double > alpha_proj_;

  std::vector< double > preselected_delta_proj_x_;
  std::vector< double > preselected_delta_proj_y_;
  std::vector< double > preselected_delta_proj_z_;
  std::vector< double > preselected_delta_proj_r_;
  std::vector< double > preselected_alpha_proj_;
  std::vector< double > preselected_upstream_proj_x_;
  std::vector< double > preselected_upstream_proj_y_;
  std::vector< double > preselected_upstream_proj_z_;
  std::vector< double > preselected_downstream_proj_x_;
  std::vector< double > preselected_downstream_proj_y_;
  std::vector< double > preselected_downstream_proj_z_;

  std::vector< double > selected_delta_proj_x_;
  std::vector< double > selected_delta_proj_y_;
  std::vector< double > selected_delta_proj_z_;
  std::vector< double > selected_delta_proj_r_;
  std::vector< double > selected_alpha_proj_;
  std::vector< double > selected_upstream_proj_x_;
  std::vector< double > selected_upstream_proj_y_;
  std::vector< double > selected_upstream_proj_z_;
  std::vector< double > selected_downstream_proj_x_;
  std::vector< double > selected_downstream_proj_y_;
  std::vector< double > selected_downstream_proj_z_;

  TVector3 tpc_face_intersection(TVector3 const& a, TVector3 const& b);

  // reset once per event
  void reset_();
};

//-----------------------------------------------------------------------
// constructor
WCTrackTPCTrackMatch::WCTrackTPCTrackMatch(fhicl::ParameterSet const& pset)
{
  // reconfigure parameters
  this->reconfigure(pset);

  // get a pointer to the geometry service provider
  geometry_ = &*(art::ServiceHandle<geo::Geometry>());

  // get a pointer to the detector properties provider
  detector_ = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // produce recob::PFParticle objects
  produces< std::vector< recob::PFParticle > >();

  // produce associations
  produces< art::Assns< recob::PFParticle, recob::Track > >();
  produces< art::Assns< recob::Track, ldp::WCTrack > >();
}

//-----------------------------------------------------------------------
// destructor
WCTrackTPCTrackMatch::~WCTrackTPCTrackMatch() {}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::beginJob()
{
  // access art's TFileService
  art::ServiceHandle< art::TFileService > tfs;

  ttree_ = tfs->make<TTree>("wctracktpctrackmatch", "wctracktpctrackmatch");

  ttree_->Branch("event",  &event_,  "event/I");
  ttree_->Branch("run",    &run_,    "run/I");
  ttree_->Branch("subrun", &subrun_, "subrun/I");

  ttree_->Branch("number_tofs", &number_tofs_, "number_tofs/I");
  ttree_->Branch("number_wctracks", &number_wctracks_, "number_wctracks/I");

  ttree_->Branch("reco_tof", &reco_tof_);
  ttree_->Branch("reco_wctrack_momentum", &reco_wctrack_momentum_);

  ttree_->Branch("wctrack_missed", &wctrack_missed_);
  ttree_->Branch("wctrack_picky", &wctrack_picky_);

  ttree_->Branch("wctrack_x", &wctrack_x_);
  ttree_->Branch("wctrack_y", &wctrack_y_);
  ttree_->Branch("wctrack_theta", &wctrack_theta_);
  ttree_->Branch("wctrack_phi", &wctrack_phi_);

  ttree_->Branch("track_key", &track_key_);
  ttree_->Branch("delta_x", &delta_x_);
  ttree_->Branch("delta_y", &delta_y_);
  ttree_->Branch("delta_z", &delta_z_);
  ttree_->Branch("delta_r", &delta_r_);
  ttree_->Branch("alpha", &alpha_);

  ttree_->Branch("preselected_delta_x", &preselected_delta_x_);
  ttree_->Branch("preselected_delta_y", &preselected_delta_y_);
  ttree_->Branch("preselected_delta_z", &preselected_delta_z_);
  ttree_->Branch("preselected_delta_r", &preselected_delta_r_);
  ttree_->Branch("preselected_alpha",   &preselected_alpha_);
  ttree_->Branch("preselected_upstream_x", &preselected_upstream_x_);
  ttree_->Branch("preselected_upstream_y", &preselected_upstream_y_);
  ttree_->Branch("preselected_upstream_z", &preselected_upstream_z_);
  ttree_->Branch("preselected_downstream_x", &preselected_downstream_x_);
  ttree_->Branch("preselected_downstream_y", &preselected_downstream_y_);
  ttree_->Branch("preselected_downstream_z", &preselected_downstream_z_);

  ttree_->Branch("selected_delta_x", &selected_delta_x_);
  ttree_->Branch("selected_delta_y", &selected_delta_y_);
  ttree_->Branch("selected_delta_z", &selected_delta_z_);
  ttree_->Branch("selected_delta_r", &selected_delta_r_);
  ttree_->Branch("selected_alpha",   &selected_alpha_);
  ttree_->Branch("selected_upstream_x", &selected_upstream_x_);
  ttree_->Branch("selected_upstream_y", &selected_upstream_y_);
  ttree_->Branch("selected_upstream_z", &selected_upstream_z_);
  ttree_->Branch("selected_downstream_x", &selected_downstream_x_);
  ttree_->Branch("selected_downstream_y", &selected_downstream_y_);
  ttree_->Branch("selected_downstream_z", &selected_downstream_z_);

  ttree_->Branch("wctrack_proj_missed", &wctrack_proj_missed_);
  ttree_->Branch("wctrack_proj_picky", &wctrack_proj_picky_);

  ttree_->Branch("wctrack_proj_x", &wctrack_proj_x_);
  ttree_->Branch("wctrack_proj_y", &wctrack_proj_y_);
  ttree_->Branch("wctrack_proj_z", &wctrack_proj_z_);
  ttree_->Branch("wctrack_proj_theta", &wctrack_proj_theta_);
  ttree_->Branch("wctrack_proj_phi", &wctrack_proj_phi_);

  ttree_->Branch("track_proj_key", &track_proj_key_);
  ttree_->Branch("delta_proj_x", &delta_proj_x_);
  ttree_->Branch("delta_proj_y", &delta_proj_y_);
  ttree_->Branch("delta_proj_z", &delta_proj_z_);
  ttree_->Branch("delta_proj_r", &delta_proj_r_);
  ttree_->Branch("alpha_proj", &alpha_proj_);

  ttree_->Branch("preselected_delta_proj_x", &preselected_delta_proj_x_);
  ttree_->Branch("preselected_delta_proj_y", &preselected_delta_proj_y_);
  ttree_->Branch("preselected_delta_proj_z", &preselected_delta_proj_z_);
  ttree_->Branch("preselected_delta_proj_r", &preselected_delta_proj_r_);
  ttree_->Branch("preselected_alpha_proj",   &preselected_alpha_proj_);
  ttree_->Branch("preselected_upstream_proj_x", &preselected_upstream_proj_x_);
  ttree_->Branch("preselected_upstream_proj_y", &preselected_upstream_proj_y_);
  ttree_->Branch("preselected_upstream_proj_z", &preselected_upstream_proj_z_);
  ttree_->Branch("preselected_downstream_proj_x", &preselected_downstream_proj_x_);
  ttree_->Branch("preselected_downstream_proj_y", &preselected_downstream_proj_y_);
  ttree_->Branch("preselected_downstream_proj_z", &preselected_downstream_proj_z_);

  ttree_->Branch("selected_delta_proj_x", &selected_delta_proj_x_);
  ttree_->Branch("selected_delta_proj_y", &selected_delta_proj_y_);
  ttree_->Branch("selected_delta_proj_z", &selected_delta_proj_z_);
  ttree_->Branch("selected_delta_proj_r", &selected_delta_proj_r_);
  ttree_->Branch("selected_alpha_proj",   &selected_alpha_proj_);
  ttree_->Branch("selected_upstream_proj_x", &selected_upstream_proj_x_);
  ttree_->Branch("selected_upstream_proj_y", &selected_upstream_proj_y_);
  ttree_->Branch("selected_upstream_proj_z", &selected_upstream_proj_z_);
  ttree_->Branch("selected_downstream_proj_x", &selected_downstream_proj_x_);
  ttree_->Branch("selected_downstream_proj_y", &selected_downstream_proj_y_);
  ttree_->Branch("selected_downstream_proj_z", &selected_downstream_proj_z_);
}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::beginRun(art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::beginSubRun(art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::endJob()
{}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::endRun(art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::endSubRun(art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::reconfigure(fhicl::ParameterSet const& pset)
{
  // read parameters from the .fcl file
  //hit_producer_label_        = pset.get< std::string >("HitLabel");
  //cluster_producer_label_    = pset.get< std::string >("ClusterLabel");
  track_producer_label_      = pset.get< std::string >("TrackLabel");
  tof_producer_label_      = pset.get< std::string >("TOFLabel");
  wctrack_producer_label_  = pset.get< std::string >("WCTrackLabel");
  alpha_cut_               = pset.get< double >("AlphaCut",      0.3);
  min_track_length_proj_z_ = pset.get< double >("MinTrackLengthZProj", 4.0);
  min_upstream_z_          = pset.get< double >("MinUpstreamZ", 2.0);
  max_upstream_z_          = pset.get< double >("MaxUpstreamZ", 6.0);
  circular_cut_x_center_   = pset.get< double >("CircularCutXCenter", 1.6);
  circular_cut_y_center_   = pset.get< double >("CircularCutYCenter", -0.17);
  circular_cut_radius_     = pset.get< double >("CircularCutRadius", 3.5);
  isThisMC_                = pset.get< bool   >("IsThisMC", false);
}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::produce(art::Event & event)
{
  //-------------------------------------------------------------------
  // point to a collection of PFParticle objects
  //-------------------------------------------------------------------

  std::vector< recob::PFParticle > pfparticle_vector;

  //std::unique_ptr< std::vector< recob::PFParticle > >
  //    pfparticle_collection(new std::vector< recob::PFParticle >);

  std::unique_ptr< art::Assns< recob::PFParticle, recob::Track > >
      pfparticle_track_assn(new art::Assns< recob::PFParticle, recob::Track >);

  std::unique_ptr< art::Assns< recob::Track, ldp::WCTrack > >
      wc_track_assn(new art::Assns< recob::Track, ldp::WCTrack >);

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
  // get hits, clusters, tracks, and vertices
  //-------------------------------------------------------------------

  //// get all the hits in the event
  //// art::ValidHandle< std::vector< recob::Hit > >
  //auto hit_handle = event.getValidHandle< std::vector< recob::Hit > >
  //    (hit_producer_label_);

  //// get all the clusters in the event
  //// art::ValidHandle< std::vector< recob::Cluster > >
  //auto cluster_handle = event.getValidHandle< std::vector< recob::Cluster > >
  //    (cluster_producer_label_);

  // get all the reconstructed tracks in the event
  // art::ValidHandle< std::vector< recob::Track > >
  auto track_handle = event.getValidHandle< std::vector< recob::Track > >
      (track_producer_label_);

  //// fill vector of hits
  //std::vector< art::Ptr< recob::Hit > > hit_vector;
  //art::fill_ptr_vector(hit_vector, hit_handle);

  // fill vector of tracks
  std::vector< art::Ptr< recob::Track > > track_vector;
  art::fill_ptr_vector(track_vector, track_handle);

  // find many hits from tracks
  const art::FindManyP< recob::Hit >
      find_many_hits(track_handle, event, track_producer_label_);

  int selected_trk_idx = -1;
  int selected_proj_trk_idx = -1;
  ////int selected_wctrk_idx = -1;

  //---------------------------------------------------------------
  // get time of flight
  //---------------------------------------------------------------

  // get all the TOF objects in the event
  art::Handle< std::vector< ldp::TOF > > tof_handle;

  // fill vector of TOF objects
  std::vector< art::Ptr< ldp::TOF > > tof_vector;

  if (event.getByLabel(tof_producer_label_, tof_handle))
  {
    art::fill_ptr_vector(tof_vector, tof_handle);
  }

  for (auto const& tof : tof_vector)
  {
    //number_tofs_ += tof->NTOF();
    for (size_t tof_idx = 0; tof_idx < tof->NTOF(); ++tof_idx)
    {
      reco_tof_.push_back(tof->SingleTOF(tof_idx));
      number_tofs_ += 1;
    }
  }

  //---------------------------------------------------------------
  // get wire chamber tracks
  //---------------------------------------------------------------

  // get all the wire chamber tracks in the event
  // art::ValidHandle< std::vector< ldp::WCTrack > >
  //auto wctrack_handle = event.getValidHandle< std::vector< ldp::WCTrack > >
  //    (wctrack_producer_label_);
  art::Handle< std::vector< ldp::WCTrack > > wctrack_handle;

  // fill vector of WC tracks
  std::vector< art::Ptr< ldp::WCTrack > > wctrack_vector;
  //art::fill_ptr_vector(wctrack_vector, wctrack_handle);

  if (event.getByLabel(wctrack_producer_label_, wctrack_handle))
  {
    art::fill_ptr_vector(wctrack_vector, wctrack_handle);
    number_wctracks_ = wctrack_vector.size();
  }

  for (auto const& wctrack : wctrack_vector)
  {
    reco_wctrack_momentum_.push_back(wctrack->Momentum());
  }

  //std::cout << "Found " << wctrack_vector.size() << " WC tracks." << std::endl;

  // if no wire chamber track found, select TPC track closest to front face
  // of the TPC as long as there is only one TPC track that starts within a
  // few centimeters of the front face

  bool wctrack_flag = false;

  if (wctrack_vector.size() > 0) wctrack_flag = true;
  //if (wctrack_vector.size() == 1) wctrack_flag = true;

  // selected WC / TPC track parameters
  double selected_delta_r = 999.;
  double selected_delta_x = 999.;
  double selected_delta_y = 999.;
  double selected_delta_z = 999.;
  double selected_alpha   = 999.;
  double selected_upstream_x = 999.;
  double selected_upstream_y = 999.;
  double selected_upstream_z = 999.;
  double selected_downstream_x = 999.;
  double selected_downstream_y = 999.;
  double selected_downstream_z = 999.;

  double selected_delta_proj_r = 999.;
  double selected_delta_proj_x = 999.;
  double selected_delta_proj_y = 999.;
  double selected_delta_proj_z = 999.;
  double selected_alpha_proj   = 999.;
  double selected_upstream_proj_x = 999.;
  double selected_upstream_proj_y = 999.;
  double selected_upstream_proj_z = 999.;
  double selected_downstream_proj_x = 999.;
  double selected_downstream_proj_y = 999.;
  double selected_downstream_proj_z = 999.;

  // loop over WC tracks
  for (auto const& wctrack : wctrack_vector)
  {
    if (!wctrack_flag) break;

    double const wctrack_x     = wctrack->XYFace(0);  // cm
    double const wctrack_y     = wctrack->XYFace(1);  // cm
    double const wctrack_theta = wctrack->Theta();
    double const wctrack_phi   = wctrack->Phi();

    int const wctrack_missed = wctrack->WCMissed();
    int const wctrack_picky = static_cast< int > (wctrack->IsPicky());

    TVector3 wctrack_dir(0, 0, 1);
    wctrack_dir.SetTheta(wctrack_theta);
    wctrack_dir.SetPhi(wctrack_phi);

    wctrack_missed_.push_back(wctrack_missed);
    wctrack_picky_.push_back(wctrack_picky);

    wctrack_x_.push_back(wctrack_x);
    wctrack_y_.push_back(wctrack_y);
    wctrack_theta_.push_back(wctrack_theta);
    wctrack_phi_.push_back(wctrack_phi);

    // loop over TPC tracks
    for (size_t trk_idx = 0; trk_idx < track_vector.size(); ++trk_idx)
    {
      // get track
      art::Ptr< recob::Track > const& track = track_vector.at(trk_idx);

      // get track length
      //double const track_length = track->Length();

      // trajectory point indices
      int upstream_traj_pt_idx = -1;
      int downstream_traj_pt_idx = -1;

      // get first trajectory point
      recob::tracking::TrajectoryPoint_t const& first_traj_pt
          = track->Trajectory().TrajectoryPoint(track->Trajectory().FirstValidPoint());

      // get last trajectory point
      recob::tracking::TrajectoryPoint_t const& last_traj_pt
          = track->Trajectory().TrajectoryPoint(track->Trajectory().LastValidPoint());

      // upstream and downstream tracjectory point directions
      recob::tracking::Vector_t upstream_traj_pt_dir;
      recob::tracking::Vector_t downstream_traj_pt_dir;

      // get upstream and downstream trajectory point indices
      if (first_traj_pt.position.Z() < last_traj_pt.position.Z())
      {
        upstream_traj_pt_idx = track->Trajectory().FirstValidPoint();
        downstream_traj_pt_idx = track->Trajectory().LastValidPoint();
        upstream_traj_pt_dir = track->Trajectory().StartDirection();
        downstream_traj_pt_dir = track->Trajectory().EndDirection();
      }
      else
      {
        upstream_traj_pt_idx = track->Trajectory().LastValidPoint();
        downstream_traj_pt_idx = track->Trajectory().FirstValidPoint();
        upstream_traj_pt_dir = -1.0 * track->Trajectory().EndDirection();
        downstream_traj_pt_dir = -1.0 * track->Trajectory().StartDirection();
      }

      // get upstream trajectory point
      recob::tracking::TrajectoryPoint_t const& upstream_traj_pt
          = track->Trajectory().TrajectoryPoint(upstream_traj_pt_idx);

      // get downstream trajectory point
      recob::tracking::TrajectoryPoint_t const& downstream_traj_pt
          = track->Trajectory().TrajectoryPoint(downstream_traj_pt_idx);

      TVector3 upstream_dir(
          upstream_traj_pt_dir.X(), upstream_traj_pt_dir.Y(),
          upstream_traj_pt_dir.Z());

      double const track_length_proj_z
          = downstream_traj_pt.position.Z() - upstream_traj_pt.position.Z();

      double const delta_x = upstream_traj_pt.position.X() - wctrack_x;
      double const delta_y = upstream_traj_pt.position.Y() - wctrack_y;
      double const delta_z = upstream_traj_pt.position.Z();

      double const alpha = wctrack_dir.Angle(upstream_dir);

      double const delta_r = std::sqrt(
                                 (delta_x - circular_cut_x_center_) *
                                 (delta_x - circular_cut_x_center_) +
                                 (delta_y - circular_cut_y_center_) *
                                 (delta_y - circular_cut_y_center_));

      track_key_.push_back(trk_idx);
      delta_x_.push_back(delta_x);
      delta_y_.push_back(delta_y);
      delta_z_.push_back(delta_z);

      alpha_.push_back(alpha);

      delta_r_.push_back(delta_r);

      double const wctrack_proj_z = upstream_traj_pt.position.Z();
      TVector3 const wctrack_proj = wctrack->ProjectionAtZ(wctrack_proj_z, isThisMC_);
      double const wctrack_proj_x = wctrack_proj.X();
      double const wctrack_proj_y = wctrack_proj.Y();

      double const wctrack_proj_theta = wctrack->DownstreamDir().Theta();
      double const wctrack_proj_phi = wctrack->DownstreamDir().Phi();

      double const delta_proj_x = upstream_traj_pt.position.X() - wctrack_proj_x;
      double const delta_proj_y = upstream_traj_pt.position.Y() - wctrack_proj_y;
      double const delta_proj_z = upstream_traj_pt.position.Z() - wctrack_proj_z;

      double const alpha_proj = wctrack->DownstreamDir().Angle(upstream_dir);

      double const delta_proj_r = std::sqrt(
                                      (delta_proj_x - circular_cut_x_center_) *
                                      (delta_proj_x - circular_cut_x_center_) +
                                      (delta_proj_y - circular_cut_y_center_) *
                                      (delta_proj_y - circular_cut_y_center_));

      wctrack_proj_missed_.push_back(wctrack_missed);
      wctrack_proj_picky_.push_back(wctrack_picky);

      wctrack_proj_x_.push_back(wctrack_proj_x);
      wctrack_proj_y_.push_back(wctrack_proj_y);
      wctrack_proj_z_.push_back(wctrack_proj_z);
      wctrack_proj_theta_.push_back(wctrack_theta);
      wctrack_proj_phi_.push_back(wctrack_phi);

      track_proj_key_.push_back(trk_idx);
      delta_proj_x_.push_back(delta_proj_x);
      delta_proj_y_.push_back(delta_proj_y);
      delta_proj_z_.push_back(delta_proj_z);
      delta_proj_r_.push_back(delta_proj_r);
      alpha_proj_.push_back(alpha_proj);

      if (delta_r < circular_cut_radius_ and
          alpha < alpha_cut_)
      {
        preselected_delta_r_.push_back(delta_r);
        preselected_delta_x_.push_back(delta_x);
        preselected_delta_y_.push_back(delta_y);
        preselected_delta_z_.push_back(delta_z);
        preselected_alpha_.push_back(alpha);
        preselected_upstream_x_.push_back(upstream_traj_pt.position.X());
        preselected_upstream_y_.push_back(upstream_traj_pt.position.Y());
        preselected_upstream_z_.push_back(upstream_traj_pt.position.Z());
        preselected_downstream_x_.push_back(downstream_traj_pt.position.X());
        preselected_downstream_y_.push_back(downstream_traj_pt.position.Y());
        preselected_downstream_z_.push_back(downstream_traj_pt.position.Z());

        if (delta_r < selected_delta_r and
            track_length_proj_z > min_track_length_proj_z_ and
            upstream_traj_pt.position.Z() > min_upstream_z_ and
            upstream_traj_pt.position.Z() < max_upstream_z_)
        {
          selected_delta_r = delta_r;
          selected_delta_x = delta_x;
          selected_delta_y = delta_y;
          selected_delta_z = delta_z;
          selected_alpha   = alpha;
          selected_trk_idx = trk_idx;
          selected_upstream_x = upstream_traj_pt.position.X();
          selected_upstream_y = upstream_traj_pt.position.Y();
          selected_upstream_z = upstream_traj_pt.position.Z();
          selected_downstream_x = downstream_traj_pt.position.X();
          selected_downstream_y = downstream_traj_pt.position.Y();
          selected_downstream_z = downstream_traj_pt.position.Z();
        }
      }

      if (delta_proj_r < circular_cut_radius_ and
          alpha_proj < alpha_cut_)
      {
        preselected_delta_proj_r_.push_back(delta_proj_r);
        preselected_delta_proj_x_.push_back(delta_proj_x);
        preselected_delta_proj_y_.push_back(delta_proj_y);
        preselected_delta_proj_z_.push_back(delta_proj_z);
        preselected_alpha_proj_.push_back(alpha_proj);
        preselected_upstream_proj_x_.push_back(upstream_traj_pt.position.X());
        preselected_upstream_proj_y_.push_back(upstream_traj_pt.position.Y());
        preselected_upstream_proj_z_.push_back(upstream_traj_pt.position.Z());
        preselected_downstream_proj_x_.push_back(downstream_traj_pt.position.X());
        preselected_downstream_proj_y_.push_back(downstream_traj_pt.position.Y());
        preselected_downstream_proj_z_.push_back(downstream_traj_pt.position.Z());

        if (delta_proj_r < selected_delta_proj_r and
            track_length_proj_z > min_track_length_proj_z_ and
            upstream_traj_pt.position.Z() > min_upstream_z_ and
            upstream_traj_pt.position.Z() < max_upstream_z_)
        {
          selected_delta_proj_r = delta_proj_r;
          selected_delta_proj_x = delta_proj_x;
          selected_delta_proj_y = delta_proj_y;
          selected_delta_proj_z = delta_proj_z;
          selected_alpha_proj   = alpha_proj;
          selected_proj_trk_idx = trk_idx;
          selected_upstream_proj_x = upstream_traj_pt.position.X();
          selected_upstream_proj_y = upstream_traj_pt.position.Y();
          selected_upstream_proj_z = upstream_traj_pt.position.Z();
          selected_downstream_proj_x = downstream_traj_pt.position.X();
          selected_downstream_proj_y = downstream_traj_pt.position.Y();
          selected_downstream_proj_z = downstream_traj_pt.position.Z();
        }
      }

    } // end loop over TPC tracks

    if (selected_trk_idx > -1)
    {
      selected_delta_r_.push_back(selected_delta_r);

      selected_delta_x_.push_back(selected_delta_x);
      selected_delta_y_.push_back(selected_delta_y);
      selected_delta_z_.push_back(selected_delta_z);

      selected_alpha_.push_back(selected_alpha);

      selected_upstream_x_.push_back(selected_upstream_x);
      selected_upstream_y_.push_back(selected_upstream_y);
      selected_upstream_z_.push_back(selected_upstream_z);

      selected_downstream_x_.push_back(selected_downstream_x);
      selected_downstream_y_.push_back(selected_downstream_y);
      selected_downstream_z_.push_back(selected_downstream_z);
    }

    if (selected_proj_trk_idx > -1)
    {
      selected_delta_proj_r_.push_back(selected_delta_proj_r);

      selected_delta_proj_x_.push_back(selected_delta_proj_x);
      selected_delta_proj_y_.push_back(selected_delta_proj_y);
      selected_delta_proj_z_.push_back(selected_delta_proj_z);

      selected_alpha_proj_.push_back(selected_alpha_proj);

      selected_upstream_proj_x_.push_back(selected_upstream_proj_x);
      selected_upstream_proj_y_.push_back(selected_upstream_proj_y);
      selected_upstream_proj_z_.push_back(selected_upstream_proj_z);

      selected_downstream_proj_x_.push_back(selected_downstream_proj_x);
      selected_downstream_proj_y_.push_back(selected_downstream_proj_y);
      selected_downstream_proj_z_.push_back(selected_downstream_proj_z);
    }

  } // end loop over WC tracks

  // fill TTree object
  ttree_->Fill();

  //if (selected_trk_idx > -1)
  //{
  //  std::vector< size_t > pfparticle_daughter_indices;
  //  pfparticle_vector.emplace_back(211, 0, 0, pfparticle_daughter_indices);

  //  art::PtrMaker< recob::PFParticle > make_pfparticle_ptr(event, *this);

  //  art::Ptr< recob::PFParticle > pfparticle_ptr
  //      = make_pfparticle_ptr(pfparticle_vector.size() - 1);

  //  pfparticle_track_assn->addSingle(
  //      pfparticle_ptr, track_vector.at(selected_trk_idx));

  //  wc_track_assn->addSingle(
  //      track_vector.at(selected_trk_idx), wctrack_vector.front());
  //}

  if (selected_proj_trk_idx > -1)
  {
    std::vector< size_t > pfparticle_daughter_indices;
    pfparticle_vector.emplace_back(211, 0, 0, pfparticle_daughter_indices);

    art::PtrMaker< recob::PFParticle > make_pfparticle_ptr(event, *this);

    art::Ptr< recob::PFParticle > pfparticle_ptr
        = make_pfparticle_ptr(pfparticle_vector.size() - 1);

    pfparticle_track_assn->addSingle(
        pfparticle_ptr, track_vector.at(selected_proj_trk_idx));

    wc_track_assn->addSingle(
        track_vector.at(selected_proj_trk_idx), wctrack_vector.front());
  }

  // convert vector to unique_ptr
  std::unique_ptr< std::vector< recob::PFParticle > >
      pfparticle_collection(
          new std::vector< recob::PFParticle >(std::move(pfparticle_vector)));

  // put collections into event
  event.put(std::move(pfparticle_collection));
  event.put(std::move(pfparticle_track_assn));
  event.put(std::move(wc_track_assn));
}

//-----------------------------------------------------------------------
TVector3 WCTrackTPCTrackMatch::tpc_face_intersection(
    TVector3 const& a, TVector3 const& b)
{
  if (a.X() == b.X()
      and a.Y() == b.Y()
      and a.Z() == b.Z())
  {
    return a;
  }

  else if (a.X() == b.X()
      or a.Y() == b.Y()
      or a.Z() == b.Z())
  {
    return TVector3();
  }

  double x = - (a.Z() / (b.Z() - a.Z())) * (b.X() - a.X()) + a.X();
  double y = - (a.Z() / (b.Z() - a.Z())) * (b.Y() - a.Y()) + a.Y();
  double z = 0;

  TVector3 vector(x, y, z);

  return vector;
}

//-----------------------------------------------------------------------
void WCTrackTPCTrackMatch::reset_()
{
  number_tofs_ = 0;
  number_wctracks_= 0;

  reco_tof_.clear();
  reco_wctrack_momentum_.clear();

  wctrack_missed_.clear();
  wctrack_picky_.clear();

  wctrack_x_.clear();
  wctrack_y_.clear();
  wctrack_theta_.clear();
  wctrack_phi_.clear();

  track_key_.clear();
  delta_x_.clear();
  delta_y_.clear();
  delta_z_.clear();
  delta_r_.clear();
  alpha_.clear();

  preselected_delta_x_.clear();
  preselected_delta_y_.clear();
  preselected_delta_z_.clear();
  preselected_delta_r_.clear();
  preselected_alpha_.clear();
  preselected_upstream_x_.clear();
  preselected_upstream_y_.clear();
  preselected_upstream_z_.clear();
  preselected_downstream_x_.clear();
  preselected_downstream_y_.clear();
  preselected_downstream_z_.clear();

  selected_delta_x_.clear();
  selected_delta_y_.clear();
  selected_delta_z_.clear();
  selected_delta_r_.clear();
  selected_alpha_.clear();
  selected_upstream_x_.clear();
  selected_upstream_y_.clear();
  selected_upstream_z_.clear();
  selected_downstream_x_.clear();
  selected_downstream_y_.clear();
  selected_downstream_z_.clear();

  wctrack_proj_missed_.clear();
  wctrack_proj_picky_.clear();

  wctrack_proj_x_.clear();
  wctrack_proj_y_.clear();
  wctrack_proj_z_.clear();
  wctrack_proj_theta_.clear();
  wctrack_proj_phi_.clear();

  track_proj_key_.clear();
  delta_proj_x_.clear();
  delta_proj_y_.clear();
  delta_proj_z_.clear();
  delta_proj_r_.clear();
  alpha_proj_.clear();

  preselected_delta_proj_x_.clear();
  preselected_delta_proj_y_.clear();
  preselected_delta_proj_z_.clear();
  preselected_delta_proj_r_.clear();
  preselected_alpha_proj_.clear();
  preselected_upstream_proj_x_.clear();
  preselected_upstream_proj_y_.clear();
  preselected_upstream_proj_z_.clear();
  preselected_downstream_proj_x_.clear();
  preselected_downstream_proj_y_.clear();
  preselected_downstream_proj_z_.clear();

  selected_delta_proj_x_.clear();
  selected_delta_proj_y_.clear();
  selected_delta_proj_z_.clear();
  selected_delta_proj_r_.clear();
  selected_alpha_proj_.clear();
  selected_upstream_proj_x_.clear();
  selected_upstream_proj_y_.clear();
  selected_upstream_proj_z_.clear();
  selected_downstream_proj_x_.clear();
  selected_downstream_proj_y_.clear();
  selected_downstream_proj_z_.clear();
}

DEFINE_ART_MODULE(WCTrackTPCTrackMatch)

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // WCTrackTPCTrackMatch_module
