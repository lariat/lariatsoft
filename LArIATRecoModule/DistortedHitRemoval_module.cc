//////////////////////////////////////////////////////////////
// Name:      DistortedHitRemoval
// Date:      6 September 2018
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef DistortedHitRemoval_Module
#define DistortedHitRemoval_Module

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
//#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
//#include "larcore/Geometry/Geometry.h"
//#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/ArtDataHelper/HitCreator.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"

// ROOT includes
//#include "TTree.h"

// C++ includes
#include <string>
#include <vector>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class DistortedHitRemoval : public art::EDProducer
{

 public:

  // standard constructor and destructor for an art module
  explicit DistortedHitRemoval(fhicl::ParameterSet const& pset);
  virtual ~DistortedHitRemoval();

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
  std::string hit_producer_label_;
  std::string track_producer_label_;
  double      z_low_;
  double      z_high_;

  // pointer to geometry provider
  //geo::GeometryCore const* geometry_;

  // pointer to detector properties
  //detinfo::DetectorProperties const* detector_;

  // pointers to TTree object
  //TTree * ttree_;

  // variables that will go into the TTree objects
  //int event_;     // number of the event being processed
  //int run_;       // number of the run being processed
  //int subrun_;    // number of the sub-run being processed

  // reset once per event
  //void reset_();
};

//-----------------------------------------------------------------------
// constructor
DistortedHitRemoval::DistortedHitRemoval(fhicl::ParameterSet const& pset)
{
  // reconfigure parameters
  this->reconfigure(pset);

  // get a pointer to the geometry service provider
  //geometry_ = &*(art::ServiceHandle<geo::Geometry>());

  // get a pointer to the detector properties provider
  //detector_ = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
  recob::HitCollectionCreator::declare_products(*this);
}

//-----------------------------------------------------------------------
// destructor
DistortedHitRemoval::~DistortedHitRemoval() {}

//-----------------------------------------------------------------------
void DistortedHitRemoval::beginJob()
{
  // access art's TFileService
  //art::ServiceHandle< art::TFileService > tfs;

  //ttree_ = tfs->make<TTree>("distortedhitremoval", "distortedhitremoval");
}

//-----------------------------------------------------------------------
void DistortedHitRemoval::beginRun(art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void DistortedHitRemoval::beginSubRun(art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void DistortedHitRemoval::endJob()
{}

//-----------------------------------------------------------------------
void DistortedHitRemoval::endRun(art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void DistortedHitRemoval::endSubRun(art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void DistortedHitRemoval::reconfigure(fhicl::ParameterSet const& pset)
{
  // read parameters from the .fcl file
  hit_producer_label_    = pset.get< std::string >("HitLabel");
  track_producer_label_  = pset.get< std::string >("TrackLabel");
  z_low_                 = pset.get< double >("ZLow", 4.0);
  z_high_                = pset.get< double >("ZHigh", 100.0);
}

//-----------------------------------------------------------------------
void DistortedHitRemoval::produce(art::Event & event)
{
  //-------------------------------------------------------------------
  // reset once per event
  //-------------------------------------------------------------------

  //this->reset_();

  //-------------------------------------------------------------------
  // get event, run, and subrun numbers
  //-------------------------------------------------------------------

  //event_  = event.id().event();
  //run_    = event.run();
  //subrun_ = event.subRun();

  //-------------------------------------------------------------------
  // get hits and spacepoints
  //-------------------------------------------------------------------

  // get all the hits in the event
  // art::ValidHandle< std::vector< recob::Hit > >
  auto const& hit_handle = event.getValidHandle< std::vector< recob::Hit > >
      (hit_producer_label_);

  // fill vector of hits
  std::vector< art::Ptr< recob::Hit > > hit_vector;
  art::fill_ptr_vector(hit_vector, hit_handle);

  // get all the reconstructed spacepoints in the event
  // art::ValidHandle< std::vector< recob::SpacePoint > >
  auto const& spacepoint_handle = event.getValidHandle< std::vector< recob::SpacePoint > >
      (track_producer_label_);

  // fill vector of spacepoints
  std::vector< art::Ptr< recob::SpacePoint > > spacepoint_vector;
  art::fill_ptr_vector(spacepoint_vector, spacepoint_handle);

  // find many hits from spacepoints
  const art::FindManyP< recob::Hit >
      find_many_hits_from_spacepoints(spacepoint_handle, event, track_producer_label_);

  //-------------------------------------------------------------------
  // find distorted hits
  //-------------------------------------------------------------------

  // initialize set for distorted hits
  std::set< recob::Hit > distorted_hits;

  // loop over spacepoints
  for (auto const& spacepoint : spacepoint_vector)
  {
    // get z-coordinate of spacepoint
    const double z = spacepoint->XYZ()[2];

    // skip if the spacepoint is within the fiducial region
    if (z > z_low_ and z < z_high_) continue;

    // std::vector< art::Ptr< recob::Hit > >
    auto const& hits = find_many_hits_from_spacepoints.at(spacepoint.key());

    // loop over associated hits
    for (auto const& hit : hits)
    {
      // store associated hit as a distorted hit
      distorted_hits.insert(*hit);
    } // end loop over associated hits

  } // end loop over spacepoints

  //-------------------------------------------------------------------
  // copy over undistorted hits
  //-------------------------------------------------------------------

  // find wire from hits
  const art::FindOneP< recob::Wire >
      wires(hit_handle, event, hit_producer_label_);

  // find raw digit from hits
  const art::FindOneP< raw::RawDigit >
      raw_digits(hit_handle, event, hit_producer_label_);

  // hit collection creator
  recob::HitCollectionCreator hit_collection_creator(
      *this, event, wires.isValid(), raw_digits.isValid());

  // loop over hits
  for (auto const& hit : hit_vector)
  {
    // skip if hit is distorted
    if (distorted_hits.find(*hit) != distorted_hits.end()) continue;

    // get associated wire and raw digit
    art::Ptr< recob::Wire >   const& wire      = wires.at(hit.key());
    art::Ptr< raw::RawDigit > const& raw_digit = raw_digits.at(hit.key());

    // copy undistorted hit, and associated wire and raw digit
    hit_collection_creator.emplace_back(*hit, wire, raw_digit);
  } // end loop over hits

  // put the hit collection and associations into the event
  hit_collection_creator.put_into(event);
}

//-----------------------------------------------------------------------
//void DistortedHitRemoval::reset_()
//{}

DEFINE_ART_MODULE(DistortedHitRemoval)

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // DistortedHitRemoval_module
