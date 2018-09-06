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
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"

// ROOT includes
#include "TTree.h"

// C++ includes
#include <cmath>
#include <memory>
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

  // reset once per event
  void reset_();
};

//-----------------------------------------------------------------------
// constructor
DistortedHitRemoval::DistortedHitRemoval(fhicl::ParameterSet const& pset)
{
  // reconfigure parameters
  this->reconfigure(pset);

  // get a pointer to the geometry service provider
  geometry_ = &*(art::ServiceHandle<geo::Geometry>());

  // get a pointer to the detector properties provider
  detector_ = lar::providerFrom<detinfo::DetectorPropertiesService>();

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
  art::ServiceHandle< art::TFileService > tfs;

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
  hit_producer_label_     = pset.get< std::string >("HitLabel");
  track_producer_label_   = pset.get< std::string >("TrackLabel");
}

//-----------------------------------------------------------------------
void DistortedHitRemoval::produce(art::Event & event)
{
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
  // get hits and tracks
  //-------------------------------------------------------------------

  // get all the hits in the event
  // art::ValidHandle< std::vector< recob::Hit > >
  auto const& hit_handle = event.getValidHandle< std::vector< recob::Hit > >
      (hit_producer_label_);

  // fill vector of hits
  std::vector< art::Ptr< recob::Hit > > hit_vector;
  art::fill_ptr_vector(hit_vector, hit_handle);

  // get all the reconstructed tracks in the event
  // art::ValidHandle< std::vector< recob::Track > >
  auto const& track_handle = event.getValidHandle< std::vector< recob::Track > >
      (track_producer_label_);

  // fill vector of tracks
  std::vector< art::Ptr< recob::Track > > track_vector;
  art::fill_ptr_vector(track_vector, track_handle);

  //-------------------------------------------------------------------
  // test
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

  // copy all hits
  for (auto const& hit : hit_vector)
  {
    art::Ptr< recob::Wire >   const& wire      = wires.at(hit.key());
    art::Ptr< raw::RawDigit > const& raw_digit = raw_digits.at(hit.key());

    hit_collection_creator.emplace_back(*hit, wire, raw_digit);
  }

  // put the hit collection and associations into the event
  hit_collection_creator.put_into(event);
}

//-----------------------------------------------------------------------
void DistortedHitRemoval::reset_()
{}

DEFINE_ART_MODULE(DistortedHitRemoval)

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // DistortedHitRemoval_module
