////////////////////////////////////////////////////////////////////////
// Class:       UpStreamTPCTrackFilter
// Module Type: filter
// File:        UpStreamTPCTrackFilter_module.cc
//
// Generated at Tue Jan  5 17:52:55 2016 by Jonathan Asaadi using artmod
// from cetpkgsupport v1_08_06.
//
// This module filters based on the number of TPC tracks (by default requiring
// only one) in the upstream portion of the TPC (by default it is requiring a
// spacepoint in the first 2cm in Z) and that have their first point
// in the active volume (0 < x < 42.5    -20 < y < 20   0 < z < 90  (cm)
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

#include <TH1F.h>

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

// #####################
// ### ROOT includes ###
// #####################
#include <TH1F.h>


class UpStreamTPCTrackFilter;

class UpStreamTPCTrackFilter : public art::EDFilter {
public:
  explicit UpStreamTPCTrackFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UpStreamTPCTrackFilter(UpStreamTPCTrackFilter const &) = delete;
  UpStreamTPCTrackFilter(UpStreamTPCTrackFilter &&) = delete;
  UpStreamTPCTrackFilter & operator = (UpStreamTPCTrackFilter const &) = delete;
  UpStreamTPCTrackFilter & operator = (UpStreamTPCTrackFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  // Declare member data here.
  std::string fTrackModuleLabel;
  double fupstreamZPosition;
  double fnTracksUpstream;

  TH1F*  fNumberTracksTot;                  
  TH1F*  fNumberTracksInUSPortion;
  TH1F*  fNumberTracksInUSPortionPassingCut;     
  
  
};

// ---------------------- Parameter Setting ---------------------
UpStreamTPCTrackFilter::UpStreamTPCTrackFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

// -------------------- FHICL Parameter Set ---------------------
void UpStreamTPCTrackFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel		= p.get< std::string >("TrackModuleLabel");
  fupstreamZPosition            = p.get< double >("upstreamZPosition", 2.0);
  fnTracksUpstream              = p.get< int >("nTracksUpstream", 1);
  
}

// ---------------------- Begin Job ---------------------------
void UpStreamTPCTrackFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fNumberTracksTot                  = tfs->make<TH1F>("NumberTracksTot"  ,"N of Tracks Totale;N Tracks;N evt"             ,30,-0.5,29.5);
  fNumberTracksInUSPortion          = tfs->make<TH1F>("NumberTracksInUSPortion"  ,"N of tracks in US portion;N Tracks;N evt"             ,30,-0.5,29.5);
  fNumberTracksInUSPortionPassingCut= tfs->make<TH1F>("NumberTracksInUSPortionPassingCut","N of tracks in US portion passing Cuts;N Tracks;N evt",30,-0.5,29.5);

  // Implementation of optional member function here.
}


// ---------------------- Event Loop ---------------------------
bool UpStreamTPCTrackFilter::filter(art::Event & evt)
{

// #####################################
// ### Getting the Track Information ###
// #####################################
art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
   
// === Filling the tracklist from the tracklistHandle ===
 if (!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return false;
 
 art::fill_ptr_vector(tracklist, trackListHandle);
 
 // === Association between SpacePoints and Tracks ===
art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);


// ### Boolian for events w/ track which ###
// ###     starts at the front face      ###
bool TrackSptsZCut = false;

// ### Counting the number of tracks in the upstream ###
// ###               for this event                  ###
int ntrksUpStream = 0;

// ### Looping over tracks ###
for(size_t i=0; i<tracklist.size();++i)
   {
   // ### Assume this track won't pass ###
   TrackSptsZCut = false;
   
   // ### Setting a temp variable for this track ###
   float tempZpoint = 100;
   
   // ### Grabbing the SpacePoints associated with this track ###
   std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
   
   // ########################################
   // ### Looping over all the SpacePoints ###
   // ########################################
   for (size_t j = 0; j<spts.size(); ++j)
      {
      // ################################################################################ 
      // ### Tracking the lowest Z point that is inside fiducial boundries of the TPC ###
      // ################################################################################
      if(spts[j]->XYZ()[2] < tempZpoint && spts[j]->XYZ()[2] > 0 &&
         spts[j]->XYZ()[2] < 90 && spts[j]->XYZ()[0] > 0 && spts[j]->XYZ()[0] < 42.5 &&
	 spts[j]->XYZ()[1] > -20 && spts[j]->XYZ()[1] < 20 )
         
	 {tempZpoint = spts[j]->XYZ()[2];}
	 
      
      // ### Only passing events with a track that has ###
      // ###  a spacepoint within the first N cm in Z  ### 
      // ###    And requiring it to be inside the TPC  ###
      if(tempZpoint < fupstreamZPosition)
      	{TrackSptsZCut = true;}
      
      }//<---End j loop
   
   // ### If this track passes then bump our counter ###
   if(TrackSptsZCut){ntrksUpStream++;}
   
   }//<---End i loop

 fNumberTracksTot->Fill(tracklist.size());
 fNumberTracksInUSPortion->Fill(ntrksUpStream);
// ### If we didn't find enough upstream tracks then return false ###
if(ntrksUpStream < fnTracksUpstream){return false;}
// ### Otherwise, keep the event ###
else{
  fNumberTracksInUSPortionPassingCut->Fill(ntrksUpStream); 
  return true;
 }



}//<---End evt loop


// ---------------------- End Job ---------------------------
void UpStreamTPCTrackFilter::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(UpStreamTPCTrackFilter)
