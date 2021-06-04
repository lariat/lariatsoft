////////////////////////////////////////////////////////////////////////
// Class:       UpStreamTPCMultiplicityFilter
// Module Type: filter
// File:        UpStreamTPCMultiplicityFilter_module.cc
//
// Generated at Wed Jan  6 09:09:31 2016 by Elena Gramellini using artmod
// from cetpkgsupport v1_08_06.
//
// This module filters events based the number of tracks present in the
// the upstream portion of the TPC.
// An event is rejected if it has more then N tracks in the first XX cm.
// N and XX are fcl parameters.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework Includes ###
// ##########################
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


class UpStreamTPCMultiplicityFilter;

class UpStreamTPCMultiplicityFilter : public art::EDFilter {
public:
  explicit UpStreamTPCMultiplicityFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UpStreamTPCMultiplicityFilter(UpStreamTPCMultiplicityFilter const &) = delete;
  UpStreamTPCMultiplicityFilter(UpStreamTPCMultiplicityFilter &&) = delete;
  UpStreamTPCMultiplicityFilter & operator = (UpStreamTPCMultiplicityFilter const &) = delete;
  UpStreamTPCMultiplicityFilter & operator = (UpStreamTPCMultiplicityFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;
  // Optional functions
  void reconfigure(fhicl::ParameterSet const & p) ;
  void beginJob   (                ) override;

private:

  // Declare member data here.
  std::string fTrackModuleLabel;
  double fupstreamZPosition;
  double fnTracksUpstream;

  TH1F*  fNumberTracksTot;                  
  TH1F*  fNumberTracksInUSPortion;
  TH1F*  fNumberTracksInUSPortionVolume;          
  TH1F*  fNumberTracksInUSPortionPassingCut;
  TH1F*  fZUSPositionTrack;
};


UpStreamTPCMultiplicityFilter::UpStreamTPCMultiplicityFilter(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

void UpStreamTPCMultiplicityFilter::beginJob()
{
  // Declare checks histograms
  art::ServiceHandle<art::TFileService> tfs;

  fNumberTracksTot                  = tfs->make<TH1F>("NumberTracksTot"                  ,"Tot N of tracks;N Tracks;N evt"                       ,30,-0.5,29.5);
  fNumberTracksInUSPortion          = tfs->make<TH1F>("NumberTracksInUSPortion"          ,"N of tracks in US portion before cut;N Tracks;N evt"             ,30,-0.5,29.5);
  fNumberTracksInUSPortionVolume    = tfs->make<TH1F>("NumberTracksInUSPortionInVolume"  ,"N of tracks in US portion in fiducial Volume;N Tracks;N evt"             ,30,-0.5,29.5);
  fNumberTracksInUSPortionPassingCut= tfs->make<TH1F>("NumberTracksInUSPortionPassingCut","N of tracks in US portion passing Cuts;N Tracks;N evt",30,-0.5,29.5);
  fZUSPositionTrack                 = tfs->make<TH1F>("ZUSPositionTrack"                 ,"Upstream Z point of the track; Upstream Z; N tracks"  ,270,0.,90.);
}

bool UpStreamTPCMultiplicityFilter::filter(art::Event & evt)
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
  bool TrackSptsVol  = false;
  // ### Counting the number of tracks in the upstream ###
  // ###               for this event                  ###
  int ntrksUpStream = 0;
  int ntrksVol      = 0;
  // ### Looping over tracks ###
  for(size_t i=0; i<tracklist.size();++i)
    {
      // ### Assume this track won't pass ###
      TrackSptsZCut = false;
      TrackSptsVol  = false;
      // ### Setting a temp variable for this track ###
      float tempZpoint = 100.;
      
      // ### Grabbing the SpacePoints associated with this track ###
      std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
     
      if (spts.size()>1) 
	{
	  if (spts[0]->XYZ()[2] > spts[1]->XYZ()[2] ) {
	    fZUSPositionTrack->Fill(spts[0]->XYZ()[2]);
	  }else fZUSPositionTrack->Fill(spts[spts.size() - 1]->XYZ()[2]);
	}
      // ########################################
      // ### Looping over all the SpacePoints ###
      // ########################################
      for (size_t j = 0; j<spts.size(); ++j)
      {
	// ################################################################################# 
	// ### Recording the lowest Z point that is inside fiducial boundries of the TPC ###
	// #################################################################################
	if(spts[j]->XYZ()[2] < tempZpoint && 
	   spts[j]->XYZ()[2] > 0   && spts[j]->XYZ()[2] < 90   && 
	   spts[j]->XYZ()[0] > 0   && spts[j]->XYZ()[0] < 42.5 &&
	   spts[j]->XYZ()[1] > -20 && spts[j]->XYZ()[1] < 20 )
	  {
	    tempZpoint = spts[j]->XYZ()[2];
	    TrackSptsVol = true;
	  }
	
	
	// ### Only passing events with a track that has ###
	// ###  a spacepoint within the first X cm in Z  ### 
	// ###    And requiring it to be inside the TPC  ###
	if(tempZpoint < fupstreamZPosition)
	  {
	    TrackSptsZCut = true; 
	    break; // Break from spacepoint loop:
	           //If you're lucky and the fist spt satisfies the condition, you don't need to look further
	  }
	
      }//<---End of loop on spacepoints
      
     
      // ### If this track passes then bump our counter ###
      if(TrackSptsZCut){ntrksUpStream++;}
      if(TrackSptsVol){ntrksVol++;}
      
    }//<---End of the loop on tracks
  
  //  Fill check histograms
  fNumberTracksTot->Fill(tracklist.size());
  fNumberTracksInUSPortion->Fill(ntrksUpStream);
  fNumberTracksInUSPortionVolume->Fill(ntrksVol);
  // ### If we found too many upstream tracks then return false ###
  if(ntrksUpStream > fnTracksUpstream) return false;
  // ### Otherwise, keep the event and save the number of tracks###
  fNumberTracksInUSPortionPassingCut->Fill(ntrksUpStream);
  return true;
}

// -------------------- FHICL Parameter Set ---------------------
void UpStreamTPCMultiplicityFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fTrackModuleLabel  = p.get< std::string >("TrackModuleLabel");
  fupstreamZPosition = p.get< double >     ("MaxUpstreamZPosition", 14.0);
  fnTracksUpstream   = p.get< int >        ("MaxNTracksUpstream"  ,    4);
  
}

DEFINE_ART_MODULE(UpStreamTPCMultiplicityFilter)
