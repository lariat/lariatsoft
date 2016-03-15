////////////////////////////////////////////////////////////////////////
// Class:       ShowerFilter
// Module Type: filter
// File:        ShowerFilter_module.cc
//
// Generated at Wed Feb 25 17:12:56 2016 by Irene Nutini using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

//This is a very *rough and preliminary* version for a Shower Filter

//It was written for the PionAnalysis  (reco and analysis in lariatsoft_v04_34_00)
//since there was still a 15% of electron induced showers contamination in our sample (especially for negative pions) after all the Selection Cuts and TPC/WC track matching  

// The filter looks at how many tracks are reconstructed for the event
// Then it counts how many of them are "short tracks"; I define as "short", tracks whose length is less than 5 cm (parameter that can be changed in the fcl file).
// We look for short tracks since in general with LineCluster and Trackers (no shower algorithms) the shower events are reconstructed as many short tracks
// If in the event there are less than 2 short tracks (parameter that can be changed in the fcl file) this is considered as a "track-like" event,
// otherwise it's considered as a "shower-like" event
// "shower-like" events are rejected

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/FindOneP.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <TH2F.h>
#include "art/Framework/Services/Optional/TFileService.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

const int kMaxTrack      = 1000;  //maximum number of tracks

class ShowerFilter;

class ShowerFilter : public art::EDFilter {
public:
  explicit ShowerFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerFilter(ShowerFilter const &) = delete;
  ShowerFilter(ShowerFilter &&) = delete;
  ShowerFilter & operator = (ShowerFilter const &) = delete;
  ShowerFilter & operator = (ShowerFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
std::string fTrackModuleLabel;

int fnShortTks;
double fShortTkLength;


  
    int ntrks;
    double trkLength[kMaxTrack]={0.};
   
//Histograms
};


ShowerFilter::ShowerFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool ShowerFilter::filter(art::Event & e)
{
  // Implementation of required member function here.
   
  // #####################################
   // ### Getting the Track Information ###	
   // #####################################
   art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
   std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks	
   // === Filling the tracklist from the tracklistHandle ===	
   if (e.getByLabel(fTrackModuleLabel,trackListHandle))
      {art::fill_ptr_vector(tracklist, trackListHandle);}
      
  //std::cout << "There are: " << tracklist.size() << " tracks produced for this event" << std::endl;   
   
   int shortTrack=0;
  // int longTrack=0;
   
   ntrks = tracklist.size();
   for (size_t p = 0; p<tracklist.size(); ++p){
   
   trkLength[p]         = tracklist[p]->Length();
   //std::cout << "Track " << p << " lenght " << tracklist[p]->Length() << std::endl;
   
   if(trkLength[p] < fShortTkLength) shortTrack++;
   //if(trkLength[p] > 20.) longTrack++;
   
   }
   
   std::cout << "There are " << shortTrack << " short reco tracks" << std::endl;
   
   if(shortTrack < fnShortTks) {
     std::cout << "This is a track-like event " << std::endl;
     return true;
   }
   
   
   else{
     std::cout << "This can be a shower-like event" << std::endl;
     //std::cout << "there are " << longTrack << " longtracks" << std::endl;
    // std::cout << "There are: " << clusterlist.size() << " clusters (blurredClu) produced for this event" << std::endl;
   // ### Saving the number of clusters in this event ###
  
    return false;
   }
}

void ShowerFilter::beginJob()
{
  // Implementation of optional member function here.
   art::ServiceHandle<art::TFileService> tfs;
  
}

void ShowerFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here
  fTrackModuleLabel = p.get< std::string  >("TrackModuleLabel","pmtrack");
  fShortTkLength = p.get< double >     ("ShortTkLength", 5.0);
  fnShortTks   = p.get< int >        ("nShortTks"  ,    2);
  // fShowerModuleLabel = p.get< std::string  >("ShowerModuleLabel");
}

DEFINE_ART_MODULE(ShowerFilter)
