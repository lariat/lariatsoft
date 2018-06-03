/////////////////////////////////////////
//
//  MuonRangeStackHits Class
//
//  gkpullia@syr.edu
// COMMENT
///////////////////////////////////////////

#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "cetlib_except/exception.h"
#include <math.h>
#include <vector>
#include <map>
#include <iostream>

namespace ldp{

  //#######################################
  MuonRangeStackHits::MuonRangeStackHits()
  {
    std::map<int,std::vector<int> > paddlemap;
    std::vector<int> blank1;
    blank1.push_back(0);
    for(int i=0; i<16;++i){
      
      paddlemap.emplace(i,blank1);
    }

  }
  
  //##########################################
  MuonRangeStackHits::MuonRangeStackHits(std::map<int, std::vector<int> > const& paddlemap,
                                         std::vector<ldp::MuRSTrack>           const& trackVector )
  {
    if( paddlemap.size() == 0 && trackVector.size() == 0 ) fIsInitializedEmpty = true;
    else{ fIsInitializedEmpty = false; }
    
    //Setting the data members with the reconstructed input
    fPaddleTimeTickMap=paddlemap;
    fMuRSTrackVector=trackVector;
    
    //Looping through the tracks and setting their penetration depths
    for( size_t iTrack = 0; iTrack < fMuRSTrackVector.size(); ++iTrack ){
      int PenPlane = fMuRSTrackVector.at(iTrack).HitVect.size();
      fMuRSTrackVector.at(iTrack).penetrationDepth = PenPlane;
      //Arbitrarily Setting to first hit (since all hits are really close together in time  anyway)
      fMuRSTrackVector.at(iTrack).arrivalTime = fMuRSTrackVector.at(iTrack).HitVect.at(0).at(2);
    }
  }
  
  

  //########################################
  std::vector<int> MuonRangeStackHits::TimeTick(int iPaddle) const
  {
    int lastpaddle=fPaddleTimeTickMap.end()->first;
    if(iPaddle > lastpaddle  ){
      throw cet::exception("MuonRangeStackHits")
      << "Requested time tick vector for paddle "
      << iPaddle
      << ".  That doesn't exist for this event.  "
      << "The last paddle number you can reference is "
      << lastpaddle;
    }
    return fPaddleTimeTickMap.find(iPaddle)->second;
  }

  //########################################
  ldp::MuRSTrack MuonRangeStackHits::GetTrack( int iTrack )
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits")
      << "Requested track number: "
      << iTrack
      << ". This doesn't exist for this event. "
      << "The last track index you can access is: "
      << fMuRSTrackVector.size() - 1;
    }
    return fMuRSTrackVector.at(iTrack);
  }

  //########################################  
  int MuonRangeStackHits::GetArrivalTime( int iTrack ) const
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits")
      << "Requested arrival time for track number: "
      << iTrack
      << ". This doesn't exist for this event. "
      << "The last track index you can access is: "
      << fMuRSTrackVector.size()-1;
    }
    return fMuRSTrackVector.at(iTrack).arrivalTime;
  }

  //########################################  
  int MuonRangeStackHits::GetPenetrationDepth( int iTrack ) const
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits")
      << "Requested penetration depth for track number: "
      << iTrack
      << ". This doesn't exist for this event. "
      << "The last track index you can access is: "
      << fMuRSTrackVector.size()-1;
    }
    return fMuRSTrackVector.at(iTrack).penetrationDepth;
  }

}// end namespace
