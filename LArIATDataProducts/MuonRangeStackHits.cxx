/////////////////////////////////////////
//
//  MuonRangeStackHits Class
//
//  gkpullia@syr.edu
// COMMENT
///////////////////////////////////////////

#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "cetlib/exception.h"
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
  MuonRangeStackHits::MuonRangeStackHits(std::map<int, std::vector<int> > paddlemap,
					 std::vector<MuRSTrack> trackVector )
  {
    //Setting the data members with the reconstructed input
    fPaddleTimeTickMap=paddlemap;
    fMuRSTrackVector=trackVector;
    
    //Looping through the tracks and setting their penetration depths
    for( size_t iTrack = 0; iTrack < fMuRSTrackVector.size(); ++iTrack ){
      int highestPlane = -2;
      for( size_t iHit = 0; iHit < fMuRSTrackVector.at(iTrack).HitVect.size(); ++iHit ){
	if( fMuRSTrackVector.at(iTrack).HitVect.at(iHit).at(0) > highestPlane ){
	  highestPlane = fMuRSTrackVector.at(iTrack).HitVect.at(iHit).at(0);
	}
      }
      fMuRSTrackVector.at(iTrack).penetrationDepth = highestPlane;
      fMuRSTrackVector.at(iTrack).arrivalTime = fMuRSTrackVector.at(iTrack).HitVect.at(0).at(2); //Arbitrarily Setting to first hit (since all hits are really close together in time  anyway)
    }
  }
  
  

  //########################################
  std::vector<int> MuonRangeStackHits::TimeTick(int iPaddle) const
  {
    int lastpaddle=fPaddleTimeTickMap.end()->first;
    if(iPaddle > lastpaddle  ){
      throw cet::exception("MuonRangeStackHits") << "Requested time tick vector for paddle "<<iPaddle<<".  That doesn't exist for this event.  The last paddle number you can reference is ""<<lastpaddle<<" <<"\n";
    }
    return fPaddleTimeTickMap.find(iPaddle)->second;
  }

  //########################################
  MuRSTrack MuonRangeStackHits::GetTrack( int iTrack )
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits") << "Requested track number: " << iTrack << ". This doesn't exist for this event. The last track index you can access is: " << fMuRSTrackVector.size()-1 << "\n";
    }
    return fMuRSTrackVector.at(iTrack);
  }

  //########################################  
  int MuonRangeStackHits::GetArrivalTime( int iTrack ) const
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits") << "Requested arrival time for track number: " << iTrack << ". This doesn't exist for this event. The last track index you can access is: " << fMuRSTrackVector.size()-1 << "\n";
    }
    return fMuRSTrackVector.at(iTrack).arrivalTime;
  }

  //########################################  
  int MuonRangeStackHits::GetPenetrationDepth( int iTrack ) const
  {
    if( size_t(iTrack) >= fMuRSTrackVector.size() ){
      throw cet::exception("MuonRangeStackHits") << "Requested penetration depth for track number: " << iTrack << ". This doesn't exist for this event. The last track index you can access is: " << fMuRSTrackVector.size()-1 << "\n";
    }
    return fMuRSTrackVector.at(iTrack).penetrationDepth;
  }


  //########################################  
  size_t MuonRangeStackHits::NTracks(){ return fMuRSTrackVector.size(); }
  
}// end namespace
