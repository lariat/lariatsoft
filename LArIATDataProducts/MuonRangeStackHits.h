//////////////////////////////////////////////////////////
// Definitions of MuonRangeStackHits object
// 
// gkpullia@syr.edu
// 
///////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H
#define LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H

#include <vector>
#include <iosfwd>
#include <string>
#include <map>
#include <iostream>

namespace ldp{

  struct MuRSTrack
  {
    std::vector<std::vector<int> > HitVect;          ///< A hit is a vector with 3 entries: 1.)
                                                     ///< Plane, 2.) Paddle, 3.) Time, and note
                                                     ///< that the combo (Plane = -1,  Paddle = -1)
                                                     ///< implies a Punchthrough hit
    int                            penetrationDepth; ///< How many planes deep did this go?
    int                            arrivalTime;      ///< At roughly what time did this track
                                                     ///< happen? (in units of the Hits' ticks)
  };

  class MuonRangeStackHits{
    
  public:
    //Constructors
    MuonRangeStackHits();    

  private:
    std::map<int, std::vector<int> > fPaddleTimeTickMap;
    std::vector<ldp::MuRSTrack>      fMuRSTrackVector;
    bool                             fIsInitializedEmpty;

#ifndef __GCCXML__

  public:
    // Non-default constructor
    MuonRangeStackHits(std::map<int, std::vector<int> > const& paddlemap,
                       std::vector<ldp::MuRSTrack>      const& trackVect );
    
    // Get Methods
    std::vector<int>       TimeTick(int iPaddle)  const; //The vector listing the time ticks when a certain iPaddle was hit.
    size_t                 NTracks();
    ldp::MuRSTrack         GetTrack(int iTrack);
    int                    GetPenetrationDepth(int iTrack) const; 
    int                    GetArrivalTime(int iTrack) const;
    bool                   WasItInitializedEmpty();

    // Set Methods
    void SetPenetrationDepth(int iTrack,
                             int depth);
    void SetArrivalTime(int iTrack, int time);

#endif
    
  };
}

#ifndef __GCCXML__

inline void   ldp::MuonRangeStackHits::SetPenetrationDepth(int iTrack, int depth)
{ fMuRSTrackVector.at(iTrack).penetrationDepth = depth; }

inline void   ldp::MuonRangeStackHits::SetArrivalTime(int iTrack, int time)
{ fMuRSTrackVector.at(iTrack).arrivalTime = time; }

inline size_t ldp::MuonRangeStackHits::NTracks()
{ return fMuRSTrackVector.size(); }

inline bool   ldp::MuonRangeStackHits::WasItInitializedEmpty()
{ return fIsInitializedEmpty; }

#endif // ___GCCXML__

#endif //LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H
	
