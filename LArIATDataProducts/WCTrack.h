////////////////////////////////////////////////////////////////////////
// $Id: WCTrack.h,v 1.00 2015/06/03 16:04:20 linehan3 Exp $
//
// Definition of wire chamber track object
//
// rlinehan@stanford.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_WCTRACK_H
#define LARIATDATAPRODUCTS_WCTRACK_H

#include <vector>
#include <iosfwd>
#include <string>


///Raw data description
namespace ldp {
  
  class WCTrack {

  public:
    WCTrack(); // Default constructor
    
  private:
    
    float fMomentum;                    //Reconstructed momentum in the XZ plane (coord system origin is at secondary target)
    float fYKink;                       //Angle difference between upstream and downstream tracks
    float fDeltaDist[3];                //Distance between upstream and downstream track ends
    float fXYFace[2];                   //X and Y position of the track on the upstream face of the TPC
    float fTheta;                       //Theta defined from the Z axis of the TPC
    float fPhi;                         //Phi defined counterclockwise from the X axis of the TPC
    
    //These are indexed by hit: each
    //hit is represented by the same
    //index in all three.
    std::vector<int> fWC;               //Defined from 1 to 4, like the Wire Chambers
    std::vector<float> fHitWire;        
    std::vector<float> fHitTime;
    
#ifndef __GCCXML__

  public:

    WCTrack( float momentum,
	     float yKink,
	     float xDist,
	     float yDist,
	     float zDist,
	     float xFace,
	     float yFace,
	     float theta,
	     float phi,
	     std::vector<int> wcVect,
	     std::vector<float> hitWireVect,
	     std::vector<float> hitTimeVect );
	     
    // Get Methods

    float               Momentum()                      const;
    float               YKink()                         const;
    float               DeltaDist(size_t i)             const;
    float               XYFace(size_t i)                const;
    float               Theta()                         const;
    float               Phi()                           const;
    int                 WC(size_t iHit)                 const;
    float               HitWire(size_t iHit)            const;
    float               HitTime(size_t iHit)            const;
    size_t              NHits()                         const;
    

#endif
  };
}

#ifndef __GCCXML__

inline float  ldp::WCTrack::Momentum() const { return fMomentum;      }
inline float  ldp::WCTrack::YKink()    const { return fYKink;         }
inline float  ldp::WCTrack::Theta()    const { return fTheta;         }
inline float  ldp::WCTrack::Phi()      const { return fPhi;           }
inline size_t ldp::WCTrack::NHits()    const { return fWC.size();     }

#endif

#endif // LARIATDATAPRODUCTS_WCTRACK_H

////////////////////////////////////////////////////////////////////////
