////////////////////////////////////////////////////////////////////////
// $Id: WCTrack.cxx,v 1.00 2015/06/03 16:04:20 linehan3 Exp $
//
// WCTrack class
//
// rlinehan@stanford.edu
//
////////////////////////////////////////////////////////////////////////

#include "LArIATDataProducts/WCTrack.h"
#include "cetlib/exception.h"

namespace ldp{

  //----------------------------------------------------------------------
  WCTrack::WCTrack() 
  {
    std::vector<int> blank1;
    std::vector<float> blank2;
    fMomentum = 0;
    fMomentum2M = 0;
    fYKink = 0;
    for( int i = 0; i < 3 ; ++i )
      fDeltaDist[i] = 0;
    for( int i = 0; i < 2 ; ++i )
      fXYFace[i] = 0;
    fTheta = 0;
    fPhi = 0;
    fWC = blank1;
    fHitWire = blank2;
    //fHitTime = blank2;
  }

  //----------------------------------------------------------------------
  WCTrack::WCTrack(float momentum,
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
		   float hitPositionVect[4][3],
		   int WCMissed,
		   float residual,
		   float WC1XMult,
		   float WC2XMult,
		   float WC3XMult,
		   float WC4XMult,
		   float WC1YMult,
		   float WC2YMult,
		   float WC3YMult,
		   float WC4YMult,		   
		   bool PickyTrackCheck,
		   float unscaledmomentum)
		   
  { 
    fMomentum = momentum;
    fYKink = yKink;
    fDeltaDist[0] = xDist; fDeltaDist[1] = yDist; fDeltaDist[2] = zDist;
    fXYFace[0] = xFace; fXYFace[1] = yFace;
    fTheta = theta;
    fPhi = phi;
    fWC = wcVect;
    fHitWire = hitWireVect;
    for(int i=0; i<4; ++i){
      for(int j=0; j<3; ++j){
    fHitPosition[i][j]= hitPositionVect[i][j];
      }
    }
    fWCMissed=WCMissed;
    fResidual=residual;
    fWC1XMult=WC1XMult;
    fWC2XMult=WC2XMult;
    fWC3XMult=WC3XMult;
    fWC4XMult=WC4XMult;
    fWC1YMult=WC1YMult;
    fWC2YMult=WC2YMult;
    fWC3YMult=WC3YMult;
    fWC4YMult=WC4YMult;    
    fPickyTrackCheck=PickyTrackCheck; 
    fUnscaledMomentum=unscaledmomentum;   
    //fHitTime = hitTimeVect;
    fDownstreamDir.SetXYZ(hitPositionVect[3][0] - hitPositionVect[2][0],
                          hitPositionVect[3][1] - hitPositionVect[2][1],
                          hitPositionVect[3][2] - hitPositionVect[2][2]);
  }

  WCTrack::WCTrack(float momentum,
                   float momentum2m,
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
		   float hitPositionVect[4][3],
		   int WCMissed,
		   float residual,
		   float WC1XMult,
		   float WC2XMult,
		   float WC3XMult,
		   float WC4XMult,
		   float WC1YMult,
		   float WC2YMult,
		   float WC3YMult,
		   float WC4YMult,
		   bool PickyTrackCheck,
		   float unscaledmomentum)
  { 
    fMomentum = momentum;
    fMomentum2M = momentum2m;
    fYKink = yKink;
    fDeltaDist[0] = xDist; fDeltaDist[1] = yDist; fDeltaDist[2] = zDist;
    fXYFace[0] = xFace; fXYFace[1] = yFace;
    fTheta = theta;
    fPhi = phi;
    fWC = wcVect;
    fHitWire = hitWireVect;
    for(int i=0; i<4; ++i){
      for(int j=0; j<3; ++j){
    fHitPosition[i][j]= hitPositionVect[i][j];
      }
    }
    fWCMissed=WCMissed;
    fResidual=residual;
    fWC1XMult=WC1XMult;
    fWC2XMult=WC2XMult;
    fWC3XMult=WC3XMult;
    fWC4XMult=WC4XMult;
    fWC1YMult=WC1YMult;
    fWC2YMult=WC2YMult;
    fWC3YMult=WC3YMult;
    fWC4YMult=WC4YMult; 
    fPickyTrackCheck=PickyTrackCheck;
    fUnscaledMomentum=unscaledmomentum;
    //fHitTime = hitTimeVect;
    fDownstreamDir.SetXYZ(hitPositionVect[3][0] - hitPositionVect[2][0],
                          hitPositionVect[3][1] - hitPositionVect[2][1],
                          hitPositionVect[3][2] - hitPositionVect[2][2]);
	
    fPz=momentum*cos(theta);
    fPy=momentum*sin(theta)*sin(phi);
    fPx=momentum*sin(theta)*cos(phi);		  
			  
			  
  }

  //--------------------------------------------------
float WCTrack::DeltaDist(size_t i) const
{
 if( i >= 3 || i < 0 ){
   throw cet::exception("WCTrack") << "illegal index requested for DeltaDist vector: "
<< i << "\n";
}
return fDeltaDist[i];
}
  
  //--------------------------------------------------
  float WCTrack::XYFace(size_t i) const
  {
    if( i >= 2 || i < 0 ){
      throw cet::exception("WCTrack") << "illegal index requested for XYFace vector: "
				      << i << "\n";
    }
    return fXYFace[i];
  }
 
  //--------------------------------------------------
  int WCTrack::WC(size_t iHit) const
  {
    if( iHit >= fWC.size() ){
      throw cet::exception("WCTrack") << "illegal index requested for WC vector: "
				      << iHit << "\n";
    }
    return fWC[iHit];
  }

  //--------------------------------------------------
  float WCTrack::HitWire(size_t iHit) const
  {
    if( iHit >= fHitWire.size() ){
      throw cet::exception("WCTrack") << "illegal index requested for HitWire vector: "
			 	      << iHit << "\n";
    }
    return fHitWire[iHit];
  }
//=====================================================  
  float WCTrack::HitPosition(int iWC, int iAx) const
  {
    if(iWC >3 ){
      throw cet::exception("WCTrack") <<"illegal WC index requested for HitPosition: "
      				      << iWC << "\n";
    }
    if(iAx >2){
      throw cet::exception("WCTrack") <<"illegal dimension index requested for HitPosition: "
      				      << iAx << "\n";
    }
    return fHitPosition[iWC][iAx];
  }

  //--------------------------------------------------
  //float WCTrack::HitTime(size_t iHit) const
 // {
   // if( iHit >= fHitTime.size() ){
     // throw cet::exception("WCTrack") << "illegal index requested for HitTime vector: "
				     // << iHit << "\n";
 //   }
  //  return fHitTime[iHit];
  //}

  //---------------------------------------------------------------------------
  TVector3 WCTrack::ProjectionAtZ(float z, bool useXYFace) const
  {
    /*--------------------------------------------------------

    Return the projected position of the downstream WCTrack
    at a given z (cm).

    --------------------------------------------------------*/

    // get position of the hit on WC4
    float const wc4x = fHitPosition[3][0];
    float const wc4y = fHitPosition[3][1];
    float const wc4z = fHitPosition[3][2];

    if (useXYFace)
    {
      // projected position at the front face of the TPC
      float const x0 = fXYFace[0];
      float const y0 = fXYFace[1];

      // get direction vector from WC4 to front face of TPC
      float const x2 = x0 - wc4x;
      float const y2 = y0 - wc4y;
      float const z2 =  0 - wc4z;

      // project to z
      float const slope = (z - wc4z) / z2;
      float const x = wc4x + slope * x2;
      float const y = wc4y + slope * y2;

      return TVector3(x, y, z);
    }

    // get direction of downstream WCTrack
    float const x1 = fDownstreamDir.X();
    float const y1 = fDownstreamDir.Y();
    float const z1 = fDownstreamDir.Z();

    // project to z
    float const slope = (z - wc4z) / z1;
    float const x = wc4x + slope * x1;
    float const y = wc4y + slope * y1;

    return TVector3(x, y, z);
  }

}
////////////////////////////////////////////////////////////////////////

