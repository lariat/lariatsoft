////////////////////////////////////////////////////////////////////////
// $Id: TOF.cxx,v 1.00 2015/06/03 16:04:20 dsmith Exp $
//
// TOF class
//
// dansmith@bu.edu
//
////////////////////////////////////////////////////////////////////////

#include "LArIATDataProducts/TOF.h"
#include "cetlib_except/exception.h"

namespace ldp{

  //----------------------------------------------------------------------
  TOF::TOF() {

    std::vector<float> blank0;
    std::vector<long> blank1;
    fTOF = blank0;
    fTimeStamp = blank1;
  }

 
  //--------------------------------------------------
  float TOF::SingleTOF(size_t iHit) const
  {
    if( iHit >= fTOF.size() ){
      throw cet::exception("TOF") << "illegal index requested for TOF vector: "
				      << iHit << "\n";
    }
    return fTOF[iHit];
  }

  //--------------------------------------------------
  long TOF::TimeStamp(size_t iHit) const
  {
    if( iHit >= fTimeStamp.size() ){
      throw cet::exception("TOF") << "illegal index requested for TimeStamp vector: "
			 	      << iHit << "\n";
    }
    return fTimeStamp[iHit];
  }

}
////////////////////////////////////////////////////////////////////////

