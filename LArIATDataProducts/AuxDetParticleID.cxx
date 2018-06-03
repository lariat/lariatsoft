////////////////////////////////////////////////////////////////////////
// $Id: WCTrack.cxx,v 1.00 2015/06/03 16:04:20 linehan3 Exp $
//
// AuxDetParticleID class
//
// rlinehan@stanford.edu
//
////////////////////////////////////////////////////////////////////////

#include "LArIATDataProducts/AuxDetParticleID.h"
#include "cetlib_except/exception.h"

namespace ldp{

  //----------------------------------------------------------------------
  AuxDetParticleID::AuxDetParticleID()
  {
    fProtonProbability = 0;
    fKaonProbability = 0;
    fPiMuProbability = 0;
    fPionProbability = 0;
    fMuonProbability = 0;
    fPDGCode = 0;

  }

  //----------------------------------------------------------------------
  AuxDetParticleID::AuxDetParticleID( float p_prob,
				      float k_prob,
				      float pimu_prob,
				      float pi_prob,
				      float mu_prob,
				      float pdgcode )
  { 
    fProtonProbability = p_prob;
    fKaonProbability = k_prob;
    fPiMuProbability = pimu_prob;
    fPionProbability = pi_prob;
    fMuonProbability = mu_prob;
    fPDGCode = pdgcode;
  }
}
   
////////////////////////////////////////////////////////////////////////

