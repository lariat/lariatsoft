////////////////////////////////////////////////////////////////////////
// $Id: WCTrack.h,v 1.00 2015/09/19 16:04:20 linehan3 Exp $
//
// Definition of Aux Det Particle ID object
//
// rlinehan@stanford.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_AUXDETPARTICLEID_H
#define LARIATDATAPRODUCTS_AUXDETPARTICLEID_H

#include <vector>
#include <iosfwd>
#include <string>


///Raw data description
namespace ldp {
  
  class AuxDetParticleID {

  public:
    AuxDetParticleID();
    
  private:    
    
    //Taken from the event
    float fProtonProbability;
    float fKaonProbability;
    float fPiMuProbability;
    float fPionProbability;
    float fMuonProbability;
   
    //Made from the largest of the above
    int fPDGCode;
    

#ifndef __GCCXML__

  public:

    AuxDetParticleID( float p_prob,
		      float k_prob,
		      float pimu_prob,
		      float pi_prob,
		      float muon_prob,
		      float pdgcode );
	     
    // Get Methods
    float           ProtonProbability()      const;
    float           KaonProbability()        const;
    float           PiMuProbability()        const;
    float           PionProbability()        const;
    float           MuonProbability()        const;
    float           PDGCode()                const;


#endif
  };
}

#ifndef __GCCXML__

inline float  ldp::AuxDetParticleID::ProtonProbability() const { return fProtonProbability; }
inline float  ldp::AuxDetParticleID::KaonProbability()   const { return fKaonProbability; }
inline float  ldp::AuxDetParticleID::PiMuProbability()   const { return fPiMuProbability; }
inline float  ldp::AuxDetParticleID::PionProbability()   const { return fPionProbability;} 
inline float  ldp::AuxDetParticleID::MuonProbability()   const { return fMuonProbability; }
inline float  ldp::AuxDetParticleID::PDGCode()           const { return fPDGCode; }

#endif

#endif // LARIATDATAPRODUCTS_WCTRACK_H
