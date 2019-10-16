////////////////////////////////////////////////////////////////////////
/// \file   AuxDetReadoutT1034.h
/// \brief  A Geant4 sensitive detector that accumulates information.
/// \author miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARIATLARG4_AUXDETREADOUTT1034_H
#define LARIATLARG4_AUXDETREADOUTT1034_H

#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

#include <vector>

// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

namespace larg4 {
  
  class AuxDetReadoutT1034 : public G4VSensitiveDetector
  {
  public:
    // Constructor.
    AuxDetReadoutT1034(std::string const& name, 
		  unsigned int       adNum,
		  unsigned int       svNum);
    
    // Destructor
    virtual ~AuxDetReadoutT1034();
    
    // Required for classes that inherit from G4VSensitiveDetector.
    //
    // Called at start and end of each event.
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);
    
    // Called to clear any accumulated information.
    virtual void clear();
    
    // The key method of this class.  It's called by Geant4 for each
    // step within the read-out geometry.  It accumulates the energy
    // in the G4Step in the ?.
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory* );
    
    // Moved here from AuxDetSimChannel
    virtual void AddParticleStep(int	inputTrackID,
				 float	inputEnergyDeposited,
				 float	inputEntryX,
				 float	inputEntryY,
				 float	inputEntryZ,
				 float	inputEntryT,
				 float	inputExitX,
				 float	inputExitY,
				 float	inputExitZ,
				 float	inputExitT,
				 float	inputExitMomentumX,
				 float	inputExitMomentumY,
				 float	inputExitMomentumZ);

    // Empty methods; they have to be defined, but they're rarely
    // used in Geant4 applications.
    virtual void DrawAll();
    virtual void PrintAll();
    
    // Independent method; returns the accumulated information
    sim::AuxDetSimChannel const GetAuxDetSimChannel() const { return fAuxDetSimChannel; };
    
  private:
    art::ServiceHandle<geo::Geometry> fGeoHandle;        ///< Handle to the Geometry service
    uint32_t                          fAuxDet;           ///< which AuxDet this AuxDetReadoutT1034 corresponds to
    uint32_t                          fAuxDetSensitive;  ///< which sensitive volume of the AuxDet this AuxDetReadoutT1034 corresponds to
    sim::AuxDetSimChannel             fAuxDetSimChannel; ///< Contains the sim::AuxDetSimChannel for this AuxDet
    std::vector<sim::AuxDetIDE>       fAuxDetIDEs;       ///< list of IDEs in one channel
};
}

#endif // LARIATLARG4_AUXDETREADOUTT1034_H
