////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNESTT1034.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using nest
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4EmSaturation.hh"

//#include "larsim/LArG4/ISCalculationNEST.h"
#include "LArIATLArG4/ISCalculationNESTT1034.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4{

  //----------------------------------------------------------------------------
  ISCalculationNESTT1034::ISCalculationNESTT1034(CLHEP::HepRandomEngine& engine)
    : fNest(0)
    , fEngine(engine)
  {
    return;
  }

  //----------------------------------------------------------------------------
  ISCalculationNESTT1034::~ISCalculationNESTT1034()
  {
    if(fNest) delete fNest;

    return;
  }

  //----------------------------------------------------------------------------
  void ISCalculationNESTT1034::Initialize()
  {
    // \todo should ideally make the yield factor passed to the NestAlg ctor a parameter
    if(!fNest) fNest = new NestAlg(1., fEngine);

    // Set the step size to small value if NEST is chosen, per Matthew Szydagis, 
    // "because without delta rays, the yields are wrong.  The ICARUS model that is 
    // in LArSoft uses a fudge factor to compensate, but NEST is "purer" -- no 
    // fudge factor. "
    fStepSize = 0.05 * CLHEP::micrometer;

    return;
  }
  
  //----------------------------------------------------------------------------
  void ISCalculationNESTT1034::Reset()
  {
    fEnergyDeposit   = 0.;
    fNumIonElectrons = 0.;
    fNumScintPhotons = 0.;

    return;
  }

  //----------------------------------------------------------------------------
  void ISCalculationNESTT1034::CalculateIonizationAndScintillation(const G4Step* step)
  {
    // get a const representation of the track for this step
    const G4Track track(*(step->GetTrack()));
    
    fNest->CalculateIonizationAndScintillation(track, *step);

    // compare the energy deposition of this step to what is in the fNest object
    if(fNest->EnergyDeposition() != step->GetTotalEnergyDeposit()/CLHEP::MeV)
      mf::LogWarning("ISCalculationNestT1034") << "NEST and G4 step depositions do not agree!\n"
					  << fNest->EnergyDeposition() << " vs " 
					  << step->GetTotalEnergyDeposit()/CLHEP::MeV;

    // Nest uses Geant units, LArSoft assumes energy is in units of MeV here
    fEnergyDeposit   = fNest->EnergyDeposition()/CLHEP::MeV;
    fNumIonElectrons = fNest->NumberIonizationElectrons();
    fNumScintPhotons = fNest->NumberScintillationPhotons();

    return;
  }

}// namespace
