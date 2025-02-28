////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationT1034.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARG4_ISCALCULATION1034_H
#define LARG4_ISCALCULATION1034_H

#include "Geant4/G4Step.hh"
#include "Geant4/G4UserLimits.hh"

//#include "larsim/LArG4/OpFastScintillation.hh"
#include "LArIATLArG4/OpFastScintillationT1034.hh"

namespace larg4{

 class ISCalculationT1034{

 public:

   ISCalculationT1034();
   virtual ~ISCalculationT1034();

   virtual void                 Initialize()                                            = 0;
   virtual void                 Reset()                            		        = 0;
   virtual void                 CalculateIonizationAndScintillation(const G4Step* step) = 0;
   virtual double               StepSizeLimit()              const 		        = 0;

   double                       EnergyDeposit()              const { return fEnergyDeposit;   }
   double       	        NumberIonizationElectrons()  const { return fNumIonElectrons; }
   double       	        NumberScintillationPhotons() const { return fNumScintPhotons; }

   //Method to get electric field
   double EFieldAtStep(double fEfield, const G4Step* step) const; //value of field with any corrections for this step  

 protected:

   double fEnergyDeposit;   ///< total energy deposited in the step
   double fNumIonElectrons; ///< number of ionization electrons for this step
   double fNumScintPhotons; ///< number of scintillation photons for this step   
   
 };
}
#endif // LARG4_ISCALCULATION_H

