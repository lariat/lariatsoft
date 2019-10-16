////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparateT1034.h
/// \brief Interface to algorithm class for a specific calculation of 
///        ionization electrons and scintillation photons assuming there
///        is no correlation between the two
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARIATLARG4_ISCALCULATIONSEPARATET1034_H
#define LARIATLARG4_ISCALCULATIONSEPARATET1034_H

#include <map>

#include "Geant4/G4EmSaturation.hh"

//#include "larsim/LArG4/ISCalculation.h"
#include "LArIATLArG4/ISCalculationT1034.h"

// forward declaration
namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

 class ISCalculationSeparateT1034 : public ISCalculationT1034 {

 public:

   ISCalculationSeparateT1034(CLHEP::HepRandomEngine&);
   virtual ~ISCalculationSeparateT1034();

   void   Initialize();
   void   Reset();
   void   CalculateIonizationAndScintillation(const G4Step* step);
   double StepSizeLimit()              const { return fStepSize;            }

 private:

   double                fStepSize;            ///< maximum step to take				  
   double                fEfield;              ///< value of electric field from LArProperties service
   double 	   	 fGeVToElectrons;      ///< conversion factor from LArProperties service	  
   double 	   	 fRecombA;             ///< from LArG4Parameters service			  
   double 	   	 fRecombk;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxA;             ///< from LArG4Parameters service			  
   double 	   	 fModBoxB;             ///< from LArG4Parameters service			  
   bool   	   	 fUseModBoxRecomb;     ///< from LArG4Parameters service			  
   bool   	   	 fScintByParticleType; ///< from LArProperties service			  
   double 	   	 fScintYieldFactor;    ///< scintillation yield factor                             
   G4EmSaturation* 	 fEMSaturation;        ///< pointer to EM saturation                            

 };
}
#endif // LARIATLARG4_ISCALCULATIONSEPARATET1034_H

