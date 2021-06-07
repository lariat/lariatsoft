////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationT1034.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

//#include "larsim/LArG4/ISCalculation.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "LArIATLArG4/ISCalculationT1034.h"

// Framework includes                                                                                                         
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace larg4 {

  //----------------------------------------------------------------------                                                    
  ISCalculationT1034::ISCalculationT1034()
  {
  }

  //----------------------------------------------------------------------                                                    
  ISCalculationT1034::~ISCalculationT1034()
  {
  }

  //......................................................................                                                    
  double ISCalculationT1034::EFieldAtStep(double fEfield, const G4Step* step) const
  {
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    if (!SCE->EnableSimEfieldSCE()) return fEfield;
    
    geo::Point_t midPoint
      { ( step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition() ) * 0.5/CLHEP::cm };
    auto EfieldDelta = fEfield * SCE->GetEfieldOffsets(midPoint);
    geo::Vector_t EfieldVec
      = { fEfield + EfieldDelta.X(), EfieldDelta.Y(), EfieldDelta.Z() };
    return EfieldVec.R();
  }

}
