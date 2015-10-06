//
//  ConditionsSummary.cpp
//  lariat-mrb
//
//  Created by Brian Rebel on 10/6/15.
//

#include <limits>

#include "cetlib/exception.h"

#include "LArIATDataProducts/ConditionsSummary.h"


namespace ldp {
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary()
  : fBeamOn       (false)
  , fMagnetCurrent(std::numeric_limits<float>::max())
  , fTPCHV        (std::numeric_limits<float>::max())
  , fMWPCVoltage  (std::vector<float>(4, std::numeric_limits<float>::max()))
  , fPMTHV        (std::numeric_limits<float>::max())
  , fSiPMHV       (std::numeric_limits<float>::max())
  {
    return;
  }
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary(bool               const& beamOn,
                                       float              const& magnetCurrent,
                                       float              const& tpchv,
                                       float              const& pmthv,
                                       float              const& sipmhv,
                                       std::vector<float> const& mwpcVoltages)
  : fBeamOn       (beamOn)
  , fMagnetCurrent(magnetCurrent)
  , fTPCHV        (tpchv)
  , fMWPCVoltage  (mwpcVoltages)
  , fPMTHV        (pmthv)
  , fSiPMHV       (sipmhv)
  {
    return;
  }

  //-------------------------------------------------------------------
  float const& ConditionsSummary::MWPCVoltage(size_t const& mwpc) const
  {
    if(mwpc > fMWPCVoltage.size() - 1)
      throw cet::exception("ConditionsSummary")
      << "requests voltage for mwpc: " << mwpc
      << " while only " << fMWPCVoltage.size()
      << " MWPCs in experiment";
    
    return fMWPCVoltage[mwpc];
  }
}
