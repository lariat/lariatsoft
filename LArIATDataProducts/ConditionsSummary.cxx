//
//  ConditionsSummary.cpp
//  lariat-mrb
//
//  Created by Brian Rebel on 10/6/15.
//

#include <limits>
#include <complex>

#include "cetlib/exception.h"

#include "LArIATDataProducts/ConditionsSummary.h"


namespace ldp {
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary()
  : fSecondaryIntensity    (std::numeric_limits<float>::max())
  , fSecondaryMomentum     (std::numeric_limits<float>::max())
  , fSecondaryPolarity     (std::numeric_limits<float>::max())
  , fMagnetCurrent         (std::numeric_limits<float>::max())
  , fMagnetPolarity        (std::numeric_limits<float>::max())
  , fTPCCathodeHV          (std::numeric_limits<float>::max())
  , fTPCCollectionV        (std::numeric_limits<float>::max())
  , fTPCInductionV         (std::numeric_limits<float>::max())
  , fTPCShieldV            (std::numeric_limits<float>::max())
  , fETLPMTHV              (std::numeric_limits<float>::max())
  , fHamamatsuPMTHV        (std::numeric_limits<float>::max())
  , fHamamatsuSiPMHV       (std::numeric_limits<float>::max())
  , fSenslSiPMHV           (std::numeric_limits<float>::max())
  , fTertiaryBeamCounters  (std::numeric_limits<float>::max())
  , fTertiaryCherenkov1    (std::numeric_limits<float>::max())
  , fTertiaryCherenkov2    (std::numeric_limits<float>::max())
  , fTertiaryCosmicCounters(std::numeric_limits<float>::max())
  , fDSTOF                 (std::numeric_limits<float>::max())
  , fUSTOF                 (std::numeric_limits<float>::max())
  , fHaloPaddle            (std::numeric_limits<float>::max())
  , fMuonRangeStack        (std::numeric_limits<float>::max())
  , fNumberMuRS            (std::numeric_limits<float>::max())
  , fPunchThrough          (std::numeric_limits<float>::max())
  , fMWPC                  (std::vector<bool>(4, false) )
  {
    return;
  }
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary(float              const& secondaryIntensity,
                                       float              const& secondaryMomentum,
                                       float              const& secondaryPolarity,
                                       float              const& magnetCurrent,
                                       float              const& magnetPolarity,
                                       float              const& tpcCathodeHV,
                                       float              const& tpcCollectionV,
                                       float              const& tpcInductionV,
                                       float              const& tpcShieldV,
                                       float              const& etlPMTHV,
                                       float              const& hamamatsuPMTHV,
                                       float              const& hamamatsuSiPMHV,
                                       float              const& senslSiPMHV,
                                       float              const& tertiaryBeamCounters,
                                       float              const& tertiaryCherenkov1,
                                       float              const& tertiaryCherenkov2,
                                       float              const& tertiaryCosmicCounters,
                                       float              const& dsTOF,
                                       float              const& usTOF,
                                       float              const& haloPaddle,
                                       float              const& muonRangeStack,
                                       float              const& numberMuRS,
                                       float              const& punchThrough,
                                       std::vector<bool>  const& mwpc)
  : fSecondaryIntensity    (secondaryIntensity    )
  , fSecondaryMomentum     (secondaryMomentum     )
  , fSecondaryPolarity     (secondaryPolarity     )
  , fMagnetCurrent         (magnetCurrent         )
  , fMagnetPolarity        (magnetPolarity        )
  , fTPCCathodeHV          (tpcCathodeHV          )
  , fTPCCollectionV        (tpcCollectionV        )
  , fTPCInductionV         (tpcInductionV         )
  , fTPCShieldV            (tpcShieldV            )
  , fETLPMTHV              (etlPMTHV              )
  , fHamamatsuPMTHV        (hamamatsuPMTHV        )
  , fHamamatsuSiPMHV       (hamamatsuSiPMHV       )
  , fSenslSiPMHV           (senslSiPMHV           )
  , fTertiaryBeamCounters  (tertiaryBeamCounters  )
  , fTertiaryCherenkov1    (tertiaryCherenkov1    )
  , fTertiaryCherenkov2    (tertiaryCherenkov2    )
  , fTertiaryCosmicCounters(tertiaryCosmicCounters)
  , fDSTOF                 (dsTOF                 )
  , fUSTOF                 (usTOF                 )
  , fHaloPaddle            (haloPaddle            )
  , fMuonRangeStack        (muonRangeStack        )
  , fNumberMuRS            (numberMuRS            )
  , fPunchThrough          (punchThrough          )
  , fMWPC                  (mwpc                  )
  {
    return;
  }

  //-------------------------------------------------------------------
  bool ConditionsSummary::MWPC(size_t const& mwpc) const
  {
    if(mwpc > fMWPC.size() - 1)
      throw cet::exception("ConditionsSummary")
      << "requests for mwpc: " << mwpc
      << " while only " << fMWPC.size()
      << " MWPCs in experiment";
    
    return fMWPC[mwpc];
  }
  
  //-------------------------------------------------------------------
  bool ConditionsSummary::BeamOn() const
  {
    return (fSecondaryIntensity > 0. && std::abs(fSecondaryMomentum) > 0.);
  }
  
} // end namespace
