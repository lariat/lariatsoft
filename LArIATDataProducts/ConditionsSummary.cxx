//
//  ConditionsSummary.cpp
//  lariat-mrb
//
//  Created by Brian Rebel on 10/6/15.
//

#include <limits>
#include <complex>

#include "cetlib_except/exception.h"

#include "LArIATDataProducts/ConditionsSummary.h"


namespace ldp {
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary()
//  : fSecondaryIntensity    (std::numeric_limits<unsigned int>::max())
//  , fSecondaryMomentum     (std::numeric_limits<unsigned int>::max())
//  , fSecondaryPolarity     (std::numeric_limits<unsigned int>::max())
//  , fMagnetCurrent         (std::numeric_limits<unsigned int>::max())
//  , fMagnetPolarity        (std::numeric_limits<unsigned int>::max())
//  , fTPCCathodeHV          (std::numeric_limits<unsigned int>::max())
//  , fTPCCollectionV        (std::numeric_limits<unsigned int>::max())
//  , fTPCInductionV         (std::numeric_limits<unsigned int>::max())
//  , fTPCShieldV            (std::numeric_limits<unsigned int>::max())
//  , fETLPMTHV              (std::numeric_limits<unsigned int>::max())
//  , fHamamatsuPMTHV        (std::numeric_limits<unsigned int>::max())
//  , fHamamatsuSiPMHV       (std::numeric_limits<unsigned int>::max())
//  , fSenslSiPMHV           (std::numeric_limits<unsigned int>::max())
//  , fTertiaryBeamCounters  (std::numeric_limits<unsigned int>::max())
//  , fTertiaryCherenkov1    (std::numeric_limits<unsigned int>::max())
//  , fTertiaryCherenkov2    (std::numeric_limits<unsigned int>::max())
//  , fTertiaryCosmicCounters(std::numeric_limits<unsigned int>::max())
//  , fDSTOF                 (std::numeric_limits<unsigned int>::max())
//  , fUSTOF                 (std::numeric_limits<unsigned int>::max())
//  , fHaloPaddle            (std::numeric_limits<unsigned int>::max())
//  , fMuonRangeStack        (std::numeric_limits<unsigned int>::max())
//  , fNumberMuRS            (std::numeric_limits<unsigned int>::max())
//  , fPunchThrough          (std::numeric_limits<unsigned int>::max())
//  , fCorrectFileFormat     (false                                   )
//  , fMWPC                  (std::vector<bool>(4, false)             )
  : fEndMC7SC1             (std::numeric_limits<unsigned int>::max())
  , fV1751CaenEnableReadout(false                                   )
  , fASICCollectionFilter  (std::numeric_limits<unsigned int>::max())
  , fASICCollectionGain    (std::numeric_limits<unsigned int>::max())
  , fASICEnableReadout     (false                                   )
  , fASICPulserOn          (false                                   )
  , fASICChannelScan       (false                                   )
  , fV1740RecordLength     (std::numeric_limits<unsigned int>::max())
  {
    return;
  }
  
  //-------------------------------------------------------------------
  ConditionsSummary::ConditionsSummary(/*size_t             const& secondaryIntensity,
                                       size_t             const& secondaryMomentum,
                                       size_t             const& secondaryPolarity,
                                       size_t             const& magnetCurrent,
                                       size_t             const& magnetPolarity,
                                       size_t             const& tpcCathodeHV,
                                       size_t             const& tpcCollectionV,
                                       size_t             const& tpcInductionV,
                                       size_t             const& tpcShieldV,
                                       size_t             const& etlPMTHV,
                                       size_t             const& hamamatsuPMTHV,
                                       size_t             const& hamamatsuSiPMHV,
                                       size_t             const& senslSiPMHV,
                                       size_t             const& tertiaryBeamCounters,
                                       size_t             const& tertiaryCherenkov1,
                                       size_t             const& tertiaryCherenkov2,
                                       size_t             const& tertiaryCosmicCounters,
                                       size_t             const& dsTOF,
                                       size_t             const& usTOF,
                                       size_t             const& haloPaddle,
                                       size_t             const& muonRangeStack,
                                       size_t             const& numberMuRS,
                                       size_t             const& punchThrough,
                                       bool               const& correctFileFormat,
                                       std::vector<bool>  const& mwpc*/
                                       size_t             const& endMC7SC1,
                                       bool               const& v1751EnableReadout,
                                       size_t             const& asicCollectionFilter,
                                       size_t             const& asicCollectionGain,
                                       bool               const& asicEnableReadout,
                                       bool               const& asicPulserOn,
                                       bool               const& asicChannelScan,
                                       size_t             const& v1740RecordLength)
//  : fSecondaryIntensity    (secondaryIntensity    )
//  , fSecondaryMomentum     (secondaryMomentum     )
//  , fSecondaryPolarity     (secondaryPolarity     )
//  , fMagnetCurrent         (magnetCurrent         )
//  , fMagnetPolarity        (magnetPolarity        )
//  , fTPCCathodeHV          (tpcCathodeHV          )
//  , fTPCCollectionV        (tpcCollectionV        )
//  , fTPCInductionV         (tpcInductionV         )
//  , fTPCShieldV            (tpcShieldV            )
//  , fETLPMTHV              (etlPMTHV              )
//  , fHamamatsuPMTHV        (hamamatsuPMTHV        )
//  , fHamamatsuSiPMHV       (hamamatsuSiPMHV       )
//  , fSenslSiPMHV           (senslSiPMHV           )
//  , fTertiaryBeamCounters  (tertiaryBeamCounters  )
//  , fTertiaryCherenkov1    (tertiaryCherenkov1    )
//  , fTertiaryCherenkov2    (tertiaryCherenkov2    )
//  , fTertiaryCosmicCounters(tertiaryCosmicCounters)
//  , fDSTOF                 (dsTOF                 )
//  , fUSTOF                 (usTOF                 )
//  , fHaloPaddle            (haloPaddle            )
//  , fMuonRangeStack        (muonRangeStack        )
//  , fNumberMuRS            (numberMuRS            )
//  , fPunchThrough          (punchThrough          )
//  , fCorrectFileFormat     (correctFileFormat     )
//  , fMWPC                  (mwpc                  )
  : fEndMC7SC1             (endMC7SC1             )
  , fV1751CaenEnableReadout(v1751EnableReadout    )
  , fASICCollectionFilter  (asicCollectionFilter  )
  , fASICCollectionGain    (asicCollectionGain    )
  , fASICEnableReadout     (asicEnableReadout     )
  , fASICPulserOn          (asicPulserOn          )
  , fASICChannelScan       (asicChannelScan       )
  , fV1740RecordLength     (v1740RecordLength     )
  {
    return;
  }

//  //-------------------------------------------------------------------
//  bool ConditionsSummary::MWPC(size_t const& mwpc) const
//  {
//    if(mwpc > fMWPC.size() - 1)
//      throw cet::exception("ConditionsSummary")
//      << "requests for mwpc: " << mwpc
//      << " while only " << fMWPC.size()
//      << " MWPCs in experiment";
//    
//    return fMWPC[mwpc];
//  }
//  
//  //-------------------------------------------------------------------
//  bool ConditionsSummary::BeamOn() const
//  {
//    return (fSecondaryIntensity > 0. && std::abs(fSecondaryMomentum) > 0.);
//  }
  
} // end namespace
