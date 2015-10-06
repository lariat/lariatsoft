//
//  ConditionsSummary.h
//  lariat-mrb
//
//  Created by Brian Rebel on 10/6/15.
//  Copyright (c) 2015 Brian Rebel. All rights reserved.
//

#ifndef LARIATDATAPRODUCTS_ConditionsSummary_h
#define LARIATDATAPRODUCTS_ConditionsSummary_h

#include <vector>

namespace ldp {
  class ConditionsSummary {
  public:
    ConditionsSummary();
    
  private:
    
    float              fSecondaryIntensity;       ///< secondary beam intensity
    float              fSecondaryMomentum;        ///< secondary beam momentum
    float              fSecondaryPolarity;        ///< secondary beam polarity, > 0 --> Positive
    float              fMagnetCurrent;            ///< value of the current supplied to the magnet
    float              fMagnetPolarity;           ///< polarity of the magnet, > 0 --> Positive
    float              fTPCCathodeHV;             ///< High voltage on the cathode
    float              fTPCCollectionV;           ///< voltage on the collection plane
    float              fTPCInductionV;            ///< voltage on the induction plane
    float              fTPCShieldV;               ///< voltage on the shield plane
    float              fETLPMTHV;                 ///< HV for the ETL PMT
    float              fHamamatsuPMTHV;           ///< HV for the Hamamatsu PMT
    float              fHamamatsuSiPMHV;          ///< HV for the Hamamatsu SiPM
    float              fSenslSiPMHV;              ///< HV for the Sensl SiPM
    float              fTertiaryBeamCounters;     ///< something about the beam counters
    float              fTertiaryCherenkov1;       ///< something about cherenkov 1
    float              fTertiaryCherenkov2;       ///< something about cherenkov 2
    float              fTertiaryCosmicCounters;   ///< something about the cosmic counters
    float              fDSTOF;                    ///< something about the downstream Time of Flight
    float              fUSTOF;                    ///< something about the downstream Time of Flight
    float              fHaloPaddle;               ///< something about the halo paddle
    float              fMuonRangeStack;           ///< something about the muon range stack
    float              fNumberMuRS;               ///< number of paddles in muon range stack?
    float              fPunchThrough;             ///< something about the punch through
    std::vector<bool>  fMWPC;                     ///< which MWPCs are on

#ifndef __GCCXML__
    
  public:
    
    ConditionsSummary(float              const& secondaryIntensity,
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
                      std::vector<bool> const& mwpc);
    
    bool         BeamOn()                        const;
    float const& SecondaryIntensity()            const;
    float const& SecondaryMomentum()             const;
    float const& SecondaryPolarity()             const;
    float const& MagnetCurrent()                 const;
    float const& MagnetPolarity()                const;
    float const& TPCCathodeHV()                  const;
    float const& TPCCollectionV()                const;
    float const& TPCInductionV()                 const;
    float const& TPCShieldV()                    const;
    bool         MWPC(size_t const& mwpc)        const;
    float const& ETLPMTHV()                      const;
    float const& HamamatsuPMTHV()                const;
    float const& HamamatsuSiPMHV()               const;
    float const& SenslSiPMHV()                   const;
    float const& TertiaryBeamCounters()          const;
    float const& TertiaryCherenkov1()            const;
    float const& TertiaryCherenkov2()            const;
    float const& TertiaryCosmicCounters()        const;
    float const& DownStreamTOF()                 const;
    float const& UpStreamTOF()                   const;
    float const& HaloPaddle()                    const;
    float const& MuonRangeStack()                const;
    float const& NumberMuRS()                    const;
    float const& PunchThrough()                  const;
    
#endif
    
  };
} // end namespace

#ifndef __GCCXML__

inline float const& ldp::ConditionsSummary::SecondaryIntensity()     const { return fSecondaryIntensity;     }
inline float const& ldp::ConditionsSummary::SecondaryMomentum()      const { return fSecondaryMomentum;      }
inline float const& ldp::ConditionsSummary::SecondaryPolarity()      const { return fSecondaryPolarity;      }
inline float const& ldp::ConditionsSummary::MagnetCurrent()          const { return fMagnetCurrent;          }
inline float const& ldp::ConditionsSummary::MagnetPolarity()         const { return fMagnetPolarity;         }
inline float const& ldp::ConditionsSummary::TPCCathodeHV()           const { return fTPCCathodeHV;           }
inline float const& ldp::ConditionsSummary::TPCCollectionV()         const { return fTPCCollectionV;         }
inline float const& ldp::ConditionsSummary::TPCInductionV()          const { return fTPCInductionV;          }
inline float const& ldp::ConditionsSummary::TPCShieldV()             const { return fTPCShieldV;             }
inline float const& ldp::ConditionsSummary::ETLPMTHV()               const { return fETLPMTHV;               }
inline float const& ldp::ConditionsSummary::HamamatsuPMTHV()         const { return fHamamatsuPMTHV;         }
inline float const& ldp::ConditionsSummary::HamamatsuSiPMHV()        const { return fHamamatsuSiPMHV;        }
inline float const& ldp::ConditionsSummary::SenslSiPMHV()            const { return fSenslSiPMHV;            }
inline float const& ldp::ConditionsSummary::TertiaryBeamCounters()   const { return fTertiaryBeamCounters;   }
inline float const& ldp::ConditionsSummary::TertiaryCherenkov1()     const { return fTertiaryCherenkov1;     }
inline float const& ldp::ConditionsSummary::TertiaryCherenkov2()     const { return fTertiaryCherenkov2;     }
inline float const& ldp::ConditionsSummary::TertiaryCosmicCounters() const { return fTertiaryCosmicCounters; }
inline float const& ldp::ConditionsSummary::DownStreamTOF()          const { return fDSTOF;                  }
inline float const& ldp::ConditionsSummary::UpStreamTOF()            const { return fUSTOF;                  }
inline float const& ldp::ConditionsSummary::HaloPaddle()             const { return fHaloPaddle;             }
inline float const& ldp::ConditionsSummary::MuonRangeStack()         const { return fMuonRangeStack;         }
inline float const& ldp::ConditionsSummary::NumberMuRS()             const { return fNumberMuRS;             }
inline float const& ldp::ConditionsSummary::PunchThrough()           const { return fPunchThrough;           }

#endif

#endif //LARIATDATAPRODUCTS_ConditionsSummary_h
