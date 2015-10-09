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
    
//    size_t             fSecondaryIntensity;       ///< secondary beam intensity
//    size_t             fSecondaryMomentum;        ///< secondary beam momentum
//    size_t             fSecondaryPolarity;        ///< secondary beam polarity, > 0 --> Positive
//    size_t             fMagnetCurrent;            ///< value of the current supplied to the magnet
//    size_t             fMagnetPolarity;           ///< polarity of the magnet, > 0 --> Positive
//    size_t             fTPCCathodeHV;             ///< High voltage on the cathode
//    size_t             fTPCCollectionV;           ///< voltage on the collection plane
//    size_t             fTPCInductionV;            ///< voltage on the induction plane
//    size_t             fTPCShieldV;               ///< voltage on the shield plane
//    size_t             fETLPMTHV;                 ///< HV for the ETL PMT
//    size_t             fHamamatsuPMTHV;           ///< HV for the Hamamatsu PMT
//    size_t             fHamamatsuSiPMHV;          ///< HV for the Hamamatsu SiPM
//    size_t             fSenslSiPMHV;              ///< HV for the Sensl SiPM
//    size_t             fTertiaryBeamCounters;     ///< something about the beam counters
//    size_t             fTertiaryCherenkov1;       ///< something about cherenkov 1
//    size_t             fTertiaryCherenkov2;       ///< something about cherenkov 2
//    size_t             fTertiaryCosmicCounters;   ///< something about the cosmic counters
//    size_t             fDSTOF;                    ///< something about the downstream Time of Flight
//    size_t             fUSTOF;                    ///< something about the downstream Time of Flight
//    size_t             fHaloPaddle;               ///< something about the halo paddle
//    size_t             fMuonRangeStack;           ///< something about the muon range stack
//    size_t             fNumberMuRS;               ///< number of paddles in muon range stack?
//    size_t             fPunchThrough;             ///< something about the punch through
//    bool               fCorrectFileFormat;        ///< ARTDAQ or other
//    std::vector<bool>  fMWPC;                     ///< which MWPCs are on
    size_t             fEndMC7SC1;                ///< MC7SC1
    bool               fV1751CaenEnableReadout;   ///< Was the 1751 Readout enabled
    size_t             fASICCollectionFilter;     ///< ASIC collection filter
    size_t             fASICCollectionGain;       ///< ASIC collection gain
    bool               fASICEnableReadout;        ///< ASIC enable Readout
    bool               fASICPulserOn;             ///< ASIC pulser on
    bool               fASICChannelScan;          ///< ASIC channel scan
    size_t             fV1740RecordLength;        ///< 1740 Record length

#ifndef __GCCXML__
    
  public:
    
    // comment out the first chunk of input arguments until we know how to
    // get them from SAM
    ConditionsSummary(/*size_t              const& secondaryIntensity,
                      size_t              const& secondaryMomentum,
                      size_t              const& secondaryPolarity,
                      size_t              const& magnetCurrent,
                      size_t              const& magnetPolarity,
                      size_t              const& tpcCathodeHV,
                      size_t              const& tpcCollectionV,
                      size_t              const& tpcInductionV,
                      size_t              const& tpcShieldV,
                      size_t              const& etlPMTHV,
                      size_t              const& hamamatsuPMTHV,
                      size_t              const& hamamatsuSiPMHV,
                      size_t              const& senslSiPMHV,
                      size_t              const& tertiaryBeamCounters,
                      size_t              const& tertiaryCherenkov1,
                      size_t              const& tertiaryCherenkov2,
                      size_t              const& tertiaryCosmicCounters,
                      size_t              const& dsTOF,
                      size_t              const& usTOF,
                      size_t              const& haloPaddle,
                      size_t              const& muonRangeStack,
                      size_t              const& numberMuRS,
                      size_t              const& punchThrough,,
                      bool                const& correctFileFormat,
                      std::vector<bool>   const& mwpc*/
                      size_t              const& endMC7SC1,
                      bool                const& v1751CaenEnableReadout,
                      size_t              const& asicCollectionFilter,
                      size_t              const& asicCollectionGain,
                      bool                const& asicEnableReadout,
                      bool                const& asicPulserOn,
                      bool                const& asicChannelScan,
                      size_t              const& v1740RecordLength);
    
    bool          BeamOn()                        const;
//    size_t const& SecondaryIntensity()            const;
//    size_t const& SecondaryMomentum()             const;
//    size_t const& SecondaryPolarity()             const;
//    size_t const& MagnetCurrent()                 const;
//    size_t const& MagnetPolarity()                const;
//    size_t const& TPCCathodeHV()                  const;
//    size_t const& TPCCollectionV()                const;
//    size_t const& TPCInductionV()                 const;
//    size_t const& TPCShieldV()                    const;
//    bool          MWPC(size_t const& mwpc)        const;
//    size_t const& ETLPMTHV()                      const;
//    size_t const& HamamatsuPMTHV()                const;
//    size_t const& HamamatsuSiPMHV()               const;
//    size_t const& SenslSiPMHV()                   const;
//    size_t const& TertiaryBeamCounters()          const;
//    size_t const& TertiaryCherenkov1()            const;
//    size_t const& TertiaryCherenkov2()            const;
//    size_t const& TertiaryCosmicCounters()        const;
//    size_t const& DownStreamTOF()                 const;
//    size_t const& UpStreamTOF()                   const;
//    size_t const& HaloPaddle()                    const;
//    size_t const& MuonRangeStack()                const;
//    size_t const& NumberMuRS()                    const;
//    size_t const& PunchThrough()                  const;
    size_t const& EndMC7SC1()                     const;
    bool   const& V1751CAENEnableReadout()        const;
    size_t const& ASICCollectionFilter()          const;
    size_t const& ASICCollectionGain()            const;
    bool   const& ASICEnableReadout()             const;
    bool   const& ASICPulserOn()                  const;
    bool   const& ASICChannelScan()               const;
    size_t const& V1740RecordLength()             const;
    
#endif
    
  };
} // end namespace

#ifndef __GCCXML__

//inline size_t const& ldp::ConditionsSummary::SecondaryIntensity()     const { return fSecondaryIntensity;     }
//inline size_t const& ldp::ConditionsSummary::SecondaryMomentum()      const { return fSecondaryMomentum;      }
//inline size_t const& ldp::ConditionsSummary::SecondaryPolarity()      const { return fSecondaryPolarity;      }
//inline size_t const& ldp::ConditionsSummary::MagnetCurrent()          const { return fMagnetCurrent;          }
//inline size_t const& ldp::ConditionsSummary::MagnetPolarity()         const { return fMagnetPolarity;         }
//inline size_t const& ldp::ConditionsSummary::TPCCathodeHV()           const { return fTPCCathodeHV;           }
//inline size_t const& ldp::ConditionsSummary::TPCCollectionV()         const { return fTPCCollectionV;         }
//inline size_t const& ldp::ConditionsSummary::TPCInductionV()          const { return fTPCInductionV;          }
//inline size_t const& ldp::ConditionsSummary::TPCShieldV()             const { return fTPCShieldV;             }
//inline size_t const& ldp::ConditionsSummary::ETLPMTHV()               const { return fETLPMTHV;               }
//inline size_t const& ldp::ConditionsSummary::HamamatsuPMTHV()         const { return fHamamatsuPMTHV;         }
//inline size_t const& ldp::ConditionsSummary::HamamatsuSiPMHV()        const { return fHamamatsuSiPMHV;        }
//inline size_t const& ldp::ConditionsSummary::SenslSiPMHV()            const { return fSenslSiPMHV;            }
//inline size_t const& ldp::ConditionsSummary::TertiaryBeamCounters()   const { return fTertiaryBeamCounters;   }
//inline size_t const& ldp::ConditionsSummary::TertiaryCherenkov1()     const { return fTertiaryCherenkov1;     }
//inline size_t const& ldp::ConditionsSummary::TertiaryCherenkov2()     const { return fTertiaryCherenkov2;     }
//inline size_t const& ldp::ConditionsSummary::TertiaryCosmicCounters() const { return fTertiaryCosmicCounters; }
//inline size_t const& ldp::ConditionsSummary::DownStreamTOF()          const { return fDSTOF;                  }
//inline size_t const& ldp::ConditionsSummary::UpStreamTOF()            const { return fUSTOF;                  }
//inline size_t const& ldp::ConditionsSummary::HaloPaddle()             const { return fHaloPaddle;             }
//inline size_t const& ldp::ConditionsSummary::MuonRangeStack()         const { return fMuonRangeStack;         }
//inline size_t const& ldp::ConditionsSummary::NumberMuRS()             const { return fNumberMuRS;             }
//inline size_t const& ldp::ConditionsSummary::PunchThrough()           const { return fPunchThrough;           }
inline size_t const& ldp::ConditionsSummary::EndMC7SC1()              const { return fEndMC7SC1;              }
inline bool   const& ldp::ConditionsSummary::V1751CAENEnableReadout() const { return fV1751CaenEnableReadout; }
inline size_t const& ldp::ConditionsSummary::ASICCollectionFilter()   const { return fASICCollectionFilter;   }
inline size_t const& ldp::ConditionsSummary::ASICCollectionGain()     const { return fASICCollectionGain;     }
inline bool   const& ldp::ConditionsSummary::ASICEnableReadout()      const { return fASICEnableReadout;      }
inline bool   const& ldp::ConditionsSummary::ASICPulserOn()           const { return fASICPulserOn;           }
inline bool   const& ldp::ConditionsSummary::ASICChannelScan()        const { return fASICChannelScan;        }
inline size_t const& ldp::ConditionsSummary::V1740RecordLength()      const { return fV1740RecordLength;      }

#endif

#endif //LARIATDATAPRODUCTS_ConditionsSummary_h
