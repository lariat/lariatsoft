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
    
    bool  fBeamOn;                   ///< is the beam on or off
    float fMagnetCurrent;            ///< value of the current supplied to the magnet
    float fTPCHV;                    ///< High voltage of the field cage
    std::vector<float> fMWPCVoltage; ///< voltage supplied to MWPCs
    float fPMTHV;                    ///< HV for the PMT
    float fSiPMHV;                   ///< HV for the SiPM

#ifndef __GCCXML__
    
  public:
    
    ConditionsSummary(bool               const& beamOn,
                      float              const& magnetCurrent,
                      float              const& tpchv,
                      float              const& pmthv,
                      float              const& sipmhv,
                      std::vector<float> const& mwpcVoltages);
    
    bool  const& BeamOn()                        const;
    float const& MagnetCurrent()                 const;
    float const& TPCHV()                         const;
    float const& MWPCVoltage(size_t const& mwpc) const;
    float const& PMTHV()                         const;
    float const& SiPMHV()                        const;
    
#endif
    
  };
}

#ifndef __GCCXML__

inline bool  const& ldp::ConditionsSummary::BeamOn()                        const { return fBeamOn;        }
inline float const& ldp::ConditionsSummary::MagnetCurrent()                 const { return fMagnetCurrent; }
inline float const& ldp::ConditionsSummary::TPCHV()                         const { return fTPCHV;         }
inline float const& ldp::ConditionsSummary::PMTHV()                         const { return fPMTHV;         }
inline float const& ldp::ConditionsSummary::SiPMHV()                        const { return fSiPMHV;        }

#endif

#endif //LARIATDATAPRODUCTS_ConditionsSummary_h
