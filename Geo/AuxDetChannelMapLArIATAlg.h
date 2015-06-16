////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapLArIATAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETCHANNELMAPLARIATALG_H
#define GEO_AUXDETCHANNELMAPLARIATALG_H

#include <vector>
#include <set>

#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geo/AuxDetChannelMapAlg.h"
#include "Geo/AuxDetGeoObjectSorterLArIAT.h"
#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

namespace geo{

  class AuxDetChannelMapLArIATAlg : public AuxDetChannelMapAlg{

  public:

    AuxDetChannelMapLArIATAlg(fhicl::ParameterSet const& p);
    
    void      Initialize( AuxDetGeometryData_t& geodata ) override;
    void      Uninitialize();
    uint32_t  PositionToAuxDetChannel(double const worldLoc[3],
				      size_t      &ad,
				      size_t      &sv) const;
    const TVector3 AuxDetChannelToPosition(uint32_t                     const& channel,
					   std::string                  const& auxDetName,
					   std::vector<geo::AuxDetGeo*> const& auxDets) const;
    
  private:
    
    geo::AuxDetGeoObjectSorterLArIAT  fSorter; ///< class to sort geo objects
  };


}
#endif // GEO_CHANNELMAPLARIATALG_H

