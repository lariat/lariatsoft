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
    
    void                     Initialize( AuxDetGeometryData_t& geodata ) override;
    void                     Uninitialize();
    

  private:
    
    geo::AuxDetGeoObjectSorterLArIAT  fSorter; ///< class to sort geo objects
  };


}
#endif // GEO_CHANNELMAPLARIATALG_H

