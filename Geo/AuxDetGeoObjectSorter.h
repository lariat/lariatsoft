////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTER_H
#define GEO_AUXDETGEOOBJECTSORTER_H

#include <vector>

#include "SimpleTypesAndConstants/geo_types.h"

#include "fhiclcpp/ParameterSet.h"

namespace geo{
  
  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  
  class AuxDetGeoObjectSorter {

  public:

    AuxDetGeoObjectSorter();
    ~AuxDetGeoObjectSorter();

    virtual void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)   const = 0;
    virtual void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)  const = 0;
    
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTER_H
